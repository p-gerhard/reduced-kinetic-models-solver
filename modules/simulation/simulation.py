#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from distutils import dir_util
import glob
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

import os
import shutil
import time

import matplotlib.pyplot as plt
import meshio
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl import _find_pyopencl_include_path

from modules.mesh import *
from .helpers import *
from .model import is_implemented, get_model_parameters


class Simulation:
    def __init__(self, parameters, mesh_filename=None, output_filename="output.xmf"):
        logging.basicConfig(
            format="[%(asctime)s] - %(levelname)s - %(message)s",
            datefmt="%m/%d/%Y %I:%M:%S %p",
        )

        is_missing = False

        # Check solver input parameter type
        ref_param = {
            "model_type": str,
            "model_name": str,
            "model_order": int,
            "tmax": float,
            "cfl": float,
            "use_muscl": bool,
        }

        for key in ref_param.keys():
            if key not in parameters.keys():
                logger.error(
                    "Required parameter {k} of type {t} is missing !".format(
                        k=key, t=ref_param[key]
                    )
                )
                is_missing = True

        if is_missing:
            exit()

        # Fill the class instance with solver parameters
        self.tmax = np.float32(parameters["tmax"])
        self.cfl = np.float32(parameters["cfl"])
        self.use_muscl = parameters["use_muscl"]

        if not (
            is_implemented(
                parameters["model_type"],
                parameters["model_name"],
                parameters["model_order"],
            )
        ):
            exit()

        self.model_type = parameters["model_type"]
        self.model_name = parameters["model_name"]
        self.model_order = np.int32(parameters["model_order"])

        # Fill the class instance model parameters
        self.m, self.dim, self.src_file, self.model_ocl_options = get_model_parameters(
            self.model_type, self.model_name, self.model_order
        )

        # Add mesh parameters
        self.mesh = MeshStructured(mesh_filename, dim=self.dim)
        self.ngrid = self.mesh.nb_cells
        self.dt = compute_dt(
            self.dim, self.cfl, self.mesh.dx, self.mesh.dy, self.mesh.dz
        )

        # Append extra model parameters to the dictionary of parameters
        self.parameters = parameters
        self.parameters.update(
            {
                "model_ocl_options": self.model_ocl_options,
                "dim": self.dim,
                "m": self.m,
                "src_file": self.src_file,
                "dt": self.dt,
                "dx": self.mesh.dx,
                "dy": self.mesh.dy,
                "dz": self.mesh.dz,
                "dz": self.dt,
                "ngrid": self.ngrid,
            }
        )

        self.base_filename = output_filename.split(".")[0]
        self.xmf_filename = "{}{}".format(self.base_filename, ".xmf")

        self.dtype = np.float32
        self.source = self.__ocl_process_source()

    def __ocl_process_source(self, print_src=False):

        with open(self.src_file, "r") as f:
            src = f.read()

        logger.info("{:<30}:".format("Simulation's parameters"))
        nb_tab = 5
        for k, v in sorted(self.parameters.items()):
            if is_num_type(v):
                if np.issubdtype(type(v), np.integer):
                    print(nb_tab * "\t" + "- {:<20} {:<12d}".format(k, v))
                if np.issubdtype(type(v), np.inexact):
                    print(nb_tab * "\t" + "- {:<20} {:<12.6f}".format(k, v))
                src = src.replace("_{}_".format(k), "({})".format(v))
            elif isinstance(v, list):
                print(nb_tab * "\t" + "- {:<20} {}".format(k, v))
            else:
                print(nb_tab * "\t" + "- {:<20} {:<12}".format(k, v))

        if print_src:
            print(src)

        return src

    def __ocl_build_and_setup(self, copy_kernel_to_ocl_path=True):
        current_file_path = os.path.dirname(os.path.realpath(__file__))
        kernel_base = os.path.abspath(
            os.path.join(current_file_path, os.path.join("..", "kernels"))
        )

        self.ocl_options = self.model_ocl_options + ["-cl-fast-relaxed-math"]
        
        if self.use_muscl:
            self.ocl_options.extend(["-D USE_MUSCL"])
        
        # Copy kernel source files to PyOpenCL's default kernel folder
        if copy_kernel_to_ocl_path:
            self.ocl_include_default_path = _find_pyopencl_include_path()
            self.dest_kernel_dir_list = build_subfolder_list(
                kernel_base, self.ocl_include_default_path
            )
            self.dest_kernel_dir_list = list(set(self.dest_kernel_dir_list))

            dir_util.copy_tree(kernel_base, self.ocl_include_default_path)

        else:
            self.ocl_options.extend(["-I", kernel_base])

        self.ocl_ctx = cl.create_some_context(interactive=True)
        self.ocl_prg = cl.Program(self.ocl_ctx, self.source)
        self.ocl_prg.build(options=self.ocl_options)

        # Remove copied files
        if copy_kernel_to_ocl_path:
            recursive_remove(self.dest_kernel_dir_list)

        self.ocl_prop = cl.command_queue_properties.PROFILING_ENABLE
        self.ocl_queue = cl.CommandQueue(self.ocl_ctx, properties=self.ocl_prop)

    def __ocl_allocate_devices_buffer(self):
        nb_cells = self.mesh.nb_cells
        buffer_size = self.m * nb_cells

        # OpenCL GPU buffers
        self.wn_gpu = cl_array.empty(self.ocl_queue, buffer_size, dtype=self.dtype)
        self.wnp1_gpu = cl_array.empty(self.ocl_queue, buffer_size, dtype=self.dtype)

        self.cells_center_gpu = cl_array.to_device(
            self.ocl_queue, self.mesh.cells_center
        )
        self.elem2elem_gpu = cl_array.to_device(self.ocl_queue, self.mesh.elem2elem)

        # Size of solution buffer
        gpu_memory = self.wn_gpu.nbytes
        gpu_memory += self.wnp1_gpu.nbytes
        gpu_memory += self.cells_center_gpu.nbytes
        gpu_memory += self.elem2elem_gpu.nbytes

        # Special case of Discrete Ordinates method
        if self.model_type == "ordinates":
            if self.dim == 2:
                self.nb_macro_field = 3
            else:
                self.nb_macro_field = 4

            self.wn_macro_gpu = cl_array.empty(
                self.ocl_queue, self.nb_macro_field * nb_cells, dtype=self.dtype
            )
            gpu_memory += self.wn_macro_gpu.nbytes

        logger.info("{:<30}: {:>6.8f} MB".format("GPU buffer's size", gpu_memory / 1e6))

    def __setup_exporter(self):

        try:
            os.mkdir(self.base_filename)
        except:
            pass

        self.xmf_filename = os.path.join(self.base_filename, self.xmf_filename)

        # Host buffer for memory copies
        if self.model_type == "ordinates":
            self.wn_macro_cpu = np.empty(
                self.nb_macro_field * self.mesh.nb_cells, dtype=self.dtype
            )

        self.wn_cpu = np.empty(self.m * self.mesh.nb_cells, dtype=self.dtype)

        # Buffer for time plotting
        self.t_tab = []
        self.w0_tot = []

        try:
            os.remove(self.xmf_filename)
        except OSError:
            pass

        try:
            os.remove(self.xmf_filename.split(".xmf")[0] + ".h5")
        except OSError:
            pass

    def __export_data(self, xdmf_writer):
        nc = self.mesh.nb_cells
        if self.model_type == "ordinates":
            self.ocl_prg.sn_compute_macro(
                self.ocl_queue,
                (self.mesh.nb_cells,),
                None,
                self.wn_gpu.data,
                self.wn_macro_gpu.data,
            ).wait()

            w = cl_array.sum(self.wn_macro_gpu[0:nc]).get()

            cd = {
                "w0": {self.mesh.cell_name: self.wn_macro_cpu[0 * nc : 1 * nc]},
                "w1": {self.mesh.cell_name: self.wn_macro_cpu[1 * nc : 2 * nc]},
                "w2": {self.mesh.cell_name: self.wn_macro_cpu[2 * nc : 3 * nc]},
            }

            if self.dim == 3:
                cd["w3"] = {self.mesh.cell_name: self.wn_macro_cpu[3 * nc : 4 * nc]}

            cl.enqueue_copy(
                self.ocl_queue, self.wn_macro_cpu, self.wn_macro_gpu.data
            ).wait()
        else:

            w = cl_array.sum(self.wn_gpu[0:nc]).get()

            cd = {
                "w0": {self.mesh.cell_name: self.wn_macro_cpu[0 * nc : 1 * nc]},
                "w1": {self.mesh.cell_name: self.wn_macro_cpu[1 * nc : 2 * nc]},
                "w2": {self.mesh.cell_name: self.wn_macro_cpu[2 * nc : 3 * nc]},
            }

            if self.dim == 3:
                cd["w3"] = {self.mesh.cell_name: self.wn_macro_cpu[3 * nc : 4 * nc]}

            cl.enqueue_copy(self.ocl_queue, self.wn_cpu, self.wn_gpu.data).wait()

        self.w0_tot.append(w)
        self.t_tab.append(self.t)

        xdmf_writer.write_data(self.t, cell_data=cd)

    def solve(self):
        self.__ocl_build_and_setup()

        self.__ocl_allocate_devices_buffer()
        self.__setup_exporter()
        self.t = 0
        ite = 0

        logger.info("{:<30}".format("OpenCL computations started"))

        t1 = time.time()
        # Init solution
        self.ocl_prg.vf_init_sol(
            self.ocl_queue,
            (self.mesh.nb_cells,),
            None,
            self.cells_center_gpu.data,
            self.wn_gpu.data,
        ).wait()

        # Setup kernel arguments
        kernel_euler_time_step = self.ocl_prg.solver_muscl_euler_time_step
        kernel_euler_time_step.set_scalar_arg_dtypes(
            [self.dtype, None, None, None, None]
        )

        with meshio.xdmf.TimeSeriesWriter(self.xmf_filename) as xdmf_writer:
            xdmf_writer.write_points_cells(
                self.mesh.nodes, [(self.mesh.cell_name, self.mesh.cells)]
            )
            while self.t <= self.tmax:
                if ite % 20 == 0:
                    self.__export_data(xdmf_writer)

                kernel_euler_time_step(
                    self.ocl_queue,
                    (self.mesh.nb_cells,),
                    None,
                    self.t,
                    self.cells_center_gpu.data,
                    self.elem2elem_gpu.data,
                    self.wn_gpu.data,
                    self.wnp1_gpu.data,
                ).wait()

                self.wn_gpu, self.wnp1_gpu = self.wnp1_gpu, self.wn_gpu
                self.t += self.dt

                print(get_ite_title(ite, self.t, 0), end="\r")
                ite += 1
        t2 = time.time()

        logger.info("{:<30}: {:>6.8f} sec.".format("Computation done in", t2 - t1))

        w0_tot_filename = os.path.join(
            self.base_filename, "{}{}".format(self.base_filename, "_w0_tot.dat")
        )

        np.savetxt(w0_tot_filename, np.stack((self.t_tab, self.w0_tot), axis=1))

        logger.info(
            "{:<30}: {:<30}".format(
                "Fields saved in", os.path.abspath(self.xmf_filename)
            )
        )

        logger.info(
            "{:<30}: {:<30}".format(
                "Total density saved in", os.path.abspath(w0_tot_filename)
            )
        )

        plt.plot(self.t_tab, self.w0_tot)
        plt.show()
