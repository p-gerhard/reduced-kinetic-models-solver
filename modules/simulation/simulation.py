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

        self.__check_required_input_parameters(parameters)

        # Fill the class instance with solver parameters
        self.tmax = np.float32(parameters["tmax"])
        self.cfl = np.float32(parameters["cfl"])
        self.use_muscl = parameters["use_muscl"]
        self.model_type = parameters["model_type"]
        self.model_name = parameters["model_name"]
        self.model_order = np.int32(parameters["model_order"])

        # Fill the object with registered model extra parameters
        (
            self.m,
            self.nb_macro_to_reconstruct,
            self.dim,
            self.src_file,
            self.model_ocl_options,
        ) = get_model_parameters(self.model_type, self.model_name, self.model_order)

        # Deal with OpenCL compilation options
        self.ocl_options = self.__ocl_build_options(parameters)

        # Read the mesh and create the connectivity
        self.mesh = MeshStructured(mesh_filename, dim=self.dim)

        # Compute dt using CFL
        self.dt = compute_dt(
            self.dim, self.cfl, self.mesh.dx, self.mesh.dy, self.mesh.dz
        )

        self.__setup_output_folder()

        # Fill global simulation dictionary with all the parameters
        self.parameters = parameters
        self.parameters.update(
            {
                "ocl_options": self.ocl_options,
                "dim": self.dim,
                "m": self.m,
                "nb_macro_to_reconstruct": self.nb_macro_to_reconstruct,
                "src_file": self.src_file,
                "dt": self.dt,
                "dx": self.mesh.dx,
                "dy": self.mesh.dy,
                "dz": self.mesh.dz,
                "ngrid": self.mesh.nb_cells,
                "xmf_filename": self.xmf_filename,
            }
        )

        # Display all simulation parameter
        logger.info("{:<30}:".format("Simulation's parameters"))

        pprint_dict(self.parameters, indent=5)

        self.source = self.__ocl_process_source()

    def __check_required_input_parameters(self, parameters):

        ref_param = {
            "model_type": str,
            "model_name": str,
            "model_order": int,
            "ocl_options": list,
            "tmax": float,
            "cfl": float,
            "use_muscl": bool,
        }
        # Check solver required input parameters and type
        if not (
            check_inclusion_and_type(
                parameters, ref_param, msg_header="Required solver parameters"
            )
        ):
            exit()

        # Check if model in implemented
        if not (
            is_implemented(
                parameters["model_type"],
                parameters["model_name"],
                parameters["model_order"],
            )
        ):
            exit()

    def __setup_output_folder(self):

        time_str = time.strftime("%Y%m%d-%H%M%S")
        self.output_folder = "{}_{}_order{}".format(
            self.model_type, self.model_name, self.model_order
        )

        fn = "{}{}".format("results", ".xmf")
        self.xmf_filename = os.path.join(self.output_folder, fn)

        fn = "{}{}".format("density_tot", ".dat")
        self.density_tot_filename = os.path.join(self.output_folder, fn)

        make_folder_empty(self.output_folder)

    def __ocl_build_options(self, parameters):
        ocl_options = []

        if self.use_muscl:
            ocl_options.extend(["-D USE_MUSCL"])

        ocl_options.extend(self.model_ocl_options)

        if "ocl_options" in parameters.keys():
            ocl_options.extend(parameters["ocl_options"])

        return ocl_options

    def __ocl_process_source(self, print_src=False):

        with open(self.src_file, "r") as f:
            src = f.read()

        # Replace "_param_" in #define PARAM _param_ inside OpenCL source
        for k, v in sorted(self.parameters.items()):
            if is_num_type(v):
                src = src.replace("_{}_".format(k), "({})".format(v))

        if print_src:
            print(src)

        return src

    def __ocl_build_and_setup(self, copy_kernel_to_ocl_path=True):
        current_file_path = os.path.dirname(os.path.realpath(__file__))
        kernel_base = os.path.abspath(
            os.path.join(current_file_path, os.path.join("..", "kernels"))
        )

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
        buffer_size = self.m * self.mesh.nb_cells

        # OpenCL GPU solution buffers
        self.wn_device = cl_array.empty(self.ocl_queue, buffer_size, dtype=np.float32)
        self.wnp1_device = cl_array.empty(self.ocl_queue, buffer_size, dtype=np.float32)

        # OpenCL cells coordinates buffer
        self.cells_center_device = cl_array.to_device(
            self.ocl_queue, self.mesh.cells_center
        )

        # OpenCL cells connectivity buffer
        self.elem2elem_device = cl_array.to_device(self.ocl_queue, self.mesh.elem2elem)

        # Size estimation of used GPU buffer
        gpu_memory = self.wn_device.nbytes
        gpu_memory += self.wnp1_device.nbytes
        gpu_memory += self.cells_center_device.nbytes
        gpu_memory += self.elem2elem_device.nbytes

        if self.nb_macro_to_reconstruct > 0:
            self.wn_macro_device = cl_array.empty(
                self.ocl_queue,
                self.nb_macro_to_reconstruct * self.mesh.nb_cells,
                dtype=np.float32,
            )
            gpu_memory += self.wn_macro_device.nbytes

        logger.info(
            "{:<30}: {:>6.8f} MB".format("Devices buffer's size", gpu_memory / 1e6)
        )

    def __allocate_host_buffer(self):
        # Host buffer for memory copies
        if self.nb_macro_to_reconstruct > 0:
            buffer_size = self.nb_macro_to_reconstruct * self.mesh.nb_cells
            self.wn_macro_host = np.empty(buffer_size, dtype=np.float32)

        # Host solution (for plotting)
        self.wn_host = np.empty(self.m * self.mesh.nb_cells, dtype=np.float32)

        # Total density and time (abscissa)
        self.t_tab = []
        self.w0_tot = []

    def __export_data(self, xdmf_writer):
        nc = self.mesh.nb_cells

        if self.nb_macro_to_reconstruct > 0:
            self.ocl_prg.sn_compute_macro(
                self.ocl_queue,
                (self.mesh.nb_cells,),
                None,
                self.wn_device.data,
                self.wn_macro_device.data,
            ).wait()

            w = cl_array.sum(self.wn_macro_device[0:nc]).get()

            cd = {
                "w0": {self.mesh.cell_name: self.wn_macro_host[0 * nc : 1 * nc]},
                "w1": {self.mesh.cell_name: self.wn_macro_host[1 * nc : 2 * nc]},
                "w2": {self.mesh.cell_name: self.wn_macro_host[2 * nc : 3 * nc]},
            }

            if self.dim == 3:
                cd["w3"] = {self.mesh.cell_name: self.wn_macro_host[3 * nc : 4 * nc]}

            cl.enqueue_copy(
                self.ocl_queue, self.wn_macro_host, self.wn_macro_device.data
            ).wait()
        else:

            w = cl_array.sum(self.wn_device[0:nc]).get()

            cd = {
                "w0": {self.mesh.cell_name: self.wn_host[0 * nc : 1 * nc]},
                "w1": {self.mesh.cell_name: self.wn_host[1 * nc : 2 * nc]},
                "w2": {self.mesh.cell_name: self.wn_host[2 * nc : 3 * nc]},
            }

            if self.dim == 3:
                cd["w3"] = {self.mesh.cell_name: self.wn_host[3 * nc : 4 * nc]}

            cl.enqueue_copy(self.ocl_queue, self.wn_host, self.wn_device.data).wait()

        self.w0_tot.append(w)
        self.t_tab.append(self.t)

        xdmf_writer.write_data(self.t, cell_data=cd)

    def solve(self):
        self.__ocl_build_and_setup()
        self.__ocl_allocate_devices_buffer()
        self.__allocate_host_buffer()

        self.t = 0
        ite = 0

        logger.info("{:<30}".format("OpenCL computations started"))

        t1 = time.time()

        # Init solution
        self.ocl_prg.vf_init_sol(
            self.ocl_queue,
            (self.mesh.nb_cells,),
            None,
            self.cells_center_device.data,
            self.wn_device.data,
        ).wait()

        # Setup kernel arguments
        kernel_euler_time_step = self.ocl_prg.solver_muscl_euler_time_step
        kernel_euler_time_step.set_scalar_arg_dtypes(
            [np.float32, None, None, None, None]
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
                    self.cells_center_device.data,
                    self.elem2elem_device.data,
                    self.wn_device.data,
                    self.wnp1_device.data,
                ).wait()

                self.wn_device, self.wnp1_device = self.wnp1_device, self.wn_device
                self.t += self.dt

                print(get_ite_title(ite, self.t, 0), end="\r")
                ite += 1
        t2 = time.time()

        logger.info("{:<30}: {:>6.8f} sec.".format("Computation done in", t2 - t1))

        np.savetxt(
            self.density_tot_filename, np.stack((self.t_tab, self.w0_tot), axis=1)
        )

        logger.info(
            "{:<30}: {:<30}".format(
                "Fields saved in", os.path.abspath(self.xmf_filename)
            )
        )

        logger.info(
            "{:<30}: {:<30}".format(
                "Total density saved in", os.path.abspath(self.density_tot_filename)
            )
        )

        plt.plot(self.t_tab, self.w0_tot)
        plt.show()
