#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from distutils import dir_util
import glob
import logging
import os
import shutil
import time


import matplotlib.pyplot as plt
import meshio
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl import _find_pyopencl_include_path, _DEFAULT_INCLUDE_OPTIONS

from modules.mesh import *


def _build_subfolder_list(base_folder, dest_folder=None):
    subfolder_list = []
    for r, d, f in os.walk(base_folder):
        for file in f:
            subfolder = os.path.relpath(r, base_folder)

            if dest_folder is not None:
                subfolder = os.path.join(dest_folder, subfolder)

            subfolder_list.append(subfolder)

    return list(set(subfolder_list))


def _recursive_remove(dir_list):
    for dir in dir_list:
        shutil.rmtree(dir, ignore_errors=False, onerror=None)


def get_ite_title(ite, t, elapsed):
    return "ite = {}, t = {:f}, elapsed (s) = {:f}".format(ite, t, elapsed)


def check_parameters(ref_param, input_param):
    is_missing = False
    for key in ref_param.keys():
        if key not in input_param.keys():
            print(
                "[Error] - Required parameter {k} of type {t} is missing".format(
                    k=key, t=ref_param[key]
                )
            )
            is_missing = True

    if is_missing:
        exit()


def is_num_type(val):
    return np.issubdtype(type(val), np.number)


def safe_assign(key, val, type_to_check):

    if type_to_check == int:
        if is_num_type(val):
            return int(val)
        else:
            raise TypeError("{} must be numeric".format(key))

    if type_to_check == float:
        if is_num_type(val):
            return float(val)
        else:
            raise TypeError("{} must be numeric".format(key))

    if type_to_check == str:
        if isinstance(val, str):
            return val
        else:
            raise TypeError("{} must be str".format(key))


class Simulation:
    def __init__(self, parameters, mesh_filename=None, output_filename="output.xmf"):
        logging.basicConfig(
            format="[%(asctime)s] - %(levelname)s - %(message)s",
            level=logging.INFO,
            datefmt="%m/%d/%Y %I:%M:%S %p",
        )
        self.__set_req_parameters(parameters)

        self.mesh = MeshStructured(mesh_filename, dim=3)

        self.base_filename = output_filename.split(".")[0]
        self.xmf_filename = "{}{}".format(self.base_filename, ".xmf")

        self.parameters = parameters

        self.__update_parameters_with_mesh_data()

        self.dtype = np.float32
        self.source = self.__ocl_process_source()

    def __set_req_parameters(self, parameters, input_mesh=False):

        # Check if required parameters are provided
        ref_param = {
            "dim": int,
            "m": int,
            "tmax": float,
            "cfl": float,
            "src_file": str,
        }

        check_parameters(ref_param, parameters)

        self.dim = safe_assign("dim", parameters["dim"], ref_param["dim"])
        self.m = safe_assign("m", parameters["m"], ref_param["m"])
        self.tmax = safe_assign("tmax", parameters["tmax"], ref_param["tmax"])
        self.cfl = safe_assign("cfl", parameters["cfl"], ref_param["cfl"])

        self.src_file = safe_assign(
            "src_file", parameters["src_file"], ref_param["src_file"]
        )

    def __update_parameters_with_mesh_data(self):
        self.ngrid = self.mesh.nb_cells
        if self.dim == 2:
            self.dt = (
                self.cfl
                * (self.mesh.dx * self.mesh.dy)
                / (2.0 * self.mesh.dx + 2.0 * self.mesh.dy)
            )

            dz = 1

        if self.dim == 3:
            self.dt = (
                self.cfl
                * (self.mesh.dx * self.mesh.dy * self.mesh.dz)
                / (
                    2.0 * self.mesh.dx * self.mesh.dy
                    + 2.0 * self.mesh.dy * self.mesh.dz
                    + 2.0 * self.mesh.dz * self.mesh.dx
                )
            )
            dz = self.mesh.dz

        self.parameters.update(
            {
                "dx": np.float32(self.mesh.dx),
                "dy": np.float32(self.mesh.dy),
                "dz": np.float32(dz),
                "dt": self.dt,
                "ngrid": self.mesh.nb_cells,
            }
        )

    def __ocl_process_source(self, print_src=False):

        with open(self.src_file, "r") as f:
            src = f.read()

        print("[info]- simulation's parameters:")

        for k, v in self.parameters.items():
            if is_num_type(v):
                if np.issubdtype(type(v), np.integer):
                    print("\t- {:<12}  {:<12d}".format(k, v))
                if np.issubdtype(type(v), np.inexact):
                    print("\t- {:<12}  {:<12.6f}".format(k, v))

                src = src.replace("_{}_".format(k), "({})".format(v))

        if print_src:
            print(src)

        return src

    def __ocl_build_and_setup(self, copy_kernel_to_ocl_path=True):
        current_file_path = os.path.dirname(os.path.realpath(__file__))
        kernel_base = os.path.abspath(
            os.path.join(current_file_path, os.path.join("..", "kernels"))
        )
        self.ocl_options = ["-cl-fast-relaxed-math"]

        if copy_kernel_to_ocl_path:
            self.ocl_include_default_path = _find_pyopencl_include_path()
            self.dest_kernel_dir_list = _build_subfolder_list(
                kernel_base, self.ocl_include_default_path
            )
            self.dest_kernel_dir_list = list(set(self.dest_kernel_dir_list))

            # shutil.copytree(
            #     kernel_base, self.ocl_include_default_path, dirs_exist_ok=True
            # )

            dir_util.copy_tree(kernel_base, self.ocl_include_default_path)

        else :
            self.ocl_options.extend(["-I", kernel_base])

        self.ocl_ctx = cl.create_some_context(interactive=True)
        self.ocl_prg = cl.Program(self.ocl_ctx, self.source)
        self.ocl_prg.build(options=self.ocl_options)

        if copy_kernel_to_ocl_path:
            _recursive_remove(self.dest_kernel_dir_list)

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

        logging.info(
            "{:<30}: {:>6.8f} MB".format("GPU buffer's size", gpu_memory / 1e6)
        )

    def __clean(self):
        _recursive_remove(self.dest_kernel_dir_list)

    def __allocate_cpu_buffer(self):
        pass

    def __setup_exporter(self):

        try:
            os.mkdir(self.base_filename)
        except:
            pass

        self.xmf_filename = os.path.join(self.base_filename, self.xmf_filename)

        # Host buffer for memory copies
        self.wn_cpu = np.empty(self.m * self.mesh.nb_cells, dtype=self.dtype)

        # Buffer for time plotting
        self.t_tab = []
        self.w0_tot = []
        nc = self.mesh.nb_cells

        self.cell_data = {
            "w0": {self.mesh.cell_name: self.wn_cpu[0 * nc : 1 * nc]},
            "w1": {self.mesh.cell_name: self.wn_cpu[1 * nc : 2 * nc]},
            "w2": {self.mesh.cell_name: self.wn_cpu[2 * nc : 3 * nc]},
            "w3": {self.mesh.cell_name: self.wn_cpu[3 * nc : 4 * nc]},
        }

        try:
            os.remove(self.xmf_filename)
        except OSError:
            pass

        try:
            os.remove(self.xmf_filename.split(".xmf")[0] + ".h5")
        except OSError:
            pass

    def __write_data(self, xdmf_writer):
        nc = self.mesh.nb_cells
        cd = {
            "w0": {self.mesh.cell_name: self.wn_cpu[0 * nc : 1 * nc]},
            "w1": {self.mesh.cell_name: self.wn_cpu[1 * nc : 2 * nc]},
            "w2": {self.mesh.cell_name: self.wn_cpu[2 * nc : 3 * nc]},
            "w3": {self.mesh.cell_name: self.wn_cpu[3 * nc : 4 * nc]},
        }

        xdmf_writer.write_data(self.t, cell_data=cd)

    def solve(self):
        self.__ocl_build_and_setup()

        self.__ocl_allocate_devices_buffer()
        self.__setup_exporter()
        self.t = 0
        ite = 0

        logging.info("{:<30}".format("OpenCL computations started"))

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

        with meshio.xdmf.TimeSeriesWriter(self.xmf_filename) as xmdf_writer:
            xmdf_writer.write_points_cells(
                self.mesh.nodes, [(self.mesh.cell_name, self.mesh.cells)]
            )
            while self.t < self.tmax:
                if ite % 20 == 0:
                    w = cl_array.sum(self.wn_gpu[0 : self.mesh.nb_cells]).get()

                    # Warning in Pn we need to multiply by sqrt(4 * pi) to get w
                    self.w0_tot.append(w)
                    self.t_tab.append(self.t)

                    cl.enqueue_copy(
                        self.ocl_queue, self.wn_cpu, self.wn_gpu.data
                    ).wait()

                    self.__write_data(xmdf_writer)

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

        logging.info("{:<30}: {:>6.8f} sec.".format("Computation done in", t2 - t1))

        w0_tot_filename = os.path.join(
            self.base_filename, "{}{}".format(self.base_filename, "_w0_tot.dat")
        )

        np.savetxt(w0_tot_filename, np.stack((self.t_tab, self.w0_tot), axis=1))

        logging.info(
            "{:<30}: {:<30}".format(
                "Fields saved in", os.path.abspath(self.xmf_filename)
            )
        )

        logging.info(
            "{:<30}: {:<30}".format(
                "Total density saved in", os.path.abspath(w0_tot_filename)
            )
        )

        plt.plot(self.t_tab, self.w0_tot)
        plt.show()
