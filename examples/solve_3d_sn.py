import sys, os

import numpy as np


sys.path.append(os.path.abspath(".."))
from modules.simulation import *

os.environ["PYOPENCL_NO_CACHE"] = "1"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"
os.environ["CUDA_CACHE_DISABLE"] = "1"


if __name__ == "__main__":

    input_fn = "../data/mesh/unit_cube_32.msh"
    output_fn = "sn_3d.xjmf"

    r = 0.05
    v = (4.0 / 3.0) * np.pi * r * r * r
    parameters = {
        "model": "ordinates",
        "m": 6,
        "dim": 3,
        "tmax": 1.0,
        "cfl": 0.05,
        "sigma": 0,
        "alpha": 1,  # 1 non-abs / 0 abs
        "beta": 0.0,  # 0 spec.   / 1 diff.
        "src_toff": 10,
        "src_x": 0.5,
        "src_y": 0.5,
        "src_z": 0.5,
        "src_r": r,
        "src_v": v,
        "src_file": "sn_3d_main.c",
    }

    simu = Simulation(parameters, mesh_filename=input_fn, output_filename=output_fn)
    simu.solve()