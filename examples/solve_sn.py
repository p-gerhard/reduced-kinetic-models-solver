import sys, os
import numpy as np

sys.path.append(os.path.abspath(".."))
from modules.simulation import *

os.environ["PYOPENCL_NO_CACHE"] = "1"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"
os.environ["CUDA_CACHE_DISABLE"] = "1"


if __name__ == "__main__":

    input_fn = "../data/mesh/unit_square_nx128_ny128.msh"
    # input_fn = "../data/mesh/unit_cube_32.msh"

    r = 0.05
    v = (4.0 / 3.0) * np.pi * r * r * r
    parameters = {
        "model_type": "ordinates",
        "model_name": "2d_uniform",
        # "model_name": "3d_lebedev",
        # "model_order": 50,
        "model_order": 64,
        "tmax": 5,
        "cfl": 0.8,
        "use_muscl": True,
        "ocl_options" : ["-cl-fast-relaxed-math"],
        "sigma": 0,
        "alpha": 1, # 1 non-abs / 0 abs
        "beta": 0,  # 0 spec.   / 1 diff.
        "src_toff": 10,
        "src_x": 0.5,
        "src_y": 0.5,
        "src_z": 0.5,
        "src_r": r,
        "src_v": v,
    }

    simu = Simulation(parameters, mesh_filename=input_fn)
    simu.solve()
