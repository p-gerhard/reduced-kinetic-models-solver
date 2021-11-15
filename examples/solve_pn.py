import sys, os
import numpy as np

sys.path.append(os.path.abspath(".."))
from modules.simulation import *

os.environ["PYOPENCL_NO_CACHE"] = "1"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"
os.environ["CUDA_CACHE_DISABLE"] = "1"


if __name__ == "__main__":

    # input_fn = "../data/mesh/unit_square_nx128_ny128.msh"
    # input_fn = "../data/mesh/unit_cube_64.msh"
    input_fn = "../data/mesh/salle_3.msh"

    r = 0.125
    v = (4.0 / 3.0) * np.pi * r * r * r
    
    simulation_parameters = {
        "model_type": "spherical_harmonics",
        "model_name": "3d",
        "model_order": 7,
        "tmax": 150,
        "cfl": 1,
        "use_muscl": True,
        "ocl_options" : ["-cl-fast-relaxed-math"],
        "sigma": 0,
        "alpha": 0.95, # 1 non-abs / 0 abs
        "beta": 1,  # 0 spec.   / 1 diff.
        "src_toff": 60,
        "src_x": 0.5,
        "src_y": 0.5,
        "src_z": 0.5,
        "src_r": r,
        "src_v": v,
    }

    postprocess_parameters = {
        "type": "paraview",
        "time": [0, 10, 150],
        "slice": [{"point": [0.5, 0.5, 0.5], "normal": [0, 0, 1]}],
        "line": [{"point_1": [0, 0.5, 0.5], "point_2": [10, 0.5, 0.5]}],
        "build_tex": False,
    }

    simu = Simulation(
        simulation_parameters,
        mesh_filename=input_fn,
        postprocess_parameters=postprocess_parameters,
    )

    simu.solve()