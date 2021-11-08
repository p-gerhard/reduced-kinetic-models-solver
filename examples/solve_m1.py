import sys, os
import numpy as np

sys.path.append(os.path.abspath(".."))
from modules.simulation import *

os.environ["PYOPENCL_NO_CACHE"] = "1"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"
os.environ["CUDA_CACHE_DISABLE"] = "1"

if __name__ == "__main__":

    # input_fn = "../data/mesh/unit_cube_64.msh"    
    # input_fn = "../data/mesh/3d_star_xy_0125.msh"
    input_fn = "../data/mesh/salle_3.msh"


    r = 0.125
    v = (4.0 / 3.0) * np.pi * r * r * r
    parameters = {
        "model_type": "moments",
        "model_name": "3d_s2",
        "model_order": 1,
        "tmax": 40,
        "cfl": 0.8,
        "use_muscl": False,
        "ocl_options": [
            "-cl-fast-relaxed-math",
            # "-D USE_KINETIC_NUM_FLUX",
            "-D USE_QUAD_LEBEDEV",
        ],
        "sigma": 0,
        "alpha": 1,  # 1 non-abs / 0 abs
        "beta": 1,  # 0 spec.   / 1 diff.
        "src_toff": 10,
        "src_x": 0.5,
        "src_y": 0.5,
        "src_z": 0.5,
        "src_r": r,
        "src_v": v,
    }

    simu = Simulation(parameters, mesh_filename=input_fn)
    simu.solve()
