import os
import numpy as np

from .helpers import check_inclusion

ORDINATES_MODELS = {
    "2d_uniform": {
        "dim": 2,
        "order": [4, 8, 16, 32, 64, 128, 256, 512],
        "src_file": "sn/sn.cl",
        "ocl_options": ["-D USE_QUAD_UNIFORM"],
    },
    "3d_lebedev": {
        "dim": 3,
        "order": [
            6,
            14,
            26,
            38,
            50,
            74,
            86,
            110,
            146,
            170,
            194,
            230,
            266,
            302,
            350,
            434,
            590,
            770,
        ],
        "src_file": "sn/sn.cl",
        "ocl_options": ["-D USE_QUAD_LEBEDEV"],
    },
}

MOMENTS_MODELS = {
    "2d_s1": {
        "type": "moments",
        "dim": 2,
        "order": [1],
        "src_file": "m1/m1.cl",
        "ocl_options": ["-D USE_M1_S1_CLOSURE", "-D IS_MOMENT_MODEL"],
    },
    "3d_s2": {
        "type": "moments",
        "dim": 3,
        "order": [1],
        "src_file": "m1/m1.cl",
        "ocl_options": ["-D USE_M1_S2_CLOSURE", "-D IS_MOMENT_MODEL"],
    },
}

IMPLEMENTED_MODELS = {"ordinates": ORDINATES_MODELS, "moments": MOMENTS_MODELS}


def get_model_parameters(model_type, model_name, model_order):

    ocl_options = IMPLEMENTED_MODELS[model_type][model_name]["ocl_options"]
    dim = IMPLEMENTED_MODELS[model_type][model_name]["dim"]
    src_file = IMPLEMENTED_MODELS[model_type][model_name]["src_file"]
    src_file = os.path.abspath(os.path.join("../modules/kernels/", src_file))
    nb_macro_to_reconstruct = 0

    if model_type == "ordinates":
        m = model_order
        if dim == 2:
            nb_macro_to_reconstruct = 3

        else:
            nb_macro_to_reconstruct = 4

    if model_type == "moments":
        if dim == 2:
            m = 3
        else:
            m = 4

    return (
        np.int32(m),
        np.int32(nb_macro_to_reconstruct),
        np.int32(dim),
        src_file,
        ocl_options,
    )


def get_model_ocl_options(model_type, model_name):
    return IMPLEMENTED_MODELS[model_type][model_name]["model_ocl_options"]


def is_implemented(model_type, model_name, model_order):
    found = True

    # Check if input model type is referenced as implemented
    found = found and check_inclusion(
        model_type, IMPLEMENTED_MODELS.keys(), msg_header="Model type"
    )

    # Check if input model name is referenced as implemented
    found = found and check_inclusion(
        model_name, IMPLEMENTED_MODELS[model_type].keys(), msg_header="Model name"
    )

    # Check if input model order is referenced as implemented
    found = found and check_inclusion(
        model_order,
        IMPLEMENTED_MODELS[model_type][model_name]["order"],
        msg_header="Model order",
    )

    return found
