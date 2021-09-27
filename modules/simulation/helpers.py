#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import logging
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

import os
import shutil
import time

import numpy as np


def compute_dt(dim, cfl, dx, dy, dz):

    if dim == 2:
        dt = cfl * (dx * dy) / (2.0 * dx + 2.0 * dy)

    if dim == 3:
        dt = cfl * (dx * dy * dz) / (2.0 * (dx * dy + dy * dz + dz * dx))

    return dt

def check_inclusion(input_param, ref_param, msg_header=""):
    is_missing = False

    # To make it iterable
    if isinstance(input_param, (str, int, float, complex, bool)):
        input_param = [input_param]

    for p in input_param:
        if p not in ref_param:
            if msg_header != "":
                error_msg = '{} "{} is missing !'.format(msg_header, p)
            else:
                error_msg = 'Parameter "{} is missing !'.format(msg_header, p)

            logger.error(error_msg)
            logger.error("Required values are : {}".format(list(ref_param)))
            is_missing = True

    if not is_missing:
        if msg_header != "":
            info_msg = "{:<30}: OK".format("{} ".format(msg_header))
        else:
            info_msg = "{:<30}: OK".format("Parameters")

        logger.info(info_msg)

    return not is_missing


def build_subfolder_list(base_folder, dest_folder=None):
    subfolder_list = []
    for r, d, f in os.walk(base_folder):
        for file in f:
            subfolder = os.path.relpath(r, base_folder)

            if dest_folder is not None:
                subfolder = os.path.join(dest_folder, subfolder)

            subfolder_list.append(subfolder)

    return list(set(subfolder_list))


def recursive_remove(dir_list):
    for dir in dir_list:
        shutil.rmtree(dir, ignore_errors=False, onerror=None)


def get_ite_title(ite, t, elapsed):
    return "ite = {}, t = {:f}, elapsed (s) = {:f}".format(ite, t, elapsed)


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

    if type_to_check == list:
        if isinstance(val, list):
            return val
        else:
            raise TypeError("{} must be list".format(key))

    if type_to_check == bool:
        if isinstance(val, bool):
            return val
        else:
            raise TypeError("{} must be list".format(key))
