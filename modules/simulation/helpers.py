#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import logging
import os
import shutil
import time

import numpy as np

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
