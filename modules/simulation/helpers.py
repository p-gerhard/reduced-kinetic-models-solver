#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

import json
import os
import shutil
import numbers
import time


import numpy as np






def check_postprocess_parameter(postprocess_parameters):
    if "time" in postprocess_parameters:
        assert all(
            isinstance(t, numbers.Number) for t in postprocess_parameters["time"]
        )

    if "line" in postprocess_parameters:
        for line in postprocess_parameters["line"]:
            assert "point_1" in line
            point_1 = line["point_1"]
            assert isinstance(point_1, list)

            len_pt1 = len(point_1)
            assert len_pt1 == 3 or len_pt1 == 2
            assert all(isinstance(coord, numbers.Number) for coord in point_1)

            assert "point_2" in line
            point_2 = line["point_2"]
            assert isinstance(point_2, list)

            len_pt2 = len(point_2)
            assert len_pt2 == 3 or len_pt2 == 2
            assert all(isinstance(coord, numbers.Number) for coord in point_2)

        if "slice" in postprocess_parameters:
            for slice in postprocess_parameters["slice"]:
                assert "point" in slice
                point = slice["point"]
                assert isinstance(point, list)

                assert len(point) == 3
                assert all(isinstance(coord, numbers.Number) for coord in point)

                assert "normal" in slice
                normal = slice["normal"]
                assert isinstance(normal, list)

                assert len(normal) == 3
                assert all(isinstance(coord, numbers.Number) for coord in normal)


def dump_dic_to_json(dic, filename):
    # Cast dic data to serializable type
    dumpable_dic = dic.copy()
    for key, value in dumpable_dic.items():
        if isinstance(value, np.floating):
            dumpable_dic[key] = float(dumpable_dic[key])

        if isinstance(value, np.integer):
            dumpable_dic[key] = int(dumpable_dic[key])

    with open(filename, "w") as file:
        json.dump(dumpable_dic, file, indent=4, sort_keys=True)


def pprint_dict(dict, indent=0):
    indent = 5
    for k, v in sorted(dict.items()):
        if is_num_type(v):
            if np.issubdtype(type(v), np.integer):
                print(indent * "\t" + "- {:<20} {:<12d}".format(k, v))
            if np.issubdtype(type(v), np.inexact):
                print(indent * "\t" + "- {:<20} {:<12.6f}".format(k, v))

        else:
            print(indent * "\t" + "- {:<20} {}".format(k, v))


def make_folder_empty(folder):
    try:
        os.mkdir(folder)
    except:
        pass

    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
        except Exception as e:
            print("Failed to delete {}. Reason: {}".format(file_path, e))


def check_inclusion_and_type(input_param, ref_param, msg_header=""):
    is_missing = False

    for key, value in ref_param.items():
        # Check inclusion
        if key not in input_param.keys():
            if msg_header != "":
                error_msg = '{} "{} is missing'.format(msg_header, key)
            else:
                error_msg = 'Parameter "{} is missing'.format(msg_header, key)

            logger.error(error_msg)
            is_missing = True

        # Check type and try to cast if possible
        else:
            if not isinstance(input_param[key], ref_param[key]):
                try:
                    input_param[key] = ref_param[key](input_param[key])
                except:
                    if msg_header != "":
                        logger.error(
                            '{msg} "{k}" of type {t_in} cannot be cast into {t_ref} !'.format(
                                msg=msg_header,
                                k=key,
                                t_in=type(input_param[key]),
                                t_ref=ref_param[key],
                            )
                        )
                    else:
                        logger.error(
                            'Parameter "{k}" of type {t_in} into {t_ref}'.format(
                                k=key, t_in=type(input_param[key]), t_ref=ref_param[key]
                            )
                        )
                    is_missing = True

    if is_missing:
        logger.error("Required values are : {}".format(list(ref_param)))

    else:
        if msg_header != "":
            info_msg = "{:<30}: OK".format("{} ".format(msg_header))
        else:
            info_msg = "{:<30}: OK".format("Parameters")

        logger.info(info_msg)

    return not is_missing


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
                error_msg = '{} "{} is missing'.format(msg_header, p)
            else:
                error_msg = 'Parameter "{} is missing'.format(msg_header, p)

            logger.error(error_msg)
            is_missing = True

    if is_missing:
        logger.error("Required values are : {}".format(list(ref_param)))
    else:
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
