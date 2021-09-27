cmake_minimum_required(VERSION 3.13)

pybind11_add_module(mesh_process SHARED src/pybind11_wrappers.cc 
                                        src/mesh_q4_process.cc
                                        src/mesh_q4_process.h
                                        src/mesh_h8_process.cc                   
                                        src/mesh_h8_process.h)                    