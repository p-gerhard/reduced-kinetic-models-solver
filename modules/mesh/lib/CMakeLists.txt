cmake_minimum_required(VERSION 3.13)

pybind11_add_module(mesh_process SHARED src/pybind11_wrappers.cc 
                                        src/mesh_process.cc
                                        src/mesh_process.h)                    