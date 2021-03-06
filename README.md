# Reduced-kinetic-models-solver
### Description


### How to Build and Install
This project is written in Python3, OpenCLv1.2 and C++11. It mainly relies on the packages [meshio](https://github.com/nschloe/meshio "meshio"), [pybind11](https://github.com/pybind/pybind11 "pybind11"), [pyopencl](https://github.com/inducer/pyopencl "PyOpenCL") and the Abseil Common libraries [abseil-cpp](https://github.com/abseil/abseil-cpp "abseil-cpp"). To use it on your host system, follow the steps:
###### Prerequisites:
1. Ensure that you have a valid OpenCL eco-system (compatible driver, headers, compiler) installed on your machine.
2. Ensure that you have a C++11 compatible compiler (eg : `g++` or `clang` ) and Cmake > 3.12 installed.
###### Init repository:
3. `git clone https://github.com/p-gerhard/reduced-kinetic-models-solver.git` -- Download the sources
4. `pip install -r requirements.txt` -- Install dependencies
5. `git submodule init` -- Initialize `abseil-cpp` and `pybind11` as git submodules in `extern` folder.
6. `git submodule update` -- Fetch all data from both projects
###### Build and install Abseil Common libraries from sources:
7. `cd ./extern/abseil-cpp/ && mkdir build && cd build`:
8. `cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE=True` -- Warning the PIC flag is important in order to link with the future shared module produced with pybind11 !
9. `make -j` -- Compile Abseil libraries
10. `sudo make install` -- System-wide installation 
###### Build the Shared C++11 lib mesh_process using pybind11:
11. `cd ./../../.. && mkdir build && cd build && cmake ..` 
12. `make` -- Create the shared librarie `mesh_process` in the `/modules/mesh/lib/` folder

### Example
To run the provided example, type:
1. `cd examples`
2. `python ./solve_3d_m1.py`
