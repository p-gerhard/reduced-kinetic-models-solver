#!/bin/bash

# Wrapper function to echo a message with a timestamp
echo_t() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] - $1"
}

echo_t Installing python3 dependencies using apt...
sudo apt install -y -q python3-numpy python3-h5py python3-meshio python3-pyopencl \
    python3-matplotlib

echo_t Initializing gitsubmodules 'abseil-cpp' and 'pybind11'...
git submodule init

echo_t Fetching last commit of 'abseil-cpp/master' and 'pybind11/master'...
git submodule update --remote --merge

echo_t Building 'abseil-cpp' target...
cd ./extern/abseil-cpp/ || exit
mkdir -p build
cd build || exit
cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE=True -DCMAKE_ABSL_PROPAGATE_CXX_STD=ON

echo_t Compiling 'abseil-cpp'...
make -j

echo_t Installing 'abseil-cpp' system-wide...
sudo make install

echo_t Building 'reduced-kinetic-models-solver' target...
cd ./../../.. || exit
mkdir -p build 
cd build || exit
cmake ..

echo_t Compiling 'reduced-kinetic-models-solver'...
make -j
