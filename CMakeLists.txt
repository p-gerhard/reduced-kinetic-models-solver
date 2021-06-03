cmake_minimum_required(VERSION 3.13)
set(CMAKE_BUILD_TYPE Release)

project(mesh_process LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../modules/mesh/lib)
find_package(absl REQUIRED)

add_subdirectory(./extern/pybind11)
add_subdirectory(./modules/mesh/lib)

target_link_libraries(mesh_process PRIVATE absl::flat_hash_map)

target_compile_options(mesh_process PRIVATE $<$<CXX_COMPILER_ID:GNU>: 
                                            -Wall 
                                            -O3>)



message("INFO: Python module created in: " ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
                     