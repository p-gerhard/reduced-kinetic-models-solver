
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "mesh_q4_process.h"
#include "mesh_h8_process.h"

void pybind11_wrapper_mesh_extract_cell_size(
	pybind11::array_t<float> np_nodes, pybind11::array_t<long> np_cells,
	pybind11::array_t<float> np_cell_size, bool is_2d)
{
	/* Unpack numpy C struct */
	pybind11::buffer_info info_nodes = np_nodes.request();
	pybind11::buffer_info info_cells = np_cells.request();
	pybind11::buffer_info info_cell_size = np_cell_size.request();

	float *nodes = static_cast<float *>(info_nodes.ptr);
	long *cells = static_cast<long *>(info_cells.ptr);
	float *cell_size = static_cast<float *>(info_cell_size.ptr);

	if (is_2d) {
		mesh_q4_extract_cell_size(nodes, cells, cell_size);
	} else {
		mesh_h8_extract_cell_size(nodes, cells, cell_size);
	}
}

void pybind11_wrapper_mesh_build_elem2elem(const long nb_cells,
										   pybind11::array_t<long> np_cells,
										   pybind11::array_t<long> np_elem2elem,
										   bool is_2d)
{
	/* Unpack numpy C struct */
	pybind11::buffer_info info_cells = np_cells.request();
	pybind11::buffer_info info_elem2elem = np_elem2elem.request();

	long *cells = static_cast<long *>(info_cells.ptr);
	long *elem2elem = static_cast<long *>(info_elem2elem.ptr);
	if (is_2d) {
		mesh_q4_build_elem2elem(nb_cells, cells, elem2elem);
	} else {
		mesh_h8_build_elem2elem(nb_cells, cells, elem2elem);
	}
}

void pybind11_wrapper_mesh_process_cells(
	const long nb_cells, pybind11::array_t<float> np_cell_size,
	pybind11::array_t<float> np_nodes, pybind11::array_t<long> np_cells,
	pybind11::array_t<float> np_cells_center, bool is_2d)
{
	/* Unpack numpy C struct */
	pybind11::buffer_info info_cell_size = np_cell_size.request();
	pybind11::buffer_info info_nodes = np_nodes.request();
	pybind11::buffer_info info_cells = np_cells.request();
	pybind11::buffer_info info_cells_center = np_cells_center.request();

	float *cell_size = static_cast<float *>(info_cell_size.ptr);
	float *nodes = static_cast<float *>(info_nodes.ptr);
	long *cells = static_cast<long *>(info_cells.ptr);
	float *cells_center = static_cast<float *>(info_cells_center.ptr);

	if (is_2d) {
		mesh_q4_process_cells(nb_cells, cell_size, nodes, cells, cells_center);
	} else {
		mesh_h8_process_cells(nb_cells, cell_size, nodes, cells, cells_center);
	}
}

void pybind11_wrapper_mesh_process_nodes(const long nb_nodes,
										 pybind11::array_t<float> np_cell_size,
										 pybind11::array_t<float> np_nodes,
										 bool is_2d)
{
	/* Unpack numpy C struct */
	pybind11::buffer_info info_cell_size = np_cell_size.request();
	pybind11::buffer_info info_nodes = np_nodes.request();

	float *cell_size = static_cast<float *>(info_cell_size.ptr);
	float *nodes = static_cast<float *>(info_nodes.ptr);
	if (is_2d) {
		mesh_q4_process_nodes(nb_nodes, cell_size, nodes);
	} else {
		mesh_h8_process_nodes(nb_nodes, cell_size, nodes);
	}
}

void pybind11_wrapper_mesh_check_elem2elem(const long nb_cells,
										   pybind11::array_t<long> np_elem2elem,
										   bool is_2d)
{
	/* Unpack numpy C struct */
	pybind11::buffer_info info_elem2elem = np_elem2elem.request();

	long *elem2elem = static_cast<long *>(info_elem2elem.ptr);
	if (is_2d) {
		mesh_q4_check_elem2elem(nb_cells, elem2elem);
	} else {
		mesh_h8_check_elem2elem(nb_cells, elem2elem);
	}
}

PYBIND11_MODULE(mesh_process, m)
{
	m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: cmake_example
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

   
    m.def("extract_cell_size", &pybind11_wrapper_mesh_extract_cell_size, 
    R"pbdoc(
        Test module.
    )pbdoc");
   
    m.def("build_elem2elem", &pybind11_wrapper_mesh_build_elem2elem, 
    R"pbdoc(
        Test module.
    )pbdoc");

    m.def("process_cells", &pybind11_wrapper_mesh_process_cells, 
    R"pbdoc(
        Test module.
    )pbdoc");

    m.def("process_nodes", &pybind11_wrapper_mesh_process_nodes, 
    R"pbdoc(
        Test module.
    )pbdoc");

	m.def("check_elem2elem", &pybind11_wrapper_mesh_check_elem2elem, 
    R"pbdoc(
        Test module.
    )pbdoc");
}
