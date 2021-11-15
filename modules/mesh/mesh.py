import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

import os
import meshio
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
import time

from ._element import *

from .lib.mesh_process import (
    extract_cell_size,
    build_elem2elem,
    process_cells,
    process_nodes,
    check_elem2elem,
)


class MeshStructured:
    def __init__(self, filename, dim, mesh_param=None):

        self.filename = filename
        assert dim in [2, 3]
        self.dim = dim

        self.elem_data = None
        is_2d = False
        
        # Load proper cell type according to the physical dimension
        if self.dim == 2:
            is_2d = True
            self.cell_name = "quad"
            self.elem_data = Q4_ELEM

        if self.dim == 3:
            self.cell_name = "hexahedron"
            self.elem_data = H8_ELEM

        self.verbose = False
        # Generate ordered 2d mesh
        if mesh_param is not None:
            assert dim == 2
            self.points, self.cells = get_ordered_2d_mesh(
                mesh_param["nx"],
                mesh_param["ny"],
                mesh_param["xmin"],
                mesh_param["xmax"],
                mesh_param["ymin"],
                mesh_param["ymax"],
            )
            self.nb_nodes, _ = np.int64(self.nodes.shape)
            self.nb_cells, _ = np.int64(self.cells.shape)
            self.dx = abs(mesh_param["xmax"] - mesh_param["xmin"]) / mesh_param["nx"]
            self.dy = abs(mesh_param["ymax"] - mesh_param["ymin"]) / mesh_param["ny"]

        # Read mesh using meshio and process it
        else:
            assert isinstance(filename, str)

            start = time.time()
            logger.info("{:<30}: {:<30}".format("Reading file", self.filename))
            self.mesh = meshio.read(filename)
            duration = time.time() - start
            logger.info("{:<30}: {:<6.8f} sec.".format("Read file in", duration))
            assert self.cell_name in self.mesh.cells_dict.keys()

            # Ref. on meshio points data
            self.cells = np.asarray(
                self.mesh.cells_dict[self.cell_name], dtype=np.int64
            )
            self.nodes = np.asarray(self.mesh.points, dtype=np.float32)

            if is_2d:
                self.nodes = np.delete(self.nodes, 2, 1)

            self.nb_nodes, phy_dim_point = self.nodes.shape
            assert phy_dim_point == self.dim

            self.nb_cells, nb_loc_node = self.cells.shape
            assert nb_loc_node == self.elem_data["NODE_PER_ELEM"]
            logger.info("{:<30}: {:<d}".format("Number of cells", self.nb_cells))

            self.cells = np.ravel(self.cells, order="C")
            self.points = np.ravel(self.nodes, order="C")

            self.cells_center = np.empty(self.dim * self.nb_cells, dtype=np.float32)
            size_elem2elem = self.elem_data["FACE_PER_ELEM"] * self.nb_cells
            self.elem2elem = np.full(size_elem2elem, -1, dtype=np.int64)
            self.cell_size = np.zeros(self.dim, dtype=np.float32)

            start = time.time()
            extract_cell_size(self.nodes, self.cells, self.cell_size, is_2d)
            duration = time.time() - start

            logger.info("{:<30}: {:>6.8f} sec.".format("Sizes extracted in", duration))

            self.dx = np.float32(self.cell_size[0])
            logger.info("{:<30}: {:>6.8f}".format("Size dx", self.cell_size[0]))
            
            self.dy = np.float32(self.cell_size[1])
            logger.info("{:<30}: {:>6.8f}".format("Size dy", self.cell_size[1]))


            if is_2d :
                 self.dz = np.float32(0.)
            else :
                self.dz = np.float32(self.cell_size[2])
                logger.info("{:<30}: {:>6.8f}".format("Size dz", self.cell_size[2]))

            start = time.time()
            process_cells(
                self.nb_cells,
                self.cell_size,
                self.nodes,
                self.cells,
                self.cells_center,
                is_2d,
            )

            duration = time.time() - start
            logger.info("{:<30}: {:>6.8f} sec.".format("Cells proceeded in", duration))
            start = time.time()

            process_nodes(self.nb_nodes, self.cell_size, self.nodes, is_2d)
            duration = time.time() - start
            logger.info("{:<30}: {:>6.8f} sec.".format("Nodes proceeded in", duration))
            start = time.time()

            build_elem2elem(self.nb_cells, self.cells, self.elem2elem, is_2d)
            duration = time.time() - start
            logger.info(
                "{:<30}: {:>6.8f} sec.".format("Connectivity built in", duration)
            )

            start = time.time()
            check_elem2elem(self.nb_cells, self.elem2elem, is_2d)
            duration = time.time() - start
            logger.info(
                "{:<30}: {:>6.8f} sec.".format("Connectivity checked in", duration)
            )

        self.points = np.reshape(self.nodes, (self.nb_nodes, self.dim))
        self.cells = np.reshape(
            self.cells, (self.nb_cells, self.elem_data["NODE_PER_ELEM"])
        )

        # self.mesh = meshio.Mesh(self.points, {self.cell_name: self.cells})
        # self.mesh.write("mesh_output.msh")