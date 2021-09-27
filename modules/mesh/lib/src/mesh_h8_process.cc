#undef NDEBUG
#include <cmath>
#include <assert.h>
#include <iostream>

#include <absl/container/flat_hash_map.h>
#include "mesh_h8_process.h"

/* GMSH hexahedron loc_nodes_id:
 *
 *         y
 *  3----------2
 *  |\     ^   |\
 *  | \    |   | \
 *  |  \   |   |  \
 *  |   7------+---6
 *  |   |  +-- |-- | -> x
 *  0---+---\--1   |
 *   \  |    \  \  |       8
 *    \ |     \  \ |
 *     \|      z  \|
 *      4----------5
 */

#define H8_ELEM_CODE 5
#define H8_PHY_DIM 3
#define H8_FACE_PER_ELEM 6
#define H8_NODE_PER_ELEM 8
#define H8_NODE_PER_FACE 4

const int h8_face_to_loc_node[H8_FACE_PER_ELEM][H8_NODE_PER_FACE] = {
	{ 1, 2, 6, 5 }, // Right
	{ 0, 4, 7, 3 }, // Left
	{ 2, 3, 7, 6 }, // Front
	{ 1, 5, 4, 0 }, // Back
	{ 4, 5, 6, 7 }, // North
	{ 0, 3, 2, 1 }, // South
};

const float h8_face_normals[H8_FACE_PER_ELEM][3] = {
	{ 1, 0, 0 }, // Right
	{ -1, 0, 0 }, // Left
	{ 0, 1, 0 }, // Front
	{ 0, -1, 0 }, // Back
	{ 0, 0, 1 }, // North
	{ 0, 0, -1 } // South
};

const int h8_face_edge_len[H8_FACE_PER_ELEM][2] = {
	{ 2, 1 }, // Right : dz, dy
	{ 2, 1 }, // Left  : dz, dy
	{ 2, 0 }, // Front : dz, dx
	{ 2, 0 }, // Back  : dz, dx
	{ 0, 1 }, // North : dx, dy
	{ 0, 1 }, // South : dx, dy
};
const int h8_lex_order[H8_NODE_PER_ELEM] = { 0, 4, 3, 7, 1, 5, 2, 6 };

#define MESH_TOL 1e-6

typedef std::tuple<long, long, long, long> key4_t;
using namespace absl;

static inline void cross_product(const float x[3], const float y[3],
								 float res[3])
{
	res[0] = x[1] * y[2] - x[2] * y[1];
	res[1] = x[2] * y[0] - x[0] * y[2];
	res[2] = x[0] * y[1] - x[1] * y[0];
}

static inline void lexsort_h8(float nodes[24], long nodes_id[8])
{
#define SWAP_NODE(i, j)                                                        \
	const float tmp_x = nodes[3 * i + 0];                                      \
	const float tmp_y = nodes[3 * i + 1];                                      \
	const float tmp_z = nodes[3 * i + 2];                                      \
	nodes[3 * i + 0] = nodes[3 * j + 0];                                       \
	nodes[3 * i + 1] = nodes[3 * j + 1];                                       \
	nodes[3 * i + 2] = nodes[3 * j + 2];                                       \
	nodes[3 * j + 0] = tmp_x;                                                  \
	nodes[3 * j + 1] = tmp_y;                                                  \
	nodes[3 * j + 2] = tmp_z;                                                  \
	const long tmp_int = nodes_id[i];                                          \
	nodes_id[i] = nodes_id[j];                                                 \
	nodes_id[j] = tmp_int;

#define REORDER(i, j)                                                          \
	if (nodes[3 * j + 0] < nodes[3 * i + 0]) {                                 \
		SWAP_NODE(i, j)                                                        \
	}                                                                          \
	if (fabs(nodes[3 * j + 0] - nodes[3 * i + 0]) < MESH_TOL &&                \
		nodes[3 * j + 1] < nodes[3 * i + 1]) {                                 \
		SWAP_NODE(i, j)                                                        \
	}                                                                          \
	if (fabs(nodes[3 * j + 0] - nodes[3 * i + 0]) < MESH_TOL &&                \
		fabs(nodes[3 * j + 1] - nodes[3 * i + 1]) < MESH_TOL &&                \
		nodes[3 * j + 2] < nodes[3 * i + 2]) {                                 \
		SWAP_NODE(i, j)                                                        \
	}

	/* BUBBLE SORT */
	REORDER(0, 1);
	REORDER(1, 2);
	REORDER(2, 3);
	REORDER(3, 4);
	REORDER(4, 5);
	REORDER(5, 6);
	REORDER(6, 7);
	REORDER(0, 1);
	REORDER(1, 2);
	REORDER(2, 3);
	REORDER(3, 4);
	REORDER(4, 5);
	REORDER(5, 6);
	REORDER(0, 1);
	REORDER(1, 2);
	REORDER(2, 3);
	REORDER(3, 4);
	REORDER(4, 5);
	REORDER(0, 1);
	REORDER(1, 2);
	REORDER(2, 3);
	REORDER(3, 4);
	REORDER(0, 1);
	REORDER(1, 2);
	REORDER(2, 3);
	REORDER(0, 1);
	REORDER(1, 2);
	REORDER(0, 1);
#undef REORDER
#undef SWAP
}

void mesh_h8_extract_cell_size(const float *nodes, long *cells,
							   float cell_size[3])
{
	long loc_cell_nodes_id[H8_NODE_PER_ELEM];
	float loc_nodes[H8_NODE_PER_ELEM * H8_PHY_DIM];
	/* Local copy of FIRST cell's node_id */
	for (int i = 0; i < H8_NODE_PER_ELEM; i++) {
		loc_cell_nodes_id[i] = cells[i];
	}

	/* Local copy of FIRST node's coordinates */
	for (int i = 0; i < H8_NODE_PER_ELEM; i++) {
		for (int k = 0; k < H8_PHY_DIM; k++) {
			long id_node = loc_cell_nodes_id[i];
			loc_nodes[H8_PHY_DIM * i + k] = nodes[H8_PHY_DIM * id_node + k];
		}
	}

	lexsort_h8(loc_nodes, loc_cell_nodes_id);

	/* Re-order */
	for (int i = 0; i < H8_NODE_PER_ELEM; i++) {
		cells[h8_lex_order[i]] = loc_cell_nodes_id[i];
	}

	for (int i = 0; i < H8_NODE_PER_ELEM; i++) {
		long id_node = cells[i];
		for (int k = 0; k < H8_PHY_DIM; k++) {
			loc_nodes[H8_PHY_DIM * i + k] = nodes[H8_PHY_DIM * id_node + k];
		}
	}

	const int dx_edge[4][2] = { { 0, 1 }, { 3, 2 }, { 4, 5 }, { 7, 6 } };
	const int dy_edge[4][2] = { { 0, 3 }, { 1, 2 }, { 4, 7 }, { 5, 6 } };
	const int dz_edge[4][2] = { { 0, 4 }, { 1, 5 }, { 3, 7 }, { 2, 6 } };

	float dx[4];
	float dy[4];
	float dz[4];

	int n0;
	int n1;

	/* Computing each H8 edge length */
	for (int k = 0; k < 4; k++) {
		/* dx */
		n0 = dx_edge[k][0];
		n1 = dx_edge[k][1];
		float ex[3] = { loc_nodes[3 * n1 + 0] - loc_nodes[3 * n0 + 0],
						loc_nodes[3 * n1 + 1] - loc_nodes[3 * n0 + 1],
						loc_nodes[3 * n1 + 2] - loc_nodes[3 * n0 + 2] };

		dx[k] = sqrt(ex[0] * ex[0] + ex[1] * ex[1] + ex[2] * ex[2]);

		/* dy */
		n0 = dy_edge[k][0];
		n1 = dy_edge[k][1];
		float ey[3] = { loc_nodes[3 * n1 + 0] - loc_nodes[3 * n0 + 0],
						loc_nodes[3 * n1 + 1] - loc_nodes[3 * n0 + 1],
						loc_nodes[3 * n1 + 2] - loc_nodes[3 * n0 + 2] };

		dy[k] = sqrt(ey[0] * ey[0] + ey[1] * ey[1] + ey[2] * ey[2]);

		/* dz */
		n0 = dz_edge[k][0];
		n1 = dz_edge[k][1];
		float ez[3] = { loc_nodes[3 * n1 + 0] - loc_nodes[3 * n0 + 0],
						loc_nodes[3 * n1 + 1] - loc_nodes[3 * n0 + 1],
						loc_nodes[3 * n1 + 2] - loc_nodes[3 * n0 + 2] };

		dz[k] = sqrt(ez[0] * ez[0] + ez[1] * ez[1] + ez[2] * ez[2]);
	}
	/* Checking edge's length */
	for (int k = 1; k < 4; k++) {
		assert(fabs(dx[k - 1] - dx[k]) < MESH_TOL);
		assert(fabs(dy[k - 1] - dy[k]) < MESH_TOL);
		assert(fabs(dz[k - 1] - dz[k]) < MESH_TOL);
	}

	cell_size[0] = dx[0];
	cell_size[1] = dy[0];
	cell_size[2] = dz[0];
}

static bool check_normals_h8(float nodes[24], const float cell_dim[3])
{
	bool check = true;

	for (int id_face = 0; id_face < H8_FACE_PER_ELEM; id_face++) {
		/* Compute face edge between node (n1) and node (n0) */
		long n0 = h8_face_to_loc_node[id_face][0];
		long n1 = h8_face_to_loc_node[id_face][1];

		const float e0[3] = { nodes[3 * n1 + 0] - nodes[3 * n0 + 0],
							  nodes[3 * n1 + 1] - nodes[3 * n0 + 1],
							  nodes[3 * n1 + 2] - nodes[3 * n0 + 2] };

		/* Compute face edge between node (n2) and node (n1) */
		long n2 = h8_face_to_loc_node[id_face][2];

		const float e1[3] = { nodes[3 * n2 + 0] - nodes[3 * n1 + 0],
							  nodes[3 * n2 + 1] - nodes[3 * n1 + 1],
							  nodes[3 * n2 + 2] - nodes[3 * n1 + 2] };

		/* Check orthogonality */
		const float e0_dot_e1 = e0[0] * e1[0] + e0[1] * e1[1] + e0[2] * e1[2];
		check = check && (fabs(e0_dot_e1) < MESH_TOL);

		/* Check both edge length */
		const float e0_ln = sqrt(e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2]);
		const float e1_ln = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);

		int edge_len_id = h8_face_edge_len[id_face][0];
		check = (check && (fabs(e0_ln - cell_dim[edge_len_id]) < MESH_TOL));

		edge_len_id = h8_face_edge_len[id_face][1];
		check = check && (fabs(e1_ln - cell_dim[edge_len_id]) < MESH_TOL);

		/* Check normal direction */
		float n[3];
		float a = (e0_ln * e1_ln);
		cross_product(e0, e1, n);

		const float n_ref[3] = { h8_face_normals[id_face][0],
								 h8_face_normals[id_face][1],
								 h8_face_normals[id_face][2] };

		check = check && (fabs(n_ref[0] - n[0] / a) < MESH_TOL);
		check = check && (fabs(n_ref[1] - n[1] / a) < MESH_TOL);
		check = check && (fabs(n_ref[2] - n[2] / a) < MESH_TOL);
	}

	return check;
}

void mesh_h8_process_cells(const long nb_cells, const float cell_size[3],
						   const float *nodes, long *cells, float *cells_center)
{
	long loc_cell_nodes_id[H8_NODE_PER_ELEM];
	float loc_nodes[H8_NODE_PER_ELEM * H8_PHY_DIM];

	for (long id_cell = 0; id_cell < nb_cells; id_cell++) {
		/* Local copy of cell's node_id */
		for (int i = 0; i < H8_NODE_PER_ELEM; i++) {
			loc_cell_nodes_id[i] = cells[H8_NODE_PER_ELEM * id_cell + i];
		}

		/* Local copy of node's coordinates */
		for (int i = 0; i < H8_NODE_PER_ELEM; i++) {
			for (int k = 0; k < H8_PHY_DIM; k++) {
				long id_node = loc_cell_nodes_id[i];
				loc_nodes[H8_PHY_DIM * i + k] = nodes[H8_PHY_DIM * id_node + k];
			}
		}

		lexsort_h8(loc_nodes, loc_cell_nodes_id);

		/* 
         * Re-order cell's node_id and  using the mapping beetween lexicographic 
         * order and desired (Gmsh) order.
         */
		for (int i = 0; i < H8_NODE_PER_ELEM; i++) {
			cells[H8_NODE_PER_ELEM * id_cell + h8_lex_order[i]] =
				loc_cell_nodes_id[i];
		}

		/* 
         * Re-order loc_nodes accordingly to cells's node_id new 
         * order and make the sum for cells centers computation.
         */
		float sum_x = 0.f;
		float sum_y = 0.f;
		float sum_z = 0.f;

		for (int i = 0; i < H8_NODE_PER_ELEM; i++) {
			long id_node = cells[H8_NODE_PER_ELEM * id_cell + i];
			for (int k = 0; k < H8_PHY_DIM; k++) {
				loc_nodes[H8_PHY_DIM * i + k] = nodes[H8_PHY_DIM * id_node + k];
			}
			sum_x += loc_nodes[H8_PHY_DIM * i + 0];
			sum_y += loc_nodes[H8_PHY_DIM * i + 1];
			sum_z += loc_nodes[H8_PHY_DIM * i + 2];
		}

		/* Check norms */
		bool check = check_normals_h8(loc_nodes, cell_size);
		assert(check);

		/* Compute cell centers using true hexahedron formula */
		const float center_1[3] = { sum_x / 8.f, sum_y / 8.f, sum_z / 8.f };

		/* Compute cell centers using displacement from first cell's node */
		const float center_2[3] = { loc_nodes[0] + 0.5f * cell_size[0],
									loc_nodes[1] + 0.5f * cell_size[1],
									loc_nodes[2] + 0.5f * cell_size[2] };

		/* Verify if both method are giving the same results */
		check = check && (fabs(center_1[0] - center_2[0]) < MESH_TOL);
		check = check && (fabs(center_1[1] - center_2[1]) < MESH_TOL);
		check = check && (fabs(center_1[2] - center_2[2]) < MESH_TOL);
		assert(check);

		cells_center[H8_PHY_DIM * id_cell + 0] = center_1[0];
		cells_center[H8_PHY_DIM * id_cell + 1] = center_1[1];
		cells_center[H8_PHY_DIM * id_cell + 2] = center_1[2];
	}
}

void mesh_h8_process_nodes(const long nb_nodes, const float cell_size[3],
						   float *nodes)
{
	/* TODO verif deplacement */
	const float half_dx = 0.5f * cell_size[0];
	const float half_dy = 0.5f * cell_size[1];
	const float half_dz = 0.5f * cell_size[2];

	float old_coord[H8_PHY_DIM];
	for (long id_node = 0; id_node < nb_nodes; id_node++) {
		old_coord[0] = nodes[H8_PHY_DIM * id_node + 0];
		old_coord[1] = nodes[H8_PHY_DIM * id_node + 1];

		const float x =
			rint(nodes[H8_PHY_DIM * id_node + 0] / half_dx) * half_dx;
		const float y =
			rint(nodes[H8_PHY_DIM * id_node + 1] / half_dy) * half_dy;

		nodes[H8_PHY_DIM * id_node + 0] = x;
		nodes[H8_PHY_DIM * id_node + 1] = y;

		assert(fabs(old_coord[0] - nodes[H8_PHY_DIM * id_node + 0]) < MESH_TOL);
		assert(fabs(old_coord[1] - nodes[H8_PHY_DIM * id_node + 1]) < MESH_TOL);

		old_coord[2] = nodes[H8_PHY_DIM * id_node + 2];
		const float z =
			rint(nodes[H8_PHY_DIM * id_node + 2] / half_dz) * half_dz;

		nodes[H8_PHY_DIM * id_node + 2] = z;
		assert(fabs(old_coord[2] - nodes[H8_PHY_DIM * id_node + 2]) < MESH_TOL);
	}
}

void mesh_h8_build_elem2elem(const long nb_cells, const long *cells,
							 long *elem2elem)
{
	/* Hashmap container : key is a tuple of 4 long and the value is 1 long */
	flat_hash_map<key4_t, long> face_map;
	key4_t key;
	key4_t key_2;
	face_map.reserve(nb_cells);
	/* Build the map using direct face orientation (node id) */
	for (long id_elem = 0; id_elem < nb_cells; id_elem++) {
		long offset = id_elem * H8_NODE_PER_ELEM;
		for (int id_face = 0; id_face < H8_FACE_PER_ELEM; id_face++) {
			key = std::make_tuple(
				cells[offset + h8_face_to_loc_node[id_face][0]],
				cells[offset + h8_face_to_loc_node[id_face][1]],
				cells[offset + h8_face_to_loc_node[id_face][2]],
				cells[offset + h8_face_to_loc_node[id_face][3]]);

			auto res = face_map.try_emplace(key, id_elem);

			if (!(res.second)) {
				std::cout << "WARNING: Duplicate cell" << std::endl;
			}
		}
	}

	/* Search corresponding faces using reverse face orientation (node id) */
	for (long id_elem = 0; id_elem < nb_cells; id_elem++) {
		long offset = id_elem * H8_NODE_PER_ELEM;
		for (int id_face = 0; id_face < H8_FACE_PER_ELEM; id_face++) {
			/* The permutation of H8 face swap the 2nd and the last element */

			key = std::make_tuple(
				cells[offset + h8_face_to_loc_node[id_face][0]],
				cells[offset + h8_face_to_loc_node[id_face][3]],
				cells[offset + h8_face_to_loc_node[id_face][2]],
				cells[offset + h8_face_to_loc_node[id_face][1]]);

			auto res = face_map.find(key);

			if (res != face_map.end()) {
				elem2elem[id_elem * H8_FACE_PER_ELEM + id_face] = res->second;
			}
		}
	}
	face_map.clear();
}

void mesh_h8_check_elem2elem(const long nb_cells, const long *elem2elem)
{
	/* Check that the neighbour of the neighbour is id_elem */
	for (long id_elem = 0; id_elem < nb_cells; id_elem++) {
		for (int id_couple = 0; id_couple < H8_FACE_PER_ELEM / 2; id_couple++) {
			long e0 = elem2elem[id_elem * H8_FACE_PER_ELEM + 2 * id_couple + 0];

			/* Avoid boundary Face */
			if (e0 != -1) {
				long e1 = elem2elem[e0 * H8_FACE_PER_ELEM + 2 * id_couple + 1];
				assert(id_elem == e1);
			}
		}
	}
}