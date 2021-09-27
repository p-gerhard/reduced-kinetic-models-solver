#undef NDEBUG
#include <cmath>
#include <assert.h>
#include <iostream>

#include <absl/container/flat_hash_map.h>
#include "mesh_q4_process.h"

/* GMSH quadrangle (Q4) loc_nodes_id:
 *
 *       v
 *       ^
 *       |
 * 3-----------2
 * |     |     |
 * |     |     |
 * |     +---- | --> u
 * |           |
 * |           |
 * 0-----------1
 */
#define Q4_ELEM_CODE 3
#define Q4_PHY_DIM 2
#define Q4_FACE_PER_ELEM 4
#define Q4_NODE_PER_ELEM 4
#define Q4_NODE_PER_FACE 2

const int q4_face_to_loc_node[Q4_FACE_PER_ELEM][Q4_NODE_PER_FACE] = {
	{ 1, 2 }, // Right
	{ 3, 0 }, // Left
	{ 2, 3 }, // North
	{ 0, 1 }, // South
};

const float q4_face_normals[Q4_FACE_PER_ELEM][3] = {
	{ 1, 0, 0 }, // Right
	{ -1, 0, 0 }, // Left
	{ 0, 1, 0 }, // North
	{ 0, -1, 0 }, // South
};

const int q4_face_edge_len[Q4_FACE_PER_ELEM] = { 1, 0, 1, 0 }; // dy, dx, dy, dx

const unsigned int q4_lex_order[Q4_NODE_PER_ELEM] = { 0, 3, 1, 2 };

#define MESH_TOL 1e-6

typedef std::tuple<long, long, long, long> key4_t;
typedef std::tuple<long, long> key2_t;
using namespace absl;

static inline void cross_product(const float x[3], const float y[3],
								 float res[3])
{
	res[0] = x[1] * y[2] - x[2] * y[1];
	res[1] = x[2] * y[0] - x[0] * y[2];
	res[2] = x[0] * y[1] - x[1] * y[0];
}

static inline void lexsort_q4(float nodes[12], long nodes_id[4])
{
#define SWAP_NODE(i, j)                                                        \
	const float tmp_x = nodes[2 * i + 0];                                      \
	const float tmp_y = nodes[2 * i + 1];                                      \
	nodes[2 * i + 0] = nodes[2 * j + 0];                                       \
	nodes[2 * i + 1] = nodes[2 * j + 1];                                       \
	nodes[2 * j + 0] = tmp_x;                                                  \
	nodes[2 * j + 1] = tmp_y;                                                  \
	const long tmp_int = nodes_id[i];                                          \
	nodes_id[i] = nodes_id[j];                                                 \
	nodes_id[j] = tmp_int;

#define REORDER(i, j)                                                          \
	if (nodes[2 * j + 0] < nodes[2 * i + 0]) {                                 \
		SWAP_NODE(i, j)                                                        \
	}                                                                          \
	if (fabs(nodes[2 * j + 0] - nodes[2 * i + 0]) < MESH_TOL &&                \
		nodes[2 * j + 1] < nodes[2 * i + 1]) {                                 \
		SWAP_NODE(i, j)                                                        \
	}

	/* BUBBLE SORT */
	REORDER(0, 1);
	REORDER(2, 3);
	REORDER(0, 2);
	REORDER(1, 3);
	REORDER(1, 2);
#undef REORDER
#undef SWAP
}

void mesh_q4_extract_cell_size(const float *nodes, long *cells,
							   float cell_size[2])
{
	long loc_cell_nodes_id[Q4_NODE_PER_ELEM];
	float loc_nodes[Q4_NODE_PER_ELEM * Q4_PHY_DIM];

	/* Local copy of FIRST cell's node_id */
	for (int i = 0; i < Q4_NODE_PER_ELEM; i++) {
		loc_cell_nodes_id[i] = cells[i];
	}

	/* Local copy of FIRST node's coordinates */
	for (int i = 0; i < Q4_NODE_PER_ELEM; i++) {
		for (int k = 0; k < Q4_PHY_DIM; k++) {
			long id_node = loc_cell_nodes_id[i];
			loc_nodes[Q4_PHY_DIM * i + k] = nodes[Q4_PHY_DIM * id_node + k];
		}
	}

	lexsort_q4(loc_nodes, loc_cell_nodes_id);

	/* Re-order */
	for (int i = 0; i < Q4_NODE_PER_ELEM; i++) {
		cells[q4_lex_order[i]] = loc_cell_nodes_id[i];
	}

	for (int i = 0; i < Q4_NODE_PER_ELEM; i++) {
		long id_node = cells[i];
		for (int k = 0; k < Q4_PHY_DIM; k++) {
			loc_nodes[Q4_PHY_DIM * i + k] = nodes[Q4_PHY_DIM * id_node + k];
		}
	}

	const int dx_edge[2][2] = { { 0, 1 }, { 3, 2 } };
	const int dy_edge[2][2] = { { 0, 3 }, { 1, 2 } };

	float dx[2];
	float dy[2];

	int n0;
	int n1;

	/* Computing each Q4 edge length */
	for (int k = 0; k < 2; k++) {
		/* dx */
		n0 = dx_edge[k][0];
		n1 = dx_edge[k][1];
		float ex[2] = {
			loc_nodes[Q4_PHY_DIM * n1 + 0] - loc_nodes[Q4_PHY_DIM * n0 + 0],
			loc_nodes[Q4_PHY_DIM * n1 + 1] - loc_nodes[Q4_PHY_DIM * n0 + 1]
		};

		dx[k] = sqrt(ex[0] * ex[0] + ex[1] * ex[1]);

		/* dy */
		n0 = dy_edge[k][0];
		n1 = dy_edge[k][1];
		float ey[2] = {
			loc_nodes[Q4_PHY_DIM * n1 + 0] - loc_nodes[Q4_PHY_DIM * n0 + 0],
			loc_nodes[Q4_PHY_DIM * n1 + 1] - loc_nodes[Q4_PHY_DIM * n0 + 1]
		};

		dy[k] = sqrt(ey[0] * ey[0] + ey[1] * ey[1]);
	}
	/* Checking edge's length */
	for (int k = 1; k < 2; k++) {
		assert(fabs(dx[k - 1] - dx[k]) < MESH_TOL);
		assert(fabs(dy[k - 1] - dy[k]) < MESH_TOL);
	}

	cell_size[0] = dx[0];
	cell_size[1] = dy[0];
}

static bool check_normals_q4(float nodes[12], const float cell_dim[2])
{
	bool check = true;

	for (int id_face = 0; id_face < Q4_FACE_PER_ELEM; id_face++) {
		/* Compute face edge between node (n1) and node (n0) */
		long n0 = q4_face_to_loc_node[id_face][0];
		long n1 = q4_face_to_loc_node[id_face][1];

		const float e[3] = {
			nodes[Q4_PHY_DIM * n1 + 0] - nodes[Q4_PHY_DIM * n0 + 0],
			nodes[Q4_PHY_DIM * n1 + 1] - nodes[Q4_PHY_DIM * n0 + 1], 0
		};

		/* Check edge length */
		const float e_ln = sqrt(e[0] * e[0] + e[1] * e[1]);

		int edge_len_id = q4_face_edge_len[id_face];
		check = (check && (fabs(e_ln - cell_dim[edge_len_id]) < MESH_TOL));

		/*@TODO SOLVE Check normal direction */
		// float n[3];
		// float a = (e_ln * e_ln);
		// cross_product(e0, e1, n);

		// const float n_ref[3] = { face_normals[id_face][0],
		// 						 face_normals[id_face][1], 0 };

		// check = check && (fabs(n_ref[0] - n[0] / a) < MESH_TOL);
		// check = check && (fabs(n_ref[1] - n[1] / a) < MESH_TOL);
	}

	return check;
}

void mesh_q4_process_cells(const long nb_cells, const float cell_size[3],
						   const float *nodes, long *cells, float *cells_center)
{
	long loc_cell_nodes_id[Q4_NODE_PER_ELEM];
	float loc_nodes[Q4_NODE_PER_ELEM * Q4_PHY_DIM];

	for (long id_cell = 0; id_cell < nb_cells; id_cell++) {
		/* Local copy of cell's node_id */
		for (int i = 0; i < Q4_NODE_PER_ELEM; i++) {
			loc_cell_nodes_id[i] = cells[Q4_NODE_PER_ELEM * id_cell + i];
		}

		/* Local copy of node's coordinates */
		for (int i = 0; i < Q4_NODE_PER_ELEM; i++) {
			for (int k = 0; k < Q4_PHY_DIM; k++) {
				long id_node = loc_cell_nodes_id[i];
				loc_nodes[Q4_PHY_DIM * i + k] = nodes[Q4_PHY_DIM * id_node + k];
			}
		}

		lexsort_q4(loc_nodes, loc_cell_nodes_id);

		/* 
         * Re-order cell's node_id and  using the mapping beetween lexicographic 
         * order and desired (Gmsh) order.
         */
		for (int i = 0; i < Q4_NODE_PER_ELEM; i++) {
			cells[Q4_NODE_PER_ELEM * id_cell + q4_lex_order[i]] =
				loc_cell_nodes_id[i];
		}

		/* 
         * Re-order loc_nodes accordingly to cells's node_id new 
         * order and make the sum for cells centers computation.
         */
		float sum_x = 0.f;
		float sum_y = 0.f;

		for (int i = 0; i < Q4_NODE_PER_ELEM; i++) {
			long id_node = cells[Q4_NODE_PER_ELEM * id_cell + i];
			for (int k = 0; k < Q4_PHY_DIM; k++) {
				loc_nodes[Q4_PHY_DIM * i + k] = nodes[Q4_PHY_DIM * id_node + k];
			}
			sum_x += loc_nodes[Q4_PHY_DIM * i + 0];
			sum_y += loc_nodes[Q4_PHY_DIM * i + 1];
		}

		/* Check norms */
		bool check = check_normals_q4(loc_nodes, cell_size);
		assert(check);

		/* Compute cell centers using true rectangle formula */
		const float center_1[2] = { sum_x / 4.f, sum_y / 4.f };

		/* Compute cell centers using displacement from first cell's node */
		const float center_2[2] = { loc_nodes[0] + 0.5f * cell_size[0],
									loc_nodes[1] + 0.5f * cell_size[1] };

		/* Verify if both method are giving the same results */
		check = check && (fabs(center_1[0] - center_2[0]) < MESH_TOL);
		check = check && (fabs(center_1[1] - center_2[1]) < MESH_TOL);
		assert(check);

		cells_center[Q4_PHY_DIM * id_cell + 0] = center_1[0];
		cells_center[Q4_PHY_DIM * id_cell + 1] = center_1[1];
	}
}

void mesh_q4_process_nodes(const long nb_nodes, const float cell_size[3],
						   float *nodes)
{
	/* TODO verif deplacement */
	const float half_dx = 0.5f * cell_size[0];
	const float half_dy = 0.5f * cell_size[1];

	float old_coord[Q4_PHY_DIM];
	for (long id_node = 0; id_node < nb_nodes; id_node++) {
		old_coord[0] = nodes[Q4_PHY_DIM * id_node + 0];
		old_coord[1] = nodes[Q4_PHY_DIM * id_node + 1];

		const float x =
			rint(nodes[Q4_PHY_DIM * id_node + 0] / half_dx) * half_dx;
		const float y =
			rint(nodes[Q4_PHY_DIM * id_node + 1] / half_dy) * half_dy;

		nodes[Q4_PHY_DIM * id_node + 0] = x;
		nodes[Q4_PHY_DIM * id_node + 1] = y;

		assert(fabs(old_coord[0] - nodes[Q4_PHY_DIM * id_node + 0]) < MESH_TOL);
		assert(fabs(old_coord[1] - nodes[Q4_PHY_DIM * id_node + 1]) < MESH_TOL);
	}
}

void mesh_q4_build_elem2elem(const long nb_cells, const long *cells,
							 long *elem2elem)
{
	/* Hashmap container : key is a tuple of 2 long and the value is 1 long */
	flat_hash_map<key2_t, long> face_map;
	key2_t key;
	key2_t key_2;
	face_map.reserve(nb_cells);
	/* Build the map using direct face orientation (node id) */
	for (long id_elem = 0; id_elem < nb_cells; id_elem++) {
		long offset = id_elem * Q4_NODE_PER_ELEM;
		for (int id_face = 0; id_face < Q4_FACE_PER_ELEM; id_face++) {
			key = std::make_tuple(
				cells[offset + q4_face_to_loc_node[id_face][0]],
				cells[offset + q4_face_to_loc_node[id_face][1]]);

			auto res = face_map.try_emplace(key, id_elem);

			if (!(res.second)) {
				std::cout << "WARNING: Duplicate cell" << std::endl;
			}
		}
	}

	/* Search corresponding faces using reverse face orientation (node id) */
	for (long id_elem = 0; id_elem < nb_cells; id_elem++) {
		long offset = id_elem * Q4_NODE_PER_ELEM;
		for (int id_face = 0; id_face < Q4_FACE_PER_ELEM; id_face++) {
			/* The permutation of H8 face swap the 2nd and the last element */

			key = std::make_tuple(
				cells[offset + q4_face_to_loc_node[id_face][1]],
				cells[offset + q4_face_to_loc_node[id_face][0]]);

			auto res = face_map.find(key);

			if (res != face_map.end()) {
				elem2elem[id_elem * Q4_FACE_PER_ELEM + id_face] = res->second;
			}
		}
	}
	face_map.clear();
}

void mesh_q4_check_elem2elem(const long nb_cells, const long *elem2elem)
{
	/* Check that the neighbour of the neighbour is id_elem */
	for (long id_elem = 0; id_elem < nb_cells; id_elem++) {
		for (int id_couple = 0; id_couple < Q4_FACE_PER_ELEM / 2; id_couple++) {
			long e0 = elem2elem[id_elem * Q4_FACE_PER_ELEM + 2 * id_couple + 0];

			/* Avoid boundary Face */
			if (e0 != -1) {
				long e1 = elem2elem[e0 * Q4_FACE_PER_ELEM + 2 * id_couple + 1];
				assert(id_elem == e1);
			}
		}
	}
}