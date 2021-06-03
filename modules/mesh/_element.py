# GMSH Quadrangle:
#
#       v
#       ^
#       |
# 3-----------2
# |     |     |
# |     |     |
# |     +---- | --> u
# |           |
# |           |
# 0-----------1

Q4_ELEM = {
    "ELEM_CODE": 3,
    "PHY_DIM": 2,
    "FACE_PER_ELEM": 4,
    "NODE_PER_ELEM": 4,
    "NODE_PER_FACE": 2,
    "FACE_TO_LOC_NODE": [
        [1, 1],  # Right
        [0, 3],  # Left
        [3, 2],  # North
        [0, 1],  # South
    ],
    "GMSH_REF_NODE_COORD": [
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
    ],
}


# GMSH Order Hexahedron:
#
#        v
# 3----------2
# |\     ^   |\
# | \    |   | \
# |  \   |   |  \
# |   7------+---6
# |   |  +-- |-- | -> u
# 0---+---\--1   |
#  \  |    \  \  |       8
#   \ |     \  \ |
#    \|      w  \|
#     4----------5






H8_ELEM = {
    "ELEM_CODE": 5,
    "PHY_DIM": 3,
    "FACE_PER_ELEM": 6,
    "NODE_PER_ELEM": 8,
    "NODE_PER_FACE": 4,
    "FACE_TO_LOC_NODE": [
        [5, 1, 6, 2],  # Right
        [4, 7, 0, 3],  # Left
        [6, 2, 7, 3],  # Front
        [4, 0, 5, 1],  # Back
        [4, 5, 7, 6],  # North
        [0, 3, 1, 2],  # South
    ],
    "GMSH_REF_NODE_COORD": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ],
}
