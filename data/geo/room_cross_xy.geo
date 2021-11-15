Mesh.Algorithm = 8;
n_long = 10;

SetFactory("OpenCASCADE");

Point(1) = {0.5, -0.5, 0, 1.0};
Point(2) = {9.5, -0.5, 0, 1.0};
Point(3) = {9.5, 0.5, 0, 1.0};
Point(4) = {0.5, 0.5, 0, 1.0};
Point(5) = {0.5, 9.5, 0, 1.0};
Point(6) = {-0.5, 9.5, 0, 1.0};
Point(7) = {-0.5, 0.5, 0, 1.0};
Point(8) = {-9.5, 0.5, 0, 1.0};
Point(9) = {-9.5, -0.5, 0, 1.0};
Point(10) = {-0.5, -0.5, 0, 1.0};
Point(11) = {-0.5, -9.5, 0, 1.0};
Point(12) = {0.5,  -9.5, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

Transfinite Line {1,3,4,6,7,9,10,12}  = n_long Using Progression 1;
Transfinite Line {2,5,8,11} = 1 Using Progression 1;

Curve Loop(1) = {6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5};
Plane Surface(1) = {1};

Transfinite Surface {2} = {8, 9, 10, 7};
Transfinite Surface {3} = {7, 4, 5, 6};
Transfinite Surface {4} = {4, 1, 2, 3};
Transfinite Surface {5} = {10, 11, 12, 1};
Transfinite Surface {6} = {7, 10, 1, 4};
Recombine Surface "*";

Extrude {0, 0, 2} {
  Surface{1}; Layers {2}; Recombine;
}

Recombine Surface "*";
Transfinite Volume "*";//+