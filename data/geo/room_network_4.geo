SetFactory("OpenCASCADE");
Mesh.Algorithm = 8;
l = 1;
main_lx = 10;
main_ly = 1;
sub_ly = 2;
sub_lx = 1;
Point(1)  = {0, 0, 0,l};
Point(2)  = {main_lx, 0, 0,l};
Point(3)  = {main_lx, main_ly, 0,l};
x_offset = 1;
Point(5)  = {main_lx-x_offset, main_ly+sub_ly, 0,l};
Point(6)  = {main_lx-x_offset-sub_lx, main_ly+sub_ly, 0,l};
x_offset = 3;
Point(9)  = {main_lx-x_offset, main_ly+sub_ly, 0,l};
Point(10)  = {main_lx-x_offset-sub_lx, main_ly+sub_ly, 0,l};
x_offset = 5;
Point(13)  = {main_lx-x_offset, main_ly+sub_ly, 0,l};
Point(14)  = {main_lx-x_offset-sub_lx, main_ly+sub_ly, 0,l};
x_offset = 7;
Point(17)  = {main_lx-x_offset, main_ly+sub_ly, 0,l};
Point(18)  = {main_lx-x_offset-sub_lx, main_ly+sub_ly, 0,l};
Point(20) = {0, main_ly, 0, l};
Point(21) = {4.25, 1, 0, 1.0};
Point(22) = {4.25, 1.25, 0, 1.0};
Point(23) = {4, 1.25, 0, 1.0};
Point(24) = {4.75, 1, 0, 1.0};
Point(25) = {4.75, 1.25, 0, 1.0};
Point(26) = {5, 1.25, 0, 1.0};
Point(27) = {8.75, 1, 0, 1.0};
Point(28) = {8.75, 1.25, 0, 1.0};
Point(29) = {9, 1.25, 0, 1.0};
Point(30) = {8.25, 1, 0, 1.0};
Point(31) = {8.25, 1.25, 0, 1.0};
Point(32) = {8, 1.25, 0, 1.0};
Point(33) = {6.75, 1, 0, 1.0};
Point(34) = {6.75, 1.25, 0, 1.0};
Point(35) = {7, 1.25, 0, 1.0};
Point(36) = {6.25, 1, 0, 1.0};
Point(37) = {6.25, 1.25, 0, 1.0};
Point(38) = {6, 1.25, 0, 1.0};
Point(40) = {10, 1, 0, 1.0};
Point(41) = {8.25, 1.25, 0, 1.0};
Point(42) = {8.25, 1, 0, 1.0};
Point(43) = {2.75, 1, 0, 1.0};
Point(44) = {2.75, 1.25, 0, 1.0};
Point(45) = {3, 1.25, 0, 1.0};
Point(47) = {2, 1.25, 0, 1.0};
Point(48) = {2.25, 1.25, -0, 1.0};
Point(49) = {2.25, 1, -0, 1.0};

Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 27};
Line(5) = {27, 28};
Line(6) = {28, 29};
Line(7) = {29, 5};
Line(8) = {5, 6};
Line(9) = {6, 32};
Line(10) = {32, 31};
Line(11) = {31, 30};
Line(12) = {30, 33};
Line(13) = {33, 34};
Line(14) = {34, 35};
Line(15) = {35, 9};
Line(16) = {9, 10};
Line(17) = {10, 38};
Line(18) = {38, 37};
Line(19) = {37, 36};
Line(20) = {36, 24};
Line(21) = {24, 25};
Line(22) = {25, 26};
Line(23) = {26, 13};
Line(24) = {13, 14};
Line(25) = {14, 23};
Line(26) = {23, 22};
Line(27) = {22, 21};
Line(28) = {21, 43};
Line(29) = {43, 44};
Line(30) = {44, 45};
Line(31) = {45, 17};
Line(32) = {17, 18};
Line(33) = {18, 47};
Line(34) = {47, 48};
Line(35) = {48, 49};
Line(36) = {49, 20};
Line(37) = {20, 1};


// Length = 0.25
Transfinite Curve {5,6,10,11,13,14,18,19,21,22,26,27,29,30,34,35} = 1 Using Progression 1;
// Length = 1
Transfinite Curve {3,8,16,24,32,37} = 5 Using Progression 1;
//Length = 1.25
Transfinite Curve {4} = 6 Using Progression 1;
// Length = 1.5
Transfinite Curve {12,20,28} = 7 Using Progression 1;
//Length = 1.75
Transfinite Curve {7,9,15,17,23,25,31,33} = 8 Using Progression 1;
// Length = 2.25
Transfinite Curve {36} = 10 Using Progression 1;
// Length = 10
Transfinite Curve {2} = 41 Using Progression 1;

Curve Loop(1) = {33, 34, 35, 36, 37, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
Plane Surface(1) = {1};

Transfinite Surface {2} = {18, 47, 45, 17};
Transfinite Surface {3} = {14, 23, 26, 13};
Transfinite Surface {4} = {10, 38, 35, 9};
Transfinite Surface {5} = {6, 32, 29, 5};
Transfinite Surface {6} = {48, 49, 43, 44};
Transfinite Surface {7} = {22, 21, 24, 25};
Transfinite Surface {8} = {37, 36, 33, 34};
Transfinite Surface {9} = {31, 30, 27, 28};
Transfinite Surface {10} = {20, 1, 2, 3};

Recombine Surface "*";

Extrude {0, 0, 1} {
    Surface{1}; Layers {4}; Recombine;
}

Recombine Surface "*";
Transfinite Volume "*";