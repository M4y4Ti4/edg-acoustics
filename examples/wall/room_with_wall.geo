Mesh.MshFileVersion = 2.2;
lc = 1.0;

Point(1)  = {  0.0,       5.98, 0.0,  lc };
Point(2)  = {  0.0,       0.0,  0.0,  lc };
Point(3)  = {  0.0,       0.0,  4.06, lc };
Point(4)  = {  0.0,       5.98, 4.06, lc };
Point(5)  = { -3.706297,  2.72, 0.0,  lc };
Point(6)  = { -3.926297,  2.72, 0.0,  lc };
Point(7)  = { -3.926297,  2.72, 4.06, lc };
Point(8)  = { -3.706297,  2.72, 4.06, lc };
Point(9)  = { -6.68,      5.98, 0.0,  lc };
Point(10) = { -6.68,      5.98, 4.06, lc };
Point(11) = { -3.706297,  0.0,  0.0,  lc };
Point(12) = { -6.68,      0.0,  0.0,  lc };
Point(13) = { -3.926297,  0.0,  0.0,  lc };
Point(14) = { -3.926297,  0.0,  4.06, lc };
Point(15) = { -3.706297,  0.0,  4.06, lc };
Point(16) = { -6.68,      0.0,  4.06, lc };

Line(1)  = {  1,  2 };
Line(2)  = {  2,  3 };
Line(3)  = {  3,  4 };
Line(4)  = {  1,  4 };
Line(5)  = {  5,  6 };
Line(6)  = {  6,  7 };
Line(7)  = {  7,  8 };
Line(8)  = {  5,  8 };
Line(9)  = {  1,  9 };
Line(10) = {  4, 10 };
Line(11) = {  9, 10 };
Line(12) = {  2, 11 };
Line(13) = {  9, 12 };
Line(14) = { 12, 13 };
Line(15) = {  6, 13 };
Line(16) = {  5, 11 };
Line(17) = { 13, 14 };
Line(18) = {  7, 14 };
Line(19) = {  3, 15 };
Line(20) = {  8, 15 };
Line(21) = { 14, 16 };
Line(22) = { 10, 16 };
Line(23) = { 14, 15 };
Line(24) = { 11, 13 };
Line(25) = { 11, 15 };
Line(26) = { 12, 16 };

Line Loop(1)  = {  1,  2,   3,  -4 };
Line Loop(2)  = {  5,  6,   7,  -8 };
Line Loop(3)  = { -9,  4,  10, -11 };
Line Loop(4)  = { -12, -1,  9,  13, 14, -15, -5, 16 };
Line Loop(5)  = { 15,  17, -18,  -6 };
Line Loop(6)  = { -3,  19, -20,  -7, 18, 21, -22, -10 };
Line Loop(7)  = { 23, -20,  -7,  18 };
Line Loop(8)  = { -24, 25, -23, -17 };
Line Loop(9)  = { 24, -15,  -5,  16 };
Line Loop(10) = { -14, 26, -21, -17 };
Line Loop(11) = { -13, 11,  22, -26 };
Line Loop(12) = { -19,  -2, 12,  25 };
Line Loop(13) = { -16,   8, 20, -25 };

Plane Surface(1)  = { 1  };
Plane Surface(2)  = { 2  };
Plane Surface(3)  = { 3  };
Plane Surface(4)  = { 4  };
Plane Surface(5)  = { 5  };
Plane Surface(6)  = { 6  };
Plane Surface(7)  = { 7  };
Plane Surface(8)  = { 8  };
Plane Surface(9)  = { 9  };
Plane Surface(10) = { 10 };
Plane Surface(11) = { 11 };
Plane Surface(12) = { 12 };
Plane Surface(13) = { 13 };

Surface Loop(1) = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
Volume(1)       = { 1 };

// Named physical groups matching BC_labels keys in wall_main.py
Physical Surface("carpet")  = { 1 };
Physical Surface("ceiling") = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
Physical Volume("air")      = { 1 };

Mesh.Algorithm                      = 6;
Mesh.Algorithm3D                    = 1;
Mesh.Optimize                       = 1;
Mesh.CharacteristicLengthFromPoints = 1;
