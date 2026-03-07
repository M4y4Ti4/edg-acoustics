Mesh.MshFileVersion = 2.2;
lc = 1.5;

// Shoebox Room: 5m (x) x 4m (y) x 3m (z)
// Matches Room.create_shoebox([5, 4, 3], ...)

// --- Points ---
// Bottom face (z = 0)
Point(1)  = { 0.0, 0.0, 0.0, lc };
Point(2)  = { 5.0, 0.0, 0.0, lc };
Point(3)  = { 5.0, 4.0, 0.0, lc };
Point(4)  = { 0.0, 4.0, 0.0, lc };

// Top face (z = 3)
Point(5)  = { 0.0, 0.0, 3.0, lc };
Point(6)  = { 5.0, 0.0, 3.0, lc };
Point(7)  = { 5.0, 4.0, 3.0, lc };
Point(8)  = { 0.0, 4.0, 3.0, lc };

// --- Lines ---
// Bottom edges
Line(1)  = { 1, 2 };
Line(2)  = { 2, 3 };
Line(3)  = { 3, 4 };
Line(4)  = { 4, 1 };

// Top edges
Line(5)  = { 5, 6 };
Line(6)  = { 6, 7 };
Line(7)  = { 7, 8 };
Line(8)  = { 8, 5 };

// Vertical edges
Line(9)  = { 1, 5 };
Line(10) = { 2, 6 };
Line(11) = { 3, 7 };
Line(12) = { 4, 8 };

// --- Line Loops & Surfaces ---

// Floor (z = 0) — carpet
Line Loop(1)    = { 1, 2, 3, 4 };
Plane Surface(1) = { 1 };

// Ceiling (z = 3) — drywall
Line Loop(2)    = { 5, 6, 7, 8 };
Plane Surface(2) = { 2 };

// Front wall (y = 0) — brick
Line Loop(3)    = { 1, 10, -5, -9 };
Plane Surface(3) = { 3 };

// Back wall (y = 4) — brick
Line Loop(4)    = { -3, 11, 7, -12 };
Plane Surface(4) = { 4 };

// Left wall (x = 0) — brick
Line Loop(5)    = { -4, 12, 8, -9 };  // corrected winding
Plane Surface(5) = { 5 };

// Right wall (x = 5) — brick
Line Loop(6)    = { 2, 11, -6, -10 };
Plane Surface(6) = { 6 };

// --- Volume ---
Surface Loop(1) = { 1, 2, 3, 4, 5, 6 };
Volume(1)       = { 1 };

// --- Physical Groups ---
// Matching materials from Python script:
//   Physical Surface(13) = carpet (floor)
//   Physical Surface(11) = hard surface (ceiling = drywall)
//   Physical Surface(14) = panel/brick walls

Physical Surface(13) = { 1 };        // floor   — carpet
Physical Surface(11) = { 2 };        // ceiling — drywall
Physical Surface(14) = { 3, 4, 5, 6 }; // walls   — brick

Physical Volume(1) = { 1 };

// --- Mesh Settings ---
Mesh.Algorithm   = 1;
Mesh.Algorithm3D = 1; // Delaunay3D
Mesh.Optimize    = 1;
Mesh.CharacteristicLengthFromPoints = 1;
