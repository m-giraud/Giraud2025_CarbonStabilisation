// SPDX-FileCopyrightInfo: Copyright (C) DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
cl1 = 1;
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {0, 1, 0, cl1};
Point(4) = {0, 0, 1, cl1};
Point(5) = {-1, 1, 1, cl1};
Line(1) = {1, 2};
Line(2) = {1, 3};
Line(3) = {1, 4};
Line(4) = {2, 3};
Line(5) = {2, 4};
Line(6) = {3, 4};
Line(7) = {1, 5};
Line(8) = {3, 5};
Line(9) = {4, 5};
Line Loop(11) = {7, -9, -3};
Plane Surface(11) = {11};
Line Loop(13) = {6, 9, -8};
Plane Surface(13) = {13};
Line Loop(15) = {8, -7, 2};
Plane Surface(15) = {15};
Line Loop(17) = {1, 5, -3};
Plane Surface(17) = {17};
Line Loop(19) = {4, 6, -5};
Plane Surface(19) = {19};
Line Loop(21) = {1, 4, -2};
Plane Surface(21) = {21};
Line Loop(23) = {2, 6, -3};
Plane Surface(23) = {23};
Surface Loop(25) = {15, 13, 11, 23};
Volume(25) = {25};
Surface Loop(27) = {21, 17, 19, 23};
Volume(27) = {27};
Physical Surface(30) = {15, 21};
Physical Surface(31) = {11, 17};
Physical Surface(32) = {13, 19};
Physical Volume(28) = {25};
Physical Volume(29) = {27};
