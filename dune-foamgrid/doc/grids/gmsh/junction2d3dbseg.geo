Point(1) = { 0, 0, 0, 10.0};
Point(2) = { 0, 0, 1, 10.0};
Point(3) = {-1, 0, 0, 10.0};
Point(4) = { 1, 0, 1, 10.0};
Point(5) = { 1, 1, 1, 10.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {1, 4};
Line(5) = {4, 2};
Line(6) = {2, 5};
Line(7) = {5, 1};

Line Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};
Line Loop(2) = {-1, 4, 5};
Plane Surface(2) = {2};
Line Loop(3) = {1, 6, 7};
Plane Surface(3) = {3};
Physical Surface(1) = {1:3};

Physical Line(42) = {2,3,4};
