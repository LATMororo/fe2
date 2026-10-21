//-----------------------------------------------------------------------
//   Dimensions in mm
//-----------------------------------------------------------------------

ls = 1.; // specimen length
hs = 1.; // specimen height
ds = 1.; // specimen depth 

//-----------------------------------------------------------------------
//   Default coordinates
//-----------------------------------------------------------------------

x1 = 0.;
y1 = 0.;
z1 = 0.;

x2 = x1 + ls;
y2 = y1;
z2 = z1;

x3 = x2;
y3 = hs;
z3 = z1;

x4 = x1;
y4 = hs;
z4 = z1;

//-----------------------------------------------------------------------
//   Characteristic element sizes
//-----------------------------------------------------------------------

he = 1.; // typical element size

//-----------------------------------------------------------------------
//   Points
//-----------------------------------------------------------------------

Point(1)  = { x1, y1,    z1, he };
Point(2)  = { x2, y2,    z2, he };
Point(3)  = { x3, y3,    z3, he };
Point(4)  = { x4, y4,    z4, he };

Point(5)  = { x1, y1, z1+ds, he };
Point(6)  = { x2, y2, z2+ds, he };
Point(7)  = { x3, y3, z3+ds, he };
Point(8)  = { x4, y4, z4+ds, he };

//-----------------------------------------------------------------------
//   Lines
//-----------------------------------------------------------------------

Line(1)   = {  1,  2 };
Line(2)   = {  2,  3 };
Line(3)   = {  3,  4 };
Line(4)   = {  4,  1 };

Line(5)   = {  5,  6 };
Line(6)   = {  6,  7 };
Line(7)   = {  7,  8 };
Line(8)   = {  8,  5 };

Line(9)   = {  2,  6 };
Line(10)  = {  7,  3 };

Line(11)  = {  1,  5 };
Line(12)  = {  8,  4 };


//-----------------------------------------------------------------------
//   Line Loops, Plane Surfaces, Surface Loop and Physical Surfaces
//-----------------------------------------------------------------------

Line Loop(1) = { -1,  -4, -3,  -2 };
Line Loop(2) = {  5,   6,  7,   8 };

Line Loop(3) = {  2, -10, -6,  -9 };
Line Loop(4) = { 11,  -8,  12,  4 };

Line Loop(5) = { 10,   3, -12, -7 };
Line Loop(6) = {  9,  -5, -11,  1 };

Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };
Plane Surface(3) = { 3 };
Plane Surface(4) = { 4 };
Plane Surface(5) = { 5 };
Plane Surface(6) = { 6 };

Surface Loop(1)  = { 1, 2, 3, 4, 5, 6 };

Volume(1)        = { 1 };

sizex = 2;
sizey = 2;
sizez = 2;

Transfinite Line { 1,  3,  5,  7 } = sizex;
Transfinite Line { 2,  4,  6,  8 } = sizey;
Transfinite Line { 9, 10, 11, 12 } = sizez;

Transfinite Surface(1) = { 4, 3, 2, 1 };
Transfinite Surface(2) = { 5, 6, 7, 8 };
Transfinite Surface(3) = { 2, 3, 7, 6 };
Transfinite Surface(4) = { 5, 8, 4, 1 };
Transfinite Surface(5) = { 3, 4, 8, 7 };
Transfinite Surface(6) = { 1, 2, 6, 5 };

Recombine Surface { 1, 2, 3, 4, 5, 6 };

Transfinite Volume {1} = { 8, 4, 3, 7, 5, 1, 2, 6 };

//Physical Surface(200)  = { 1, 2, 3, 4, 5, 6 };

//Mesh.SecondOrderLinear = 1;
//Mesh.ElementOrder = 2;

Mesh.SecondOrderLinear = 1;

Physical Volume(100)   = { 1 };

Mesh 3;
Mesh.MshFileVersion = 2.2;





