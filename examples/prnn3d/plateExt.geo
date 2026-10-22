
//-----------------------------------------------------------------------
//   This *.geo script creats a 1/4 plate with a hole (centered), which 
//   is meshed with 3D linear 4-node elements (i.e., TET4). The main 
//   idea is to use such model as a macro mesh in a FE^2 simulation,
//   where each TET4 element has a single integration point which in
//   turn associated with an RVE.
//
//   Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
//   Date: July 2026
// 
//   Changelog:
//
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
//   Dimensions in mm
//-----------------------------------------------------------------------

ls = 2.;   // specimen length
hs = 1.;   // specimen height
ds = 0.05; // specimen depth

rc = 0.25; // radius

//-----------------------------------------------------------------------
//   Default coordinates
//-----------------------------------------------------------------------

x1 = rc;
y1 = 0.;
z1 = 0.;

x2 = ls * 0.5;
y2 = y1;
z2 = z1;

x3 = x2;
y3 = hs * 0.5;
z3 = z1;

x4 = 0.;
y4 = y3;
z4 = z1;

x5 = 0.;
y5 = rc;
z5 = z1;

x6 = 0.;
y6 = y1;
z6 = 0.;

//-----------------------------------------------------------------------
//   Characteristic element sizes
//-----------------------------------------------------------------------

// 70 HEX-8 elems...
//he = 0.20; // typical element size
//fe = 0.07; // finest mesh

// 52 HEX-8 elems...
// he = 0.22; // typical element size
// fe = 0.07; // finest mesh

// 275 HEX-8 elems...
//he = 0.100;  // typical element size
//fe = 0.025; // finest mesh

// 55 TET4 elems...
//he = 0.17; // typical element size
//fe = 0.09; // finest mesh

// 159 TET4 elems...
he = 0.18; // typical element size
fe = 0.18; // finest mesh

// 159 TET4 elems...
// he = 0.05; // typical element size
// fe = 0.05; // finest mesh

//-----------------------------------------------------------------------
//   Points
//-----------------------------------------------------------------------

Point(1)  = { x1, y1,    z1, fe };
Point(2)  = { x2, y2,    z2, he };
Point(3)  = { x3, y3,    z3, he };
Point(4)  = { x4, y4,    z4, fe };
Point(5)  = { x5, y5,    z5, fe };

Point(6)  = { x6, y6,    z6, fe };

//-----------------------------------------------------------------------
//   Lines
//-----------------------------------------------------------------------

Line(1)    = {  1,  2 };
Line(2)    = {  2,  3 };
Line(3)    = {  3,  4 };
Line(4)    = {  4,  5 };

Circle(5)  = { 5, 6, 1 };

//-----------------------------------------------------------------------
//   Curve Loops, Plane Surfaces, Extrude, and Physical Volume
//-----------------------------------------------------------------------

Curve Loop(1) = { 1, 2, 3, 4, 5 };

Plane Surface(1) = { 1 };

// To create quadrangles instead of triangles
//Recombine Surface{1};

// Extrude in z direction 
e() = Extrude {0,0,ds} {
  Surface{1}; Layers{ {1}, {1} }; 
};

// ... and extrude with hexes (quadrangles) 
//e() = Extrude {0,0,ds} {
//  Surface{1}; Layers{ {1}, {1} }; Recombine;
//};

// Whole domain (exclude other elements, e.g., lines, plane elements)
Physical Volume(1) = {e(1)};

// Create a 3D mesh
Mesh 3;
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 1;

// File format
Mesh.MshFileVersion = 2.2;

Save "plateExt.msh";
