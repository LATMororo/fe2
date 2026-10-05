#------------------------------------------------------------------------
#    Create a rectangular domain meshed with 4-node finite element (Q4)
#    by using Gmsh Python API functionalities.
#
#    Reference: Gmsh documentations
#    https://gmsh.info/doc/texinfo/gmsh.html#t1. Accessed: May 2023
#
#    Author: Luiz Ant. T. Mororo
#    Date: Aug 2025
#------------------------------------------------------------------------

import gmsh
import sys

# Before using any functions in the Python API, Gmsh must be initialized.

gmsh.initialize()

gmsh.model.add('macro')

# Parametric data.

b = 1.     # length
h = 1.     # width
lc = 0.025 # typical element size
fi = 0.05

rc = 0.125 # radius of inclusions

xc = b * 0.5
yc = h * 0.5

n_nodes_b = 21 # number of nodes along the length
n_nodes_h = 21 # number of nodes along the width 

# Alias for 'gmsh.model.geo'.
# NOTE: one can make use of 'gmsh.model.occ' instead.

factory = gmsh.model.geo

#------------------------------------------------------------------------
#    Points
#------------------------------------------------------------------------

factory.addPoint(0,      0, 0, lc, 1)
# it could be done like, p1 = factory.addPoint(0, 0, 0, lc)
factory.addPoint(b,      0, 0, lc, 2)
factory.addPoint(b,      h, 0, lc, 3)
factory.addPoint(0,      h, 0, lc, 4)

factory.addPoint(xc, yc-rc, 0, lc, 5) 
factory.addPoint(xc,    yc, 0, lc, 6) # inclusion center
factory.addPoint(xc, yc+rc, 0, lc, 7) 

#------------------------------------------------------------------------
#    Lines (or Curves)
#------------------------------------------------------------------------

factory.addLine(1, 2, 1)
# it could also be done like: factory.addLine(p1, 2, 1)
factory.addLine(2, 3, 2)
factory.addLine(3, 4, 3)
factory.addLine(4, 1, 4)

factory.addCircleArc(5, 6, 7, 5)
factory.addCircleArc(7, 6, 5, 6)

#factory.mesh.setTransfiniteCurve(1, n_nodes_b) 
#factory.mesh.setTransfiniteCurve(2, n_nodes_h) 
#factory.mesh.setTransfiniteCurve(3, n_nodes_b) 
#factory.mesh.setTransfiniteCurve(4, n_nodes_h) 

#factory.mesh.setTransfiniteCurve(5, n_nodes_c) 
#factory.mesh.setTransfiniteCurve(6, n_nodes_c) 



#------------------------------------------------------------------------
#    Surfaces
#------------------------------------------------------------------------

factory.addCurveLoop([1, 2, 3, 4, 5, 6], 1)
factory.addCurveLoop([5, 6], 2)


factory.addPlaneSurface([1], 1)
factory.addPlaneSurface([2], 2)

#factory.mesh.setTransfiniteSurface(1)
#factory.mesh.setTransfiniteSurface(2)


gmsh.model.geo.synchronize()


gmsh.model.addPhysicalGroup(2, [1], 1)
gmsh.model.addPhysicalGroup(2, [2], 2)

#------------------------------------------------------------------------
# Generating a 2D mesh and saving it ...
#------------------------------------------------------------------------

gmsh.model.mesh.setPeriodic(1, [1], [2], [1, b, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1])
gmsh.model.mesh.setPeriodic(1, [3], [4], [1, 0, 0, 0,  0, 1, h, 0,  0, 0, 1, 0,  0, 0, 0, 1])

#gmsh.model.mesh.setRecombine(2, 1)
gmsh.model.mesh.generate(2)

gmsh.option.setNumber("Mesh.MshFileVersion", 2)

gmsh.write('rve.msh')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Finalise when it is done using the Gmsh Python API.

gmsh.finalize()
