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

b = 12.   # length
h = 4.    # width
lc = 0.5  # typical element size

n_nodes_b = 13 # number of nodes along the length
n_nodes_h = 5  # number of nodes along the width 

# Alias for 'gmsh.model.geo'.
# NOTE: one can make use of 'gmsh.model.occ' instead.

factory = gmsh.model.geo

#------------------------------------------------------------------------
# Create a rectangle b x h size
#------------------------------------------------------------------------

factory.addPoint(0, 0, 0, lc, 1)
factory.addPoint(b, 0, 0, lc, 2)
factory.addPoint(b, h, 0, lc, 3)
factory.addPoint(0, h, 0, lc, 4)

factory.addLine(1, 2, 1)
factory.addLine(2, 3, 2)
factory.addLine(3, 4, 3)
factory.addLine(4, 1, 4)

# ... and make sure they have the same number of nodes on each edge.

factory.mesh.setTransfiniteCurve(1, n_nodes_b) 
factory.mesh.setTransfiniteCurve(2, n_nodes_h) 
factory.mesh.setTransfiniteCurve(3, n_nodes_b) 
factory.mesh.setTransfiniteCurve(4, n_nodes_h) 

#------------------------------------------------------------------------
# Define curve loops 
#------------------------------------------------------------------------

# Curve loop of the matrix domain.

factory.addCurveLoop([1, 2, 3, 4], 1) 

#------------------------------------------------------------------------
# Define plane surfaces 
#------------------------------------------------------------------------

factory.addPlaneSurface([1], 1)

factory.mesh.setTransfiniteSurface(1)

#------------------------------------------------------------------------
# Synchronize Gmsh API and create physical groups 
#------------------------------------------------------------------------

factory.synchronize()

gmsh.model.addPhysicalGroup(2, [1], 1)

#------------------------------------------------------------------------
# Generating a 2D mesh and saving it ...
#------------------------------------------------------------------------

gmsh.model.mesh.setRecombine(2, 1)
gmsh.model.mesh.generate(2)

gmsh.option.setNumber("Mesh.MshFileVersion", 2)

gmsh.write('macro.msh')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Finalise when it is done using the Gmsh Python API.

gmsh.finalize()
