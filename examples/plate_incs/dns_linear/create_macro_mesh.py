#------------------------------------------------------------------------
#    Create a mesh with an array of inclusions embedded in a matrix by
#    using Gmsh Python API functionalities.
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

gmsh.model.add('rec_inclusions')

# Parametric data.

b = 12.   # length
h = 4.    # width
lc = 0.05 # typical element size
fi = 0.05

rc = 0.125 # radius of inclusions

n_inc_b = 12  # number of inclusions along the length
n_inc_h = 4   # number of inclusions along the width 

n_nodes_b = int(b//lc)
n_nodes_h = int(h//lc) 

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

factory.addPoint(0,h/2,0, lc, 5)

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
# Create an array of inclusions.
#------------------------------------------------------------------------

dx = 0.5 * b / n_inc_b
dy = 0.5 * h / n_inc_h

print('delta_x', dx)
print('delta_y', dx)


circ_tags = [] # tags of circle arcs that form an inclusion

for iy in range(n_inc_h):
    for ix in range(n_inc_b):
        xc = 2*ix*dx+dx
        yc = 2*iy*dy+dy
        print(2*ix*dx+dx, 2*iy*dy+dy)

        # Create points for two circle arcs.

        pt_tags = []

        pt_tags.append(factory.addPoint(xc, yc   , 0, fi))
        pt_tags.append(factory.addPoint(xc, yc-rc, 0, fi))
        pt_tags.append(factory.addPoint(xc, yc+rc, 0, fi))

        # Create the circle arcs.

        circ_tags.append(factory.addCircleArc(pt_tags[1], pt_tags[0], pt_tags[2]))
        circ_tags.append(factory.addCircleArc(pt_tags[2], pt_tags[0], pt_tags[1]))

        #factory.mesh.setTransfiniteCurve(circ_tags[-1], n_nodes_c) 
        #factory.mesh.setTransfiniteCurve(circ_tags[-2], n_nodes_c) 

#print(circ_tags)

#------------------------------------------------------------------------
# Define curve loops 
#------------------------------------------------------------------------

# Curve loop of the matrix domain.

factory.addCurveLoop([1, 2, 3, 4] + circ_tags, 1) 

# Curve loop of all the inclusions.

inc_loop_tags = []

for i in range(len(circ_tags)//2):
    inc_loop_tags.append(factory.addCurveLoop([circ_tags[2*i], circ_tags[2*i+1]]))

#print(inc_loop_tags)

#------------------------------------------------------------------------
# Define plane surfaces 
#------------------------------------------------------------------------

# Plane surface for the matrix.

factory.addPlaneSurface([1], 1)

# Plane surface for the inclusions.

inc_surf_tags = []

for isurf in inc_loop_tags:
    inc_surf_tags.append(factory.addPlaneSurface([isurf]))

#------------------------------------------------------------------------
# Synchronize Gmsh API and create physical groups 
#------------------------------------------------------------------------

factory.synchronize()

gmsh.model.addPhysicalGroup(2, [1], 1)
gmsh.model.addPhysicalGroup(2, inc_surf_tags, 2)

#------------------------------------------------------------------------
# Generating a 2D mesh and saving it ...
#------------------------------------------------------------------------

gmsh.model.mesh.generate(2)

gmsh.option.setNumber("Mesh.MshFileVersion", 2)

gmsh.write('rec_inc.msh')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Finalise when it is done using the Gmsh Python API.

gmsh.finalize()
