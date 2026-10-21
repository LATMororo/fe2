
from pvplotter.mesh.PVMesh       import readFromGmsh 
from pvplotter.plotters.JivePlot import JivePlot

import pyvista as pv

# Read gmsh mesh and convert it to PyVista/Vtk format.

micro_mesh = readFromGmsh('rve_vf20.msh')

# Set a PyVista plotter

pl = pv.Plotter()

# ... and a JivePlot.

jpl = JivePlot(pl)

jpl.setMesh(micro_mesh) # set a mesh

jpl.addPlot('outIP.0.31.mech',     # file/path 
            tableType='NodeTable', # table type: NodeTable or ElementTable 
            dataType='strain',     # data type: strain, stress, ...
            comp='strain_xx',            # component: 'strain_xx', 'mag', ...
            scale=5.)             # scale deformation

jpl.clim([1.e-3, 5.e-2])

pl.add_scalar_bar(title='strain')
pl.add_axes()

pl.show()
