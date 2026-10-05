
from pvplotter.mesh.PVMesh       import readFromGmsh 
from pvplotter.plotters.JivePlot import JivePlot

import pyvista as pv

# Read gmsh mesh and convert it to PyVista/Vtk format.

dns_mesh   = readFromGmsh('dns_j2/rec_inc.msh')
macro_mesh = readFromGmsh('fe2_j2/macro.msh')
micro_mesh = readFromGmsh('fe2_j2/rve.msh')

# Read DNS, macro, and micro strain fields

dns_table   = JivePlot.findTable('dns_j2/macro_out.40.mech', 
                                  tableType='NodeTable', dataType='strain')

macro_table = JivePlot.findTable('fe2_j2/macro_out.40.mech', 
                                  tableType='NodeTable', dataType='strain')

micro_table = JivePlot.findTable('fe2_j2/outIP.25.40.mech', 
                                  tableType='NodeTable', dataType='strain')

# ... and get only the component strain_xx.

dns_mesh.point_data['dns_strain']     = [float(d[1]) for d in dns_table[1]]
macro_mesh.point_data['macro_strain'] = [float(d[1]) for d in macro_table[1]]
micro_mesh.point_data['micro_strain'] = [float(d[1]) for d in micro_table[1]]

# PyVista plotter.

pl = pv.Plotter(shape='1/2')

pl.subplot(2)
mactor = pl.add_mesh(dns_mesh, 
                     show_edges=True, 
                     scalars='dns_strain', 
                     clim=[0., 6.5e-3])

pl.camera_position='xy'


pl.subplot(1)
mactor = pl.add_mesh(micro_mesh, 
                     show_edges=True,
                     scalars='micro_strain', # RVE strain
                     clim=[0., 6.5e-3])

pl.camera_position='xy'


pl.subplot(0)
mactor = pl.add_mesh(macro_mesh, 
                     show_edges=True, 
                     scalars='macro_strain', 
                     clim=[0., 6.5e-3])

pl.camera_position='xy'

pl.show()
