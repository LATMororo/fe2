#!/usr/bin/python3

#------------------------------------------------------------------------
#    Simple GmshModel Python API functionalities for generating squared
#    RVEs with embedded sphere inclusions.
#
#    Reference: GmshModel documentations
#    https://gmshmodel.readthedocs.io/en/latest/index.html
#    Accessed: July 2025
#
#    Author: Luiz Ant. T. Mororo, luiz.mororo@ifce.edu.br
#    Date:   July 2025
#    
#    Changelog:
#    
#    Author: Viny P. Barros, pereira@alu.ufc.br
#    Date: Jan 2026 
#
#    Organize the original script into functions, as well as add attempts
#    functionalities and control seeds.
#
#    Author: Luiz Ant. T. Mororo, luiz.mororo@ifce.edu.br
#    Date: June 2026
#   
#    Code clean-up and add more comments.
#------------------------------------------------------------------------


import sys
import math
import gmsh
import gmshModel
import numpy as np
import random 


#========================================================================
#    Main function (deal with Gmsh, loop fail-fast and mesh)
#========================================================================


def _gen_random_rve(a, r, n_incs, vol_frac):
    """
    It attempts to generate a random RVE.
    """

    gmAPI = gmshModel.Model
    mesh_size = 0.15 * a  
    
    initParameters={                                                           
        "inclusionSets": [r, n_incs],                                  
        "inclusionType": "Sphere",                                             
        "size": [a, a, a],                                                     
        "origin": [0, 0, 0],                                                   
        "periodicityFlags": [1, 1, 1],                                         
        "domainGroup": "matrix",                                               
        "inclusionGroup": "inclusions",                                        
        "gmshConfigChanges": {
            "General.Terminal": 0,                             
            "Mesh.CharacteristicLengthExtendFromBoundary": 0,
            "Mesh.SecondOrderLinear": 0,
            "Mesh.ElementOrder": 2,
            "Geometry.MatchMeshTolerance": 1.e-06, 
            "Geometry.Tolerance": 1.e-06,
            "Mesh.MshFileVersion": 2.2             
        }
    }

    modelingParameters = {                                                     
        "placementOptions": {
            "maxAttempts": 20000,      
            "minRelDistBnd": 0.15,     
            "minRelDistInc": 0.15,     
        }
    }

    meshingParameters={                                                        
        "threads": None,                                                       
        "refinementOptions": {
            "maxMeshSize": mesh_size,                            
            "inclusionRefinement": True,                     
            "interInclusionRefinement": False,                
            "elementsPerCircumference": 15,  
            "elementsBetweenInclusions": 3,                  
            "inclusionRefinementWidth": 3,                   
            "transitionElements": "auto",                    
            "aspectRatio": 1.                    
        }
    }

    MAX_ATTEMPTS = 20
    tentativa = 0
    attempt = 0
    sucess = False

    while attempt < MAX_ATTEMPTS:
        attempt += 1
        print(f"\n--- Attempt {attempt}/{MAX_ATTEMPTS} ---")
        
        testRVE = gmAPI.RandomInclusionRVE(**initParameters)
        
        try:
            # 1. Geometry settings
            testRVE.defineGeometricObjects(**modelingParameters)
            testRVE.addGeometricObjectsToGmshModel()
            
            # --- FAIL-FAST ---
            gmsh.model.occ.synchronize()
            num_volumes_3d = len(gmsh.model.occ.getEntities(3))
            num_spheres = num_volumes_3d - 1 # minus matrix/cube 
            
            if num_spheres < n_incs:
                print(f"[WARN] Only {num_spheres}/{n_incs} inclusions were inserted.")
                print("Changing seeds...")
                testRVE.close()
                continue 
                
            print(f"[OK] {n_incs}/{n_incs} embedded inclusions. " \
                                  "Initializing boolean operation...")
            
            # 2. Boolean operations
            testRVE.defineBooleanOperations()                                          
            testRVE.performBooleanOperationsForGmshModel()                             

            # ... and define DEFAULT physical groups
            testRVE.definePhysicalGroups()                                             
            
            # ... and then, remove 'boundary' group
            # NOTE: this is necessary for better reading task of JIVE's 
            #       module 'GmshInputModule'
            del testRVE.groups['boundary']
            testRVE.physicalGroups.pop()
            
            testRVE.addPhysicalGroupsToGmshModel()                                     
            #NOTE: in order to remove `$PhysicalNames' names, run:
            testRVE.gmshAPI.removePhysicalName('matrix')
            testRVE.gmshAPI.removePhysicalName('inclusions')

            # 3. Set up periodicity constraints
            testRVE.setupPeriodicity()                                                 
            
            # 4. Mesh generation
            print("Valid geometry: meshing...") 
            testRVE.createMesh(**meshingParameters)

            # ... and show RVE's mesh
            testRVE.visualizeMesh()
            
            # 5. Save mesh file 
            fname = f'rve_vf{int(vol_frac)}.msh2'
            testRVE.saveMesh(fname)
            print(f"[SUCESS] RVE mesh saved as: {fname}")
            
            testRVE.close()
            sucess = True
            break 
            
        except Exception as e:
            print(f"[MESH ERROR] Gmsh error: {e}")
            print("Geometry overlapping during boolean operation. Re-init model...")
            testRVE.close()

    if not sucess:
        print(f"\n[FATAL ERROR] Unable of generating Vf={vol_frac}% after {MAX_ATTEMPTS} attempts")


#========================================================================
#    N INCLUSIONS + FIXED RADIUS
#========================================================================


def gen_random_fixed_radius(radius, n_incs, vol_fracs):
    """
    Keep radii of inclusions fixed, and compute RVE's edges (a).
    """
    print(f"\n{'='*60}\n RANDOM GENERATION: FIXED RADIUS (r = {radius:.5f})\n{'='*60}")
    
    sphere_vol = (4.0 / 3.0) * math.pi * (radius**3) # vol. of a single sphere
    total_incs_vol = n_incs * sphere_vol # total vol. of all inclusions

    for vol_frac in vol_fracs:
        a = (total_incs_vol / (vol_frac / 100.0)) ** (1.0 / 3.0)

        print(f'\n>>> VFrac: {vol_frac}% | Computed RVE size (a): {a:.5e} mm')
        _gen_random_rve(a, radius, n_incs, vol_frac)


#========================================================================
#    N INCLUSIONS + FIXED RVE EDGE
#========================================================================


def gen_random_fixed_edge(edge, n_incs, vol_fracs):
    """
    Keep a fixed RVE (fixed edge), and adjust radii of inclusions based on
    volume fraction. NOTE: all the inclusions have the same radius.
    """
    print(f"\n{'='*60}\n RANDOM GENERATION: FIXED EDGE (a = {edge:.2f})\n{'='*60}")
    
    total_rve_vol = edge ** 3.0

    for vol_frac in vol_fracs:
        # Volume of a single inclusion
        single_inc_vol = (total_rve_vol * (vol_frac / 100.0)) / n_incs
        
        # Compute the radius of a single inclusions (sphere)
        r = ((3.0 * single_inc_vol) / (4.0 * math.pi)) ** (1.0 / 3.0)

        print(f'\n>>> Vfrac: {vol_frac}% | Computed radius (r): {r:.5e} mm')
        _gen_random_rve(edge, r, n_incs, vol_frac)


#========================================================================
#    MAIN FUNCTION/CALL
#========================================================================


if __name__ == "__main__":

    seed = 5
    random.seed(seed)
    np.random.seed(seed)
    
    vol_fracs = [10, 20, 30] # [%]
    vol_fracs = [20]
    n_incs    = 5 # nbr. of spheres
    
    #-----------------------------------------------------------
    # Fixed radius - a radius is given
    #-----------------------------------------------------------
    radius = 0.005 / 2.0 # [mm]
    gen_random_fixed_radius(radius, n_incs, vol_fracs)
    
    #-----------------------------------------------------------
    # Fixed edge - an edge is given, which usually equals to one
    #-----------------------------------------------------------
    edge = 1.0
    #gen_random_fixed_edge(edge, n_incs, vol_fracs)
