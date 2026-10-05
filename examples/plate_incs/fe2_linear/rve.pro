
userInput =
{
  modules = [ "mesh", "pbcgroups" ];

  mesh = 
  {
    type = "GmshInput";
    file = "rve.msh";

    eGrpIdx = [ 1, 2 ];
  };

  pbcgroups = 
  {
    type = "PBCGroupInput";
  };
};


model =
{
  type = "Matrix";

  matrix0 = 
  {
    type      = "FEM";
    symmetric = true;
  };

  model = 
  {
    type = "Multi";

    models = [ "matrix", "inc" ];

    matrix = 
    {
      type = "Solid";

      elements = "gmsh1";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };

      material = 
      {
        type = "Isotropic";
        rank = 2;
        
        state = "PlaneStrain";
        
        young	= 79500.;
        poisson	= 0.33;
      };

      thickness = 1.;
    };

    inc = matrix;
    inc.elements       = "gmsh2";
    inc.material.young   = 410000.;
    inc.material.poisson = 0.2;
  };
};


extra =
{
  modules = [ "solver" ];

  solver =
  {
    type = "MicroPBC";

    // Solver for equilibrium problem

    esolver =
    {
      type = "Nonlin"; // one could use `Linsolve'

      solver = 
      {
        type = "GMRES";
        precon.type = "ILUd";
      };
    };

    // Solver for micro-problem, i.e., homogenization

    msolver =
    {
      type = "SkylineLU";
    };
  };
};

