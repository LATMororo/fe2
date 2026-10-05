
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
        type = "J2";
        rank = 2;
        
        state = "PlaneStrain";
        
        young	= 79500.;
        poisson	= 0.33;

        yield    = "90.0740 * ( 1. - exp(-113.0268 * x)) + 444.1184"; // exp.
      	//yield    = "486.3960 + 794.0059 * x"; // linear
      };

      thickness = 1.;
    };

    inc = 
    {
      type = "Solid";

      elements = "gmsh2";

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
        
        young   = 410000.;
        poisson = 0.2;
      };

      thickness = 1.;
    };
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
      type = "Nonlin"; 

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

