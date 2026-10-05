
log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};


userInput =
{
  modules = [ "mesh" ];

  mesh = 
  {
    type = "GmshInput";
    file = "rec_inc.msh";

    nodeGroups = [ "left", "right", "mid" ];

    eGrpIdx     = [ 1, 2 ];

    left.xtype  = "min";
    right.xtype = "max";

    mid.xtype   = "min";
    mid.ytype   = "mid";
    mid.eps     = 1.e-7;
  };
};


control =
{
  pause    = 0.2;
  runWhile = "i < 11";
  fgMode   = false;
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

    models = [ "matrix", "inc", "lodi", "diri" ];

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
    inc.elements         = "gmsh2";
    inc.material.young   = 410000.;
    inc.material.poisson = 0.2;

    lodi =
    {
      type  = "LoadDisp";
      group = "right";
    };

    diri =
    {
      type = "Dirichlet";

      stateIncr = 0.002;
      initState = 0.0;
        
      nodeGroups = [ "left", "mid", "right" ];
      dofs       = [   "dx",  "dy",    "dx" ];
      loaded     = 2;
    };
  };
};


extra =
{
  modules = [ "solver", "output", "sample" ];

  solver =
  {
    type = "Nonlin"; // `Linsolve would be enough'

    solver =
    {
      type = "GMRES";

      precon =
      {
        type = "ILUd";
      };
    };
  };

  output =
  {
    type = "Output";

    file = "$(CASE_NAME)_out.%i.mech";

    vectors = [ "state = disp" ];

    tables  = [ "nodes/strain", "nodes/stress" ]; 
  };

  sample =
  {
    type = "Sample";

    file = "lodi.out";

    separator = ",";

    dataSets  = [ "i", "model.model.lodi.disp[0]", 
                       "model.model.lodi.load[0]" ];
  };
};
