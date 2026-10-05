
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
    file = "macro.msh";

    nodeGroups = [ "left", "right", "mid" ];

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
  runWhile = "i < 55";
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

    models = [ "bulk", "lodi", "diri" ];

    bulk = 
    {
      type = "Solid";

      elements = "all";

      // Compare to the corresponding `linear' problem, let's increase
      // the nbr. of threads. It requires more `iterations'!
      threads  = 8;

      shape =
      {
        type      = "Quad4";
        intScheme = "Gauss1 * Gauss1";
      };

      material = 
      {
        type = "FE2";  
        rank = 2;

        state = "PlaneStrain";

        outIPs = [ 1, 2, 25, 26 ];

        micro = 
        {
          include "rve.pro";

          output =
          {
            vectors = [ "state = disp" ];
            tables  = [ "nodes/strain" ];
          };
        };
      };

      thickness = 1.;
    };

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
    type = "Nonlin"; 

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
