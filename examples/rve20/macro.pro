
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
    file = "cube.msh";
    
    dim = 3;

    nodeGroups = [ "left", "right", "corner" ];

    left.xtype  = "min";
    right.xtype = "max";

    corner.xval  = 0.0;
    corner.yval  = 0.0;
    corner.zval  = 0.0;
  };
};


control =
{
  pause    = 0.2;
  runWhile = "i <= 30";
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

      shape =
      {
        type      = "Hex8";
        intScheme = "Gauss1 * Gauss1 * Gauss1";
      };

      material = 
      {
        type = "FE2";  
        rank = 3;

        outIPs  = [ 0 ];

	    sampIPs = [ 0 ];

        micro = 
        {
          include "rve.pro";

          output =
          {
            vectors = [ "state = disp" ];
            tables  = [ "nodes/strain" ];
          };

          sample =
          {
            separator = ",";

            dataSets = [ "i", "micro.extra.solver.epsM[0]", 
                              "micro.extra.solver.epsM[1]",
                              "micro.extra.solver.epsM[2]",
                              "micro.extra.solver.epsM[3]",
                              "micro.extra.solver.epsM[4]",
                              "micro.extra.solver.epsM[5]",
                              "micro.extra.solver.sigM[0]",
                              "micro.extra.solver.sigM[1]",
                              "micro.extra.solver.sigM[2]",
                              "micro.extra.solver.sigM[3]",
                              "micro.extra.solver.sigM[4]",
                              "micro.extra.solver.sigM[5]" ];
          };
        };
      };
    };

    lodi =
    {
      type  = "LoadDisp";
      group = "right";
    };

    diri =
    {
      type = "Dirichlet";

      stateIncr = 0.001;
      initState = 0.0;
        
      nodeGroups = [ "left", "right", "corner", "corner" ];
      dofs       = [   "dx",   "dx" ,     "dy",     "dz" ];
      loaded     = 1;
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

    tables  = [ "elements/strain", "elements/stress" ]; 
  };

  sample =
  {
    type = "Sample";

    file = "lodi.out";

    separator = ",";

    dataSets = [ "i", "model.model.lodi.disp[0]", 
                      "model.model.lodi.load[0]" ];
  };
};
