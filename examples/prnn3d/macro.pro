
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
    file = "plateExt.msh";
    
    dim = 3;

    nodeGroups = [ "left", "right", "bottom", "zface" ];

    left.xtype   = "min";
    right.xtype  = "max";
    bottom.ytype = "min";

    //zface.xtype  = "min";
    zface.ztype  = "min";
  };
};


control =
{
  pause    = 0.2;
  runWhile = "i <= 25";
  fgMode   = false;
};


model =
{
  type = "Matrix";

  matrix0 = 
  {
    type      = "FEM";
    symmetric = false;
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
        type      = "Tet4";
        intScheme = "Gauss1";
      };

      material = 
      {
        type = "PRNN";
        rank = 3;
        weights = "prnn_composite_loading_5_dense.pth";
        module = "pyprnn.checkpoint";
        loader = "load_model";
        //pythonPath = "../../python";

        material =
        {
          type        = "J2";
          rank        = 3;
          young       = 79500.; 
          poisson     = 0.33; 
          yield       = "90.0740 * ( 1. - exp(-113.0268 * x)) + 444.1184";
          rmTolerance = 1.e-6; 
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

      stateIncr = 0.0005;
      initState = 0.0;
        
      nodeGroups = [ "left", "right", "bottom", "zface" ];
      dofs       = [   "dx",   "dx" ,     "dy",    "dz" ];

      //nodeGroups = [ "left", "right", "bottom" ];
      //dofs       = [   "dx",   "dx" ,     "dy" ];
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

    solver.type = "GMRES";
    solver.precon.type = "SparseLU";

    //solver.type = "SparseLU";
  };

  /*
  graph =
  {
    type = "Graph";
    dataSets = "loadDisp";
    loadDisp =
    {
      key = "Load-displacement curve";
      xData = "model.model.lodi.disp[0]";
      yData = "model.model.lodi.load[0]";
    };
  };
  */

  output =
  {
    type = "Output";

    file = "$(CASE_NAME)_out.%i.mech";

    vectors = [ "state = disp" ];

    tables  = [ "elements/strain", "elements/stress", "nodes/strain" ]; 
  };

  sample =
  {
    type = "Sample";

    file = "lodi.out";

    separator = ",";
    dataSets  = [ "i", "model.model.lodi.disp[0]", "model.model.lodi.load[0]" ];
  };
};
