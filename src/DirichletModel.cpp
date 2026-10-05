
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This class implements a model for setting Dirichlet boundary 
 *  conditions (BCs) for certain group of nodes. 
 *  
 *  It inherits and follows what has been implemented F. P. van der Meer
 *  his collegues at TU Delft. Clean-up, refactored member functions in 
 *  order to meet new Jem/Jive (version 3.0) data types (e.g., `idx_t')
 *  have been accomplished.
 *
 *  However, the state increment is accomplished by means of a function,
 *  which can also be defined by the user. In this regard, such function
 *  can be written in term of coordinates (x, y, and z - they depend on
 *  the dimensionality of the problem), TIME_STEP (i) and TIME (t);
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date:   April 2026
 */


#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/util/PropertyException.h>


#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/Globdat.h>
#include <jive/model/ModelFactory.h>
#include <jive/util/FuncUtils.h>
#include <jive/util/Constraints.h>
#include <jive/util/DofSpace.h>


#include "DirichletModel.h"


JEM_DEFINE_CLASS ( DirichletModel );


using namespace jem::literals;


//=======================================================================
//   class DirichletModel 
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  DirichletModel::TYPE_NAME        = "Dirichlet";

const char*  DirichletModel::NODE_GROUPS_PROP = "nodeGroups";
const char*  DirichletModel::DOFS_PROP        = "dofs";
const char*  DirichletModel::STATE_FUNCS_PROP = "stateFuncs";
const char*  DirichletModel::LOADED_PROP      = "loaded";
const char*  DirichletModel::STATE_INCR_PROP  = "stateIncr";
const char*  DirichletModel::INIT_STATE_PROP  = "initState";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


DirichletModel::DirichletModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  stateFuncs_ = nullptr;
  ngroups_    = 0_idx;
}


DirichletModel::~DirichletModel

  ()

{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void DirichletModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jem::Array;

  using jive::util::FuncUtils;


  Properties  myProps = props.findProps ( myName_ );

  // Get node groups on which BCs will be set.

  myProps.get ( nodeGroups_, NODE_GROUPS_PROP );

  ngroups_ = nodeGroups_.size ();

  if ( !( ngroups_ > 0_idx ) )
  {
    throw jem::IllegalInputException ( 
      CLASS_NAME, "an empty list ``nodeGroups'!"
    );
  }

  // Get DOFs associated with the node groups. 

  myProps.get ( dofTypes_, DOFS_PROP );

  // The length of DOFs must be equal to the length of node groups. 

  if ( dofTypes_.size() != ngroups_ )
  {
    throw jem::IllegalInputException ( 
      CLASS_NAME, 
      "`dofTypes' must have the same length as `nodeGroups'!" 
    );
  }

  // Define function arguments and a temp. function array.

  const idx_t rank = NodeSet::get( globdat, getContext() ).rank ();
  String      args = "i, t";

  if      ( rank == 1_idx )
  {
    args = args + ", x";
  }
  else if ( rank == 2_idx )
  {
    args = args + ", x, y";
  }
  else 
  {
    args = args + ", x, y, z";
  }
  
  // Get state functions associated with the node groups. The length of
  // functions must be equal to the length of the node groups. Otherwise,
  // check whether 'loaded' parameter is set.

  Array< Ref<Function> > funcs ( ngroups_ );
  funcs = nullptr;

  if (  myProps.contains( STATE_FUNCS_PROP ) && 
       !myProps.contains( LOADED_PROP      ) )
  {
    FuncUtils::configFuncs ( funcs, args,
                             STATE_FUNCS_PROP,
                             myProps, globdat );

    if ( funcs.size() > 0_idx && funcs.size() != ngroups_ )
    {
      throw jem::IllegalInputException ( 
        CLASS_NAME,
        "`stateFuncs' must have the same length as `nodeGroups'!" 
      );
    }
  }
  else
  {
    // Find which group will be 'loaded'.
    // NOTE: `loaded' is an indice in `nodeGroups'.

    idx_t loaded = 0_idx;

    if ( ! myProps.find( loaded, LOADED_PROP, 0_idx, ngroups_-1_idx ) )
    {
      jem::System::warn() << CLASS_NAME 
                          << ": 'loading' the first group `" 
                          << nodeGroups_[loaded] << "' with dof type `" 
                          << dofTypes_[loaded] << "'!\n";
    }

    // For a loaded setting, first, check if there is an initial state
    // and/or an state increment.

    Array< Ref<Function> > func1 ( 1 );
    func1[0] = nullptr;

    double stateIncr = 0.0;
    double initState = 0.0;

    // If there is an initial state

    myProps.find ( stateIncr, STATE_INCR_PROP );

    // ... cancel the first increment.

    if ( myProps.find( initState, INIT_STATE_PROP ) )
    {
      initState -= stateIncr;
    }

    // Define a function to handle such situation based on 
    // TIME_STEP only, i.e., `i' function argument.

    func1[0] = FuncUtils::newFunc ( 
      args, 
      String( String( initState ) + " + " + String( stateIncr ) + "*i" ), 
      globdat 
    ); 

    // Second, for a loaded setting, check if there is a single function 
    // for such node group.

    if ( myProps.contains( STATE_FUNCS_PROP ) )
    {
      // First, try to get properties from a string, otherwise, try to 
      // get from a 'list' of functions (which must have only one 
      // element), e.g., `stateFuncs = ['0.001*i']'.

      try
      {
        FuncUtils::configFunc ( func1[0], args,
                                STATE_FUNCS_PROP,
                                myProps, globdat );
      }
      catch ( const jem::util::PropertyException& ex ) 
      {
        FuncUtils::configFuncs ( func1, args,
                                 STATE_FUNCS_PROP,
                                 myProps, globdat );
      }
    }

    // Update the function list.

    const String zeroFunc = String ( "0.0" );

    for ( idx_t ig = 0; ig < ngroups_; ig++ )
    {
      if ( ig == loaded )
      {
        funcs[ig] = func1[0];

        continue;
      }

      // Constant functions that always return zero.

      funcs[ig] = FuncUtils::newFunc ( args, zeroFunc, globdat );
    }
  }

  // Everything is ok, then, commit.

  stateFuncs_.swap ( funcs );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DirichletModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  using jive::util::FuncUtils;


  Properties myConf = conf.makeProps ( myName_ );

  myConf.set ( NODE_GROUPS_PROP, nodeGroups_ );
  myConf.set ( DOFS_PROP       , dofTypes_   );

  FuncUtils::getConfig ( myConf, stateFuncs_, STATE_FUNCS_PROP );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool DirichletModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;


  // Initialization.

  if ( action == Actions::INIT )
  {
    init_ ( globdat );

    return true;
  }

  // Advance to next time step.

  if ( action == Actions::ADVANCE )
  {
    globdat.set ( "var.accepted", true );

    //advance_    ( globdat );

    return true;
  }

  // Apply constraints.

  if ( action == Actions::GET_CONSTRAINTS )
  {
    applyConstraints_ ( params, globdat );

    return true;
  }
  
  // Proceed to next time step.

  if ( action == Actions::COMMIT )
  {
    //commit_ ( params, globdat );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   makeNew 
//-----------------------------------------------------------------------


Ref<Model> DirichletModel::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  return jem::newInstance<Self> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declare 
//-----------------------------------------------------------------------


void DirichletModel::declare

  ()

{
  using jive::model::ModelFactory;


  ModelFactory::declare ( TYPE_NAME , & makeNew );
  ModelFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   init_ 
//-----------------------------------------------------------------------


void DirichletModel::init_ 

  ( const Properties&  globdat )

{
  // Get nodes, then DOFs of nodes, and constraints of DOFs from globdat.

  const String context = getContext ();

  nodes_ = NodeSet::get     ( globdat, context );
  dofs_  = DofSpace::get    ( nodes_.getData(), globdat, context );
  cons_  = Constraints::get ( dofs_, globdat );    
}


//-----------------------------------------------------------------------
//   advance_
//-----------------------------------------------------------------------


void DirichletModel::advance_

  ( const Properties&  globdat )

{}


//-----------------------------------------------------------------------
//   applyConstraints_
//-----------------------------------------------------------------------


void DirichletModel::applyConstraints_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::IdxVector;
  using jive::Vector;

  using jive::fem::NodeGroup;
  using jive::model::Globdat;


  const String           ctxt = getContext  ();
  const idx_t            rank = nodes_.rank ();

  Assignable<NodeGroup>  group;

  Vector                 icoords ( rank ); 

  IdxVector              inodes;

  double                 args[5];
  double                 itime;
  double                 val;

  idx_t                  istep = -1_idx;
  idx_t                  inode = -1_idx;
  idx_t                  nn    = -1_idx;
  idx_t                  itype = -1_idx;
  idx_t                  idof  = -1_idx;

  // Get current TIME_STEP and/or TIME from globdat.

  globdat.find ( istep, Globdat::TIME_STEP );
  globdat.find ( itime, Globdat::TIME      );

  args[0] = (double) istep;
  args[1] = itime;

  // Loop over all the node groups.

  for ( idx_t ig = 0; ig < ngroups_; ig++ )
  {
    group  = NodeGroup::get ( nodeGroups_[ig], nodes_, globdat, ctxt );

    nn     = group.size     ();

    inodes . resize ( nn );
    inodes = group.getIndices ();

    itype  = dofs_->findType  ( dofTypes_[ig] );

    JEM_PRECHECK2 ( itype >= 0_idx, "invalid dof type!" );

    // Apply constraint.

    for ( idx_t in = 0; in < nn; in++ )
    {
      // Get node index, its attached dof index and coordinates.

      inode = inodes[in];

      idof  = dofs_->getDofIndex ( inode,   itype );
      nodes_.getNodeCoords       ( icoords, inode );

      // Add nodal coordinates into 'args'.

      for ( idx_t ir = 0; ir < rank; ir++ )
      {
        args[ir + 2] = icoords[ir];
      }

      // Add constraint.

      val = stateFuncs_[ig]->getValue ( args );

      cons_->addConstraint ( idof, val );
    }
  }

  // Compress for more efficient storage.

  cons_->compress ();
}


//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------


void DirichletModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{}


//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------


void DirichletModel::setStepSize_

  ( const Properties&  params )

{}
