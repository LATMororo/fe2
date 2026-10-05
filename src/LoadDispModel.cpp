
/*
 *  Model that prints load-displacement data for nodegroup to file:
 *  The sum of the nodal loads and the average of the nodal displacements
 *
 *  Adapted from: MonitorModel, 
 *  FPM, 29-1-2008
 * 
 *  Data is computed in GET_MATRIX0 and GET_INT_VECTOR
 *         and optionally written to file in COMMIT
 *         (rather using the SampleModule to write to file recommended)
 *
 *  Changelog:
 *
 *  Modified: Luiz Ant. T. Mororo, l.a.taumaturgomororo@tudelft.nl
 *  Date:     July, 2019
 *
 *  The class has been changed in order to enable parallel computing.
 *
 *  Modified: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date: Fev 2025
 *
 *  Clean-up code, and refactored member functions for better
 *  performance, and to meet new Jem/Jive (version 3.0) data types, for 
 *  instance, `idx_t'. 
 */


#include <jem/base/array/operators.h>
#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/System.h>
#include <jem/mp/Buffer.h>
#include <jem/mp/Context.h>
#include <jem/util/ArrayBuffer.h>
#include <jem/util/Properties.h>


#include <jive/model/Actions.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/StateVector.h>
#include <jive/mp/Globdat.h>
#include <jive/util/utilities.h>


#include "LoadDispModel.h"


JEM_DEFINE_CLASS ( LoadDispModel );


using namespace jem::literals;



//=======================================================================
//   class LoadDispModel
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  LoadDispModel::TYPE_NAME      = "LoadDisp";

const char*  LoadDispModel::NODES_PROP     = "nodes";
const char*  LoadDispModel::GROUP_PROP     = "group";
const char*  LoadDispModel::FILE_NAME_PROP = "file";
const char*  LoadDispModel::TYPES_PROP     = "types";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


LoadDispModel::LoadDispModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  const String context = getContext ();

  nodes_     = NodeSet::get     ( globdat, context ); 
  dofs_      = DofSpace::get    ( nodes_.getData (), globdat, context );
  nn_        = -1_idx;
  nnOld_     = -1_idx;
  rank_      = nodes_.rank      ();
  typeCount_ = dofs_->typeCount ();
  groupName_ = "";
  vex_       = nullptr;
  mpx_       = nullptr;
}


LoadDispModel::~LoadDispModel 

  ()

{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void LoadDispModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jem::util::ArrayBuffer;

  using jive::util::Globdat;


  Properties    myProps = props.findProps ( myName_ );

  const String  context = getContext ();

  IdxVector     nodeIDs;

  Properties    myVars = Globdat::getVariables ( myName_, globdat );

  if ( myProps.find( nodeIDs, NODES_PROP ) )
  {
    // Get node and dof numbers specified with 'nodes = IdxVector'.

    ArrayBuffer<idx_t> inbuff;
    
    nn_ = nodeIDs.size ();

    for ( idx_t i = 0; i < nn_; i++ )
    {
      idx_t k = nodes_.findNode ( nodeIDs[i] );

      if ( k < 0_idx ) continue;

      inbuff.pushBack ( k );
    }

    if ( inbuff.size() > 0 )
    {
      inodes_.ref    ( inbuff.toArray() );

      idofs_ .resize ( inodes_.size(), types_.size() );

      for ( idx_t i = 0; i < inodes_.size(); i++ )
      {
        for ( idx_t j = 0; j < types_.size(); j++ )
        {
          idofs_(i,j) = dofs_->getDofIndex ( inodes_[i], types_[j] );
        }
      }
    }
  }

  if ( myProps.find( groupName_, GROUP_PROP ) )
  {
    // Get node and dof numbers specified with 'group = String'.
    
    ngroup_ = NodeGroup::get ( groupName_, nodes_, globdat, context );
  }
  
  if ( myProps.find( typeNames_, TYPES_PROP ) )
  {
    types_.resize ( typeNames_.size() );
  }
  else
  {
    types_    .resize ( rank_ );
    typeNames_.resize ( rank_ );

    for ( idx_t i = 0; i < rank_; i++ )
    {
      typeNames_[i] = dofs_->getTypeName ( i );
    }
  }

  for ( idx_t i = 0; i < typeNames_.size(); i++ )
  {
    types_[i] = dofs_->getTypeIndex ( typeNames_[i] );
  }

  // Get the number of processores and set isParallel_.

  isParallel_ = 
      ( jive::mp::Globdat::procCount( globdat ) > 1 ? true : false );

  // Get the VectorExchanger and MPContext if they are required.

  if ( isParallel_ )
  {
    if ( ! vex_ )
    {
      vex_ = VectorExchanger::find ( dofs_, globdat );
    }

    if ( ! vex_ )
    {
      throw jem::IllegalInputException 
                              ( CLASS_NAME, "VectorExchanger is NULL!" );
    }

    if ( ! mpx_ )
    {
      //mpx_ = jive::mp::Globdat::getMPContext ( globdat );
      mpx_ = vex_->getMPContext ();
    }

    if ( ! mpx_                                       || 
         mpx_->size() <= 1                            ||
         jive::mp::Globdat::procCount( globdat ) <= 1  )

    {
      throw jem::IllegalInputException ( 
                CLASS_NAME, String::format (
                "MPContext is not found! %d", 
                mpx_->size() 
              ) 
            );
    }
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void LoadDispModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  if ( inodes_.size() )
  {
    myConf.set ( NODES_PROP, inodes_    );
  }

  myConf.set   ( GROUP_PROP, groupName_ );

  myConf.set   ( TYPES_PROP, typeNames_ );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool LoadDispModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jem::ALL;
  using jem::select;
  using jem::sum;

  using jive::util::Globdat;
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;


  if ( nn_ < 0_idx && ! nodes_ )
  {
    return false;
  }

  Vector  load ( types_.size() );
  Vector  disp ( types_.size() );

  if ( action == Actions::GET_MATRIX0   ||
       action == Actions::GET_INT_VECTOR )
  {
    Properties  myVars = Globdat::getVariables ( myName_, globdat );

    // Get global data.

    Vector      fint;
    Vector      state;

    params.get       ( fint , ActionParams::INT_VECTOR );
    StateVector::get ( state, dofs_, globdat           );

    if ( groupName_.size() > 0_idx )
    {
      updateIDofs_ ();
    }

    // Exchange data if it is a parallel computing.

    if ( isParallel_ )
    {
      gatherData_ ( disp, load, state, fint, globdat );

      // Store in globdat.
      // NOTE: it is only meaningful for rank 0 (root).

      myVars.set ( "load" , load );
      myVars.set ( "disp" , disp );

      return true;
    }

    // Compute cumulative load and average displacement.

    for ( idx_t i = 0; i < types_.size(); i++ )
    {
      load[i] = sum ( select( fint , idofs_(ALL,i) ) );
      disp[i] = sum ( select( state, idofs_(ALL,i) ) ) / (double)nn_;
    }

    // Store in globdat.

    myVars.set ( "load" , load );
    myVars.set ( "disp" , disp );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   makeNew 
//-----------------------------------------------------------------------


Ref<Model> LoadDispModel::makeNew

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


void LoadDispModel::declare

  ()

{
  using jive::model::ModelFactory;

  ModelFactory::declare ( TYPE_NAME , & makeNew );
  ModelFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   updateIDofs_
//-----------------------------------------------------------------------


void LoadDispModel::updateIDofs_

  ()

{
  using jem::util::ArrayBuffer;

  using jive::mp::BorderSet;


  nn_ = ngroup_.size ();

  inodes_.resize ( nn_ );
  idofs_ .resize ( nn_ , types_.size() );

  inodes_ = ngroup_.getIndices ();

  for ( idx_t i = 0; i < nn_; i++ )
  {
    for ( idx_t j = 0; j < types_.size(); j++ )
    {
      idofs_(i,j) = dofs_->findDofIndex ( inodes_[i], types_[j] );
    }
  }

  // Update idofborders_ variable.

  if ( isParallel_ && nn_ != nnOld_ )
  {
    ArrayBuffer<idx_t> idofbuff;

    BorderSet          recvb   ( vex_->getRecvBorders() );

    IdxVector          iitems  ( recvb.maxBorderSize () );

    int                rbCount = recvb.size ();
    int                iitemCount;

    for ( idx_t irb = 0; irb < rbCount; irb++ )
    {
      iitemCount = recvb.getBorderItems ( iitems, irb );

      for ( idx_t i = 0; i < iitemCount; i++ )
      {
        if ( jem::testany( iitems[i] == inodes_ ) )
	{
	  for ( idx_t j = 0; j < types_.size(); j++ )
	  {
	    idofbuff.pushBack 
	              ( dofs_->findDofIndex( iitems[i], types_[j] ) );
	  }
	}
      }
    }

    idofborders_.ref ( jive::util::makeUnique( idofbuff.toArray() ) );

    nnOld_ = nn_;
  }
}


//-----------------------------------------------------------------------
//   gatherData_ 
//-----------------------------------------------------------------------


void LoadDispModel::gatherData_

  ( const Vector&      disp,
    const Vector&      load,
    const Vector&      state,
    const Vector&      fint,
    const Properties&  globdat ) const

{
  using jem::ALL;

  using jem::mp::SendBuffer;
  using jem::mp::RecvBuffer;

  using jive::util::Globdat;


  JEM_ASSERT ( state.size() == fint.size() );


  const int root = 0;

  Vector    lfint   ( fint  .size() );
  Vector    lstate  ( fint  .size() );

  Vector    weights ( fint  .size() );
  Vector    ldisp   ( types_.size() );
  Vector    lload   ( types_.size() );

  // Scatter weights.

  weights = 1.;

  vex_->scatter ( weights );

  // Compute 'averaged' vectors.

  lstate = state / weights;

  if ( inodes_.size() > 0 )
  {
    if ( vex_->hasOverlap() )
    {
      lfint = fint;

      if ( idofborders_.size() > 0 )
      {
        lfint[ idofborders_ ] = 0.;
      }
    }
    else
    {
      lfint.ref ( fint );
    }
  }

  // Compute cumulative load and displacement. Skip processors
  // without inodes_.

  if ( inodes_.size () > 0 )
  {
    for ( idx_t i = 0; i < types_.size(); i++ )
    {
      lload[i] = jem::sum ( select( lfint , idofs_(ALL,i) ) );
      ldisp[i] = jem::sum ( select( lstate, idofs_(ALL,i) ) );
    }
  }
  else
  {
    lload = 0.;
    ldisp = 0.;
  }

  // Reduce operation.
  // NOTE: the outcome of this operation is stored into rank 0 (root).
  // NOTE: SampleModule only takes rank 0 into account.

  disp = 0.;
  load = 0.;

  mpx_->reduce ( RecvBuffer( disp .addr(), disp .size() ), 
                 SendBuffer( ldisp.addr(), ldisp.size() ),
		 root,
		 jem::mp::SUM );

  mpx_->reduce ( RecvBuffer( load .addr(), load .size() ), 
                 SendBuffer( lload.addr(), lload.size() ),
		 root,
		 jem::mp::SUM );

  // Compute averaged displacement.

  if ( mpx_->myRank() == 0 )
  {
    int nn;

    if ( groupName_.size() > 0 )
    {
      Properties  myVars = Globdat::getVariables ( groupName_, globdat );

      myVars.get  ( nn, groupName_ );
    }
    else
    {
      nn = nn_;
    }

    disp = disp / (double)nn;
  }
}
