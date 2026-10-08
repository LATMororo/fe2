
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  Simple class for stress analysis. It has been intended to be used 
 *  with continuum finite element meshes. This implementation has been 
 *  adapted from the one given in JIVE's tutorial and course by Erik Jan 
 *  Lingen (Dynaflow Reaserch Group).
 *
 *  The main tasks of such model are:
 *
 *  Assemble internal force vector and stiffness matrix necessary for 
 *  linear and nonlinear analyses, and post-processing results after
 *  COMMIT action.
 *
 *  NOTE: columns of tables for post-processing operations are defined
 *        here; however, each material could define its own set of 
 *        columns. For instance, materials with plasticity might define
 *        a column for plastic strain field.
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date:   Fev 2025
 *  
 *  Changelog:
 *
 *  Modified: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br, luiztaumaturgo@gmail.com
 *  Date: Jul 2026
 *
 *  Compute average strains and stresses in the elements.
 *  
 */


#include <thread>


#include <jem/base/array/operators.h>
#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/Monitor.h>
#include <jem/base/System.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/util/PropertyException.h>


#include <jive/algebra/MatrixBuilder.h>
#include <jive/fem/ElementGroup.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/StateVector.h>
#include <jive/util/FuncUtils.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/XTable.h>
#include <jive/util/utilities.h>


#include "SolidModel.h"


JEM_DEFINE_CLASS ( SolidModel );


using namespace jem::literals;


using jem::System;

using jem::numeric::matmul;
using jem::numeric::MatmulChain;


typedef MatmulChain<double, 1> MChain1;
typedef MatmulChain<double, 3> MChain3;


//=======================================================================
//   class SolidModel::Worker_
//=======================================================================


class SolidModel::Worker_ : public jem::mt::WorkPool::Task
{
 public:

  explicit                  Worker_ 

    ( idx_t                   ieStart,
      idx_t                   ieEnd,
      const Ref<SolidModel>&  model );

  void                      setParams

    ( Ref<MatrixBuilder>&     mbld,
      const Vector&           force,
      const Vector&           disp );

  virtual void              run                   ()       override;  


 protected:

  virtual                  ~Worker_               ();


 private:

  jem::Monitor              monitor_;

  Ref<SolidModel>           model_;

  Ref<MatrixBuilder>        mbld_;

  Ref<jem::numeric::Function>
                            thickFunc_; 

  IdxVector                 dofTypes_;    

  Vector                    force_;
  Vector                    disp_;

  FEMUtils::ShapeGradsFunc  getShapeGrads_;

  idx_t                     ieStart_;
  idx_t                     ieEnd_;

  idx_t                     rank_;      
  idx_t                     strCount_;    

  idx_t                     nodeCount_;   
  idx_t                     ipCount_;     
  idx_t                     dofCount_;    

};


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


SolidModel::Worker_::Worker_

  ( idx_t                   ieStart,
    idx_t                   ieEnd,
    const Ref<SolidModel>&  model ) : 

    ieStart_ ( ieStart ), 
    ieEnd_   ( ieEnd   )

{
  model_     = model;

  rank_      = model_->rank_;
  strCount_  = model_->strCount_;
  nodeCount_ = model_->nodeCount_;
  ipCount_   = model_->ipCount_;
  dofCount_  = model_->dofCount_;

  dofTypes_.resize ( rank_ );
  dofTypes_  = model_->dofTypes_;

  thickFunc_ = model_->thickFunc_;

  getShapeGrads_ = FEMUtils::getShapeGradsFunc ( rank_ );

  mbld_      = nullptr;
  //force_     = nullptr;
  //disp_      = nullptr;
}


SolidModel::Worker_::~Worker_

  ()

{
  mbld_  = nullptr;
  //force_ = nullptr;
  //disp_  = nullptr;
}


//-----------------------------------------------------------------------
//   setParams
//-----------------------------------------------------------------------


void SolidModel::Worker_::setParams 

  ( Ref<MatrixBuilder>&  mbld,
    const Vector&        fint,
    const Vector&        disp ) 

{
  mbld_  = mbld;
  force_.ref ( fint );
  disp_ .ref ( disp );
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------


void SolidModel::Worker_::run

  ()

{
  using jem::ALL;
  using jem::Lock;
  using jem::Monitor;

  using jive::Cubix;


  const SolidModel& m = * model_;


  Cubix       grads      ( rank_, nodeCount_, ipCount_ );

  Matrix      stiff      ( strCount_, strCount_  );
  Matrix      coords     ( rank_,     nodeCount_ );
  Matrix      elemMat    ( dofCount_, dofCount_  );
  Matrix      b          ( strCount_, dofCount_  );
  Matrix      bt         = b.transpose ();
  Matrix      coordIp    ( rank_,     ipCount_   );

  Vector      elemForce  ( dofCount_  );
  Vector      elemDisp   ( dofCount_  );

  Vector      strain     ( strCount_  );
  Vector      stress     ( strCount_  );

  Vector      ipWeights  ( ipCount_   );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  MChain1     mc1;
  MChain3     mc3;

  idx_t       ipoint = ipCount_ * m.ielems_[ieStart_];

  for ( idx_t ie = ieStart_; ie < ieEnd_; ie++ )
  {
    {
      // Lock all this block so that 'global' (or 'shared') quantities
      // can be accessed by each thread without causing deadlocks and/or
      // race conditions.

      Lock<Monitor> lock ( monitor_ );

      // Get the global element index.

      const idx_t ielem = m.ielems_[ie];

      // Get the element nodes, coordinates and DOF indices.

      m.elems_.getElemNodes  ( inodes, ielem  );
      m.nodes_.getSomeCoords ( coords, inodes );
      m.dofs_->getDofIndices ( idofs , inodes, dofTypes_ );

      // Get the element shape function gradients and integration
      // point weights.

      m.shape_->getShapeGradients ( grads, ipWeights, coords );

      // Get the global coordinates of integration points.
    
      m.shape_->getGlobalIntegrationPoints ( coordIp, coords );

      // Get the current solution for the element.

      elemDisp = jem::select ( disp_, idofs );
    }

    // Assemble the element matrix and internal force vector.

    elemMat   = 0.0;
    elemForce = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Update the integration point weights. For 2D: multiply ipWeights 
      // with thickness.
      // NOTE: for 3D problems, thickness = 1.0.
      
      ipWeights[ip] *= thickFunc_->eval
                                ( coordIp(0, ip), coordIp(1, ip) );

      // Compute the B-matrix for this integration point.

      getShapeGrads_ ( b, grads(ALL,ALL,ip) );

      // Compute the strain: {strain} = [B] * {elemDisp}.

      matmul ( strain, b, elemDisp );

      // Compute the stress and constitutive matrix.

      m.material_->update ( stress, stiff, strain, ipoint++ );

      // Compute the element force vector and stiffness matrix

      elemForce += ipWeights[ip] * mc1.matmul ( bt, stress );
      elemMat   += ipWeights[ip] * mc3.matmul ( bt, stiff, b );
    }

    // Assemble: add element quantities to the global ones.
    // Lock again in order to avoid race conditions.

    {
      Lock<Monitor> lock ( monitor_ );
    
      // 1. Stiffness matrix.

      mbld_->addBlock ( idofs, idofs, elemMat );

      // 2. Internal force vector. 
    
      jem::select     ( force_, idofs ) += elemForce;
    }
  }
}


//=======================================================================
//   class SolidModel
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  SolidModel::TYPE_NAME      = "Solid";

const char*  SolidModel::SHAPE_PROP     = "shape";
const char*  SolidModel::MATERIAL_PROP  = "material";
const char*  SolidModel::THICKNESS_PROP = "thickness";
const char*  SolidModel::THREADS_PROP   = "threads";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


SolidModel::SolidModel

   ( const String&      name,
     const Properties&  conf,
     const Properties&  props,
     const Properties&  globdat ) : Super ( name )

{
  using jive::fem::ElementGroup;
  using jive::geom::IShapeFactory;
  using jive::util::XDofSpace;
  using jive::util::FuncUtils;

  using FEMUtils::STRAIN_COUNTS;
  using FEMUtils::DOF_TYPE_NAMES;


  // Get properties for this model.

  Properties myConf    = conf.makeProps ( myName_ );
  Properties myProps   = props.getProps ( myName_ );

  const String context = getContext     ();

  // Get elements and nodes from the globdat.

  const ElementGroup egroup = ElementGroup::get 
                                   ( myConf, myProps, globdat, context );

  elems_  = egroup.getElements ();
  nodes_  = elems_.getNodes    ();

  // Get element indices w.r.t. such ElementGroup 

  ielems_.ref ( egroup.getIndices() );

  // ... and make sure that TLS elements are sorted and unique.

  jem::sort              ( ielems_ );
  jive::util::makeUnique ( ielems_ );

  // Set extra members based on mesh settings.

  rank_     = nodes_.rank  ();
  numElem_  = ielems_.size ();

  strCount_ = STRAIN_COUNTS[ rank_ ];

  // Make sure that the number of spatial dimensions (the rank of the
  // mesh) is valid.

  if ( rank_ < 1_idx || rank_ > 3_idx )
  {
    throw jem::IllegalInputException (
      CLASS_NAME,
      String::format (
        "invalid node rank: %d (should be 1, 2 or 3)", rank_
      )
    );
  }

  // Create the internal shape functions

  shape_ = IShapeFactory::newInstance ( 
                jive::util::joinNames( myName_, SHAPE_PROP ), conf, props 
              );

  // ... and set more members.

  nodeCount_  = shape_->nodeCount   ();
  ipCount_    = shape_->ipointCount ();
  dofCount_   = rank_ * nodeCount_;

  // Make sure that the rank of the shape matches the rank of the mesh.

  if ( shape_->globalRank() != rank_ )
  {
    throw jem::IllegalInputException (
      CLASS_NAME,
      String::format (
        "shape has invalid rank: %d (should be %d)",
        shape_->globalRank(),
        rank_
      )
    );
  }

  elems_.checkSomeElements ( context, ielems_, nodeCount_ );

  // Create DOF space.
  
  Ref<XDofSpace> xdofs = XDofSpace::get ( nodes_.getData(), globdat );

  // 1. Map dof types.

  dofTypes_.resize ( rank_ );

  for (idx_t i = 0; i < rank_; i++) 
  {
    dofTypes_[i] = xdofs->addType ( DOF_TYPE_NAMES[i] );
  }

  // 2. Set dofs to the nodes.

  xdofs->addDofs ( elems_.getUniqueNodesOf( ielems_ ), dofTypes_ );

  // 3. Set DofSpace member object.

  dofs_ = xdofs;

  // Initialize the material.

  material_ = newMaterial   ( MATERIAL_PROP, myConf, myProps, globdat );

  JEM_ASSERT2 ( rank_ == material_->rank(), "material rank mismatch!" );

  // Set thickness member, which can be a 'constant' or a 'function'.

  // NOTE: for 2D cases, it could be the cross-section 'area'.
  
  String thickVars = "x,y"; 

  double thickness = 1.;

  // First, check if thickness is constant ('double'). Otherwise, check 
  // if thickness is a function ('string').

  if ( rank_ < 3_idx )
  {
    try 
    {
      myProps.find ( thickness, THICKNESS_PROP );

      thickFunc_ = FuncUtils::newFunc ( thickVars, String( thickness ) );

      myConf.set   ( THICKNESS_PROP, thickness );
    }
    catch ( const jem::util::PropertyException& ex ) 
    {
      String thickFunc = "";

      myProps.find ( thickFunc, THICKNESS_PROP );

      thickFunc_ = FuncUtils::newFunc ( thickVars, thickFunc );

      myConf.set   ( THICKNESS_PROP, thickFunc );
    }
  }
  else
  {
    thickFunc_ = FuncUtils::newFunc ( thickVars, String( thickness ) );
  }

  // Get shape function gradient.

  getShapeGrads = FEMUtils::getShapeGradsFunc ( rank_ );
}


SolidModel::~SolidModel

  ()

{
  if ( pool_ )
  {
    pool_->cancelAll   ();
    pool_->killWorkers ();
  }
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void SolidModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties myProps = props.findProps ( myName_ );

  // Let material handle its own 'configure'
  
  material_->configure   ( myProps, globdat );

  // ... also, allocate new material points that have state quantities,
  // for instance, history and plastic strain.

  material_->allocPoints ( numElem_ * ipCount_ );

  // Check if it is required to set a number of threads/workers.

  pool_ = nullptr;

  if ( myProps.contains( THREADS_PROP ) )
  {
    int maxThread = (int) std::thread::hardware_concurrency (); 

    if ( maxThread == 0 )
    {
      System::warn() << CLASS_NAME << ": total number of available" 
                                   << " threads is not well definedd!\n";
    }

    int nth;

    myProps.get ( nth, THREADS_PROP );

    if ( nth >= maxThread )
    {
      throw jem::IllegalInputException ( 
        CLASS_NAME, "setting more workers thant available threads!" 
      );
    }

    if ( nth == 1 )
    {
      nth += 1;

      System::warn() << CLASS_NAME << ": concurrency tasks are not"
                                   << " going to be performde with"
                                   << " 2 threads!\n";
    }

    if ( nth <= 1 )
    {
      System::warn() << CLASS_NAME << ": concurrency tasks are not"
                                   << " going to be performed!\n";
    }

    JEM_ASSERT2 ( nth <= maxThread, "outnumber available threads!" );

    if ( nth >= 2 )
    {
      pool_ = jem::newInstance<jem::mt::WorkPool> ( nth );

      jobs_   .reserve ( nth );
      workers_.reserve ( nth );

      const idx_t chunkSize = numElem_ / nth;
      const idx_t remainder = numElem_ % nth;

      idx_t currStart = 0_idx;

      idx_t currChunkSize, currEnd;

      for ( idx_t it = 0; it < nth; it++ )
      {
        currChunkSize = chunkSize + ( it < remainder ? 1_idx : 0_idx );
        currEnd       = currStart + currChunkSize;

        workers_.insert ( 
	      it, jem::newInstance<Worker_> ( currStart, currEnd, this ) 
	    );

        jobs_.insert ( 
          it, pool_->newJob ( workers_.getAs<Worker_>( it ) )
        );

        currStart = currEnd;
      }
    }
  }
} 


//-----------------------------------------------------------------------
//   getConfig 
//-----------------------------------------------------------------------


void SolidModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const 

{
  Properties myConf = conf.makeProps ( myName_ );

  // Let material handle its own 'getConfig'.

  material_->getConfig ( myConf, globdat );
  
  if ( ! pool_ )
  {
    myConf.set ( THREADS_PROP, 0 );
  }
  else
  {
    myConf.set ( THREADS_PROP, pool_->maxWorkers() );
  }
}


//-----------------------------------------------------------------------
//   takeAction 
//-----------------------------------------------------------------------


bool SolidModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat ) 

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;


  if ( action == Actions::INIT )
  {
//    jem::System::out() << CLASS_NAME << ": INIT action\n";
    //material_->allocPoints ( 1_idx );
    
    return true;
  }

  if ( action == Actions::ADVANCE )
  {
//    jem::System::out() << CLASS_NAME << ": ADVANCE . . . .\n";

    return true;
  }

  if ( action == Actions::GET_MATRIX0 )
  {
//    jem::System::out() << CLASS_NAME << ": GET_MATRIX0 . . . .\n";

    Ref<MatrixBuilder> mbld;
    Vector             fint;
    Vector             disp;

    // Get the current displacement field.

    StateVector::get ( disp, dofs_, globdat );

    // Get the (global) matrix builder and the (global) internal force
    // vector

    params.get ( mbld, ActionParams::MATRIX0    );
    params.get ( fint, ActionParams::INT_VECTOR );

    // ... add model contributions in such quantities.

    if ( ! pool_ )
    {
      assemble_  ( *mbld, fint, disp ); 
    }
    else
    {
      const idx_t njobs = jobs_.size ();

      for ( idx_t iw = 0; iw < njobs; iw++ )
      {
        workers_.getAs<Worker_>( iw )->setParams ( mbld, fint, disp );
      }

      for ( idx_t ij = 0; ij < njobs; ij++ )
      {
        jobs_.getAs<jem::mt::WorkPool::Job>( ij )->start ();
      }

      pool_->waitAll ();
    }

    return true;
  }

  if ( action == Actions::GET_TABLE )
  {
//    jem::System::out() << CLASS_NAME << ": GET_TABLE. . . .\n";

    const bool hasTable = getTable_ ( params, globdat );

//    jem::System::out() << CLASS_NAME << "    : " <<  hasTable << ". . . .\n";

    return hasTable; 
  }

  if ( action == Actions::COMMIT )
  {
//    jem::System::out() << CLASS_NAME << ": COMMIT . . . .\n";

    material_->commit ();

    return true;
  }

  if ( action == Actions::CANCEL )
  {
//    jem::System::out() << CLASS_NAME << ": CANCEL . . . .\n";

    material_->cancel ();

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   makeNew 
//-----------------------------------------------------------------------


Ref<Model> SolidModel::makeNew

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


void SolidModel::declare

  ()

{
  using jive::model::ModelFactory;

  ModelFactory::declare ( TYPE_NAME , & makeNew );
  ModelFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   assemble_
//-----------------------------------------------------------------------


void SolidModel::assemble_

  ( MatrixBuilder&  mbld,
    const Vector&   fint,
    const Vector&   disp ) const

{
  using jem::ALL;

  using jive::Cubix;
  

  Cubix       grads      ( rank_, nodeCount_, ipCount_ );

  Matrix      stiff      ( strCount_, strCount_  );
  Matrix      coords     ( rank_,     nodeCount_ );
  Matrix      elemMat    ( dofCount_, dofCount_  );
  Matrix      b          ( strCount_, dofCount_  );
  Matrix      bt         = b.transpose ();
  Matrix      coordIp    ( rank_,     ipCount_   );

  Vector      elemForce  ( dofCount_  );
  Vector      elemDisp   ( dofCount_  );

  Vector      strain     ( strCount_  );
  Vector      stress     ( strCount_  );

  Vector      ipWeights  ( ipCount_   );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  MChain1     mc1;
  MChain3     mc3;

  idx_t       ipoint = 0_idx;


  // Loop over all the elements.

  for (idx_t ie = 0; ie < numElem_; ie++) 
  {
    // Get the global element index.

    const idx_t ielem = ielems_[ie];

    // Get the element nodes, coordinates and DOF indices.

    elems_.getElemNodes  ( inodes, ielem  );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs , inodes, dofTypes_ );

    // Get the element shape function gradients and integration
    // point weights.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the global coordinates of integration points.
    
    shape_->getGlobalIntegrationPoints ( coordIp, coords );

    // Get the current solution for the element.

    elemDisp = jem::select ( disp, idofs );

    // Assemble the element matrix and internal force vector.

    elemMat   = 0.0;
    elemForce = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Update the integration point weights. For 2D: multiply ipWeights 
      // with thickness.
      // NOTE: for 3D problems, thickness = 1.0.
      
      ipWeights[ip] *= thickFunc_->eval
                                ( coordIp(0, ip), coordIp(1, ip) );

      // Compute the B-matrix for this integration point.

      getShapeGrads ( b, grads(ALL,ALL,ip) );

      // Compute the strain: {strain} = [B] * {elemDisp}.

      matmul ( strain, b, elemDisp );

      // Compute the stress and constitutive matrix.

      // FIX: passing ip to update function will cause an error for 
      //      other kind of materials rather than IsotropicMaterial.

      material_->update ( stress, stiff, strain, ipoint++ );

      // Compute the element force vector and stiffness matrix

      elemForce += ipWeights[ip] * mc1.matmul ( bt, stress );
      elemMat   += ipWeights[ip] * mc3.matmul ( bt, stiff, b );
    }

    // Assemble: add element quantities to the global ones.
    
    // 1. Stiffness matrix.

    mbld.addBlock ( idofs, idofs, elemMat );

    // 2. Internal force vector. 
    
    jem::select   ( fint, idofs ) += elemForce;
  }
}


//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool SolidModel::getTable_

  ( const Properties&  params,
    const Properties&  globdat ) const

{
  using jive::model::ActionParams;
  using jive::model::StateVector;


  Ref<XTable> table;

  Vector      weights;
  Vector      disp;

  String      name;

  // Get the table, the name of the table, and the table row weights
  // from the action parameters.

  params.get ( table  , ActionParams::TABLE         );
  params.get ( name   , ActionParams::TABLE_NAME    );
  params.get ( weights, ActionParams::TABLE_WEIGHTS );

  StateVector::get ( disp, dofs_, globdat );

  // Nodal strains.

  if ( name == "strain" && table->getRowItems() == nodes_.getData() )
  {
    setStrainTable_ ( *table, weights, disp );

    return true;
  }

  // Nodal stresses.

  if ( name == "stress" && table->getRowItems() == nodes_.getData() )
  {
    setStressTable_ ( *table, weights, disp );

    return true;
  }

  // Element strains.

  if ( name == "strain" && table->getRowItems() == elems_.getData() )
  {
    setElemStrainTable_ ( *table, disp );

    return true;
  }

  // Element stresses.

  if ( name == "stress" && table->getRowItems() == elems_.getData() )
  {
    setElemStressTable_ ( *table, disp );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   setStrainTable_
//-----------------------------------------------------------------------


void SolidModel::setStrainTable_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  disp ) const

{
  using jem::ALL;

  using jive::Cubix;


  Cubix       grads      ( rank_, nodeCount_, ipCount_ );

  Matrix      coords     ( rank_,     nodeCount_ );
  Matrix      b          ( strCount_, dofCount_  );

  Matrix      ndStrain   ( nodeCount_, strCount_ );

  Matrix      N          = shape_->getShapeFunctions ();

  Vector      ipStrain   ( strCount_  );
  Vector      elemDisp   ( dofCount_  );
  Vector      ndWeights  ( nodeCount_ );
  Vector      ipWeights  ( ipCount_   );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );
  IdxVector   jcols      ( strCount_  );


  // Add the columns for the strain components to the table.

  switch ( strCount_ )
  {
    case 1:
      
      jcols[0] = table.addColumn ( "strain_xx" );

      break;
    
    case 3:

      jcols[0] = table.addColumn ( "strain_xx" );
      jcols[1] = table.addColumn ( "strain_yy" );
      jcols[2] = table.addColumn ( "strain_xy" );

      break;
    
    case 6:

      jcols[0] = table.addColumn ( "strain_xx" );
      jcols[1] = table.addColumn ( "strain_yy" );
      jcols[2] = table.addColumn ( "strain_zz" );
      jcols[3] = table.addColumn ( "strain_xy" );
      jcols[4] = table.addColumn ( "strain_yz" );
      jcols[5] = table.addColumn ( "strain_xz" );

      break;

    default:

      throw jem::IllegalInputException ( 
        CLASS_NAME, 
        String::format ( 
        "unexpected number of strain components %d", strCount_ ) 
      );
  }

  // Loop over all elements. 

  for ( idx_t ie = 0; ie < numElem_; ie++ )
  {
    const idx_t ielem = ielems_[ie];

    // Get the element nodes, coordinates and DOF indices.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    
    // Get displacement field for this element.

    elemDisp = jem::select ( disp, idofs );

    // Get the element shape function gradients, and 
    // integration point weights.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Loop over all the integration points.

    ndStrain  = 0.0;
    ndWeights = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Compute strain field at a given integration point.

      ipStrain = 0.0;

      getShapeGrads ( b, grads( ALL,ALL,ip ) );

      matmul        ( ipStrain, b, elemDisp  );

      // Extrapolate the integration point strain to the nodes by
      // using shape functions.

      ndStrain  += matmul ( N( ALL, ip), ipStrain );

      ndWeights += N      ( ALL, ip );
    }

    // Update the table and the table weights.

    table.addBlock ( inodes, jcols, ndStrain );

    jem::select    ( weights, inodes ) += ndWeights;
  }
}


//-----------------------------------------------------------------------
//   setStressTable_
//-----------------------------------------------------------------------


void SolidModel::setStressTable_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  disp ) const

{
  using jem::ALL;

  using jive::Cubix;


  Cubix       grads      ( rank_, nodeCount_, ipCount_ );

  Matrix      coords     ( rank_,     nodeCount_ );
  Matrix      b          ( strCount_, dofCount_  );

  Matrix      ndStress   ( nodeCount_, strCount_ );

  Matrix      N          = shape_->getShapeFunctions ();

  Vector      ipStress   ( strCount_  );
  Vector      ipStrain   ( strCount_  );
  Vector      elemDisp   ( dofCount_  );
  Vector      ndWeights  ( nodeCount_ );
  Vector      ipWeights  ( ipCount_   );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );
  IdxVector   jcols      ( strCount_  );

  idx_t       ipoint = 0_idx;


  // Add the columns for the strain components to the table.

  switch ( strCount_ )
  {
    case 1:
      
      jcols[0] = table.addColumn ( "stress_xx" );

      break;
    
    case 3:

      jcols[0] = table.addColumn ( "stress_xx" );
      jcols[1] = table.addColumn ( "stress_yy" );
      jcols[2] = table.addColumn ( "stress_xy" );

      break;
    
    case 6:

      jcols[0] = table.addColumn ( "stress_xx" );
      jcols[1] = table.addColumn ( "stress_yy" );
      jcols[2] = table.addColumn ( "stress_zz" );
      jcols[3] = table.addColumn ( "stress_xy" );
      jcols[4] = table.addColumn ( "stress_yz" );
      jcols[5] = table.addColumn ( "stress_xz" );

      break;

    default:

      throw jem::IllegalInputException ( 
        CLASS_NAME, 
        String::format ( 
        "unexpected number of stress components %d", strCount_ ) 
      );
  }

  // Loop over all elements. 

  for ( idx_t ie = 0; ie < numElem_; ie++ )
  {
    const idx_t ielem = ielems_[ie];

    // Get the element nodes, coordinates and DOF indices.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    
    // Get displacement field for this element.

    elemDisp = jem::select ( disp, idofs );

    // Get the element shape function gradients, and 
    // integration point weights.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Loop over all the integration points.

    ndStress  = 0.0;
    ndWeights = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Compute strain and stress fields at a given integration point.

      ipStress = 0.0;
      ipStrain = 0.0;

      getShapeGrads     ( b, grads( ALL,ALL,ip ) );

      matmul            ( ipStrain, b, elemDisp  );

      material_->update ( ipStress, ipStrain, ipoint++ );

      // Extrapolate the integration point stress to the nodes by
      // using shape functions.

      ndStress  += matmul ( N( ALL, ip), ipStress );

      ndWeights += N      ( ALL, ip );
    }

    // Update the table and the table weights.

    table.addBlock ( inodes, jcols, ndStress );

    jem::select    ( weights, inodes ) += ndWeights;
  }
}


//-----------------------------------------------------------------------
//   setElemStrainTable_
//-----------------------------------------------------------------------


void SolidModel::setElemStrainTable_

  ( XTable&        table,
    const Vector&  disp ) const

{
  using jem::ALL;

  using jive::Cubix;


  Cubix       grads      ( rank_, nodeCount_, ipCount_ );

  Matrix      coords     ( rank_,     nodeCount_ );
  Matrix      b          ( strCount_, dofCount_  );

  Vector      elStrain   ( strCount_  );
  Vector      elemDisp   ( dofCount_  );
  Vector      ipWeights  ( ipCount_   );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );
  IdxVector   jcols      ( strCount_  );


  // Add the columns for the strain components to the table.

  switch ( strCount_ )
  {
    case 1:
      
      jcols[0] = table.addColumn ( "strain_xx" );

      break;
    
    case 3:

      jcols[0] = table.addColumn ( "strain_xx" );
      jcols[1] = table.addColumn ( "strain_yy" );
      jcols[2] = table.addColumn ( "strain_xy" );

      break;
    
    case 6:

      jcols[0] = table.addColumn ( "strain_xx" );
      jcols[1] = table.addColumn ( "strain_yy" );
      jcols[2] = table.addColumn ( "strain_zz" );
      jcols[3] = table.addColumn ( "strain_xy" );
      jcols[4] = table.addColumn ( "strain_yz" );
      jcols[5] = table.addColumn ( "strain_xz" );

      break;

    default:

      throw jem::IllegalInputException ( 
        CLASS_NAME, 
        String::format ( 
        "unexpected number of strain components %d", strCount_ ) 
      );
  }

  // Loop over all elements. 

  for ( idx_t ie = 0; ie < numElem_; ie++ )
  {
    const idx_t ielem = ielems_[ie];

    // Get the element nodes, coordinates and DOF indices.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    
    // Get displacement field for this element.

    elemDisp = jem::select ( disp, idofs );

    // Get the element shape function gradients, and 
    // integration point weights.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Loop over all the integration points.

    elStrain = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Compute strain field at a given integration point.

      getShapeGrads ( b, grads( ALL,ALL,ip ) );

      // Sum all the strain fields of such element.

      elStrain += matmul ( b, elemDisp  );
    }

    // Average the strain field.
    
    elStrain /= static_cast<double> ( ipCount_ );

    // Update the table.

    table.setRowValues ( ielem, jcols, elStrain );
  }
}


//-----------------------------------------------------------------------
//   setElemStressTable_
//-----------------------------------------------------------------------


void SolidModel::setElemStressTable_

  ( XTable&        table,
    const Vector&  disp ) const

{
  using jem::ALL;

  using jive::Cubix;


  Cubix       grads      ( rank_, nodeCount_, ipCount_ );

  Matrix      coords     ( rank_,     nodeCount_ );
  Matrix      b          ( strCount_, dofCount_  );

  Vector      ipStress   ( strCount_  );
  Vector      ipStrain   ( strCount_  );
  Vector      elStress   ( strCount_  );
  Vector      elemDisp   ( dofCount_  );
  Vector      ipWeights  ( ipCount_   );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );
  IdxVector   jcols      ( strCount_  );

  idx_t       ipoint = 0_idx;


  // Add the columns for the strain components to the table.

  switch ( strCount_ )
  {
    case 1:
      
      jcols[0] = table.addColumn ( "stress_xx" );

      break;
    
    case 3:

      jcols[0] = table.addColumn ( "stress_xx" );
      jcols[1] = table.addColumn ( "stress_yy" );
      jcols[2] = table.addColumn ( "stress_xy" );

      break;
    
    case 6:

      jcols[0] = table.addColumn ( "stress_xx" );
      jcols[1] = table.addColumn ( "stress_yy" );
      jcols[2] = table.addColumn ( "stress_zz" );
      jcols[3] = table.addColumn ( "stress_xy" );
      jcols[4] = table.addColumn ( "stress_yz" );
      jcols[5] = table.addColumn ( "stress_xz" );

      break;

    default:

      throw jem::IllegalInputException ( 
        CLASS_NAME, 
        String::format ( 
        "unexpected number of stress components %d", strCount_ ) 
      );
  }

  // Loop over all elements. 

  for ( idx_t ie = 0; ie < numElem_; ie++ )
  {
    const idx_t ielem = ielems_[ie];

    // Get the element nodes, coordinates and DOF indices.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    
    // Get displacement field for this element.

    elemDisp = jem::select ( disp, idofs );

    // Get the element shape function gradients, and 
    // integration point weights.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Loop over all the integration points.

    elStress = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Compute strain and stress fields at a given integration point.

      ipStress = 0.0;
      ipStrain = 0.0;

      getShapeGrads     ( b, grads( ALL,ALL,ip ) );

      matmul            ( ipStrain, b, elemDisp  );

      material_->update ( ipStress, ipStrain, ipoint++ );

      // Update the element stress vector.

      elStress += ipStress;
    }

    // Average the stress field.

    elStress /= static_cast<double> ( ipCount_ );

    // Update the table.

    table.setRowValues ( ielem, jcols, elStress );
  }
} 
