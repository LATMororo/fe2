
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This solver applies PBC on a RVE and computes its homogenized state
 *  quantities [1-2], i.e., the macro constitutive law, namely, macro 
 *  stress field and (tangent) stiffness matrix. PBC are evaluated in 
 *  two modes:
 *
 *  1) `standalone_ == true', the user can run a simple RVE without the 
 *                          whole FE^2 framework, and a macro strain 
 *                          field has to be provided.
 *
 *  2) `standalone_ == false', it is coupled in a FE^2 framework, in 
 *                           which a macro strain is given and attached
 *                           to an integration point on a upper level.
 *
 *  It follows what has been used implemented by TU Delft's EM group 
 *  [3-5], and part of this implementation has been adapted from 
 *  Rocha [5].
 *
 *  [1] R. J. M. Smit, W. A. M. Brekelmans, and H. E. H. Meijer. 
 *      Prediction of the mechanical behavior of nonlinear heterogeneous 
 *      systems by multi-level finite element modeling. Computer Methods 
 *      in Applied Mechanics and Engineering 155 (1) (1998), pp. 181–192.
 *
 *  [2] MultiFreedom Constraints I. 
 *      https://quickfem.com/wp-content/uploads/IFEM.Ch08.pdf. 
 *      Accessed: 02-01-2023.
 *
 *  [3] V. P. Nguyen, O. Lloberas-Valls, M. Stroeven, and L. J. Sluys. 
 *      Computational homogenization for multiscale crack modeling.
 *      Implementational and computational aspects. Int. J. Numer. Meth. 
 *      Engng, 89:192-226, 2011.
 *
 *  [4] V. P. Nguyen, Multiscale failure modelling of quasi-brittle 
 *      materials. PhD dissertation. Delft University of Technology, 
 *      2011.
 *
 *  [5] I. B. C. M. Rocha. Numerical and Experimental Investigation of 
 *      Hygrothermal Aging in Laminated Composites. Ph.D. dissertation. 
 *      Delft University of Technology, 2019.
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date:   Sep 2025
 *
 *  Changelog:
 *
 *  Functions have been defined in order to control macro strain field
 *  scale in a standalone mode. It allows the user to define a 
 *  function(s) that describes the increment/behaviour of a macroscopic
 *  strain component of a given RVE.
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br, luiztaumaturgo@gmail.com
 *  Date:   May 2026
 */


#include <jem/base/array/operators.h>
#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/io/FileWriter.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>


#include <jive/algebra/AbstractMatrix.h>
#include <jive/algebra/ConstrainedMatrix.h>
#include <jive/algebra/StdVectorSpace.h>
#include <jive/app/ModuleFactory.h>
#include <jive/fem/FEMatrixBuilder.h>
#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/Model.h>
#include <jive/solver/Solver.h>
#include <jive/solver/SolverParams.h>
#include <jive/util/FuncUtils.h>
#include <jive/util/Globdat.h>
#include <jive/util/utilities.h>


#include "MicroPBCSolverModule.h"
#include "PBCGroupInputModule.h"
#include "FEMUtils.h"


JEM_DEFINE_CLASS ( MicroPBCSolverModule );


using namespace jem::literals;


using jem::System;


using jive::algebra::ConstrainedMatrix;
using jive::fem::FEMatrixBuilder;
using jive::model::Model;
using jive::solver::Solver;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::DofSpace;
using jive::util::FuncUtils;


//=======================================================================
//   class MicroPBCSolverModule::RunMicroData_ 
//=======================================================================


class MicroPBCSolverModule::RunMicroData_ : public jem::Object
{
 public:

  explicit                  RunMicroData_

    ( const String&           ctxt );

  void                      init

    ( const Properties&       globdat );
  
  void                      initMicroSolver

    ( const String&           name, 
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  void                      setConstraints        

    ( const Vector&           strain )               const;

  void                      calcState             

    ( MicroState&             mstate, 
      const Properties&       globdat )              const;


 protected:

  virtual                  ~RunMicroData_         ();

 
 private:

  void                      calcInvVol_           ();
  void                      calcBndHMatrix_       ();
  void                      getCtrlAndRestDofs_   ();

  void                      setPeriodicCons_      

    ( Constraints&            c )                     const;

  void                      setCtrlNodeCons_       

    ( const Vector&           strain )                const;


 public:

  Ref<Solver>               msolver;
  Ref<FEMatrixBuilder>      mBuilder;
  Ref<ConstrainedMatrix>    consMat;

  Assignable<ElementSet>    elems;
  Assignable<NodeSet>       nodes;

  Ref<Model>                model;

  Ref<DofSpace>             dofs;
  Ref<Constraints>          cons;
  Ref<Constraints>          cons2;

  Matrix                    H; // boundary matrix

  Vector                    dxRVE;

  IdxVector                 ctrlNodes;
  IdxVector                 ctrlDofs; 
  IdxVector                 restDofs;
  IdxVector                 bndNodes[6];
  IdxVector                 dofTypes;  // dof type indices

  const String              context;

  idx_t                     rank;
  idx_t                     strCount;

  double                    v0Inv;

};


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


MicroPBCSolverModule::RunMicroData_::RunMicroData_ 

  ( const String&  ctxt ) : context ( ctxt )

{}


MicroPBCSolverModule::RunMicroData_::~RunMicroData_ 

  ()

{}


//-----------------------------------------------------------------------
//   MicroPBCSolverModule::RunMicroData_::init
//-----------------------------------------------------------------------


void MicroPBCSolverModule::RunMicroData_::init

  ( const Properties&  globdat )

{
  using jive::algebra::AbstractMatrix;
  using jive::fem::NodeGroup;

  using FEMUtils::STRAIN_COUNTS;
  using FEMUtils::DOF_TYPE_NAMES;


  // From globdat, get the primary Model, ElementSet (and its 
  // corresponding NodeSet), DofSpace, and Constraints.

  model = Model       :: get ( globdat, context );
  elems = ElementSet  :: get ( globdat, context );
  dofs  = DofSpace    :: get ( globdat, context );
  cons  = Constraints :: get ( dofs,    globdat );

  nodes = elems.getNodes     ();

  // Get node rank and the number of strain components.

  rank     = nodes.rank ();
  strCount = STRAIN_COUNTS[ rank ];

  // Map dof types.

  JEM_ASSERT      ( rank == dofs->typeCount() );

  dofTypes.resize ( rank );

  for ( idx_t i = 0; i < rank; i++ )
  {
    idx_t k = dofs->findType ( DOF_TYPE_NAMES[i] );

    if ( k < 0_idx )
    {
      throw jem::IllegalInputException ( 
        context, 
        String::format ( "Dof tyep `%s' not found!", DOF_TYPE_NAMES[i] ) 
      );
    }

    dofTypes[i] = k;
  }

  // Get boundary nodes.

  for ( idx_t i = 0; i < 2*rank; i++ )
  {
    NodeGroup bnd = NodeGroup::get
             ( PBCGroupInputModule::EDGES[i], nodes, globdat, context );

    IdxVector inodes = bnd.getIndices ();

    JEM_ASSERT ( inodes.size() > 0_idx );

//    jem::System::out() << "EDGES: " << PBCGroupInputModule::EDGES[i] << " " << inodes << "\n";

    bndNodes[i].resize ( inodes.size() );
    bndNodes[i] = inodes;
  }

  // Get corner (controlling) nodes.

  ctrlNodes.resize ( rank + 1_idx );

  for ( idx_t i = 0; i < rank+1; i++ )
  {
    NodeGroup ctrl = NodeGroup::get
           ( PBCGroupInputModule::CORNERS[i], nodes, globdat, context );

    IdxVector inodes = ctrl.getIndices ();

    JEM_PRECHECK ( inodes.size() == 1_idx );

    ctrlNodes[i] = inodes[0];
  }

//  jem::System::out() << "ctrlNodes = " << ctrlNodes << "\n";

  // Compute RVE volume. 
  
  calcInvVol_ ();

  // Compute boundary matrix H.

  calcBndHMatrix_ ();

  // Determine dofs attached to ctrl. nodes and remaining ones.

  getCtrlAndRestDofs_ ();

  // Build a Constraints object for ConstrainedMatrix object to be used 
  // in homogenization process.

  // NOTE: 'cons' object has periodic constraints and prescribed 
  //       displacements at controlling nodes. The object 'cons' can have 
  //       non-zero values at controlling nodes and zero ones for setting 
  //       homogeneous constraints during the homogenization process 
  //       (micro problem). 
  //       On the other hand, 'cons2' object has only periodic 
  //       constraints to be used in the ConstrainedMatrix object in such 
  //       homogenized procedure (i.e., micro problem). 

  cons2 = jem::newInstance<Constraints> ( dofs ); 

  setPeriodicCons_ ( *cons2 );

/*
  Ref<jem::io::PrintWriter> prn = jem::newInstance<jem::io::PrintWriter> ( & jem::System::out() );
  dofs->printTo( *prn );
  cons->printTo( *prn );
  cons2->printTo( *prn );
*/

  // Build AbstractMatrix and solver for micro problem. 

  mBuilder = jem::newInstance<FEMatrixBuilder> 
                                         ( "microBuilder", elems, dofs ); 

  Ref<AbstractMatrix> microMatrix = mBuilder->getMatrix ();

/*
  solver = jem::newInstance<jive::solver::SkylineLU> 
                                    ( "microSolver", microMatrix, cons );
				    */

 /* 
  using jive::util::joinNames;
  msolver  = newSolver ( joinNames( myName_, MSOLVER_PROP ), 
                              conf, props, globdat );
*/

  consMat = jem::newInstance<ConstrainedMatrix>   ( microMatrix, cons2 );
}


//-----------------------------------------------------------------------
//   MicroPBCSolverModule::RunMicroData_::initMicroSolver
//-----------------------------------------------------------------------


void MicroPBCSolverModule::RunMicroData_::initMicroSolver

  ( const String&      name, 
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  using jem::newInstance;

  using jive::algebra::VectorSpace;
  using jive::algebra::StdVectorSpace;
  using jive::solver::newSolver;
  using jive::solver::SolverParams;


  // Build a VectorSpace object if required.

  Ref<VectorSpace> vspace = newInstance<StdVectorSpace> ( dofs );

  // SolverParams. Read preconditioners from property file if required.
  // For instance, the makeNew() function from GMRES solver tries to set
  // a preconditioner from property file.

  // Properties out = SolverParams::newInstance 
  //
  //        ( AbstractMatrix, VectorSpace, Constraints, Preconditioner );

  Properties solverParams = SolverParams::newInstance   

                        ( mBuilder->getMatrix(), vspace, cons, nullptr );

//  jem::System::out() << "NAME  " << name  << "\n";
//  jem::System::out() << "PROPS " << props.contains( name ) << "\n";

  msolver = newSolver ( name, conf, props, solverParams, globdat );

  msolver->configure  ( props );
  msolver->getConfig  ( conf  );
}


//-----------------------------------------------------------------------
//   MicroPBCSolverModule::RunMicroData_::setConstraints
//-----------------------------------------------------------------------


void MicroPBCSolverModule::RunMicroData_::setConstraints        

  ( const Vector&  strain ) const

{
  setCtrlNodeCons_ ( strain );
  setPeriodicCons_ ( *cons  );
}


//-----------------------------------------------------------------------
//   MicroPBCSolverModule::RunMicroData_::calcState
//-----------------------------------------------------------------------


void MicroPBCSolverModule::RunMicroData_::calcState

  ( MicroState&        mstate, 
    const Properties&  globdat ) const

{
  using jem::ALL;

  using jem::numeric::matmul;
  using jive::model::Actions;
  using jive::model::ActionParams;


  // Assemble the tangent micro stiffness matrix and compute the internal
  // force vector.
  
  const idx_t dofCount = dofs->dofCount ();

  Properties  microParams;

  Matrix      unit    ( strCount, strCount );

  Vector      strIncr ( strCount );
  Vector      fint    ( dofCount );

  unit = 0.;
  for ( idx_t i = 0; i < strCount; i++ ) unit(i,i) = 1.0;

  // Set homogeneous constraints for the micro model.

  //Ref<jem::io::PrintWriter> prn = jem::newInstance<jem::io::PrintWriter> ( & jem::System::out() );
  //*prn << "HOMOGENEOS. . .\n";

  Vector epsM ( strCount );
  epsM = 0.0;
  setConstraints ( epsM );

/*
  *prn << "cons . . .\n";
  cons->printTo ( *prn );
  *prn << ". . .\n";

  *prn << "cons2 . . .\n";
  cons2->printTo ( *prn );
  *prn << ". . .\n";
*/

  mBuilder->setToZero ();

  fint = 0.;

  microParams.set   ( ActionParams::MATRIX0   , mBuilder );
  microParams.set   ( ActionParams::INT_VECTOR, fint     );

  model->takeAction ( Actions::GET_MATRIX0, microParams, globdat );

  mBuilder->updateMatrix ();

  // Assemble micro strain-stress relations column by column.

  const idx_t bndDofCount = H.size ( 1 );

  //Matrix stiff  ( strCount, strCount );
  //Vector stress ( strCount    );

  Vector fi     ( bndDofCount ); 
  Vector hi     ( bndDofCount ); 
  Vector qi     ( bndDofCount );
  Vector rhs1   ( dofCount    );
  Vector rhs2   ( dofCount    );
  Vector sol    ( dofCount    );
  Vector tmp    ( dofCount    );

  Matrix Ht    = H.transpose ();

  for ( idx_t i = 0; i < strCount; i++ )
  {
    fi = matmul ( Ht, unit(ALL, i) );

    // Assemble r1 = {0 fi}^T vector.

    rhs1 = 0.;
    jem::select ( rhs1, ctrlDofs ) = fi;

    // Compute {gi hi}^T = \bar{K} * r1 vector.

    consMat->matmul  ( tmp, rhs1 );

    hi = jem::select ( tmp, ctrlDofs );

    // Assemble r2 = {gi 0}^T vector.

    rhs2 = 0.;
    jem::select ( rhs2, restDofs ) = jem::select ( tmp, restDofs );

    // Compute {ui 0}^T = \bar{K}^-1 * r2 vector. 

    //solver->solve ( sol, rhs2 );
    msolver->solve ( sol, rhs2 );

    // Compute {pi qi}^T = \bar{K} * {ui 0}^T vector.

    // tmp = 0.;
    consMat-> matmul ( tmp, sol );

    // Assemble/retrive {qi} vector. 

    qi = jem::select ( tmp, ctrlDofs );

    // Compute ith column of macroscopic tangent stiffness.

    mstate.stiff(ALL, i) = matmul ( H, std::move( Vector( hi - qi ) ) );
  }

  // Homogenized tangent moduli matrix and Homogenized stress.

  mstate.stiff *= v0Inv;

  jive::util::evalMasterDofs ( fint, *cons2 );

  Vector fc ( jem::select( fint, ctrlDofs ) );


  mstate.stress  = matmul  ( H, fc ); 
  mstate.stress *= v0Inv;


  mBuilder->setToZero ();
  mBuilder->clear     ();
}


//-----------------------------------------------------------------------
//   MicroPBCSolverModule::RunMicroData_::calcInvVol_
//-----------------------------------------------------------------------


void MicroPBCSolverModule::RunMicroData_::calcInvVol_ 

  ()

{
  // Compute area (or volume) based on corner node coordinates.

  // NOTE: for 2D cases, it is necessary to take into account thickness
  //       for the 'actual' volume.

  dxRVE.resize ( rank );

  Vector x0 ( rank );
  Vector x1 ( rank );

  nodes.getNodeCoords ( x0, ctrlNodes[0] );

  for ( idx_t i = 0; i < rank; i++ )
  {
    nodes.getNodeCoords ( x1, ctrlNodes[i+1] );

    dxRVE[i] = x1[i] - x0[i];
  }

  v0Inv = 1. / jem::product ( dxRVE );
}


//-----------------------------------------------------------------------
//   MicroPBCSolverModule::RunMicroData_::calcBndHMatrix_
//-----------------------------------------------------------------------


void MicroPBCSolverModule::RunMicroData_::calcBndHMatrix_

  ()

{
  // compute H_   0x 0y  xx  xy  yx  yy
  //          xx [ 0, 0, dx,  0,  0,  0 ]
  // 2D: H_ = yy [ 0, 0,  0,  0,  0, dy ]
  //          xy [ 0, 0,  0, dx, dy,  0 ]/2

  //              0x 0y 0z  xx  xy  xz  yx  yy  yz  zx  zy  zz 
  //          xx [ 0, 0, 0, dx,  0,  0,  0,  0,  0,  0,  0,  0 ]
  //          yy [ 0, 0, 0,  0,  0,  0,  0, dy,  0,  0,  0,  0 ]
  // 3D: H_ = zz [ 0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0, dz ]
  //          xy [ 0, 0, 0,  0, dx,  0, dy,  0,  0,  0,  0,  0 ]/2
  //          yz [ 0, 0, 0,  0,  0,  0,  0,  0, dy,  0, dz,  0 ]/2
  //          xz [ 0, 0, 0,  0,  0, dx,  0,  0,  0, dz,  0,  0 ]/2

  const idx_t nCtrlDofs = rank * ctrlNodes.size ();

  H.resize ( strCount, nCtrlDofs ); 

  H = 0.;

  if      ( rank == 2_idx )
  {
    H(0,2) = dxRVE[0];
    H(1,5) = dxRVE[1];
    H(2,3) = dxRVE[0] * .5;
    H(2,4) = dxRVE[1] * .5;
  }
  else if ( rank == 3_idx )
  {
    H(0,3)  = dxRVE[0];
    H(1,7)  = dxRVE[1];
    H(2,11) = dxRVE[2];

    H(3,4)  = dxRVE[0] * .5;
    H(3,6)  = dxRVE[1] * .5;

    H(4,8)  = dxRVE[1] * .5;
    H(4,10) = dxRVE[2] * .5;

    H(5,5)  = dxRVE[0] * .5;
    H(5,9)  = dxRVE[2] * .5;
  }
  else 
  {
    throw jem::IllegalInputException ( 
      context, 
      "At the moment, MicroPBCSolverModule only works for 2D and 3D problems!"
    );
  }
}


//-----------------------------------------------------------------------
//   MicroPBCSolverModule::RunMicroData_::getCtrlAndRestDofs_
//-----------------------------------------------------------------------


void MicroPBCSolverModule::RunMicroData_::getCtrlAndRestDofs_ 

  () 

{
  using jive::BoolVector;


  const idx_t dofCount  = dofs->dofCount        ();
  const idx_t nCtrlDofs = rank * ctrlNodes.size ();

  // First, determine dofs associated to ctrl. nodes.

  ctrlDofs.resize ( nCtrlDofs );

  const idx_t k = dofs->collectDofIndices 
                                       ( ctrlDofs, ctrlNodes, dofTypes );

  JEM_PRECHECK ( k == ctrlDofs.size() ); 

  // Then, determine the remaining dofs.

  BoolVector  dofMask ( dofCount );

  dofMask = true;

  jem::select ( dofMask, ctrlDofs ) = false;

  restDofs.resize ( jem::count( dofMask ) );

  idx_t i, j;

  for ( i = j = 0; i < dofCount; i++ )
  {
    if ( dofMask[i] )
    {
      restDofs[j++] = i;
    }
  }

  JEM_ASSERT ( j == restDofs.size() );
}


//-----------------------------------------------------------------------
//   MicroPBCSolverModule::RunMicroData_::setPeriodicCons_
//-----------------------------------------------------------------------


void MicroPBCSolverModule::RunMicroData_::setPeriodicCons_

  ( Constraints&  c ) const

{
  IdxVector idofs ( 2 );
  Vector    coefs ( 2 );

  coefs = 1.0;

  // Loop over faces.

  for ( idx_t ix = 0; ix < rank; ix++ )
  {
    IdxVector inodes ( bndNodes[2*ix]   );
    IdxVector jnodes ( bndNodes[2*ix+1] );
    idx_t     imaster = ctrlNodes[ ix+1 ];

/*
    jem::System::out() << "master = " << inodes << jem::io::endl;
    jem::System::out() << "slaves = " << jnodes << jem::io::endl;
    jem::System::out() << "ctrl   = " << imaster << jem::io::endl;
*/

    // Loop over the nodes of boundaries.

    for ( idx_t ine = 0; ine < inodes.size(); ine++ )
    {
      // Skip control points.

      if ( jnodes[ine] != imaster )
      {
        // Loop over dof types.

        for ( idx_t jx = 0; jx < rank; jx++ )
        {
          idofs[0] = dofs->getDofIndex ( inodes[ine], dofTypes[jx] );
          idofs[1] = dofs->getDofIndex ( imaster,     dofTypes[jx] );

          // Dofs of dependent (slave) node j.

          idx_t jdof = dofs->getDofIndex ( jnodes[ine], dofTypes[jx] );

          // Add Constraints.

	  idx_t jitem, item, iitem, doftype;

	  dofs->decodeDofIndex( jitem, doftype, jdof );
	  dofs->decodeDofIndex(  item, doftype, idofs[0] );
	  dofs->decodeDofIndex( iitem, doftype, idofs[1] );

	  ++jitem;
	  ++item;
	  ++iitem;

	  //jem::System::out()<< FEMUtils::DOF_TYPE_NAMES[jx] << "[" << jitem<< "] = " << item << " + " << iitem << "\n";

          c.addConstraint ( jdof, idofs, coefs );
        }
      }
    }
  }

  c.compress ();
}


//-----------------------------------------------------------------------
//   MicroPBCSolverModule::RunMicroData_::setCtrlNodeCons_
//-----------------------------------------------------------------------


void MicroPBCSolverModule::RunMicroData_::setCtrlNodeCons_

  ( const Vector&  strain ) const

{
  Vector vals ( ctrlDofs.size() );

  jem::numeric::matmul ( vals, H.transpose(), strain );

  // Prescribe displacement on the three control points (nodes).
  
  for ( idx_t i = 0; i < ctrlDofs.size(); i++ )
  {
    cons->addConstraint ( ctrlDofs[i], vals[i] );
  }

  cons->compress ();
}


//=======================================================================
//   class MicroPBCSolverModule 
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  MicroPBCSolverModule::TYPE_NAME        = "MicroPBC";

const char*  MicroPBCSolverModule::ESOLVER_PROP     = "esolver";
const char*  MicroPBCSolverModule::MSOLVER_PROP     = "msolver";
const char*  MicroPBCSolverModule::EPS_M_PROP       = "epsM";
const char*  MicroPBCSolverModule::EPS_M_INCR_PROP  = "epsMIncr";
const char*  MicroPBCSolverModule::EPS_M_FUNCS_PROP = "epsMFuncs";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


MicroPBCSolverModule::MicroPBCSolverModule

  ( const String&  name ) : Super ( name )

{
  solver_     = nullptr;
  out_        = nullptr;
  strCount_   = 0_idx;
  standalone_ = false;
}


MicroPBCSolverModule::~MicroPBCSolverModule

  ()

{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status MicroPBCSolverModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) 

{
  using jive::implict::newSolverModule;
  using jive::util::joinNames;

  using FEMUtils::STRAIN_COUNTS;


  // Init an equilibrium SolverModule.

  solver_ = newSolverModule ( joinNames( myName_, ESOLVER_PROP ), 
                              conf, props, globdat );

  solver_->configure ( props, globdat );
  solver_->getConfig ( conf , globdat );

  Status res = solver_->init ( conf, props, globdat );

  // Set a new micro run data object.

  Ref<RunMicroData_> newData = 
                        jem::newInstance<RunMicroData_> ( getContext() );

  newData->init            ( globdat );
  newData->initMicroSolver ( joinNames( myName_, MSOLVER_PROP ), 
                             conf, props, globdat );

  // Everything Ok, then, commit.

  runMData_.swap ( newData );

  // If a standalone mode is set, check if the given macro strain is 
  // compatible with problem settings.

  const idx_t rank = NodeSet::get( globdat, getContext() ).rank ();

  strCount_        = STRAIN_COUNTS[ rank ];

  if ( standalone_ )
  {
    if ( epsM_.size() != strCount_ || epsMFuncs_.size() != strCount_ )
    {
      throw jem::IllegalInputException ( 
        CLASS_NAME, 
        "macro strain vector is incompatible with problem dimension!"
      ); 
    }
  }

  // Resize current macro stress vector. It is used mostly for output 
  // purpose.
  
  sigM_.resize ( strCount_ );
  sigM_ = 0.0;

  return res;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------

/*
Module::Status MicroPBCSolverModule::run

  ( const Properties&  globdat )

{
  runMData_->epsM[2] = 0.1;

  runMData_->setCtrlNodeCons ();
  runMData_->setPeriodicCons ();

  Ref<jem::io::PrintWriter> prn = jem::newInstance<jem::io::PrintWriter> ( & jem::System::out() );

  *prn << "cons . . .\n";
  runMData_->cons->printTo ( *prn );
  *prn << ". . .\n";

  *prn << "cons2 . . .\n";
  runMData_->cons2->printTo ( *prn );
  *prn << ". . .\n";

  Status res = solver_->run  ( globdat );

 Vector disp (runMData_->dofs->dofCount());

  runMData_->calcState       ( globdat );

  // Postprocess data.

  return res;
}
*/


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void MicroPBCSolverModule::shutdown

  ( const Properties&  globdat )

{
  solver_->shutdown ( globdat );
}


//-----------------------------------------------------------------------
//   configure 
//-----------------------------------------------------------------------

 
void MicroPBCSolverModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties myProps = props.getProps  ( myName_ );

  standalone_ = false; 

  // Define a temp. function array and function arguments.

  Array< Ref<Function> > funcs;
  String args = "i, t";

  // First, check if EPS_M_PROP has been defined.

  if ( myProps.find( epsM_, EPS_M_PROP ) )
  {
    const idx_t strCount = epsM_.size ();

    epsM0_.resize ( strCount );

    epsM0_ = epsM_;

    standalone_ = true; // standalone mode

    Vector epsMIncr;

    if ( ! myProps.find( epsMIncr, EPS_M_INCR_PROP ) )
    {
      epsMIncr.resize ( strCount );

      epsMIncr = 0.0;
    }

    JEM_PRECHECK2 ( strCount <= epsMIncr.size(), "mismatch size!" );

    // Define functions associated to this setting.

    funcs.resize ( strCount );
    funcs = nullptr;

    for ( idx_t icomp = 0; icomp < strCount; icomp++ )
    {
      String ifuncExp = String( epsM0_[icomp] )   + " + "   +
                        String( epsMIncr[icomp] ) + "*i + " + 
                        String( epsMIncr[icomp] ) + "*t"; 

      funcs[icomp] = FuncUtils::newFunc ( args, ifuncExp, globdat ); 
    }
  }

  if ( myProps.contains( EPS_M_FUNCS_PROP ) )
  {
    standalone_ = true; // standalone mode

    FuncUtils::configFuncs ( funcs, args,
                             EPS_M_FUNCS_PROP,
                             myProps, globdat );
  
    const idx_t strCount = funcs.size ();

    epsM_.resize  ( strCount );
    epsM_ = 0.0;

    epsM0_.resize ( strCount );
    epsM0_ = 0.0;
  }

  // Everything if OK, then, commit.

  epsMFuncs_.swap ( funcs );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void MicroPBCSolverModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties myConf = conf.makeProps ( myName_ );

  myConf.set ( "standalone", standalone_ );

  if ( standalone_ )
  {
    FuncUtils::getConfig ( myConf, epsMFuncs_, EPS_M_FUNCS_PROP );
  }
}


//-----------------------------------------------------------------------
//   advance
//-----------------------------------------------------------------------


void MicroPBCSolverModule::advance

  ( const Properties&  globdat )

{
  if ( standalone_ )
  {
    // Scale macro strain field.

    try
    {
      FuncUtils::resolve ( epsMFuncs_, globdat );
      scaleMacroStrain_  ( globdat );
    }
    catch ( jem::Exception& ex )
    {
      ex.setContext ( getContext() );
      throw;
    }

    runMData_->setConstraints ( epsM_ );

/*
  Ref<jem::io::PrintWriter> prn = jem::newInstance<jem::io::PrintWriter> ( & jem::System::out() );
  
  *prn << "cons . . .\n";
  runMData_->cons->printTo ( *prn );
  *prn << ". . .\n";

  *prn << "cons2 . . .\n";
  runMData_->cons2->printTo ( *prn );
  *prn << ". . .\n";
*/
  }

  solver_->advance ( globdat );
}


//-----------------------------------------------------------------------
//   solve
//-----------------------------------------------------------------------


void MicroPBCSolverModule::solve

  ( const Properties&  info,
    const Properties&  globdat )

{
  solver_->solve ( info, globdat );
}


//-----------------------------------------------------------------------
//   cancel
//-----------------------------------------------------------------------


void MicroPBCSolverModule::cancel

  ( const Properties&  globdat )

{
  solver_->cancel ( globdat );
}


//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------


bool MicroPBCSolverModule::commit

  ( const Properties&  globdat )

{
  bool ok = solver_->commit ( globdat );

  if ( ok && standalone_ )
  {
    Ref<MicroState> mstate = jem::newInstance<MicroState> ();

    mstate->stiff .resize ( strCount_, strCount_ );
    mstate->stress.resize ( strCount_ );

    runMData_->calcState  ( *mstate, globdat );

    // Update initial macro strain field.

    epsM0_ = epsM_;

    // Write output.

    writeStates_ ( *mstate, globdat );

    // Update current macro stress field.

    sigM_ = mstate->stress;
  }

  // Store the current macro strain and stress fields in globdat variable 
  // so that one can use it with SampleModule, for instance.

  Properties  myVars = 

                   jive::util::Globdat::getVariables ( myName_, globdat );

  myVars.set ( "epsM", epsM_ );
  myVars.set ( "sigM", sigM_ );

  return ok;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void MicroPBCSolverModule::setPrecision

  ( double  eps )

{
  solver_->setPrecision ( eps );
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double MicroPBCSolverModule::getPrecision 

  () const

{
  return solver_->getPrecision ();
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool MicroPBCSolverModule::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using MicroUtils::MicroActions;
  using MicroUtils::MicroActionParams;


  if ( action == MicroActions::GET_MACRO_STRAIN )
  {
    if ( standalone_ )
    {
      throw jem::IllegalInputException ( 
        CLASS_NAME, 
	"it only works in a multiscale framework, "
	"not in a `standalone' mode!" 
      );
    }

    params.get  ( epsM_, MicroActionParams::MACRO_STRAIN );

    runMData_->setConstraints ( epsM_ );

    return true;
  }

  if ( action == MicroActions::GET_MICRO_STATE )
  {
    const RunMicroData_&  md = * runMData_; 

    Ref<MicroState>       mstate;

    params.get            ( mstate, MicroActionParams::MICRO_STATE );

    mstate->stiff .resize ( md.strCount, md.strCount );
    mstate->stress.resize ( md.strCount );

    mstate->stiff  = 0.;
    mstate->stress = 0.;

    md.calcState ( *mstate, globdat );

    // Update macro stress field.

    sigM_ = mstate->stress;

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   makeNew 
//-----------------------------------------------------------------------


Ref<Module> MicroPBCSolverModule::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  return jem::newInstance<Self> ( name );
}


//-----------------------------------------------------------------------
//   declare 
//-----------------------------------------------------------------------


void MicroPBCSolverModule::declare

  ()

{
  using jive::app::ModuleFactory;
  
  ModuleFactory::declare ( TYPE_NAME , & makeNew );
  ModuleFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   scaleMacroStrain_
//-----------------------------------------------------------------------


void MicroPBCSolverModule::scaleMacroStrain_

  ( const Properties&  globdat )

{
  using jive::util::Globdat;


  double args[5];

  double t = 0.0;

  idx_t  j = 0;

  globdat.find ( j, Globdat::TIME_STEP );
  globdat.find ( t, Globdat::TIME      );

  args[0] = (double) j;
  args[1] = t;

  for ( idx_t icomp = 0; icomp < strCount_; icomp++ )
  {
    epsM_[icomp] = epsMFuncs_[icomp]->getValue ( args );
  }

  jem::System::out() << CLASS_NAME 
                     << ": updated macro strain = " << epsM_ << "\n";
}


//-----------------------------------------------------------------------
//   getWriter_
//-----------------------------------------------------------------------


Ref<PrintWriter>  MicroPBCSolverModule::getWriter_

  () const

{
  using jem::newInstance;

  using jem::util::StringUtils;
  using jem::io::FileWriter;

  using jive::StringVector;


  // File for output.

  StringVector fname ( 2 );

  fname[0] = "MicroPBC";
  fname[1] = "out";

  return newInstance<PrintWriter> ( 

    newInstance<FileWriter> ( StringUtils::join( fname, "." ) ) 

  );
} 


//-----------------------------------------------------------------------
//   writeStates_
//-----------------------------------------------------------------------


void  MicroPBCSolverModule::writeStates_

  ( const MicroState&  mstate, 
    const Properties&  globdat ) const

{
  using jem::io::endl;

  using jive::util::Globdat;


  // Get current time step.

  idx_t istep;

  globdat.get ( istep, Globdat::TIME_STEP );

  // Get writer.

  if ( out_ == nullptr )
  {
    const_cast<Self*> ( this )->out_ = getWriter_ ();

    out_->nformat.setFractionDigits ( 8    );
    out_->nformat.setScientific     ( true );
    out_->nformat.setFloatWidth     ( 16   );
  }

  *out_ << "TIME_STEP: " << istep << endl;

  *out_ << "stress:"     << endl;
  *out_ << mstate.stress << endl;

  *out_ << "stiff:"      << endl;
  *out_ << mstate.stiff  << endl;

  *out_ << endl;

  out_->flush ();
}

