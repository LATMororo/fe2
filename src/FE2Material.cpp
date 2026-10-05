
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This material model provides an interface between macro- and micro- 
 *  models of a multiscale analysis [1-3]. For each integration point, a 
 *  chain/models of FE modules and models are set. 
 *
 *  This framework follows what has been implemented by Iuri Rocha and 
 *  his collegues at TU Delft [3-5]. Unlike their code, micro solvers are
 *  provided separately, which are responsible for solving the micro BVP 
 *  (Boundary Value Problem) subjected to appropriate BCs obtained from 
 *  macro strain field, and for homogenizing state quantities, i.e., 
 *  macro stress and stiffness.
 *
 *  [1] R. J. M. Smit, W. A. M. Brekelmans, and H. E. H. Meijer. 
 *      Prediction of the mechanical behavior of nonlinear heterogeneous 
 *      systems by multi-level finite element modeling. Computer Methods 
 *      in Applied Mechanics and Engineering 155 (1) (1998), pp. 181–192.
 *
 *  [2] D. Peric, E. A. de Souza Neto, R. A. Feijo, M. Partovi, and 
 *      A. J. C. Molina. On micro-to-macro transitions for multi-scale 
 *      analysis of non-linear heterogeneous materials: unified 
 *      variational basis and finite element implementation. Int. J. 
 *      Numer. Meth. Engng, 87 (1-5) (2011), pp. 149–170.
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
 *  Modified: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br, luiztaumaturgo@gmail.com
 *  Date: Jul 2026
 *
 *  Allow the user to set a SampleModule settings for some integration 
 *  points, i.e., some RVEs (or all of them).
 * 
 */


#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/array/operators.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/RuntimeException.h>
#include <jem/base/System.h>
#include <jem/io/FileWriter.h>
#include <jem/io/PrintWriter.h>
#include <jem/io/Logger.h>
#include <jem/io/PatternLogger.h>


#include <jive/app/ChainModule.h>
#include <jive/app/UserconfModule.h>
#include <jive/fem/InitModule.h>
#include <jive/fem/NodeSet.h>
#include <jive/fem/ShapeModule.h>
#include <jive/util/Globdat.h>
#include <jive/app/OutputModule.h>
#include <jive/app/SampleModule.h>


#include "FE2Material.h"
#include "MaterialFactory.h"
#include "MicroUtils.h"


JEM_DEFINE_CLASS ( FE2Material );


using namespace jem::literals;


using MicroUtils::MicroState;


//=======================================================================
//   class FE2Material::MicroChain_ 
//=======================================================================


class FE2Material::MicroChain_ : public Object
{
 public:

  explicit                  MicroChain_ 

    ( idx_t                   ip,
      const Properties&       conf,
      const Properties&       props );

  idx_t                     rank                  () const;

  void                      initOutput

    ( const Properties&       conf,
      const Properties&       props );

  void                      initSample

    ( const Properties&       conf,
      const Properties&       props );

  void                      update 
    
    ( const Vector&           stress,
      const Matrix&           stiff,
      const Vector&           strain )               const;

  void                      advance               () const;
  void                      commit                ();


 protected:

  virtual                  ~MicroChain_           ();
 

 public:

  Ref<jive::app::ChainModule> 
                            chain;

  Ref<jive::app::OutputModule>
                            out;

  Ref<jive::app::SampleModule>
                            sample;

  Ref<MicroState>           state;

  Properties                globdat;
  Properties                conf;
  Properties                props;

  String                    myID;
 

 private:

  bool                      commit_;

};


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


FE2Material::MicroChain_::MicroChain_ 

  ( idx_t              ip,
    const Properties&  conf,
    const Properties&  props )

{
  using jem::newInstance;

  using jive::app::ChainModule;
  using jive::app::UserconfModule;
  using jive::fem::InitModule;
  using jive::fem::NodeSet;
  using jive::fem::ShapeModule;
  using jive::util::Globdat;

  
  myID = String ( ip );

  // Create a new globdat for the new micro model at a given integration
  // point, i.e., 'ip'

  globdat = Globdat::newInstance ( "microPoint" + myID ); 

  // ... also, create a new runtime variable for it.

  Globdat::getVariables ( globdat );

  Globdat::initStep     ( globdat );

  // Create a chain of modules w.r.t. the micro problem.

  chain = newInstance<ChainModule> ( "micro.chain" );

  chain->pushBack ( newInstance<UserconfModule>( "micro.userInput" ) );
  chain->pushBack ( newInstance<ShapeModule>   (                   ) );
  chain->pushBack ( newInstance<InitModule>    ( "micro.init"      ) );
  chain->pushBack ( newInstance<UserconfModule>( "micro.extra"     ) );

  chain->configure ( props, globdat );
  chain->getConfig ( conf , globdat );
  chain->init      ( conf , props, globdat );

  // Set a MicroState object for the given integration point.

  state = newInstance<MicroState> ();

  // Set commit member class.

  commit_ = true;

  // Default members.
  
  out    = nullptr;
  sample = nullptr;
}


FE2Material::MicroChain_::~MicroChain_ 

  ()

{
  if ( out )
  {
    out->shutdown ( globdat );
  }

  if ( sample )
  {
    sample->shutdown ( globdat );
  }
}


//-----------------------------------------------------------------------
//   FE2Material::MicroChain_::rank 
//-----------------------------------------------------------------------


idx_t FE2Material::MicroChain_::rank 

  () const

{
  const String context = "MicroChain constructor (" + myID + ")" ;

  return jive::fem::NodeSet::get( globdat, context ).rank ();
}


//-----------------------------------------------------------------------
//   FE2Material::MicroChain_::initOutput
//-----------------------------------------------------------------------


void FE2Material::MicroChain_::initOutput 

  ( const Properties&  conf,
    const Properties&  props )

{
  out = jem::newInstance<jive::app::OutputModule> ( "micro.output" );

  props.set ( "micro.output.file", "outIP." + myID + ".%i.mech" );

  out->configure ( props, globdat );
  out->getConfig ( conf , globdat );
  out->init      ( conf , props, globdat );
}


//-----------------------------------------------------------------------
//   FE2Material::MicroChain_::initSample
//-----------------------------------------------------------------------


void FE2Material::MicroChain_::initSample 

  ( const Properties&  conf,
    const Properties&  props )

{
  sample = jem::newInstance<jive::app::SampleModule> ( "micro.sample" );

  props.set ( "micro.sample.file", "samp." + myID + ".out" );

  sample->configure ( props, globdat );
  sample->getConfig ( conf , globdat );
  sample->init      ( conf , props, globdat );
}


//-----------------------------------------------------------------------
//   FE2Material::MicroChain_::update
//-----------------------------------------------------------------------


void FE2Material::MicroChain_::update

  ( const Vector&  stress,
    const Matrix&  stiff,
    const Vector&  strain ) const

{
  using MicroUtils::MicroActions;
  using MicroUtils::MicroActionParams;


  Properties params;

  bool       ok = false;

  
  // Set a macro strain field.

  params.set         ( MicroActionParams::MACRO_STRAIN, strain );

  ok = chain->takeAction 
                     ( MicroActions::GET_MACRO_STRAIN, params, globdat );

  if ( ! ok ) 
  {
    throw jem::RuntimeException ( 
      CLASS_NAME,  
      String::format ( 
      "Macro strain field has not been set at the "
      "integration point (%d)", myID )  
    );
  }

  // Advance material. Basically, update micro/chain TIME_STEP and 
  // OLD_TIME_STEP.

  advance ();

  // Solve the micro model.

  try
  {
    chain->run ( globdat );
  }
  catch ( const jem::Exception& ex )
  {
    jem::System::out() << "Micromodel " << myID << " with strain " 
                       << strain << " did not converge\n";
    
    throw;
  }

  // Compute homogenized/macro state: stress and stiffness matrix.

  params.clear ();

  params.set   ( MicroActionParams::MICRO_STATE, state );

  ok = chain->takeAction
               ( MicroActions::GET_MICRO_STATE, params, globdat );

  if ( ! ok ) 
  {
    throw jem::RuntimeException ( 
      CLASS_NAME,  
      String::format ( 
      "Homogenized stress and stiffness quantities have been not "
      "obtained at the integration point (%d)", myID )  
    );
  }

  stress = 0.0;
  stiff  = 0.0;

  stress = state->stress;
  stiff  = state->stiff;

  // Update commit member class. Macro model might require more 
  // iterations, and then, micro model might be requested to update its 
  // own model; however, micro model/chain is not actually committed.
  // This is necessary to avoid advancing micro TIME_STEP and 
  // OLD_TIME_STEP when macro model is not actually committed.

  const_cast<FE2Material::MicroChain_*>( this )->commit_ = false;
}


//-----------------------------------------------------------------------
//   FE2Material::MicroChain_::advance
//-----------------------------------------------------------------------


void FE2Material::MicroChain_::advance

  () const 

{
  using jive::util::Globdat;


  if ( commit_ )
  {
    Globdat::advanceStep ( globdat );
  }
  else
  {
    idx_t istep, istep0;

    globdat.get ( istep , Globdat::TIME_STEP     );
    globdat.get ( istep0, Globdat::OLD_TIME_STEP );

    if ( istep == 1_idx )
    {
      istep  = 1_idx;
      istep0 = 0_idx;
    }
    else
    {
      --istep;
      --istep0;
    }

    globdat.set ( Globdat::TIME_STEP    , istep  );
    globdat.set ( Globdat::OLD_TIME_STEP, istep0 );
  }
}


//-----------------------------------------------------------------------
//   FE2Material::MicroChain_::commit
//-----------------------------------------------------------------------


void FE2Material::MicroChain_::commit 

  ()

{
  if ( out )
  {
    out->run ( globdat );
  }

  if ( sample )
  {
    sample->run ( globdat );
  }

  jive::util::Globdat::commitStep ( globdat );

  // Update commit member class.

  commit_ = true;
}



//=======================================================================
//   class FE2Material 
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  FE2Material::TYPE_NAME     = "FE2";

const char*  FE2Material::STATE_PROP    = "state";
const char*  FE2Material::OUT_IPS_PROP  = "outIPs";
const char*  FE2Material::SAMP_IPS_PROP = "sampIPs";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


FE2Material::FE2Material 

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) : 

    Super ( name )

{
  // Get material Properties.

  props_ = props.getProps ( myName_ );
  conf_  = conf.makeProps ( myName_ );

  // Get "dim" (or "rank" of material) parameter.

  rank_ = -1_idx;

  props_.get ( rank_, RANK_PROP, 1_idx, 3_idx ); 
  conf_ .set ( RANK_PROP, rank_ );

  // Set macro the number of strain components.

  strCount_ = FEMUtils::STRAIN_COUNTS[ rank_ ];
}


FE2Material::~FE2Material 

  ()

{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void FE2Material::configure 

  ( const Properties&  props,
    const Properties&  globdat )

{
  if ( ! props.contains( myName_ ) )
  {
    props_ = props.getProps ( myName_ );
  }

  // Get the rank/dim and state of material.

  if ( rank_ == 2_idx )
  {
    props_.get( stateString_, STATE_PROP );

    if      ( stateString_ == "PlaneStrain" )
    {
      state_ = ProblemType::PlaneStrain;
    }
    else if ( stateString_ == "PlaneStress" )
    {
      state_ = ProblemType::PlaneStress;
    }
  }

  if ( rank_ == 3_idx )
  {
    state_ = ProblemType::Solid;
     
    stateString_ = "Solid";
  }

  // Get integration-point list for output.

  props_.find ( outIPs_, OUT_IPS_PROP );

  // Get integration-point list for sample output.

  props_.find ( sampIPs_, SAMP_IPS_PROP );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void FE2Material::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const 

{
  if ( ! conf_.contains( myName_ ) )
  {
    const_cast<Self*>( this )->conf_ = conf.makeProps ( myName_ );
  }

  conf_.set ( STATE_PROP   , stateString_ );
  conf_.set ( OUT_IPS_PROP , outIPs_      );
  conf_.set ( SAMP_IPS_PROP, sampIPs_     );
}


//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------


void FE2Material::commit 

  ()

{
  // Commit all micro models/chains.

  for ( idx_t ip = 0; ip < microChains_.size(); ip++ )
  {
    microChains_.getAs<MicroChain_>( ip )->commit ();
  }
}


//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------


void FE2Material::allocPoints 

  ( idx_t  count )

{
  using jem::newInstance;
  using jem::System;

  using jem::io::FileWriter;
  using jem::io::Logger;
  using jem::io::PatternLogger;
  using jem::io::PrintWriter;
  using jem::io::Writer;
  using jem::util::ObjFlex;

  
  JEM_ASSERT ( count > 0_idx );


  Ref<Writer> macroOut = & System::out ();
  Ref<Writer> microOut = newInstance<PrintWriter> 
                              ( newInstance<FileWriter>( "micro.log" ) );


  Ref<Logger> macroLog = newInstance<PatternLogger> 
                                                  ( macroOut, "*.info" );
  Ref<Logger> microLog = newInstance<PatternLogger> 
                                                  ( microOut, "*.info" );

  macroOut->flush       ();

  System::setLogger     ( microLog );
  System::setOutStream  ( microOut );
  System::setWarnStream ( microOut );

  // Create micro models/chains.

  ObjFlex newChains;

  microRank_     = -1_idx;

  const bool hasOutIp  = outIPs_.size() > 0_idx ? true : false;
  const bool hasOutput = props_.contains ( "micro.output" );

  const bool hasSampIp = sampIPs_.size() > 0_idx ? true : false;
  const bool hasSample = props_.contains ( "micro.sample" );

  // Check if an OutputModule setting has been defined.

  if ( hasOutIp && ! hasOutput )
  {
    throw jem::IllegalInputException ( 
      CLASS_NAME, "OutputModule has not been set!" 
    );
  }

  // Check if a SampleModule setting has been defined.

  if ( hasSampIp && ! hasSample )
  {
    throw jem::IllegalInputException ( 
      CLASS_NAME, "SampleModule has not been set!" 
    );
  }

  // Loop over all integration points and build micro chains/models.

  for ( idx_t ip = 0; ip < count; ip++ )
  {
    Ref<MicroChain_> newChain = newInstance<MicroChain_> 
                                                   ( ip, conf_, props_ );

    // Check if there is a unique 'rank' for all  micro models.

    if ( ip == 0_idx ) 
    {
      microRank_ = newChain->rank (); 
    }
    else
    {
      if ( microRank_ != newChain->rank() )
      {
        throw jem::IllegalInputException 
                   ( CLASS_NAME, "Region with different micro models!" );
      }
    }

    JEM_ASSERT ( 
      microRank_ == rank_ || ( microRank_ == 3_idx && rank_ == 2_idx ) 
    );

    // Check if it is required to output results.

    // 1. Write results for all integration points.

    if ( ! hasOutIp && hasOutput )
    {
      newChain->initOutput ( conf_, props_ );
    }

    // 2. Write results only of a set of integration points in 
    // outIPs_ member.

    if ( hasOutIp && jem::testany( outIPs_ == ( int ) ip ) )
    {
      newChain->initOutput ( conf_, props_ );
    }

    // Check if it is required to output sample stuff.

    // 1. Write sample results for all integration points.

    if ( ! hasSampIp && hasSample )
    {
      newChain->initSample ( conf_, props_ );
    }

    // 2. Write results only of a set of integration points in 
    // sampIPs_ member.

    if ( hasSampIp && jem::testany( sampIPs_ == ( int ) ip ) )
    {
      newChain->initSample ( conf_, props_ );
    }

    newChains.insert ( ip, std::move( newChain ) );
  }

  // Everything Ok, then, commit it.

  newChains.trimToSize ();

  microChains_.swap    ( newChains );

  // Define permutation vector if macro and micro ranks are different.

  if ( microRank_ > rank_ )
  {
    iperm23_.resize ( 3 );

    iperm23_[0] = 0_idx;  // micro_xx == macro_xx
    iperm23_[1] = 1_idx;  // micro_yy == macro_yy
    iperm23_[2] = 3_idx;  // micro_xy == macro_xy
  }

  // Define microStrCount_ member.

  microStrCount_ = FEMUtils::STRAIN_COUNTS[ microRank_ ];

  microLog->flush ();

  System::setLogger     ( macroLog );
  System::setOutStream  ( macroOut );
  System::setWarnStream ( macroOut );
}


//-----------------------------------------------------------------------
//   rank
//-----------------------------------------------------------------------


idx_t FE2Material::rank

  () const noexcept

{
  return rank_;
}


//-----------------------------------------------------------------------
//   getProblemType
//-----------------------------------------------------------------------


ProblemType FE2Material::getProblemType 

  () const noexcept

{
  return state_;
}


//-----------------------------------------------------------------------
//   elasticUpdate
//-----------------------------------------------------------------------


void FE2Material::elasticUpdate

  ( const Vector&  stress,
    const Matrix&  stiff,
    const Vector&  strain ) const

{
  throw jem::RuntimeException ( 
    CLASS_NAME, "elasticUpdate() is not supported!" 
  );
}


Vector FE2Material::elasticUpdate

  ( const Vector&  strain ) const

{
  throw jem::RuntimeException ( 
    CLASS_NAME, "elasticUpdate() is not supported!" 
  );
}


//-----------------------------------------------------------------------
//   update 
//-----------------------------------------------------------------------


void FE2Material::update 

  ( const Vector&  stress,
    const Matrix&  stiff,
    const Vector&  strain,
    idx_t          ip ) const 

{
  Matrix microStiff   ( microStrCount_, microStrCount_ );
  Vector microStress  ( microStrCount_ );
  Vector macroStrain  ( microStrCount_ );


  jem::io::Writer&  out = jem::System::out ();


  // Get micro model/chain attached to a given integration point on the 
  // macro level.

  const MicroChain_& mc = * microChains_.getAs<MicroChain_> ( ip );

  // Add missing strain componentes in case of problems using a 3D RVE 
  // on 2D macro domains.

  completeStrain_ ( macroStrain, strain );

  // Compute homogenized quantities: stress and stiffness.

  print ( out, "\n", CLASS_NAME, ": homogenizing i.p. `", ip, "'\n\n" );

  mc.update ( microStress, microStiff, macroStrain );

  microStressToMacroStress_ ( stress, microStress  );
  microStiffToMacroStiff_   ( stiff , microStiff   );
}


void FE2Material::update 

  ( const Vector&  stress,
    const Vector&  strain,
    idx_t          ip ) const 

{
  // NOTE: unlike the same function update() (see above), which solves a
  //       micro BVP, this one only gets stress at the converged state;
  //       therefore, it is intended to be used in post-processing tasks.

  // Get micro model/chain attached to a given integration point on the 
  // macro level.

  const MicroChain_& mc = * microChains_.getAs<MicroChain_> ( ip );

  // Get micro stress.

  microStressToMacroStress_ ( stress, mc.state->stress );
}


//-----------------------------------------------------------------------
//   getStiffMat
//-----------------------------------------------------------------------


Matrix FE2Material::getStiffMat           

  () const 

{
  throw jem::RuntimeException ( 
    CLASS_NAME, "getStiffMat() is not supported!" 
  );
}


//-----------------------------------------------------------------------
//   fill3DStress
//-----------------------------------------------------------------------


Tuple<double,6> FE2Material::fill3DStress

  ( const Vector&  v3 ) const

{
  throw jem::RuntimeException ( 
    CLASS_NAME, "fill3DStress() is not supported!" 
  );
}


//-----------------------------------------------------------------------
//   fill3DStrain
//-----------------------------------------------------------------------


Tuple<double,6> FE2Material::fill3DStrain

  ( const Vector&  v3 ) const

{
  throw jem::RuntimeException ( 
    CLASS_NAME, "fill3DStrain() is not supported!" 
  );
}


//-----------------------------------------------------------------------
//   makeNew 
//-----------------------------------------------------------------------

 
Ref<Material> FE2Material::makeNew

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


void FE2Material::declare

  ()

{
  MaterialFactory::declare ( TYPE_NAME,  & makeNew );
  MaterialFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   completeStrain_
//-----------------------------------------------------------------------


void FE2Material::completeStrain_

  ( const Vector&  microStrain,
    const Vector&  macroStrain ) const
{
  // Fill in (add missing strain components) the micro strain vector in 
  // case micro-rank = 3 and macro-rank = 2.

  if ( rank_ == microRank_ )
  {
    microStrain = macroStrain;
  }
  else
  {
    microStrain = 0.0;
    microStrain[iperm23_] = macroStrain;
  }
}


//-----------------------------------------------------------------------
//   microStressToMacroStress_ 
//-----------------------------------------------------------------------


void FE2Material::microStressToMacroStress_

  ( const Vector&  macroStress,
    const Vector&  microStress ) const

{
  // Get plane stress components in case micro-rank = 3 and 
  // macro-rank = 2.

  if ( rank_ == microRank_ )
  {
    macroStress = microStress;
  }
  else
  {
    macroStress = microStress[iperm23_];
  }
}


//-----------------------------------------------------------------------
//   microStiffToMacroStiff_
//-----------------------------------------------------------------------


void FE2Material::microStiffToMacroStiff_

  ( const Matrix&  macroStiff,
    const Matrix&  microStiff ) const

{
  // Get plane stiff components in case micro-rank = 3 and 
  // macro-rank = 2.

  if ( rank_ == microRank_ )
  {
    macroStiff = microStiff;
  }
  else
  {
    macroStiff = microStiff ( iperm23_, iperm23_ );
  }
}
