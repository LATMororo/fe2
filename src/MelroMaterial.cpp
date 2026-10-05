
/*
 *  Copyright (C) 2014 TU Delft. All rights reserved.
 *  
 *  This class implements the material model for polymers
 *  from Melro et al. (2013)
 *  
 *  Author: F.P. van der Meer, f.p.vandermeer@tudelft.nl
 *  Date:   October 2014
 *
 *  Changelog:
 *
 *  Modified: Luiz Antonio T. Mororo, October 2017
 *  Date:     October 2017
 *
 *  This class is a child material member of DamTLSMaterial
 *  classes. 
 *
 *  Modified: Luiz Antonio T. Mororo, luiz.mororo@ifce.edu.br
 *  Date:     June, 2025
 *  
 *  Clean-up code, and refactored member functions for better 
 *  performance, storage scheme, and to meet new Jem/Jive (version 3.0) 
 *  data types, for instance, `idx_t'. Besides, desperate mode tasks, 
 *  used in oscillating equilibrium point, was removed 
 */


#include <jem/base/ClassTemplate.h>
#include <jem/base/Float.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/utilities.h>


#include <jive/util/FuncUtils.h>


#include "MelroMaterial.h"
#include "MaterialFactory.h"


JEM_DEFINE_CLASS ( MelroMaterial );


using namespace jem::literals;


using jem::numeric::Function;


//=======================================================================
//   class MelroMaterial::Hist_ 
//=======================================================================


//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


MelroMaterial::Hist_::Hist_

  () : epspeq ( 0. ), dissipation ( 0. )

{
  epsp = 0.;
  loading = false;
}

//-----------------------------------------------------------------------
//   MelroMaterial::Hist_::toVector
//-----------------------------------------------------------------------


inline void MelroMaterial::Hist_::toVector

 ( const Vector&  vec ) const

{
  vec[0] = epsp[0];
  vec[1] = epsp[1];
  vec[2] = epsp[2];
  vec[3] = epsp[3];
  vec[4] = epsp[4];
  vec[5] = epsp[5];
  vec[6] = epspeq;
  vec[7] = dissipation;
  vec[8] = loading;
}


//=======================================================================
//   class MelroMaterial::YieldFunc_
//=======================================================================


class MelroMaterial::YieldFunc_ : public jem::Object
{
 public:

                            YieldFunc_  
      
    ( double                  young,
      double                  poisson );

  void                      configure

    ( const Properties&       props,
      const Properties&       globdat );
  
  void                      getConfig 

    ( const Properties&       conf )                 const;

  bool                      isPlastic

    ( const Invariants&       inv,
      double                  epspeq0);

  double                    findRoot 
      
    ( double                  dgam0 );

  void                      findBounds

    ( double&                 gmin,
      double&                 gmax,
      double&                 fmin,
      double&                 fmax );

  double                    findRoot

    ( double                  gmin,
      double                  gmax,
      double                  fmin,
      double                  fmax );

  void                      getTangentParameters 

    ( double&                 betam, 
      double&                 pbm, 
      double&                 hfac, 
      double&                 beta, 
      double&                 phi, 
      double&                 rho, 
      double&                 chi, 
      double&                 psi, 
      double&                 xi, 
      double&                 omega, 
      double                  dgam );

  void                      getTangentParameters 

    ( double&                 beta, 
      double&                 phi, 
      double&                 rho, 
      double&                 chi, 
      double&                 psi ); 

  void                      getPGradientParameters

    ( double&                 phiF,
      double&                 rhoF,
      double&                 chiF,
      double&                 psiF )              const;

  inline double             getAlpha              () const noexcept;
  inline double             getZetaP              () const noexcept;
  inline double             getZetaS              () const noexcept;
  inline double             getEpspeq             () const noexcept;

    //inline double         getK             () const;
    //inline double         getG             () const;

 protected:

  virtual                  ~YieldFunc_            ();


 private:

  double                    evalTrial_            () const;

  double                    eval_

    ( double                  dgam );

  void                      setDGam_ 
      
    ( double                  dgam );

  void                      resetDGam_            ();

  double                    estimateRoot_

    ( double                  gmin,
      double                  gmax,
      double                  fmin,
      double                  fmax )                 const;

  void                      improveBounds_

    ( double&                 gmin,
      double&                 gmax,
      double&                 fmin,
      double&                 fmax );


 private:

  Ref<Function>             sigmaT_;
  Ref<Function>             sigmaC_;

  double                    G_;
  double                    K_;
  double                    alpha_;
  double                    alpha2_;
  double                    nuPFac_;
  double                    Ka_;

  double                    rmTol_;
  double                    poissonP_;

  double                    epspeq0_;
  double                    j2tr_;
  double                    i1tr_;
  double                    j2fac_;
  double                    i1fac_;

  double                    zetaP_;
  double                    zetaS_;
  double                    zS2_;
  double                    zP2_;
  double                    depspeq_;
  double                    epspeq_;
  double                    sigC_;
  double                    sigT_;
  double                    sigct_;
  double                    HC_;
  double                    HT_;
  double                    gjzs3_;
  double                    kaizp2_;
  double                    dcdg_;

  idx_t                     maxIter_;

};


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


MelroMaterial::YieldFunc_::YieldFunc_

  ( double  young,
    double  poisson )

{
  G_ = young / 2. / ( 1. + poisson );
  K_ = young / 3. / ( 1. - 2. * poisson );

  rmTol_    = 1.e-10;
  poissonP_ = 0;
  maxIter_  = 10_idx;
}


MelroMaterial::YieldFunc_::~YieldFunc_ 

  ()

{}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::configure
//-----------------------------------------------------------------------


void MelroMaterial::YieldFunc_::configure 

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jive::util::FuncUtils;


  props.find  ( rmTol_   , RM_TOL_PROP          );
  props.get   ( poissonP_, PLASTIC_POISSON_PROP );
  props.find  ( maxIter_ , RM_MAX_ITER_PROP     );

  alpha_  = ( 4.5 - 9. * poissonP_ ) / ( 1. + poissonP_ );
  alpha2_ = alpha_ * alpha_;
  Ka_     = K_ * alpha_;
  nuPFac_ = std::sqrt ( 1. / ( 1. + 2. * poissonP_ * poissonP_ ) );

  // Get hardening curves (tension and compression).

  String args = "x";
  props.find ( args, "args" );

  sigmaT_ = FuncUtils::newFunc ( args, SIGMA_T_PROP, props, globdat );
  sigmaC_ = FuncUtils::newFunc ( args, SIGMA_C_PROP, props, globdat );

  FuncUtils::resolve           ( *sigmaT_, globdat );
  FuncUtils::resolve           ( *sigmaC_, globdat );
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::getConfig
//-----------------------------------------------------------------------


void MelroMaterial::YieldFunc_::getConfig

  ( const Properties&  conf ) const

{
  using jive::util::FuncUtils;


  conf.set ( RM_TOL_PROP         , rmTol_    );
  conf.set ( PLASTIC_POISSON_PROP, poissonP_ );
  conf.set ( RM_MAX_ITER_PROP    , maxIter_  );

  FuncUtils::getConfig ( conf, sigmaT_, SIGMA_T_PROP );
  FuncUtils::getConfig ( conf, sigmaC_, SIGMA_C_PROP );
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::isPlastic
//-----------------------------------------------------------------------


bool MelroMaterial::YieldFunc_::isPlastic

  ( const Invariants&  inv,
    double             epspeq0 )

{
  // Update all variables that are constant inside the return 
  // mapping scheme.

  epspeq0_ = epspeq0;

  j2tr_  = inv.getJ2 ();
  i1tr_  = inv.getFirstInvariant ();
  j2fac_ = 18.*j2tr_;
  i1fac_ = 4.*alpha2_*i1tr_*i1tr_ / 27.;

  // A = 18 J2tr / zetaS^2 + (4/27) * (a * I1tr / zetaP)^2 

  // Evaluate failure criterion for the trial state.

  const double crit = evalTrial_    ();

  const double sigT = sigmaT_->eval ( epspeq0_ );
  const double sigC = sigmaC_->eval ( epspeq0_ );

  return (crit/sigC/sigT) > -rmTol_;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::findRoot
//-----------------------------------------------------------------------


double MelroMaterial::YieldFunc_::findRoot 

  ( double  dgam0 )

{
  double dgam     = dgam0;
  double oldddgam = -1.;
  double oldcrit  = evalTrial_ ();

  for ( idx_t irm = 0; irm < maxIter_; )
  {
    double crit = eval_ ( dgam );

    double normalized = jem::numeric::abs( crit ) / ( sigC_*sigT_ );

    if ( normalized < rmTol_ )
    {
      // System::out() << " return mapping converged in " << irm 
      // << " iterations " << normalized << endl;
      break;
    }
    if ( ++irm == maxIter_ )
    {
      throw jem::Exception ( JEM_FUNC, "no convergence" );
    }
    if ( jem::Float::isNaN( dcdg_ ) )
    {
      throw jem::Exception ( JEM_FUNC, "nan" );
    }
    // System::out() << " return mapping iteration " << irm << " dgam " <<
    // dgam << ", criterion: " << crit << ", dcdg " << dcdg << endl;

    double ddgam = crit / dcdg_;

    if ( ddgam * oldddgam < 0. )
    {
      // Divergence detection.

      if ( oldcrit * crit < 0.  && 
           jem::numeric::abs( ddgam ) > jem::numeric::abs( oldddgam ) )
      {
        // There might be an inflection point around the root
        // use linear interpolation rather than linearization.

        // System::out() << "divergence prevention" << endl;
        ddgam = -oldddgam * crit / ( crit - oldcrit );
      }
    }

    dgam -= ddgam;

    oldddgam = ddgam;
    oldcrit  = crit;
  }

  if ( dgam < -jem::Limits<double>::EPSILON )
  {
    throw jem::Exception ( JEM_FUNC, "negative dgam" );
  }

  return dgam;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::findBounds
//-----------------------------------------------------------------------


void  MelroMaterial::YieldFunc_::findBounds

  ( double&  gmin,
    double&  gmax,
    double&  fmin,
    double&  fmax )

{
  gmin = 0.;
  gmax = -1.;
  fmax = -1.;

  fmin = evalTrial_ ();

  JEM_ASSERT ( fmin > 0. );

  double dgam = jem::Limits<double>::EPSILON;

  while ( gmax < 0. )
  {
    double crit = eval_ ( dgam );

    if ( crit > 0. )
    {
      gmin = dgam;
      fmin = crit;
    }
    else
    {
      gmax = dgam;
      fmax = crit;
    }
    dgam *= 10.;
  }
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::findRoot (with bounds)
//-----------------------------------------------------------------------


double MelroMaterial::YieldFunc_::findRoot

  ( double  gmin,
    double  gmax,
    double  fmin,
    double  fmax )

{
  double dgam;

  try 
  {
    double dgam0 = estimateRoot_ ( gmin, gmax, fmin, fmax );

    dgam = findRoot ( dgam0 );
  }
  catch ( const jem::Exception& ex )
  {
    improveBounds_  ( gmin, gmax, fmin, fmax );

    dgam = findRoot ( gmin, gmax, fmin, fmax );
  }

  return dgam;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::getTangentParameters
//-----------------------------------------------------------------------
 

void MelroMaterial::YieldFunc_::getTangentParameters 

  ( double&  betam, 
    double&  pbm, 
    double&  hfac, 
    double&  beta, 
    double&  phi, 
    double&  rho, 
    double&  chi, 
    double&  psi, 
    double&  xi, 
    double&  omega, 
    double   dgam )

{
  setDGam_ ( dgam );

  betam = 6. * G_ / zetaS_;
  pbm = 2.*Ka_/3./zetaP_ - betam/3.;

  if ( jem::numeric::abs( depspeq_ ) > 1.e-20 )
  {
    hfac = 2. * ( ( HC_ - HT_ ) * i1tr_ / zetaP_ - (sigC_*HT_+sigT_*HC_) ) 
                * dgam*dgam*nuPFac_*nuPFac_/depspeq_;
  }
  else
  {
    hfac =  0.;
  }

  double eta  = -dcdg_;

  beta        = 2. * G_ / zetaS_;
  phi         = K_ / zetaP_ * ( 1. - 4. * kaizp2_ * sigct_ / eta );
  rho         = 36. * K_ * G_ * sigct_ / eta / zS2_ / zetaP_;
  chi         = 72. * G_ * G_ / eta / zS2_ / zS2_;
  psi         = 8. * G_ * kaizp2_  / eta / zS2_;
  xi          = 2. * kaizp2_ / eta / 3.;
  omega       = 6. * G_ / zS2_ / eta;
}


void MelroMaterial::YieldFunc_::getTangentParameters

  ( double&  beta, 
    double&  phi, 
    double&  rho, 
    double&  chi, 
    double&  psi )

{
  // Fast version of getTangentParameters for dgam = 0.

  resetDGam_ ();

  double eta = -dcdg_;

  beta = 2. * G_;
  phi  = K_ * ( 1. - 4. * kaizp2_ * sigct_ / eta );
  rho = 36. * K_ * G_ * sigct_ / eta;
  chi = 72. * G_ * G_ / eta;
  psi = 8. * G_ * kaizp2_ / eta;
}


void MelroMaterial::YieldFunc_::getPGradientParameters

  ( double&  phiF,
    double&  rhoF,
    double&  chiF,
    double&  psiF ) const

{
  double eta = -dcdg_;

  phiF = 4. * kaizp2_ * sigct_ / 3. / eta;
  rhoF = 18. * K_ * sigct_ / eta;
  chiF = 36. * G_ / eta;
  psiF = 8. * G_ * alpha_ * i1tr_ / 3. / eta;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::getAlpha 
//-----------------------------------------------------------------------


inline double MelroMaterial::YieldFunc_::getAlpha

  () const noexcept

{
  return alpha_;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::getZetaP
//-----------------------------------------------------------------------


inline double MelroMaterial::YieldFunc_::getZetaP

  () const noexcept

{
  return zetaP_;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::getZetaS
//-----------------------------------------------------------------------


inline double MelroMaterial::YieldFunc_::getZetaS

  () const noexcept

{
  return zetaS_;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::getEpspeq
//-----------------------------------------------------------------------


inline double MelroMaterial::YieldFunc_::getEpspeq

  () const noexcept

{
  return epspeq_;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::evalTrial_
//-----------------------------------------------------------------------


double MelroMaterial::YieldFunc_::evalTrial_ 

  () const

{
  // Fast version of eval_ function, without computing all 'dgam' 
  // dependent variables.

  double sigC = sigmaC_->eval ( epspeq0_ );
  double sigT = sigmaT_->eval ( epspeq0_ );
  double sigct = sigC - sigT;

  return 6.*j2tr_ + 2.*i1tr_*sigct - 2.*sigC*sigT;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::eval_
//-----------------------------------------------------------------------


double MelroMaterial::YieldFunc_::eval_

  ( double  dgam )

{
  setDGam_ ( dgam );

  return 6.*j2tr_/zS2_ + 2.*i1tr_*sigct_/zetaP_ - 2.*sigC_*sigT_;
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::setDGam_
//-----------------------------------------------------------------------


void MelroMaterial::YieldFunc_::setDGam_

  ( double  dgam ) 

{
  // Update all variables that are a function of 'dgam'.

  zetaS_ = 1. + 6. * G_  * dgam;
  zetaP_ = 1. + 2. * Ka_ * dgam;

  zS2_   = zetaS_*zetaS_;
  zP2_   = zetaP_*zetaP_;

  const double j2fac = 18.*j2tr_;
  const double i1fac = 4.*alpha2_*i1tr_*i1tr_ / 27.;
  const double sqrtA = std::sqrt ( j2fac / zS2_ + i1fac / zP2_ );

  depspeq_ = nuPFac_ * dgam * sqrtA;
  epspeq_  = epspeq0_ + depspeq_;

  sigC_    = sigmaC_->eval ( epspeq_ );
  sigT_    = sigmaT_->eval ( epspeq_ );
  sigct_   = sigC_ - sigT_;

  HC_      = sigmaC_->deriv ( epspeq_ );
  HT_      = sigmaT_->deriv ( epspeq_ );

  gjzs3_   = G_  * j2tr_ / zS2_ / zetaS_;
  kaizp2_  = Ka_ * i1tr_ / zP2_;

  // NOTE: adding an additional I1tr in the final term of depedg
  //       (Melro's paper is incorrect here).

  // Eq. (14) in Frans' article.

  double depedg = nuPFac_ * ( sqrtA - .5 * dgam / sqrtA * 
      ( 216.*gjzs3_ + 16.*alpha2_*kaizp2_*i1tr_/27./zetaP_ ) );

  double dsigCe = HC_ * depedg;
  double dsigTe = HT_ * depedg;

  // Eq. (12) in Frans' article.

  dcdg_ = 2.*i1tr_*(dsigCe-dsigTe)/zetaP_ - 4.*kaizp2_*sigct_
          - 72.*gjzs3_ - 2.*(sigC_*dsigTe+sigT_*dsigCe);
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::resetDGam_
//-----------------------------------------------------------------------


void MelroMaterial::YieldFunc_::resetDGam_ 

  ()

{
  // Update all variables that are a function of 'dgam', using dgam = 0.

  zetaS_ = 1.;
  zetaP_ = 1.;

  zS2_   = 1.;
  zP2_   = 1.;

  const double j2fac = 18.*j2tr_;
  const double i1fac = 4.*alpha2_*i1tr_*i1tr_ / 27.;
  const double sqrtA = std::sqrt ( j2fac + i1fac );

  depspeq_ = 0.;
  epspeq_  = epspeq0_;

  sigC_    = sigmaC_->eval ( epspeq_ );
  sigT_    = sigmaT_->eval ( epspeq_ );
  sigct_   = sigC_ - sigT_;

  HC_      = sigmaC_->deriv ( epspeq_ );
  HT_      = sigmaT_->deriv ( epspeq_ );

  gjzs3_   = G_  * j2tr_;
  kaizp2_  = Ka_ * i1tr_;

  double depedg = nuPFac_ * sqrtA;

  double dsigCe = HC_ * depedg;
  double dsigTe = HT_ * depedg;

  dcdg_ = 2.*i1tr_*(dsigCe-dsigTe) - 4.*kaizp2_*sigct_
          - 72.*gjzs3_ - 2.*(sigC_*dsigTe+sigT_*dsigCe);
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::estimateRoot_
//-----------------------------------------------------------------------


double MelroMaterial::YieldFunc_::estimateRoot_

  ( double  gmin,
    double  gmax,
    double  fmin,
    double  fmax ) const

{
  // linear interpolation

  // fmin  o___ 
  //           ---___
  // phi=0 |----------o-----------|------
  //                      ---___
  // fmax  -                    --o
  //     gmin        ret        gmax

  return gmin + fmin * (gmax-gmin) / (fmin-fmax);
}


//-----------------------------------------------------------------------
//   MelroMaterial::YieldFunc_::improveBounds_
//-----------------------------------------------------------------------


void  MelroMaterial::YieldFunc_::improveBounds_

  ( double&  gmin,
    double&  gmax,
    double&  fmin,
    double&  fmax )

{
  const idx_t  np = 10_idx;
  const double dg = ( gmax - gmin ) / double( np );
  double dgam     = gmin;

  while ( dgam < gmax )
  {
    dgam += dg;
    double crit = eval_ ( dgam );

    if ( crit > 0. )
    { 
      gmin = dgam;
      fmin = crit;
    }
    else
    {
      gmax = dgam;
      fmax = crit;

      break;
    }
  }
}


//=======================================================================
//   class MelroMaterial 
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  MelroMaterial::TYPE_NAME            = "Melro";

const char*  MelroMaterial::PLASTIC_POISSON_PROP = "poissonP";
const char*  MelroMaterial::RM_TOL_PROP          = "rmTolerance";
const char*  MelroMaterial::RM_MAX_ITER_PROP     = "rmMaxIter";
const char*  MelroMaterial::SIGMA_T_PROP         = "sigmaT";
const char*  MelroMaterial::SIGMA_C_PROP         = "sigmaC";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


MelroMaterial::MelroMaterial 

  ( const String&      name, 
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

    : Super ( name, globdat, 3_idx )

{
  // Get material Properties.

  Properties myProps = props.getProps ( myName_ );
  Properties myConf  = conf.makeProps ( myName_ );

  // Get "dim" (or "rank" of Melro material) parameter.

  melroRank_ = -1_idx;

  myProps.get ( melroRank_, RANK_PROP, 2_idx, 3_idx ); 
  myConf .set ( RANK_PROP , melroRank_);
}


MelroMaterial::~MelroMaterial 

  ()

{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void MelroMaterial::configure 

  ( const Properties&  props,
    const Properties&  globdat )

{
  // First, let Super material handle its own configure action.

  Super::configure ( props, globdat );

  // Then, MelroMaterial deals with its onw configure member function.

  Properties myProps = props.getProps ( myName_ );

  if ( melroRank_== 2_idx )
  {
    myProps.get( stateString_, STATE_PROP );

    if      ( stateString_ == "PlaneStrain" )
    {
      state_ = ProblemType::PlaneStrain;
    }
    else if ( stateString_ == "PlaneStress" )
    {
      throw jem::IllegalInputException ( 
        CLASS_NAME, 
        "plasticity for Plane Stress is not implemented yet!" 
      );
    }

    // Compute the 3x3 matrix for 2D problems.

    comp2DStiffMat_ ();
  }

  if ( melroRank_ == 3_idx )
  {
    state_ = ProblemType::Solid;
     
    stateString_ = "Solid";
  }

  // Set a YieldFunc_ member

  y_ = jem::newInstance<YieldFunc_> ( young(), poisson() );

  // ... and configure it.

  y_->configure ( myProps, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void MelroMaterial::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const 

{
  Super::getConfig ( conf, globdat );

  Properties myConf = conf.makeProps ( myName_ );

  y_->getConfig ( myConf );
}


//-----------------------------------------------------------------------
//   commit 
//-----------------------------------------------------------------------


void MelroMaterial::commit 

  () 

{
  newHist_.swap ( preHist_ );

  latestHist_ = & preHist_;
}


//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------


void MelroMaterial::allocPoints

  ( idx_t  count )

{
  for ( idx_t i = 0; i < count; i++ )
  {
    preHist_.pushBack ( Hist_() );
    newHist_.pushBack ( Hist_() );
  }

  latestHist_ = & preHist_;
}


//-----------------------------------------------------------------------
//   rank
//-----------------------------------------------------------------------


idx_t MelroMaterial::rank

  () const noexcept

{
  return melroRank_;
}


//-----------------------------------------------------------------------
//   getProblemType
//-----------------------------------------------------------------------


ProblemType MelroMaterial::getProblemType 

  () const noexcept

{
  return state_;
}


//-----------------------------------------------------------------------
//   elasticUpdate
//-----------------------------------------------------------------------


void MelroMaterial::elasticUpdate

  ( const Vector&  stress,
    const Matrix&  stiff,
    const Vector&  strain ) const

{
  if ( melroRank_ == 3_idx )
  {
    Super::elasticUpdate ( stress, stiff, strain );
  }
  else
  {
    stiff = stiff2D_;

    jem::numeric::matmul ( stress, stiff, strain );
  }
}


//-----------------------------------------------------------------------
//   elasticUpdate
//-----------------------------------------------------------------------


Vector MelroMaterial::elasticUpdate

  ( const Vector&  strain ) const

{
  if ( melroRank_ == 3_idx )
  {
    return Super::elasticUpdate ( strain );
  }
  else
  {
    return jem::numeric::matmul ( stiff2D_, strain );
  }
}


//-----------------------------------------------------------------------
//   update 
//-----------------------------------------------------------------------


void MelroMaterial::update 

  ( const Vector&  stress,
    const Matrix&  stiff,
    const Vector&  strain,
    idx_t          ip ) const 

{
  Ref<Self> self = const_cast<Self*> ( this ); // remove constness for
                                               // updt. history-related
                                               // member variables


  m66       dmat;
  Vec6      eps   ( fill3DStrain( strain ) ); 
  Vec6      sig;
  Vec6      sigtr;

  double    dgam;
  double    epspeq0 = preHist_[ip].epspeq;
  bool      loading = false;

  sigtr = v6_to_t6_ ( 
      Super::elasticUpdate ( 
          t6_to_v6_ ( std::move( Vec6( eps - preHist_[ip].epsp ) ) 
        ) 
      ) 
    ); // elastic stress field {sig} = [D] {eps - epsp}

  StressInvariants inv ( sigtr );

  if ( y_->isPlastic( inv, epspeq0 ) )
  {
    try
    {
      // Classical return mapping scheme starting with dgam = 0.

      dgam = y_->findRoot ( 0. );
    }
    catch ( const jem::Exception& ex )
    {
      handleException_ ( ex, strain, ip );

      // More robust return mapping scheme.

      double gmin, gmax, fmin, fmax;

      y_->findBounds      ( gmin, gmax, fmin, fmax );

      dgam = y_->findRoot ( gmin, gmax, fmin, fmax );
    }

    // Compute stress.

    const double ptr    = inv.getFirstInvariant() / 3.;
    const double p      = ptr    / y_->getZetaP ();

    Vec6         sigDtr = deviatoric_ ( sigtr, ptr );
    Vec6         sigD   = sigDtr / y_->getZetaS();

    sig     = sigD;
    sig[0] += p;
    sig[1] += p;
    sig[2] += p;
    
    // Compute plastic strain increment. 

    const double alphaP23 = y_->getAlpha() * p * 2. / 3.;

    Vec6 mvec;

    mvec[0]  = 3. * sigD[0] + alphaP23;
    mvec[1]  = 3. * sigD[1] + alphaP23;
    mvec[2]  = 3. * sigD[2] + alphaP23;
    mvec[3]  = 6. * sigD[3];  // engng strain factor 2
    mvec[4]  = 6. * sigD[4];  
    mvec[5]  = 6. * sigD[5]; 

    Vec6 depsp = dgam * mvec;

    // Consistent tangent.

    double betam, pbm, hfac, beta, phi, rho, chi, psi, xi, omega;

    y_->getTangentParameters 
      ( betam, pbm, hfac, beta, phi, rho, chi, psi, xi, omega, dgam );

    m66  dmde = aI4_plus_bII_ ( betam, pbm );

    Vec6 hvec = hfac * jem::numeric::matmul ( dmde, mvec );

    dmat = aI4_plus_bII_ ( beta, phi-beta/3. );

    for ( idx_t i = 0; i < 6; ++i )
    {
      for ( idx_t j = 0; j < 6; ++j )
      {
        dmat(i,j) -= chi   * sigDtr[i] * sigDtr[j];
        dmat(i,j) -= omega * sigDtr[i] * hvec[j];
      }
      for ( idx_t j = 0; j < 3; ++j )
      {
        dmat(i,j) -= rho * sigDtr[i];
        dmat(j,i) -= psi * sigDtr[i];
        dmat(j,i) -= xi  * hvec[i];
      }
    }

    // Update history.

    double dG = jem::dot ( sig, depsp );

    self->newHist_[ip].epsp        = preHist_[ip].epsp + depsp;
    self->newHist_[ip].epspeq      = y_->getEpspeq ();
    self->newHist_[ip].dissipation = preHist_[ip].dissipation + dG;

    loading = true;
  }
  else
  {
    Matrix stiff3D ( 6, 6 );

    stiff3D = Super::getStiffMat ();
    sig     = sigtr;

    for ( idx_t i = 0; i < 6; i++ )
    {
      for ( idx_t j = 0; j < 6; j++ )
      {
        dmat(i,j) = stiff3D(i,j);
      }
    }

    self->newHist_[ip].epsp        = preHist_[ip].epsp;
    self->newHist_[ip].epspeq      = preHist_[ip].epspeq;
    self->newHist_[ip].dissipation = preHist_[ip].dissipation;
  }

  self->newHist_[ip].loading = loading;
  self->latestHist_          = & (self->newHist_);

  reduce3DVector_ ( stress, sig  );
  reduce3DMatrix_ ( stiff , dmat );
}


void MelroMaterial::update 

  ( const Vector&  stress,
    const Vector&  strain,
    idx_t          ip ) const 

{
  // NOTE: this member function is suposed to be used in functions that 
  //       write data for post-processing tasks at a converged state.
  //       Therefore, here, it does not execute a return mapping scheme
  //       again; it just computes stresses as:
  //       {sig} = [D] {eps - epsp} 'elastic' strain/stress field(s)

  // Compute elastic strain field.

  Vec6 eps    ( fill3DStrain( strain ) ); // total   strain

  Vec6 epsp = (*latestHist_)[ip].epsp;    // plastic strain

  Vec6 epse = eps - epsp;                 // elastic strain

  // Get 3D elastic stiffness matrix.

  Matrix stiff3D ( 6, 6 );
  stiff3D = Super::getStiffMat ();

  // Compute stress field.

  m66 tstiff;

  for ( idx_t i = 0; i < 6; i++ )
  {
    for ( idx_t j = 0; j < 6; j++ )
    {
      tstiff(i,j) = stiff3D(i,j);
    }
  }

  Vec6 sig = jem::numeric::matmul ( tstiff, epse );

  reduce3DVector_ ( stress, sig );
}


//-----------------------------------------------------------------------
//   getStiffMat
//-----------------------------------------------------------------------


Matrix MelroMaterial::getStiffMat           

  () const 

{
  if ( melroRank_ == 3_idx )
  {
    return Super::getStiffMat ();
  }
  else
  {
    return stiff2D_;
  }
}


//-----------------------------------------------------------------------
//   makeNew 
//-----------------------------------------------------------------------

 
Ref<Material> MelroMaterial::makeNew

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


void MelroMaterial::declare

  ()

{
  MaterialFactory::declare ( TYPE_NAME,  & makeNew );
  MaterialFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   comp2DStiffMat_
//-----------------------------------------------------------------------


void MelroMaterial::comp2DStiffMat_ 

  ()

{
  const idx_t  n  = FEMUtils::STRAIN_COUNTS[melroRank_];

  const double E  = young   ();
  const double nu = poisson ();

  stiff2D_.resize ( n, n );
  stiff2D_ = 0.;

  if      ( melroRank_ == 2_idx && state_ == ProblemType::PlaneStrain )
  {
    const double a = E / ((1. + nu) * (1. - 2. * nu));
    const double b = .5 * (1. - 2. * nu);
    const double c = 1. - nu;

    stiff2D_ (0,0) = a * c;
    stiff2D_ (0,1) = stiff2D_ (1,0) = a * nu;
    stiff2D_ (1,1) = a * c;
    stiff2D_ (2,2) = a * b;
  }
  else if ( melroRank_ == 2_idx && state_ == ProblemType::PlaneStress )
  {
    const double a = E / (1. - nu * nu);

    stiff2D_ (0,0) = a;
    stiff2D_ (0,1) = stiff2D_ (1,0) = a * nu;
    stiff2D_ (1,1) = a;
    stiff2D_ (2,2) = a * .5 * (1. - nu);
  }
  else
  {
    throw jem::IllegalInputException ( 
      CLASS_NAME, "unexpected 2D material for local rank!" 
    );
  }
}


//-----------------------------------------------------------------------
//   handleException_
//-----------------------------------------------------------------------


void MelroMaterial::handleException_ 

  ( const jem::Exception&  ex,
    const Vector&          strain,
    idx_t                  ip ) const

{
  using jem::io::endl;
  using jem::io::PrintWriter;


  Ref<PrintWriter> out = jem::newInstance<PrintWriter>
                                     ( & jem::System::debug( "melro" ) );

  out->nformat.setFractionDigits ( 10 );

  if      ( ex.what().equals( "no convergence" ) )
  {
    *out << "No convergence in return mapping algorithm " << endl;
  }
  else if ( ex.what().equals( "nan" ) )
  {
    *out << "NaN detected in return mapping algorithm " << endl;
  }
  else if ( ex.what().equals( "negative dgam" ) )
  {
    *out << "Negative increment found" << endl;
  }
  else
  {
    *out << "Caught unkown exception: " << ex.what() << endl;
  }

  *out << "epsp    " << preHist_[ip].epsp   << endl;
  *out << "epspeq0 " << preHist_[ip].epspeq << endl;
  *out << "strain  " << strain << endl << endl;
}


//-----------------------------------------------------------------------
//   aI4_plus_bII_
//-----------------------------------------------------------------------


Tuple<double,6,6> MelroMaterial::aI4_plus_bII_ 

  ( double  a,
    double  b ) const

{
  // Make 6x6 matrix with a*I_4^s+b*II in Voigt notation.

  m66 ret;
  ret = 0.;

  ret(0,0) += a+b;  ret(0,1) += b;    ret(0,2) += b;
  ret(1,0) += b;    ret(1,1) += a+b;  ret(1,2) += b;
  ret(2,0) += b;    ret(2,1) += b;    ret(2,2) += a+b;
  ret(3,3) += .5*a;
  ret(4,4) += .5*a;
  ret(5,5) += .5*a;

  return ret;
}


//-----------------------------------------------------------------------
//   reduce3DVector_
//-----------------------------------------------------------------------

 
void MelroMaterial::reduce3DVector_
    
  ( const Vector&  v,
    const Vec6&    t ) const

{
  // Reduce a full 3D tuple to a 2D or 3D vector (works the same for both 
  // stress and strain)

  const idx_t n = v.size ();

  if       ( n == 3_idx )
  {
    v[0] = t[0];
    v[1] = t[1];
    v[2] = t[3];
  }
  else if ( n == 4_idx )
  {
    Tuple<idx_t,4> aximap ( 0, 1, 3, 2 );

    for ( idx_t i = 0; i < 4; ++i )
    {
      v[i] = t[aximap[i]];
    }
  }
  else
  {
    v[0] = t[0];
    v[1] = t[1];
    v[2] = t[2];
    v[3] = t[3];
    v[4] = t[4];
    v[5] = t[5];
  }
}


//-----------------------------------------------------------------------
//   reduce3DMatrix_
//-----------------------------------------------------------------------


void MelroMaterial::reduce3DMatrix_

  ( const Matrix&  m,
    const m66&     t ) const

{
  using jem::numeric::LUSolver;


  const idx_t n = m.size ( 0 );


  if      ( n == 3_idx )
  {
    if ( state_ == ProblemType::PlaneStrain )
    {
      select2DMatrix_ ( m, t );
    }
    else
    {
      double d;

      m66    tmp66 = t;

      LUSolver::invert ( tmp66, d     );
      select2DMatrix_  ( m    , tmp66 );
      LUSolver::invert ( m    , d     );
    }
  }
  else if ( n == 4_idx )
  {
    selectAxiMatrix_ ( m, t );
  }
  else
  {
    for ( idx_t i = 0; i < 6; i++ )
    {
      for ( idx_t j = 0; j < 6; j++ )
      {
        m(i,j) = t(i,j);
      }
    }
  }
}


//-----------------------------------------------------------------------
//   select2DMatrix_
//-----------------------------------------------------------------------


void MelroMaterial::select2DMatrix_

  ( const Matrix&  m,
    const m66&     t ) const

{
  // Selecting [xx, yy, xy] components from [xx, yy, zz, xy, yz, xz].
 
  for ( idx_t i = 0; i < 2; i++ )
  {
    m(i,2) = t(i,3_idx);
    m(2,i) = t(3_idx,i);
 
    for ( idx_t j = 0; j < 2; j++ )
    {
      m(i,j) = t(i,j);
    }
  }

  m(2_idx,2_idx) = t(3_idx,3_idx);
}


//-----------------------------------------------------------------------
//   selectAxiMatrix_ 
//-----------------------------------------------------------------------


void MelroMaterial::selectAxiMatrix_

  ( const Matrix&  m,
    const m66&     t ) const

{
  // Selecting [xx, yy, xy, zz] components from [xx, yy, zz, xy, yz, xz].
 
  Tuple<idx_t,4> aximap ( 0, 1, 3, 2 );
 
  for ( idx_t i = 0; i < 4; ++i )
  {
    for ( idx_t j = 0; j < 4; ++j )
    {
      m(i,j) = t ( aximap[i], aximap[j] );
    }
  }
}

