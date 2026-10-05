
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This concrete class implements a material based on Hooke's law for 
 *  2D (plane stress and plane strain) and 3D problems. 
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date: Jan 2025
 */


#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/numeric/algebra/matmul.h>


#include "IsotropicMaterial.h"
#include "MaterialFactory.h"


JEM_DEFINE_CLASS ( IsotropicMaterial );


using namespace jem::literals;


using FEMUtils::STRAIN_COUNTS;


typedef Array<double,1>     Vector;
typedef Array<double,2>     Matrix;


//=======================================================================
//   class IsotropicMaterial 
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  IsotropicMaterial::TYPE_NAME    = "Isotropic";

const char*  IsotropicMaterial::YOUNG_PROP   = "young";
const char*  IsotropicMaterial::POISSON_PROP = "poisson";
const char*  IsotropicMaterial::AREA_PROP    = "area";
const char*  IsotropicMaterial::RHO_PROP     = "rho";
const char*  IsotropicMaterial::STATE_PROP   = "state";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


IsotropicMaterial::IsotropicMaterial

  ( const String&      name,
    const Properties&  globdat,
    idx_t              rank ) : Super ( name ), rank_ ( rank )

{
  JEM_ASSERT ( rank >= 1_idx && rank <= 3_idx );

  // Set (default values) and initialize some members.

  young_   = 1.0;
  poisson_ = 1.0;
  rho_     = 1.0;

  stiff_  .resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] ); 
  stiff_  = 0.0;
}


IsotropicMaterial::IsotropicMaterial

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

  : Super ( name )

{
  // Get material Properties.

  Properties myProps = props.getProps ( myName_ );
  Properties myConf  = conf.makeProps ( myName_ );

  // Get "dim" (or "rank" of material) parameter.

  rank_ = -1_idx;

  myProps.get ( rank_, RANK_PROP, 1_idx, 3_idx ); 
  myConf .set ( RANK_PROP, rank_ );

  // Set (default values) and initialize some members.

  young_   = 1.0;
  poisson_ = 1.0;
  rho_     = 1.0;

  stiff_  .resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] ); 
  stiff_  = 0.0;
}


IsotropicMaterial::~IsotropicMaterial

  ()

{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void IsotropicMaterial::configure 

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jem::maxOf;


  Properties myProps = props.getProps ( myName_ );

  myProps.get  ( young_,   YOUNG_PROP,   0.0, maxOf( young_ ) );
  myProps.get  ( poisson_, POISSON_PROP, 0.0, 0.5             );

  myProps.find ( rho_, RHO_PROP, 0.0, maxOf( rho_ ) );

  // NOTE: the stateString_ and state_ could be handled in a different 
  //       way. Perhaps, one could use HashMap, i.e., 
  //       map(ProblemType, stateString). 

  if ( rank_ == 1_idx )
  {
    myProps.get ( area_, AREA_PROP );

    stateString_ = "Bar";
  }

  if ( rank_ == 2_idx )
  {
    myProps.get( stateString_, STATE_PROP );

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

  // Compute the elastic moduli, only once time.

  computeStiffMat_ ();
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void IsotropicMaterial::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const 

{
  Properties myConf = conf.makeProps ( myName_ );

  myConf.set ( YOUNG_PROP  , young_       );
  myConf.set ( POISSON_PROP, poisson_     );
  myConf.set ( STATE_PROP  , stateString_ );

  if ( rank_ == 1_idx )
  {
    myConf.set ( AREA_PROP, area_ );
  }
}


//-----------------------------------------------------------------------
//   rank
//-----------------------------------------------------------------------


idx_t IsotropicMaterial::rank

  () const noexcept

{
  return rank_;
}


//-----------------------------------------------------------------------
//   getProblemType
//-----------------------------------------------------------------------


ProblemType IsotropicMaterial::getProblemType 

  () const noexcept

{
  return state_;
}


//-----------------------------------------------------------------------
//   elasticUpdate
//-----------------------------------------------------------------------


void IsotropicMaterial::elasticUpdate

  ( const Vector&  stress,
    const Matrix&  stiff,
    const Vector&  strain ) const

{
  // Compute the elastic moduli.

  stiff = stiff_;   

  // Compute the stress vector.

  jem::numeric::matmul ( stress, stiff_ , strain );
}


//-----------------------------------------------------------------------
//   elasticUpdate
//-----------------------------------------------------------------------


Vector IsotropicMaterial::elasticUpdate

  ( const Vector&  strain ) const

{
  return jem::numeric::matmul ( stiff_, strain );
}


//-----------------------------------------------------------------------
//   update 
//-----------------------------------------------------------------------


void IsotropicMaterial::update 

  ( const Vector&  stress,
    const Matrix&  stiff,
    const Vector&  strain,
    idx_t          ip ) const 

{
  elasticUpdate ( stress, stiff, strain );
}


void IsotropicMaterial::update 

  ( const Vector&  stress,
    const Vector&  strain,
    idx_t          ip ) const 

{
  stress = elasticUpdate ( strain );
}


//-----------------------------------------------------------------------
//   getStiffMat
//-----------------------------------------------------------------------


Matrix IsotropicMaterial::getStiffMat           

  () const 

{
  return stiff_;
}


//-----------------------------------------------------------------------
//   fill3DStress
//-----------------------------------------------------------------------


Tuple<double,6> IsotropicMaterial::fill3DStress

  ( const Vector&  v3 ) const

{
  JEM_ASSERT ( v3.size() >= 1 && v3.size() <= 3 );

  if ( v3.size() == 3 )
  {
    double sig_zz = state_ == ProblemType::PlaneStress
                  ? 0.
                  : poisson_ * ( v3[0] + v3[1] );

    return FEMUtils::fillFrom2D ( v3, sig_zz );
  }
  else
  {
    return FEMUtils::fillFrom3D ( v3 );
  }
}


//-----------------------------------------------------------------------
//   fill3DStrain
//-----------------------------------------------------------------------


Tuple<double,6> IsotropicMaterial::fill3DStrain

  ( const Vector&  v3 ) const

{
  JEM_ASSERT ( v3.size() >= 1_idx && v3.size() <= 6_idx );

  if ( v3.size() == 3_idx )
  {
    double eps_zz = state_ == ProblemType::PlaneStress
                  ? -poisson_ / (1.-poisson_) * (v3[0]+v3[1])
                  : 0.;

    return FEMUtils::fillFrom2D ( v3, eps_zz );
  }
  else
  {
    return FEMUtils::fillFrom3D ( v3 );
  }
}


//-----------------------------------------------------------------------
//   young
//-----------------------------------------------------------------------


double IsotropicMaterial::young 

  () const noexcept

{
  return young_;
}


//-----------------------------------------------------------------------
//   poisson
//-----------------------------------------------------------------------


double IsotropicMaterial::poisson

  () const noexcept

{
  return poisson_;
}


//-----------------------------------------------------------------------
//   makeNew 
//-----------------------------------------------------------------------

 
Ref<Material> IsotropicMaterial::makeNew

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


void IsotropicMaterial::declare

  ()

{
  MaterialFactory::declare ( TYPE_NAME,  & makeNew );
  MaterialFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   computeStiffMat_
//-----------------------------------------------------------------------


void IsotropicMaterial::computeStiffMat_

  () 

{
  const int     n  = STRAIN_COUNTS[rank_];

  const double  e  = young_;
  const double  nu = poisson_;


  if      ( rank_ == 1 )
  {
    stiff_ (0,0) = e * area_;
  }
  else if ( rank_ == 3 )
  {
    const double  a = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double  b = 0.5 * (1.0 - 2.0 * nu);
    const double  c = 1.0 - nu;

    stiff_ (0,0) = a * c;
    stiff_ (0,1) = stiff_ (0,2) = a * nu;
    stiff_ (1,1) = stiff_ (0,0);
    stiff_ (1,2) = stiff_ (0,1);
    stiff_ (2,2) = stiff_ (0,0);
    stiff_ (3,3) = a * b;
    stiff_ (4,4) = stiff_ (3,3);
    stiff_ (5,5) = stiff_ (3,3);

    // Copy lower triangle of the stress-strain matrix.

    for ( int i = 0; i < n; i++ )
    {
      for ( int j = 0; j < i; j++ )
      {
        stiff_ (i,j) = stiff_ (j,i);
      }
    }
  }
  else if ( state_ == ProblemType::PlaneStrain )
  {
    const double  a = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double  b = 0.5 * (1.0 - 2.0 * nu);
    const double  c = 1.0 - nu;

    stiff_ (0,0) = a * c;
    stiff_ (0,1) = stiff_ (1,0) = a * nu;
    stiff_ (1,1) = a * c;
    stiff_ (2,2) = a * b;
  }
  else if ( state_ == ProblemType::PlaneStress )
  {
    const double  a = e / (1.0 - nu * nu);

    stiff_ (0,0) = a;
    stiff_ (0,1) = stiff_ (1,0) = a * nu;
    stiff_ (1,1) = a;
    stiff_ (2,2) = a * 0.5 * (1.0 - nu);
  }
  else
  {
    // TODO: a better way could be designed to handle problem types.

    throw jem::IllegalInputException ( 
      CLASS_NAME, 
      "unexpected problem type: `" + stateString_ + 
      "'; available types: `Bar', `PlaneStrain', `PlaneStress', `Solid'" 
    );
  }
}
