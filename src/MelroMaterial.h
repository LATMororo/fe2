
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


#ifndef MELRO_MATERIAL_H 
#define MELRO_MATERIAL_H 


#include <jem/util/Flex.h>


#include "IsotropicMaterial.h"
#include "Invariants.h"


//=======================================================================
//   class MelroMaterial 
//=======================================================================


class MelroMaterial : public IsotropicMaterial 
{
 public:

  JEM_DECLARE_CLASS       ( MelroMaterial, IsotropicMaterial );

  static const char*        TYPE_NAME;

  static const char*        PLASTIC_POISSON_PROP;
  static const char*        RM_TOL_PROP;
  static const char*        RM_MAX_ITER_PROP;
  static const char*        SIGMA_T_PROP;
  static const char*        SIGMA_C_PROP;

  typedef Tuple<double,6,6> m66;


  explicit                  MelroMaterial 

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat )                    override;

  virtual void              getConfig 

    ( const Properties&       conf,
      const Properties&       globdat )              const override; 

  virtual void              commit                ()       override;

  virtual void              allocPoints

    ( idx_t                   count )                      override;

  virtual idx_t             rank                  () const noexcept
                                                           override;

  virtual ProblemType       getProblemType        () const noexcept 
                                                           override;
  virtual void              elasticUpdate

    ( const Vector&           stress,
      const Matrix&           stiff,
      const Vector&           strain )               const override;

  virtual Vector            elasticUpdate

    ( const Vector&           strain )               const override;

  virtual void              update 

    ( const Vector&           stress,
      const Matrix&           stiff,
      const Vector&           strain,
      idx_t                   ip )                   const override;

  virtual void              update 

    ( const Vector&           stress,
      const Vector&           strain,
      idx_t                   ip )                   const override;

  virtual Matrix            getStiffMat           () const override;

  static Ref<Material>      makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare               ();

  
 protected:

  virtual                  ~MelroMaterial         ();
 

 private:

  class                     YieldFunc_;

  class                     Hist_
  {
   public:

                            Hist_                 ();

    void                    toVector 
    
      ( const                 Vector& vec )          const;


    Vec6                    epsp;    // plastic strain

    double                  epspeq;  // equivalent plastic strain
    double                  dissipation;
    bool                    loading;

  };

  void                      comp2DStiffMat_       (); 

  void                      handleException_ 

    ( const jem::Exception&   ex,
      const Vector&           strain,
      idx_t                   ip )                   const;

  m66                       aI4_plus_bII_ 

    ( double                  a,
      double                  b )                    const;

  void                      reduce3DVector_
    
    ( const Vector&           v,
      const Vec6&             t )                    const; 

  void                      reduce3DMatrix_
    
    ( const Matrix&           m,
      const m66&              t )                    const;

  void                      select2DMatrix_
    
    ( const Matrix&           m,
      const m66&              t )                    const;

  void                      selectAxiMatrix_
    
    ( const Matrix&           m,
      const m66&              t )                    const;

  inline void               t6_to_v6_

    ( const Vector&           v6,
      const Vec6&             t6 )                   const;

  inline Vector             t6_to_v6_

    ( const Vec6&             t6 )                   const;

  inline void               v6_to_t6_

    (       Vec6&             t6,
      const Vector&           v6 )                   const;

  inline Vec6               v6_to_t6_

    ( const Vector&           v6 )                   const;

  inline Vec6               deviatoric_

    ( const Vec6&             full,
      double                  p )                    const;
 
  
 private:

  Ref<YieldFunc_>           y_;

  jem::util::Flex<Hist_>    preHist_;
  jem::util::Flex<Hist_>    newHist_;
  jem::util::Flex<Hist_>*   latestHist_;

  Matrix                    stiff2D_;

  String                    stateString_;

  idx_t                     melroRank_;

  ProblemType               state_;

};


//#######################################################################
//   Implementation
//#######################################################################


//-----------------------------------------------------------------------
//   t6_to_v6_ 
//-----------------------------------------------------------------------


inline void MelroMaterial::t6_to_v6_

  ( const Vector&  v6,
    const Vec6&    t6 ) const

{
  v6[0] = t6[0]; v6[1] = t6[1]; v6[2] = t6[2]; 
  v6[3] = t6[3]; v6[4] = t6[4]; v6[5] = t6[5]; 
}


inline Vector MelroMaterial::t6_to_v6_ 

  ( const Vec6&  t6 ) const

{
  Vector v6  ( 6 );
  t6_to_v6_  ( v6, t6 );

  return v6;
}


//-----------------------------------------------------------------------
//   v6_to_t6_ 
//-----------------------------------------------------------------------


inline void MelroMaterial::v6_to_t6_

  (       Vec6&    t6,
    const Vector&  v6 ) const

{
  t6[0] = v6[0]; t6[1] = v6[1]; t6[2] = v6[2]; 
  t6[3] = v6[3]; t6[4] = v6[4]; t6[5] = v6[5]; 
}


inline Vec6 MelroMaterial::v6_to_t6_

  ( const Vector&  v6 ) const

{
  Vec6 t6;
  v6_to_t6_ ( t6, v6 );

  return t6;
}


//-----------------------------------------------------------------------
//   deviatoric_ 
//-----------------------------------------------------------------------


Vec6 MelroMaterial::deviatoric_

  ( const Vec6&  full,
    double       p ) const

{
  Vec6 ret = full;

  ret[0] -= p;
  ret[1] -= p;
  ret[2] -= p;

  return ret;
}


#endif
