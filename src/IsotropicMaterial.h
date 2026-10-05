
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This concrete class implements a material based on Hooke's law for 
 *  2D (plane stress and plane strain) and 3D problems. 
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date: Jan 2025
 */


#ifndef ISOTROPIC_MATERIAL_H
#define ISOTROPIC_MATERIAL_H


#include "Material.h"


//=======================================================================
//   class IsotropicMaterial
//=======================================================================


class IsotropicMaterial : public Material
{
 public:

  JEM_DECLARE_CLASS       ( IsotropicMaterial, Material );

  static const char*        TYPE_NAME;

  static const char*        YOUNG_PROP;
  static const char*        POISSON_PROP;
  static const char*        AREA_PROP;
  static const char*        RHO_PROP;
  static const char*        STATE_PROP;


  explicit                  IsotropicMaterial

    ( const String&           name,
      const Properties&       globdat,
      idx_t                   rank );

  explicit                  IsotropicMaterial

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

  virtual Tuple<double,6>   fill3DStrain

    ( const Vector&           v3 )                   const override;

  virtual Tuple<double,6>   fill3DStress

    ( const Vector&           v3 )                   const override;

  virtual double            young                 () const noexcept;
  virtual double            poisson               () const noexcept;

  static Ref<Material>      makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare               ();


 protected:

  virtual                  ~IsotropicMaterial     ();


 private:

  void                      computeStiffMat_      ();


 private:

  Matrix                    stiff_;

  String                    stateString_;

  double                    young_;
  double                    poisson_;
  double                    area_;
  double                    rho_;

  idx_t                     rank_;

  ProblemType               state_;

};


#endif
