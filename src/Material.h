
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This abstract class implements an interface for continuum mechanics 
 *  constitutive models.
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date: Jan 2025
 */


#ifndef MATERIAL_H
#define MATERIAL_H


#include <jem/base/NamedObject.h>
#include <jem/util/Properties.h>


#include "FEMUtils.h"


using jem::Ref;
using jem::String;


using jem::util::Properties;


using FEMUtils::ProblemType;


//=======================================================================
//   class Material 
//=======================================================================


class Material : public jem::NamedObject
{
 public:

  JEM_DECLARE_CLASS       ( Material, NamedObject );

  static const char*        RANK_PROP;

  typedef Array<double,1>   Vector;
  typedef Array<double,2>   Matrix;

  
  virtual String            getContext            () const override;

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig 

    ( const Properties&       conf,
      const Properties&       globdat )              const; 

  virtual void              commit                ();

  virtual void              cancel                ();

  virtual void              allocPoints

    ( idx_t                   count );

  virtual idx_t             rank                  () const noexcept = 0;
  virtual ProblemType       getProblemType        () const noexcept = 0;

  virtual void              elasticUpdate

    ( const Vector&           stress,
      const Matrix&           stiff,
      const Vector&           strain )               const = 0;

  virtual Vector            elasticUpdate

    ( const Vector&           strain )               const = 0;

  virtual void              update 

    ( const Vector&           stress,
      const Matrix&           stiff,
      const Vector&           strain,
      idx_t                   ip )                   const = 0;

  virtual void              update 

    ( const Vector&           stress,
      const Vector&           strain,
      idx_t                   ip )                   const = 0;

  virtual Matrix            getStiffMat           () const = 0;

  virtual Tuple<double,6>   fill3DStrain

    ( const Vector&           v3 )                   const = 0;

  virtual Tuple<double,6>   fill3DStress

    ( const Vector&           v3 )                   const = 0;


 protected:

  explicit                  Material

    ( const String&           name = "" );

  virtual                  ~Material              ();

};


//=======================================================================
//   related functions
//=======================================================================


Ref<Material>               newMaterial

  ( const String&             name,
    const Properties&         conf,
    const Properties&         props,
    const Properties&         globdat );


#endif
