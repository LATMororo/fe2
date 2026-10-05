
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This class implements a material factory for creation of new 
 *  continuum materials (constitutive models).
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date: Jan 2025
 */


#ifndef MATERIAL_FACTORY_H
#define MATERIAL_FACTORY_H


#include <jem/util/Properties.h>


#include <jive/util/Factory.h>


using jem::idx_t;
using jem::Ref;
using jem::String;


using jem::util::Properties;


using jive::StringVector;


//=======================================================================
//   forward declaration
//=======================================================================


class Material;


//=======================================================================
//   class MaterialFactory 
//=======================================================================


class MaterialFactory : public jive::util::Factory
{
 public:

  typedef MaterialFactory   Self; 
  typedef Factory           Super;


  typedef Ref<Material>    (*Constructor)

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare

    ( const String&           type,
      const Constructor       ctor );

  static bool               exists

    ( const String&           type );

  static StringVector       listKnownTypes        ();

  static Ref<Material>      newInstance

    ( const String&           type, 
      const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static Ref<Material>      newInstance

    ( const String&           name, 
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

 private:

  class                     CtorMap_;

};


#endif
