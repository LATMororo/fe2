
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This class implements a material factory for creation of new 
 *  continuum materials (constitutive models).
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date: Jan 2025
 */


#include <jive/util/CtorMap.h>

#include <jem/base/System.h>


#include "MaterialFactory.h"
#include "Material.h"


//=======================================================================
//   class MaterialFactory::CtorMap_
//=======================================================================


class MaterialFactory::CtorMap_ : public 
                                        jive::util::CtorMap<Constructor>

{};


//=======================================================================
//   class MaterialFactory 
//=======================================================================


//-----------------------------------------------------------------------
//   declare
//-----------------------------------------------------------------------


void MaterialFactory::declare

  ( const String&  type,
    Constructor    ctor )

{
  JEM_PRECHECK ( ctor );

  CtorMap_::insert ( type, ctor );
}


//-----------------------------------------------------------------------
//   exists
//-----------------------------------------------------------------------


bool MaterialFactory::exists 

  ( const String&  type )

{
  return CtorMap_::contains ( type );
}


//-----------------------------------------------------------------------
//   listKnownTypes
//-----------------------------------------------------------------------


StringVector MaterialFactory::listKnownTypes 

  ()

{
  return CtorMap_::listKnownTypes ();
}


//-----------------------------------------------------------------------
//   newInstance
//-----------------------------------------------------------------------


Ref<Material> MaterialFactory::newInstance

  ( const String&      type,
    const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  Constructor ctor = CtorMap_::find ( type );

  if ( ctor )
  {
    return ctor ( name, conf, props, globdat );
  }
  else
  {
    return nullptr;
  }
}


Ref<Material> MaterialFactory::newInstance

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  Properties  subConf  = conf.makeProps ( name );
  Properties  subProps = props.getProps ( name );

  Ref<Material>  mat;
  String         type;

  // TYPE_PROP = "type" is defined in jive/util/Factory.h.

  subProps.get ( type, TYPE_PROP ); 
  subConf .set ( TYPE_PROP, type );

  mat = newInstance ( type, name, conf, props, globdat );

  if ( ! mat )
  {
    noSuchTypeError (
      subProps.getContext ( TYPE_PROP ),
      type, "material",
      CtorMap_::listKnownTypes ()
    );
  }

  return mat;
}

