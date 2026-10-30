
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This abstract class implements an interface for continuum mechanics (CM)
 *  constitutive models.
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date: Jan 2025
 */


#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/IllegalOperationException.h>
#include <jem/util/None.h>


#include "Material.h"
#include "MaterialFactory.h"


JEM_DEFINE_CLASS ( Material );


//=======================================================================
//   class Material
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  Material::RANK_PROP = "rank"; 


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


Material::Material

  ( const String&  name ) : Super ( name )

{}


Material::~Material

  ()

{}


//-----------------------------------------------------------------------
//   getContext 
//-----------------------------------------------------------------------


String Material::getContext 

  () const

{
  return NamedObject::makeContext ( "material", myName_ );
}


//-----------------------------------------------------------------------
//   configure 
//-----------------------------------------------------------------------


void Material::configure

  ( const Properties&  props,
    const Properties&  globdat )

{}


//-----------------------------------------------------------------------
//   getConfig 
//-----------------------------------------------------------------------


void Material::getConfig 

  ( const Properties&  conf,
    const Properties&  globdat ) const 

{}


//-----------------------------------------------------------------------
//   commit 
//-----------------------------------------------------------------------


void Material::commit

  ()

{}


//-----------------------------------------------------------------------
//   cancel 
//-----------------------------------------------------------------------


void Material::cancel

  ()

{}


//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------


void Material::allocPoints 

  ( idx_t  count )

{}


//-----------------------------------------------------------------------
//   historyCount
//-----------------------------------------------------------------------


idx_t Material::historyCount

  () const

{
  return -1;
}


//-----------------------------------------------------------------------
//   getHistory & setHistory
//-----------------------------------------------------------------------


void Material::getHistory

  ( const Vector&  hist,
    idx_t          ip ) const

{
  throw jem::IllegalOperationException (
    getContext (),
    "material does not expose its history variables"
  );
}


void Material::setHistory

  ( const Vector&  hist,
    idx_t          ip )

{
  throw jem::IllegalOperationException (
    getContext (),
    "material does not expose its history variables"
  );
}


//=======================================================================
//   related functions
//=======================================================================


Ref<Material> newMaterial

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  Ref<jem::Object> obj;

  props.get ( obj, name );

  if ( jem::util::isNone( obj ) ) 
  {
    throw jem::IllegalInputException (  
      JEM_FUNC, "No Material `material' property has been defined!" 
    );
  }

  // Get the correct material object.

  return MaterialFactory::newInstance ( name, conf, props, globdat );
}
