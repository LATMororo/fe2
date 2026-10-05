
/*
 *  Copyright (C) 2014 TU Delft. All rights reserved.
 *  
 *  Frans van der Meer, September 2014
 *  
 *  Module to generate NodeGroups and ElementGroups from input file
 *  used as child Module in GmshInputModule
 *
 *  Changelog:
 *
 *  Modified: Luiz Ant. T. Mororo, l.a.taumaturgomororo@tudelft.nl
 *  Date:     July 2020
 * 
 *  The class was changed in order to enable parallel computing.
 *
 *  Modified: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date:     Fev 2025
 * 
 *  Clean-up code, remove parallel functionalities, and refactored member 
 *  functions for better performance, storage scheme, and to meet new 
 *  Jem/Jive (version 3.0) data types, for instance, idx_t. 
 *
 *  NGroup and EGroup are no longer members of this class. They are 
 *  actually helper classes that are defined and implemented in the
 *  corresponding *.cpp file.
 *
 *  Besides, functions that do not store NodeGroup and ElementGroup are
 *  designed. They just find nodes and elements.
 */


#ifndef GROUP_INPUT_MODULE_H
#define GROUP_INPUT_MODULE_H


#include <jive/app/Module.h>
#include <jive/fem/ElementSet.h>


using jem::idx_t;
using jem::Ref;
using jem::String;


using jem::util::Properties;


using jive::app::Module;
using jive::fem::ElementSet;
using jive::fem::NodeSet;


//=======================================================================
//   class GroupInputModule 
//=======================================================================


class GroupInputModule : public Module
{
 public:

  JEM_DECLARE_CLASS       ( GroupInputModule, Module );

  static const char*        TYPE_NAME;

  static const char*        NODE_GRPS_PROP;
  static const char*        ELEM_GRPS_PROP;


  explicit                  GroupInputModule 

    ( const String&           name = "groupInput" );

  virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat )                    override;

  jive::IdxVector           findNodes

    ( const NodeSet&          nodes,
      const Properties&       conf,
      const Properties&       props,
      const String&           name )                 const;

  jive::IdxVector           findElems

    ( const ElementSet&       elems,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat, 
      const String&           name )                 const;
      
  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare               ();
 

 protected:

  virtual                  ~GroupInputModule      ();

};


#endif
