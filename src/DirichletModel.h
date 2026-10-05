
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This class implements a model for setting Dirichlet boundary 
 *  conditions (BCs) for certain group of nodes. 
 *  
 *  It inherits and follows what has been implemented F. P. van der Meer
 *  his collegues at TU Delft. Clean-up, refactored member functions in 
 *  order to meet new Jem/Jive (version 3.0) data types (e.g., `idx_t')
 *  have been accomplished.
 *
 *  However, the state increment is accomplished by means of a function,
 *  which can also be defined by the user. In this regard, such function
 *  can be written in term of coordinates (x, y, and z - they depend on
 *  the dimensionality of the problem), TIME_STEP (i) and TIME (t);
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date:   April 2026
 */


#ifndef DIRICHLET_MODEL_H
#define DIRICHLET_MODEL_H


#include <jive/model/Model.h>
#include <jive/fem/NodeSet.h>
#include <jive/util/Assignable.h>


using jem::idx_t;
using jem::Ref;
using jem::String;


using jem::util::Properties;
using jem::numeric::Function;


using jive::model::Model;
using jive::fem::NodeSet;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::DofSpace;


//=======================================================================
//   class DirichletModel
//=======================================================================


class DirichletModel : public Model
{
 public:

  JEM_DECLARE_CLASS       ( DirichletModel, Model );

  static const char*        TYPE_NAME; 

  static const char*        NODE_GROUPS_PROP;
  static const char*        DOFS_PROP;
  static const char*        STATE_FUNCS_PROP;
  static const char*        LOADED_PROP;
  static const char*        STATE_INCR_PROP;
  static const char*        INIT_STATE_PROP;


                            DirichletModel

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

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat )                    override;

  static Ref<Model>         makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare               ();


 protected:

  virtual                  ~DirichletModel        ();

 
 private:

  void                      init_

    ( const Properties&       globdat );

  void                      advance_

    ( const Properties&       globdat );  

  void                      applyConstraints_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      commit_

    ( const Properties&       params,
      const Properties&       globdat );
      
  void                      setStepSize_

    ( const Properties&       params );      


 private:

  Assignable<NodeSet>       nodes_;
  Ref<DofSpace>             dofs_;
  Ref<Constraints>          cons_;

  jem::Array
       < Ref<Function> >    stateFuncs_;

  jive::StringVector        nodeGroups_;
  jive::StringVector        dofTypes_;

  idx_t                     ngroups_;

};


#endif
