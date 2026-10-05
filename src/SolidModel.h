
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  Simple class for stress analysis. It has been intended to be used 
 *  with continuum finite element meshes. This implementation has been 
 *  adapted from the one given in JIVE's tutorial and course by Erik Jan 
 *  Lingen (Dynaflow Reaserch Group).
 *
 *  The main tasks of such model are:
 *
 *  Assemble internal force vector and stiffness matrix necessary for 
 *  linear and nonlinear analyses, and post-processing results after
 *  COMMIT action.
 *
 *  NOTE: columns of tables for post-processing operations are defined
 *        here; however, each material could define its own set of 
 *        columns. For instance, materials with plasticity might define
 *        a column for plastic strain field.
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date:   Fev 2025
 *  
 *  Changelog:
 *
 *  Modified: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br, luiztaumaturgo@gmail.com
 *  Date: Jul 2026
 *
 *  Compute average strains and stresses in the elements.
 *  
 */


#ifndef SOLID_MODEL_H
#define SOLID_MODEL_H


#include <jem/mt/WorkPool.h>
#include <jem/util/ObjFlex.h>


#include <jive/fem/ElementSet.h>
#include <jive/geom/InternalShape.h>
#include <jive/model/Model.h>
#include <jive/util/Assignable.h>


#include "Material.h"
#include "FEMUtils.h"


using jem::util::Properties;


using jive::algebra::MatrixBuilder;
using jive::fem::ElementSet;
using jive::fem::NodeSet;
using jive::geom::IShape;
using jive::model::Model;
using jive::util::Assignable;
using jive::util::DofSpace;
using jive::util::XTable;


//=======================================================================
//   class SolidModel
//=======================================================================


class SolidModel : public Model
{
 public:

  JEM_DECLARE_CLASS       ( SolidModel, Model );

  static const char*        TYPE_NAME; 

  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;
  static const char*        THICKNESS_PROP;
  static const char*        THREADS_PROP;


                            SolidModel 
                         
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

  virtual                  ~SolidModel            ();
 

 private:

  void                      setIPMap_             () const;

  idx_t                     getIPIdx_

    ( idx_t                   ielem,
      idx_t                   ip )                   const;

  void                      assemble_

    ( MatrixBuilder&          mbld,
      const Vector&           force,
      const Vector&           disp )                 const;

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat )              const;

  void                      setStrainTable_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           disp )                 const;

  void                      setStressTable_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           disp )                 const;

  void                      setElemStrainTable_

    ( XTable&                 table,
      const Vector&           disp )                 const;

  void                      setElemStressTable_

    ( XTable&                 table,
      const Vector&           disp )                 const;


 private:
  
  class                     Worker_;       // friend class to help in  
  friend class              Worker_;       // assemb. global sys. of eq.
                                           // concurrently
  Assignable<ElementSet>    elems_;
  Assignable<NodeSet>       nodes_;

  Ref<jem::mt::WorkPool>    pool_;
  jem::util::ObjFlex        jobs_;
  jem::util::ObjFlex        workers_;

  Ref<IShape>               shape_;
  Ref<DofSpace>             dofs_;

  Ref<jem::numeric::Function>
                            thickFunc_;    // user-defined functions for 
                                           // computing thickness (rank < 3)

  Ref<Material>             material_;

  IdxVector                 ielems_;       // element indices
  IdxVector                 dofTypes_;     // dof type indices

  FEMUtils::ShapeGradsFunc  getShapeGrads; // shape function gradients

  idx_t                     rank_;         // dimension of the mesh/problem
  idx_t                     numElem_;      // nr of elements in the mesh
  idx_t                     strCount_;     // nr of strain components

  idx_t                     nodeCount_;    // nodes per element
  idx_t                     ipCount_;      // integration point per element
  idx_t                     dofCount_;     // nr of dofs per element

};


#endif
