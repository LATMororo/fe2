
/*
 *  Model that prints load-displacement data for nodegroup to file:
 *  The sum of the nodal loads and the average of the nodal displacements
 *
 *  Adapted from: MonitorModel, 
 *  FPM, 29-1-2008
 * 
 *  Data is computed in GET_MATRIX0 and GET_INT_VECTOR
 *         and optionally written to file in COMMIT
 *         (rather using the SampleModule to write to file recommended)
 *
 *  Changelog:
 *
 *  Modified: Luiz Ant. T. Mororo, l.a.taumaturgomororo@tudelft.nl
 *  Date:     July, 2019
 *
 *  The class has been changed in order to enable parallel computing.
 *
 *  Modified: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date: Fev 2025
 *
 *  Clean-up code, and refactored member functions for better
 *  performance, and to meet new Jem/Jive (version 3.0) data types, for 
 *  instance, `idx_t'. 
 */


#ifndef LOAD_DISP_MODEL_H 
#define LOAD_DISP_MODEL_H


#include <jive/fem/NodeGroup.h>
#include <jive/model/Model.h>
#include <jive/mp/VectorExchanger.h>
#include <jive/util/Assignable.h>
#include <jive/util/DofSpace.h>


using jem::idx_t;
using jem::Ref;
using jem::String;


using jem::util::Properties;


using jive::IdxMatrix;
using jive::IdxVector;
using jive::MPContext;
using jive::StringVector;
using jive::Vector;


using jive::fem::NodeGroup;
using jive::fem::NodeSet;
using jive::model::Model;
using jive::mp::VectorExchanger;
using jive::util::Assignable;
using jive::util::DofSpace;


//=======================================================================
//   class LoadDispModel
//=======================================================================


class LoadDispModel : public Model
{
 public:

  JEM_DECLARE_CLASS       ( LoadDispModel, Model );

  static const char*        TYPE_NAME; 

  static const char*        NODES_PROP;
  static const char*        GROUP_PROP;
  static const char*        FILE_NAME_PROP;
  static const char*        TYPES_PROP;

                            LoadDispModel

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

  virtual                  ~LoadDispModel         ();


 private:

  void                      updateIDofs_          ();

  void                      gatherData_    
  
    ( const Vector&           disp,
      const Vector&           load,
      const Vector&           state,
      const Vector&           fint,
      const Properties&       globdat )              const;

 private:

  Ref<DofSpace>             dofs_;
  Assignable<NodeSet>       nodes_;
  Assignable<NodeGroup>     ngroup_;
  Ref<VectorExchanger>      vex_;
  Ref<MPContext>            mpx_;

  IdxMatrix                 idofs_;

  IdxVector                 inodes_;
  IdxVector                 idofborders_;
  IdxVector                 types_;

  StringVector              typeNames_;

  String                    groupName_;

  idx_t                     nn_;
  idx_t                     nnOld_;
  idx_t                     rank_;
  idx_t                     typeCount_;

  bool                      isParallel_;
};

#endif


