
#ifndef PYTHON_MATERIAL_BRIDGE_H
#define PYTHON_MATERIAL_BRIDGE_H


#include <string>
#include <vector>
#include <memory> // smart pointers


#include "Material.h"


//==================================================================================================
//   forward declarations
//==================================================================================================


class PythonMaterialFactory;


//==================================================================================================
//   Material response
//==================================================================================================

// Plain C++ data across the Python boundary; Jem objects remain on a Jem-managed executor thread, 
// including when called by PyTorch workers.


struct MaterialResponse
{
  std::vector<double>                     stress;
  std::vector<double>                     stiff;

};


// Gradients of one material step (strain_t, history_{t-1}) -> (stress_t, history_t).

struct StepGradients
{
  std::vector<double>                     strain;
  std::vector<double>                     history;

};


//==================================================================================================
//   class PythonMaterialBatch
//==================================================================================================


class PythonMaterialBatch 
{
 public:
                                         ~PythonMaterialBatch   ();

  idx_t                                   count                 () const noexcept;
  idx_t                                   components            () const noexcept;

  // Number of history variables per point; negative if the material does not expose them.

  idx_t                                   historyCount          () const noexcept;

  MaterialResponse                        update 

    ( const std::vector<double>&            strain );

  void                                    commit                ();
  void                                    cancel                ();

  // History variables of all points after the latest update.

  std::vector<double>                     history               ();

  // Vector-Jacobian product of a single step, given the history before that step.

  StepGradients                           stepVjp 

    ( const std::vector<double>&            history,
      const std::vector<double>&            strain,
      const std::vector<double>&            stressAdjoint,
      const std::vector<double>&            historyAdjoint,
      double                                step );

  // Computes the gradient of stress of given earlier strain fields by replaying them. Fallback for
  // materials that do not expose their history.

  std::vector<double>                     historyVjp 

    ( const std::vector<double>&            history,
      idx_t                                 steps,
      const std::vector<double>&            strain,
      const std::vector<double>&            adjoint, 
      double                                step );


 private:

  friend class                            PythonMaterialFactory;
  class                                   State_;

  explicit                                PythonMaterialBatch 

    ( State_*                               state, 
      idx_t                                 count, 
      idx_t                                 components,
      idx_t                                 historyCount );

  State_*                                 state_;
  idx_t                                   count_;
  idx_t                                   components_;
  idx_t                                   historyCount_;

};


//==================================================================================================
//   class PythonMaterialFactory 
//==================================================================================================


class PythonMaterialFactory
{
 public:

  explicit                                PythonMaterialFactory  

    ( const std::string&                    props );

                                          PythonMaterialFactory 

    ( const Properties&                     matProps,
      const Properties&                     globdat );

  std::shared_ptr<PythonMaterialBatch>    create 

    ( idx_t                                 count )                 const;
  
  static std::vector<std::string>         types                 ();

                                         ~PythonMaterialFactory () = default;

 private:

  // Class to handle existing Jive materials, i.e., materials that have been implemented in pure 
  // Jive.

  class                                   Specification_; 
  std::shared_ptr<Specification_>         spec_;

};


#endif
