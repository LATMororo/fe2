
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


//==================================================================================================
//   class PythonMaterialBatch
//==================================================================================================


class PythonMaterialBatch 
{
 public:
                                         ~PythonMaterialBatch   ();

  idx_t                                   count                 () const noexcept;
  idx_t                                   components            () const noexcept;

  MaterialResponse                        update 

    ( const std::vector<double>&            strain );

  void                                    commit                ();
  void                                    cancel                ();

  // Computes the gradient of stress of given earlier strain fields. 

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
      idx_t                                 components );

  State_*                                 state_;
  idx_t                                   count_;
  idx_t                                   components_;

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
