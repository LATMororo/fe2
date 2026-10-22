
#include <cstring>


#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>


#include "PythonMaterialBridge.h"
#include "PythonMaterialBindings.h"


namespace py = pybind11;

using Input  = py::array_t<double, py::array::c_style|py::array::forcecast>;


namespace {


//==================================================================================================
//   Helper functions
//==================================================================================================


static std::vector<double> values_

  ( const Input&  array )

{
  return std::vector<double>(array.data(),array.data()+array.size());
}

static py::array_t<double> array_ 

  ( const std::vector<double>&          data,
    std::initializer_list<py::ssize_t>  shape )

{
  py::array_t<double> result ( shape );

  if ( !data.empty() ) 
  {
    std::memcpy ( result.mutable_data(), data.data(), data.size()*sizeof( double ) );
  }

  return result;
}


} // End of unnamed namespace


//==================================================================================================
//   bindMaterials functions
//==================================================================================================

static void bindMaterials 

  ( py::module_&  module )

{
  module.doc () = "Generic Jem/Jive material factory, trial updates and history derivatives";

  // Use `std::shared_ptr<>' as holder type for PythonMaterialBatch and PythonMaterialFactory.
  // Usage example: `from _fe2_material_native import MaterialFactory'.

  // PythonMaterialBatch.

  py::class_<PythonMaterialBatch, std::shared_ptr<PythonMaterialBatch>> ( module, "MaterialBatch" )
    // Binding for the function 'count': getter function
    .def_property_readonly ( "count",
      [] ( const PythonMaterialBatch&  m ) { 
        return static_cast<int>( m.count() ); 
      } 
    )
    // Binding for the function 'components': getter function 
    .def_property_readonly ( "components",
      [] ( const PythonMaterialBatch&  m ) {
        return static_cast<int>( m.components() );
      }
    )
    // Binding for the function 'update': update material => internal variables, e.g., plasticity
    .def ( "update", 
      [] ( PythonMaterialBatch&  mat, const Input&  strain ) {
        if ( strain.ndim()  != 2           || 
            strain.shape(0) != mat.count() || 
            strain.shape(1) != mat.components() )
        {
          throw py::value_error ( "Expected strain [material points, components]!" );
        }

        MaterialResponse response = mat.update ( values_( strain) );

        return py::make_tuple ( 
          array_( response.stress, { mat.count(), mat.components() } ),
          array_( response.stiff,  { mat.count(), mat.components(), mat.components() } )
        );
    })
    // Bindings for the functions 'commit' and 'cancel': they handle internal variables
    .def( "commit", &PythonMaterialBatch::commit )
    .def( "cancel", &PythonMaterialBatch::cancel )
    // Bindings the vector-Jacobian product
    .def( "history_vjp",
      [] ( PythonMaterialBatch&  mat, 
           const Input&          history, 
           const Input&          strain,
           const Input&          adjoint, 
           double                step ) {

        if ( history.ndim()   != 3           || 
             history.shape(1) != mat.count() || 
             history.shape(2) != mat.components() )
        {
          throw py::value_error ( "Expected history [steps, material points, components]!" );
        }

        std::vector<double> result = mat.historyVjp ( 
          values_( history ), 
          history.shape(0), 
          values_( strain ), 
          values_( adjoint ), 
          step 
        );

        return array_ ( result, { history.shape(0), mat.count(), mat.components() } );
    }, py::arg( "history" ), py::arg( "strain" ), py::arg( "adjoint" ), py::arg( "step" ) = 1.e-7 );

  // MaterialFactory.
 
  py::class_
    <PythonMaterialFactory, std::shared_ptr<PythonMaterialFactory>> ( module, "MaterialFactory" )
    // Binding for Python constructor __init__
    .def( py::init<const std::string&>(), py::arg( "props" ) )
    // Binding for __call__ dunder method for creating 'material points'
    .def( "__call__", [] ( const PythonMaterialFactory& factory, int count ) { 
        return factory.create( count ); 
      } 
    )
    // Binding for the function type: list of known material types of Jive
    .def_static( "types", &PythonMaterialFactory::types );
} 


#ifdef FE2_PYTHON_EXTENSION
PYBIND11_MODULE ( _fe2_material_native, module )
{
  // The public Python wrapper ensures this library is loaded on the main thread: Jem's static 
  // initializers allocate before this function runs. Bootstrap the executor before callbacks from 
  // temporary Python workers.
   
  PythonMaterialFactory::types ();
  bindMaterials                ( module );
}
#else
PYBIND11_EMBEDDED_MODULE ( _fe2_material, module) 
{ 
  bindMaterials ( module ); 
}
#endif
