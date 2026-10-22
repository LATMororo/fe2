
#include <mutex>
#include <vector>


#include <pybind11/embed.h>
#include <pybind11/numpy.h>


#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/IllegalOperationException.h>


#include "PRNNMaterial.h"
#include "PythonMaterialBridge.h"
#include "MaterialFactory.h"


JEM_DEFINE_CLASS ( PRNNMaterial );


namespace py = pybind11;


using namespace jem::literals;


namespace
{

//==================================================================================================
//   Helper functions
//==================================================================================================


std::string text_ 

  ( const String&  value ) 

{
  return std::string ( value.addr(), value.size() );
}


void initializePython_ 

  ()

{
  static std::once_flag once;

  std::call_once(once, [] {
    if ( !Py_IsInitialized() )
    {
      PyConfig config;
      PyConfig_InitPythonConfig ( &config );

      PyStatus status = PyConfig_SetBytesString ( 
        &config, &config.program_name, FE2_PYTHON_EXECUTABLE 
      );

      if ( PyStatus_Exception( status ) ) 
      {
        throw std::runtime_error ( "Cannot configure embedded Python executable!" );
      }

      py::initialize_interpreter ( &config );
      PyEval_SaveThread          ();
    }
  });
  // Interpreter lives for the host process. Never finalize PyTorch underneath another material or 
  // a running Jem worker.
}

py::array_t<double> input_ 

  ( const Vector&  strain )
{
  py::array_t<double> result ( static_cast<py::ssize_t>( strain.size() ) );

  for ( idx_t i = 0; i < strain.size(); i++ ) result.mutable_data()[i] = strain[i];

  return result;
}


}


//==================================================================================================
//   class PRNNMaterial::Impl_
//==================================================================================================


class __attribute__((visibility("hidden"))) PRNNMaterial::Impl_
{
 public:

  std::shared_ptr<PythonMaterialFactory>  factory;

  py::object                              loadModel;
  py::object                              runtime;

  std::vector<py::object>                 models;

  Properties                              materialProps;

  String                                  module     = "pyprnn.checkpoint"; 
  String                                  loader     = "load_model";
  String                                  pythonPath = "";
  String                                  stateName;
  String                                  weights;

  idx_t                                   rank       = 0_idx;
  idx_t                                   components = 0_idx;

  ProblemType                             state      = ProblemType::PlaneStrain;

  bool                                    configured = false;

};


//==================================================================================================
//   class PRNNMaterial
//==================================================================================================


//--------------------------------------------------------------------------------------------------
//   static data
//--------------------------------------------------------------------------------------------------


const char* PRNNMaterial::TYPE_NAME        = "PRNN";

const char* PRNNMaterial::WEIGHTS_PROP     = "weights";
const char* PRNNMaterial::MODULE_PROP      = "module";
const char* PRNNMaterial::LOADER_PROP      = "loader";
const char* PRNNMaterial::PYTHON_PATH_PROP = "pythonPath";
const char* PRNNMaterial::MATERIAL_PROP    = "material";
const char* PRNNMaterial::STATE_PROP       = "state";


//--------------------------------------------------------------------------------------------------
//   constructor & destructor
//--------------------------------------------------------------------------------------------------


PRNNMaterial::PRNNMaterial 

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) : Super( name ), impl_( new Impl_ )

{
  // Read and set material rank

  props.getProps( name ).get( impl_->rank, RANK_PROP, idx_t( 2 ), idx_t( 3 ) );

  conf.makeProps( name ).set( RANK_PROP, impl_->rank );

  // ... andi set strain/stress components based on material rank.

  impl_->components = FEMUtils::STRAIN_COUNTS[ impl_->rank ];
}


PRNNMaterial::~PRNNMaterial 

  ()

{
  if ( Py_IsInitialized() ) 
  {
    py::gil_scoped_acquire gil;
    impl_.reset ();
  }
}


//--------------------------------------------------------------------------------------------------
//   configure
//--------------------------------------------------------------------------------------------------


void PRNNMaterial::configure 

  ( const Properties&  props,
    const Properties&  globdat )

{
  if ( !impl_->models.empty() )
  {
    throw jem::IllegalOperationException (
      getContext(), "Cannot reconfigure allocated Python PRNN history!"
    );
  }
  
  Properties myProps = props.getProps ( myName_ );

  myProps.get  ( impl_->weights   , WEIGHTS_PROP     );
  myProps.find ( impl_->module    , MODULE_PROP      );
  myProps.find ( impl_->loader    , LOADER_PROP      );
  myProps.find ( impl_->pythonPath, PYTHON_PATH_PROP );

  // Read material state.

  impl_->stateName = "Solid";

  if ( impl_->rank == 2_idx ) 
  {
    myProps.get ( impl_->stateName, STATE_PROP );
  }

  if      ( impl_->rank == 2_idx && impl_->stateName == "PlaneStress" ) 
  {
    impl_->state = ProblemType::PlaneStress;
  }
  else if ( impl_->rank == 2_idx && impl_->stateName == "PlaneStrain" ) 
  {
    impl_->state = ProblemType::PlaneStrain;
  }
  else if ( impl_->rank == 3_idx && impl_->stateName == "Solid" ) 
  {
    impl_->state = ProblemType::Solid;
  }
  else 
  {
    throw jem::IllegalInputException ( getContext(), "Invalid macro rank/state!" );
  }

  impl_->materialProps = myProps.getProps ( MATERIAL_PROP );
  impl_->factory       = std::make_shared<PythonMaterialFactory> ( impl_->materialProps, globdat );

  initializePython_ ();

  py::gil_scoped_acquire gil;
  try
  {
    auto path = py::module_::import( "sys" ).attr( "path" );
    path.attr( "insert" )( 0, FE2_PYTHON_SOURCE );

    if ( impl_->pythonPath.size() ) 
    {
      path.attr( "insert" )( 0, text_( impl_->pythonPath ) );
    }

    py::module_::import ( "_fe2_material" );

    impl_->runtime    = py::module_::import ( "pyprnn.runtime" );

    impl_->loadModel  = py::module_::import ( 
      text_ ( impl_->module ).c_str() ).attr( text_( impl_->loader ).c_str() 
    );

    impl_->configured = true;
  }
  catch ( const py::error_already_set&  error ) 
  {
    throw jem::IllegalInputException ( getContext(), String( error.what() ) );
  }
}


//--------------------------------------------------------------------------------------------------
//   getConfig
//--------------------------------------------------------------------------------------------------


void PRNNMaterial::getConfig 

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties myConf = conf.makeProps ( myName_ );

  myConf.set ( RANK_PROP       , impl_->rank       );
  myConf.set ( STATE_PROP      , impl_->stateName  );
  myConf.set ( WEIGHTS_PROP    , impl_->weights    );
  myConf.set ( MODULE_PROP     , impl_->module     );
  myConf.set ( LOADER_PROP     , impl_->loader     );
  myConf.set ( PYTHON_PATH_PROP, impl_->pythonPath );

  myConf.makeProps( MATERIAL_PROP ).mergeWith( impl_->materialProps );
}


//--------------------------------------------------------------------------------------------------
//   commit
//--------------------------------------------------------------------------------------------------


void PRNNMaterial::commit 

  ()

{
  py::gil_scoped_acquire gil;
  try 
  {
    for ( auto& model : impl_->models ) 
    {
      model.attr ( "commit" )();
    }
  }
  catch ( const py::error_already_set&  error ) 
  {
    throw jem::IllegalOperationException ( getContext(), String( error.what() ) );
  }
}


//--------------------------------------------------------------------------------------------------
//   cancel
//--------------------------------------------------------------------------------------------------


void PRNNMaterial::cancel 

  ()

{
  py::gil_scoped_acquire gil;
  try 
  {
    for ( auto& model : impl_->models ) 
      model.attr ("cancel")();
  }
  catch ( const py::error_already_set&  error ) 
  {
    throw jem::IllegalOperationException ( getContext(), String( error.what() ) );
  }
}


//--------------------------------------------------------------------------------------------------
//   allocPoints
//--------------------------------------------------------------------------------------------------


void PRNNMaterial::allocPoints 

  ( idx_t  count )

{
  JEM_PRECHECK ( impl_->configured && count >= 0_idx );

  py::gil_scoped_acquire gil;
  try
  {
    for (idx_t i = 0; i < count; i++ )
    {
      impl_->models.push_back( 
        impl_->loadModel( text_( impl_->weights ), py::cast( impl_->factory ) ) 
      );
    }
  }
  catch ( const py::error_already_set&  error ) 
  {
    throw jem::IllegalInputException ( getContext(), String( error.what() ) );
  }
}



//--------------------------------------------------------------------------------------------------
//   rank
//--------------------------------------------------------------------------------------------------


idx_t PRNNMaterial::rank 

  () const noexcept 

{
  return impl_->rank;
}


//--------------------------------------------------------------------------------------------------
//   getProblemType
//--------------------------------------------------------------------------------------------------


ProblemType PRNNMaterial::getProblemType 

  () const noexcept 

{
  return impl_->state;
}


//--------------------------------------------------------------------------------------------------
//   elasticUpdate
//--------------------------------------------------------------------------------------------------


void PRNNMaterial::elasticUpdate 

  ( const Vector&  stress,
    const Matrix&  stiff,
    const Vector&  strain ) const

{
  throw jem::IllegalOperationException ( 
    getContext(), "PRNN has no separate elastic constitutive law!" 
  );
}


Vector PRNNMaterial::elasticUpdate 

  ( const Vector&  strain ) const

{
  throw jem::IllegalOperationException (
    getContext(), "PRNN has no separate elastic constitutive law!"
  );
}


//--------------------------------------------------------------------------------------------------
//   update
//--------------------------------------------------------------------------------------------------


void PRNNMaterial::update 

  ( const Vector&  stress,
    const Matrix&  stiff,
    const Vector&  strain,
    idx_t          ip ) const

{
  JEM_PRECHECK ( 
    ip >= 0_idx                                     && 
    ip < static_cast<idx_t>( impl_->models.size() ) && 
    strain.size() == impl_->components              &&
    stress.size() == impl_->components              && 
    stiff.size(0) == impl_->components              && 
    stiff.size(1) == impl_->components 
  );

  py::gil_scoped_acquire gil;
  try
  {
    auto result = impl_->runtime.attr( "evaluate_with_tangent" )( 
      impl_->models[ip], input_( strain ) ).cast<py::tuple>();

    auto sig    = result[0].cast<py::array_t<double>>();
    auto c      = result[1].cast<py::array_t<double>>();

    if ( sig.ndim()   != 1                 || 
         sig.shape(0) != impl_->components || 
         c.ndim()     != 2                 ||
         c.shape(0)   != impl_->components || 
         c.shape(1)   != impl_->components )
    {
      throw jem::IllegalInputException ( getContext(), "Invalid Python response dimensions!" );
    }

    auto s = sig.unchecked<1> ();
    auto k = c.unchecked<2>   ();

    for ( idx_t i = 0; i < impl_->components; i++ )
    {
      stress[i] = s(i);

      for ( idx_t j = 0; j < impl_->components; j++ ) 
      {
        stiff(i, j) = k(i, j);
      }
    }
  }
  catch ( const py::error_already_set&  error )
  {
    impl_->models[ip].attr( "cancel" )();

    throw jem::IllegalOperationException ( getContext(), String( error.what() ) );
  }
}


void PRNNMaterial::update 

  ( const Vector&  stress,
    const Vector&  strain,
    idx_t          ip ) const

{
  JEM_PRECHECK ( 
    ip >= 0_idx                                     && 
    ip < static_cast<idx_t>( impl_->models.size() ) && 
    strain.size() == impl_->components              && 
    stress.size() == impl_->components
  );

  py::gil_scoped_acquire gil;
  try
  {
    auto result = impl_->runtime.attr( "observe" )( 
        impl_->models[ip], input_( strain ) ).cast<py::array_t<double>>();

    if ( result.ndim() != 1 || result.shape(0) != impl_->components )
    {
      throw jem::IllegalInputException ( getContext(), "Invalid Python stress dimensions!" );
    }

    auto s = result.unchecked<1> ();

    for ( idx_t i = 0; i < impl_->components; i++ ) 
    {
      stress[i] = s(i);
    }
  }
  catch ( const py::error_already_set&  error ) 
  {
    throw jem::IllegalOperationException ( getContext(), String( error.what() ) );
  }
}


//--------------------------------------------------------------------------------------------------
//   getStiffMat
//--------------------------------------------------------------------------------------------------


Matrix PRNNMaterial::getStiffMat 

  () const

{
  throw jem::IllegalOperationException (
    getContext(), "Use the tangent returned by PRNN update!"
  );
}


//--------------------------------------------------------------------------------------------------
//   fill3DStrain
//--------------------------------------------------------------------------------------------------


Tuple<double,6> PRNNMaterial::fill3DStrain 

  ( const Vector&  v ) const

{
  if ( impl_->rank == 3_idx ) 
  {
    return FEMUtils::fillFrom3D ( v );
  }

  if ( impl_->state == ProblemType::PlaneStrain ) 
  {
    return FEMUtils::fillFrom2D ( v, 0.0 );
  }

  throw jem::IllegalOperationException ( 
    getContext(), "Out-of-plane strain is not a network output!"
  );
}


//--------------------------------------------------------------------------------------------------
//   fill3DStress
//--------------------------------------------------------------------------------------------------


Tuple<double,6> PRNNMaterial::fill3DStress 

  ( const Vector&  v ) const

{
  if ( impl_->rank == 3_idx ) 
  {
    return FEMUtils::fillFrom3D ( v );
  }

  if ( impl_->state == ProblemType::PlaneStress ) 
  {
    return FEMUtils::fillFrom2D ( v, 0.0 );
  }

  throw jem::IllegalOperationException ( 
    getContext(), "Out-of-plane stress is not a network output!"
  );
}


//--------------------------------------------------------------------------------------------------
//   makeNew
//--------------------------------------------------------------------------------------------------


Ref<Material> PRNNMaterial::makeNew 

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  return jem::newInstance<Self> ( name, conf, props, globdat );
}


//--------------------------------------------------------------------------------------------------
//   declare
//--------------------------------------------------------------------------------------------------


void PRNNMaterial::declare 

  () 

{
  MaterialFactory::declare ( TYPE_NAME,  & makeNew );
  MaterialFactory::declare ( CLASS_NAME, & makeNew );
}
