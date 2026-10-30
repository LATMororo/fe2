
#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>


#include <jem/base/Exception.h>
#include <jem/base/Thread.h>
#include <jem/io/StringReader.h>
#include <jem/util/Tokenizer.h>


#include "declare.h"
#include "MaterialFactory.h"
#include "PythonMaterialBridge.h"


using namespace jem::literals;


namespace {


//==================================================================================================
//   Thread safe class 
//==================================================================================================


class Executor_ : public jem::Thread 
{
 public:
  
  virtual void                            run                   ()       override;

  template<class F> 
    auto                                  call 

    ( F                                     work ) -> decltype( work() );

 private:

  std::mutex                              mutex_;
  std::condition_variable                 ready_;

  std::queue<std::function<void()>>       queue_;

};


//==================================================================================================
//   Helper functions
//==================================================================================================


void require_ 

  ( bool         condition, 
    const char*  message )

{
  if ( !condition ) throw std::invalid_argument ( message );
}


std::string text_ 

  ( const String&  value )

{
  return std::string ( value.addr(), value.size() );
}


Executor_& executor_ 

  ()

{
  // Shared for process lifetime. No Python objects live on this worker, so a Python caller may hold 
  // the GIL while waiting for the constitutive result.

  static Executor_* instance = [] {
    jem::Thread::getCurrent ();

    auto worker = jem::newInstance<Executor_> ();

    worker->start ();

    auto retained = new Ref<Executor_> ( worker );

    return retained->get ();
  }();

  return *instance;
}


void myJiveDeclares_

  ()

{
  // NOTE: FE2Material sets roughly a second chain of a FE analysis. Maybe, it will require more
  //       models and modules from JIVE.
 
  declareModules   ();
  declareModels    ();
  declareMaterials ();
}


//==================================================================================================
//   Executor_ class 
//==================================================================================================


//--------------------------------------------------------------------------------------------------
//   run
//--------------------------------------------------------------------------------------------------


void Executor_::run 

  ()

{
  for (;;)
  {
    std::function<void()> work;

    {
      std::unique_lock<std::mutex> lock ( mutex_ );

      ready_.wait ( lock, [this] { 
          return !queue_.empty(); 
        } 
      );

      work = std::move ( queue_.front() );
      queue_.pop       ();
    } // lock is going to be released after this block
    
    work ();
  }
}


//--------------------------------------------------------------------------------------------------
//   call
//--------------------------------------------------------------------------------------------------


template<class F> 

  auto Executor_::call 

  ( F  work ) -> decltype( work() )

{
  if ( jem::Thread::getCurrent() == this ) return work ();

  using Result = decltype ( work() );

  auto task = std::make_shared<std::packaged_task<Result()>> ( [ work = std::move( work ) ] {
      try 
      { 
        return work(); 
      }
      catch ( const jem::Exception&  error ) 
      { 
        throw std::runtime_error ( text_( error.where() ) +": "+ text_( error.what() ) ); 
      }
    }
  );

  auto result = task->get_future ();

  // Push work into queue_; thefore, lock to avoid race condition.
  {
   std::lock_guard<std::mutex> lock ( mutex_ );
   
   queue_.push( [task] { (*task)(); } );
  }

  ready_.notify_one ();

  return result.get ();
}


} // End of unnamed namespace


//==================================================================================================
//   class PythonMaterialBatch::State_
//==================================================================================================


class PythonMaterialBatch::State_
{
 public:

  MaterialResponse                        evaluate 

    ( const Ref<Material>&                  model, 
      const std::vector<double>&            strain )               const;
 

  std::vector<double>                     replay 

    ( const std::vector<double>&            history, 
      idx_t                                 steps,
      const std::vector<double>&            strain )               const;

  std::function<Ref<Material>()>          newJiveMaterial;

  Ref<Material>                           material;
  Ref<Material>                           scratch; // for derivatives; created on first use

  idx_t                                   count;
  idx_t                                   components;
  idx_t                                   historyCount;

  bool                                    pending = false;

};


//--------------------------------------------------------------------------------------------------
//   evaluate
//--------------------------------------------------------------------------------------------------


MaterialResponse PythonMaterialBatch::State_::evaluate 

  ( const Ref<Material>&        model, 
    const std::vector<double>&  strain ) const

{
  // Do some checks.

  require_ ( 
    strain.size() == static_cast<size_t>( count * components ), "Wrong material strain dimensions!" 
  );


  require_ ( // check inf strain components 
    std::all_of( strain.begin(), strain.end(), []( double x ) { return std::isfinite( x ); } ),
    "Strain must be finite!" 
  );

  Matrix           stiff   ( components, components );

  Vector           eps     ( components );
  Vector           sig     ( components );

  MaterialResponse result;

  result.stress.resize ( count * components );
  result.stiff.resize  ( count * components * components );

  // Loop over all 'material' points.

  for ( idx_t p = 0; p < count; p++ )
  {
    // Loop over all components.
    
    for ( idx_t i = 0; i < components; i++) eps[i] = strain[ p*components+i ];

    // Update/Compute stress vector and stiffness matrix for a given point.

    model->update ( sig, stiff, eps, p );

    for (idx_t i = 0; i < components; i++ )
    {
      result.stress[ p*components+i ] = sig[i];

      for ( idx_t j = 0; j < components; j++ )
      {
        result.stiff[ (p*components+i)*components+j ] = stiff( i, j );
      }
    }
  }

  // Check if there is any infinite value for both stress and stiffness.

  require_ ( 
    std::all_of( result.stress.begin(), result.stress.end(), 
      []( double x ){ return std::isfinite( x ); } 
    ), 
    "Nonfinite material stress!"
  );

  require_ ( 
    std::all_of( result.stiff.begin(),result.stiff.end(),
      []( double x ){ return std::isfinite( x ); } 
    ), 
    "Nonfinite stiffness matrix!"
  );

  return result;
}


//--------------------------------------------------------------------------------------------------
//   replay
//--------------------------------------------------------------------------------------------------


std::vector<double> PythonMaterialBatch::State_::replay 

  ( const std::vector<double>&  history, 
    idx_t                       steps,
    const std::vector<double>&  strain ) const

{
  Ref<Material> model = newJiveMaterial ();
  const idx_t   width = count * components;

  for ( idx_t t = 0; t < steps; t++ )
  {
    std::vector<double> eps ( history.begin() + t*width, history.begin() + ( t + 1 )*width );

    evaluate      ( model, eps );
    model->commit ();
  }

  return evaluate( model, strain ).stress;
}


//==================================================================================================
//   class PythonMaterialFactory::Specification_
//==================================================================================================


class PythonMaterialFactory::Specification_
{
 public:

  static const char*                      MATERIAL;
  static const char*                      TYPE;

  Ref<Material>                           newJiveMaterial 

    ( idx_t                                 count )                const;

  Properties                              props;
  Properties                              globdat;

};


//--------------------------------------------------------------------------------------------------
//   static data
//--------------------------------------------------------------------------------------------------


const char*  PythonMaterialFactory::Specification_::MATERIAL = "material";
const char*  PythonMaterialFactory::Specification_::TYPE     = "type";


//--------------------------------------------------------------------------------------------------
//   newJiveMaterial
//--------------------------------------------------------------------------------------------------


Ref<Material> PythonMaterialFactory::Specification_::newJiveMaterial 

  ( idx_t  count ) const

{
  Properties conf;
  String     type;

  // Get material type.

  props.getProps ( MATERIAL ).get( type, TYPE );

  // PRNN material cannot be used.

  require_ ( type != "PRNN", "A PRNN cannot recursively be its own physical material layer" );

  // Use MaterialFactory to get a new material

  Ref<Material> material = newMaterial ( MATERIAL, conf, props, globdat );

  // ... and configure and allocate 'material points'.
  
  material->configure ( props, globdat );

  if ( count > 0_idx ) material->allocPoints ( count );

  return material;
}


//==================================================================================================
//   class PythonMaterialFactory 
//==================================================================================================


//--------------------------------------------------------------------------------------------------
//   constructor
//--------------------------------------------------------------------------------------------------


PythonMaterialFactory::PythonMaterialFactory

  ( const std::string&  props )

{
  auto raw = executor_().call( [&] {
      myJiveDeclares_ ();

      auto spec      = std::make_unique<Specification_>        ();
    
      auto reader    = jem::newInstance<jem::io::StringReader> ( String( props.c_str() ) );

      auto tokenizer = jem::newInstance<jem::util::Tokenizer>  ( reader );

      spec->props.parseFrom ( *tokenizer, "Python material configuration" );
      spec->newJiveMaterial ( 0_idx ); // validate settings without creating online history
    
      return spec.release   ();
    }
  );

  spec_ = std::shared_ptr<Specification_> ( raw, [] ( Specification_* p ) {
      executor_().call( [p] { delete p; } );
    }
  );
}


PythonMaterialFactory::PythonMaterialFactory

  ( const Properties&  matProps,
    const Properties&  globdat )

{
  auto raw=executor_().call ( [&] {
    myJiveDeclares_ ();

    auto spec = std::make_unique<Specification_> ();

    spec->props.set ( Specification_::MATERIAL, matProps );

    spec->globdat = globdat;

    spec->newJiveMaterial ( 0_idx );

    return spec.release   ();
  });

  spec_ = std::shared_ptr<Specification_> ( raw, [] ( Specification_* p ) { 
      executor_().call( [p] { delete p; } ); 
    } 
  );
}


//--------------------------------------------------------------------------------------------------
//   create
//--------------------------------------------------------------------------------------------------


std::shared_ptr<PythonMaterialBatch> PythonMaterialFactory::create 

  ( idx_t  count ) const

{
  require_ ( count > 0_idx, "Material point count must be positive" );

  auto state = executor_().call( [&] {
    auto value = std::make_unique<PythonMaterialBatch::State_> ();

    // Instantiate a new material, and set material points.

    value->newJiveMaterial = [ spec = spec_, count ] { return spec->newJiveMaterial( count ); };

    value->material        = value->newJiveMaterial ();
    value->count           = count;

    // Cehck material rank.

    const idx_t rank = value->material->rank  ();

    require_ ( rank >= 1_idx && rank <= 3_idx, "Unsupported material rank!" );

    value->components   = FEMUtils::STRAIN_COUNTS[ rank ];
    value->historyCount = value->material->historyCount ();

    return value.release ();
  });

  return std::shared_ptr<PythonMaterialBatch> ( 
    new PythonMaterialBatch( state, count, state->components, state->historyCount ) 
  );
}


//--------------------------------------------------------------------------------------------------
//   types
//--------------------------------------------------------------------------------------------------


std::vector<std::string> PythonMaterialFactory::types 

  ()

{
  return executor_().call ( [] {
      myJiveDeclares_ ();

      StringVector             names = MaterialFactory::listKnownTypes ();
      std::vector<std::string> result;

      for ( idx_t i = 0; i < names.size(); i++ ) 
      {
        result.push_back( text_( names[i] ) );
      }

      return result;
    }
  );
}


//==================================================================================================
//   class PythonMaterialBatch
//==================================================================================================


//--------------------------------------------------------------------------------------------------
//   constructor & destructor
//--------------------------------------------------------------------------------------------------


PythonMaterialBatch::PythonMaterialBatch 

  ( State_*  state, 
    idx_t    count, 
    idx_t    components,
    idx_t    historyCount )

  : state_( state ), count_( count ), components_( components ), historyCount_( historyCount ) 

{}

PythonMaterialBatch::~PythonMaterialBatch 

  () 

{ 
  executor_().call ( [this] { delete state_; } ); 
}


//--------------------------------------------------------------------------------------------------
//   count
//--------------------------------------------------------------------------------------------------


idx_t PythonMaterialBatch::count 

  () const noexcept 

{
  return count_;
}


//--------------------------------------------------------------------------------------------------
//   components
//--------------------------------------------------------------------------------------------------


idx_t PythonMaterialBatch::components 

  () const noexcept 

{
  return components_;
}


//--------------------------------------------------------------------------------------------------
//   historyCount
//--------------------------------------------------------------------------------------------------


idx_t PythonMaterialBatch::historyCount 

  () const noexcept 

{
  return historyCount_;
}


//--------------------------------------------------------------------------------------------------
//   update 
//--------------------------------------------------------------------------------------------------


MaterialResponse PythonMaterialBatch::update 

  ( const std::vector<double>&  strain )

{
  return executor_().call( [&] {
    try
    {
      MaterialResponse result = state_->evaluate ( state_->material, strain );
      state_->pending         = true;

      return result;
    }
    catch (...) 
    {
      state_->material->cancel ();
      state_->pending = false;

      throw;
    }
  });
}


//--------------------------------------------------------------------------------------------------
//   commit
//--------------------------------------------------------------------------------------------------


void PythonMaterialBatch::commit 

  ()

{
  executor_().call( [this] {
    if ( state_->pending ) 
    {
      state_->material->commit ();
    }

    state_->pending = false;
  });
}


//--------------------------------------------------------------------------------------------------
//   cancel
//--------------------------------------------------------------------------------------------------


void PythonMaterialBatch::cancel 

  ()

{
  executor_().call( [this] { 
    state_->material->cancel (); 
    state_->pending = false;
  });
}


//--------------------------------------------------------------------------------------------------
//   historyVjp - Vector-Jacobian Product
//--------------------------------------------------------------------------------------------------
//
// It is used during training. The native material tangent gives the derivative of current stress 
// with respect to current strain, holding accepted history fixed. Training also needs the 
// dependence on 'earlier strains': for instance, earlier plastic loading changes the stress during 
// later unloading.
//
// Given past encoded strains H and the incoming stress gradient g_sigma, it computes:
//
//   
//   g_H = (\frac{ \partial \sigma_t }{ \partial H })^T g_sigma
//   
//
// It approximates this by perturbing each past strain component, replaying the accepted history 
// through fresh native material instances, and evaluating the resulting current stress. The live 
// material state stays intact. PyTorch then propagates these gradients through the earlier encoder 
// operations.
//
// This is expensive, but avoids requiring analytical history derivatives from every native 
// material. Online Newton tangents do not use this replay; they use the native tangent combined 
// with derivatives of the neural layers.


std::vector<double> PythonMaterialBatch::historyVjp 

  ( const std::vector<double>&  history,
    idx_t                       steps, 
    const std::vector<double>&  strain, 
    const std::vector<double>&  adjoint, 
    double                      step )
{
  require_ ( 
    steps >= 0.0 && history.size() == static_cast<size_t>( steps*count_*components_ ), 
    "Wrong material history dimensions!" 
  );

  require_ ( 
    adjoint.size() == static_cast<size_t>( count_*components_ ), 
    "Wrong stress adjoint dimensions!" 
  );

  require_ ( 
    step > 0.0 && std::isfinite( step ), "Finite-difference step must be positive and finite!" 
  );

  return executor_().call( [&] {
    // Isolated replicas: differentiation never advances the live material.

    std::vector<double> result    ( history.size() );
    std::vector<double> perturbed = history;

    for ( size_t i = 0; i < history.size(); i++ )
    {
      const double h = step * std::max( 1.0, std::abs( history[i] ) );

      perturbed[i]              = history[i] + h;
      std::vector<double> plus  = state_->replay ( perturbed, steps, strain );

      perturbed[i]              = history[i] - h;
      std::vector<double> minus = state_->replay ( perturbed, steps, strain );

      perturbed[i] = history[i];

      for ( size_t j = 0; j < adjoint.size(); j++ ) 
      {
        result[i] += ( plus[j] - minus[j] ) * adjoint[j] / ( 2 * h );
      }
    }

    return result;
  });
}


//--------------------------------------------------------------------------------------------------
//   history
//--------------------------------------------------------------------------------------------------


std::vector<double> PythonMaterialBatch::history 

  ()

{
  require_ ( historyCount_ >= 0_idx, "Material does not expose its history variables!" );

  return executor_().call( [this] {
    std::vector<double> result ( static_cast<size_t>( count_*historyCount_ ) );
    Vector              hist   ( historyCount_ );

    for ( idx_t p = 0; p < count_; p++ )
    {
      state_->material->getHistory ( hist, p );

      for ( idx_t k = 0; k < historyCount_; k++ ) result[ p*historyCount_+k ] = hist[k];
    }

    return result;
  });
}


//--------------------------------------------------------------------------------------------------
//   stepVjp - Vector-Jacobian Product of a single step
//--------------------------------------------------------------------------------------------------
//
// Training treats each material step as a function of the current strain and the history before
// the step, (eps_t, alpha_{t-1}) -> (sigma_t, alpha_t). PyTorch then chains these steps through
// time, as with any recurrent cell, so the cost of a backward pass is linear in the sequence 
// length. Given the adjoints g_sigma and g_alpha, it returns:
//
//   g_eps   = ( d alpha_t / d eps_t )^T g_alpha
//   g_alpha = ( d sigma_t / d alpha_{t-1} )^T g_sigma + ( d alpha_t / d alpha_{t-1} )^T g_alpha
//
// The term ( d sigma_t / d eps_t )^T g_sigma is left to the caller, which has the native tangent.
// Derivatives are central differences of one step at one point, evaluated on a scratch material 
// so the live material state stays intact. Points are independent, so each perturbation only 
// touches its own point.


StepGradients PythonMaterialBatch::stepVjp 

  ( const std::vector<double>&  history,
    const std::vector<double>&  strain, 
    const std::vector<double>&  stressAdjoint, 
    const std::vector<double>&  historyAdjoint, 
    double                      step )
{
  const idx_t nc = components_;
  const idx_t nh = historyCount_;

  require_ ( nh >= 0_idx, "Material does not expose its history variables!" );

  require_ ( 
    history.size() == static_cast<size_t>( count_*nh ) && 
    historyAdjoint.size() == history.size(),
    "Wrong material history dimensions!" 
  );

  require_ ( 
    strain.size() == static_cast<size_t>( count_*nc ) && 
    stressAdjoint.size() == strain.size(),
    "Wrong material strain or stress adjoint dimensions!" 
  );

  require_ ( 
    step > 0.0 && std::isfinite( step ), "Finite-difference step must be positive and finite!" 
  );

  return executor_().call( [&] {
    if ( state_->scratch == nullptr ) state_->scratch = state_->newJiveMaterial ();

    Material&     model = *state_->scratch;
    StepGradients result;

    result.strain .assign ( strain .size(), 0.0 );
    result.history.assign ( history.size(), 0.0 );

    Matrix stiff ( nc, nc );
    Vector eps   ( nc ), sigPlus ( nc ), sigMinus ( nc );
    Vector alpha ( nh ), newPlus ( nh ), newMinus ( nh );

    for ( idx_t p = 0; p < count_; p++ )
    {
      for ( idx_t i = 0; i < nc; i++ ) eps[i]   = strain [ p*nc+i ];
      for ( idx_t k = 0; k < nh; k++ ) alpha[k] = history[ p*nh+k ];

      // Perturb strain components first, then history variables.

      for ( idx_t i = 0; i < nc + nh; i++ )
      {
        double& x     = i < nc ? eps[i] : alpha[i-nc];
        const double x0 = x;
        const double h  = step * std::max( 1.0, std::abs( x0 ) );

        x = x0 + h;
        model.setHistory ( alpha, p );
        model.update     ( sigPlus, stiff, eps, p );
        model.getHistory ( newPlus, p );

        x = x0 - h;
        model.setHistory ( alpha, p );
        model.update     ( sigMinus, stiff, eps, p );
        model.getHistory ( newMinus, p );

        x = x0;

        double g = 0.0;

        for ( idx_t k = 0; k < nh; k++ ) 
        {
          g += ( newPlus[k] - newMinus[k] ) * historyAdjoint[ p*nh+k ];
        }

        if ( i < nc )
        {
          result.strain[ p*nc+i ] = g / ( 2*h );
        }
        else
        {
          for ( idx_t j = 0; j < nc; j++ ) 
          {
            g += ( sigPlus[j] - sigMinus[j] ) * stressAdjoint[ p*nc+j ];
          }

          result.history[ p*nh+i-nc ] = g / ( 2*h );
        }
      }
    }

    // Leave the scratch material without a pending trial.

    model.cancel ();

    return result;
  });
}
