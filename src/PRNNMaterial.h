
#ifndef PYTHON_PRNN_MATERIAL_H
#define PYTHON_PRNN_MATERIAL_H


#include <memory>


#include "Material.h"


//==================================================================================================
//   class PRNNMaterial 
//==================================================================================================

// Thin online driver. The Python model owns its architecture and material context; Jive owns the 
// accepted/rejected increment lifecycle.


class PRNNMaterial : public Material
{
 public:

  JEM_DECLARE_CLASS                     ( PRNNMaterial, Material );

  static const char*                      TYPE_NAME;

  static const char*                      WEIGHTS_PROP;
  static const char*                      MODULE_PROP;
  static const char*                      LOADER_PROP;
  static const char*                      PYTHON_PATH_PROP;
  static const char*                      MATERIAL_PROP;
  static const char*                      STATE_PROP;

                                          PRNNMaterial 

    ( const String&                         name, 
      const Properties&                     conf,
      const Properties&                     props, 
      const Properties&                     globdat );

  void                                    configure 

    ( const Properties&                     props, 
      const Properties&                     globdat )                    override;

  void                                    getConfig 

    ( const Properties&                     conf, 
      const Properties&                     globdat )              const override;

  void                                    commit                ()       override;
  void                                    cancel                ()       override;

  void                                    allocPoints 

    ( idx_t                                 count )                      override;

  idx_t                                   rank                  () const noexcept override;
  ProblemType                             getProblemType        () const noexcept override;

  void                                    elasticUpdate 

    ( const Vector&                         stress, 
      const Matrix&                         stiff,
      const Vector&                         strain )               const override;

  Vector                                  elasticUpdate 

    ( const Vector&                         strain )               const override;

  void                                    update 

    ( const Vector&                         stress, 
      const Matrix&                         stiff,
      const Vector&                         strain, 
      idx_t                                 ip )                   const override;

  void                                    update 

    ( const Vector&                         stress, 
      const Vector&                         strain, 
      idx_t                                 ip )                   const override;

  Matrix                                  getStiffMat           () const override;

  Tuple<double,6>                         fill3DStrain 

    ( const Vector&                         strain )               const override;

  Tuple<double,6>                         fill3DStress 

    ( const Vector&                         stress )               const override;

  static Ref<Material>                    makeNew 

    ( const String&                         name, 
      const Properties&                     conf,
      const Properties&                     props, 
      const Properties&                     globdat );

  static void                             declare               ();


 protected:

  virtual                                ~PRNNMaterial          ();


 private:

  class                                   Impl_;
  std::unique_ptr<Impl_>                  impl_;

};


#endif
