
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This material model provides an interface between macro- and micro- 
 *  models of a multiscale analysis [1-3]. For each integration point, a 
 *  chain/models of FE modules and models are set. 
 *
 *  This framework follows what has been implemented by Iuri Rocha and 
 *  his collegues at TU Delft [3-5]. Unlike their code, micro solvers are
 *  provided separately, which are responsible for solving the micro BVP 
 *  (Boundary Value Problem) subjected to appropriate BCs obtained from 
 *  macro strain field, and for homogenizing state quantities, i.e., 
 *  macro stress and stiffness.
 *
 *  [1] R. J. M. Smit, W. A. M. Brekelmans, and H. E. H. Meijer. 
 *      Prediction of the mechanical behavior of nonlinear heterogeneous 
 *      systems by multi-level finite element modeling. Computer Methods 
 *      in Applied Mechanics and Engineering 155 (1) (1998), pp. 181–192.
 *
 *  [2] D. Peric, E. A. de Souza Neto, R. A. Feijo, M. Partovi, and 
 *      A. J. C. Molina. On micro-to-macro transitions for multi-scale 
 *      analysis of non-linear heterogeneous materials: unified 
 *      variational basis and finite element implementation. Int. J. 
 *      Numer. Meth. Engng, 87 (1-5) (2011), pp. 149–170.
 *
 *  [3] V. P. Nguyen, O. Lloberas-Valls, M. Stroeven, and L. J. Sluys. 
 *      Computational homogenization for multiscale crack modeling.
 *      Implementational and computational aspects. Int. J. Numer. Meth. 
 *      Engng, 89:192-226, 2011.
 *
 *  [4] V. P. Nguyen, Multiscale failure modelling of quasi-brittle 
 *      materials. PhD dissertation. Delft University of Technology, 
 *      2011.
 *
 *  [5] I. B. C. M. Rocha. Numerical and Experimental Investigation of 
 *      Hygrothermal Aging in Laminated Composites. Ph.D. dissertation. 
 *      Delft University of Technology, 2019.
 *  
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date:   Sep 2025
 *
 *  Changelog:
 *
 *  Modified: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br, luiztaumaturgo@gmail.com
 *  Date: Jul 2026
 *
 *  Allow the user to set a SampleModule settings for some integration 
 *  points, i.e., some RVEs (or all of them).
 * 
 */


#ifndef FE2_MATERIAL_H 
#define FE2_MATERIAL_H  


#include <jem/util/ObjFlex.h>


#include "Material.h"


//=======================================================================
//   class FE2Material
//=======================================================================


class FE2Material : public Material
{
 public:

  JEM_DECLARE_CLASS       ( FE2Material, Material );

  static const char*        TYPE_NAME;

  static const char*        STATE_PROP;
  static const char*        OUT_IPS_PROP;
  static const char*        SAMP_IPS_PROP;


  explicit                  FE2Material

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

  virtual void              commit                ();

  virtual void              allocPoints

    ( idx_t                   count )                      override;

  virtual idx_t             rank                  () const noexcept
                                                           override;

  virtual ProblemType       getProblemType        () const noexcept 
                                                           override;

  virtual void              elasticUpdate

    ( const Vector&           stress,
      const Matrix&           stiff,
      const Vector&           strain )               const override;

  virtual Vector            elasticUpdate

    ( const Vector&           strain )               const override;

  virtual void              update 

    ( const Vector&           stress,
      const Matrix&           stiff,
      const Vector&           strain,
      idx_t                   ip )                   const override;

  virtual void              update 

    ( const Vector&           stress,
      const Vector&           strain,
      idx_t                   ip )                   const override;

  virtual Matrix            getStiffMat           () const override;

  virtual Tuple<double,6>   fill3DStrain

    ( const Vector&           v3 )                   const override;

  virtual Tuple<double,6>   fill3DStress

    ( const Vector&           v3 )                   const override;

  static Ref<Material>      makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare               ();


 private:

  virtual                  ~FE2Material           ();


 private:

  class                     MicroChain_;

  void                      completeStrain_

    ( const Vector&           microStrain,
      const Vector&           macroStrain )          const;
  
  void                      microStressToMacroStress_

    ( const Vector&           macroStress,
      const Vector&           microStress )          const;

  void                      microStiffToMacroStiff_

    ( const Matrix&           macroStiff,
      const Matrix&           microStiff )           const;


 private:
  
  jem::util::ObjFlex        microChains_;

  Properties                conf_;
  Properties                props_;

  IdxVector                 iperm23_;
  jive::IntVector           outIPs_;
  jive::IntVector           sampIPs_;

  String                    stateString_;

  ProblemType               state_;

  idx_t                     rank_;         // macro rank
  idx_t                     strCount_;     // macro strain components

  idx_t                     microRank_;
  idx_t                     microStrCount_;

};


#endif
