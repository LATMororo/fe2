
/*
 *  Copyright (C) 2025 UFC-LMCV/IFCE. All rights reserved.
 *
 *  This solver applies PBC on a RVE and computes its homogenized state
 *  quantities [1-2], i.e., the macro constitutive law, namely, macro 
 *  stress field and (tangent) stiffness matrix. PBC are evaluated in 
 *  two modes:
 *
 *  1) `standalone_ == true', the user can run a simple RVE without the 
 *                          whole FE^2 framework, and a macro strain 
 *                          field has to be provided.
 *
 *  2) `standalone_ == false', it is coupled in a FE^2 framework, in 
 *                           which a macro strain is given and attached
 *                           to an integration point on a upper level.
 *
 *  It follows what has been used implemented by TU Delft's EM group 
 *  [3-5], and part of this implementation has been adapted from 
 *  Rocha [5].
 *
 *  [1] R. J. M. Smit, W. A. M. Brekelmans, and H. E. H. Meijer. 
 *      Prediction of the mechanical behavior of nonlinear heterogeneous 
 *      systems by multi-level finite element modeling. Computer Methods 
 *      in Applied Mechanics and Engineering 155 (1) (1998), pp. 181–192.
 *
 *  [2] MultiFreedom Constraints I. 
 *      https://quickfem.com/wp-content/uploads/IFEM.Ch08.pdf. 
 *      Accessed: 02-01-2023.
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
 *  Functions have been defined in order to control macro strain field
 *  scale in a standalone mode. It allows the user to define a 
 *  function(s) that describes the increment/behaviour of a macroscopic
 *  strain component of a given RVE.
 *
 *  Author: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br, luiztaumaturgo@gmail.com
 *  Date:   May 2026
 */


#ifndef MICRO_PBC_SOLVER_MODULE_H
#define MICRO_PBC_SOLVER_MODULE_H


#include <jive/implict/SolverModule.h>


#include "MicroUtils.h"


using jem::idx_t;
using jem::String;
using jem::Ref;


using jem::io::PrintWriter;
using jem::util::Properties;
using jem::numeric::Function;


using MicroUtils::MicroState;


//=======================================================================
//   class MicroPBCSolverModule 
//=======================================================================


class MicroPBCSolverModule : public jive::implict::SolverModule 
{
 public:

  JEM_DECLARE_CLASS       ( MicroPBCSolverModule, SolverModule );

  static const char*        TYPE_NAME;

  static const char*        ESOLVER_PROP;
  static const char*        MSOLVER_PROP;
  static const char*        EPS_M_PROP;
  static const char*        EPS_M_INCR_PROP;
  static const char*        EPS_M_FUNCS_PROP;


  explicit                  MicroPBCSolverModule 

    ( const String&           name = "pbc" );

 virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat )                    override;
/*
  virtual Status            run

    ( const Properties&       globdat )                    override;
*/

  virtual void              shutdown

    ( const Properties&       globdat )                    override;

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat )                    override;

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )              const override;

  virtual void              advance

    ( const Properties&       globdat )                    override;

  virtual void              solve

    ( const Properties&       info,
      const Properties&       globdat )                    override; 

  virtual void              cancel

    ( const Properties&       globdat )                    override; 

  virtual bool              commit

    ( const Properties&       globdat )                    override; 

  virtual void              setPrecision

    ( double                  eps )                        override; 

  virtual double            getPrecision          () const override; 

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat )                    override;

  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare               ();


 protected:

  virtual                  ~MicroPBCSolverModule  ();

 
 private:

  void                      scaleMacroStrain_     

    ( const Properties&       globdat );

  Ref<PrintWriter>          getWriter_            () const;
  
  void                      writeStates_

    ( const MicroState&       mstate, 
      const Properties&       globdat )              const;


 private:
  
  class                     RunMicroData_;

  Ref<RunMicroData_>        runMData_;

  Ref<SolverModule>         solver_;

  jem::Array
       < Ref<Function> >    epsMFuncs_;

  Ref<PrintWriter>          out_;

  jive::Vector              sigM_;
  jive::Vector              epsM_;
  jive::Vector              epsM0_;

  idx_t                     strCount_;

  bool                      standalone_;

};


#endif
