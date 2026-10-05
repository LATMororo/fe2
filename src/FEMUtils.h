
#ifndef FEM_UTILS_H 
#define FEM_UTILS_H


//#include <jem/base/array/Array.h>
//#include <jem/base/tuple/Tuple.h>

#include <jive/Array.h>


using jem::Array;
using jem::idx_t;
using jem::Tuple;


using jive::IdxVector;
using jive::Matrix;
using jive::Vector;


namespace FEMUtils
{
  //=====================================================================
  //   General constants
  //=====================================================================


  // An integer array that maps the number of spatial dimensions (1, 2,
  // or 3) to the number of strain/stress components.


  extern const idx_t        STRAIN_COUNTS[4];


  // An string array that maps the number of spatial dimensions to the 
  // type of DOFs ('dx', 'dy', 'dz').


  extern const char*        DOF_TYPE_NAMES[3]; 


  // An enum that maps the types of problems: plane strain, plane stress,
  // and solid (3D problems).


  enum class ProblemType : int 
  { 
                              Bar = 0,
                              PlaneStrain, 
                              PlaneStress,
                              Solid
  }; 


  //=====================================================================
  //   Helper function: fillFrom3D & fillFrom2D
  //=====================================================================


  Tuple<double,6>           fillFrom3D

    ( const Array<double,1>&  v6 ); 


  Tuple<double,6>           fillFrom2D

    ( const Array<double,1>&  v3,
      double                  zzval ); 


  //=====================================================================
  //   typedefs
  //=====================================================================


  // A pointer to a function that computes the spatial derivatives of
  // the interpolation matrix. This is the so-called B-matrix.

  typedef void              (*ShapeGradsFunc)

    ( const Matrix&           b,
      const Matrix&           g );

  // Similar stuff, from the shape functions N_i, compute the shape 
  // matrix N.

  typedef void              (*ShapeFuncsFunc)

    ( const Matrix&           sfuncs,
      const Vector&           n );

  // Similar stuff, from computing grad(U).

  typedef void              (*GradUFunc)

    ( const Matrix&           gradU,
      const Matrix&           gradN,
      const Vector&           elemDisp );


  //=====================================================================
  //   public functions
  //=====================================================================


  //---------------------------------------------------------------------
  // These functions compute the B-matrix given an interpolation matrix
  //---------------------------------------------------------------------

  void                      get1DShapeGrads

    ( const Matrix&           b,
      const Matrix&           g );

  void                      get2DShapeGrads

    ( const Matrix&           b,
      const Matrix&           g );

  void                      get3DShapeGrads

    ( const Matrix&           b,
      const Matrix&           g );


  // A function that returns a pointer to a function that computes the
  // B-matrix given the number of spatial dimensions.


  ShapeGradsFunc            getShapeGradsFunc

    ( idx_t                   rank );


  //---------------------------------------------------------------------
  //   matrix of shape functions  N
  //---------------------------------------------------------------------

  void                      get1DShapeFuncs

    ( const Matrix&           sfuncs,
      const Vector&           n );

  void                      get2DShapeFuncs

    ( const Matrix&           sfuncs,
      const Vector&           n );

  void                      get3DShapeFuncs

    ( const Matrix&           sfuncs,
      const Vector&           n );


  // A function that returns a pointer to a function that computes the
  // B-matrix given the number of spatial dimensions.


  ShapeFuncsFunc            getShapeFuncsFunc

    ( idx_t                   rank );


  //---------------------------------------------------------------------
  //   computing grad(u) from grad(N) and elemDisp
  //---------------------------------------------------------------------

  // gradU = du_i/dx_j,  size: ( rank,rank )
  // gradN = dN_i/dx_j,  size: ( rank,nodeCount )
  // elemDisp = u        numbering: (0x,0y,1x,1y...nx,ny)

  void                      get1DGradU

    ( const Matrix&           gradU,
      const Matrix&           gradN,
      const Vector&           elemDisp );

  void                      get2DGradU

    ( const Matrix&           gradU,
      const Matrix&           gradN,
      const Vector&           elemDisp );

  void                      get3DGradU

    ( const Matrix&           gradU,
      const Matrix&           gradN,
      const Vector&           elemDisp );

  // A function that returns a pointer to a function that computes the
  // B-matrix given the number of spatial dimensions.

  GradUFunc                 getGradUFunc

    ( idx_t                   rank );

}


#endif
