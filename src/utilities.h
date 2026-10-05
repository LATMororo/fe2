
#ifndef UTILITIES_H
#define UTILITIES_H


#include <jive/Array.h>


using jem::idx_t;
using jem::Tuple;


using jive::Matrix;
using jive::Vector;


namespace utilities
{
  // Solving linear equation of one variable, a*x + b = 0.

  void                      solveLinearEqua

    ( double&                 ans,
      double                  a,
      double                  b );


  // Solving quadratic equation, a*x^2 + b*x + c = 0.

  void                      solveQuadEqua

    ( Vector&                 ans,
      double                  a,
      double                  b,
      double                  c );


  // Solving cubic equation, a*x^3 + b*x^2 + c*x + d = 0. Only for
  // functions with real roots, and it returns the number of roots.

  idx_t                     solveCubicEqua

    ( Tuple<double,3>&        ans,
      double                  a,
      double                  b,
      double                  c,
      double                  d );


  // Invert 2x2 matrix (no size check performed!).

  void                      invert2x2

    ( const Matrix&           mat );


  // Compute the ramp function, <x> = 1/2(x+x).

  double                    evalMcAuley

    ( double                  x );

  double                    evalHeaviside

    ( double                  x );

}

#endif
