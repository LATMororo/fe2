
#include <jem/base/Error.h>
#include <jem/numeric/utilities.h>


#include "utilities.h"


//-----------------------------------------------------------------------
//   constants
//-----------------------------------------------------------------------


const double PI       = 3.14159265;
const double RAD_120  = 0.666666667 * PI;
const double EPS      = jem::Limits<double>::EPSILON;


//-----------------------------------------------------------------------
//   solveLinearEqua
//-----------------------------------------------------------------------


void utilities::solveLinearEqua

  ( double&  ans,
    double   a,
    double   b )
{
  ans = -b/a;
}

//-----------------------------------------------------------------------
//   solveQuadEqua
//-----------------------------------------------------------------------


void utilities::solveQuadEqua

  ( Vector&  ans,
    double   a,
    double   b,
    double   c )
{
  if ( jem::numeric::abs(a) < EPS )
  {
    ans.resize( 1 );
    ans[0] = - c / b;
  }
  else
  {
    double d = b * b - 4.0 * a * c;
    
    if ( d < 0 )
    {
      throw jem::Error (
        JEM_FUNC,
        "imaginary principal values !!! "
      );  
    }
    else
    {
      ans.resize( 2 );
      ans[0] = 0.5 * (-b + sqrt( d )) / a;
      ans[1] = 0.5 * (-b - sqrt( d )) / a;
    }
  }
}


//-----------------------------------------------------------------------
//   solveCubicEqua
//-----------------------------------------------------------------------


idx_t utilities::solveCubicEqua

  ( Tuple<double,3>&  ans, 
    double            a,
    double            b,
    double            c,
    double            d)

{  
  ans = NAN;

  if ( jem::numeric::abs( a ) < EPS )
  {
    if ( jem::numeric::abs( b ) < EPS )
    {
      ans[0] = - d / c;

      return 1;
    }
    else
    {
      double D = c * c - 4.0 * b * d;
      
      if ( D < 0.0 )
      {        
        return 0;
      }
      else
      {
        D      = sqrt(D);

        ans[0] = 0.5 * (-c + D) / b;
        ans[1] = 0.5 * (-c - D) / b;
       
        return 2;
      }
    }
  }
  else
  {
    double kk  = 1.0 / a;
    
    double aa  = b * kk;
    double bb  = c * kk;
    double cc  = d * kk;

    double aa3 = aa / 3.0;

    double p, q, r;
    double phi, help;
   
    q  = ( aa * aa - 3.0 * bb ) / 9.0;
    r  = ( 2.0 * aa * aa * aa - 9.0 * aa * bb + 27.0 * cc ) / 54.0;
 
    help = r / sqrt( q * q * q );
   
    if ( jem::numeric::abs( help ) > 1.0 )
    {
      help = ( help < 0 ) ? -1.0 : 1.0; // prevent rounding errors
    }
   
    phi = acos ( help );
    p   = sqrt ( q    );

    ans[0] = -2.0 * p * cos ( phi / 3.0 )           - aa3;
    ans[1] = -2.0 * p * cos ( phi / 3.0 + RAD_120 ) - aa3;
    ans[2] = -2.0 * p * cos ( phi / 3.0 - RAD_120 ) - aa3;

    return 3;
  }
}
  
//-----------------------------------------------------------------------
//   invert2x2
//-----------------------------------------------------------------------

void utilities::invert2x2

  ( const Matrix&  m )

{
  // no check for size and singularity!

  double d = m(0,0) * m(1,1) - m(0,1) * m(1,0);
  double s = 1. / d;

  double t = m(0,0);

  m(0,0) =  m(1,1) * s;
  m(0,1) = -m(0,1) * s;
  m(1,0) = -m(1,0) * s;
  m(1,1) =  t      * s;
}


//-----------------------------------------------------------------------
//   evalMcAuley
//-----------------------------------------------------------------------


double utilities::evalMcAuley

  ( double  x )
{
  return (x > 0 ? x : 0.0);
}


//-----------------------------------------------------------------------
//    evalHeaviside
//-----------------------------------------------------------------------


double utilities::evalHeaviside

  ( double  x )
{
  return (x < 0 ? -1.0 : 1.0);
} 

