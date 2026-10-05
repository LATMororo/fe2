
#include <jem/base/array/utilities.h>


#include "FEMUtils.h"


using namespace jem::literals;


//=======================================================================
//   General constants
//=======================================================================


const idx_t  FEMUtils::STRAIN_COUNTS[4]  = { 0, 1, 3, 6 };
const char*  FEMUtils::DOF_TYPE_NAMES[3] = { "dx", "dy", "dz" }; 


//=======================================================================
//   Helper function: fillFrom3D & fillFrom2D
//=======================================================================


Tuple<double,6> FEMUtils::fillFrom3D

  ( const Vector&  v6 ) 

{
  Tuple<double,6> t6;

  t6[0] = v6[0];
  t6[1] = v6[1];
  t6[2] = v6[2];
  t6[3] = v6[3];
  t6[4] = v6[4];
  t6[5] = v6[5];

  return  t6;
}


Tuple<double,6> FEMUtils::fillFrom2D

  ( const Vector&  v3,
    double         zzval ) 

{
  Tuple<double,6> t6;

  t6[0] = v3[0];
  t6[1] = v3[1];
  t6[2] = zzval;
  t6[3] = v3[2];
  t6[4] = 0.;
  t6[5] = 0.;

  return  t6;
}


//=======================================================================
//   public functions
//=======================================================================


//-----------------------------------------------------------------------
//   These functions compute the B-matrix given an interpolation matrix
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
//   get1DShapeGrads
//-----------------------------------------------------------------------


void FEMUtils::get1DShapeGrads

  ( const Matrix&  b,
    const Matrix&  g )

{
  JEM_ASSERT ( b.size(0) == 1 &&
               g.size(0) == 1 &&
               b.size(1) == g.size(1) );

  b = g;
}


//-----------------------------------------------------------------------
//   get2DShapeGrads
//-----------------------------------------------------------------------


void FEMUtils::get2DShapeGrads

  ( const Matrix&  b,
    const Matrix&  g )

{
  JEM_ASSERT ( b.size(0) == 3 &&
               g.size(0) == 2 &&
               b.size(1) == 2 * g.size(1) );

  const idx_t  nodeCount = g.size (1);


  b = 0.0;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 2 * inode;

    b(0,i + 0) = g(0,inode);
    b(1,i + 1) = g(1,inode);

    b(2,i + 0) = g(1,inode);
    b(2,i + 1) = g(0,inode);
  }
}


//-----------------------------------------------------------------------
//   get3DShapeGrads
//-----------------------------------------------------------------------


void FEMUtils::get3DShapeGrads

  ( const Matrix&  b,
    const Matrix&  g )

{
  JEM_ASSERT ( b.size(0) == 6 &&
               g.size(0) == 3 &&
               b.size(1) == 3 * g.size(1) );

  const idx_t  nodeCount = g.size (1);


  b = 0.0;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 3 * inode;

    b(0,i + 0) = g(0,inode);
    b(1,i + 1) = g(1,inode);
    b(2,i + 2) = g(2,inode);

    b(3,i + 0) = g(1,inode);
    b(3,i + 1) = g(0,inode);

    b(4,i + 1) = g(2,inode);
    b(4,i + 2) = g(1,inode);

    b(5,i + 2) = g(0,inode);
    b(5,i + 0) = g(2,inode);
  }
}


//-----------------------------------------------------------------------
//   getShapeGradsFunc
//-----------------------------------------------------------------------


FEMUtils::ShapeGradsFunc FEMUtils::getShapeGradsFunc 

  ( const idx_t  rank )

{
  JEM_ASSERT ( rank >= 1_idx && rank <= 3_idx );


  if      ( rank == 1_idx )
  {
    return & FEMUtils::get1DShapeGrads;
  }
  else if ( rank == 2_idx )
  {
    return & FEMUtils::get2DShapeGrads;
  }
  else
  {
    return & FEMUtils::get3DShapeGrads;
  }
}


//-----------------------------------------------------------------------
//   matrix of shape functions  N
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
//   get1DShapeFuncs
//-----------------------------------------------------------------------


void FEMUtils::get1DShapeFuncs

  ( const Matrix&  sfuncs,
    const Vector&  n )
{
  sfuncs(0,0) = n[0];
}


//-----------------------------------------------------------------------
//   get2DShapeFuncs
//-----------------------------------------------------------------------


void FEMUtils::get2DShapeFuncs

  ( const Matrix&  s,
    const Vector&  n )
{
  JEM_ASSERT ( s.size(0) == 2 &&
               s.size(1) == 2 * n.size() );

  const idx_t  nodeCount = n.size ();

  s = 0.0;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 2 * inode;

    s(0,i + 0) = n[inode];
    s(1,i + 1) = n[inode];
  }
}


//-----------------------------------------------------------------------
//   get3DShapeFuncs
//-----------------------------------------------------------------------


void FEMUtils::get3DShapeFuncs

  ( const Matrix&  s,
    const Vector&  n )
{
  JEM_ASSERT ( s.size(0) == 3 &&
               s.size(1) == 3 * n.size() );

  const idx_t  nodeCount = n.size ();

  s = 0.0;

  for ( idx_t inode = 0; inode < nodeCount; inode++ )
  {
    idx_t  i = 3 * inode;

    s(0,i + 0) = n[inode];
    s(1,i + 1) = n[inode];
    s(2,i + 2) = n[inode];
  }
}


//-----------------------------------------------------------------------
//   getShapeFuncsFunc
//-----------------------------------------------------------------------


FEMUtils::ShapeFuncsFunc FEMUtils::getShapeFuncsFunc

  ( idx_t  rank )

{
  JEM_ASSERT ( rank >= 1_idx && rank <= 3_idx );


  if      ( rank == 1_idx )
  {
    return & FEMUtils::get1DShapeFuncs;
  }
  else if ( rank == 2_idx )
  {
    return & FEMUtils::get2DShapeFuncs;
  }
  else
  {
    return & FEMUtils::get3DShapeFuncs;
  }
}


//-----------------------------------------------------------------------
//   computing grad(u) from grad(N) and elemDisp
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
//   get1DGradU
//-----------------------------------------------------------------------


void FEMUtils::get1DGradU

  ( const Matrix&  gradU,
    const Matrix&  gradN,
    const Vector&  elemDisp )

{
  using jem::dot;
  using jem::ALL;

  gradU(0,0) = jem::dot ( elemDisp, gradN(0,ALL ) );
}


//-----------------------------------------------------------------------
//   get2DGradU
//-----------------------------------------------------------------------


void FEMUtils::get2DGradU

  ( const Matrix&  gradU,
    const Matrix&  gradN,
    const Vector&  elemDisp )

{
  // gradU = du_i/dx_j,  size: ( rank,rank )
  // gradN = dN_i/dx_j,  size: ( rank,nodeCount )
  // elemDisp = u        numbering: (0x,0y,1x,1y...nx,ny)

  using jem::ALL;
  using jem::END;
  using jem::slice;

  Vector ux   ( elemDisp[slice(0,END,2)] );
  Vector uy   ( elemDisp[slice(1,END,2)] );
  Vector dNdx ( gradN(0,ALL)             );
  Vector dNdy ( gradN(1,ALL)             );

  gradU(0,0) = dot ( ux, dNdx );
  gradU(0,1) = dot ( ux, dNdy );
  gradU(1,0) = dot ( uy, dNdx );
  gradU(1,1) = dot ( uy, dNdy );
}


//-----------------------------------------------------------------------
//   get3DGradU
//-----------------------------------------------------------------------


void FEMUtils::get3DGradU

  ( const Matrix&  gradU,
    const Matrix&  gradN,
    const Vector&  elemDisp )

{
  using jem::ALL;
  using jem::END;
  using jem::slice;

  Vector ux   ( elemDisp[slice(0,END,3)] );
  Vector uy   ( elemDisp[slice(1,END,3)] );
  Vector uz   ( elemDisp[slice(2,END,3)] );
  Vector dNdx ( gradN(0,ALL)             );
  Vector dNdy ( gradN(1,ALL)             );
  Vector dNdz ( gradN(2,ALL)             );

  gradU(0,0) = dot ( ux, dNdx );
  gradU(0,1) = dot ( ux, dNdy );
  gradU(0,2) = dot ( ux, dNdy );
  gradU(1,0) = dot ( uy, dNdx );
  gradU(1,1) = dot ( uy, dNdy );
  gradU(1,2) = dot ( uy, dNdy );
  gradU(2,0) = dot ( uy, dNdx );
  gradU(2,1) = dot ( uy, dNdy );
  gradU(2,2) = dot ( uy, dNdy );
}


//-----------------------------------------------------------------------
//   getGradUFunc
//-----------------------------------------------------------------------


FEMUtils::GradUFunc FEMUtils::getGradUFunc 

  ( idx_t  rank )

{
  JEM_ASSERT ( rank >= 1_idx && rank <= 3_idx );


  if      ( rank == 1_idx )
  {
    return & FEMUtils::get1DGradU;
  }
  else if ( rank == 2_idx )
  {
    return & FEMUtils::get2DGradU;
  }
  else
  {
    return & FEMUtils::get3DGradU;
  }
}

