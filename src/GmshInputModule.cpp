
/*
 *  Copyright (C) 2014 TU Delft. All rights reserved.
 *  
 *  Frans van der Meer, September 2014
 *  
 *  Module to generate NodeGroups and ElementGroups from input file
 *  used as child Module in GmshInputModule
 *
 *  Changelog:
 *
 *  Modified: Luiz Ant. T. Mororo, l.a.taumaturgomororo@tudelft.nl
 *  Date:     July 2020
 * 
 *  The class was changed in order to enable parallel computing.
 *
 *  Modified: Luiz A. T. Mororo, luiz.mororo@ifce.edu.br
 *  Date:     Fev 2025
 * 
 *  Clean-up code, remove parallel functionalities, and refactored member 
 *  functions for better performance, storage scheme, and to meet new 
 *  Jem/Jive (version 3.0) data types, for instance, idx_t. 
 *
 *  NGroup and EGroup are no longer members of this class. They are 
 *  actually helper classes that are defined and implemented in the
 *  corresponding *.cpp file.
 *
 *  Besides, functions that not store NodeGroup and ElementGroup are
 *  designed. They just find nodes and elements.
 */


#include <jem/base/array/operators.h>
#include <jem/base/array/utilities.h>
#include <jem/base/System.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/io/FileReader.h>
#include <jem/numeric/utilities.h>
#include <jem/util/ArrayBuffer.h>
#include <jem/util/Properties.h>
#include <jem/util/SparseArray.h>


#include <jive/app/ModuleFactory.h>
#include <jive/fem/ElementGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/Model.h>
#include <jive/util/StorageMode.h>


#include "GmshInputModule.h"
#include "GroupInputModule.h"


JEM_DEFINE_CLASS ( GmshInputModule );


using namespace jem::literals;


using jive::IdxVector;
using jive::Vector;

using jive::util::StorageMode;


//=======================================================================
//   class GmshInputModule 
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  GmshInputModule::TYPE_NAME       = "GmshInput";

const char*  GmshInputModule::FILE_NAME_PROP  = "file";
const char*  GmshInputModule::DIM_PROP        = "dim";
const char*  GmshInputModule::EL_GRP_IDX_PROP = "eGrpIdx";
const char*  GmshInputModule::STG_MODE_PROP   = "storage";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


GmshInputModule::GmshInputModule

  ( const String&  name ) : Super ( name )

{
  rank_ = 2;
}


GmshInputModule::~GmshInputModule

  ()

{}


//-----------------------------------------------------------------------
//   init 
//-----------------------------------------------------------------------


Module::Status GmshInputModule::init 

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  Properties myProps = props.findProps ( myName_ );
  Properties myConf  = conf.makeProps  ( myName_ );

  // Read filename (if not specified: module inactive).

  myProps.get ( fileName_, FILE_NAME_PROP );
  myConf .set ( FILE_NAME_PROP, fileName_ );

  // Read dimension.

  myProps.find ( rank_, DIM_PROP );
  myConf .set  ( DIM_PROP, rank_ );

  // Read element group indices.

  IntVector eGrpIdx;

  myProps.find ( eGrpIdx, EL_GRP_IDX_PROP );
  myConf .set  ( EL_GRP_IDX_PROP, eGrpIdx );

  // Read storage mode.

  StorageMode stgMode = StorageMode::DEFAULT_STORAGE;

  if ( ! getStorageMode( stgMode, myConf, myProps ) )
  {
    stgMode = StorageMode::DEFAULT_STORAGE;
  }

  // Init sets.

  initItemSets_ ( stgMode, globdat );

  // Read mesh.

  readMesh_ ( fileName_, eGrpIdx, globdat );

  // Ceck higher order elements.

  checkHighOrderMesh_ ();

  // Read possible groups. Make use of GroupInputModule.

  Ref<GroupInputModule> groupInput = 
                          jem::newInstance<GroupInputModule> ( myName_ );

  groupInput->init ( conf, props, globdat );

  return DONE;
}


//-----------------------------------------------------------------------
//   makeNew 
//-----------------------------------------------------------------------


Ref<Module> GmshInputModule::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  return jem::newInstance<Self> ( name );
}


//-----------------------------------------------------------------------
//   declare 
//-----------------------------------------------------------------------


void GmshInputModule::declare

  ()

{
  using jive::app::ModuleFactory;
  
  ModuleFactory::declare ( GmshInputModule::TYPE_NAME , & makeNew );
  ModuleFactory::declare ( GmshInputModule::CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   initItemSets_
//-----------------------------------------------------------------------


void GmshInputModule::initItemSets_

  ( StorageMode        mode, 
    const Properties&  globdat )

{
  using jem::IllegalInputException;

  using jive::fem::ElementSet;
  using jive::fem::NodeSet;
  using jive::fem::newXElementSet;
  using jive::fem::newXNodeSet;


  // Create a new NodeSet.

  if ( NodeSet::find( globdat ) )
  {
    throw IllegalInputException ( 
      CLASS_NAME, "pre-existing NodeSet!" 
    );
  }

  xnodes_ = newXNodeSet ( mode ); 

  // Create a new ElementSet.

  if ( ElementSet::find( globdat ) )
  {
    throw IllegalInputException ( 
      CLASS_NAME, "pre-existing ElementSet!" 
    );
  }

  xelems_ = newXElementSet ( xnodes_, mode );
}

//-----------------------------------------------------------------------
//   readMesh_
//-----------------------------------------------------------------------


void GmshInputModule::readMesh_ 

  ( const String&      fileName,
    const IntVector&   eGrpIdx, 
    const Properties&  globdat )

{
  using jem::BEGIN;

  using jem::io::FileReader;
  using jem::util::ArrayBuffer;
  using jem::util::SparseArray;

  using jive::fem::ElementGroup;


  jem::io::Writer&  out  = jem::System::out  ();
  jem::io::Writer&  warn = jem::System::warn ();

  SparseArray< ArrayBuffer<idx_t> > egroups;
  Assignable < ElementGroup >       newGroup;

  const String ctxt = getContext ();

  Vector       coords ( 3 );

  IdxVector    perm ( 20 ); // max # of nodes per element
  IdxVector    inodes; 
  IdxVector    groupElems;

  String       line;
  String       groupName;

  idx_t        nn;
  idx_t        numNodes;
  idx_t        nodeID;
  idx_t        numElems;

  idx_t        iel;
  idx_t        eltype; 
  idx_t        nrtags; 
  idx_t        igroup; 


  // Open gmsh file.

  Ref<FileReader> file = jem::newInstance<FileReader> ( fileName );
  
  print ( out, ctxt, 
          ": reading input mesh from file `", fileName, "'.\n");

  // Jump over the first four lines of the input file.

  file->readLine ();
  file->readLine ();
  file->readLine ();

  line = file->readLine().stripWhite ();

  if ( line != "$Nodes" )
  {
    throw jem::IllegalInputException ( 
      CLASS_NAME, 
      "Unknown file format. Are any nodes specified in this Gmsh-file?" 
    );
  }

  // Read nodes.

  numNodes = file->parseInt ();

  xnodes_.reserve ( numNodes );

  print ( out, "  ... Adding ", numNodes, " nodes.\n" );

  for ( idx_t in = 0; in < numNodes; ++in)
  {
    nodeID    = file->parseInt   ();

    coords[0] = file->parseFloat ();
    coords[1] = file->parseFloat ();
    coords[2] = file->parseFloat ();

    xnodes_.addNode ( nodeID, coords[jem::slice(BEGIN,rank_)] );
  }

  xnodes_.store ( globdat ); 

  // Read current line and jump over next two ones.

  file->readLine ();
  file->readLine ();

  line = file->readLine().stripWhite ();

  if ( line != "$Elements" )
  {
    throw jem::IllegalInputException ( 
      CLASS_NAME, 
      "Unknown file format. "
      "Are any elements specified in this Gmsh-file?" 
    );
  }

  // Read elements.

  numElems = file->parseInt ();

  xelems_.reserve ( numElems );

  print ( out, "  ... Adding ", numElems, " elements.\n" );

  for ( idx_t ie = 0; ie < numElems; ie++ )
  {
    iel    = file->parseInt ();
    eltype = file->parseInt ();
    nrtags = file->parseInt ();
    igroup = file->parseInt ();

    for ( idx_t j = 1; j < nrtags; j++ )
    {
      file->parseInt ();
    }

    switch ( eltype )
    {
      case  1_idx:   nn = 2_idx;   break; // 2-node line
      case  2_idx:   nn = 3_idx;   break; // T3
      case  3_idx:   nn = 4_idx;   break; // Q4
      case  4_idx:   nn = 4_idx;   break; // TET4 
      case  5_idx:   nn = 8_idx;   break; // HEX8
      case  6_idx:   nn = 6_idx;   break; // WEDGE6/PRISM6
      case  8_idx:   nn = 3_idx;   break; // 3-node line
      case  9_idx:   nn = 6_idx;   break; // T6
      case 10_idx:   nn = 9_idx;   break; // Q9
      case 11_idx:   nn = 10_idx;  break; // TET10
      case 16_idx:   nn = 8_idx;   break; // Q8

      default:

      throw jem::IllegalInputException ( 
        JEM_FUNC, 
        String::format ( 
          "unknown element type: ( iel = %i, eltype = %i )", 
          iel, eltype 
        ) 
      );
    }

    // Get permutation vector (for quadratic elements).

    if      ( eltype == 9_idx )
    {
      perm[0] = 0_idx;
      perm[1] = 2_idx;
      perm[2] = 4_idx;
      perm[3] = 1_idx;
      perm[4] = 3_idx;
      perm[5] = 5_idx;
    }
    else if ( eltype == 11_idx )
    {
    /*
      perm[0] = 0_idx;
      perm[1] = 4_idx;
      perm[2] = 1_idx;
      perm[3] = 5_idx; 
      perm[4] = 2_idx;
      perm[5] = 6_idx;
      perm[6] = 8_idx;
      perm[7] = 7_idx;
      perm[8] = 9_idx;
      perm[9] = 3_idx;
      */

      perm[0] = 0_idx;
      perm[1] = 2_idx;
      perm[2] = 4_idx;
      perm[3] = 9_idx; 
      perm[4] = 1_idx;
      perm[5] = 3_idx;
      perm[6] = 5_idx;
      perm[7] = 6_idx;
      perm[8] = 8_idx;
      perm[9] = 7_idx;
    }
    else if ( eltype == 16_idx )
    {
      perm[0] = 0_idx;
      perm[1] = 2_idx;
      perm[2] = 4_idx;
      perm[3] = 6_idx; 
      perm[4] = 1_idx;
      perm[5] = 3_idx;
      perm[6] = 5_idx;
      perm[7] = 7_idx;
    }
    else
    {
      perm[jem::slice(0,nn)] = jem::iarray ( nn );
    }

    // Read element nodes

    inodes.resize ( nn );

    for ( idx_t in = 0; in < nn; in++ )
    {
      inodes[ perm[in] ] = xnodes_.findNode ( file->parseInt() );
    }

    // ... and elements.

    xelems_.addElement ( inodes );

    // Check if it is necessary to create element groups.

    if ( eGrpIdx.size() > 0 )
    {
      egroups[igroup].pushBack ( ie );
    }
  }

  xelems_.store ( globdat );

  // Store ElemGroup if required.

  if ( eGrpIdx.size() > 0 )
  {
    for ( idx_t igrp = 0; igrp < eGrpIdx.size(); igrp++ )
    {
      idx_t iegroup = (idx_t) eGrpIdx[igrp];

      ArrayBuffer<idx_t>* egroup = egroups.find ( iegroup );

      if ( ! egroup )
      {
        print ( warn, ctxt, ":  ... ElemGroup `gmsh", 
                iegroup, "' not found.\n" );
	
	continue;
      }

      groupElems.ref ( egroup->toArray() );

      groupName = "gmsh" + String ( iegroup ); 

      print ( out, "  ... store ElementGroup `", groupName, "' with ", 
              groupElems.size(), " elements.\n" );

      newGroup = newElementGroup ( groupElems, xelems_ );
      newGroup .store            ( groupName , globdat );
    }
  }

  // Close gmsh file.

  file->close ();
}


//-----------------------------------------------------------------------
//   checkHighOrderMesh_
//-----------------------------------------------------------------------


void GmshInputModule::checkHighOrderMesh_ 

  () const

{
  using jem::ALL;

  using jive::Matrix;


  // When ( maxx > tiny ) or ( maxy > tiny ) there are curved element 
  // boundaries. 

  const idx_t nn = xelems_.maxElemNodeCount (); 
  const idx_t ne = xelems_.size             ();

  Matrix      coords ( rank_, nn );

  IdxVector   inodes ( nn    );

  Vector      dd     ( rank_ );
  Vector      c0;
  Vector      c1;
  Vector      cm;

  double      maxx = -jem::maxOf<double> ();
  double      maxy = -jem::maxOf<double> ();

  if ( nn == 6 || nn == 8 )
  {
    idx_t nb = nn / 2;

    for ( idx_t ie = 0; ie < ne; ie++ )
    {
      xelems_.getElemNodes  ( inodes, ie     );
      xnodes_.getSomeCoords ( coords, inodes );

      for ( idx_t ib = 0; ib < nb; ++ib )
      {
        c0.ref ( coords( ALL, ib * 2 ) );
        c1.ref ( coords( ALL, ( ib * 2 + 2 ) % nn ) );
        cm.ref ( coords( ALL, ib * 2 + 1) );

        dd = (c1-cm) - (cm-c0); 

        maxx = jem::max ( maxx, jem::numeric::abs(dd[0]) );
        maxy = jem::max ( maxy, jem::numeric::abs(dd[1]) );
      }
    }
  }

  //if ( jem::max( maxx, maxy ) > jem::Limits<double>::EPSILON )
  if ( jem::max( maxx, maxy ) > 1.e-07 )
  {
    jem::System::warn() << getContext() 
      << ": the mesh seems to have curved element boundaries.\n"
      << "Consider setting Mesh.SecondOrderLinear = 1 in .geo-file.\n";
  }
}
