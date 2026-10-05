
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
 *  Besides, functions that do not store NodeGroup and ElementGroup are
 *  designed. They just find nodes and elements.
 */


#include <jem/base/array/utilities.h>
#include <jem/base/array/operators.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/numeric/utilities.h>
#include <jem/util/ArrayBuffer.h>
#include <jem/util/Properties.h>


#include <jive/app/ModuleFactory.h>
#include <jive/fem/ElementGroup.h>
#include <jive/fem/NodeGroup.h>


#include "GroupInputModule.h"


JEM_DEFINE_CLASS ( GroupInputModule );


using namespace jem::literals;


using jem::Array;
using jem::Tuple;


using jive::IdxVector;
using jive::Matrix;
using jive::StringVector;
using jive::Vector;


//=======================================================================
//   class NGroup 
//=======================================================================


/*
 *  Class that generates NodeGroups from input in .pro-file
 *
 *  Use xtype/ytype/ztype = "min"/"max",
 *   or xval/yval/zval = value,
 *   to specify NodeGroup location
 *  Specify 1, 2 or 3 coordinates for plane, line or point. 
 */


class NGroup 
{
 public:

  static const double       PI;
  static const char*        X_NAMES[3];
  static const char*        X_TYPES[4];
  static const char*        NODES_PROP;
  static const char*        ALL_PROP;
  static const char*        EPS_PROP;
  static const char*        X0_PROP;
  static const char*        Y0_PROP;
  static const char*        RADIUS_PROP;
  static const char*        ANGLE_PROP;

  enum  class               Check : int { NONE = 0, VALUE, BOUNDS };


                            NGroup                ();

  void                      configure

    ( const Properties&       conf,
      const Properties&       props,
      const String&           name );

  void                      makeGroup

    ( const NodeSet&          nodes,
      const Properties&       globdat );

  IdxVector                 findNodes 
    
    ( const NodeSet&          nodes, 
      const Properties&       conf,
      const Properties&       props,
      const String&           name );


 private:

  void                      findNodes_

    ( const NodeSet&          nodes );

  void                      store_

    ( const NodeSet&          nodes,
      const Properties&       globdat );


 private:

  // Every entry of doX_ is a switch that determines how a 
  // node is checked in that dimension (x,y,z).
  // If ( all(doX_==0) ) a line is specified.

  Matrix                    xbounds_;

  IdxVector                 inodes_;

  Array<Check>              doX_;  

  Tuple<Vector,3>           xvals_;
  Tuple<String,3>           xtype_;

  StringVector              knownTypes_;
  String                    myName_;

  double                    eps_;
  double                    angle_;
  double                    radius_;

  idx_t                     rank_;

  bool                      all_;

};


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const double  NGroup::PI         = 3.14159265358979323846;
const char*   NGroup::X_NAMES[3] = { "x", "y", "z"};
const char*   NGroup::X_TYPES[4] = { "type", "val", "vals", "bounds" }; 
const char*   NGroup::NODES_PROP = "nodes";
const char*   NGroup::ALL_PROP   = "all";
const char*   NGroup::EPS_PROP   = "eps";
const char*   NGroup::X0_PROP    = "x0";
const char*   NGroup::Y0_PROP    = "y0";
const char*   NGroup::RADIUS_PROP= "radius";
const char*   NGroup::ANGLE_PROP = "angle";


//-----------------------------------------------------------------------
//   constructor 
//-----------------------------------------------------------------------


NGroup::NGroup

  () 

{
  xbounds_ .resize ( 3, 2 );  // ( x/y/z , min/max )

  doX_     .resize ( 3 );

  xvals_[0].resize ( 1 );
  xvals_[1].resize ( 1 );
  xvals_[2].resize ( 1 );

  xtype_ = "";

  knownTypes_.resize ( 5 );

  knownTypes_[0] = "min";
  knownTypes_[1] = "max";
  knownTypes_[2] = "firstHalf";
  knownTypes_[3] = "secondHalf"; 
  knownTypes_[4] = "mid";

  //eps_    = jem::Limits<double>::EPSILON;
  eps_    = 1.e-7;
  angle_  = 0.;
  radius_ = -1.;

  rank_   = 3_idx;
  all_    = false;
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void NGroup::configure 

  ( const Properties&  conf,
    const Properties&  props,
    const String&      name )

{
  // Set group name.

  myName_ = name;

  // Get properties.

  Properties  myConf  = conf .makeProps ( myName_ );
  Properties  myProps = props. getProps ( myName_ );

  // Set nodes.

  if ( myProps.find( inodes_, NODES_PROP ) )
  {
    myConf.set ( NODES_PROP, inodes_ );
    doX_ = Check::NONE;

    return;
  }

  myProps.find ( all_, ALL_PROP );
  myConf.set   ( ALL_PROP, all_ );

  if ( all_ )
  {
    return;
  }

  // Set tolerance.

  myProps.find ( eps_, EPS_PROP );
  myConf.set   ( EPS_PROP, eps_ );

  // Find coordinate types or values.

  for ( idx_t ir = 0; ir < rank_; ir++ )
  {
    Vector bounds;

    String thiss = X_NAMES[ir];

    String xtypename ( thiss + X_TYPES[0] );
    String xvalname  ( thiss + X_TYPES[1] );
    String xvalsname ( thiss + X_TYPES[2] );
    String xinname   ( thiss + X_TYPES[3] );

    doX_[ir] = Check::VALUE;

    if      ( myProps.find( xtype_[ir], xtypename ) )
    {
      myConf.set ( xtypename, xtype_[ir] );

      if ( ! jem::testany( knownTypes_ == xtype_[ir] ) )
      {
        throw jem::IllegalInputException ( 
          props.getContext ( myName_ ), 
          String::format   ( "unknown xtype: `%s'", xtype_[ir] ) 
        );
      }
    }
    else if ( myProps.find( xvals_[ir][0], xvalname ) )
    {
      myConf.set ( xvalname, xvals_[ir][0] );
    }
    else if ( myProps.find( xvals_[ir], xvalsname ) )
    {
      myConf.set ( xvalsname, xvals_[ir] );
    }
    else if ( myProps.find( bounds, xinname ) )
    {
      myConf.set ( xinname, bounds );

      if ( bounds.size() == 2 )
      {
        xbounds_(ir,0) = bounds[0] - eps_;
        xbounds_(ir,1) = bounds[1] + eps_;

        doX_[ir] = Check::BOUNDS;
      }
      else
      {
        throw jem::IllegalInputException ( 
          props.getContext ( myName_ ), 
          String::format   ( "`%s' does not have length 2!", xinname ) 
        );
      }
    }
    else
    {
      doX_[ir] = Check::NONE;
    }
  }

  if ( ! jem::testany( doX_ ) )
  {
    // line specified with x0, y0 and angle
    // or circle with x0, y0 and radius

    myProps.get ( xvals_[0][0], X0_PROP );
    myProps.get ( xvals_[1][0], Y0_PROP );

    myConf.set  ( X0_PROP, xvals_[0][0] );
    myConf.set  ( Y0_PROP, xvals_[1][0] );

    if ( myProps.find( radius_, RADIUS_PROP, 0., jem::maxOf( radius_ ) ) )
    {
      myConf.set ( RADIUS_PROP, radius_ );
    }
    else
    {
      myProps.get ( angle_, ANGLE_PROP );
      myConf.set  ( ANGLE_PROP, angle_ );

      angle_ *= PI / 180.;
    }
  }
}


//-----------------------------------------------------------------------
//   makeGroup
//-----------------------------------------------------------------------


void NGroup::makeGroup

  ( const NodeSet&     nodes,
    const Properties&  globdat )

{
  findNodes_ ( nodes );
  store_     ( nodes, globdat );
}


//-----------------------------------------------------------------------
//   findNodes
//-----------------------------------------------------------------------


IdxVector NGroup::findNodes 
    
  ( const NodeSet&     nodes, 
    const Properties&  conf,
    const Properties&  props,
    const String&      name )

{
  configure  ( conf, props, name );
  findNodes_ ( nodes );

  return inodes_;
}


//-----------------------------------------------------------------------
//   findNodes_
//-----------------------------------------------------------------------


void NGroup::findNodes_

  ( const NodeSet&  nodes )

{
  using jem::util::ArrayBuffer;


  Vector              xcoords ( nodes.size() );
  Vector              coords  ( rank_        );

  ArrayBuffer<idx_t>  nodebuf;

  bool                include;

  if      ( all_ )
  {
    // Add all nodes.

    nodebuf.reserve ( nodes.size() );

    for ( idx_t in = 0; in < nodes.size(); ++in )
    {
      nodebuf.pushBack ( in );
    }
  }
  else if ( jem::testany( doX_ ) ) 
  {
    // Find coordinate if specified with type.

    for ( idx_t ir = 0; ir < nodes.rank(); ir++ )
    {
      if ( jem::testany( xtype_[ir] == knownTypes_ ) )
      {
        nodes.getData()->getXCoords ( xcoords, ir );

        double minX = jem::min ( xcoords );
        double maxX = jem::max ( xcoords );

        if      ( xtype_[ir] == knownTypes_[0] ) // min
        {
          xvals_[ir] = minX;
        }
        else if ( xtype_[ir] == knownTypes_[1] ) // max
        {
          xvals_[ir] = maxX;
        }
        else if ( xtype_[ir] == knownTypes_[2] ) // firstHalf
        {
          xbounds_(ir,0) = minX - eps_;
          xbounds_(ir,1) = minX - eps_ + ( maxX - minX ) / 2.;

          doX_[ir] = Check::BOUNDS;
        }
        else if ( xtype_[ir] == knownTypes_[3] ) // secondHalf
        {
          xbounds_(ir,0) = minX + eps_ + ( maxX - minX ) / 2.;
          xbounds_(ir,1) = maxX + eps_;

          doX_[ir] = Check::BOUNDS;
        }
        else if ( xtype_[ir] == knownTypes_[4] ) // mid
        {
          xvals_[ir] = ( maxX + minX ) / 2.;
        }
      }
    }

    // Loop over nodes, node is included unless it violates the criterion
    // specified for one of the dimensions.

    for ( idx_t in = 0; in < nodes.size(); in++ )
    {
      include = true;

      nodes.getNodeCoords ( coords, in );

      // Loop over dimensions.

      for ( idx_t ir = 0; ir < nodes.rank(); ir++ )
      {
        if ( doX_[ir] == Check::VALUE )
        {
          // Check if coordinate is close to one of specified values.

          bool thisx = false;

          for ( idx_t ival = 0; ival < xvals_[ir].size(); ival++ )
          {
            thisx |= jem::numeric::abs( xvals_[ir][ival] - coords[ir] ) 
                          < eps_ ;
          }

          include &= thisx;
        }
        else if ( doX_[ir] == Check::BOUNDS )
        {
          // check if coordinate is between specified bounds

          include &= ( coords[ir] > xbounds_(ir,0) &&
              coords[ir] < xbounds_(ir,1) );
        }
      }

      if ( include )
      {
        nodebuf.pushBack ( in );
      }
    }
  }
  else if ( radius_ > 0. )
  {
    for ( idx_t in = 0; in < nodes.size(); in++ )
    {
      nodes.getNodeCoords ( coords, in );

      double dx = coords[0] - xvals_[0][0];
      double dy = coords[1] - xvals_[1][0];

      double dist = std::sqrt ( dx * dx + dy * dy );

      if ( jem::numeric::abs( dist - radius_ ) < eps_ )
      {
        nodebuf.pushBack ( in );
      }
    }
  }
  else
  {
    // Find nodes specified with x0, y0, and angle.

    for ( idx_t in = 0; in < nodes.size(); in++ )
    {
      nodes.getNodeCoords ( coords, in );

      double dx = coords[0] - xvals_[0][0];
      double dy = coords[1] - xvals_[1][0];

      double angle = std::atan2 ( dy, dx );
      double dist  = std::sqrt  ( dx * dx + dy * dy ); 

      if ( jem::numeric::abs( angle - angle_ ) < eps_ || 
           jem::numeric::abs( angle - angle_ - PI ) < eps_ ||
           dist < eps_ )
      {
        nodebuf.pushBack ( in );
      }
    }
  }

  inodes_.ref ( nodebuf.toArray() );
}


//-----------------------------------------------------------------------
//   store
//-----------------------------------------------------------------------


void NGroup::store_

  ( const NodeSet&     nodes,
    const Properties&  globdat )

{
  using jive::fem::NodeGroup;


  NodeGroup group = newNodeGroup ( inodes_, nodes );

  group.store ( myName_, globdat );

  jem::System::out() << "  ... Created NodeGroup `" << myName_ << 
    "' with " << inodes_.size() << " nodes.\n";

  if ( inodes_.size() < 1 )
  {
    jem::System::warn() << "NodeGroup `" << myName_ << "' is empty.\n";
  }
}


//=======================================================================
//   class EGroup 
//=======================================================================

/*
 *  Class that generates ElementGroups from input in .pro-file
 *
 *  Use xtype/ytype/ztype = "min"/"max",
 *   or xval/yval/zval = value,
 *   to specify ElementGroup location.
 *  Specify 1, 2 or 3 coordinates for plane, line or point. 
 */


class EGroup
{
 public:

  static const char*        X_NAMES[3];
  static const char*        X_TYPES[5];
  static const char*        FROM_GRP_PROP;
  static const char*        ELEMS_PROP; 
  static const char*        ALL_PROP;
  static const char*        EPS_PROP; 

  enum  class               Check : int { NONE = 0, VALUE, BOUNDS };


                            EGroup                ();

  void                      configure

    ( const Properties&       conf,
      const Properties&       props,
      const String&           name );

  void                      makeGroup

      ( const ElementSet&     elems,
        const Properties&     globdat );

  IdxVector                 findElems
    
    ( const ElementSet&       elems, 
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat,
      const String&           name );


 private:

  void                      findElems_

    ( const ElementSet&       elems,
      const Properties&       globdat );

  void                      store_ 

    ( const ElementSet&       elems,
      const Properties&       globdat );


 private:

  // Every entry of doX_ is a switch that determines how an 
  // element is checked in that dimension (x,y,z).

  Matrix                    xbounds_;

  Array<Check>              doX_;  

  IdxVector                 ielems_;

  Tuple<Vector,3>           xvals_;
  Tuple<String,3>           xtype_;
  Tuple<bool,3>             completely_;

  StringVector              knownTypes_;
  String                    myName_;
  String                    parent_;

  double                    eps_;
  idx_t                     rank_;
  bool                      all_;

};


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  EGroup::X_NAMES[3]    = { "x", "y", "z" };
const char*  EGroup::X_TYPES[5]    = 
                         { "type", "val", "vals", "bounds", "complete" }; 
const char*  EGroup::FROM_GRP_PROP = "fromGroup";
const char*  EGroup::ELEMS_PROP    = "elems";
const char*  EGroup::ALL_PROP      = "all";
const char*  EGroup::EPS_PROP      = "eps";


//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


EGroup::EGroup

  ()

{
  xbounds_ .resize ( 3, 2 );  // ( x/y/z , min/max )

  doX_     .resize ( 3 );

  xvals_[0].resize ( 1 );
  xvals_[1].resize ( 1 );
  xvals_[2].resize ( 1 );

  xtype_      = "";
  completely_ = false;

  knownTypes_.resize ( 5 );

  knownTypes_[0] = "min";
  knownTypes_[1] = "max";
  knownTypes_[2] = "firstHalf";
  knownTypes_[3] = "secondHalf"; 
  knownTypes_[4] = "mid";

  parent_ = "all";

  eps_    = jem::Limits<double>::EPSILON;
  eps_    = 10.e-12;
  rank_   = 3_idx;
  all_    = false;
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void EGroup::configure

  ( const Properties&  conf,
    const Properties&  props,
    const String&      name )

{
  // Set group name.

  myName_ = name;

  // Get properties.

  Properties  myConf  = conf .makeProps ( myName_ );
  Properties  myProps = props. getProps ( myName_ );

  // Set elements.

  myProps.find ( parent_, FROM_GRP_PROP );
  myConf.set   ( FROM_GRP_PROP, parent_ );

  if ( myProps.find( ielems_, ELEMS_PROP ) )
  {
    myConf.set ( ELEMS_PROP, ielems_ );
    doX_ = Check::NONE;

    return;
  }

  myProps.find ( all_, ALL_PROP );
  myConf.set   ( ALL_PROP, all_ );

  if ( all_ )
  {
    return;
  }

  // Set tolerance.

  myProps.find ( eps_, EPS_PROP );
  myConf.set   ( EPS_PROP, eps_ );

  // Find coordinate types or values.

  for ( idx_t ir = 0; ir < rank_; ir++ )
  {
    Vector bounds;

    String thiss = X_NAMES[ir];

    String xtypename ( thiss + X_TYPES[0] );
    String xvalname  ( thiss + X_TYPES[1] );
    String xvalsname ( thiss + X_TYPES[2] );
    String xinname   ( thiss + X_TYPES[3] );
    String xcomplete ( thiss + X_TYPES[4] );

    doX_[ir] = Check::VALUE;

    if      ( myProps.find( xtype_[ir], xtypename ) )
    {
      myConf.set ( xtypename, xtype_[ir]);

      if ( ! jem::testany( knownTypes_ == xtype_[ir] ) )
      {
        throw jem::IllegalInputException ( 
          props.getContext ( myName_ ), 
          String::format   ( "unknown xtype: `%s'", xtype_[ir] ) 
        );
      }
    }
    else if ( myProps.find( xvals_[ir][0], xvalname ) )
    {
      myConf.set ( xvalname, xvals_[ir][0] );
    }
    else if ( myProps.find( xvals_[ir], xvalsname ) )
    {
      myConf.set ( xvalsname, xvals_[ir] );
    }
    else if ( myProps.find( bounds, xinname ) )
    {
      myConf.set ( xinname, bounds );

      if ( bounds.size() == 2 )
      {
        xbounds_(ir,0) = bounds[0] - eps_;
        xbounds_(ir,1) = bounds[1] + eps_;

        doX_[ir] = Check::BOUNDS;

        myProps.find ( completely_[ir], xcomplete );
        myConf.set   ( xcomplete, completely_[ir] );
      }
      else
      {
        throw jem::IllegalInputException ( 
          props.getContext ( myName_ ), 
          String::format   ( "`%s' does not have length 2!", xinname ) 
        );
      }
    }
    else
    {
      doX_[ir] = Check::NONE;
    }
  }
}


//-----------------------------------------------------------------------
//   makeGroup
//-----------------------------------------------------------------------


void EGroup::makeGroup

  ( const ElementSet&  elems,
    const Properties&  globdat )

{
  findElems_ ( elems, globdat );
  store_     ( elems, globdat );
}


//-----------------------------------------------------------------------
//   findElems
//-----------------------------------------------------------------------


IdxVector EGroup::findElems
    
  ( const ElementSet&  elems, 
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat,
    const String&      name )

{
  configure  ( conf, props, name );
  findElems_ ( elems, globdat );

  return ielems_;
}


//-----------------------------------------------------------------------
//   findElems_
//-----------------------------------------------------------------------


void EGroup::findElems_
 
  ( const ElementSet&  elems,
    const Properties&  globdat )

{
  using jem::util::ArrayBuffer;

  using jive::fem::ElementGroup;


  const idx_t  nodeCount = elems.maxElemNodeCount ();

  NodeSet             nodes   ( elems.getNodes() );

  Matrix              coords  ( nodes.rank(), nodeCount );

  Vector              xcoords ( elems.size() );

  IdxVector           inodes  ( nodeCount    );
  IdxVector           ieparent;

  ArrayBuffer<idx_t>  elembuf;

  bool                include;

  if ( all_ )
  {
    // Add all elements.

    elembuf.reserve ( elems.size() );

    for ( idx_t in = 0; in < elems.size(); in++ )
    {
      elembuf.pushBack ( in );
    }
  }
  else if ( jem::testany( doX_ ) )
  {
    // Find coordinate if specified with type.

    for ( idx_t ir = 0; ir < nodes.rank(); ir++ )
    {
      if ( jem::testany( xtype_[ir] == knownTypes_ ) )
      {
        nodes.getData()->getXCoords ( xcoords, ir );

        double minX = jem::min ( xcoords );
        double maxX = jem::max ( xcoords );

        if      ( xtype_[ir] == knownTypes_[0] ) // min
        {
          xvals_[ir] = minX;
        }
        else if ( xtype_[ir] == knownTypes_[1] ) // max
        {
          xvals_[ir] = maxX;
        }
        else if ( xtype_[ir] == knownTypes_[2] ) // firstHalf
        {
          xbounds_(ir,0) = minX - eps_;
          xbounds_(ir,1) = minX - eps_ + ( maxX - minX ) / 2.;

          doX_[ir] = Check::BOUNDS;
        }
        else if ( xtype_[ir] == knownTypes_[3] ) // secondHalf
        {
          xbounds_(ir,0) = minX + eps_ + ( maxX - minX ) / 2.;
          xbounds_(ir,1) = maxX + eps_;

          doX_[ir] = Check::BOUNDS;
        }
        else if ( xtype_[ir] == knownTypes_[4] ) // mid
        {
          xvals_[ir] = ( maxX + minX ) / 2.;
        }
      }
    }

    // Loop over elems, elem is included unless it violates the criterion
    // specified for one of the dimensions.

    if ( parent_ == "all" )
    {
      ieparent.ref ( IdxVector( jem::iarray ( elems.size() ) ) );
    }
    else
    {
      ElementGroup egroup ( 
        ElementGroup::get ( 
          parent_, elems, globdat, 
          "GroupInputModule search for parent" 
        ) 
      );

      ieparent.ref ( egroup.getIndices() );
    }

    for ( idx_t ie = 0; ie < ieparent.size(); ie++ )
    {
      include     = true;
      idx_t ielem = ieparent[ie];

      elems.getElemNodes  ( inodes, ielem  );
      nodes.getSomeCoords ( coords, inodes );

      // Loop over dimensions.

      for ( idx_t ir = 0; ir < nodes.rank(); ir++ )
      {
        double ixmin = jem::min( coords(ir,jem::ALL) );
        double ixmax = jem::max( coords(ir,jem::ALL) );

        if    ( doX_[ir] == Check::VALUE )
        {
          // Check if one of specified values is inside element.

          bool thisx = false;

          for ( idx_t ival = 0; ival < xvals_[ir].size(); ival++ )
          {
            thisx |= ( ixmin <= xvals_[ir][ival] && 
                       ixmax >= xvals_[ir][ival] );
          }

          include &= thisx;
        }
        else if ( doX_[ir] == Check::BOUNDS )
        {
          if ( completely_[ir] )
          {
            // Check if element is completely between specified bounds.

            include &= ( ixmin > xbounds_(ir,0) &&
                         ixmax < xbounds_(ir,1) );
          }
          else
          {
            // Check if element is partially between specified bounds.

            include &= ( ixmin < xbounds_(ir,1) &&
                         ixmax > xbounds_(ir,0) );
          }
        }
      }
      
      if ( include )
      {
        elembuf.pushBack ( ielem );
      }
    }
  }

  ielems_.ref ( elembuf.toArray() );
}
 

//-----------------------------------------------------------------------
//   store_
//-----------------------------------------------------------------------


void EGroup::store_
 
  ( const ElementSet&  elems,
    const Properties&  globdat )

{
  using jive::fem::ElementGroup;


  ElementGroup group = newElementGroup ( ielems_, elems );

  group.store ( myName_, globdat );

  jem::System::out() << "  ... Created ElementGroup `" << myName_ << 
    "' with " << ielems_.size() << " elements.\n";
  
  if ( ielems_.size() < 1 )
  {
    jem::System::warn() << "ElementGroup `" << myName_ << "' is empty.\n";
  }
}
 

//=======================================================================
//   class GroupInputModule
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  GroupInputModule::TYPE_NAME      = "groupInput";

const char*  GroupInputModule::NODE_GRPS_PROP = "nodeGroups";
const char*  GroupInputModule::ELEM_GRPS_PROP = "elemGroups";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


GroupInputModule::GroupInputModule

  ( const String&  name ) : Super ( name )

{}


GroupInputModule::~GroupInputModule

  ()

{}


//-----------------------------------------------------------------------
//   init 
//-----------------------------------------------------------------------


Module::Status GroupInputModule::init 

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  Properties  myConf  = conf .makeProps ( myName_ );
  Properties  myProps = props.findProps ( myName_ );

  // Read names of groups that are to be created.

  StringVector nGroupNames;
  StringVector eGroupNames;

  if ( myProps.find( nGroupNames, NODE_GRPS_PROP ) )
  {
    myConf.set ( NODE_GRPS_PROP, nGroupNames );
  }

  if ( myProps.find( eGroupNames, ELEM_GRPS_PROP ) )
  {
    myConf.set ( ELEM_GRPS_PROP, eGroupNames );
  }

  // Make NodeGroups.

  const NodeSet nodes = NodeSet::get 
                                 ( globdat, props.getContext( myName_ ) );

  for ( idx_t ing = 0; ing < nGroupNames.size(); ing++ )
  {
    NGroup ngrpInp;

    ngrpInp.configure ( myConf, myProps, nGroupNames[ing] );
    ngrpInp.makeGroup ( nodes , globdat ); 
  }

  // Make ElementGroups.

  const ElementSet& elems = ElementSet::get 
                                ( globdat, props.getContext( myName_ ) );

  for ( idx_t ieg = 0; ieg < eGroupNames.size(); ieg++ )
  {
    EGroup egrpInp;

    egrpInp.configure ( myConf, myProps, eGroupNames[ieg] );
    egrpInp.makeGroup ( elems , globdat ); 
  }

  return DONE;
} 


//-----------------------------------------------------------------------
//   findNodes
//-----------------------------------------------------------------------


IdxVector GroupInputModule::findNodes

  ( const NodeSet&     nodes,
    const Properties&  conf,
    const Properties&  props,
    const String&      name )  const

{
  Properties  myConf  = conf .makeProps ( myName_ );
  Properties  myProps = props.findProps ( myName_ );
 
  NGroup ngrpInp;

  return ngrpInp.findNodes ( nodes, myConf, myProps, name );
}


//-----------------------------------------------------------------------
//   findElems
//-----------------------------------------------------------------------


IdxVector GroupInputModule::findElems

  ( const ElementSet&  elems,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat,
    const String&      name )  const

{
  Properties  myConf  = conf .makeProps ( myName_ );
  Properties  myProps = props.findProps ( myName_ );
 
  EGroup egrpInp;

  return egrpInp.findElems ( elems, myConf, myProps, globdat, name ); 
}


//-----------------------------------------------------------------------
//   makeNew 
//-----------------------------------------------------------------------


Ref<Module> GroupInputModule::makeNew

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


void GroupInputModule::declare

  ()

{
  using jive::app::ModuleFactory;
  
  ModuleFactory::declare ( GroupInputModule::TYPE_NAME , & makeNew );
  ModuleFactory::declare ( GroupInputModule::CLASS_NAME, & makeNew );
}
