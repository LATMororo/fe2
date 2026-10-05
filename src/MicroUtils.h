
#ifndef MICRO_UTILS_H
#define MICRO_UTILS_H


#include <jem/base/array/Array.h>
#include <jem/base/Object.h>


//#include <jive/Array.h>


//=======================================================================
//   namespace MicroUtils 
//=======================================================================


namespace MicroUtils
{ 
  //=====================================================================
  //   Micro model actions and related params 
  //=====================================================================


  class MicroActions
  {
   public:
   
    static const char*      GET_MACRO_STRAIN; // update macro strain 
    static const char*      GET_MICRO_STATE;  // update micro state

  };

  class MicroActionParams
  {
   public:

    static const char*      MACRO_STRAIN;
    static const char*      MICRO_STATE;

  };


  //=====================================================================
  //   Micro model state
  //=====================================================================


  class MicroState : public jem::Object
  {
   public:

    jem::Array<double,2>    stiff;
    jem::Array<double>      stress;

  };

}


#endif
