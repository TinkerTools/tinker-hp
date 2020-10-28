#include <stdint.h>
#include <stdio.h>
#include "thrust_cache.h"

namespace THRUST_WRAPPER{
   typedef char value_type;
   ThrustCache::ThrustCache()
      : nbytes(0),ext(0)
      , ptr(nullptr)
   {}


   value_type* ThrustCache::attach_memory( void** ptr_arg, ptrdiff_t numbyte)
   {
      if (numbyte > nbytes) {
         nbytes = numbyte;
         ext    = 1;
         ptr    = reinterpret_cast<value_type*>(*ptr_arg);
         //printf( "attaching memory  %lu %lu %ld %d", nbytes, numbyte, (long int)ptr, ext );
      }
      return ptr;
   }


   void ThrustCache::deallocate(value_type*, size_t)
   {
      // does not do anything
      return;
   }


   void ThrustCache::clear()
   {
      ptr = nullptr;
      nbytes = 0;
   }


   ThrustCache thrust_cache;


   extern "C" {
      void thrust_cache_dealloc()
      {
         thrust_cache.clear();
      }


      void thrust_cache_alloc(void** ptr,intptr_t numbytes)
      {
         thrust_cache.attach_memory(ptr,(ptrdiff_t)numbytes);
      }
   }
}
