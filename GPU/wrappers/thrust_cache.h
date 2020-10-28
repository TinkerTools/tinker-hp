#pragma once
#include <thrust/system/cuda/vector.h>
#include <cstddef>
#include <iostream>
#include <stdio.h>

namespace THRUST_WRAPPER{
class ThrustCache
{
public:
   typedef char value_type;

   ThrustCache();
   ~ThrustCache(){};

   value_type* allocate(ptrdiff_t num_bytes){
      //printf(" calling allocator thrust with %lu \n at disposal %lu \n ",num_bytes, nbytes);
      if (num_bytes>nbytes){
         if (ext) {
            std::cout << "Allocated thrust cache memory is not enough to perform the routine;" <<std::endl;
            std::cout << "Increase initial thrust cache" <<std::endl <<"Switching to thrust regular cache allocator" <<std::endl;

            // no allocation of the right size exists
            // create a new one with cuda::malloc
            // try if cuda::malloc can't satisfy the request
            try { ptr  = thrust::cuda::malloc<char>(num_bytes).get(); }
            catch(std::runtime_error &e) { throw; }
            ext = 0;
         }
         else {
            // transform the pointer to cuda::pointer before calling cuda::free
            thrust::cuda::free(thrust::cuda::pointer<char>(ptr)); 
            // reallocate pointer with th appropriate suze
            try { ptr  = thrust::cuda::malloc<char>(num_bytes).get(); }
            catch(std::runtime_error &e) { throw; }
         }
         return ptr;
      }
      else return ptr;
   }
   value_type* attach_memory(void**,ptrdiff_t);
   void deallocate(value_type*, size_t);
   void clear();

private:
   size_t nbytes;
   int ext;
   value_type* ptr;
};
extern ThrustCache thrust_cache;
}

