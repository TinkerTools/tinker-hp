#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/remove.h>
#include "thrust_cache.h"


extern "C" {
namespace THRUST_WRAPPER{
   // sorting wrappers
   void sort_int_wrapper    ( int*    data, int N )
   {
      thrust::device_ptr<int>    dev_ptr(data);
      thrust::sort(dev_ptr, dev_ptr+N);
   }
   void sort_int_async_wrapper ( int* data, int N, cudaStream_t s)
   {
      thrust::device_ptr<int>    dev_ptr(data);
      thrust::sort(thrust::cuda::par(thrust_cache).on(s), dev_ptr, dev_ptr+N);
   }
   void sort_float_wrapper  ( float*  data, int N )
   {
      thrust::device_ptr<float>  dev_ptr(data);
      thrust::sort(dev_ptr, dev_ptr+N);
   }
   void sort_double_wrapper ( double* data, int N )
   {
      thrust::device_ptr<double> dev_ptr(data);
      thrust::sort(dev_ptr, dev_ptr+N);
   }

   // inclusive scan wrappers
   void inclusive_scan_int_wrapper ( int* in, int N, int* out )
   {
      thrust::device_ptr<int> dev_in(in);
      thrust::device_ptr<int> dev_out(out);
      thrust::inclusive_scan (dev_in,dev_in+N,dev_out);
   }

   // exclusive scan wrappers
   void exclusive_scan_int_wrapper ( int* in, int N, int* out )
   {
      thrust::device_ptr<int> dev_in(in);
      thrust::device_ptr<int> dev_out(out);
      thrust::exclusive_scan (dev_in,dev_in+N,dev_out);
   }
   void exclusive_scan_int_async_wrapper ( int* in, int N, int* out, cudaStream_t s)
   {
      thrust::device_ptr<int> dev_in(in);
      thrust::device_ptr<int> dev_out(out);
      thrust::exclusive_scan (thrust::cuda::par(thrust_cache).on(s),dev_in,dev_in+N,dev_out);
   }

   // Stable sort wrappers
   void stable_sort_by_key_int_wrapper ( int* in, int N, int* values )
   {
      thrust::device_ptr<int> dev_in(in);
      thrust::device_ptr<int> dev_values(values);
      thrust::stable_sort_by_key (dev_in,dev_in+N,dev_values);
   }
   void stable_sort_by_key_int_async_wrapper ( int* in, int N, int* values, cudaStream_t s )
   {
      thrust::device_ptr<int> dev_in(in);
      thrust::device_ptr<int> dev_values(values);
      thrust::stable_sort_by_key (thrust::cuda::par(thrust_cache).on(s),dev_in,dev_in+N,dev_values);
   }

   // Remove wrappers
   void remove_zero_wrapper (int* in, int N)
   {
      thrust::device_ptr<int> dev_in(in);
      thrust::remove(dev_in,dev_in+N,0);
   }
   void remove_zero_async_wrapper (int* in, int N, cudaStream_t s)
   {
      thrust::device_ptr<int> dev_in(in);
      thrust::remove(thrust::cuda::par(thrust_cache).on(s),dev_in,dev_in+N,0);
   }
   void remove_async_wrapper (int* in, int N, int value, int* new_last, cudaStream_t s)
   {
      thrust::device_ptr<int> dev_in(in);
      thrust::device_ptr<int> dev_last(
            thrust::remove(thrust::cuda::par(thrust_cache).on(s),dev_in,dev_in+N,value));
      *new_last = dev_last - dev_in;
   }
}
}
