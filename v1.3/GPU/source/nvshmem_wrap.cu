#include "utils.h"
#include "mpi.h"
#include "nvshmem.h"
#include "nvshmemx.h"

EXTERN_C_BEG

nvshmemx_init_attr_t attr;

__global__ static void zero_put (char* ptr, const size_t n ){
   for (int i=threadIdx.x + blockDim.x*blockIdx.x; i<n; i+=blockDim.x*gridDim.x)
      ptr[i] = 0x00;
}

void nvshmem_init_f ( int *mype, int *npes, int *comm ){
   MPI_Comm c_comm = MPI_Comm_f2c( (MPI_Fint) *comm );

   attr.mpi_comm = &c_comm;
   nvshmemx_init_attr (NVSHMEMX_INIT_WITH_MPI_COMM, &attr);
   *mype = nvshmem_my_pe();
   *npes = nvshmem_n_pes();

   if (*mype==0) printf("\n *****  Using NVShmem Feature *****\n\n");
}

void nvshmem_finalize_f() { nvshmem_finalize(); }

void* nvshmem_malloc_f( size_t size ){ return nvshmem_malloc(size); }

void* nvshmem_calloc_f( size_t size ){
   /*  Allocation */
   void* ptr = nvshmem_malloc(size);
   /*  initiate with zero */
   zero_put<<<(size>>8),1<<8,0>>>((char*)ptr,size);
   return ptr;
}

void* nvshmem_prt_f( void* ptr, int pe ){ return nvshmem_ptr(ptr, pe); }

EXTERN_C_END
