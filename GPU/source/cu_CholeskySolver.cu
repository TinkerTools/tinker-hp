#include "utils.h"

#if (defined(SINGLE)||defined(MIXED))
#   define cuPOTRF_buffSize cusolverDnSpotrf_bufferSize
#   if (CUDART_VERSION>10010)
#     define cuGESV_buffSize cusolverDnSSgesv_bufferSize
#     define cuGESV           cusolverDnSSgesv
#   else
#     define cuGETRF_buffsize cusolverDnSgetrf_bufferSize
#     define cuGETRF          cusolverDnSgetrf
#     define cuGETRS          cusolverDnSgetrs
#   endif
#   define cuPOTRF          cusolverDnSpotrf
#   define cuPOTRS          cusolverDnSpotrs
#else
#   define cuPOTRF_buffSize cusolverDnDpotrf_bufferSize
#   if (CUDART_VERSION>10010)
#     define cuGESV_buffSize  cusolverDnDDgesv_bufferSize
#     define cuGESV           cusolverDnDDgesv
#   else
#     define cuGETRF_buffsize cusolverDnDgetrf_bufferSize
#     define cuGETRF         cusolverDnDgetrf
#     define cuGETRS         cusolverDnDgetrs
#   endif
#   define cuPOTRF          cusolverDnDpotrf
#   define cuPOTRS          cusolverDnDpotrs
#endif

extern const int rank;
extern const int tinkerdebug;

/* ---------
   Cu Solver global environnement
   ---------
*/
cusolverDnHandle_t cuCholHandle = NULL;
const cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
real* d_workSpace=NULL;
size_t s_workSpaceSize;
int info;
int* d_info;

__global__ void CheckcuSolverInfo (int* d_info, int line, int rank) {
  if (*d_info != 0) printf (" Error info %d with cuSolver in " __FILE__ " :%d \n", *d_info, line);
}

EXTERN_C_BEG

void initcuSolverHandle(cudaStream_t stream){
   gpuErrchkSolver( cusolverDnCreate(&cuCholHandle) )

   gpuErrchkSolver( cusolverDnSetStream(cuCholHandle, stream) )

   gpuErrchk( cudaMalloc (&d_info, sizeof(int)) )

   if (rank==0) printf ("\n *** Using CuSolver Library ***\n\n" );
}

/* ----------------
   Reallocation procedure based on MOD_utilgpu.f reallocate_acc
   ---------------- */
void device_reallocate(void** array, const size_t bytesSize, size_t& PrevSize){
   if ( !(*array) ){
      gpuErrchk( cudaMalloc(array,bytesSize) )
      PrevSize = bytesSize;
      //printf(" device_reallocate size %lu \n", PrevSize);
   }
   else {
      if (bytesSize > PrevSize) {
         gpuErrchk( cudaFree( *array ) )
         gpuErrchk( cudaMalloc(array,bytesSize) )
         PrevSize = bytesSize;
         /*printf(" device_reallocate size %lu \n", PrevSize);*/
      }
   }
}

void cuPOTRF_Wrapper(const int n, real* A, const int lda, cudaStream_t stream){
   int Lwork=0;
   cusolverStatus_t status1;

   gpuErrchk( cuPOTRF_buffSize(cuCholHandle, uplo, n, A, lda, &Lwork) )

   device_reallocate((void**)&d_workSpace, (size_t)Lwork*sizeof(real), s_workSpaceSize);

   status1 = cuPOTRF(cuCholHandle,uplo, n, A, lda, d_workSpace, Lwork, d_info);
   if (status1!=CUSOLVER_STATUS_SUCCESS) printf( "Cholesky Factorisation on device failed with Error %d \n",status1 );

   if (tinkerdebug) {
      CheckcuSolverInfo<<<1,1,0,stream>>>(d_info, __LINE__, rank);
      gpuErrchk( cudaGetLastError() )
   }

}

void cuPOTRS_Wrapper(const int n, real* A, const int lda, real* B, const int ldb, cudaStream_t stream){
   cusolverStatus_t status1;

   status1 = cuPOTRS(cuCholHandle, uplo, n, 1, A, lda, B, ldb, d_info);
   if (status1!=CUSOLVER_STATUS_SUCCESS) printf( "Error %d solving Linear system \n",status1);

   if (tinkerdebug) {
      CheckcuSolverInfo<<<1,1,0,stream>>>(d_info, __LINE__, rank);
      gpuErrchk( cudaGetLastError() )
   }
}

void cuGESV_Wrapper(const int n, const int nrhs, real* A, const int lda, int* Ipiv, real* B, const int ldb, cudaStream_t stream){
#if (CUDART_VERSION>10010)
   size_t lwork_bytes=0;
   int iter;
   gpuErrchkSolver( cuGESV_buffSize(cuCholHandle, n, nrhs, A, lda, Ipiv, B, ldb, B, ldb, d_workSpace, &lwork_bytes) )
   device_reallocate((void**)&d_workSpace, (size_t)lwork_bytes, s_workSpaceSize);

   gpuErrchkSolver( cuGESV(cuCholHandle, n, nrhs, A, lda, Ipiv, B, ldb, B, ldb, d_workSpace, lwork_bytes, &iter, d_info) )
#else
   int Lwork=0;
   gpuErrchkSolver( cuGETRF_buffsize(cuCholHandle, n, n, A, lda, &Lwork) )
   //printf(" LU solve n %d nrhs %d lda %d ldb %d Lwork %d\n",n,nrhs,lda,ldb, Lwork);
   device_reallocate((void**)&d_workSpace, (size_t)Lwork*sizeof(real), s_workSpaceSize);

   gpuErrchkSolver( cuGETRF(cuCholHandle, n, n, A, lda, d_workSpace, Ipiv, d_info) )
   gpuErrchkSolver( cuGETRS(cuCholHandle, CUBLAS_OP_N, n, 1, A, lda, Ipiv, B, ldb, d_info) )
#endif
   if (tinkerdebug) {
      CheckcuSolverInfo<<<1,1,0,stream>>>(d_info, __LINE__, rank);
      gpuErrchk( cudaGetLastError() )
   }
}

void destroycuSolverHandle(){
   gpuErrchkSolver( cusolverDnDestroy(cuCholHandle) )
}
EXTERN_C_END
