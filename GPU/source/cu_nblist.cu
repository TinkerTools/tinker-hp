
#include "utils.h"
#include "image.h"

extern const real eps_cell;
extern const real _xcell;
extern const real _ycell;
extern const real _zcell;
extern const real _box34;
extern const int octahedron;

#define FILTER_PARAMS_GEN \
        const int* cell_glob, const int* block_scan, const real* x, const real* y, const real* z, const int* matb_lst,       \
        const int nlocnlb, const int nblock, const int matsize, const int nblocknlb_pair, const real cutbuff2,               \
        int* lst, int* nlst
#define FILTER_PARAMS     \
        FILTER_PARAMS_GEN , TINKER_IMAGE_PARAMS
#define FILTER1_PARAMS    \
        FILTER_PARAMS_GEN , TINKER_IMAGE1_PARAMS

#define FILTER_ARGS       \
        cell_glob, block_scan, x, y, z, \
        matb_lst, nlocnlb, nblock, matsize, nblocknlb_pair, cutbuff2,          \
        lst, nlst , TINKER_IMAGE_ARGS

__global__
void cu_filter_lst_sparse (FILTER1_PARAMS)
{
   const int ithread = threadIdx.x + (blockIdx.x)*blockDim.x;
   const int iwarp   =     ithread / WARP_SIZE;
   const int nwarp   = (blockDim.x*gridDim.x) / WARP_SIZE;
   const int ilane   = threadIdx.x & (WARP_SIZE-1);

   int ii,j;

   //if (ithread==0) printf("%d %d %d %d %lf %lf %lf",nlocnlb,nblock,matsize,nblocknlb_pair,cutbuff2,xcell,ycell2);

   for ( ii=iwarp; ii<nblocknlb_pair; ii+=nwarp ){
      int iblock = lst[2*ii]-1;
      int idx   = (iblock)*WARP_SIZE + ilane;
      int iscan = block_scan[iblock]*WARP_SIZE + 2*nblocknlb_pair;
      real xi   = x[idx];
      real yi   = y[idx];
      real zi   = z[idx];

      int kdx   = (lst[2*ii+1]-1)*WARP_SIZE + ilane;
      //int kglob = cell_glob[kdx];
      real xk   = x[kdx];
      real yk   = y[kdx];
      real zk   = z[kdx];

      int jth   = WARP_SIZE-1;
      //int bset  = 0;
      int nbits = 0;
      int oldnb = 0;

      if (idx==kdx) lst [iscan + ilane] = kdx+1;
      else {
         for ( j=0;j<WARP_SIZE;j++ ){
            int srclane = j;
            real xpos = xi - __shfl_sync(ALL_LANES,xk,srclane);
            real ypos = yi - __shfl_sync(ALL_LANES,yk,srclane);
            real zpos = zi - __shfl_sync(ALL_LANES,zk,srclane);
            Image (xpos,ypos,zpos);
            int ilane_set = (xpos*xpos+ypos*ypos+zpos*zpos<=cutbuff2) ? 1 : 0;
            int accept = __ballot_sync(ALL_LANES,ilane_set);
            if (accept) {
               if (nbits==ilane) jth = j;
               nbits += 1;
            }
         }

         if (ilane==0) oldnb = atomicAdd(nlst+iblock,nbits);

         int kdx_ = __shfl_sync(ALL_LANES,kdx,jth);
         oldnb    = __shfl_sync(ALL_LANES,oldnb,0);

         if (ilane<nbits) lst[iscan+oldnb+ilane] = kdx_+1;
      }
   }
}

EXTERN_C_BEG
void filter_lst_sparse (FILTER_PARAMS,cudaStream_t st){

   const int gS  = 320;
   const int bS  = 4*WARP_SIZE;
   const int sh  = 0;

   cu_filter_lst_sparse<<<gS,bS,sh,st>>> (FILTER_ARGS);
   cudaError_t ierrSync = cudaGetLastError();
   //cudaError_t ierrSync = cudaDeviceSynchronize();
   if (ierrSync != cudaSuccess) {
      printf("nblist C filter kernel error: %d \n %s",ierrSync,
             cudaGetErrorString(ierrSync));
   }
   return;
}
EXTERN_C_END
