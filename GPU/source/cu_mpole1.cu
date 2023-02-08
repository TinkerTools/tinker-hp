#include "utils.h"
#include "image.h"


extern const int nproc;
extern const int nrank;
extern const int tinkerdebug;
extern const cudaDeviceProp devProp;
extern const real eps_cell;
extern const real _xcell;
extern const real _ycell;
extern const real _zcell;
extern const real _ixcell;
extern const real _iycell;
extern const real _izcell;
extern const real _box34;
extern const int octahedron;


template<int do_correct>
__device__
void mpole1_couple(const real r2, const real& xr, const real& yr, const real& zr, const rpole_elt& ip, const rpole_elt& kp, const real c_mscale, const real& aewald, const real& f, const real& alsq2n, const real& alsq2,
                   real& e, real3& frc, real3m& frc_i, real3m& frc_k, real3& ttmi, real3& ttmk){
#pragma acc routine
   real de;
   real r,invr,invr2;
   real rr1,rr3,rr5,rr7,rr9,rr11;
   real bn0,bn1,bn2,bn3,bn4,bn5;
   real alsqt,ralpha,exp2a;
   real dikx,diky,dikz;
   real qrix,qriy,qriz;
   real qrkx,qrky,qrkz;
   real qrixr,qriyr,qrizr;
   real qrkxr,qrkyr,qrkzr;
   real qrrx,qrry,qrrz;
   real qikrx,qikry,qikrz;
   real qkirx,qkiry,qkirz;
   real diqkx,diqky,diqkz;
   real dkqix,dkqiy,dkqiz;
   real dqiqkx,dqiqky,dqiqkz;
   real dri,drk,qrri,qrrk;
   real diqrk,dkqri;
   real qik,qrrik;
   real term1,term2,term3;
   real term4,term5,term6;

   /* get reciprocal distance terms for this interaction */

   r      = f_sqrt(r2);
   invr   = 1.0/r;
   invr2  = 1.0/r2;
   rr1    = c_mscale*f * invr;
   rr3    =         rr1 * invr2;
   rr5    = 3.0   * rr1 * invr2 * invr2;
   rr7    = 15.0  * rr1 * invr2 * invr2 * invr2;
   rr9    = 105.0 * rr3 * invr2 * invr2 * invr2;
   rr11   = 63.0  * rr7 * invr2 * invr2;

   /* calculate the real space Ewald error function terms */

   if (do_correct) {
      bn0 = 0.0; bn1 = 0.0; bn2 = 0.0;
      bn3 = 0.0; bn4 = 0.0; bn5 = 0.0;  }
   else {
      ralpha = aewald * r;
      alsqt  = alsq2*alsq2;
      exp2a  = f_exp(-ralpha*ralpha);
      bn0    = f_erfc(ralpha);
      bn0    = bn0*invr;

      bn1    = (    bn0+alsq2            *alsq2n*exp2a) * invr2;
      bn2    = (3.0*bn1+alsqt            *alsq2n*exp2a) * invr2;
      bn3    = (5.0*bn2+alsqt*alsq2      *alsq2n*exp2a) * invr2;
      bn4    = (7.0*bn3+alsqt*alsqt      *alsq2n*exp2a) * invr2;
      bn5    = (9.0*bn4+alsqt*alsqt*alsq2*alsq2n*exp2a) * invr2;

      bn0    = f * bn0;
      bn1    = f * bn1;
      bn2    = f * bn2;
      bn3    = f * bn3;
      bn4    = f * bn4;
      bn5    = f * bn5;
   }

  /* modify distances to account for Ewald and exclusions */

   rr1    = bn0 - rr1;
   rr3    = bn1 - rr3;
   rr5    = bn2 - rr5;
   rr7    = bn3 - rr7;
   rr9    = bn4 - rr9;
   rr11   = bn5 - rr11;

   /* intermediates involving moments and distance separation */

   dikx   = ip.dy*kp.dz - ip.dz*kp.dy;
   diky   = ip.dz*kp.dx - ip.dx*kp.dz;
   dikz   = ip.dx*kp.dy - ip.dy*kp.dx;

   dri    = ip.dx*xr  + ip.dy*yr  + ip.dz*zr;
   drk    = kp.dx*xr  + kp.dy*yr  + kp.dz*zr;
   qrix   = ip.qxx*xr + ip.qxy*yr + ip.qxz*zr;
   qriy   = ip.qxy*xr + ip.qyy*yr + ip.qyz*zr;
   qriz   = ip.qxz*xr + ip.qyz*yr + ip.qzz*zr;
   qrkx   = kp.qxx*xr + kp.qxy*yr + kp.qxz*zr;
   qrky   = kp.qxy*xr + kp.qyy*yr + kp.qyz*zr;
   qrkz   = kp.qxz*xr + kp.qyz*yr + kp.qzz*zr;
   qrri   = qrix*xr + qriy*yr + qriz*zr;
   qrrk   = qrkx*xr + qrky*yr + qrkz*zr;
   qrrik  = qrix*qrkx + qriy*qrky + qriz*qrkz;

   qik    = 2.0*(ip.qxy*kp.qxy + ip.qxz*kp.qxz + ip.qyz*kp.qyz)
               + ip.qxx*kp.qxx + ip.qyy*kp.qyy + ip.qzz*kp.qzz;

   qrixr  = qriz*yr - qriy*zr;
   qriyr  = qrix*zr - qriz*xr;
   qrizr  = qriy*xr - qrix*yr;
   qrkxr  = qrkz*yr - qrky*zr;
   qrkyr  = qrkx*zr - qrkz*xr;
   qrkzr  = qrky*xr - qrkx*yr;

   qrrx   = qrky*qriz - qrkz*qriy;
   qrry   = qrkz*qrix - qrkx*qriz;
   qrrz   = qrkx*qriy - qrky*qrix;

   qikrx  = ip.qxx*qrkx + ip.qxy*qrky + ip.qxz*qrkz;
   qikry  = ip.qxy*qrkx + ip.qyy*qrky + ip.qyz*qrkz;
   qikrz  = ip.qxz*qrkx + ip.qyz*qrky + ip.qzz*qrkz;
   qkirx  = kp.qxx*qrix + kp.qxy*qriy + kp.qxz*qriz;
   qkiry  = kp.qxy*qrix + kp.qyy*qriy + kp.qyz*qriz;
   qkirz  = kp.qxz*qrix + kp.qyz*qriy + kp.qzz*qriz;

   diqkx  = ip.dx*kp.qxx  + ip.dy*kp.qxy + ip.dz*kp.qxz;
   diqky  = ip.dx*kp.qxy  + ip.dy*kp.qyy + ip.dz*kp.qyz;
   diqkz  = ip.dx*kp.qxz  + ip.dy*kp.qyz + ip.dz*kp.qzz;
   dkqix  = kp.dx*ip.qxx  + kp.dy*ip.qxy + kp.dz*ip.qxz;
   dkqiy  = kp.dx*ip.qxy  + kp.dy*ip.qyy + kp.dz*ip.qyz;
   dkqiz  = kp.dx*ip.qxz  + kp.dy*ip.qyz + kp.dz*ip.qzz;
   diqrk  = ip.dx*qrkx    + ip.dy*qrky   + ip.dz*qrkz;
   dkqri  = kp.dx*qrix    + kp.dy*qriy   + kp.dz*qriz;

   dqiqkx = ip.dy*qrkz - ip.dz*qrky +  kp.dy*qriz - kp.dz*qriy
          - 2* (ip.qxy*kp.qxz + ip.qyy*kp.qyz + ip.qyz*kp.qzz
              - ip.qxz*kp.qxy - ip.qyz*kp.qyy - ip.qzz*kp.qyz);
   dqiqky = ip.dz*qrkx - ip.dx*qrkz + kp.dz*qrix - kp.dx*qriz
          - 2*(ip.qxz*kp.qxx + ip.qyz*kp.qxy + ip.qzz*kp.qxz
             - ip.qxx*kp.qxz - ip.qxy*kp.qyz - ip.qxz*kp.qzz);
   dqiqkz = ip.dx*qrky - ip.dy*qrkx + kp.dx*qriy - kp.dy*qrix
          - 2*(ip.qxx*kp.qxy + ip.qxy*kp.qyy + ip.qxz*kp.qyz
             - ip.qxy*kp.qxx - ip.qyy*kp.qxy - ip.qyz*kp.qxz);

   /* calculate intermediate terms for multipole energy */

   term1  = ip.c*kp.c;
   term2  = kp.c*dri   - ip.c*drk   + ip.dx*kp.dx + ip.dy*kp.dy + ip.dz*kp.dz;
   term3  = ip.c*qrrk  + kp.c*qrri  - dri*drk + 2*(dkqri - diqrk + qik);
   term4  = dri*qrrk   - drk*qrri   - 4.0*qrrik;
   term5  = qrri*qrrk;

   /* compute the energy contributions for this interaction */
   e     += term1*rr1 + term2*rr3 + term3*rr5 + term4*rr7 + term5*rr9;
   de     = term1*rr3 + term2*rr5 + term3*rr7 + term4*rr9 + term5*rr11;

   /* calculate intermediate terms for force and torque */
   term1  = -kp.c*rr3 + drk*rr5 - qrrk*rr7;
   term2  =  ip.c*rr3 + dri*rr5 + qrri*rr7;
   term3  = 2.0 * rr5;
   term4  = 2.0 * (-kp.c*rr5+drk*rr7-qrrk*rr9);
   term5  = 2.0 * (-ip.c*rr5-dri*rr7-qrri*rr9);
   term6  = 4.0 * rr7;

   /* compute the force components for this interaction */
   frc.x  = de*xr 
          + term1*ip.dx         + term2*kp.dx
          + term3*(diqkx-dkqix) + term4*qrix
          + term5*qrkx          + term6*(qikrx+qkirx);
   frc.y  = de*yr 
          + term1*ip.dy         + term2*kp.dy
          + term3*(diqky-dkqiy) + term4*qriy
          + term5*qrky          + term6*(qikry+qkiry);
   frc.z  = de*zr 
          + term1*ip.dz         + term2*kp.dz
          + term3*(diqkz-dkqiz) + term4*qriz
          + term5*qrkz          + term6*(qikrz+qkirz);

   /* increment the compute the torque components for this interaction */

   ttmi.x += -rr3*dikx + term1*(ip.dy*zr - ip.dz*yr)  + term3*(dqiqkx+dkqiz*yr  - dkqiy*zr) - term4*qrixr - term6*(qikrz*yr  - qikry*zr + qrrx);
   ttmi.y += -rr3*diky + term1*(ip.dz*xr - ip.dx*zr)  + term3*(dqiqky+dkqix*zr  - dkqiz*xr) - term4*qriyr - term6*(qikrx*zr  - qikrz*xr + qrry);
   ttmi.z += -rr3*dikz + term1*(ip.dx*yr - ip.dy*xr)  + term3*(dqiqkz+dkqiy*xr  - dkqix*yr) - term4*qrizr - term6*(qikry*xr  - qikrx*yr + qrrz);
   ttmk.x +=  rr3*dikx + term2*(kp.dy*zr - kp.dz*yr)  - term3*(dqiqkx+diqkz*yr  - diqky*zr) - term5*qrkxr - term6*(qkirz*yr  - qkiry*zr - qrrx);
   ttmk.y +=  rr3*diky + term2*(kp.dz*xr - kp.dx*zr)  - term3*(dqiqky+diqkx*zr  - diqkz*xr) - term5*qrkyr - term6*(qkirx*zr  - qkirz*xr - qrry);
   ttmk.z +=  rr3*dikz + term2*(kp.dx*yr - kp.dy*xr)  - term3*(dqiqkz+diqky*xr  - diqkx*yr) - term5*qrkzr - term6*(qkiry*xr  - qkirx*yr - qrrz);

   /* store in large container for mixed precision */
   frc_i.x += frc.x;
   frc_i.y += frc.y;
   frc_i.z += frc.z;
   frc_k.x -= frc.x;
   frc_k.y -= frc.y;
   frc_k.z -= frc.z;
}


#define MPOLE1_PARAMS_GEN                                              \
        const int*restrict ipole, const int*restrict pglob, const int*restrict loc, const int*restrict ieblst, const int*restrict eblst   \
        , const real*restrict x, const real*restrict y, const real*restrict z, const real (*restrict rpole)[13]  \
        , realm (*restrict dem)[3], real (*restrict tem)[3]                                                      \
        , realm *restrict em_buffer, real *restrict vir_buffer                                                   \
        , const int npolelocnlb, const int npolelocnlb_pair, const int npolebloc, const int n                    \
        , const real off2, const real f, const real alsq2, const real alsq2n, const real aewald
#define MPOLE1_PARAMS                                       \
        MPOLE1_PARAMS_GEN , MIDPOINTIMAGE_PARAMS  
#define MPOLE1_PARAMS1                                      \
        MPOLE1_PARAMS_GEN , MIDPOINTIMAGE1_PARAMS

#define MPOLE1_ARGS                                         \
          ipole, pglob, loc, ieblst, eblst, x, y, z, rpole  \
        , dem, tem, em_buffer, vir_buffer                   \
        , npolelocnlb, npolelocnlb_pair, npolebloc, n       \
        , off2, f, alsq2, alsq2n, aewald                    \
        , MIDPOINTIMAGE_ARGS, nproc

#define GROUPS_PARAMS \
        const int ngrp, const int use_group \
        , real (*restrict wgrp) \
        , int (*restrict grplist)
#define GROUPS_ARGS \
        ngrp, use_group, wgrp, grplist

__global__
void cu_emreal1c_core (MPOLE1_PARAMS1,const int nproc,GROUPS_PARAMS){

   const int ithread = threadIdx.x + blockIdx.x*blockDim.x;
   const int iwarp   =               ithread / WARP_SIZE;
   const int nwarp   =  blockDim.x*gridDim.x / WARP_SIZE;
   const int ilane   = threadIdx.x & (WARP_SIZE-1);
   const int beg     = 0;
   int accept_mid    = 1;

   int klane,srclane;
   int ii,j;
   int iblock,idx,kdx;
   int iipole,iglob,i,kpole,kglob,kbis,kglob_;
   int do_pair,same_block;
   //int iga,igb,igb_;
   real xk_,yk_,zk_,d2;
   real em_;
   real vir_[6];
   rpole_elt ip;
   __shared__ rpole_elt kp[BLOCK_DIM];
   real3 posi,pos;
   real3 ttmi,frc;
   real3m frc_i;
   __shared__ real3  posk[BLOCK_DIM],ttmk[BLOCK_DIM];
   __shared__ real3m frc_k[BLOCK_DIM];
   real scale_;

   //if (ithread==0) printf( " %i %i %i %i %i %i %i " r_Format r_Format r_Format "\n", nwarp,RED_BUFF_SIZE,nproc,npolelocnlb,npolebloc,n,npolelocnlb_pair,off2,alsq2,aewald);

   for ( ii=iwarp; ii<npolelocnlb_pair; ii+=nwarp ){
      /*  Load atom block i parameters */
      iblock  = ieblst[ii];
      if (iblock==0) continue;
      idx     = (iblock-1)*WARP_SIZE + ilane;
      iipole  = pglob[idx] -1;
      iglob   = ipole[idx] -1;
      i       = loc  [idx] -1;
      posi.x  = x[idx];
      posi.y  = y[idx];
      posi.z  = z[idx];
      ip.c    = rpole[iipole][0];
      ip.dx   = rpole[iipole][1];
      ip.dy   = rpole[iipole][2];
      ip.dz   = rpole[iipole][3];
      ip.qxx  = rpole[iipole][4];
      ip.qxy  = rpole[iipole][5];
      ip.qxz  = rpole[iipole][6];
      ip.qyy  = rpole[iipole][8];
      ip.qyz  = rpole[iipole][9];
      ip.qzz  = rpole[iipole][12];

      /*  Load atom block k parameters */
      kdx     = eblst[beg+ ii*WARP_SIZE+ ilane] -1;
      kpole   = pglob[kdx] -1;
      kglob   = ipole[kdx] -1;
      kbis    = loc  [kdx] -1;
      posk[threadIdx.x].x  = x[kdx];
      posk[threadIdx.x].y  = y[kdx];
      posk[threadIdx.x].z  = z[kdx];
      kp[threadIdx.x].c    = rpole[kpole][0];
      kp[threadIdx.x].dx   = rpole[kpole][1];
      kp[threadIdx.x].dy   = rpole[kpole][2];
      kp[threadIdx.x].dz   = rpole[kpole][3];
      kp[threadIdx.x].qxx  = rpole[kpole][4];
      kp[threadIdx.x].qxy  = rpole[kpole][5];
      kp[threadIdx.x].qxz  = rpole[kpole][6];
      kp[threadIdx.x].qyy  = rpole[kpole][8];
      kp[threadIdx.x].qyz  = rpole[kpole][9];
      kp[threadIdx.x].qzz  = rpole[kpole][12];

      /* set compute Data to]0 */
      em_ = 0.0;
      vir_[0]=0.0; vir_[1]=0.0; vir_[2]=0.0;
      vir_[3]=0.0; vir_[4]=0.0; vir_[5]=0.0;

      frc_i.x = 0.0;
      frc_i.y = 0.0;
      frc_i.z = 0.0;
      frc_k[threadIdx.x].x = 0.0;
      frc_k[threadIdx.x].y = 0.0;
      frc_k[threadIdx.x].z = 0.0;
      ttmi.x = 0.0;
      ttmi.y = 0.0;
      ttmi.z = 0.0;
      ttmk[threadIdx.x].x = 0.0;
      ttmk[threadIdx.x].y = 0.0;
      ttmk[threadIdx.x].z = 0.0;

      same_block = ( idx!=kdx )? 0:1 ;

      //if(use_group){
      //  iga=grplist[iglob];
      //  igb=grplist[kglob];
      //}

      for ( j=0; j<WARP_SIZE; j++ ){
         srclane  = (ilane+j) & (WARP_SIZE-1);
         klane    = threadIdx.x-ilane + srclane;
         kglob_   =      __shfl_sync(ALL_LANES,kglob ,srclane);
         /*kp_.c    =      __shfl_sync(ALL_LANES,kp.c  ,srclane);
         kp_.dx   =      __shfl_sync(ALL_LANES,kp.dx ,srclane);
         kp_.dy   =      __shfl_sync(ALL_LANES,kp.dy ,srclane);
         kp_.dz   =      __shfl_sync(ALL_LANES,kp.dz ,srclane);
         kp_.qxx  =      __shfl_sync(ALL_LANES,kp.qxx,srclane);
         kp_.qxy  =      __shfl_sync(ALL_LANES,kp.qxy,srclane);
         kp_.qxz  =      __shfl_sync(ALL_LANES,kp.qxz,srclane);
         kp_.qyy  =      __shfl_sync(ALL_LANES,kp.qyy,srclane);
         kp_.qyz  =      __shfl_sync(ALL_LANES,kp.qyz,srclane);
         kp_.qzz  =      __shfl_sync(ALL_LANES,kp.qzz,srclane);*/
         //if (use_group) {
         // igb_ = __shfl_sync(ALL_LANES,igb,klane);
         // scale_ =  wgrp[iga*(ngrp+1)+igb_];
         //} else {
          scale_ = 1.0;
         //}

         if (nproc>1) {
            xk_   = posk[klane].x;
            yk_   = posk[klane].y;
            zk_   = posk[klane].z;
            pos.x = posi.x - xk_;
            pos.y = posi.y - yk_;
            pos.z = posi.z - zk_;
            accept_mid = Midpointimage(xk_,yk_,zk_,pos.x,pos.y,pos.z);
            pos.x=-pos.x; pos.y=-pos.y; pos.z=-pos.z; 
         }
         else {
            pos.x = posk[klane].x - posi.x;
            pos.y = posk[klane].y - posi.y;
            pos.z = posk[klane].z - posi.z;
            Image(pos.x,pos.y,pos.z);
         }
         d2      = pos.x*pos.x + pos.y*pos.y + pos.z*pos.z;
         do_pair = (same_block)? (iglob < kglob_):1;

         if (do_pair && d2<=off2 && accept_mid) {
            /* Compute one interaction
               This interaction is symetrical we thus don't need to switch comput */
            mpole1_couple<0>(d2,pos.x,pos.y,pos.z,ip,kp[klane],1.0-scale_,aewald,f,alsq2n,alsq2,
                             em_,frc,frc_i,frc_k[klane],ttmi,ttmk[klane]);
            /* increment local virial */
            vir_[0] -= pos.x * frc.x;
            vir_[1] -= pos.y * frc.x;
            vir_[2] -= pos.z * frc.x;
            vir_[3] -= pos.y * frc.y;
            vir_[4] -= pos.z * frc.y;
            vir_[5] -= pos.z * frc.z;
         }
      }

      /* Update energy */

      atomicAdd ( &em_buffer[ithread & (RED_BUFF_SIZE-1)], em_ );

      /* increment force-based gradient on both sites */
 
      atomicAdd( &dem[i   ][0],frc_i.x  );
      atomicAdd( &dem[i   ][1],frc_i.y  );
      atomicAdd( &dem[i   ][2],frc_i.z  );
      atomicAdd( &dem[kbis][0],frc_k[threadIdx.x].x );
      atomicAdd( &dem[kbis][1],frc_k[threadIdx.x].y );
      atomicAdd( &dem[kbis][2],frc_k[threadIdx.x].z );

      /* increment force-based torques on both sites */

      atomicAdd( &tem[i   ][0],ttmi.x );
      atomicAdd( &tem[i   ][1],ttmi.y );
      atomicAdd( &tem[i   ][2],ttmi.z );
      atomicAdd( &tem[kbis][0],ttmk[threadIdx.x].x );
      atomicAdd( &tem[kbis][1],ttmk[threadIdx.x].y );
      atomicAdd( &tem[kbis][2],ttmk[threadIdx.x].z );

      /* increment the virial due to pairwise Cartesian forces */

      atomicAdd( &vir_buffer[                  (ithread & (RED_BUFF_SIZE-1))],vir_[0] );
      atomicAdd( &vir_buffer[  RED_BUFF_SIZE + (ithread & (RED_BUFF_SIZE-1))],vir_[1] );
      atomicAdd( &vir_buffer[2*RED_BUFF_SIZE + (ithread & (RED_BUFF_SIZE-1))],vir_[2] );
      atomicAdd( &vir_buffer[3*RED_BUFF_SIZE + (ithread & (RED_BUFF_SIZE-1))],vir_[3] );
      atomicAdd( &vir_buffer[4*RED_BUFF_SIZE + (ithread & (RED_BUFF_SIZE-1))],vir_[4] );
      atomicAdd( &vir_buffer[5*RED_BUFF_SIZE + (ithread & (RED_BUFF_SIZE-1))],vir_[5] );
   }
}



int first_call_emreal = 1;
int gS_emreal         = 120;


EXTERN_C_BEG
void cu_emreal1c(MPOLE1_PARAMS1,cudaStream_t st,GROUPS_PARAMS){
   const int sh = 0;
   cudaError_t ierrSync;

   if (first_call_emreal){
      first_call_emreal=0;
      cudaKernelMaxGridSize(gS_emreal,cu_emreal1c_core,BLOCK_DIM,0)  /* This a Macro Function */
      if (tinkerdebug) printf (" gridSize emreal1c   %d \n", gS_emreal);
   }
   gS_emreal = npolelocnlb_pair/4;

   cu_emreal1c_core<<<gS_emreal,BLOCK_DIM,sh,st>>> (MPOLE1_ARGS,GROUPS_ARGS);
   if (tinkerdebug) ierrSync = cudaDeviceSynchronize();
   else             ierrSync = cudaGetLastError();
   if (ierrSync != cudaSuccess)
      printf("emreal1c_core C kernel error: %d \n  %s",ierrSync, cudaGetErrorString(ierrSync));
   return;
}
EXTERN_C_END
