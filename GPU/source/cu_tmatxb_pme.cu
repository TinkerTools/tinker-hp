#include "utils.h"
#include "image.h"

#define TMATXB_PARAMS_GEN                                              \
        const int*restrict ipole, const int*restrict pglob, const int*restrict ploc, const int*restrict ieblst, const int*restrict eblst   \
        , const real*restrict x, const real*restrict y, const real*restrict z, const real*restrict pdamp, const real*restrict thole, const real*restrict polarity   \
        , const real (*restrict mu)[6] , real (*restrict efi)[6]                                               \
        , const int npolelocnlb, const int npolelocnlb_pair, const int npolebloc, const int n, const int nproc, const int balanced \
        , const real cut2, const real alsq2, const real alsq2n, const real aewald
#define TMATXB_PARAMS                                               \
        TMATXB_PARAMS_GEN , MIDPOINTIMAGE_PARAMS
#define TMATXB_PARAMS1                                              \
        TMATXB_PARAMS_GEN , MIDPOINTIMAGE1_PARAMS

#define EFLD0_PARAMS_GEN                                              \
        const int*restrict ipole, const int*restrict pglob, const int*restrict ploc, const int*restrict ieblst, const int*restrict eblst   \
        , const real*restrict x, const real*restrict y, const real*restrict z, const real*restrict pdamp, const real*restrict thole, const real*restrict polarity, const real (*restrict rpole)[13]  \
        , real (*restrict efi)[6]                                               \
        , const int npolelocnlb, const int npolelocnlb_pair, const int npolebloc, const int n, const int nproc, const int balanced \
        , const real cut2, const real alsq2, const real alsq2n, const real aewald
#define OTFDC_EFLD0_PARAMS_GEN                                              \
        const int*restrict ipole, const int*restrict pglob, const int*restrict ploc, const int*restrict grplst, const int*restrict atmofst, const int*restrict npergrp, const int*restrict kofst, const int*restrict ieblst, const int*restrict eblst   \
        , const real*restrict x, const real*restrict y, const real*restrict z, const real*restrict pdamp, const real*restrict thole, const real*restrict polarity, const real (*restrict rpole)[13]  \
        , real (*restrict efi)[6], real* restrict zmat                            \
        , const int npolelocnlb, const int npolelocnlb_pair, const int npolebloc, const int n, const int nproc \
        , const real cut2, const real alsq2, const real alsq2n, const real aewald

#define EFLD0_PARAMS                                               \
        EFLD0_PARAMS_GEN , MIDPOINTIMAGE_PARAMS
#define OTFDC_EFLD0_PARAMS                                         \
        OTFDC_EFLD0_PARAMS_GEN , MIDPOINTIMAGE_PARAMS
#define EFLD0_PARAMS1                                              \
        EFLD0_PARAMS_GEN , MIDPOINTIMAGE1_PARAMS
#define OTFDC_EFLD0_PARAMS1                                        \
        OTFDC_EFLD0_PARAMS_GEN , MIDPOINTIMAGE1_PARAMS

#define TMATXB_ARGS                                                 \
          ipole,pglob,ploc,ieblst,eblst,x,y,z,pdamp,thole,polarity  \
        , mu,efi                                                    \
        , npolelocnlb,npolelocnlb_pair,npolebloc,n,nproc,balanced   \
        , cut2,alsq2,alsq2n,aewald                                  \
        , MIDPOINTIMAGE_ARGS

#define EFLD0_ARGS                                                       \
          ipole,pglob,ploc,ieblst,eblst,x,y,z,pdamp,thole,polarity,rpole \
        , efi                                                            \
        , npolelocnlb,npolelocnlb_pair,npolebloc,n,nproc,balanced        \
        , cut2,alsq2,alsq2n,aewald                                       \
        , MIDPOINTIMAGE_ARGS
#define OTFDC_EFLD0_ARGS                                                 \
          ipole,pglob,ploc,grplst,atmofst,npergrp,kofst,ieblst,eblst     \
        , x,y,z,pdamp,thole,polarity,rpole,efi,zmat                      \
        , npolelocnlb,npolelocnlb_pair,npolebloc,n,nproc                 \
        , cut2,alsq2,alsq2n,aewald                                       \
        , MIDPOINTIMAGE_ARGS

__device__ inline
void tmatxb_couple
        ( real& d2, const real3& dist, const real6& dpui, const real6& dpuk, const real& sdamp, const real& pgamma, const real& aewald, const real& alsq2, const real& alsq2n, const real& uscale,
                       real3& fid, real3& fip, real3& fkd, real3& fkp ){
   real d1,ralpha,exp2a;
   real bn0,bn1,bn2;
   real sdamp1,expdamp1;
   real rr3,rr5,rr3_bn1,rr5_bn2,sc3,sc5;
   real duir,dukr,puir,pukr;  /* Scalar products duir = (du(i).r) */

   /* compute the distances and the scaling factors according to Thole's model. */
   d1      = f_sqrt(d2);
   d2      = 1/d2;

   ralpha  = aewald * d1;
   d1      = 1/d1;
   exp2a   = f_exp(-ralpha*ralpha);
   bn0     = f_erfc(ralpha);  // becomes a macro with USE_ERFC_HASTINGS
   bn0     = bn0 * d1;
   bn1     = (     bn0 +         alsq2 * alsq2n * exp2a) * d2;
   bn2     = ( 3 * bn1 + alsq2 * alsq2 * alsq2n * exp2a) * d2;

   if( sdamp == 0.0 ) {
     sdamp1  = -100.0;
     sc3      =   1 - f_exp(sdamp1) * uscale;
     sc5      =   1 - f_exp(sdamp1) * uscale * (1 - sdamp1);
   }
   else {
     sdamp1 = - pgamma / ((d1*sdamp)*(d1*sdamp)*(d1*sdamp));
     if (sdamp1 > -50.0) {
       expdamp1 = f_exp(sdamp1);
       sc3      =   1 - expdamp1 * uscale;
       sc5      =   1 - expdamp1 * uscale * (1 - sdamp1);
     }
     else {
       sc3     = 1;
       sc5     = 1;
     }
   }

   /* compute the field */
   rr3     =     (1 - sc3) * (d1 * d2);
   rr5     = 3 * (1 - sc5) * (d1 * d2 * d2);
   rr3_bn1 = rr3 - bn1;
   rr5_bn2 = rr5 - bn2;

   duir    = dpui.x *dist.x + dpui.y *dist.y + dpui.z *dist.z;
   dukr    = dpuk.x *dist.x + dpuk.y *dist.y + dpuk.z *dist.z;

   puir    = dpui.xx *dist.x + dpui.yy *dist.y + dpui.zz *dist.z;
   pukr    = dpuk.xx *dist.x + dpuk.yy *dist.y + dpuk.zz *dist.z;

   fid.x  += -rr3_bn1*dpuk.x  + rr5_bn2*dukr*dist.x;
   fid.y  += -rr3_bn1*dpuk.y  + rr5_bn2*dukr*dist.y;
   fid.z  += -rr3_bn1*dpuk.z  + rr5_bn2*dukr*dist.z;

   fip.x  += -rr3_bn1*dpuk.xx + rr5_bn2*pukr*dist.x;
   fip.y  += -rr3_bn1*dpuk.yy + rr5_bn2*pukr*dist.y;
   fip.z  += -rr3_bn1*dpuk.zz + rr5_bn2*pukr*dist.z;

   fkd.x  += -rr3_bn1*dpui.x  + rr5_bn2*duir*dist.x;
   fkd.y  += -rr3_bn1*dpui.y  + rr5_bn2*duir*dist.y;
   fkd.z  += -rr3_bn1*dpui.z  + rr5_bn2*duir*dist.z;

   fkp.x  += -rr3_bn1*dpui.xx + rr5_bn2*puir*dist.x;
   fkp.y  += -rr3_bn1*dpui.yy + rr5_bn2*puir*dist.y;
   fkp.z  += -rr3_bn1*dpui.zz + rr5_bn2*puir*dist.z;
}

__device__ inline
void efld0_couple(const real d2,const real3& pos, const rpole_elt& ip, const rpole_elt& kp, const real& alsq2, const real& alsq2n,
         const real& aewald, const real& damp, const real& pgamma, const real& dscale, const real& pscale,
        real3& fid,real3& fip,real3& fkd,real3& fkp,real& d1, real& sc3, real& sc5, real& bn1, real& bn2, const int do_correct) {

   real exp2a;
   real invd1,invd2,invd3,invd5,invd7;
   real sc7;
   real drr3,drr5,drr7,prr3,prr5,prr7;
   real dir,qirr,dkr,qkrr;
   real qirx,qiry,qirz,qkrx,qkry,qkrz;
   real fkmx,fkmy,fkmz,fimx,fimy,fimz;
   real invdamp,expdamp1,damp1;
   real ralpha,bn0,bn3;
   const real one=1.;
   const real two=2.;

   damp1   = -100.0;
   invdamp = 1/damp;
   invd2   = 1/d2;
   d1      = f_sqrt(d2);
   invd1   = 1/d1;

   sc3     = one;
   sc5     = one;
   sc7     = one;

   invd3   = invd1  * invd2;
   invd5   = invd3  * invd2;
   invd7   = invd5  * invd2;

   if (damp!=0.0) damp1 = - pgamma*(d1*invdamp)*(d1*invdamp)*(d1*invdamp);

   if (damp1 > -50.0)  {
      expdamp1  = f_exp(damp1);
      sc3  = one - expdamp1;
      sc5  = one - expdamp1*(one - damp1);
      sc7  = one - expdamp1*(one - damp1 + 0.6*damp1*damp1);
   }

   if (do_correct) {
      /*  [dp]scale equal to 1-[dp]scale in this case */
      drr3    =      sc3*dscale * invd3;
      drr5    =  3 * sc5*dscale * invd5;
      drr7    = 15 * sc7*dscale * invd7;

      prr3    =      sc3*pscale * invd3;
      prr5    =  3 * sc5*pscale * invd5;
      prr7    = 15 * sc7*pscale * invd7;
   }
   else {

      /* calculate the error function damping terms */
      ralpha  = aewald * d1;
      exp2a   = f_exp( -ralpha*ralpha );
      bn0     = f_erfc(ralpha);  // becomes a macro with USE_ERFC_HASTINGS
 
      bn0     =    bn0                                     * invd1;
      bn1     = (  bn0  + alsq2             *alsq2n*exp2a) * invd2;
      bn2     = (3*bn1  + alsq2*alsq2       *alsq2n*exp2a) * invd2;
      bn3     = (5*bn2  + alsq2*alsq2*alsq2 *alsq2n*exp2a) * invd2;

      drr3    =      (one - sc3*dscale) * invd3;
      drr5    =  3 * (one - sc5*dscale) * invd5;
      drr7    = 15 * (one - sc7*dscale) * invd7;
      
      prr3    =      (one - sc3*pscale) * invd3;
      prr5    =  3 * (one - sc5*pscale) * invd5;
      prr7    = 15 * (one - sc7*pscale) * invd7;
   }
 
   /* compute some intermediate quantities */
   dir     = ip.dx *pos.x + ip.dy *pos.y + ip.dz *pos.z;
   qirx    = ip.qxx*pos.x + ip.qxy*pos.y + ip.qxz*pos.z;
   qiry    = ip.qxy*pos.x + ip.qyy*pos.y + ip.qyz*pos.z;
   qirz    = ip.qxz*pos.x + ip.qyz*pos.y + ip.qzz*pos.z;
   qirr    =   qirx*pos.x +   qiry*pos.y +   qirz*pos.z;

   dkr     = kp.dx *pos.x + kp.dy *pos.y +  kp.dz *pos.z;
   qkrx    = kp.qxx*pos.x + kp.qxy*pos.y +  kp.qxz*pos.z;
   qkry    = kp.qxy*pos.x + kp.qyy*pos.y +  kp.qyz*pos.z;
   qkrz    = kp.qxz*pos.x + kp.qyz*pos.y +  kp.qzz*pos.z;
   qkrr    =   qkrx*pos.x +   qkry*pos.y +    qkrz*pos.z;

   if (do_correct) {
      fimx = 0.0; fimy = 0.0; fimz = 0.0;
      fkmx = 0.0; fkmy = 0.0; fkmz = 0.0;
   }
   else {
      fimx = -( bn1*kp.c  - bn2*dkr + bn3*qkrr )*pos.x - bn1*kp.dx + two*bn2*qkrx;
      fimy = -( bn1*kp.c  - bn2*dkr + bn3*qkrr )*pos.y - bn1*kp.dy + two*bn2*qkry;
      fimz = -( bn1*kp.c  - bn2*dkr + bn3*qkrr )*pos.z - bn1*kp.dz + two*bn2*qkrz;
      fkmx =  ( bn1*ip.c  + bn2*dir + bn3*qirr )*pos.x - bn1*ip.dx - two*bn2*qirx;
      fkmy =  ( bn1*ip.c  + bn2*dir + bn3*qirr )*pos.y - bn1*ip.dy - two*bn2*qiry;
      fkmz =  ( bn1*ip.c  + bn2*dir + bn3*qirr )*pos.z - bn1*ip.dz - two*bn2*qirz;
   }

   fid.x  += fimx + ( drr3*kp.c  - drr5*dkr + drr7*qkrr )*pos.x +  drr3*kp.dx - two*drr5*qkrx;
   fid.y  += fimy + ( drr3*kp.c  - drr5*dkr + drr7*qkrr )*pos.y +  drr3*kp.dy - two*drr5*qkry;
   fid.z  += fimz + ( drr3*kp.c  - drr5*dkr + drr7*qkrr )*pos.z +  drr3*kp.dz - two*drr5*qkrz;
   fip.x  += fimx + ( prr3*kp.c  - prr5*dkr + prr7*qkrr )*pos.x +  prr3*kp.dx - two*prr5*qkrx;
   fip.y  += fimy + ( prr3*kp.c  - prr5*dkr + prr7*qkrr )*pos.y +  prr3*kp.dy - two*prr5*qkry;
   fip.z  += fimz + ( prr3*kp.c  - prr5*dkr + prr7*qkrr )*pos.z +  prr3*kp.dz - two*prr5*qkrz;

   fkd.x  += fkmx - ( drr3*ip.c  + drr5*dir + drr7*qirr )*pos.x +  drr3*ip.dx + two*drr5*qirx;
   fkd.y  += fkmy - ( drr3*ip.c  + drr5*dir + drr7*qirr )*pos.y +  drr3*ip.dy + two*drr5*qiry;
   fkd.z  += fkmz - ( drr3*ip.c  + drr5*dir + drr7*qirr )*pos.z +  drr3*ip.dz + two*drr5*qirz;
   fkp.x  += fkmx - ( prr3*ip.c  + prr5*dir + prr7*qirr )*pos.x +  prr3*ip.dx + two*prr5*qirx;
   fkp.y  += fkmy - ( prr3*ip.c  + prr5*dir + prr7*qirr )*pos.y +  prr3*ip.dy + two*prr5*qiry;
   fkp.z  += fkmz - ( prr3*ip.c  + prr5*dir + prr7*qirr )*pos.z +  prr3*ip.dz + two*prr5*qirz;

}

__global__ void check_loc( const int*restrict pglob, const int*restrict ipole, const int*restrict ploc, const int npolelocnlb, const int nbloc, const int rank ){
   for ( int ii=threadIdx.x + blockIdx.x*blockDim.x;ii<npolelocnlb;ii+=blockDim.x*gridDim.x ){
       const int iipole = pglob[ii];
       const int iploc  = ploc [ii];
       if (iploc==0 || iploc>nbloc ) printf("check_ploc(%d) pole(%d) glob(%d) Idx(%d) rank(%d)\n",iploc,ipole[ii],iipole,ii+1,rank);
       //iploc = ploc_s(iipole)
       //if (iploc==0 || iploc.gt.nbloc) printf("out ploc %d %d %d %d",iploc,iipole,i,rank);
   }
}

__global__
void cu_efld0_direct_core (EFLD0_PARAMS1){

   const int ithread = threadIdx.x + blockIdx.x*blockDim.x;
   const int iwarp   =               ithread / WARP_SIZE;
   const int nwarp   =  blockDim.x*gridDim.x / WARP_SIZE;
   const int ilane   = threadIdx.x & (WARP_SIZE-1);
   int accept_mid    = 1;

   int klane,srclane;
   int ii,j;
   int iblock,idx,kdx;
   int iipole,iglob,iploc,kpole,kglob,kploc,kglob_;
   int do_pair,same_block;
   real xk_,yk_,zk_,d2;
   real ipdp,ipgm,kpdp,kpgm,pdp,pgm;
   rpole_elt ip;
   real3 posi,pos;
   real3 fid,fip;
   __shared__ real3 posk[BLOCK_DIM],fkd[BLOCK_DIM],fkp[BLOCK_DIM];
   __shared__ rpole_elt kp[BLOCK_DIM];
   //__shared__ int ncalc[4];
   //__shared__ int cont;

   //if (ithread==0) printf( " %i %i %i %i %i %i %i " r_Format r_Format r_Format "\n", nwarp,nproc,npolelocnlb,npolebloc,n,npolelocnlb_pair,cut2,alsq2,aewald);

   for ( ii=iwarp; ii<npolelocnlb_pair; ii+=nwarp ){
      /*  Load atom block i parameters */
      iblock  = ieblst[ii];
      if (iblock==0) continue;
      idx     = (iblock-1)*WARP_SIZE + ilane;
      iipole  = pglob[idx] -1;
      iglob   = ipole[idx] -1;
      iploc   = ploc [idx] -1;
      posi.x  = x[idx];
      posi.y  = y[idx];
      posi.z  = z[idx];
      ipdp    = pdamp[iipole];
      ipgm    = thole[iipole];
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
      kdx     = eblst[ii*WARP_SIZE+ ilane] -1;
      kpole   = pglob[kdx] -1;
      kglob   = ipole[kdx] -1;
      kploc   = ploc [kdx] -1;
      posk[threadIdx.x].x  = x[kdx];
      posk[threadIdx.x].y  = y[kdx];
      posk[threadIdx.x].z  = z[kdx];
      kpdp    = pdamp[kpole];
      kpgm    = thole[kpole];
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
      //if (ilane==1) ncalc[threadIdx.x/WARP_SIZE]=0;

      /* set compute Data to]0 */
      fid.x   = 0;
      fid.y   = 0;
      fid.z   = 0;
      fip.x   = 0;
      fip.y   = 0;
      fip.z   = 0;
      fkd[threadIdx.x].x = 0;
      fkd[threadIdx.x].y = 0;
      fkd[threadIdx.x].z = 0;
      fkp[threadIdx.x].x = 0;
      fkp[threadIdx.x].y = 0;
      fkp[threadIdx.x].z = 0;
      //cont=0;

      same_block = ( idx!=kdx )? 0:1 ;

      for ( j=0; j<WARP_SIZE; j++ ){
         srclane  = (ilane+j) & (WARP_SIZE-1);
         klane    = threadIdx.x-ilane + srclane;
         kglob_   =      __shfl_sync(ALL_LANES,kglob ,srclane);
         pdp      = ipdp*__shfl_sync(ALL_LANES,kpdp  ,srclane);
         pgm      =      __shfl_sync(ALL_LANES,kpgm  ,srclane);
         if (ipgm<pgm) pgm = ipgm;

         if (nproc>1 && balanced ) {
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

         if (do_pair && d2<=cut2 && accept_mid) {
            //atomicAdd( &cont, 1);
            /* Compute one interaction
               Since the interaction is not symetrical we need to switch comput when necessary */
            real d,sc3,sc5,bn1,bn2;
            //if (iblock<500) atomicAdd( &ncalc[threadIdx.x/WARP_SIZE],1 );
            if (iglob<kglob_)
               efld0_couple(d2,pos,ip,kp[klane],alsq2,alsq2n,aewald,pdp,pgm,1.0,1.0,
                            fid,fip,fkd[klane],fkp[klane],d,sc3,sc5,bn1,bn2,0);
            else{
               pos.x = -pos.x; pos.y=-pos.y; pos.z=-pos.z;
               efld0_couple(d2,pos,kp[klane],ip,alsq2,alsq2n,aewald,pdp,pgm,1.0,1.0,
                            fkd[klane],fkp[klane],fid,fip,d,sc3,sc5,bn1,bn2,0);
            }
         }
      }

      //if (ilane==0) printf("%d  %d  %d \n",ii,iblock,cont);
      /* increment electric field for each atoms */
      atomicAdd( &efi[iploc][0],fid.x );
      atomicAdd( &efi[iploc][1],fid.y );
      atomicAdd( &efi[iploc][2],fid.z );
      atomicAdd( &efi[iploc][3],fip.x );
      atomicAdd( &efi[iploc][4],fip.y );
      atomicAdd( &efi[iploc][5],fip.z );
      atomicAdd( &efi[kploc][0],fkd[threadIdx.x].x );
      atomicAdd( &efi[kploc][1],fkd[threadIdx.x].y );
      atomicAdd( &efi[kploc][2],fkd[threadIdx.x].z );
      atomicAdd( &efi[kploc][3],fkp[threadIdx.x].x );
      atomicAdd( &efi[kploc][4],fkp[threadIdx.x].y );
      atomicAdd( &efi[kploc][5],fkp[threadIdx.x].z );
      __syncwarp(ALL_LANES);
      //if (ilane==1 && iblock<500) printf("ii %8d %8d    %8d \n",ii+1,iblock,ncalc[threadIdx.x/WARP_SIZE]);
   }
}

__global__
void cu_otfdc_efld0_direct_core (OTFDC_EFLD0_PARAMS1){

   const int ithread = threadIdx.x + blockIdx.x*blockDim.x;
   const int iwarp   =               ithread / WARP_SIZE;
   const int nwarp   =  blockDim.x*gridDim.x / WARP_SIZE;
   const int ilane   = threadIdx.x & (WARP_SIZE-1);
   int accept_mid    = 1;

   int klane,srclane;
   int ii,j;
   int iblock,idx,kdx,li,atii,maxrow,kofi;
   int iipole,iglob,iploc,kpole,kglob,kploc,kglob_;
   int do_pair,same_block;
   real xk_,yk_,zk_,d2;
   real ipdp,ipgm,kpdp,kpgm,pdp,pgm;
   rpole_elt ip;
   real3 posi,pos;
   real3 fid,fip;
   __shared__ int lk[BLOCK_DIM];
   __shared__ real3 posk[BLOCK_DIM],fkd[BLOCK_DIM],fkp[BLOCK_DIM];
   __shared__ rpole_elt kp[BLOCK_DIM];

   //if (ithread==0) printf( " %i %i %i %i %i %i %i " r_Format r_Format r_Format "\n", nwarp,nproc,npolelocnlb,npolebloc,n,npolelocnlb_pair,cut2,alsq2,aewald);

   for ( ii=iwarp; ii<npolelocnlb_pair; ii+=nwarp ){
      /*  Load atom block i parameters */
      iblock  = ieblst[ii];
      if (iblock==0) continue;
      idx     = (iblock-1)*WARP_SIZE + ilane;
      iipole  = pglob[idx] -1;
      iglob   = ipole[idx] -1;
      iploc   = ploc [idx] -1;
      posi.x  = x[idx];
      posi.y  = y[idx];
      posi.z  = z[idx];
      // gather zmat data 
      li      =   grplst[iglob];
      atii    = (atmofst[iglob]-1)*3;
      maxrow  =  npergrp[li-1]*3;
      kofi    =    kofst[li-1];

      ipdp    = pdamp[iipole];
      ipgm    = thole[iipole];
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
      kdx     = eblst[ii*WARP_SIZE+ ilane] -1;
      kpole   = pglob[kdx] -1;
      kglob   = ipole[kdx] -1;
      kploc   = ploc [kdx] -1;
      lk  [threadIdx.x]    = grplst[kglob];
      posk[threadIdx.x].x  = x[kdx];
      posk[threadIdx.x].y  = y[kdx];
      posk[threadIdx.x].z  = z[kdx];
      kpdp    = pdamp[kpole];
      kpgm    = thole[kpole];
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
      fid.x   = 0;
      fid.y   = 0;
      fid.z   = 0;
      fip.x   = 0;
      fip.y   = 0;
      fip.z   = 0;
      fkd[threadIdx.x].x = 0;
      fkd[threadIdx.x].y = 0;
      fkd[threadIdx.x].z = 0;
      fkp[threadIdx.x].x = 0;
      fkp[threadIdx.x].y = 0;
      fkp[threadIdx.x].z = 0;

      same_block = ( idx!=kdx )? 0:1 ;

      for ( j=0; j<WARP_SIZE; j++ ){
         srclane  = (ilane+j) & (WARP_SIZE-1);
         klane    = threadIdx.x-ilane + srclane;
         kglob_   =      __shfl_sync(ALL_LANES,kglob ,srclane);
         pdp      = ipdp*__shfl_sync(ALL_LANES,kpdp  ,srclane);
         pgm      =      __shfl_sync(ALL_LANES,kpgm  ,srclane);
         if (ipgm<pgm) pgm = ipgm;

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

         if (do_pair && d2<=cut2 && accept_mid) {
            /* Compute one interaction
               Since the interaction is not symetrical we need to switch comput when necessary */
            real d,sc3,sc5,bn1,bn2;
            if (iglob<kglob_)
               efld0_couple(d2,pos,ip,kp[klane],alsq2,alsq2n,aewald,pdp,pgm,1.0,1.0,
                            fid,fip,fkd[klane],fkp[klane],d,sc3,sc5,bn1,bn2,0);
            else{
               pos.x = -pos.x; pos.y=-pos.y; pos.z=-pos.z;
               efld0_couple(d2,pos,kp[klane],ip,alsq2,alsq2n,aewald,pdp,pgm,1.0,1.0,
                            fkd[klane],fkp[klane],fid,fip,d,sc3,sc5,bn1,bn2,0);
            }
            if (li==lk[klane] && li!=-1) {

               bn1    = bn1 -     (1.0 - sc3)/ (d*d2);
               bn2    = bn2 - 3.0*(1.0 - sc5)/ (d*d2*d2);
               int atkk = (atmofst[kglob_] - 1)*3;
               int cofst1,cofst2,cofst3,rofst;
               if (atii < atkk) {
                 cofst1 = atii + 1;
                 cofst2 = atii + 2;
                 cofst3 = atii + 3;
                 rofst  = atkk;
               } else {
                 cofst1 = atkk + 1;
                 cofst2 = atkk + 2;
                 cofst3 = atkk + 3;
                 rofst  = atii;
               }

               cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2;
               cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2;
               cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2;
               zmat[rofst+0+cofst1+kofi] =  bn1 - bn2*pos.x*pos.x;
               zmat[rofst+1+cofst1+kofi] =      - bn2*pos.x*pos.y;
               zmat[rofst+2+cofst1+kofi] =      - bn2*pos.x*pos.z;
               zmat[rofst+0+cofst2+kofi] =      - bn2*pos.x*pos.y;
               zmat[rofst+1+cofst2+kofi] =  bn1 - bn2*pos.y*pos.y;
               zmat[rofst+2+cofst2+kofi] =      - bn2*pos.y*pos.z;
               zmat[rofst+0+cofst3+kofi] =      - bn2*pos.x*pos.z;
               zmat[rofst+1+cofst3+kofi] =      - bn2*pos.y*pos.z;
               zmat[rofst+2+cofst3+kofi] =  bn1 - bn2*pos.z*pos.z;

            }
         }
      }

      /* increment electric field for each atoms */
      atomicAdd( &efi[iploc][0],fid.x );
      atomicAdd( &efi[iploc][1],fid.y );
      atomicAdd( &efi[iploc][2],fid.z );
      atomicAdd( &efi[iploc][3],fip.x );
      atomicAdd( &efi[iploc][4],fip.y );
      atomicAdd( &efi[iploc][5],fip.z );
      atomicAdd( &efi[kploc][0],fkd[threadIdx.x].x );
      atomicAdd( &efi[kploc][1],fkd[threadIdx.x].y );
      atomicAdd( &efi[kploc][2],fkd[threadIdx.x].z );
      atomicAdd( &efi[kploc][3],fkp[threadIdx.x].x );
      atomicAdd( &efi[kploc][4],fkp[threadIdx.x].y );
      atomicAdd( &efi[kploc][5],fkp[threadIdx.x].z );
   }
}


__global__
void cu_tmatxb_pme_core (TMATXB_PARAMS1){

   const int ithread = threadIdx.x + blockIdx.x*blockDim.x;
   const int iwarp   =               ithread / WARP_SIZE;
   const int nwarp   =  blockDim.x*gridDim.x / WARP_SIZE;
   const int ilane   = threadIdx.x & (WARP_SIZE-1);
   int accept_mid    = 1;

   int klane,srclane;
   int ii,j;
   int iblock,idx,kdx;
   int iipole,iglob,iploc,kpole,kploc;
   int do_pair,same_block;
   real xk_,yk_,zk_,d2;
   real ipdp,ipgm,pdp,pgm;
   real6 dpuk_;
   real3 pos;
   __shared__ int kglob[BLOCK_DIM], ikstat[BLOCK_DIM];
   __shared__ real kpdp[BLOCK_DIM],kpgm[BLOCK_DIM];
   __shared__ real3 posk[BLOCK_DIM],posi[BLOCK_DIM];
   __shared__ real3 fkd[BLOCK_DIM],fkp[BLOCK_DIM];
   __shared__ real3 fid[BLOCK_DIM],fip[BLOCK_DIM];
   __shared__ real6 dpui[BLOCK_DIM],dpuk[BLOCK_DIM];


   //if (ithread==0) printf( " %i %i %i %i %i %i %i " r_Format r_Format r_Format "\n", nwarp,nproc,npolelocnlb,npolebloc,n,npolelocnlb_pair,cut2,alsq2,aewald);

   for ( ii=iwarp; ii<npolelocnlb_pair; ii+=nwarp ){
      /*  Load atom block i parameters */
      iblock  = ieblst[ii];
      if (iblock==0) continue;
      idx     = (iblock-1)*WARP_SIZE + ilane;
      iipole  = pglob[idx] -1;
      iglob   = ipole[idx] -1;
      iploc   = ploc [idx] -1;
      posi[threadIdx.x].x  = x[idx];
      posi[threadIdx.x].y  = y[idx];
      posi[threadIdx.x].z  = z[idx];
      ipdp    = pdamp[iipole];
      ipgm    = thole[iipole];
      dpui[threadIdx.x].x  = mu[iploc][0];
      dpui[threadIdx.x].y  = mu[iploc][1];
      dpui[threadIdx.x].z  = mu[iploc][2];
      dpui[threadIdx.x].xx = mu[iploc][3];
      dpui[threadIdx.x].yy = mu[iploc][4];
      dpui[threadIdx.x].zz = mu[iploc][5];

      /*  Load atom block k parameters */
      kdx     = eblst[ii*WARP_SIZE+ ilane] -1;
      kpole   = pglob[kdx] -1;
      kglob[threadIdx.x]   = ipole[kdx] -1;
      kploc   = ploc [kdx] -1;
      posk[threadIdx.x].x  = x[kdx];
      posk[threadIdx.x].y  = y[kdx];
      posk[threadIdx.x].z  = z[kdx];
      kpdp[threadIdx.x]    = pdamp[kpole];
      kpgm[threadIdx.x]    = thole[kpole];
      dpuk[threadIdx.x].x  = mu[kploc][0];
      dpuk[threadIdx.x].y  = mu[kploc][1];
      dpuk[threadIdx.x].z  = mu[kploc][2];
      dpuk[threadIdx.x].xx = mu[kploc][3];
      dpuk[threadIdx.x].yy = mu[kploc][4];
      dpuk[threadIdx.x].zz = mu[kploc][5];

      /* set compute Data to]0 */
      ikstat[threadIdx.x]   = 0;
      fid   [threadIdx.x].x = 0;
      fid   [threadIdx.x].y = 0;
      fid   [threadIdx.x].z = 0;
      fip   [threadIdx.x].x = 0;
      fip   [threadIdx.x].y = 0;
      fip   [threadIdx.x].z = 0;
      fkd   [threadIdx.x].x = 0;
      fkd   [threadIdx.x].y = 0;
      fkd   [threadIdx.x].z = 0;
      fkp   [threadIdx.x].x = 0;
      fkp   [threadIdx.x].y = 0;
      fkp   [threadIdx.x].z = 0;

      same_block = ( idx!=kdx )? 0:1 ;

      #pragma unroll
      for ( int i=0; i<2; i++ ){
         srclane     = (ilane+i) & (WARP_SIZE-1);
         int iilane  = threadIdx.x-ilane + srclane;
      for ( j=0; j<WARP_SIZE; j++ ){
         if (atomicOr( &ikstat[iilane],1<<j ) & 1<<j) continue;

         srclane  = (ilane+j) & (WARP_SIZE-1);
         klane    = threadIdx.x-ilane + srclane;
         dpuk_.x  = dpuk[klane].x ;
         dpuk_.y  = dpuk[klane].y ;
         dpuk_.z  = dpuk[klane].z ;
         dpuk_.xx = dpuk[klane].xx;
         dpuk_.yy = dpuk[klane].yy;
         dpuk_.zz = dpuk[klane].zz;
         pdp      = ipdp*kpdp[klane];
         pgm      =      kpgm[klane];
         if (ipgm<pgm) pgm = ipgm;

         if (nproc>1) {
            xk_   = posk[klane].x;
            yk_   = posk[klane].y;
            zk_   = posk[klane].z;
            pos.x = posi[iilane].x - xk_;
            pos.y = posi[iilane].y - yk_;
            pos.z = posi[iilane].z - zk_;
            accept_mid = Midpointimage(xk_,yk_,zk_,pos.x,pos.y,pos.z);
         }
         else {
            pos.x = posi[iilane].x - posk[klane].x;
            pos.y = posi[iilane].y - posk[klane].y;
            pos.z = posi[iilane].z - posk[klane].z;
            Image(pos.x,pos.y,pos.z);
         }
         d2      = pos.x*pos.x + pos.y*pos.y + pos.z*pos.z;
         do_pair = (same_block)? iglob < kglob[klane] : 1 ;

         if (do_pair && d2<=cut2 && accept_mid) {
             /* compute one interaction */
             tmatxb_couple(d2,pos,dpui[iilane],dpuk_,pdp,pgm,aewald,alsq2,alsq2n,1.
                          ,fid[iilane],fip[iilane],fkd[klane],fkp[klane]);
         }
      }
      }
      __syncwarp(ALL_LANES);
      /*if (ii==0&&ilane==3) printf (" %5d %5d %5i %7d " r10_Format "\n" ,
         ii,iglob,klane,kdx,dpui.x,dpuk_.x,dpui.y,dpuk_.y,
         fid.x,fid.y,fip.x); */

      /* increment electric field for each atoms */
      atomicAdd( &efi[iploc][0],fid[threadIdx.x].x );
      atomicAdd( &efi[iploc][1],fid[threadIdx.x].y );
      atomicAdd( &efi[iploc][2],fid[threadIdx.x].z );
      atomicAdd( &efi[iploc][3],fip[threadIdx.x].x );
      atomicAdd( &efi[iploc][4],fip[threadIdx.x].y );
      atomicAdd( &efi[iploc][5],fip[threadIdx.x].z );
      atomicAdd( &efi[kploc][0],fkd[threadIdx.x].x );
      atomicAdd( &efi[kploc][1],fkd[threadIdx.x].y );
      atomicAdd( &efi[kploc][2],fkd[threadIdx.x].z );
      atomicAdd( &efi[kploc][3],fkp[threadIdx.x].x );
      atomicAdd( &efi[kploc][4],fkp[threadIdx.x].y );
      atomicAdd( &efi[kploc][5],fkp[threadIdx.x].z );
   }
}


EXTERN_C_BEG

int nproc      = 1;
int rank       = 0;
int devicenum  = 0;
int tinkerdebug= 0;
cudaDeviceProp devProp;

void C_init_env (int devicenum_=0, int nproc_=1, int rank_=0, int tinkerdebug_=0){
   cudaGetDeviceProperties(&devProp,devicenum_);
   devicenum   = devicenum_;
   nproc       = nproc_;
   rank        = rank_;
   tinkerdebug = tinkerdebug_;
}

real _xcell;
real _ycell;
real _zcell;
real _ixcell;
real _iycell;
real _izcell;
real eps_cell;
real _box34;
int octahedron;

void C_get_cell( real xcell_, real ycell_, real zcell_, real eps_cell_, int octa_, real box34_ ){
   _xcell   = xcell_;
   _ycell   = ycell_;
   _zcell   = zcell_;
   _ixcell  = (real) 1.0/ (double)xcell_;
   _iycell  = (real) 1.0/ (double)ycell_;
   _izcell  = (real) 1.0/ (double)zcell_;
   eps_cell = eps_cell_;
   octahedron = octa_;
   _box34   = box34_;
}

int dynamic_gS        = 1;
int first_call_efld0  = 1;
int gS_efld    = 160;
int gS_loc     = 160;
const int maxBlock = 1<<16;


void cu_efld0_direct(EFLD0_PARAMS,cudaStream_t st){
   const int sh = 0;
   cudaError_t ierrSync;

   if (first_call_efld0){
      first_call_efld0=0;
      cudaKernelMaxGridSize(gS_efld,cu_efld0_direct_core,BLOCK_DIM,0)  /* This a Macro Function */
      if(rank==0 && tinkerdebug&1) {
         printf (" gridSize efld0     %d \n", gS_efld);
         printf (" balanced computation     %d \n ", balanced);
      }
      if (nproc>1) {
         cudaKernelMaxGridSize(gS_loc,check_loc,BLOCK_DIM,0)  /* This a Macro Function */
         check_loc<<<gS_loc,BLOCK_DIM,sh,st>>>(pglob,ipole,ploc,npolelocnlb,npolebloc,rank);
         ierrSync = tinkerdebug ? cudaDeviceSynchronize() : cudaGetLastError();
         if (ierrSync != cudaSuccess) printf("check_loc kernel error: %d ( %s )\n",ierrSync,cudaGetErrorString(ierrSync));
      }
   }

   if (dynamic_gS) gS_efld= (npolelocnlb_pair>>2 < maxBlock) ? npolelocnlb_pair>>2 : maxBlock ;
   cu_efld0_direct_core<<<gS_efld,BLOCK_DIM,sh,st>>> (EFLD0_ARGS);
   ierrSync = tinkerdebug ? cudaDeviceSynchronize() : cudaGetLastError();
   if (ierrSync != cudaSuccess) printf("efld0_direct_core C kernel error: %d ( %s )\n",ierrSync,cudaGetErrorString(ierrSync));

   return;
}

int first_call_otfdc_efld0  = 1;
int gS_otfdc_efld           = 160;

void cu_otfdc_efld0_direct(OTFDC_EFLD0_PARAMS,cudaStream_t st){
   const int sh = 0;

   if (first_call_otfdc_efld0){
      first_call_otfdc_efld0=0;
      cudaKernelMaxGridSize(gS_otfdc_efld,cu_otfdc_efld0_direct_core,BLOCK_DIM,0)  /* This a Macro Function */
      if (rank==0 && tinkerdebug&1) {
         printf (" gridSize oftdc_efld0     %d \n", gS_efld);
      }
   }
   if (dynamic_gS) gS_otfdc_efld= npolelocnlb_pair/4;

   cu_otfdc_efld0_direct_core<<<gS_otfdc_efld,BLOCK_DIM,sh,st>>> (OTFDC_EFLD0_ARGS);
   if  (tinkerdebug) gpuErrchk( cudaDeviceSynchronize() )
   else              gpuErrchk( cudaGetLastError() )
   return;
}

int first_call_tmatxb = 1;
int gS_tmat           = 160;

void cu_tmatxb_pme(TMATXB_PARAMS,cudaStream_t st){
   //int gS = 160;
   const int sh = 0;

   if (first_call_tmatxb){
      first_call_tmatxb = 0;
      cudaKernelMaxGridSize(gS_tmat,cu_tmatxb_pme_core,BLOCK_DIM,0)  /* This a Macro Function */
      if (rank==0 && tinkerdebug&1) printf (" gridSize tmatxb_cu %d \n", gS_tmat);
   }
   if (dynamic_gS) gS_tmat= npolelocnlb_pair/8;

   cu_tmatxb_pme_core<<<gS_tmat,BLOCK_DIM,sh,st>>> (TMATXB_ARGS);
   cudaError_t ierrSync;
   if(tinkerdebug) ierrSync = cudaDeviceSynchronize();
   else            ierrSync = cudaGetLastError();
   if (ierrSync != cudaSuccess)
      printf("tmatxb_pme_core C kernel error: %d \n  %s",ierrSync, cudaGetErrorString(ierrSync));
   return;
}
EXTERN_C_END
