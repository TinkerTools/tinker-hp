c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1cu" : driver for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates on device
c
c
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_precision.h"
#include "tinker_cudart.h"
#include "tinker_types.h"
      module empole1cu
        use utilcu   ,only: nproc,BLOCK_DIM
        use utilgpu  ,only: BLOCK_SIZE,RED_BUFF_SIZE,WARP_SIZE
        use tinTypes ,only: real3,real6,mdyn3_r,rpole_elt 
        use tinheader,only: zeror,oner
        use sizes    ,only: maxgrp
        implicit none
        private
#include "atomicOp.h.f"

        public:: emreal1_kcu,emreal1shr_kcu,emreal1lgr_kcu
     &         , emreal3_kcu
     &         , emreal_scaling_cu,emrealShlg_scaling_cu
        contains

#include "convert.f.inc"
#include "image.f.inc"
#include "midpointimage.f.inc"
#include "groups.inc.f"
#include "atomicOp.inc.f"
#include "pair_mpole1.f.inc"

        M_subroutine mdyn3r_zr(r3)
        type(mdyn3_r),intent(out)::r3
        r3%x=0;r3%y=0;r3%z=0;
        end subroutine
        M_subroutine real3_zr(r3)
        type(real3),intent(out)::r3
        r3%x=0;r3%y=0;r3%z=0;
        end subroutine

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_lambdadyn__)
#define __sufx__ 1_kcu
#include "empolecu.tpl.f"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_lambdadyn__+__use_shortRange__)
#define __sufx__ 1shr_kcu
#include "empolecu.tpl.f"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_grd__+__use_ene__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_groups__+__use_lambdadyn__+__use_longRange__)
#define __sufx__ 1lgr_kcu
#include "empolecu.tpl.f"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

#define __tver__ (__use_ene__+__use_act__)
#define __tfea__ (__use_mpi__+__use_groups__)
#define __sufx__ 3_kcu
#include "empolecu.tpl.f"
#undef  __tver__
#undef  __tfea__
#undef  __sufx__

        attributes(global) subroutine emreal_scaling_cu
     &            (mcorrect_ik,mcorrect_scale,ipole,loc,locp,x,y,z,rpole
     &            ,dem,tem,em_buff,vir_buff
     &            ,n,nbloc,n_mscale,sizemc,extract
     &            ,off2,f,aewald,alsq2,alsq2n)
        implicit none
        integer  ,value,intent(in):: n,nbloc,n_mscale,sizemc
        logical  ,value,intent(in):: extract
        real(t_p),value,intent(in):: f,aewald,alsq2n,alsq2,off2
        integer  ,device,intent(in)::mcorrect_ik(sizemc,2)
     &           ,loc(n),locp(n),ipole(n)
        real(t_p),device,intent(in):: x(n),y(n),z(n),rpole(13,n)
     &           ,mcorrect_scale(sizemc)
        real(t_p),device,intent(inout)::tem(3,*)
     &           ,vir_buff(6*RED_BUFF_SIZE)
        ener_rtyp,device,intent(inout)::em_buff(RED_BUFF_SIZE)
        mdyn_rtyp,device,intent(inout)::dem(3,nbloc)
        integer   i,j,k,iglob,kglob,kbis,ithread
        integer   ii,kk,iipole,kkpole,lot
        ener_rtyp e
        real(t_p) r2,rstat
        real(t_p) xr,yr,zr
        type(rpole_elt) ip,kp
        type(mdyn3_r) frc_r
        type(real3) ttmi,ttmk,frc
        real(t_p) mscale

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        do ii = ithread, n_mscale, blockDim%x*gridDim%x
           iipole = mcorrect_ik(ii,1)
           kkpole = mcorrect_ik(ii,2)
           mscale = mcorrect_scale(ii)
           iglob  = ipole(iipole)
           kglob  = ipole(kkpole)
           i      = loc(iglob)
           kbis   = loc(kglob)

           xr     = x(kglob) - x(iglob)
           yr     = y(kglob) - y(iglob)
           zr     = z(kglob) - z(iglob)
           call image_inl(xr,yr,zr)
           r2     = xr*xr + yr*yr + zr*zr
           if (r2.gt.off2) cycle       !apply cutoff

           ip%c   = rpole(01,iipole)
           ip%dx  = rpole(02,iipole)
           ip%dy  = rpole(03,iipole)
           ip%dz  = rpole(04,iipole)
           ip%qxx = rpole(05,iipole)
           ip%qxy = rpole(06,iipole)
           ip%qxz = rpole(07,iipole)
           ip%qyy = rpole(09,iipole)
           ip%qyz = rpole(10,iipole)
           ip%qzz = rpole(13,iipole)
           kp%c   = rpole(01,kkpole)
           kp%dx  = rpole(02,kkpole)
           kp%dy  = rpole(03,kkpole)
           kp%dz  = rpole(04,kkpole)
           kp%qxx = rpole(05,kkpole)
           kp%qxy = rpole(06,kkpole)
           kp%qxz = rpole(07,kkpole)
           kp%qyy = rpole(09,kkpole)
           kp%qyz = rpole(10,kkpole)
           kp%qzz = rpole(13,kkpole)

c          mscale = 0.0
c20        continue
           e = 0
           call mdyn3r_zr(frc_r)
           call  real3_zr(ttmi)
           call  real3_zr(ttmk)
           ! compute mpole one interaction
           call mpole1_couple(r2,xr,yr,zr,ip,kp,mscale,
     &                        aewald,f,alsq2n,alsq2,
     &                        e,frc,frc_r,ttmi,ttmk,.true.)

           lot   = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
           ! update energy
           rstat = atomicAdd( em_buff(lot), e )
           ! increment the virial due to pairwise Cartesian forces c
           rstat = atomicAdd(vir_buff(0*RED_BUFF_SIZE+lot),-xr*frc%x)
           rstat = atomicAdd(vir_buff(1*RED_BUFF_SIZE+lot),-yr*frc%x)
           rstat = atomicAdd(vir_buff(2*RED_BUFF_SIZE+lot),-zr*frc%x)
           rstat = atomicAdd(vir_buff(3*RED_BUFF_SIZE+lot),-yr*frc%y)
           rstat = atomicAdd(vir_buff(4*RED_BUFF_SIZE+lot),-zr*frc%y)
           rstat = atomicAdd(vir_buff(5*RED_BUFF_SIZE+lot),-zr*frc%z)
           ! increment force-based gradient and torque on first site
           rstat = atomicAdd( dem(1,i)    ,-frc_r%x )
           rstat = atomicAdd( dem(2,i)    ,-frc_r%y )
           rstat = atomicAdd( dem(3,i)    ,-frc_r%z )
           rstat = atomicAdd( dem(1,kbis) , frc_r%x )
           rstat = atomicAdd( dem(2,kbis) , frc_r%y )
           rstat = atomicAdd( dem(3,kbis) , frc_r%z )
           ! increment force-based gradient and torque on second site
           if (extract) then
              i      = locp(iglob)
              kbis   = locp(kglob)
           end if
           rstat = atomicAdd( tem(1,i)    , ttmi%x )
           rstat = atomicAdd( tem(2,i)    , ttmi%y )
           rstat = atomicAdd( tem(3,i)    , ttmi%z )
           rstat = atomicAdd( tem(1,kbis) , ttmk%x )
           rstat = atomicAdd( tem(2,kbis) , ttmk%y )
           rstat = atomicAdd( tem(3,kbis) , ttmk%z )

        end do
        end subroutine

        attributes(global) subroutine emrealShlg_scaling_cu
     &            (mcorrect_ik,mcorrect_scale,ipole,loc,locp,x,y,z,rpole
     &            ,dem,tem,em_buff,vir_buff
     &            ,n,nbloc,n_mscale,sizemc,extract
     &            ,r_cut,loff2,off2,sheal,f,aewald,alsq2,alsq2n)
        implicit none
        integer  ,value,intent(in):: n,nbloc,n_mscale,sizemc
        logical  ,value,intent(in):: extract
        real(t_p),value,intent(in):: f,aewald,alsq2n,alsq2
     &           ,off2,loff2,r_cut,sheal
        integer  ,device,intent(in)::mcorrect_ik(sizemc,2)
     &           ,loc(n),locp(n),ipole(n)
        real(t_p),device,intent(in):: x(n),y(n),z(n),rpole(13,n)
     &           ,mcorrect_scale(sizemc)
        real(t_p),device,intent(inout)::tem(3,*)
     &           ,vir_buff(6*RED_BUFF_SIZE)
        ener_rtyp,device,intent(inout)::em_buff(RED_BUFF_SIZE)
        mdyn_rtyp,device,intent(inout)::dem(3,nbloc)
        integer   i,j,k,iglob,kglob,kbis,ithread
        integer   ii,kk,iipole,kkpole,lot,ver,fea
        ener_rtyp e
        integer(1) mutik
        logical ugrp,ulamdyn
        real(t_p) r2,fgrp,elambda,delambdae_
        real(t_p) xr,yr,zr,mscale
        type(rpole_elt) ip,kp
        type(mdyn3_r) frc_r
        type(real3) ttmi,ttmk,frc
        real(t_p) rstat
        parameter(ugrp=.false.,ulamdyn=.false.
     &           ,fgrp=1.0
     &           ,ver =__use_ene__+__use_grd__+__use_vir__+__use_sca__
     &           ,fea =__use_shortRange__)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        do ii = ithread, n_mscale, blockDim%x*gridDim%x
           iipole = mcorrect_ik(ii,1)
           kkpole = mcorrect_ik(ii,2)
           mscale = mcorrect_scale(ii)
           iglob  = ipole(iipole)
           kglob  = ipole(kkpole)
           i      = loc(iglob)
           kbis   = loc(kglob)

           xr     = x(kglob) - x(iglob)
           yr     = y(kglob) - y(iglob)
           zr     = z(kglob) - z(iglob)
           call image_inl(xr,yr,zr)
           r2     = xr*xr + yr*yr + zr*zr
           if (r2.gt.off2) cycle       !apply cutoff

           ip%c   = rpole(01,iipole)
           ip%dx  = rpole(02,iipole)
           ip%dy  = rpole(03,iipole)
           ip%dz  = rpole(04,iipole)
           ip%qxx = rpole(05,iipole)
           ip%qxy = rpole(06,iipole)
           ip%qxz = rpole(07,iipole)
           ip%qyy = rpole(09,iipole)
           ip%qyz = rpole(10,iipole)
           ip%qzz = rpole(13,iipole)
           kp%c   = rpole(01,kkpole)
           kp%dx  = rpole(02,kkpole)
           kp%dy  = rpole(03,kkpole)
           kp%dz  = rpole(04,kkpole)
           kp%qxx = rpole(05,kkpole)
           kp%qxy = rpole(06,kkpole)
           kp%qxz = rpole(07,kkpole)
           kp%qyy = rpole(09,kkpole)
           kp%qyz = rpole(10,kkpole)
           kp%qzz = rpole(13,kkpole)

c          mscale = 0.0
c20        continue
           e = 0
           call  real3_zr(ttmi)
           call  real3_zr(ttmk)
           ! compute mpole one interaction
           call duo_mpole(r2,xr,yr,zr,ip,kp,mscale
     &                   ,r_cut,sheal,aewald,f,alsq2n,alsq2
     &                   ,ugrp,fgrp,ulamdyn,mutik,elambda
     &                   ,delambdae_,e,frc,ttmi,ttmk
     &                   ,ver,fea)

           lot   = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
           ! update energy
           rstat = atomicAdd( em_buff(lot), e )
           ! increment the virial due to pairwise Cartesian forces c
           rstat = atomicAdd(vir_buff(0*RED_BUFF_SIZE+lot),-xr*frc%x)
           rstat = atomicAdd(vir_buff(1*RED_BUFF_SIZE+lot),-yr*frc%x)
           rstat = atomicAdd(vir_buff(2*RED_BUFF_SIZE+lot),-zr*frc%x)
           rstat = atomicAdd(vir_buff(3*RED_BUFF_SIZE+lot),-yr*frc%y)
           rstat = atomicAdd(vir_buff(4*RED_BUFF_SIZE+lot),-zr*frc%y)
           rstat = atomicAdd(vir_buff(5*RED_BUFF_SIZE+lot),-zr*frc%z)
           ! increment force-based gradient and torque on first site
           call atomic_add_f1( dem(1,i)    , frc%x )
           call atomic_add_f1( dem(2,i)    , frc%y )
           call atomic_add_f1( dem(3,i)    , frc%z )
           call atomic_add_f1( dem(1,kbis) ,-frc%x )
           call atomic_add_f1( dem(2,kbis) ,-frc%y )
           call atomic_add_f1( dem(3,kbis) ,-frc%z )
           ! increment force-based gradient and torque on second site
           if (extract) then
              i      = locp(iglob)
              kbis   = locp(kglob)
           end if
           rstat = atomicAdd( tem(1,i)    , ttmi%x )
           rstat = atomicAdd( tem(2,i)    , ttmi%y )
           rstat = atomicAdd( tem(3,i)    , ttmi%z )
           rstat = atomicAdd( tem(1,kbis) , ttmk%x )
           rstat = atomicAdd( tem(2,kbis) , ttmk%y )
           rstat = atomicAdd( tem(3,kbis) , ttmk%z )

        end do
        end subroutine
      end module
#endif
