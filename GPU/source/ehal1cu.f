c
c     sorbonne university
c     washington university in saint louis
c     university of texas at austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1cu  --  buffered 14-7 energy & derivatives##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1" calculates the buffered 14-7 van der waals energy and
c     its first derivatives with respect to cartesian coordinates
c
c     *_t texture variable destined to be attached to their target
#ifdef _CUDA
#define TINKER_CUF
#include  "tinker_macro.h"
      module ehal1cu
        use cudafor
        use neigh    ,only: i_bl=>inside_block,n_bl=>normal_block
     &               ,d_bl=>disc_block
        use sizes    ,only: maxvalue,maxclass,maxgrp
        use tinheader,only: ti_p
        use utilcu   ,only: nproc,ngrp,all_lanes,VDW_BLOCK_DIM
     &               ,skipvdw12,cu_update_skipvdw12,use_virial
     &               ,vcouple,xcell2,ycell2,zcell2,f_abs
#if TINKER_SINGLE_PREC + TINKER_MIXED_PREC
     &               ,f_sqrt
#endif
        use utilgpu  ,only: BLOCK_SIZE,RED_BUFF_SIZE
        use vdw      ,only: vdw_skipvdw12=>skipvdw12
        use vdwpot   ,only: ARITHMETIC_RL,GEOMETRIC_RL,CUBIC_MEAN_RL
     &               ,HARMONIC_RL,HHG_RL

        integer(1),private:: one1,two1
        integer  ,pointer,device::ired_t(:)
     &           ,loc_ired_t(:),vblst_t(:),ivblst_t(:)
     &           ,jvdw_t(:),i12_p(:,:)
        real(t_p),pointer,texture::radmin_t(:,:),epsilon_t(:,:)
        real(t_p),pointer,device::
     &            kred_t(:),xred_t(:),yred_t(:),zred_t(:)
        parameter(one1=1, two1=2)
#include "atomicOp.h.f"
        contains

#include "convert.f.inc"
#include "image.f.inc"
#include "midpointimage.f.inc"
#include "pair_ehal.inc.f"
#include "groups.inc.f"
#include "atomicOp.inc.f"

        M_subroutine 
     &        get_rad(radi,radk,rv,radrule)
        implicit none
        integer radrule
        real(t_p),radi,radk,rv
        if (radrule.eq.ARITHMETIC_RL) then
        rv = (radi+radk)
        else if (radrule.eq.GEOMETRIC_RL) then
        rv = 2.0*f_sqrt(radi*radk)
        else if (radrule.eq.CUBIC_MEAN_RL) then
        rv = 2.0*(radi**3+radk**3)/(radi**2+radk**2)
        end if
        end subroutine
        M_subroutine 
     &        get_eps(epsi,epsk,ep,epsrule)
        implicit none
        integer epsrule
        real(t_p) epsi,epsk,ep
        if      (epsrule.eq.ARITHMETIC_RL) then
        ep = 0.5*(epsi+epsk)
        else if (epsrule.eq.GEOMETRIC_RL) then
        ep = f_sqrt(epsi*epsk)
        else if (epsrule.eq.HHG_RL) then
        ep = 4.0*(epsi*epsk)/(f_sqrt(epsi)+f_sqrt(epsk))**2
        end if
        end subroutine

c       All available features for hal kernel
c#define __tfea__ (__use_mpi__+__use_softcore__+__radepsOpt__+__use_shortRange__+__use_longRange__+__use_groups__)
c#define __tver__ (__use_ene__+__use_grd__+__use_vir__+__use_act__)

#define __tver__ (__use_ene__+__use_grd__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_softcore__+__use_groups__+__use_lambdadyn__)
#define __sufx__ 1_kcu
#include "ehalcu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ (__use_grd__)
#define __tfea__ (0)
#define __sufx__ 1_v0_kcu
#include "ehalcu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ (__use_ene__+__use_grd__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_softcore__+__use_groups__+__use_lambdadyn__+__use_shortRange__)
#define __sufx__ 1short_kcu
#include "ehalcu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ (__use_ene__+__use_grd__+__use_vir__)
#define __tfea__ (__use_mpi__+__use_softcore__+__use_groups__+__use_lambdadyn__+__use_longRange__)
#define __sufx__ 1long_kcu
#include "ehalcu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

#define __tver__ (__use_ene__+__use_act__)
#define __tfea__ (__use_mpi__+__use_softcore__+__use_groups__)
#define __sufx__ 3_kcu
#include "ehalcu.tpl.f"
#undef __tver__
#undef __tfea__
#undef __sufx__

#if 1
        attributes(global)
     &  subroutine ehal1short_cu_deb
     &           (xred,yred,zred,sgl_id,slc_id,loc_ired
     &           ,ivblst,vblst,jvdw,epsilon,radmin
     &           ,ired,kred,grplist,wgrp,dev,ev_buff,vir_buff
     &           ,nvdwlocnlb_pair,n,nbloc,nvdwlocnl,nvdwlocnlb
     &           ,nvdwclass,rinv
     &           ,c0,c1,c2,c3,c4,c5,cut2,off2,off
     &           ,scexp,vlambda,scalpha,mut
     &           ,shortheal,ghal,dhal,use_short,use_group
     &           ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &           ,int1,int2,rank
#ifdef TINKER_DEBUG
     &           ,inter
#endif
     &           )

        implicit none
        integer  ,value,intent(in):: nvdwlocnlb_pair,n,nbloc
     &           ,nvdwlocnl,nvdwlocnlb,nvdwclass,int1,int2,rank
        real(t_p),value,intent(in):: c0,c1,c2,c3,c4,c5
     &           ,cut2,ghal,dhal,off2,off,shortheal
     &           ,scexp,vlambda,scalpha
        integer(1),device :: mut(n)
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,rinv
        ener_rtyp,device :: ev_buff (RED_BUFF_SIZE)
        real(t_p),device :: vir_buff(RED_BUFF_SIZE)
        integer  ,device,intent(in)::sgl_id(nvdwlocnlb)
     &           ,loc_ired(nvdwlocnlb),ivblst(nvdwlocnlb_pair)
     &           ,vblst(nvdwlocnlb_pair*(BLOCK_SIZE))
     &           ,slc_id(nvdwlocnlb)
        integer  ,device,intent(in)::ired(n),jvdw(nvdwlocnlb)
        real(t_p),device,intent(in)::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass),kred(n)
        real(t_p),device,intent(in):: xred(nvdwlocnlb)
     &           ,yred(nvdwlocnlb),zred(nvdwlocnlb)
        mdyn_rtyp,device,intent(inout)::dev(3,nbloc)
        logical  ,value,intent(in)::use_short
        logical, value, intent(in) :: use_group
        integer, device, intent(in) :: grplist(*)
        real(t_p), device, intent(in) :: wgrp(ngrp+1,ngrp+1)
#ifdef TINKER_DEBUG
        integer,device:: inter(*)
#endif

        integer iga,igb,igb_
        integer ithread,iwarp,nwarp,ilane,srclane
        integer dstlane,klane,iblock
        integer idx,kdx,kdx_,kglob_,ii,j,i,iglob,it
        integer kbis,ist,k
        integer,shared,dimension(VDW_BLOCK_DIM):: kt,kglob
        integer ai12(maxvalue)
        real(t_p) e,istat
        ener_rtyp ev_
        real(t_p) xi,yi,zi,xk_,yk_,zk_,xpos,ypos,zpos
        real(t_p),shared,dimension(VDW_BLOCK_DIM)::xk,yk,zk
        real(t_p) rik2,rv2,eps2,scale_
        real(t_p) dedx,dedy,dedz,devx,devy,devz
        mdyn_rtyp,shared,dimension(VDW_BLOCK_DIM)::gxk,gyk,gzk
        mdyn_rtyp gxi,gyi,gzi
        real(t_p) vxx_,vxy_,vxz_,vyy_,vyz_,vzz_
        logical do_pair,same_block,accept_mid,ik12
        integer(1) muti,mutik
        integer(1),shared::mutk(VDW_BLOCK_DIM)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, nvdwlocnlb_pair-1, nwarp

           ! Load atom block k neighbor parameter
           kdx    = vblst( ii*warpsize + ilane )
           kglob(threadIdx%x)  = sgl_id(kdx)
           kbis   = slc_id(kdx)
           kt(threadIdx%x) = jvdw  (kdx)
           xk(threadIdx%x) = xred  (kdx)
           yk(threadIdx%x) = yred  (kdx)
           zk(threadIdx%x) = zred  (kdx)
#ifndef TINKER_NO_MUTATE
           mutk(threadIdx%x) = mut(kglob(threadIdx%x))
#endif
           call syncwarp(ALL_LANES)

           ! Set Data to compute to zero
           ev_    = 0
           gxi               = 0
           gyi               = 0
           gzi               = 0
           gxk(threadIdx%x)  = 0
           gyk(threadIdx%x)  = 0
           gzk(threadIdx%x)  = 0
           vxx_   = 0
           vxy_   = 0
           vxz_   = 0
           vyy_   = 0
           vyz_   = 0
           vzz_   = 0

           ! Load atom block i parameters
           iblock = ivblst(ii+1)
           if (iblock.eq.0) cycle
           idx    = (iblock-1)*warpsize + ilane 
           iglob  = sgl_id(idx)
           i      = slc_id (idx)
           it     = jvdw  (idx)
           xi     = xred  (idx)
           yi     = yred  (idx)
           zi     = zred  (idx)
#ifndef TINKER_NO_MUTATE
           muti   = mut(iglob)
#endif
           if (skipvdw12) then
             do k = 1,maxvalue
                ai12(k) = i12_p(k,iglob)
             end do
           end if

           same_block = (idx.ne.kdx)

           ! load group factors
            if (use_group) then
              iga=grplist(iglob)
              igb=grplist(kglob(threadIdx%x))
            endif

           ! Interact block i with block k
           do j = 0,warpsize-1
              srclane = iand( ilane+j-1,warpsize-1 ) + 1
              klane   = threadIdx%x-ilane + srclane
              kglob_  = kglob(klane)
#ifdef TINKER_DEBUG
              kdx_    = __shfl(kdx  ,srclane)
#endif
#ifndef TINKER_NO_MUTATE
              mutik   = muti + mutk(klane)
#endif
              if(use_group) then
                igb_ = __shfl(igb,klane)
                scale_ = wgrp(iga+1,igb_+1)
              else
                scale_  = 1.0_ti_p
              endif

              if (skipvdw12) then
                 ! Look for 1-2 bond interaction
                 ! and skip it
                 ik12 =.false.
                 do k = 1, maxvalue
                    if (ai12(k).eq.kglob_) ik12=.true.
                 end do
                 if (ik12) cycle
              end if

              if (nproc.gt.1) then
                 xk_   = xk(klane)
                 yk_   = yk(klane)
                 zk_   = zk(klane)
                 xpos  = xi - xk_
                 ypos  = yi - yk_
                 zpos  = zi - zk_
                 call midpointimage_inl(xk_,yk_,zk_,xpos,ypos,zpos)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 xpos    = xi - xk(klane)
                 ypos    = yi - yk(klane)
                 zpos    = zi - zk(klane)
                 call image_inl(xpos,ypos,zpos)
              end if

              dedx = 0.0; dedy = 0.0; dedz = 0.0;

              rik2  = xpos**2 + ypos**2 + zpos**2

              do_pair = merge(.true.,iglob.lt.kglob_,same_block)
              if (do_pair.and.rik2<=off2
     &           .and.accept_mid) then

#ifndef TINKER_NO_MUTATE
                 ! Annihilate
                 if (vcouple.and.mutik.eq.two1) mutik=one1
#endif

                 rv2   =  radmin_t (kt(klane),it)
                 eps2  = epsilon_t (kt(klane),it)

c                call ehal1_couple_short(xpos,ypos,zpos,rik2,rv2,eps2
c    &                     ,1.0_ti_p,cut2,off
c    &                     ,scexp,vlambda,scalpha,mutik
c    &                     ,shortheal,ghal,dhal,e,dedx,dedy,dedz)

                 call ehal1_couple(xpos,ypos,zpos,rik2,rv2,eps2,1.0_ti_p
     &                         ,cut2,rinv,off,ghal,dhal
     &                         ,scexp,vlambda,scalpha,mutik
     &                         ,e,dedx,dedy,dedz
     &                         ,__use_ene__,__use_softcore__)

                 ev_   = ev_ + tp2enr(e)

#ifdef TINKER_DEBUG
                 if (iglob<kglob(klane)) then
                    ist = Atomicadd(inter(iglob),1)
                 else
                    ist = Atomicadd(inter(kglob(klane)),1)
                 end if
#endif
                 if ((iglob.eq.int1.or.kglob_.eq.int1)
     &              .and.abs(dedz).gt.1e-7) then
                    print*,min(iglob,kglob_),max(iglob,kglob_)
     &                    ,real(rik2,4),real(dedz,4),'i1',rank
                 end if
                 if ((iglob.eq.int2.or.kglob_.eq.int2)
     &              .and.abs(dedz).gt.1e-7) then
                    print*,min(iglob,kglob_),max(iglob,kglob_)
     &                    ,real(rik2,4),real(dedz,4),'i2',rank
                 end if

                 if (use_virial) then
                 vxx_  = vxx_  + xpos * dedx
                 vxy_  = vxy_  + ypos * dedx
                 vxz_  = vxz_  + zpos * dedx
                 vyy_  = vyy_  + ypos * dedy
                 vyz_  = vyz_  + zpos * dedy
                 vzz_  = vzz_  + zpos * dedz
                 end if

                 !dstlane = iand( ilane-1+warpsize-j, warpsize-1 ) + 1

                 ! Resolve gradient
                 gxi = gxi + tp2mdr(dedx)
                 gyi = gyi + tp2mdr(dedy)
                 gzi = gzi + tp2mdr(dedz)

                 gxk(klane) = gxk(klane) + tp2mdr(dedx)
                 gyk(klane) = gyk(klane) + tp2mdr(dedy)
                 gzk(klane) = gzk(klane) + tp2mdr(dedz)

              end if

           end do

           it = iand(ithread-1,RED_BUFF_SIZE-1) + 1
           !increment the van der Waals energy
           istat = atomicAdd( ev_buff(it) ,ev_ )
 
c          if (use_virial) then
c          ! Increment virial term of van der Waals
c          istat = atomicAdd( vir_buff(0*RED_BUFF_SIZE+it),vxx_ )
c          istat = atomicAdd( vir_buff(1*RED_BUFF_SIZE+it),vxy_ )
c          istat = atomicAdd( vir_buff(2*RED_BUFF_SIZE+it),vxz_ )
c          istat = atomicAdd( vir_buff(3*RED_BUFF_SIZE+it),vyy_ )
c          istat = atomicAdd( vir_buff(4*RED_BUFF_SIZE+it),vyz_ )
c          istat = atomicAdd( vir_buff(5*RED_BUFF_SIZE+it),vzz_ )
c          end if

c          !increment the van der Waals derivatives
c          if (idx.le.nvdwlocnl) then
c             istat   = atomicAdd (dev(1,i),gxi)
c             istat   = atomicAdd (dev(2,i),gyi)
c             istat   = atomicAdd (dev(3,i),gzi)
c          end if

c          call syncwarp(ALL_LANES)
c          if (kdx.le.nvdwlocnl) then
c             istat   = atomicSub (dev(1,kbis),gxk(threadIdx%x))
c             istat   = atomicSub (dev(2,kbis),gyk(threadIdx%x))
c             istat   = atomicSub (dev(3,kbis),gzk(threadIdx%x))
c          end if

        end do
        end subroutine
#endif

        subroutine attach_vdwcu_data(i12,kred,radmin,epsilon
     &                   ,n,nvdwclass)
        implicit none
        integer  ,intent(in):: n,nvdwclass
        integer  ,intent(in),device,target:: i12(maxvalue,n)
        real(t_p),intent(in),device,target:: kred(n)
        real(t_p),device,target::radmin(nvdwclass,nvdwclass)
     &           ,epsilon(nvdwclass,nvdwclass)
        kred_t    => kred
        i12_p     => i12
        radmin_t  => radmin
        epsilon_t => epsilon
        end subroutine

      end module

      subroutine ehal1c_kernel
      use ehal1cu
      implicit none
      end subroutine
#else
      ! For Portability issue
      subroutine void_ehalcu
      end
#endif
