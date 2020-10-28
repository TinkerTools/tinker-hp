c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "echargecu" : driver for calculation of the point charge
c     energy and derivatives with respect to Cartesian coordinates on device
c
c
#ifndef TINKER_CUF
#define TINKER_CUF
#include "tinker_precision.h"
#include "tinker_cudart.h"
      module echargecu
        use utilcu  ,only: nproc,ndir,BLOCK_DIM
        use utilgpu ,only: real3,real6,real3_red,rpole_elt
     &              ,BLOCK_SIZE,RED_BUFF_SIZE,WARP_SIZE
        contains

#include "image.f.inc"
#include "midpointimage.f.inc"
#include "switch_respa.f.inc"
#include "pair_charge.f.inc"

        attributes(global) subroutine ecreal1d_core_cu
     &        ( iion, cglob, loc, ieblst, eblst
     &        , nionlocnlb, nionlocnlb_pair, nionbloc, n
     &        , x, y, z, pchg
     &        , off2, f, aewald, ebuffer
     &        , dec, ec_buff, vir_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &        , inter
#endif
     &        )
        implicit none
        integer,value,intent(in)::nionlocnlb,nionbloc,n
     &         ,nionlocnlb_pair
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
     &           ,off2,aewald,f,ebuffer
        integer,device,intent(in)::iion(*),cglob(*),loc(*)
     &                ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),pchg(*)
        real(t_p),device:: vir_buff(*)
        real(r_p),device:: dec(3,*),ec_buff(*)
#ifdef TINKER_DEBUG
        integer  ,device:: inter(*)
#endif

        integer ithread,iwarp,nwarp,ilane,klane,srclane
        integer ii,j,i,kbis
        integer iblock,idx,kdx
        integer iichg,iglob,icloc,kchg,kglob,kcloc
#ifdef TINKER_DEBUG
        integer kglob_,kdx_,istat
#endif
        integer location
        integer,parameter::no_scaling=0
        real(t_p) xk_,yk_,zk_,d2,fi
        real(t_p) ec_
        real(t_p) rstat,zero
        type(real3) posi,pos
        type(real3) frc
        type(real3_red) frc_i
        real(t_p)      ,shared::   fk(BLOCK_DIM)
        type(real3_red),shared::frc_k(BLOCK_DIM)
        type(real3)    ,shared:: posk(BLOCK_DIM)
        real(t_p) vir_(6)
        logical do_pair,same_block,accept_mid
        parameter(zero=0.0)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, nionlocnlb_pair, nwarp
           iblock = ieblst(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iichg   = cglob(idx)
           iglob   = iion (idx)
           i       = loc  (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           fi      = f*pchg(iichg)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kchg    = cglob(kdx)
           kglob   = iion (kdx)
           kbis    = loc  (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           fk  (threadIdx%x)    = pchg(kchg)

           ! zero data to compute
           frc_i%x = 0.0;
           frc_i%y = 0.0;
           frc_i%z = 0.0;
           frc_k(threadIdx%x)%x = 0.0;
           frc_k(threadIdx%x)%y = 0.0;
           frc_k(threadIdx%x)%z = 0.0;

           !* set compute Data to 0
           ec_ = 0.0
           vir_(1)=0.0; vir_(2)=0.0; vir_(3)=0.0;
           vir_(4)=0.0; vir_(5)=0.0; vir_(6)=0.0;

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_     = __shfl(kdx ,srclane)
              kglob_   = __shfl(kglob ,srclane)
#endif
              if (ndir.gt.1) then
                 xk_   = posk(klane)%x
                 yk_   = posk(klane)%y
                 zk_   = posk(klane)%z
                 pos%x = posi%x - xk_
                 pos%y = posi%y - yk_
                 pos%z = posi%z - zk_
                 call midpointimage_inl(xk_,yk_,zk_,pos%x,pos%y,pos%z)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 pos%x = posi%x - posk(klane)%x 
                 pos%y = posi%y - posk(klane)%y 
                 pos%z = posi%z - posk(klane)%z 
                 call image_inl(pos%x,pos%y,pos%z)
              end if

              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                        ,same_block)

              if (do_pair.and.d2<=off2.and.accept_mid) then
                 ! compute one interaction
                 call charge_couple(d2,pos%x,pos%y,pos%z,ebuffer
     &                             ,fi*fk(klane),aewald,1.0_ti_p
     &                             ,ec_,frc,frc_k(klane),no_scaling)

                 vir_(1) = vir_(1) + pos%x * frc%x
                 vir_(2) = vir_(2) + pos%y * frc%x
                 vir_(3) = vir_(3) + pos%z * frc%x
                 vir_(4) = vir_(4) + pos%y * frc%y
                 vir_(5) = vir_(5) + pos%z * frc%y
                 vir_(6) = vir_(6) + pos%z * frc%z

                 frc_i%x = frc_i%x + frc%x
                 frc_i%y = frc_i%y + frc%y
                 frc_i%z = frc_i%z + frc%z

#ifdef TINKER_DEBUG
                 if (iglob<kglob_) then
                    istat=AtomicAdd(inter(iglob),1)
                 else
                    istat=AtomicAdd(inter(kglob_),1)
                 end if
#endif
              end if
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

           ! Update energy buffer
           rstat = atomicAdd(ec_buff(location), real(ec_,r_p))

           ! Update forces
           rstat = atomicAdd( dec(1,i   ), frc_i%x )
           rstat = atomicAdd( dec(2,i   ), frc_i%y )
           rstat = atomicAdd( dec(3,i   ), frc_i%z )
           rstat = atomicAdd( dec(1,kbis), frc_k(threadIdx%x)%x )
           rstat = atomicAdd( dec(2,kbis), frc_k(threadIdx%x)%y )
           rstat = atomicAdd( dec(3,kbis), frc_k(threadIdx%x)%z )

           ! Update virial buffer
           rstat = atomicAdd(vir_buff(0*RED_BUFF_SIZE+location),vir_(1))
           rstat = atomicAdd(vir_buff(1*RED_BUFF_SIZE+location),vir_(2))
           rstat = atomicAdd(vir_buff(2*RED_BUFF_SIZE+location),vir_(3))
           rstat = atomicAdd(vir_buff(3*RED_BUFF_SIZE+location),vir_(4))
           rstat = atomicAdd(vir_buff(4*RED_BUFF_SIZE+location),vir_(5))
           rstat = atomicAdd(vir_buff(5*RED_BUFF_SIZE+location),vir_(6))

        end do

        end subroutine

        attributes(global) subroutine ecreal3d_core_cu
     &        ( iion, cglob, loc, ieblst, eblst
     &        , nionlocnlb, nionlocnlb_pair, nionbloc, n
     &        , x, y, z, pchg
     &        , off2, f, aewald, ebuffer
     &        , ec_buff, nec_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &        )
        implicit none
        integer,value,intent(in)::nionlocnlb,nionbloc,n
     &         ,nionlocnlb_pair
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
     &           ,off2,aewald,f,ebuffer
        integer,device,intent(in)::iion(*),cglob(*),loc(*)
     &                ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),pchg(*)
        integer  ,device:: nec_buff(*)
        real(r_p),device::  ec_buff(*)
#ifdef TINKER_DEBUG
#endif

        integer ithread,iwarp,nwarp,ilane,klane,srclane
        integer ii,j,i,kbis
        integer iblock,idx,kdx
        integer iichg,iglob,icloc,kchg,kglob,kcloc
#ifdef TINKER_DEBUG
        integer kglob_,kdx_,istat
#endif
        integer location
        integer,parameter::no_scaling=0
        integer nec_,istat
        real(t_p) xk_,yk_,zk_,d2,fi
        real(t_p) ec_
        real(t_p) rstat,zero
        type(real3) posi,pos
        real(t_p)      ,shared::   fk(BLOCK_DIM)
        type(real3)    ,shared:: posk(BLOCK_DIM)
        logical do_pair,same_block,accept_mid
        parameter(zero=0.0)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, nionlocnlb_pair, nwarp
           iblock = ieblst(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iichg   = cglob(idx)
           iglob   = iion (idx)
           i       = loc  (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           fi      = f*pchg(iichg)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kchg    = cglob(kdx)
           kglob   = iion (kdx)
           kbis    = loc  (kdx)
           posk(threadIdx%x)%x = x(kdx)
           posk(threadIdx%x)%y = y(kdx)
           posk(threadIdx%x)%z = z(kdx)
           fk  (threadIdx%x)   = pchg(kchg)

           !* set compute Data to 0
           ec_  = 0.0
           nec_ = 0

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_     = __shfl(kdx ,srclane)
              kglob_   = __shfl(kglob ,srclane)
#endif
              if (ndir.gt.1) then
                 xk_   = posk(klane)%x
                 yk_   = posk(klane)%y
                 zk_   = posk(klane)%z
                 pos%x = posi%x - xk_
                 pos%y = posi%y - yk_
                 pos%z = posi%z - zk_
                 call midpointimage_inl(xk_,yk_,zk_,pos%x,pos%y,pos%z)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    accept_mid = .false.
                 else
                    accept_mid = .true.
                 end if
              else
                 pos%x = posi%x - posk(klane)%x 
                 pos%y = posi%y - posk(klane)%y 
                 pos%z = posi%z - posk(klane)%z 
                 call image_inl(pos%x,pos%y,pos%z)
              end if

              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                        ,same_block)

              if (do_pair.and.d2<=off2.and.accept_mid) then
                 ! compute one interaction
                 call charge3_couple(d2,pos%x,pos%y,pos%z,ebuffer
     &                             ,fi*fk(klane),aewald,1.0_ti_p
     &                             ,ec_,no_scaling)
                 nec_ = nec_ + 1

#ifdef TINKER_DEBUG
#endif
              end if
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

           ! Update energy buffer
           rstat = atomicAdd( ec_buff(location), real(ec_,r_p))
           istat = atomicAdd(nec_buff(location), nec_)
        end do

        end subroutine

      end module
#endif
