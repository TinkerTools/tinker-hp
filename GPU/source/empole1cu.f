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
        use utilcu  ,only: nproc,BLOCK_DIM
        use utilgpu ,only: BLOCK_SIZE,RED_BUFF_SIZE,WARP_SIZE
        use tinTypes,only: real3,real6,mdyn3_r,rpole_elt 

        contains

#include "convert.f.inc"
#include "image.f.inc"
#include "midpointimage.f.inc"
#include "pair_mpole1.f.inc"

        attributes(global) subroutine emreal1c_core_cu
     &        ( ipole, pglob, loc, ieblst, eblst
     &        , npolelocnlb, npolelocnlb_pair, npolebloc, n
     &        , x, y, z, rpole
     &        , off2, f, alsq2, alsq2n, aewald
     &        , dem, tem, em_buff, vir_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &        , inter
#endif
     &        )
        implicit none
        integer,value,intent(in)::npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
     &           ,off2,alsq2,alsq2n,aewald,f
        integer,device,intent(in)::ipole(*),pglob(*),loc(*)
     &                ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),rpole(13,*)
        real(t_p),device:: tem(3,*),vir_buff(*)
        ener_rtyp,device:: dem(3,*)
        ener_rtyp,device:: em_buff(*)
#ifdef TINKER_DEBUG
        integer  ,device:: inter(*)
#endif

        integer ithread,iwarp,nwarp,ilane,klane,srclane
        integer ii,j,i,kbis
        integer iblock,idx,kdx
        integer iipole,iglob,iploc,kpole,kglob,kploc
#ifdef TINKER_DEBUG
        integer kglob_,kdx_,istat
#endif
        integer location
        real(t_p) xk_,yk_,zk_,d2
        ener_rtyp em_
        real(t_p) rstat,zero
        type(real3) posi,pos
        type(real3) frc
        type(mdyn3_r) frc_i
        type(mdyn3_r),shared::frc_k(BLOCK_DIM)
        type(real3),shared::posk(BLOCK_DIM)
        type(real3) ttmi
        type(real3),shared:: ttmk(BLOCK_DIM)
        real(t_p) vir_(6)
        type(rpole_elt) ip,kp,kp_
        logical do_pair,same_block,accept_mid
        parameter(zero=0.0)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, npolelocnlb_pair, nwarp
           iblock = ieblst(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           i       = loc  (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           ip%c    = rpole( 1, iipole)
           ip%dx   = rpole( 2, iipole)
           ip%dy   = rpole( 3, iipole)
           ip%dz   = rpole( 4, iipole)
           ip%qxx  = rpole( 5, iipole)
           ip%qxy  = rpole( 6, iipole)
           ip%qxz  = rpole( 7, iipole)
           ip%qyy  = rpole( 9, iipole)
           ip%qyz  = rpole(10, iipole)
           ip%qzz  = rpole(13, iipole)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob   = ipole(kdx)
           kbis    = loc  (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kp%c    = rpole( 1, kpole)
           kp%dx   = rpole( 2, kpole)
           kp%dy   = rpole( 3, kpole)
           kp%dz   = rpole( 4, kpole)
           kp%qxx  = rpole( 5, kpole)
           kp%qxy  = rpole( 6, kpole)
           kp%qxz  = rpole( 7, kpole)
           kp%qyy  = rpole( 9, kpole)
           kp%qyz  = rpole(10, kpole)
           kp%qzz  = rpole(13, kpole)

           ! zero data to compute

           frc_i%x = 0;
           frc_i%y = 0;
           frc_i%z = 0;
           frc_k(threadIdx%x)%x = 0;
           frc_k(threadIdx%x)%y = 0;
           frc_k(threadIdx%x)%z = 0;
           ttmi%x  = 0.0;
           ttmi%y  = 0.0;
           ttmi%z  = 0.0;
           ttmk(threadIdx%x)%x = 0.0;
           ttmk(threadIdx%x)%y = 0.0;
           ttmk(threadIdx%x)%z = 0.0;

           !* set compute Data to 0
           em_ = 0
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
              kp_%c    = __shfl(kp%c   ,klane)
              kp_%dx   = __shfl(kp%dx  ,klane)
              kp_%dy   = __shfl(kp%dy  ,klane)
              kp_%dz   = __shfl(kp%dz  ,klane)
              kp_%qxx  = __shfl(kp%qxx ,klane)
              kp_%qxy  = __shfl(kp%qxy ,klane)
              kp_%qxz  = __shfl(kp%qxz ,klane)
              kp_%qyy  = __shfl(kp%qyy ,klane)
              kp_%qyz  = __shfl(kp%qyz ,klane)
              kp_%qzz  = __shfl(kp%qzz ,klane)

              if (nproc.gt.1) then
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
                    pos%x=-pos%x; pos%y=-pos%y; pos%z=-pos%z;
                 end if
              else
                 pos%x = posk(klane)%x - posi%x
                 pos%y = posk(klane)%y - posi%y
                 pos%z = posk(klane)%z - posi%z
                 call image_inl(pos%x,pos%y,pos%z)
              end if

              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                        ,same_block)

              if (do_pair.and.d2<=off2.and.accept_mid) then
                 ! compute one interaction
                 call mpole1_couple(d2,pos%x,pos%y,pos%z,ip,kp_,zero
     &                             ,aewald,f,alsq2n,alsq2
     &                            ,em_,frc,frc_k(klane),ttmi,ttmk(klane)
     &                             ,.false.)

                 vir_(1) = vir_(1) - pos%x * frc%x
                 vir_(2) = vir_(2) - pos%y * frc%x
                 vir_(3) = vir_(3) - pos%z * frc%x
                 vir_(4) = vir_(4) - pos%y * frc%y
                 vir_(5) = vir_(5) - pos%z * frc%y
                 vir_(6) = vir_(6) - pos%z * frc%z

                 frc_i%x = frc_i%x + tp2mdr(frc%x)
                 frc_i%y = frc_i%y + tp2mdr(frc%y)
                 frc_i%z = frc_i%z + tp2mdr(frc%z)

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
           rstat = atomicAdd(em_buff(location), em_)

           ! Update forces
           rstat = atomicAdd( dem(1,i   ), frc_i%x )
           rstat = atomicAdd( dem(2,i   ), frc_i%y )
           rstat = atomicAdd( dem(3,i   ), frc_i%z )
           rstat = atomicAdd( dem(1,kbis), frc_k(threadIdx%x)%x )
           rstat = atomicAdd( dem(2,kbis), frc_k(threadIdx%x)%y )
           rstat = atomicAdd( dem(3,kbis), frc_k(threadIdx%x)%z )

           ! Update torque
           rstat = atomicAdd( tem(1,i)   ,ttmi%x )
           rstat = atomicAdd( tem(2,i)   ,ttmi%y )
           rstat = atomicAdd( tem(3,i)   ,ttmi%z )
           rstat = atomicAdd( tem(1,kbis),ttmk(threadIdx%x)%x )
           rstat = atomicAdd( tem(2,kbis),ttmk(threadIdx%x)%y )
           rstat = atomicAdd( tem(3,kbis),ttmk(threadIdx%x)%z )

           ! Update virial buffer
           rstat = atomicAdd(vir_buff(0*RED_BUFF_SIZE+location),vir_(1))
           rstat = atomicAdd(vir_buff(1*RED_BUFF_SIZE+location),vir_(2))
           rstat = atomicAdd(vir_buff(2*RED_BUFF_SIZE+location),vir_(3))
           rstat = atomicAdd(vir_buff(3*RED_BUFF_SIZE+location),vir_(4))
           rstat = atomicAdd(vir_buff(4*RED_BUFF_SIZE+location),vir_(5))
           rstat = atomicAdd(vir_buff(5*RED_BUFF_SIZE+location),vir_(6))

        end do

        end subroutine

        attributes(global) subroutine emrealshortlong1c_core_cu
     &        ( ipole, pglob, loc, ieblst, eblst
     &        , npolelocnlb, npolelocnlb_pair, npolebloc, n
     &        , x, y, z, rpole
     &        , shortheal, r_cut, sh_cut2, off2, f
     &        , alsq2, alsq2n, aewald, mode
     &        , dem, tem, em_buff, vir_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
     &        , inter
#endif
     &        )
        implicit none
        integer,value,intent(in)::npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair,mode
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,sh_cut2,shortheal
     &           ,r_cut,off2,alsq2,alsq2n,aewald,f
        integer,device,intent(in)::ipole(*),pglob(*),loc(*)
     &                ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),rpole(13,*)
        real(t_p),device:: tem(3,*),vir_buff(*)
        mdyn_rtyp,device:: dem(3,*)
        ener_rtyp,device:: em_buff(*)
#ifdef TINKER_DEBUG
        integer  ,device:: inter(*)
#endif

        integer ithread,iwarp,nwarp,ilane,klane,srclane
        integer ii,j,i,kbis
        integer iblock,idx,kdx
        integer iipole,iglob,iploc,kpole,kglob,kploc
#ifdef TINKER_DEBUG
        integer kglob_,kdx_,istat
#endif
        integer location
        real(t_p) xk_,yk_,zk_,d2
        ener_rtyp em_
        real(t_p) rstat,zero
        type(real3) posi,pos
        type(real3) frc
        type(mdyn3_r) frc_i
        type(mdyn3_r),shared::frc_k(BLOCK_DIM)
        type(real3) ,shared::posk(BLOCK_DIM)
        type(real3) ttmi
        type(real3) ,shared:: ttmk(BLOCK_DIM)
        real(t_p) vir_(6)
        type(rpole_elt) ip,kp,kp_
        logical do_pair,same_block,accept_mid
        parameter(zero=0.0)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, npolelocnlb_pair, nwarp
           iblock = ieblst(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           i       = loc  (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           ip%c    = rpole( 1, iipole)
           ip%dx   = rpole( 2, iipole)
           ip%dy   = rpole( 3, iipole)
           ip%dz   = rpole( 4, iipole)
           ip%qxx  = rpole( 5, iipole)
           ip%qxy  = rpole( 6, iipole)
           ip%qxz  = rpole( 7, iipole)
           ip%qyy  = rpole( 9, iipole)
           ip%qyz  = rpole(10, iipole)
           ip%qzz  = rpole(13, iipole)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob   = ipole(kdx)
           kbis    = loc  (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kp%c    = rpole( 1, kpole)
           kp%dx   = rpole( 2, kpole)
           kp%dy   = rpole( 3, kpole)
           kp%dz   = rpole( 4, kpole)
           kp%qxx  = rpole( 5, kpole)
           kp%qxy  = rpole( 6, kpole)
           kp%qxz  = rpole( 7, kpole)
           kp%qyy  = rpole( 9, kpole)
           kp%qyz  = rpole(10, kpole)
           kp%qzz  = rpole(13, kpole)

           ! zero data to compute

           frc_i%x = 0;
           frc_i%y = 0;
           frc_i%z = 0;
           frc_k(threadIdx%x)%x = 0;
           frc_k(threadIdx%x)%y = 0;
           frc_k(threadIdx%x)%z = 0;
           ttmi%x  = 0.0;
           ttmi%y  = 0.0;
           ttmi%z  = 0.0;
           ttmk(threadIdx%x)%x = 0.0;
           ttmk(threadIdx%x)%y = 0.0;
           ttmk(threadIdx%x)%z = 0.0;

           !* set compute Data to 0
           em_ = 0
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
              kp_%c    = __shfl(kp%c   ,klane)
              kp_%dx   = __shfl(kp%dx  ,klane)
              kp_%dy   = __shfl(kp%dy  ,klane)
              kp_%dz   = __shfl(kp%dz  ,klane)
              kp_%qxx  = __shfl(kp%qxx ,klane)
              kp_%qxy  = __shfl(kp%qxy ,klane)
              kp_%qxz  = __shfl(kp%qxz ,klane)
              kp_%qyy  = __shfl(kp%qyy ,klane)
              kp_%qyz  = __shfl(kp%qyz ,klane)
              kp_%qzz  = __shfl(kp%qzz ,klane)

              if (nproc.gt.1) then
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
                    pos%x=-pos%x; pos%y=-pos%y; pos%z=-pos%z;
                 end if
              else
                 pos%x = posk(klane)%x - posi%x
                 pos%y = posk(klane)%y - posi%y
                 pos%z = posk(klane)%z - posi%z
                 call image_inl(pos%x,pos%y,pos%z)
              end if

              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                        ,same_block)

              if (do_pair.and.sh_cut2<=d2.and.d2<=off2.and.accept_mid)
     &           then
                 ! compute one interaction
                 call mpole1_couple_shortlong
     &                          (d2,pos%x,pos%y,pos%z,ip,kp_,zero
     &                          ,r_cut,shortheal,aewald,f,alsq2n,alsq2
     &                          ,em_,frc,frc_k(klane),ttmi,ttmk(klane)
     &                          ,.false.,mode)

                 vir_(1) = vir_(1) - pos%x * frc%x
                 vir_(2) = vir_(2) - pos%y * frc%x
                 vir_(3) = vir_(3) - pos%z * frc%x
                 vir_(4) = vir_(4) - pos%y * frc%y
                 vir_(5) = vir_(5) - pos%z * frc%y
                 vir_(6) = vir_(6) - pos%z * frc%z

                 frc_i%x = frc_i%x + tp2mdr(frc%x)
                 frc_i%y = frc_i%y + tp2mdr(frc%y)
                 frc_i%z = frc_i%z + tp2mdr(frc%z)

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
           rstat = atomicAdd(em_buff(location), em_)

           ! Update forces
           rstat = atomicAdd( dem(1,i   ), frc_i%x )
           rstat = atomicAdd( dem(2,i   ), frc_i%y )
           rstat = atomicAdd( dem(3,i   ), frc_i%z )
           rstat = atomicAdd( dem(1,kbis), frc_k(threadIdx%x)%x )
           rstat = atomicAdd( dem(2,kbis), frc_k(threadIdx%x)%y )
           rstat = atomicAdd( dem(3,kbis), frc_k(threadIdx%x)%z )

           ! Update torque
           rstat = atomicAdd( tem(1,i)   ,ttmi%x )
           rstat = atomicAdd( tem(2,i)   ,ttmi%y )
           rstat = atomicAdd( tem(3,i)   ,ttmi%z )
           rstat = atomicAdd( tem(1,kbis),ttmk(threadIdx%x)%x )
           rstat = atomicAdd( tem(2,kbis),ttmk(threadIdx%x)%y )
           rstat = atomicAdd( tem(3,kbis),ttmk(threadIdx%x)%z )

           ! Update virial buffer
           rstat = atomicAdd(vir_buff(0*RED_BUFF_SIZE+location),vir_(1))
           rstat = atomicAdd(vir_buff(1*RED_BUFF_SIZE+location),vir_(2))
           rstat = atomicAdd(vir_buff(2*RED_BUFF_SIZE+location),vir_(3))
           rstat = atomicAdd(vir_buff(3*RED_BUFF_SIZE+location),vir_(4))
           rstat = atomicAdd(vir_buff(4*RED_BUFF_SIZE+location),vir_(5))
           rstat = atomicAdd(vir_buff(5*RED_BUFF_SIZE+location),vir_(6))

        end do

        end subroutine


        attributes(global) subroutine emreal3_cu
     &        ( ipole, pglob, loc, ieblst, eblst
     &        , npolelocnlb, npolelocnlb_pair, npolebloc, n
     &        , x, y, z, rpole
     &        , off2, f, alsq2, alsq2n, aewald
     &        , em_buff, nem_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &        )
        implicit none
        integer,value,intent(in)::npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
     &           ,off2,alsq2,alsq2n,aewald,f
        integer,device,intent(in)::ipole(*),pglob(*),loc(*)
     &                ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),rpole(13,*)
        integer  ,device:: nem_buff(*)
        ener_rtyp,device::  em_buff(*)
#ifdef TINKER_DEBUG
#endif

        integer ithread,iwarp,nwarp,ilane,klane,srclane
        integer ii,j,i,kbis
        integer iblock,idx,kdx
        integer iipole,iglob,iploc,kpole,kglob,kploc
#ifdef TINKER_DEBUG
        integer kglob_,kdx_
#endif
        integer location,istat
        integer nem_
        real(t_p) xk_,yk_,zk_,d2
        ener_rtyp em_
        real(t_p) rstat,zero
        type(real3) posi,pos
        type(real3),shared::posk(BLOCK_DIM)
        type(rpole_elt) ip,kp,kp_
        logical do_pair,same_block,accept_mid
        parameter(zero=0.0)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, npolelocnlb_pair, nwarp
           iblock = ieblst(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           i       = loc  (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           ip%c    = rpole( 1, iipole)
           ip%dx   = rpole( 2, iipole)
           ip%dy   = rpole( 3, iipole)
           ip%dz   = rpole( 4, iipole)
           ip%qxx  = rpole( 5, iipole)
           ip%qxy  = rpole( 6, iipole)
           ip%qxz  = rpole( 7, iipole)
           ip%qyy  = rpole( 9, iipole)
           ip%qyz  = rpole(10, iipole)
           ip%qzz  = rpole(13, iipole)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob   = ipole(kdx)
           kbis    = loc  (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kp%c    = rpole( 1, kpole)
           kp%dx   = rpole( 2, kpole)
           kp%dy   = rpole( 3, kpole)
           kp%dz   = rpole( 4, kpole)
           kp%qxx  = rpole( 5, kpole)
           kp%qxy  = rpole( 6, kpole)
           kp%qxz  = rpole( 7, kpole)
           kp%qyy  = rpole( 9, kpole)
           kp%qyz  = rpole(10, kpole)
           kp%qzz  = rpole(13, kpole)

           !* set compute Data to 0
           em_ = 0
           nem_ = 0

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_     = __shfl(kdx ,srclane)
              kglob_   = __shfl(kglob ,srclane)
#endif
              kp_%c    = __shfl(kp%c   ,klane)
              kp_%dx   = __shfl(kp%dx  ,klane)
              kp_%dy   = __shfl(kp%dy  ,klane)
              kp_%dz   = __shfl(kp%dz  ,klane)
              kp_%qxx  = __shfl(kp%qxx ,klane)
              kp_%qxy  = __shfl(kp%qxy ,klane)
              kp_%qxz  = __shfl(kp%qxz ,klane)
              kp_%qyy  = __shfl(kp%qyy ,klane)
              kp_%qyz  = __shfl(kp%qyz ,klane)
              kp_%qzz  = __shfl(kp%qzz ,klane)

              if (nproc.gt.1) then
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
                    pos%x=-pos%x; pos%y=-pos%y; pos%z=-pos%z;
                 end if
              else
                 pos%x = posk(klane)%x - posi%x
                 pos%y = posk(klane)%y - posi%y
                 pos%z = posk(klane)%z - posi%z
                 call image_inl(pos%x,pos%y,pos%z)
              end if

              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                        ,same_block)

              if (do_pair.and.d2<=off2.and.accept_mid) then
                 ! compute one interaction
                 call mpole3_couple(d2,pos%x,pos%y,pos%z,ip,kp_,zero
     &                             ,aewald,f,alsq2n,alsq2
     &                             ,em_,.false.)

                 nem_ = nem_ + 1
#ifdef TINKER_DEBUG
#endif
              end if
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

           ! Update energy & interaction number buffer
           rstat = atomicAdd( em_buff(location),  em_)
           istat = atomicAdd(nem_buff(location), nem_)
        end do

        end subroutine

        attributes(global) subroutine emrealshortlong3_cu
     &        ( ipole, pglob, loc, ieblst, eblst
     &        , npolelocnlb, npolelocnlb_pair, npolebloc, n
     &        , x, y, z, rpole
     &        , shortheal, r_cut, sh_cut2, off2, f
     &        , alsq2, alsq2n, aewald, mode
     &        , em_buff, nem_buff
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#ifdef TINKER_DEBUG
#endif
     &        )
        implicit none
        integer,value,intent(in)::npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair,mode
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,sh_cut2,shortheal
     &           ,r_cut,off2,alsq2,alsq2n,aewald,f
        integer,device,intent(in)::ipole(*),pglob(*),loc(*)
     &                ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),rpole(13,*)
        integer  ,device:: nem_buff(*)
        ener_rtyp,device::  em_buff(*)
#ifdef TINKER_DEBUG
#endif

        integer ithread,iwarp,nwarp,ilane,klane,srclane
        integer ii,j,i,kbis
        integer iblock,idx,kdx
        integer iipole,iglob,iploc,kpole,kglob,kploc
#ifdef TINKER_DEBUG
        integer kglob_,kdx_
#endif
        integer location,istat,nem_
        real(t_p) xk_,yk_,zk_,d2
        ener_rtyp em_
        real(t_p) rstat,zero
        type(real3) posi,pos
        type(real3),shared::posk(BLOCK_DIM)
        type(rpole_elt) ip,kp,kp_
        logical do_pair,same_block,accept_mid
        parameter(zero=0.0)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = ((ithread-1) / warpsize) + 1
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        accept_mid = .true.

        do ii = iwarp, npolelocnlb_pair, nwarp
           iblock = ieblst(ii)
           if (iblock==0) cycle
           idx     = (iblock-1)*WARP_SIZE + ilane;
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           i       = loc  (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           ip%c    = rpole( 1, iipole)
           ip%dx   = rpole( 2, iipole)
           ip%dy   = rpole( 3, iipole)
           ip%dz   = rpole( 4, iipole)
           ip%qxx  = rpole( 5, iipole)
           ip%qxy  = rpole( 6, iipole)
           ip%qxz  = rpole( 7, iipole)
           ip%qyy  = rpole( 9, iipole)
           ip%qyz  = rpole(10, iipole)
           ip%qzz  = rpole(13, iipole)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob   = ipole(kdx)
           kbis    = loc  (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kp%c    = rpole( 1, kpole)
           kp%dx   = rpole( 2, kpole)
           kp%dy   = rpole( 3, kpole)
           kp%dz   = rpole( 4, kpole)
           kp%qxx  = rpole( 5, kpole)
           kp%qxy  = rpole( 6, kpole)
           kp%qxz  = rpole( 7, kpole)
           kp%qyy  = rpole( 9, kpole)
           kp%qyz  = rpole(10, kpole)
           kp%qzz  = rpole(13, kpole)

           !* set compute Data to 0
           em_  = 0
           nem_ = 0

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
#ifdef TINKER_DEBUG
              kdx_     = __shfl(kdx ,srclane)
              kglob_   = __shfl(kglob ,srclane)
#endif
              kp_%c    = __shfl(kp%c   ,klane)
              kp_%dx   = __shfl(kp%dx  ,klane)
              kp_%dy   = __shfl(kp%dy  ,klane)
              kp_%dz   = __shfl(kp%dz  ,klane)
              kp_%qxx  = __shfl(kp%qxx ,klane)
              kp_%qxy  = __shfl(kp%qxy ,klane)
              kp_%qxz  = __shfl(kp%qxz ,klane)
              kp_%qyy  = __shfl(kp%qyy ,klane)
              kp_%qyz  = __shfl(kp%qyz ,klane)
              kp_%qzz  = __shfl(kp%qzz ,klane)

              if (nproc.gt.1) then
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
                    pos%x=-pos%x; pos%y=-pos%y; pos%z=-pos%z;
                 end if
              else
                 pos%x = posk(klane)%x - posi%x
                 pos%y = posk(klane)%y - posi%y
                 pos%z = posk(klane)%z - posi%z
                 call image_inl(pos%x,pos%y,pos%z)
              end if

              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                        ,same_block)

              if (do_pair.and.sh_cut2<=d2.and.d2<=off2.and.accept_mid)
     &           then
                 ! compute one interaction
                 call mpole3_couple_shortlong
     &                          (d2,pos%x,pos%y,pos%z,ip,kp_,zero
     &                          ,r_cut,shortheal,aewald,f,alsq2n,alsq2
     &                          ,em_,.false.,mode)

                 nem_ = nem_ + 1
#ifdef TINKER_DEBUG
#endif
              end if
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

           ! Update energy & interaction number buffer
           rstat = atomicAdd( em_buff(location),  em_)
           istat = atomicAdd(nem_buff(location), nem_)
        end do

        end subroutine

      end module
#endif
