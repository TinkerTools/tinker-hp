c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar1  --  polarization energy & derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar1cu" calculates the induced dipole polarization energy
c     and derivatives with respect to Cartesian coordinates in a CUDA Fortran kernel
c
c
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_macro.h"
#include "tinker_cudart.h"

      module epolar1cu
        use utilcu  ,only: nproc,ndir,ngrp,BLOCK_DIM,ALL_LANES
        use utilgpu ,only: BLOCK_SIZE,RED_BUFF_SIZE
     &              ,WARP_SIZE
        use tinTypes,only: real3,real6,mdyn3_r,rpole_elt

        contains

#include "image.f.inc"
#include "midpointimage.f.inc"
#include "convert.f.inc"
#include "pair_polar.inc.f"

        attributes(global) subroutine epreal1c_core_cu
     &        ( ipole, pglob, loc, ploc, ieblst, eblst
     &        , x, y, z, rpole, pdamp, dirdamp, thole, uind, uinp
     &        , pot, dep, trq, ep_buff, vir_buff
     &        , npolelocnlb, npolelocnlb_pair, npolebloc, n, u_cflx
     &        , u_ddamp, off2, f, alsq2, alsq2n, aewald
     &        ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend )
        implicit none
        integer  ,value,intent(in)::npolelocnlb,npolebloc,n
     &           ,npolelocnlb_pair
        logical  ,value,intent(in):: u_cflx, u_ddamp
        real(t_p),value:: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
     &           ,off2,alsq2,alsq2n,aewald,f
        integer  ,device,intent(in)::ipole(*),pglob(*),loc(*),ploc(*)
     &           ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),rpole(13,*)
     &           ,pdamp(*),dirdamp(*),thole(*),uind(3,*),uinp(3,*)
        real(t_p),device:: trq(3,*),vir_buff(*),pot(*)
        mdyn_rtyp,device:: dep(3,*)
        ener_rtyp,device:: ep_buff(*)

        integer ithread,iwarp,nwarp,ilane,klane,istat,srclane
        integer beg,ii,j,i,kbis
        integer iblock,idx,kdx,kdx_
        integer iipole,iglob,iploc,kpole,kploc
        integer location
        integer,shared::kglob(BLOCK_DIM)
        real(t_p) xk_,yk_,zk_,d2
        real(t_p) ep_,poti
        real(t_p) ipdp,ipgm,ipgm1,pdp,pgm,pgm1,rstat
        real(t_p), parameter :: eps = 1.e-5
        type(real6) dpui
        type(real3) posi,pos
        type(real3) frc
        type(mdyn3_r) frc_i
        type(mdyn3_r),shared::frc_k(BLOCK_DIM)
        type(real3)  ,shared:: posk(BLOCK_DIM)
        real(t_p)    ,shared:: kpdp(BLOCK_DIM),kpgm(BLOCK_DIM)
     &               ,kpgm1(BLOCK_DIM),potk(BLOCK_DIM)
        type(real3) trqi
        type(real3),shared:: trqk(BLOCK_DIM)
        type(real6),shared:: dpuk(BLOCK_DIM)
        real(t_p) vir_(6)
        type(rpole_elt),shared:: kp(BLOCK_DIM)
        type(rpole_elt) ip,kp_
        logical do_pair,same_block,accept_mid

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
           iploc   = ploc (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           ipdp    = pdamp(iipole)
           ipgm    = thole(iipole)
           ipgm1   = dirdamp(iipole)
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
           dpui%x  = uind ( 1, iipole)
           dpui%y  = uind ( 2, iipole)
           dpui%z  = uind ( 3, iipole)
           dpui%xx = uinp ( 1, iipole)
           dpui%yy = uinp ( 2, iipole)
           dpui%zz = uinp ( 3, iipole)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob(threadIdx%x)  = ipole(kdx)
           kbis    = loc  (kdx)
           kploc   = ploc (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kpdp(threadIdx%x)    = pdamp(kpole)
           kpgm(threadIdx%x)    = thole(kpole)
           kpgm1(threadIdx%x)   = dirdamp(kpole)
           kp(threadIdx%x)%c    = rpole( 1, kpole)
           kp(threadIdx%x)%dx   = rpole( 2, kpole)
           kp(threadIdx%x)%dy   = rpole( 3, kpole)
           kp(threadIdx%x)%dz   = rpole( 4, kpole)
           kp(threadIdx%x)%qxx  = rpole( 5, kpole)
           kp(threadIdx%x)%qxy  = rpole( 6, kpole)
           kp(threadIdx%x)%qxz  = rpole( 7, kpole)
           kp(threadIdx%x)%qyy  = rpole( 9, kpole)
           kp(threadIdx%x)%qyz  = rpole(10, kpole)
           kp(threadIdx%x)%qzz  = rpole(13, kpole)

           dpuk(threadIdx%x)%x  = uind ( 1, kpole)
           dpuk(threadIdx%x)%y  = uind ( 2, kpole)
           dpuk(threadIdx%x)%z  = uind ( 3, kpole)
           dpuk(threadIdx%x)%xx = uinp ( 1, kpole)
           dpuk(threadIdx%x)%yy = uinp ( 2, kpole)
           dpuk(threadIdx%x)%zz = uinp ( 3, kpole)

           call syncwarp(ALL_LANES)
           frc_i%x = 0;
           frc_i%y = 0;
           frc_i%z = 0;
           frc_k(threadIdx%x)%x = 0;
           frc_k(threadIdx%x)%y = 0;
           frc_k(threadIdx%x)%z = 0;
           trqi%x  = 0.0_ti_p;
           trqi%y  = 0.0_ti_p;
           trqi%z  = 0.0_ti_p;
           trqk(threadIdx%x)%x = 0.0_ti_p;
           trqk(threadIdx%x)%y = 0.0_ti_p;
           trqk(threadIdx%x)%z = 0.0_ti_p;

           if (u_cflx) then
              poti = 0.0
              potk(threadIdx%x)=0.0
           end if

           !* set compute Data to 0
           ep_ = 0
           vir_(1)=0.0; vir_(2)=0.0; vir_(3)=0.0;
           vir_(4)=0.0; vir_(5)=0.0; vir_(6)=0.0;

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
              kp_%c    = kp(klane)%c
              kp_%dx   = kp(klane)%dx
              kp_%dy   = kp(klane)%dy
              kp_%dz   = kp(klane)%dz
              kp_%qxx  = kp(klane)%qxx
              kp_%qxy  = kp(klane)%qxy
              kp_%qxz  = kp(klane)%qxz
              kp_%qyy  = kp(klane)%qyy
              kp_%qyz  = kp(klane)%qyz
              kp_%qzz  = kp(klane)%qzz
              pdp      = ipdp*kpdp(klane)
              pgm      = min(ipgm,kpgm(klane))
              if (u_ddamp) then
                 pgm1 = min(ipgm1,kpgm1(klane))
                 if (pgm1.ne.0.0) pgm1 = max(ipgm1,kpgm1(klane))
              end if

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
              do_pair = merge(.true.,iglob.lt.kglob(klane)
     &                       ,same_block)
              if (do_pair.and.d2<=off2.and.accept_mid
     &          .and.d2.gt.eps) then
                 ! compute one interaction
                 call epolar1_couple(dpui,ip,dpuk(klane),kp_,d2,pos
     &                       ,aewald,alsq2,alsq2n,pgm,pgm1,pdp,f
     &                       ,1.0_ti_p,1.0_ti_p,1.0_ti_p,u_ddamp
     &                       ,u_cflx,poti,potk(klane),ep_,frc
     &                       ,frc_k(klane),trqi,trqk(klane),.false.)

                 vir_(1) = vir_(1) + pos%x * frc%x
                 vir_(2) = vir_(2) + pos%y * frc%x
                 vir_(3) = vir_(3) + pos%z * frc%x
                 vir_(4) = vir_(4) + pos%y * frc%y
                 vir_(5) = vir_(5) + pos%z * frc%y
                 vir_(6) = vir_(6) + pos%z * frc%z

                 frc_i%x = frc_i%x - tp2mdr(frc%x)
                 frc_i%y = frc_i%y - tp2mdr(frc%y)
                 frc_i%z = frc_i%z - tp2mdr(frc%z)
              end if
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

           rstat = atomicAdd(ep_buff(location), tp2enr(ep_))

           rstat = atomicAdd( dep(1,i   ), frc_i%x )
           rstat = atomicAdd( dep(2,i   ), frc_i%y )
           rstat = atomicAdd( dep(3,i   ), frc_i%z )
           rstat = atomicAdd( trq(1,iploc),trqi%x )
           rstat = atomicAdd( trq(2,iploc),trqi%y )
           rstat = atomicAdd( trq(3,iploc),trqi%z )
           call syncwarp(ALL_LANES)
           rstat = atomicAdd( dep(1,kbis), frc_k(threadIdx%x)%x )
           rstat = atomicAdd( dep(2,kbis), frc_k(threadIdx%x)%y )
           rstat = atomicAdd( dep(3,kbis), frc_k(threadIdx%x)%z )
           rstat = atomicAdd( trq(1,kploc),trqk(threadIdx%x)%x )
           rstat = atomicAdd( trq(2,kploc),trqk(threadIdx%x)%y )
           rstat = atomicAdd( trq(3,kploc),trqk(threadIdx%x)%z )

           if (u_cflx) then
              rstat = atomicAdd( pot(i),poti )
              rstat = atomicAdd( pot(kbis),potk(threadIdx%x) )
           end if

           rstat = atomicAdd(vir_buff(0*RED_BUFF_SIZE+location),vir_(1))
           rstat = atomicAdd(vir_buff(1*RED_BUFF_SIZE+location),vir_(2))
           rstat = atomicAdd(vir_buff(2*RED_BUFF_SIZE+location),vir_(3))
           rstat = atomicAdd(vir_buff(3*RED_BUFF_SIZE+location),vir_(4))
           rstat = atomicAdd(vir_buff(4*RED_BUFF_SIZE+location),vir_(5))
           rstat = atomicAdd(vir_buff(5*RED_BUFF_SIZE+location),vir_(6))

        end do

        end subroutine

        attributes(global) subroutine mpreal1c_core_cu
     &        ( ipole,pglob,loc,ploc,ieblst,eblst
     &        , x,y,z,rpole,pdamp,thole,uind,uinp,use_group,grplist,wgrp
     &        , dep, trq, ep_buff, vir_buff
     &        , npolelocnlb, npolelocnlb_pair,npolebloc,n,mode
     &        , off2,f,alsq2,alsq2n,aewald,r_cut,s_cut2,shortheal
     &        , p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
        implicit none
        integer,value,intent(in)::npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair,mode
        logical,value,intent(in)::use_group
        real(t_p),value:: p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &           ,off2,alsq2,alsq2n,aewald,f,shortheal,r_cut,s_cut2
        integer,device,intent(in)::ipole(*),pglob(*),loc(*),ploc(*)
     &                ,ieblst(*),eblst(*),grplist(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),rpole(13,*)
     &                  ,pdamp(*),thole(*),uind(3,*),uinp(3,*)
     &                  ,wgrp(ngrp+1,ngrp+1)

        real(t_p),device:: trq(3,*),vir_buff(*)
        mdyn_rtyp,device:: dep(3,*)
        ener_rtyp,device:: ep_buff(*)

        integer ithread,iwarp,nwarp,ilane,klane,istat,srclane
        integer beg,ii,j,i,kbis
        integer iblock,idx,kdx,kdx_
        integer iipole,iglob,iploc,kpole,kploc
        integer location
        integer,shared::kglob(BLOCK_DIM)
        integer iga,igb_,igb
        !integer, shared :: igb(BLOCK_DIM)
        real(t_p) xk_,yk_,zk_,d2
        real(t_p) ep_,scale_
        real(t_p) ipdp,ipgm,pdp,pgm,rstat
        type(real6) dpui
        type(real3) posi,pos
        type(real3) frc
        type(mdyn3_r) frc_i
        type(mdyn3_r),shared::frc_k(BLOCK_DIM)
        type(real3)  ,shared:: posk(BLOCK_DIM)
        real(t_p)    ,shared:: kpdp(BLOCK_DIM),kpgm(BLOCK_DIM)
        type(real3) trqi,trqk_
        type(real3),shared:: trqk(BLOCK_DIM)
        type(real6),shared:: dpuk(BLOCK_DIM)
        real(t_p) vir_(6)
        type(rpole_elt),shared:: kp(BLOCK_DIM)
        type(rpole_elt) ip
        logical do_pair,same_block,accept_mid

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
           iploc   = ploc (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           ipdp    = pdamp(iipole)
           ipgm    = thole(iipole)
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
           dpui%x  = uind ( 1, iipole)
           dpui%y  = uind ( 2, iipole)
           dpui%z  = uind ( 3, iipole)
           dpui%xx = uinp ( 1, iipole)
           dpui%yy = uinp ( 2, iipole)
           dpui%zz = uinp ( 3, iipole)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob(threadIdx%x)  = ipole(kdx)
           kbis    = loc  (kdx)
           kploc   = ploc (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kpdp(threadIdx%x)    = pdamp(kpole)
           kpgm(threadIdx%x)    = thole(kpole)
           kp(threadIdx%x)%c    = rpole( 1, kpole)
           kp(threadIdx%x)%dx   = rpole( 2, kpole)
           kp(threadIdx%x)%dy   = rpole( 3, kpole)
           kp(threadIdx%x)%dz   = rpole( 4, kpole)
           kp(threadIdx%x)%qxx  = rpole( 5, kpole)
           kp(threadIdx%x)%qxy  = rpole( 6, kpole)
           kp(threadIdx%x)%qxz  = rpole( 7, kpole)
           kp(threadIdx%x)%qyy  = rpole( 9, kpole)
           kp(threadIdx%x)%qyz  = rpole(10, kpole)
           kp(threadIdx%x)%qzz  = rpole(13, kpole)

           dpuk(threadIdx%x)%x  = uind ( 1, kpole)
           dpuk(threadIdx%x)%y  = uind ( 2, kpole)
           dpuk(threadIdx%x)%z  = uind ( 3, kpole)
           dpuk(threadIdx%x)%xx = uinp ( 1, kpole)
           dpuk(threadIdx%x)%yy = uinp ( 2, kpole)
           dpuk(threadIdx%x)%zz = uinp ( 3, kpole)

           !if(use_group) then
           !   iga = grplist(iglob)
           !   !igb(threadIdx%x) = grplist(kglob(threadIdx%x))
           !   igb = grplist(kglob(threadIdx%x))
           !endif

           call syncwarp(ALL_LANES)
           frc_i%x = 0;
           frc_i%y = 0;
           frc_i%z = 0;
           frc_k(threadIdx%x)%x = 0;
           frc_k(threadIdx%x)%y = 0;
           frc_k(threadIdx%x)%z = 0;
           trqi%x  = 0.0_ti_p;
           trqi%y  = 0.0_ti_p;
           trqi%z  = 0.0_ti_p;
           trqk(threadIdx%x)%x = 0.0_ti_p;
           trqk(threadIdx%x)%y = 0.0_ti_p;
           trqk(threadIdx%x)%z = 0.0_ti_p;

           !* set compute Data to 0
           ep_ = 0
           vir_(1)=0.0; vir_(2)=0.0; vir_(3)=0.0;
           vir_(4)=0.0; vir_(5)=0.0; vir_(6)=0.0;

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
              pdp      = ipdp*kpdp(klane)
              pgm      = min(ipgm,kpgm(klane))
              !if(use_group) then
              !  igb_ = __shfl(igb,klane)
              !  !igb_ = igb(klane)
              !  scale_ = wgrp(iga+1,igb_+1)
              !else
                scale_  = 1.0_ti_p
              !endif

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
              do_pair = merge(.true.,iglob.lt.kglob(klane)
     &                       ,same_block)
              if (do_pair.and.s_cut2<=d2.and.d2<=off2.and.accept_mid)
     &           then
                 ! compute one interaction
                 call mpolar1_couple_comp
     &               (dpui,ip,dpuk(klane),kp(klane),d2,pos,
     &                aewald,alsq2,alsq2n,pgm,pdp,f,r_cut,shortheal,
     &                      1.0_ti_p-scale_,1.0_ti_p,1.0_ti_p,1.0_ti_p,
     &                ep_,frc,frc_k(klane),trqi,trqk(klane),.false.,
     &                mode,iglob,kglob(klane))

                 vir_(1) = vir_(1) - pos%x * frc%x
                 vir_(2) = vir_(2) - pos%y * frc%x
                 vir_(3) = vir_(3) - pos%z * frc%x
                 vir_(4) = vir_(4) - pos%y * frc%y
                 vir_(5) = vir_(5) - pos%z * frc%y
                 vir_(6) = vir_(6) - pos%z * frc%z

                 frc_i%x = frc_i%x + tp2mdr(frc%x)
                 frc_i%y = frc_i%y + tp2mdr(frc%y)
                 frc_i%z = frc_i%z + tp2mdr(frc%z)
              end if
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

           rstat = atomicAdd(ep_buff(location), tp2enr(ep_))

           rstat = atomicAdd(vir_buff(0*RED_BUFF_SIZE+location),vir_(1))
           rstat = atomicAdd(vir_buff(1*RED_BUFF_SIZE+location),vir_(2))
           rstat = atomicAdd(vir_buff(2*RED_BUFF_SIZE+location),vir_(3))
           rstat = atomicAdd(vir_buff(3*RED_BUFF_SIZE+location),vir_(4))
           rstat = atomicAdd(vir_buff(4*RED_BUFF_SIZE+location),vir_(5))
           rstat = atomicAdd(vir_buff(5*RED_BUFF_SIZE+location),vir_(6))

           rstat = atomicAdd( dep(1,i   ), frc_i%x )
           rstat = atomicAdd( dep(2,i   ), frc_i%y )
           rstat = atomicAdd( dep(3,i   ), frc_i%z )
           rstat = atomicAdd( trq(1,iploc),trqi%x )
           rstat = atomicAdd( trq(2,iploc),trqi%y )
           rstat = atomicAdd( trq(3,iploc),trqi%z )
           call syncwarp(ALL_LANES)
           rstat = atomicAdd( dep(1,kbis), frc_k(threadIdx%x)%x )
           rstat = atomicAdd( dep(2,kbis), frc_k(threadIdx%x)%y )
           rstat = atomicAdd( dep(3,kbis), frc_k(threadIdx%x)%z )
           rstat = atomicAdd( trq(1,kploc),trqk(threadIdx%x)%x )
           rstat = atomicAdd( trq(2,kploc),trqk(threadIdx%x)%y )
           rstat = atomicAdd( trq(3,kploc),trqk(threadIdx%x)%z )

        end do
        end subroutine

        attributes(global) subroutine epreal3_cu
     &        ( ipole, pglob, loc, ploc, ieblst, eblst
     &        , x, y, z, rpole, pdamp, thole, uind, uinp
     &        , ep_buff, nep_buff
     &        , npolelocnlb, npolelocnlb_pair, npolebloc, n
     &        , off2, f, alsq2, alsq2n, aewald, off, shortheal
     &        , use_short, use_dirdamp
     &        ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend )
        implicit none
        integer,value,intent(in)::npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair
        logical,value,intent(in)::use_short,use_dirdamp
        real(t_p),value:: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend
     &           ,off2,alsq2,alsq2n,aewald,f,shortheal,off
        integer,device,intent(in)::ipole(*),pglob(*),loc(*),ploc(*)
     &                ,ieblst(*),eblst(*)
        real(t_p),device,intent(in):: x(*),y(*),z(*),rpole(13,*)
     &                  ,pdamp(*),thole(*),uind(3,*),uinp(3,*)

        integer  ,device:: nep_buff(*)
        mdyn_rtyp,device:: ep_buff(*)

        integer ithread,iwarp,nwarp,ilane,klane,istat,srclane
        integer beg,ii,j,i,kbis
        integer iblock,idx,kdx,kdx_
        integer iipole,iglob,iploc,kpole,kploc
        integer location,nep_
        integer,shared::kglob(BLOCK_DIM)
        real(t_p) xk_,yk_,zk_,d2
        ener_rtyp ep_
        real(t_p) ipdp,ipgm,pdp,pgm,rstat
        type(real6) dpui
        type(real3) posi,pos
        type(real3),shared:: posk(BLOCK_DIM)
        real(t_p)  ,shared:: kpdp(BLOCK_DIM),kpgm(BLOCK_DIM)
        type(real6),shared:: dpuk(BLOCK_DIM)
        type(rpole_elt) ip,kp,kp_
        logical do_pair,same_block,accept_mid

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
           iploc   = ploc (idx)
           posi%x  =     x(idx)
           posi%y  =     y(idx)
           posi%z  =     z(idx)
           ipdp    = pdamp(iipole)
           ipgm    = thole(iipole)
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
           dpui%x  = uind ( 1, iipole)
           dpui%y  = uind ( 2, iipole)
           dpui%z  = uind ( 3, iipole)
           dpui%xx = uinp ( 1, iipole)
           dpui%yy = uinp ( 2, iipole)
           dpui%zz = uinp ( 3, iipole)

           !  Load atom block k parameters
           kdx     = eblst( (ii-1)*WARP_SIZE+ ilane )
           kpole   = pglob(kdx)
           kglob(threadIdx%x)   = ipole(kdx)
           kbis    = loc  (kdx)
           kploc   = ploc (kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kpdp(threadIdx%x)    = pdamp(kpole)
           kpgm(threadIdx%x)    = thole(kpole)
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

           dpuk(threadIdx%x)%x  = uind ( 1, kpole)
           dpuk(threadIdx%x)%y  = uind ( 2, kpole)
           dpuk(threadIdx%x)%z  = uind ( 3, kpole)
           dpuk(threadIdx%x)%xx = uinp ( 1, kpole)
           dpuk(threadIdx%x)%yy = uinp ( 2, kpole)
           dpuk(threadIdx%x)%zz = uinp ( 3, kpole)

           !* set compute Data to 0
           ep_ = 0
           nep_= 0

           same_block = (idx.ne.kdx)
 
           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
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
              pdp      = ipdp*kpdp(klane)
              pgm      = min(ipgm,kpgm(klane))
              if (pgm.eq.0.0) pgm = max(ipgm,kpgm(klane))

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
                    pos%x=-pos%x; pos%y=-pos%y; pos%z=-pos%z;
                 end if
              else
                 pos%x = posk(klane)%x - posi%x
                 pos%y = posk(klane)%y - posi%y
                 pos%z = posk(klane)%z - posi%z
                 call image_inl(pos%x,pos%y,pos%z)
              end if
              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = merge(.true.,iglob.lt.kglob(klane)
     &                       ,same_block)
              if (do_pair.and.d2<=off2.and.accept_mid) then
                 ! compute one interaction
                 call epolar3_couple(dpui,ip,dpuk(klane),kp_,d2,pos
     &                       ,aewald,alsq2,alsq2n,pgm,pdp,use_dirdamp
     &                       ,f,off,shortheal,1.0_ti_p
     &                       ,ep_,use_short,.false.)
                 nep_ = nep_ + 1

              end if
           end do

           location = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1

           rstat = atomicAdd( ep_buff(location),  ep_)
           istat = atomicAdd(nep_buff(location), nep_)
        end do

        end subroutine

        ! Loop on scale interaction for correction
        attributes(global) subroutine epreal1c_scaling_cu
     &             (dpucorrect_ik,dpucorrect_scale,ipole,loc
     &             ,polelocnl,x,y,z,pdamp,dirdamp,thole,rpole,uind,uinp
     &             ,pot,dep,trq,ep_buff,vir_buff,n,n_dpuscale
     &             ,nbloc,u_cflx,u_ddamp,off2,aewald,alsq2,alsq2n,f)
        implicit none

        integer  ,device,intent(in)::dpucorrect_ik(*),ipole(n),loc(n)
     &           ,polelocnl(n)
        real(t_p),device,intent(in)::dpucorrect_scale(*),x(n),y(n),z(n)
     &           ,pdamp(n),thole(n),dirdamp(*)
     &           ,rpole(13,n),uind(3,n),uinp(3,n)
        integer  ,value,intent(in)::n,n_dpuscale,nbloc
        logical  ,value,intent(in)::u_cflx,u_ddamp
        real(t_p),value,intent(in)::aewald,f,alsq2,alsq2n,off2
        real(t_p),device,intent(inout)::trq(3,*),pot(*)
     &           ,vir_buff(6*RED_BUFF_SIZE)
        mdyn_rtyp,device,intent(inout)::dep(3,nbloc)
        ener_rtyp,device,intent(inout)::ep_buff(RED_BUFF_SIZE)

        integer   i,k,iglob,kglob,iploc,kploc,ithread
        integer   ii,iipole,kpole,j,lot
        real(t_p) r2,pgamma,pgama1,damp,rstat
        real(t_p) e
        real(t_p) pdi,pti,ddi
        real(t_p) pscale,dscale,uscale,poti,potk
        type(rpole_elt) ip,kp
        type(real6)     dpui,dpuk
        type(real3)     pos,posi,trqi,trqk,frc
        type(mdyn3_r)   frc_r
        integer iga,igb

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        do ii = ithread, n_dpuscale, blockDim%x*gridDim%x
           iipole   = dpucorrect_ik(2*(ii-1)+1)
           kpole    = dpucorrect_ik(2*(ii-1)+2)

           iglob    = ipole(iipole)
           kglob    = ipole(kpole)
           i        = loc  (iglob)
           k        = loc  (kglob)
           iploc    = polelocnl(iipole)
           kploc    = polelocnl(kpole)

           pos%x    = x(kglob) - x(iglob)
           pos%y    = y(kglob) - y(iglob)
           pos%z    = z(kglob) - z(iglob)

           call image_inl(pos%x,pos%y,pos%z)
           ! cutoff
           r2       = pos%x**2 + pos%y**2 + pos%z**2
           if (r2>off2) cycle

           pdi      = pdamp(iipole)
           pti      = thole(iipole)
           ddi      = dirdamp(iipole)

           ip%c     = rpole( 1,iipole)
           ip%dx    = rpole( 2,iipole)
           ip%dy    = rpole( 3,iipole)
           ip%dz    = rpole( 4,iipole)
           ip%qxx   = rpole( 5,iipole)
           ip%qxy   = rpole( 6,iipole)
           ip%qxz   = rpole( 7,iipole)
           ip%qyy   = rpole( 9,iipole)
           ip%qyz   = rpole(10,iipole)
           ip%qzz   = rpole(13,iipole)

           dpui%x   = uind ( 1,iipole)
           dpui%y   = uind ( 2,iipole)
           dpui%z   = uind ( 3,iipole)
           dpui%xx  = uinp ( 1,iipole)
           dpui%yy  = uinp ( 2,iipole)
           dpui%zz  = uinp ( 3,iipole)

           dscale   = dpucorrect_scale(3*(ii-1)+1)
           pscale   = dpucorrect_scale(3*(ii-1)+2)
           uscale   = dpucorrect_scale(3*(ii-1)+3)

           kp%c     = rpole( 1, kpole)
           kp%dx    = rpole( 2, kpole)
           kp%dy    = rpole( 3, kpole)
           kp%dz    = rpole( 4, kpole)
           kp%qxx   = rpole( 5, kpole)
           kp%qxy   = rpole( 6, kpole)
           kp%qxz   = rpole( 7, kpole)
           kp%qyy   = rpole( 9, kpole)
           kp%qyz   = rpole(10, kpole)
           kp%qzz   = rpole(13, kpole)

           dpuk%x   = uind ( 1, kpole)
           dpuk%y   = uind ( 2, kpole)
           dpuk%z   = uind ( 3, kpole)
           dpuk%xx  = uinp ( 1, kpole)
           dpuk%yy  = uinp ( 2, kpole)
           dpuk%zz  = uinp ( 3, kpole)

           pgamma   = min( pti,thole(kpole) )
           damp     = pdi*pdamp(kpole)
           if (u_ddamp) then
              pgama1   = min( ddi,dirdamp(kpole) )
              if (pgama1.eq.0.0) pgama1   = max( ddi,dirdamp(kpole) )
           end if

           ! zero outputs
           e=0;
           frc_r%x=0; frc_r%y=0; frc_r%z=0;
            trqk%x=0;  trqk%y=0;  trqk%z=0;
            trqi%x=0;  trqi%y=0;  trqi%z=0;
             poti=0.0; potk=0.0

           ! Compute polar interaction
           call epolar1_couple(dpui,ip,dpuk,kp,r2,pos
     &                 ,aewald,alsq2,alsq2n,pgamma,pgama1,damp,f
     &                 ,dscale,pscale,uscale,u_ddamp
     &                 ,u_cflx,poti,potk
     &                 ,e,frc,frc_r,trqi,trqk,.true.)

           ! Increment energy
           lot   = iand( ithread-1,RED_BUFF_SIZE-1 ) + 1
           rstat = atomicAdd( ep_buff(lot), tp2enr(e) )
           ! Increment gradient and virial due to Cartesian forces
           rstat = atomicAdd( dep(1,i),-frc_r%x )
           rstat = atomicAdd( dep(2,i),-frc_r%y )
           rstat = atomicAdd( dep(3,i),-frc_r%z )
           rstat = atomicAdd( dep(1,k), frc_r%x )
           rstat = atomicAdd( dep(2,k), frc_r%y )
           rstat = atomicAdd( dep(3,k), frc_r%z )

           if (u_cflx) then
              rstat = atomicAdd( pot(i),poti )
              rstat = atomicAdd( pot(k),potk )
           end if

           rstat = atomicAdd(vir_buff(0*RED_BUFF_SIZE+lot),pos%x*frc%x)
           rstat = atomicAdd(vir_buff(1*RED_BUFF_SIZE+lot),pos%y*frc%x)
           rstat = atomicAdd(vir_buff(2*RED_BUFF_SIZE+lot),pos%z*frc%x)
           rstat = atomicAdd(vir_buff(3*RED_BUFF_SIZE+lot),pos%y*frc%y)
           rstat = atomicAdd(vir_buff(4*RED_BUFF_SIZE+lot),pos%z*frc%y)
           rstat = atomicAdd(vir_buff(5*RED_BUFF_SIZE+lot),pos%z*frc%z)
           ! Increment torque
           rstat = atomicAdd( trq(1,iploc),trqi%x )
           rstat = atomicAdd( trq(2,iploc),trqi%y )
           rstat = atomicAdd( trq(3,iploc),trqi%z )
           rstat = atomicAdd( trq(1,kploc),trqk%x )
           rstat = atomicAdd( trq(2,kploc),trqk%y )
           rstat = atomicAdd( trq(3,kploc),trqk%z )
        end do
        end subroutine

      end module
#endif
