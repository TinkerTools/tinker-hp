c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_precision.h"
#include "tinker_cudart.h"
#include "tinker_types.h"

      module tmatxb_pmecu
        use utilcu  ,only: nproc,BLOCK_DIM,ALL_LANES
     &              ,bal_comp=>balanced_comput
        use utilgpu ,only: BLOCK_SIZE
        use tintypes,only: rpole_elt,real3,real6,real7
        use neigh   ,only: i_bl=>inside_block,d_bl=>disc_block

        contains

#include "image.f.inc"
#include "midpointimage.f.inc"
#include "pair_tmatxb.f.inc"
#define __tfea__ (__use_mpi__)

        attributes(global) subroutine tmatxb_pme_core_cu
     &                    (ipole,pglob,b_stat,ploc,ieblst,eblst,x,y,z
     &                    ,pdamp,thole,polarity,mu,efi
     &              ,npolelocnl,npolelocnlb,npolelocnlb_pair,npolebloc,n
     &                    ,cut2,alsq2,alsq2n,aewald
     &                    ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
        implicit none
        integer,value,intent(in):: npolelocnl,npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,cut2,alsq2,alsq2n,aewald
        integer,device,intent(in)::ipole(npolelocnlb),pglob(npolelocnlb)
     &         ,ploc(npolelocnlb),ieblst(npolelocnlb_pair)
     &         ,eblst(npolelocnlb_pair*(BLOCK_SIZE+2)),b_stat(*)
        real(t_p),device,intent(in)::pdamp(n),thole(n),polarity(n)
     &           ,x(npolelocnlb),y(npolelocnlb),z(npolelocnlb)
        real(t_p),device,intent(in):: mu(3,2,npolebloc)
        real(t_p),device,intent(inout):: efi(npolebloc,3,2)

        integer ithread,iwarp,nwarp,ilane,klane,nlane,istat,srclane
        integer beg,ii,j
        integer iblock,idx,kdx,kdx_
        integer iipole,iglob,kglob,iploc,kpole,kploc
        real(t_p) xk_,yk_,zk_,d2
        real(t_p) ipdp,ipgm,pdp,pgm,rstat
        type(real6) dpui
        type(real3),shared:: posi(BLOCK_DIM)
        type(real3) pos,posk
        type(real3) fid,fip
        type(real3) fkd,fkp
        real(t_p)  ,shared:: kpdp(BLOCK_DIM),kpgm(BLOCK_DIM)
        type(real7),shared:: dpuk(BLOCK_DIM)
        logical do_pair

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        nlane   = merge(1,ilane+1,ilane.eq.warpsize)
        !beg     = 2*npolelocnlb_pair

c       if (ithread.eq.1) print*,npolelocnlb,npolebloc,n,cut2,alsq2
c    &     ,aewald,npolelocnlb_pair

        do ii = iwarp, npolelocnlb_pair-1, nwarp

           ! Load atom block k parameters
           kdx     = eblst(ii*warpsize + ilane)
           kpole   = pglob(kdx)
           kglob   = ipole(kdx)
           kploc   = ploc(kdx)
           posk%x  = x(kdx)
           posk%y  = y(kdx)
           posk%z  = z(kdx)
           kpdp(threadIdx%x)    = pdamp(kpole)
           kpgm(threadIdx%x)    = thole(kpole)
           dpuk(threadIdx%x)%x  = mu(1,1,kploc)
           dpuk(threadIdx%x)%y  = mu(2,1,kploc)
           dpuk(threadIdx%x)%z  = mu(3,1,kploc)
           dpuk(threadIdx%x)%xx = mu(1,2,kploc)
           dpuk(threadIdx%x)%yy = mu(2,2,kploc)
           dpuk(threadIdx%x)%zz = mu(3,2,kploc)

           ! Load atom block i parameters
           iblock  = ieblst(ii+1)
           if (iblock.eq.0) cycle
           idx     = (iblock-1)*warpsize + ilane
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           iploc   = ploc(idx)
           posi(threadIdx%x)%x  = x(idx)
           posi(threadIdx%x)%y  = y(idx)
           posi(threadIdx%x)%z  = z(idx)
           ipdp    = pdamp(iipole)
           ipgm    = thole(iipole)
           dpui%x  = mu(1,1,iploc)
           dpui%y  = mu(2,1,iploc)
           dpui%z  = mu(3,1,iploc)
           dpui%xx = mu(1,2,iploc)
           dpui%yy = mu(2,2,iploc)
           dpui%zz = mu(3,2,iploc)

           !set compute Data to 0
           fid%x   = 0
           fid%y   = 0
           fid%z   = 0
           fip%x   = 0
           fip%y   = 0
           fip%z   = 0
           fkd%x = 0
           fkd%y = 0
           fkd%z = 0
           fkp%x = 0
           fkp%y = 0
           fkp%z = 0

           do j = 0,warpsize-1
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
              !kdx_     = __shfl(kdx ,srclane)
              pdp      = ipdp*kpdp(klane)
              pgm      = min(ipgm,kpgm(klane))

              do_pair = merge(iglob.lt.kglob,.true.,idx.eq.kdx)
#if __tfea__ & __use_mpi__
              if (nproc.gt.1.and.bal_comp) then
                 block; real(t_p) xk_,yk_,zk_
                 xk_   = posk%x
                 yk_   = posk%y
                 zk_   = posk%z
                 pos%x = posi(threadIdx%x)%x - xk_
                 pos%y = posi(threadIdx%x)%y - yk_
                 pos%z = posi(threadIdx%x)%z - zk_
                 call midpointimage_inl(xk_,yk_,zk_,pos%x,pos%y,pos%z)
                 if ((zk_.lt.p_zbeg).or.(zk_.ge.p_zend)
     &           .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &           .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                    do_pair = .false.
                 end if
                 end block
              else
                 pos%x = posi(threadIdx%x)%x - posk%x
                 pos%y = posi(threadIdx%x)%y - posk%y
                 pos%z = posi(threadIdx%x)%z - posk%z
                 call image_inl(pos%x,pos%y,pos%z)
              end if
#else
              pos%x = posi(threadIdx%x)%x - posk%x
              pos%y = posi(threadIdx%x)%y - posk%y
              pos%z = posi(threadIdx%x)%z - posk%z
              call image_inl(pos%x,pos%y,pos%z)
#endif
              d2      = pos%x**2 + pos%y**2 + pos%z**2

              do_pair = do_pair.and.d2<=cut2
              if (do_pair) then
                  ! compute one interaction
                  call duo_tmatxb(d2,pos,dpui,dpuk(klane)
     &                    ,pdp,pgm,aewald,alsq2,alsq2n,real(1,t_p)
     &                    ,fid,fip,fkd,fkp,.false.)
              end if

              kglob  = __shfl(  kglob,nlane )
              posk%x = __shfl( posk%x,nlane )
              posk%y = __shfl( posk%y,nlane )
              posk%z = __shfl( posk%z,nlane )
               fkd%x = __shfl(  fkd%x,nlane )
               fkd%y = __shfl(  fkd%y,nlane )
               fkd%z = __shfl(  fkd%z,nlane )
               fkp%x = __shfl(  fkp%x,nlane )
               fkp%y = __shfl(  fkp%y,nlane )
               fkp%z = __shfl(  fkp%z,nlane )
 
           end do

           ! increment electric field for each atoms
           if (idx.le.npolelocnl) then
           rstat = atomicAdd( efi(idx,1,1),fid%x )
           rstat = atomicAdd( efi(idx,2,1),fid%y )
           rstat = atomicAdd( efi(idx,3,1),fid%z )
           rstat = atomicAdd( efi(idx,1,2),fip%x )
           rstat = atomicAdd( efi(idx,2,2),fip%y )
           rstat = atomicAdd( efi(idx,3,2),fip%z )
           end if
           if (kdx.le.npolelocnl) then
           rstat = atomicAdd( efi(kdx,1,1),fkd%x )
           rstat = atomicAdd( efi(kdx,2,1),fkd%y )
           rstat = atomicAdd( efi(kdx,3,1),fkd%z )
           rstat = atomicAdd( efi(kdx,1,2),fkp%x )
           rstat = atomicAdd( efi(kdx,2,2),fkp%y )
           rstat = atomicAdd( efi(kdx,3,2),fkp%z )
           end if
        end do
        end subroutine

        attributes(global) subroutine tmatxb_pme_core_v4_cu
     &                    (ipole,pglob,ploc,ieblst,eblst,x,y,z
     &                    ,pdamp,thole,polarity
     &                    ,mu,mu3,mu4,efi,efi3,efi4
     &                    ,npolelocnlb,npolelocnlb_pair,npolebloc,n
     &                    ,cut2,alsq2,alsq2n,aewald
     &                    ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
        implicit none
        integer,value,intent(in)::npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,cut2,alsq2,alsq2n,aewald
        integer,device,intent(in)::ipole(npolelocnlb),pglob(npolelocnlb)
     &         ,ploc(npolelocnlb),ieblst(npolelocnlb_pair)
     &         ,eblst(npolelocnlb_pair*(BLOCK_SIZE+2))
        real(t_p),device,intent(in)::pdamp(n),thole(n),polarity(n)
     &           ,x(npolelocnlb),y(npolelocnlb),z(npolelocnlb)
        real(t_p),device,intent(in):: mu(3,2,npolebloc),mu3(3,npolebloc)
     &           ,mu4(3,npolebloc)
        real(t_p),device,intent(inout):: efi(3,2,npolebloc)
     &           ,efi3(3,npolebloc),efi4(3,npolebloc)

        integer ithread,iwarp,nwarp,ilane,klane,istat,srclane
        integer beg,ii,j
        integer iblock,idx,kdx,kdx_
        integer iipole,iglob,iploc,kpole,kploc
        integer,shared:: kglob(BLOCK_DIM)
        real(t_p) xk_,yk_,zk_,d2
        real(t_p) ipdp,ipgm,pdp,pgm,rstat
        type(real6) dpui,equi,dpuk_
        type(real3) posi,pos
        type(real3) fid,fip,fie,fiq
        type(real3),shared::fkd(BLOCK_DIM),fkp(BLOCK_DIM)
     &             ,fke(BLOCK_DIM),fkq(BLOCK_DIM)
     &             ,posk(BLOCK_DIM)
        real(t_p)  ,shared:: kpdp(BLOCK_DIM),kpgm(BLOCK_DIM)
        type(real6),shared:: dpuk(BLOCK_DIM),equk(BLOCK_DIM)
        logical do_pair,same_block,accept_mid

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        !beg     = 2*npolelocnlb_pair
        accept_mid = .true.

c       if (ithread.eq.1) print*,npolelocnlb,npolebloc,n,cut2,alsq2
c    &     ,aewald,npolelocnlb_pair

        do ii = iwarp, npolelocnlb_pair-1, nwarp

           ! Load atom block k parameters
           kdx     = eblst(ii*warpsize + ilane)
           kpole   = pglob(kdx)
           kglob(threadIdx%x)   = ipole(kdx)
           kploc   = ploc(kdx)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kpdp(threadIdx%x)    = pdamp(kpole)
           kpgm(threadIdx%x)    = thole(kpole)
           dpuk(threadIdx%x)%x  = mu(1,1,kploc)
           dpuk(threadIdx%x)%y  = mu(2,1,kploc)
           dpuk(threadIdx%x)%z  = mu(3,1,kploc)
           dpuk(threadIdx%x)%xx = mu(1,2,kploc)
           dpuk(threadIdx%x)%yy = mu(2,2,kploc)
           dpuk(threadIdx%x)%zz = mu(3,2,kploc)
           equk(threadIdx%x)%x  = mu3(1,kploc)
           equk(threadIdx%x)%y  = mu3(2,kploc)
           equk(threadIdx%x)%z  = mu3(3,kploc)
           equk(threadIdx%x)%xx = mu4(1,kploc)
           equk(threadIdx%x)%yy = mu4(2,kploc)
           equk(threadIdx%x)%zz = mu4(3,kploc)
           call syncwarp(ALL_LANES)

           ! Load atom block i parameters
           iblock  = ieblst(ii+1)
           if (iblock.eq.0) cycle
           idx     = (iblock-1)*warpsize + ilane
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           iploc   = ploc(idx)
           posi%x  = x(idx)
           posi%y  = y(idx)
           posi%z  = z(idx)
           ipdp    = pdamp(iipole)
           ipgm    = thole(iipole)
           dpui%x  = mu(1,1,iploc)
           dpui%y  = mu(2,1,iploc)
           dpui%z  = mu(3,1,iploc)
           dpui%xx = mu(1,2,iploc)
           dpui%yy = mu(2,2,iploc)
           dpui%zz = mu(3,2,iploc)
           equi%x  = mu3(1,iploc)
           equi%y  = mu3(2,iploc)
           equi%z  = mu3(3,iploc)
           equi%xx = mu4(1,iploc)
           equi%yy = mu4(2,iploc)
           equi%zz = mu4(3,iploc)

           !set compute Data to 0
           fid%x   = 0
           fid%y   = 0
           fid%z   = 0
           fip%x   = 0
           fip%y   = 0
           fip%z   = 0
           fie%x   = 0
           fie%y   = 0
           fie%z   = 0
           fiq%x   = 0
           fiq%y   = 0
           fiq%z   = 0
           fkd(threadIdx%x)%x = 0
           fkd(threadIdx%x)%y = 0
           fkd(threadIdx%x)%z = 0
           fkp(threadIdx%x)%x = 0
           fkp(threadIdx%x)%y = 0
           fkp(threadIdx%x)%z = 0
           fke(threadIdx%x)%x = 0
           fke(threadIdx%x)%y = 0
           fke(threadIdx%x)%z = 0
           fkq(threadIdx%x)%x = 0
           fkq(threadIdx%x)%y = 0
           fkq(threadIdx%x)%z = 0

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
              !kdx_     = __shfl(kdx ,srclane)
              pdp      = ipdp*kpdp(klane)
              pgm      = min(ipgm,kpgm(klane))

              if (nproc.gt.1.and.bal_comp) then
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
              do_pair = merge(.true.,iglob.lt.kglob(klane)
     &                       ,same_block)
              if (do_pair.and.d2<=cut2.and.accept_mid) then
                  ! compute one interaction
                  call
     &           tmatxb4_couple(d2,pos,dpui,dpuk(klane),equi,equk(klane)
     &                    ,pdp,pgm,aewald,alsq2,alsq2n,real(1,t_p)
     &                    ,fid,fip,fkd(klane),fkp(klane)
     &                    ,fie,fiq,fke(klane),fkq(klane),.false.)
              end if
 
           end do

           ! increment electric field for each atoms
           rstat = atomicAdd( efi(1,1,iploc),fid%x )
           rstat = atomicAdd( efi(2,1,iploc),fid%y )
           rstat = atomicAdd( efi(3,1,iploc),fid%z )
           rstat = atomicAdd( efi(1,2,iploc),fip%x )
           rstat = atomicAdd( efi(2,2,iploc),fip%y )
           rstat = atomicAdd( efi(3,2,iploc),fip%z )
           rstat = atomicAdd( efi3(1,iploc),fie%x )
           rstat = atomicAdd( efi3(2,iploc),fie%y )
           rstat = atomicAdd( efi3(3,iploc),fie%z )
           rstat = atomicAdd( efi4(1,iploc),fiq%x )
           rstat = atomicAdd( efi4(2,iploc),fiq%y )
           rstat = atomicAdd( efi4(3,iploc),fiq%z )
           call syncwarp(ALL_LANES)
           rstat = atomicAdd( efi(1,1,kploc),fkd(threadIdx%x)%x )
           rstat = atomicAdd( efi(2,1,kploc),fkd(threadIdx%x)%y )
           rstat = atomicAdd( efi(3,1,kploc),fkd(threadIdx%x)%z )
           rstat = atomicAdd( efi(1,2,kploc),fkp(threadIdx%x)%x )
           rstat = atomicAdd( efi(2,2,kploc),fkp(threadIdx%x)%y )
           rstat = atomicAdd( efi(3,2,kploc),fkp(threadIdx%x)%z )
           rstat = atomicAdd( efi3(1,kploc),fke(threadIdx%x)%x )
           rstat = atomicAdd( efi3(2,kploc),fke(threadIdx%x)%y )
           rstat = atomicAdd( efi3(3,kploc),fke(threadIdx%x)%z )
           rstat = atomicAdd( efi4(1,kploc),fkq(threadIdx%x)%x )
           rstat = atomicAdd( efi4(2,kploc),fkq(threadIdx%x)%y )
           rstat = atomicAdd( efi4(3,kploc),fkq(threadIdx%x)%z )
        end do
        end subroutine


        attributes(global) subroutine otfdc_tmatxb_pme_core_cu
     &                    (ipole,pglob,ploc,grplst,ieblst,eblst,x,y,z
     &                    ,pdamp,thole,polarity,mu,efi
     &                    ,npolelocnlb,npolelocnlb_pair,npolebloc,n
     &                    ,cut2,alsq2,alsq2n,aewald
     &                    ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
        implicit none
        integer,value,intent(in)::npolelocnlb,npolebloc,n
     &         ,npolelocnlb_pair
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,cut2,alsq2,alsq2n,aewald
        integer,device,intent(in)::ipole(npolelocnlb),pglob(npolelocnlb)
     &         ,ploc(npolelocnlb),ieblst(npolelocnlb_pair)
     &         ,grplst(n)
     &         ,eblst(npolelocnlb_pair*(BLOCK_SIZE+2))
        real(t_p),device,intent(in)::pdamp(n),thole(n),polarity(n)
     &           ,x(npolelocnlb),y(npolelocnlb),z(npolelocnlb)
        real(t_p),device,intent(in):: mu(3,2,npolebloc)
        real(t_p),device,intent(inout):: efi(3,2,npolebloc)

        integer ithread,iwarp,nwarp,ilane,klane,istat,srclane
        integer beg,ii,j,li
        integer iblock,idx,kdx,kdx_
        integer iipole,iglob,iploc,kpole,kglob,kploc
        real(t_p) xk_,yk_,zk_,d2
        real(t_p) ipdp,ipgm,pdp,pgm,rstat
        type(real6) dpui,dpuk_
        type(real3) posi,pos
        type(real3) fid,fip
        type(real3),shared::fkd(BLOCK_DIM),fkp(BLOCK_DIM)
     &             ,posk(BLOCK_DIM)
        integer    ,shared:: lk(BLOCK_DIM)
        real(t_p)  ,shared:: kpdp(BLOCK_DIM),kpgm(BLOCK_DIM)
        type(real6),shared:: dpuk(BLOCK_DIM)
        logical do_pair,same_block,accept_mid

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        !beg     = 2*npolelocnlb_pair
        accept_mid = .true.

c       if (ithread.eq.1) print*,npolelocnlb,npolebloc,n,cut2,alsq2
c    &     ,aewald,npolelocnlb_pair

        do ii = iwarp, npolelocnlb_pair-1, nwarp
           ! Load atom block i parameters
           iblock  = ieblst(ii+1)
           if (iblock.eq.0) cycle
           idx     = (iblock-1)*warpsize + ilane
           iipole  = pglob(idx)
           iglob   = ipole(idx)
           iploc   = ploc(idx)
           posi%x  = x(idx)
           posi%y  = y(idx)
           posi%z  = z(idx)
           ipdp    = pdamp(iipole)
           ipgm    = thole(iipole)
           li      = grplst(iglob)
           dpui%x  = mu(1,1,iploc)
           dpui%y  = mu(2,1,iploc)
           dpui%z  = mu(3,1,iploc)
           dpui%xx = mu(1,2,iploc)
           dpui%yy = mu(2,2,iploc)
           dpui%zz = mu(3,2,iploc)
           call syncwarp(ALL_LANES)

           ! Load atom block k parameters
           kdx     = eblst(ii*warpsize + ilane)
           kpole   = pglob(kdx)
           kglob   = ipole(kdx)
           kploc   = ploc(kdx)
           lk  (threadIdx%x)    = grplst(kglob)
           posk(threadIdx%x)%x  = x(kdx)
           posk(threadIdx%x)%y  = y(kdx)
           posk(threadIdx%x)%z  = z(kdx)
           kpdp(threadIdx%x)    = pdamp(kpole)
           kpgm(threadIdx%x)    = thole(kpole)
           dpuk(threadIdx%x)%x  = mu(1,1,kploc)
           dpuk(threadIdx%x)%y  = mu(2,1,kploc)
           dpuk(threadIdx%x)%z  = mu(3,1,kploc)
           dpuk(threadIdx%x)%xx = mu(1,2,kploc)
           dpuk(threadIdx%x)%yy = mu(2,2,kploc)
           dpuk(threadIdx%x)%zz = mu(3,2,kploc)

           !set compute Data to 0
           fid%x   = 0
           fid%y   = 0
           fid%z   = 0
           fip%x   = 0
           fip%y   = 0
           fip%z   = 0
           fkd(threadIdx%x)%x = 0
           fkd(threadIdx%x)%y = 0
           fkd(threadIdx%x)%z = 0
           fkp(threadIdx%x)%x = 0
           fkp(threadIdx%x)%y = 0
           fkp(threadIdx%x)%z = 0

           same_block = (idx.ne.kdx)

           do j = 1,warpsize
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
              !kdx_     = __shfl(kdx ,srclane)
              pdp      = ipdp*kpdp(klane)
              pgm      = min(ipgm,kpgm(klane))

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
                 end if
              else
                 pos%x = posi%x - posk(klane)%x
                 pos%y = posi%y - posk(klane)%y
                 pos%z = posi%z - posk(klane)%z
                 call image_inl(pos%x,pos%y,pos%z)
              end if
              d2      = pos%x**2 + pos%y**2 + pos%z**2
              do_pair = (merge(.true.,iglob.lt.__shfl(kglob,srclane)
     &                   ,same_block).and.(li.ne.lk(klane).or.li.lt.0))
              if (do_pair.and.d2<=cut2.and.accept_mid) then
                  ! compute one interaction
                  call tmatxb_couple(d2,pos,dpui,dpuk(klane)
     &                    ,pdp,pgm,aewald,alsq2,alsq2n,real(1,t_p)
     &                    ,fid,fip,fkd(klane),fkp(klane),.false.)
              end if
 
           end do

           ! increment electric field for each atoms
           rstat = atomicAdd( efi(1,1,iploc),fid%x )
           rstat = atomicAdd( efi(2,1,iploc),fid%y )
           rstat = atomicAdd( efi(3,1,iploc),fid%z )
           rstat = atomicAdd( efi(1,2,iploc),fip%x )
           rstat = atomicAdd( efi(2,2,iploc),fip%y )
           rstat = atomicAdd( efi(3,2,iploc),fip%z )
           call syncwarp(ALL_LANES)
           rstat = atomicAdd( efi(1,1,kploc),fkd(threadIdx%x)%x )
           rstat = atomicAdd( efi(2,1,kploc),fkd(threadIdx%x)%y )
           rstat = atomicAdd( efi(3,1,kploc),fkd(threadIdx%x)%z )
           rstat = atomicAdd( efi(1,2,kploc),fkp(threadIdx%x)%x )
           rstat = atomicAdd( efi(2,2,kploc),fkp(threadIdx%x)%y )
           rstat = atomicAdd( efi(3,2,kploc),fkp(threadIdx%x)%z )
        end do
        end subroutine

        attributes(global) subroutine efld0_direct_scaling_cu
     &            (dpcorrect_ik,dpcorrect_scale,poleloc,ipole,pdamp
     &            ,thole,x,y,z,rpole,ef,n_dpscale,n,npolebloc
     &            ,cut2,aewald,alsq2,alsq2n)
        implicit none
        integer,device,intent(in)::dpcorrect_ik(*),poleloc(n),ipole(n)
        real(t_p),device,intent(in):: dpcorrect_scale(*),pdamp(n)
     &           ,thole(n),rpole(13,n),x(*),y(*),z(*)
        real(t_p),device,intent(inout)::ef(3,2,*)
        integer,value,intent(in)::n,npolebloc,n_dpscale
        real(t_p),value,intent(in)::cut2,aewald,alsq2,alsq2n

        integer ii,ithread,iipole,kpole,iploc,iglob,kbis,kglob
        real(t_p) pdi,pti,d2,damp,thole1,pgamma,d,bn1,bn2,sc3,sc5
     &           ,pscale,dscale,rstat
        type(real3) fid,fip,fkd,fkp,pos
        type(rpole_elt) ip,kp

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x

        do ii = ithread , n_dpscale, blockDim%x*gridDim%x
           iipole = dpcorrect_ik(2*(ii-1)+1)
           kpole  = dpcorrect_ik(2*(ii-1)+2)

           dscale = dpcorrect_scale(2*(ii-1)+1)
           pscale = dpcorrect_scale(2*(ii-1)+2)

           iploc  = poleloc   (iipole)
           iglob  = ipole     (iipole)

           kbis   = poleloc(kpole)
           kglob  = ipole  (kpole)

           pdi    = pdamp(iipole)
           pti    = thole(iipole)

           if (.not.bal_comp.and.
     &        (iploc.lt.1.or.iploc.gt.npolebloc.or.
     &          kbis.lt.1.or.kbis.gt.npolebloc)) cycle

           pos%x  = x(kglob) - x(iglob)
           pos%y  = y(kglob) - y(iglob)
           pos%z  = z(kglob) - z(iglob)
           call image_inl(pos%x,pos%y,pos%z)
           d2     = pos%x**2 + pos%y**2 + pos%z**2
           if (d2.gt.cut2) cycle

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

           kp%c   = rpole(01, kpole)
           kp%dx  = rpole(02, kpole)
           kp%dy  = rpole(03, kpole)
           kp%dz  = rpole(04, kpole)
           kp%qxx = rpole(05, kpole)
           kp%qxy = rpole(06, kpole)
           kp%qxz = rpole(07, kpole)
           kp%qyy = rpole(09, kpole)
           kp%qyz = rpole(10, kpole)
           kp%qzz = rpole(13, kpole)

           thole1 = thole(kpole)
           damp   = pdi * pdamp(kpole)
           pgamma = min( pti,thole1 )

           call efld0_couple(d2,pos,ip,kp,alsq2,alsq2n,
     &                aewald,damp,pgamma,dscale,pscale,
     &                fid,fip,fkd,fkp,d,bn1,bn2,sc3,sc5,.true.)

           if (dscale.ne.0.0_ti_p) then
              rstat = atomicAdd(ef(1,1,iploc), fid%x)
              rstat = atomicAdd(ef(2,1,iploc), fid%y)
              rstat = atomicAdd(ef(3,1,iploc), fid%z)
              rstat = atomicAdd(ef(1,1,kbis ), fkd%x)
              rstat = atomicAdd(ef(2,1,kbis ), fkd%y)
              rstat = atomicAdd(ef(3,1,kbis ), fkd%z)
           end if
           if (pscale.ne.0.0_ti_p) then
              rstat = atomicAdd(ef(1,2,iploc), fip%x)
              rstat = atomicAdd(ef(2,2,iploc), fip%y)
              rstat = atomicAdd(ef(3,2,iploc), fip%z)
              rstat = atomicAdd(ef(1,2,kbis ), fkp%x)
              rstat = atomicAdd(ef(2,2,kbis ), fkp%y)
              rstat = atomicAdd(ef(3,2,kbis ), fkp%z)
           end if
           end do
        end subroutine
      end module
#else
      ! For Portability issue
      subroutine void_tmatxbcu
      end
#endif
