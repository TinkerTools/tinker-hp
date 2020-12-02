c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_precision.h"
#include "tinker_cudart.h"

      module tmatxb_pmecu
        use utilcu  ,only: nproc,BLOCK_DIM,ALL_LANES
        use utilgpu ,only: real3,real6,BLOCK_SIZE

        contains

#include "image.f.inc"
#include "midpointimage.f.inc"
#include "pair_tmatxb.f.inc"

        attributes(global) subroutine tmatxb_pme_core_cu
     &                    (ipole,pglob,ploc,ieblst,eblst,x,y,z
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
     &         ,eblst(npolelocnlb_pair*(BLOCK_SIZE+2))
        real(t_p),device,intent(in)::pdamp(n),thole(n),polarity(n)
     &           ,x(npolelocnlb),y(npolelocnlb),z(npolelocnlb)
        real(t_p),device,intent(in):: mu(3,2,npolebloc)
        real(t_p),device,intent(inout):: efi(3,2,npolebloc)

        integer ithread,iwarp,nwarp,ilane,klane,istat,srclane
        integer beg,ii,j
        integer iblock,idx,kdx,kdx_
        integer iipole,iglob,iploc,kpole,kploc
        integer,shared:: kglob(BLOCK_DIM)
        real(t_p) xk_,yk_,zk_,d2
        real(t_p) ipdp,ipgm,pdp,pgm,rstat
        type(real6) dpui,dpuk_
        type(real3) posi,pos
        type(real3) fid,fip
        type(real3),shared::fkd(BLOCK_DIM),fkp(BLOCK_DIM)
     &             ,posk(BLOCK_DIM)
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
              do_pair = merge(.true.,iglob.lt.kglob(klane)
     &                       ,same_block)
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
      end module
#else
      ! For Portability issue
      subroutine void_tmatxbcu
      end
#endif
