#ifdef _CUDA
#define TINKER_CUF
#include "tinker_macro.h"
#include "tinker_cudart.h"

      module tmatxb_pme_cpencu
        use utilcu  ,only: nproc,BLOCK_DIM,ALL_LANES
     &              ,bal_comp=>balanced_comput
        use utilgpu ,only: BLOCK_SIZE
        use tintypes,only: rpole_elt,real3,real6,real7
        use neigh   ,only: i_bl=>inside_block,d_bl=>disc_block
        private
        public :: tmatxb_pme_cpen_kcu

        contains

#include "image.f.inc"
#include "midpointimage.f.inc"
#include "pair_tmatxb.f.inc"
#define __tfea__ (__use_mpi__)

        attributes(global) subroutine tmatxb_pme_cpen_kcu
     &            (ipole,pglob,b_stat,ploc,ieblst,eblst,x,y,z
     &            ,palpha,polarity,mu,efi
     &            ,na,nab,nbap,nabloc,n,pentyp,icall
     &            ,cut2,alsq2,alsq2n,aewald
     &            ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &            ,scal_ik,scal_val,poleloc,ipole_,x_,y_,z_,n_scale,efi_
     &            )
        implicit none
        integer,value,intent(in):: na,nab,nabloc,n
     &         ,nbap,pentyp,icall,n_scale
        real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &           ,p_zbeg,p_zend,cut2,alsq2,alsq2n,aewald
        integer,device,intent(in)::ipole(nab),pglob(nab)
     &         ,ploc(nab),ieblst(nbap)
     &         ,eblst(nbap*(BLOCK_SIZE+2)),b_stat(*)
     &         ,scal_ik(*),poleloc(*),ipole_(*)
        real(t_p),device,intent(in)::palpha(n),polarity(n)
     &           ,x(nab),y(nab),z(nab)
     &           ,scal_val(*),x_(*),y_(*),z_(*)
        real(t_p),device,intent(in):: mu(3,2,nabloc)
        real(t_p),device,intent(inout):: efi(nabloc,3,2),efi_(3,2,*)

        integer ithread,iwarp,nwarp,ilane,klane,nlane,istat,srclane
        integer beg,ii,j,ver1,ver,fea
        integer iblock,idx,kdx,kdx_
        integer iipole,iglob,kglob,iploc,kpole,kploc
        real(t_p) xk_,yk_,zk_,d2
        real(t_p) alphai,alphak,rstat
        type(real6) dpui
        type(real3),shared:: posi(BLOCK_DIM)
        type(real3) pos,posk
        type(real3) fid,fip
        type(real3) fkd,fkp
        type(real7),shared:: dpuk(BLOCK_DIM)
        logical do_pair
        parameter( ver=__use_grd__
     &           ,ver1=ver+__use_sca__,fea=__use_mpi__)

        ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
        iwarp   = (ithread-1) / warpsize
        nwarp   = blockDim%x*gridDim%x / warpsize
        ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
        nlane   = merge(1,ilane+1,ilane.eq.warpsize)
        !beg     = 2*nbap
c
c     Correct Scaling Interactions
c
        if (icall.eq.1) then
        do ii = ithread, n_scale, blockDim%x*gridDim%x; block
        real(t_p) uscal
           iipole  = scal_ik(2*(ii-1)+1)
           kpole   = scal_ik(2*(ii-1)+2)
           uscal   = scal_val(ii)
           
           iglob   = ipole_ (iipole)
           kglob   = ipole_ (kpole)
           iploc   = poleloc(iipole)
           kploc   = poleloc(kpole)
           !skip atom if it is not polarizable
           !FIXME do we need to test on polarity
           if (polarity(iipole) == 0) cycle
           if (iploc.eq.0.or.iploc.gt.nabloc) cycle
           if (kploc.eq.0.or.kploc.gt.nabloc) cycle
           
           dpui%x  = mu(1,1,iploc)
           dpui%y  = mu(2,1,iploc)
           dpui%z  = mu(3,1,iploc)
           dpui%xx = mu(1,2,iploc)
           dpui%yy = mu(2,2,iploc)
           dpui%zz = mu(3,2,iploc)
           
           dpuk%x  = mu(1,1,kploc)
           dpuk%y  = mu(2,1,kploc)
           dpuk%z  = mu(3,1,kploc)
           dpuk%xx = mu(1,2,kploc)
           dpuk%yy = mu(2,2,kploc)
           dpuk%zz = mu(3,2,kploc)
           
           pos%x   = x_(kglob) - x_(iglob)
           pos%y   = y_(kglob) - y_(iglob)
           pos%z   = z_(kglob) - z_(iglob)
           call image_inl(pos%x,pos%y,pos%z)
           d2  = pos%x*pos%x + pos%y*pos%y + pos%z*pos%z
           if (d2 > cut2) cycle
           
           alphai   = palpha(iipole)
           alphak   = palpha( kpole)

           call duo_tmatxb_cpen(d2,pos,dpui,dpuk(threadIdx%x)
     &             ,pentyp,alphai,alphak,aewald,uscal
     &             ,fid,fip,fkd,fkp,ver1,fea)

           ! increment electric field for each atoms
           rstat= atomicAdd( efi_(1,1,iploc),fid%x )
           rstat= atomicAdd( efi_(2,1,iploc),fid%y )
           rstat= atomicAdd( efi_(3,1,iploc),fid%z )
           rstat= atomicAdd( efi_(1,2,iploc),fip%x )
           rstat= atomicAdd( efi_(2,2,iploc),fip%y )
           rstat= atomicAdd( efi_(3,2,iploc),fip%z )
           rstat= atomicAdd( efi_(1,1,kploc),fkd%x )
           rstat= atomicAdd( efi_(2,1,kploc),fkd%y )
           rstat= atomicAdd( efi_(3,1,kploc),fkd%z )
           rstat= atomicAdd( efi_(1,2,kploc),fkp%x )
           rstat= atomicAdd( efi_(2,2,kploc),fkp%y )
           rstat= atomicAdd( efi_(3,2,kploc),fkp%z )
        end block; end do
        end if
c
c       Compute pairwise interactions
c
        do ii = iwarp, nbap-1, nwarp

           ! Load atom block k parameters
           kdx     = eblst(ii*warpsize + ilane)
           kpole   = pglob(kdx)
           kglob   = ipole(kdx)
           kploc   = ploc(kdx)
           posk%x  = x(kdx)
           posk%y  = y(kdx)
           posk%z  = z(kdx)
           alphak  = palpha(kpole)
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
           alphai  = palpha(iipole)
           dpui%x  = mu(1,1,iploc)
           dpui%y  = mu(2,1,iploc)
           dpui%z  = mu(3,1,iploc)
           dpui%xx = mu(1,2,iploc)
           dpui%yy = mu(2,2,iploc)
           dpui%zz = mu(3,2,iploc)

           !set compute Data to 0
           fid%x   = 0; fkd%x = 0
           fid%y   = 0; fkd%y = 0
           fid%z   = 0; fkd%z = 0
           fip%x   = 0; fkp%x = 0
           fip%y   = 0; fkp%y = 0
           fip%z   = 0; fkp%z = 0

           do j = 0,warpsize-1
              srclane  = iand( ilane+j-1,warpsize-1 ) + 1
              klane    = threadIdx%x-ilane + srclane
              !kdx_     = __shfl(kdx ,srclane)

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
                  call duo_tmatxb_cpen2(d2,pos,dpui,dpuk(klane)
     &                    ,pentyp,alphai,alphak,aewald,1.0
     &                    ,fid,fip,fkd,fkp,ver,fea)
              end if

               kglob = __shfl(  kglob,nlane )
              alphak = __shfl( alphak,nlane )
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
           if (idx.le.na) then
              rstat = atomicAdd( efi(idx,1,1),fid%x )
              rstat = atomicAdd( efi(idx,2,1),fid%y )
              rstat = atomicAdd( efi(idx,3,1),fid%z )
              rstat = atomicAdd( efi(idx,1,2),fip%x )
              rstat = atomicAdd( efi(idx,2,2),fip%y )
              rstat = atomicAdd( efi(idx,3,2),fip%z )
           end if
           if (kdx.le.na) then
              rstat = atomicAdd( efi(kdx,1,1),fkd%x )
              rstat = atomicAdd( efi(kdx,2,1),fkd%y )
              rstat = atomicAdd( efi(kdx,3,1),fkd%z )
              rstat = atomicAdd( efi(kdx,1,2),fkp%x )
              rstat = atomicAdd( efi(kdx,2,2),fkp%y )
              rstat = atomicAdd( efi(kdx,3,2),fkp%z )
           end if
        end do
        end subroutine

      end module
#endif
