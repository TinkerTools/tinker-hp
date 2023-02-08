c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#ifdef _CUDA
#define TINKER_CUF
#include "tinker_macro.h"
#include "tinker_cudart.h"

      module efld0_cpencu
      use utilcu  ,only: nproc,BLOCK_DIM,ALL_LANES
     &            ,bal_comp=>balanced_comput
      use utilgpu ,only: BLOCK_SIZE
      use tintypes,only: rpole_elt,real3,real6,real7
      use neigh   ,only: i_bl=>inside_block,d_bl=>disc_block
      private

      public :: efld0_cpen_kcu
      contains

#include "image.f.inc"
#include "midpointimage.f.inc"
#include "pair_efld_cp.inc.f"

      attributes(global) subroutine efld0_cpen_kcu
     &          (ipole,pglob,ploc,bapl,abpl
     &          ,x,y,z,pcore,pval,palpha,rpole
     &          ,na,nab,nbap,npolebloc,n,pentyp,cut2,aewald
     &          ,efi
     &          ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
     &          ,scal_ik,scal_val,poleloc,ipole_,x_,y_,z_,n_dpscale)
      implicit none
      integer  ,value,intent(in):: na,nab,npolebloc,n,nbap,pentyp
     &         ,n_dpscale
      real(t_p),value,intent(in):: p_xbeg,p_xend,p_ybeg,p_yend
     &         ,p_zbeg,p_zend,cut2,aewald
      integer  ,device,intent(in)::ipole(nab),pglob(nab)
     &         ,ploc(nab),bapl(nbap)
     &         ,abpl(nbap*(BLOCK_SIZE+2))
     &         ,scal_ik(*),poleloc(*),ipole_(*)
      real(t_p),device,intent(in)::pcore(n),pval(n),palpha(n)
     &         ,x(nab),y(nab),z(nab),scal_val(*),x_(*),y_(*),z_(*)
     &         ,rpole(13,*)
      real(t_p),device,intent(inout):: efi(3,2,npolebloc)

      integer ithread,iwarp,nwarp,ilane,klane,nlane,istat,srclane
      integer beg,ii,j,iblock,idx,kdx,kdx_,ver,ver1,fea
      integer iipole,iglob,kglob,iploc,kpole,kploc
      real(t_p) xk_,yk_,zk_,d2,rstat
      real(t_p) corei,vali,alphai
      real(t_p)  ,shared:: corek(BLOCK_DIM),valk(BLOCK_DIM)
     &           ,alphak(BLOCK_DIM)
      type(real3),shared:: posi(BLOCK_DIM)
      type(real3) pos,posk
      type(rpole_elt) :: ip
      type(rpole_elt),shared :: kp(BLOCK_DIM)
      real(t_p) fid(3),fip(3),fkd(3),fkp(3)
      logical do_pair
      parameter( ver=0
     &         ,ver1=__use_grd__+__use_sca__
     &         ,fea=0 )

      ithread = threadIdx%x + (blockIdx%x-1)*blockDim%x
      iwarp   = (ithread-1) / warpsize
      nwarp   = blockDim%x*gridDim%x / warpsize
      ilane   = iand( threadIdx%x-1,warpsize-1 ) + 1
      nlane   = merge(1,ilane+1,ilane.eq.warpsize)
c
c     Correct Scaling Interactions
c
      do ii = ithread , n_dpscale, blockDim%x*gridDim%x; block
         real(t_p) dscale,pscale
         iipole = scal_ik(2*(ii-1)+1)
         kpole  = scal_ik(2*(ii-1)+2)

         dscale = scal_val(2*(ii-1)+1)
         pscale = scal_val(2*(ii-1)+2)

         iploc  = poleloc   (iipole)
         iglob  = ipole_    (iipole)

         kploc  = poleloc(kpole)
         kglob  = ipole_ (kpole)

         !pdi    = pdamp(iipole)
         !pti    = thole(iipole)
         corei  = pcore(iipole)
         vali   = pval(iipole)
         alphai = palpha(iipole)

         if (.not.bal_comp.and.
     &      (iploc.lt.1.or.iploc.gt.npolebloc.or.
     &       kploc.lt.1.or.kploc.gt.npolebloc)) cycle

         pos%x  = x_(kglob) - x_(iglob)
         pos%y  = y_(kglob) - y_(iglob)
         pos%z  = z_(kglob) - z_(iglob)
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

         kp(threadIdx%x)%c   = rpole(01, kpole)
         kp(threadIdx%x)%dx  = rpole(02, kpole)
         kp(threadIdx%x)%dy  = rpole(03, kpole)
         kp(threadIdx%x)%dz  = rpole(04, kpole)
         kp(threadIdx%x)%qxx = rpole(05, kpole)
         kp(threadIdx%x)%qxy = rpole(06, kpole)
         kp(threadIdx%x)%qxz = rpole(07, kpole)
         kp(threadIdx%x)%qyy = rpole(09, kpole)
         kp(threadIdx%x)%qyz = rpole(10, kpole)
         kp(threadIdx%x)%qzz = rpole(13, kpole)

         corek (threadIdx%x) = pcore (kpole)
         valk  (threadIdx%x) = pval  (kpole)
         alphak(threadIdx%x) = palpha(kpole)
         !thole1 = thole(kpole)
         !damp   = pdi * pdamp(kpole)
         !pgamma = min( pti,thole1 )

         call duo_efld0Cpen (ip,kp(threadIdx%x),d2,pos,pscale,dscale
     &           ,aewald,pentyp,corei,corek(threadIdx%x)
     &           ,vali,valk(threadIdx%x),alphai,alphak(threadIdx%x)
     &           ,fid,fip,fkd,fkp,ver1,fea)

         if (dscale.ne.0.0_ti_p) then
            rstat = atomicAdd(efi(1,1,iploc), fid(1))
            rstat = atomicAdd(efi(2,1,iploc), fid(2))
            rstat = atomicAdd(efi(3,1,iploc), fid(3))
            rstat = atomicAdd(efi(1,1,kploc), fkd(1))
            rstat = atomicAdd(efi(2,1,kploc), fkd(2))
            rstat = atomicAdd(efi(3,1,kploc), fkd(3))
         end if
         if (pscale.ne.0.0_ti_p) then
            rstat = atomicAdd(efi(1,2,iploc), fip(1))
            rstat = atomicAdd(efi(2,2,iploc), fip(2))
            rstat = atomicAdd(efi(3,2,iploc), fip(3))
            rstat = atomicAdd(efi(1,2,kploc), fkp(1))
            rstat = atomicAdd(efi(2,2,kploc), fkp(2))
            rstat = atomicAdd(efi(3,2,kploc), fkp(3))
         end if
      end block; end do
c
c       Compute pairwise interaction
c
      do ii = iwarp, nbap-1, nwarp

         ! Load atom block k parameters
         kdx     = abpl(ii*warpsize + ilane)
         kpole   = pglob(kdx)
         kglob   = ipole(kdx)
         kploc   = ploc(kdx)
         posk%x  = x(kdx)
         posk%y  = y(kdx)
         posk%z  = z(kdx)
         kp(threadIdx%x)%c   = rpole(01,kpole)
         kp(threadIdx%x)%dx  = rpole(02,kpole)
         kp(threadIdx%x)%dy  = rpole(03,kpole)
         kp(threadIdx%x)%dz  = rpole(04,kpole)
         kp(threadIdx%x)%qxx = rpole(05,kpole)
         kp(threadIdx%x)%qxy = rpole(06,kpole)
         kp(threadIdx%x)%qxz = rpole(07,kpole)
         kp(threadIdx%x)%qyy = rpole(09,kpole)
         kp(threadIdx%x)%qyz = rpole(10,kpole)
         kp(threadIdx%x)%qzz = rpole(13,kpole)
         corek(threadIdx%x)    = pcore(kpole)
         valk (threadIdx%x)    = pval (kpole)
         alphak(threadIdx%x)   = palpha(kpole)

         ! Load atom block i parameters
         iblock  = bapl(ii+1)
         if (iblock.eq.0) cycle
         idx     = (iblock-1)*warpsize + ilane
         iipole  = pglob(idx)
         iglob   = ipole(idx)
         iploc   = ploc(idx)
         posi(threadIdx%x)%x = x(idx)
         posi(threadIdx%x)%y = y(idx)
         posi(threadIdx%x)%z = z(idx)
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
         corei   = pcore (iipole)
         vali    = pval  (iipole)
         alphai  = palpha(iipole)

         !set compute Data to 0
         fid(1)  = 0; fip(1)  = 0
         fid(2)  = 0; fip(2)  = 0
         fid(3)  = 0; fip(3)  = 0
         fkd(1)  = 0; fkp(1)  = 0
         fkd(2)  = 0; fkp(2)  = 0
         fkd(3)  = 0; fkp(3)  = 0

         do j = 0,warpsize-1
            srclane  = iand( ilane+j-1,warpsize-1 ) + 1
            klane    = threadIdx%x-ilane + srclane
            !kdx_     = __shfl(kdx ,srclane)
            !pdp      = ipdp*kpdp(klane)
            !pgm      = min(ipgm,kpgm(klane))

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
     &         .or.(yk_.lt.p_ybeg).or.(yk_.ge.p_yend)
     &         .or.(xk_.lt.p_xbeg).or.(xk_.ge.p_xend)) then
                  do_pair = .false.
               end if
               pos%x=-pos%x; pos%y=-pos%y; pos%z=-pos%z;
               end block
            else
               pos%x = posk%x - posi(threadIdx%x)%x
               pos%y = posk%y - posi(threadIdx%x)%y
               pos%z = posk%z - posi(threadIdx%x)%z
               call image_inl(pos%x,pos%y,pos%z)
            end if
#else
            pos%x =  posk%x - posi(threadIdx%x)%x
            pos%y =  posk%y - posi(threadIdx%x)%y
            pos%z =  posk%z - posi(threadIdx%x)%z
            call image_inl(pos%x,pos%y,pos%z)
#endif
            d2      = pos%x**2 + pos%y**2 + pos%z**2

            do_pair = do_pair.and.d2<=cut2
            if (do_pair) then
                ! compute one interaction
                call duo_efld0Cpen2 (ip,kp(klane),d2,pos,1.0,1.0,aewald
     &                  ,pentyp,corei,corek(klane),vali,valk(klane)
     &                  ,alphai,alphak(klane)
     &                  ,fid,fip,fkd,fkp,ver,fea)
            end if

            kglob  = __shfl(  kglob,nlane )
            posk%x = __shfl( posk%x,nlane )
            posk%y = __shfl( posk%y,nlane )
            posk%z = __shfl( posk%z,nlane )
             fkd(1)= __shfl(  fkd(1),nlane )
             fkd(2)= __shfl(  fkd(2),nlane )
             fkd(3)= __shfl(  fkd(3),nlane )
             fkp(1)= __shfl(  fkp(1),nlane )
             fkp(2)= __shfl(  fkp(2),nlane )
             fkp(3)= __shfl(  fkp(3),nlane )
         end do

         ! increment electric field for each atoms
         if (idx.le.na) then
         rstat = atomicAdd( efi(1,1,iploc),fid(1) )
         rstat = atomicAdd( efi(2,1,iploc),fid(2) )
         rstat = atomicAdd( efi(3,1,iploc),fid(3) )
         rstat = atomicAdd( efi(1,2,iploc),fip(1) )
         rstat = atomicAdd( efi(2,2,iploc),fip(2) )
         rstat = atomicAdd( efi(3,2,iploc),fip(3) )
         end if
         if (kdx.le.na) then
         rstat = atomicAdd( efi(1,1,kploc),fkd(1) )
         rstat = atomicAdd( efi(2,1,kploc),fkd(2) )
         rstat = atomicAdd( efi(3,1,kploc),fkd(3) )
         rstat = atomicAdd( efi(1,2,kploc),fkp(1) )
         rstat = atomicAdd( efi(2,2,kploc),fkp(2) )
         rstat = atomicAdd( efi(3,2,kploc),fkp(3) )
         end if
      end do
      end subroutine

      end module
#endif
