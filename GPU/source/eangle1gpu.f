c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangle1  --  angle bend energy and derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangle1" calculates the angle bending potential energy and
c     the first derivatives with respect to Cartesian coordinates;
c     projected in-plane angles at trigonal centers, special linear
c     or Fourier angle bending terms are optionally used
c
c
c
c     Comment
c     Zhi Wang, July 1, 2019
c     
c     The original Tinker implementation has following expressions
c     xip = xib + xt * delta
c     xap = xia - xip
c     
c     And they were reorganized to
c     xap = xia - xib - xt * delta
c     for higher accuracy in the single precision mode.
c     
c     Consider an example where
c     xia = 33.553368, xib = 34.768604
c     xt = 0.33142909, delta = 0.0044494048,
c     the later expression gives a better numerical result.
c
#include "tinker_macro.h"
      module eangle1gpu_inl
#include "atomicOp.h.f"
        contains
#include "ker_angle.inc.f"
      end module

      subroutine eangle1gpu_(dea,deW1aMD,atmType,afld,s)
      use angle     ,only: nangleloc,iang,anat,ak
      use angpot
      use atmlst
      use atoms     ,only: n,x,y,z
      use bound
      !use deriv
      use domdec
      use energi
      use eangle1gpu_inl
      use group
      use inform,only: deb_Path,minmaxone
      use math
      use potent,only:use_amd_wat1
      use usage
      use virial
      use timestat,only:timer_enter,timer_exit,timer_eangle
      use tinheader
      use mamd
      implicit none
      integer  ,intent(in)   :: atmType(:),s
      real(t_p),intent(in)   :: afld(:)
      real(r_p),intent(inout):: dea(:,:),deW1aMD(3,s)
      integer i,ia,ib,ic,id,iangle,angtypii,tfea,tver
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      logical proceed
      real(t_p) ideal,force,e,fgrp
      parameter(
     &    tfea=__use_gamd__+__use_groups__+__use_polymer__,
     &    tver=__use_grd__+__use_ene__+__use_vir__
     &    )

      if(deb_Path) write(*,*) 'eangle1gpu'
      call timer_enter( timer_eangle )
c
c     calculate the bond angle bending energy term
c
!$acc parallel loop async
#ifdef USE_NVSHMEM_CUDA
!$acc&     present(x,y,z,loc,use,angleglob,anat,
#else
!$acc&     present(x,y,z,loc,use,iang,angleglob,anat,
#endif
!$acc&     ak,afld,angtypI,grplist,wgrp,atmType,aMDwattype,deW1aMD,
!$acc&     dea,ea,eW1aMD,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz,
!$acc&     fmat_ps,dfmat_ps)
!$acc&     reduction(+:ea,eW1aMD,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&     private(fgrp)
      do iangle = 1, nangleloc
         i   = angleglob(iangle)
#ifdef USE_NVSHMEM_CUDA
         ipe =     (i-1)/nangle_pe
         ind = mod((i-1),nangle_pe) +1
         ia  = d_iang(ipe)%pel(1,ind)
         ib  = d_iang(ipe)%pel(2,ind)
         ic  = d_iang(ipe)%pel(3,ind)
         id  = d_iang(ipe)%pel(4,ind)
#else
         ia  = iang(1,i)
         ib  = iang(2,i)
         ic  = iang(3,i)
         id  = iang(4,i)
#endif
         ideal    = anat(i)
         force    = ak(i)
         angtypii = angtypI(i)
c
c     decide whether to compute the current interaction
c
         if (angtypii .eq. ANG_IN_PLANE) then
            if(use_group)
     &         call groups4_inl(fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
            proceed = (use(ia).or. use(ib).or. use(ic).or. use(id))
         else
            if (use_group)
     &         call groups3_inl (fgrp,ia,ib,ic,ngrp,grplist,wgrp)
            proceed = (use(ia).or. use(ib).or. use(ic))
         end if
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            if (angtypii .ne. ANG_IN_PLANE) then
              !compute the bond angle bending energy and gradient
               call ker_angle(i,ia,ib,ic,loc,ideal,force
     &                 ,angunit,cang,pang,sang,qang,fgrp
     &                 ,use_group,use_polymer,angtypii,x,y,z,afld
     &                 ,ea,e,dea,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz
     &                 ,g_vzz,tver,tfea,fmat_ps,dfmat_ps)
            else
               call ker_angle_plan(i,ia,ib,ic,id,loc,ideal,force
     &                 ,angunit,cang,pang,sang,qang,fgrp
     &                 ,use_group,use_polymer,angtypii,x,y,z,afld
     &                 ,use_amd_wat1,atmType,aMDwattype,eW1aMD,ea,e
     &                 ,dea,deW1aMD,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                 ,tver,tfea)
            end if
         end if
      end do
      call timer_exit( timer_eangle )
      end

      subroutine eangle1gpu
      use angle
      use atoms
      use deriv
      implicit none
      interface
      subroutine eangle1gpu_(dea,deW1aMD,atmType,afld,s)
       integer  ,intent(in)   :: atmType(:),s
       real(t_p),intent(in)   :: afld(:)
       real(r_p),intent(inout):: dea(:,:),deW1aMD(3,s)
      end subroutine
      end interface
      call eangle1gpu_(dea,deW1aMD,type,afld,max(1,size(deW1aMD,2)))
      end subroutine
