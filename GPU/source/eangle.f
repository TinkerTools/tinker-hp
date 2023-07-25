c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine eangle  --  angle bending potential energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "eangle" calculates the angle bending potential energy;
c     projected in-plane angles at trigonal centers, special
c     linear or Fourier angle bending terms are optionally used
c
c
#include "tinker_macro.h"
      module eanglegpu_inl
#include "atomicOp.h.f"
        contains
#include "ker_angle.inc.f"
      end module

      subroutine eangle_(dea,deW1aMD,atmType,afld,lena)
      use analyz
      use angle    ,only: nangleloc,iang,anat,ak
      use angpot
      use atmlst
      use atmtyp
      use atoms    ,only: n,x,y,z
      use bound
      use domdec
      use energi
      use eanglegpu_inl
      use group
      use inform
      use iounit
      use math
      use mamd
      use potent  ,only: use_amd_wat1
      use usage
      use tinheader
      use timestat,only:timer_enter,timer_exit,timer_eangle
      use virial
      implicit none
      integer  ,intent(in)   :: lena
      integer  ,intent(in)   :: atmType(n)
      real(t_p),intent(inout):: afld(lena)
      real(r_p),intent(inout):: dea(1),deW1aMD(1)
      integer   i,ia,ib,ic,id,iangle,angtypii,tfea,tver
      logical   proceed,header,huge
#ifdef USE_NVSHMEM_CUDA
      integer   ipe,ind
#endif
      real(t_p) ideal,force,e,fgrp
      character*9 label
      parameter(
     &    tfea=__use_groups__+__use_polymer__,
     &    tver=__use_ene__
     &    )

      if(deb_Path) write(*,*) 'eanglegpu'
      call timer_enter( timer_eangle )
c
c     zero out the angle bending energy and partitioning terms
c
      ea  = 0.0_re_p
c     aea = 0.0_ti_p
      header = rank.eq.0
c
c     calculate the bond angle bending energy term
c
!$acc parallel loop 
#ifdef USE_NVSHMEM_CUDA
!$acc&         present(x,y,z,loc,use,aea,angleglob
#else
!$acc&         present(x,y,z,loc,use,iang,aea,angleglob
#endif
!$acc&    ,anat,ak,afld,angtypI,grplist,wgrp,atmType,aMDwattype,deW1aMD
!$acc&    ,dea,ea,eW1aMD,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
!$acc&    ,fmat_ps,dfmat_ps)
!$acc&         present(ea) async reduction(+:ea)
      do iangle = 1, nangleloc
         i     = angleglob(iangle)
#ifdef USE_NVSHMEM_CUDA
         ipe   =     (i-1)/nangle_pe
         ind   = mod((i-1),nangle_pe) +1
         ia    = d_iang(ipe)%pel(1,ind)
         ib    = d_iang(ipe)%pel(2,ind)
         ic    = d_iang(ipe)%pel(3,ind)
#else
         ia    = iang(1,i)
         ib    = iang(2,i)
         ic    = iang(3,i)
#endif
         ideal = anat(i)
         force = ak(i)
         angtypii = angtypI(i)
c
c     decide whether to compute the current interaction
c
         if (angtypii .eq. ANG_IN_PLANE) then
            id   = iang(4,i)
            if(use_group)
     &         call groups4_inl(fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
            proceed = (use(ia).or. use(ib).or. use(ic).or. use(id))
         else
            if(use_group)
     &         call groups3_inl(fgrp,ia,ib,ic,ngrp,grplist,wgrp)
            proceed = (use(ia).or. use(ib).or. use(ic))
         end if
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
c
c     compute the bond angle bending energy
c
            if (angtypii .ne. ANG_IN_PLANE) then
 
               call ker_angle(i,ia,ib,ic,loc,ideal,force
     &                 ,angunit,cang,pang,sang,qang,fgrp
     &                 ,use_group,use_polymer,angtypii,x,y,z,afld
     &                 ,ea,e,dea,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz
     &                 ,g_vzz,tver,tfea,fmat_ps,dfmat_ps)
c
c     compute the projected in-plane angle bend energy
c
            else
               call ker_angle_plan(i,ia,ib,ic,id,loc,ideal,force
     &                 ,angunit,cang,pang,sang,qang,fgrp
     &                 ,use_group,use_polymer,angtypii,x,y,z,afld
     &                 ,use_amd_wat1,atmType,aMDwattype,eW1aMD,ea,e
     &                 ,deW1aMD,dea,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &                 ,tver,tfea)
            end if
         end if
      end do
      call timer_exit( timer_eangle )
      end

      subroutine eangle
      use angle
      use atoms
      use deriv
      implicit none
      call eangle_(dea,deW1aMD,type,afld,size(afld))
      end subroutine
