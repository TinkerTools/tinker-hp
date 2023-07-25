c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eurey1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eurey1" calculates the Urey-Bradley interaction energy and
c     its first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module eurey1gpu_inl
#include "atomicOp.h.f"
      contains
#include "ker_urey.inc.f"
      end module

      subroutine eurey1gpu
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use eurey1gpu_inl
      use group
      use inform    ,only: deb_Path
      use tinheader ,only: ti_p,re_p
      use tinTypes  ,only: real3
      use urey
      use urypot
      use usage
      use virial
      use timestat
      implicit none
      integer i,ia,ic,iurey,grp,tver,tfea,ureytypii
      real(t_p) ideal,force,fgrp,e
      type(real3) ded
      logical proceed
      parameter(
     &         grp=__use_groups__,
     &        tver=__use_grd__+__use_ene__+__use_vir__,
     &        tfea=__use_groups__+__use_polymer__
     &         )

      if(deb_Path) write(*,*) 'eurey1gpu'
      call timer_enter( timer_eurey1 )
c
c     calculate the Urey-Bradley 1-3 energy and first derivatives
c
!$acc parallel loop present(ureyglob,iury,loc,x,y,z,ul,uk,grplist,wgrp
!$acc&         ,use,deub,eub,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
!$acc&         ,ureytypI) async
!$acc&         reduction(+:eub,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         private(ded,fgrp)
      do iurey = 1, nureyloc
         i     = ureyglob(iurey)
         ureytypii = ureytypI(i)
         ia    = iury(1,i)
         ic    = iury(3,i)
         ideal = ul(i)
         force = uk(i)
c
c     decide whether to compute the current interaction
c
         proceed = (use(ia) .or. use(ic))

         if (IAND(tfea,grp).NE.0.and.use_group)
     &      call groups2_inl (fgrp,ia,ic,ngrp,grplist,wgrp)
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            call ker_urey(i,ia,ic,nbloc,loc,ideal,force,ureyunit
     &               ,cury,qury,fgrp,ureytypii,use_group,use_polymer
     &               ,x,y,z
     &               ,eub,e,ded,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &               ,tver,tfea)

           !increment the total Urey-Bradley first derivatives
            call atomic_add( deub(1,ia),ded%x )
            call atomic_add( deub(2,ia),ded%y )
            call atomic_add( deub(3,ia),ded%z )
            call atomic_sub( deub(1,ic),ded%x )
            call atomic_sub( deub(2,ic),ded%y )
            call atomic_sub( deub(3,ic),ded%z )
         end if
      end do
c
      call timer_exit( timer_eurey1 )
      end
