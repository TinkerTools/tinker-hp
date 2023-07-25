c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine eurey  --  Urey-Bradley potential energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "eurey" calculates the Urey-Bradley 1-3 interaction energy
c
c
#include "tinker_macro.h"
      module eurey_inl
#include "atomicOp.h.f"
      contains
#include "ker_urey.inc.f"
      end module

      subroutine eurey
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use eurey_inl
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
     &        tver=__use_ene__,
     &        tfea=__use_groups__+__use_polymer__
     &         )

      if(deb_Path) write(*,*) 'eurey'
c
c     zero out the Urey-Bradley interaction energy
c
      eub = 0.0_re_p
c
c     calculate the Urey-Bradley 1-3 energy term
c
!$acc parallel loop present(ureyglob,iury,loc,x,y,z,ul,uk,grplist,wgrp
!$acc&         ,use,deub,eub,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
!$acc&         ,ureytypI) async
!$acc&         reduction(+:eub)
!$acc&         private(fgrp)
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
     &              ,cury,qury,fgrp,ureytypii,use_group,use_polymer
     &              ,x,y,z
     &              ,eub,e,ded,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &              ,tver,tfea)
         end if
      end do
      end
