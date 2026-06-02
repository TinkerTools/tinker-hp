!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###########################################################
!     ##                                                       ##
!     ##  subroutine eurey  --  Urey-Bradley potential energy  ##
!     ##                                                       ##
!     ###########################################################
!
!
!     "eurey" calculates the Urey-Bradley 1-3 interaction energy
!
!
#include "tinker_macro.h"
module eurey_inl
#include "atomicOp.h.f90"
contains
#include "ker_urey.inc.f90"
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
   parameter(&
      &grp=__use_groups__,&
      &tver=__use_ene__,&
      &tfea=__use_groups__+__use_polymer__&
      &)

   if(deb_Path) write(*,*) 'eurey'
!
!     zero out the Urey-Bradley interaction energy
!
   eub = 0.0_re_p
!
!     calculate the Urey-Bradley 1-3 energy term
!
!$acc parallel loop present(ureyglob,iury,loc,x,y,z,ul,uk,grplist,wgrp &
!$acc         ,use,deub,eub,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz &
!$acc         ,ureytypI) async &
!$acc         reduction(+:eub) &
!$acc         private(fgrp)
   do iurey = 1, nureyloc
      i     = ureyglob(iurey)
      ureytypii = ureytypI(i)
      ia    = iury(1,i)
      ic    = iury(3,i)
      ideal = ul(i)
      force = uk(i)
!
!     decide whether to compute the current interaction
!
      proceed = (use(ia) .or. use(ic))

      if (IAND(tfea,grp).NE.0.and.use_group)&
         &call groups2_inl (fgrp,ia,ic,ngrp,grplist,wgrp)
!
!     compute the value of the 1-3 distance deviation
!
      if (proceed) then
         call ker_urey(i,ia,ic,nbloc,loc,ideal,force,ureyunit&
            &,cury,qury,fgrp,ureytypii,use_group,use_polymer&
            &,x,y,z&
            &,eub,e,ded,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz&
            &,tver,tfea)
      end if
   end do
end
