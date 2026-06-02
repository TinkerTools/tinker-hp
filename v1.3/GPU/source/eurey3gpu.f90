!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine eurey3  --  Urey-Bradley energy & analysis  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "eurey3" calculates the Urey-Bradley energy; also
!     partitions the energy among the atoms
!
!
#include "tinker_macro.h"
module eurey3gpu_inl
#include "atomicOp.h.f90"
contains
#include "ker_urey.inc.f90"
end module

subroutine eurey3gpu
   use action
   use analyz
   use atmlst
   use atoms
   use bound
   use deriv
   use domdec
   use energi
   use eurey3gpu_inl
   use group
   use inform    ,only: verbose,debug,deb_Path
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
   logical proceed,header,huge
   parameter(&
      &grp=__use_groups__,&
      &tver=__use_ene__,&
      &tfea=__use_groups__+__use_polymer__&
      &)

   if(deb_Path) write(*,*) 'eurey3gpu'
!
!     zero out the Urey-Bradley interaction energy
!
   eub    = 0.0_re_p
   neub   = 0
   header = rank.eq.0
!
!     calculate the Urey-Bradley 1-3 energy term
!
!$acc parallel loop present(ureyglob,iury,loc,x,y,z,ul,uk,grplist,wgrp &
!$acc         ,use,eub,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz &
!$acc         ,ureytypI) async &
!$acc         reduction(+:eub,neub) &
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

      if (use_group.and.IAND(tfea,grp).NE.0)&
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
         if(e.ne.0.0) neub = neub+1
#if 0
         ib = loc(ib)
         ic = loc(ic)
         aub(ib) = aub(ib) + 0.5_ti_p*e
         aub(ic) = aub(ic) + 0.5_ti_p*e
!
!     print a message if the energy of this interaction is large
!
         huge = (e .gt. 5.0_ti_p)
         if (debug .or. (verbose.and.huge)) then
            if (header) then
               header = .false.
               write (iout,10)
10             format (/,' Individual Urey-Bradley Interactions :',&
                  &//,' Type',18x,'Atom Names',18x,'Ideal',&
                  &4x,'Actual',6x,'Energy',/)
            end if
            write (iout,20)  ia,name(ia),ib,name(ib),&
               &ic,name(ic),ideal,rac,e
20          format (' UreyBrad',2x,i7,'-',a3,i7,'-',a3,&
               &i7,'-',a3,2x,2f10.4,f12.4)
         end if
#endif
      end if
   end do
end
