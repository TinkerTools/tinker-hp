!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine estrbnd  --  stretch-bend cross term energy  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "estrbnd" calculates the stretch-bend potential energy
!
!
#include "tinker_macro.h"
module estrbnd_inl
#include "atomicOp.h.f90"
contains
#include "ker_strbnd.inc.f90"
end module

subroutine estrbnd_(isb,typeA,anat,bl,deba,deaMD)
   use angle     ,only: iang
   use angpot    ,only: stbnunit
   use atmlst    ,only: strbndglob
   use atoms     ,only: x,y,z
   use bound
   use domdec
   use energi
   use estrbnd_inl
   use group
   use inform
   use iounit
   use mamd
   use nvshmem
   use potent    ,only: use_amd_wat1
   use strbnd    ,only: nstrbndloc,sbk
   use tinheader ,only: ti_p,re_p
   use usage
   use sizes     ,only: tinkerdebug
   use virial
   implicit none
   integer  ,intent(in):: isb(:,:),typeA(:)
   real(t_p),intent(in):: anat(:),bl(:)
   real(r_p),intent(out):: deba(1),deaMD(1)

   integer   grp,tver,tfea
   integer   i,ia,ib,ic,istrbnd,istrbndloc,ipe,ind
   logical   proceed
   real(t_p) e,fgrp,force1,force2
   parameter(&
      &grp=__use_groups__,&
      &tver=__use_ene__,&
      &tfea=__use_groups__+__use_polymer__+__use_gamd__&
      &)
!
!     zero out the energy component and partitioning terms
!
   if(deb_Path) write(*,*) 'estrbnd'
   eba    = 0
!
!     calculate the stretch-bend energy term
!
!$acc parallel loop present(eba,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz, &
!$acc     eW1aMD,deba,deaMD,isb,typeA,anat,bl,iang,sbk,strbndglob) &
!$acc     default(present) async &
!$acc     reduction(+:eba)
   do istrbndloc = 1, nstrbndloc
      istrbnd = strbndglob(istrbndloc)
      i       = isb(1,istrbnd)
#ifdef USE_NVSHMEM
      ipe     =     (i-1)/nangle_pe
      ind     = mod((i-1),nangle_pe) +1
      ia      = d_iang(ipe)%pel(1,ind)
      ib      = d_iang(ipe)%pel(2,ind)
      ic      = d_iang(ipe)%pel(3,ind)
#else
      ia      = iang(1,i)
      ib      = iang(2,i)
      ic      = iang(3,i)
#endif
      force1  = sbk(1,istrbnd)
      force2  = sbk(2,istrbnd)
      if (use_group.and.IAND(tfea,grp).NE.0)&
         &call groups3_inl(fgrp,ia,ib,ic,ngrp,grplist,wgrp)
!
!     decide whether to compute the current interaction
!
      proceed = (use(ia).or. use(ib).or. use(ic))
!
!     get the coordinates of the atoms in the angle
!
      if (proceed) then
         call ker_strbnd(i,ia,ib,ic,istrbnd,loc,isb,typeA,aMDwattype&
            &,use_polymer,use_group,use_amd_wat1&
            &,stbnunit,fgrp,force1,force2,anat,bl,x,y,z&
            &,eba,eW1aMD,e,deba,deaMD&
            &,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz&
            &,tver,tfea)
      end if
   end do
end

subroutine estrbnd
   use atoms
   use angle
   use bond
   use deriv
   use strbnd
   implicit none
   interface
      subroutine estrbnd_(isb,typeA,anat,bl,deba,deW1aMD)
         integer  ,intent(in):: isb(:,:),typeA(:)
         real(t_p),intent(in):: anat(:),bl(:)
         real(r_p),intent(out):: deba(1),deW1aMD(1)
      end subroutine
   end interface
   call estrbnd_(isb,type,anat,bl,deba,deW1aMD)
end subroutine
