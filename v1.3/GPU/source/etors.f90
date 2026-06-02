!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ########################################################
!     ##                                                    ##
!     ##  subroutine etors  --  torsional potential energy  ##
!     ##                                                    ##
!     ########################################################
!
!
!     "etors" calculates the torsional potential energy
!
!
#include "tinker_macro.h"
module etors_inl
#include "atomicOp.h.f90"
contains
#include "ker_tors.inc.f90"
end module

subroutine etors
   implicit none
   call etors0a
end
!
!
!     #########################################################
!     ##                                                     ##
!     ##  subroutine etors0a  --  standard torsional energy  ##
!     ##                                                     ##
!     #########################################################
!
!
!     "etors0a" calculates the torsional potential energy
!     using a standard sum of Fourier terms
!
!
subroutine etors0a_(tors1,tors2,tors3,tors4,tors5,tors6,det)
   use atmlst
   use atoms
   use bound
   use domdec    ,only: loc
   use energi
   use etors_inl
   use group
   use inform
   use tinheader
   use torpot
   use tors      ,only: ntorsloc,itors
   use usage
   use virial
   implicit none
   real(t_p),intent(in):: tors1(:,:),tors2(:,:),tors3(:,:),tors4(:,:)&
      &,tors5(:,:),tors6(:,:)
   real(r_p),intent(inout):: det(1)
   integer itor,i,ia,ib,ic,id,tver,tfea,grp
#ifdef USE_NVSHMEM_CUDA
   integer ipe,ind
#endif
   real(t_p) e,fgrp
   logical proceed
   parameter(&
      &grp=__use_groups__,&
      &tver=__use_ene__,&
      &tfea=__use_groups__+__use_polymer__+__use_gamd__&
      &)
!
   if(deb_Path) write(*,*) 'etors'
!
!     zero out the torsional potential energy
!
   et = 0.0_ti_p
!
!     calculate the torsional angle energy term
!
!$acc parallel loop default(present) present(et) async &
!$acc         reduction(+:et)
   do itor = 1, ntorsloc
      i = torsglob(itor)
#ifdef USE_NVSHMEM_CUDA
      ipe = (i-1)/ntors_pe
      ind = mod((i-1),ntors_pe) +1
      ia  = d_itors(ipe)%pel(1,ind)
      ib  = d_itors(ipe)%pel(2,ind)
      ic  = d_itors(ipe)%pel(3,ind)
      id  = d_itors(ipe)%pel(4,ind)
#else
      ia  = itors(1,i)
      ib  = itors(2,i)
      ic  = itors(3,i)
      id  = itors(4,i)
#endif
      if (use_group.and.IAND(tfea,grp).NE.0)&
         &call groups4_inl (fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
!
!     decide whether to compute the current interaction
!
      proceed = (use(ia).or. use(ib).or. use(ic).or. use(id))
!
!     compute the value of the torsional angle
!
      if (proceed) then
         call ker_tors(i,ia,ib,ic,id,loc&
            &,use_group,use_polymer&
            &,torsunit,fgrp,tors1,tors2,tors3,tors4,tors5,tors6&
            &,x,y,z&
            &,et,e,det,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz&
            &,tver,tfea)
      end if
   end do
end

subroutine etors0a
   use deriv   ,only: det
   use tors
   implicit none
   interface
      subroutine etors0a_(tors1,tors2,tors3,tors4,tors5,tors6,det)
         real(t_p),intent(in):: tors1(:,:),tors2(:,:),tors3(:,:),tors4(:,:)&
            &,tors5(:,:),tors6(:,:)
         real(r_p),intent(inout)::  det(1)
      end subroutine
   end interface
   call etors0a_(tors1,tors2,tors3,tors4,tors5,tors6,det)
end subroutine
