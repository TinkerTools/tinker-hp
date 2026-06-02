!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine eopbend1  --  out-of-plane energy and derivs  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "eopbend1" computes the out-of-plane bend potential energy and
!     first derivatives at trigonal centers via a Wilson-Decius-Cross
!     or Allinger angle
!
!
#include "tinker_macro.h"
module eopbend1gpu_inl
#include "atomicOp.h.f90"
contains
#include "ker_opbend.inc.f90"

   subroutine eopbend1gpu_(deopb)
      use angle
      use angpot
      use atmlst
      use atomsMirror
      use bound
      use domdec
      use energi
      use group
      use inform    ,only: deb_Path
      use math
      use opbend
      use tinheader ,only:ti_p,re_p
      use timestat
      use usage
      use virial
      implicit none
      real(r_p),intent(inout):: deopb(:,:)
      integer i,iopbend,iopbendloc
      integer ia,ib,ic,id
      integer tver,tfea
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,force,fgrp
      logical   proceed
      parameter(&
         &tver=__use_grd__+__use_ene__+__use_vir__,&
         &tfea=__use_groups__+__use_polymer__&
         &)
!
!     calculate the out-of-plane bending energy and derivatives
!
      if(deb_Path) write(*,*) 'eopbend1gpu'
      call timer_enter( timer_eopbend1 )

!$acc parallel loop async &
#ifdef USE_NVSHMEM_CUDA
!$acc      present(loc,use,opbendglob,opbk,iopb &
#else
!$acc      present(loc,use,opbendglob,iang,opbk,iopb &
#endif
!$acc     ,deopb,eopb,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz) &
!$acc      reduction(+:eopb,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do iopbendloc = 1, nopbendloc
         iopbend = opbendglob(iopbendloc)
         i       = iopb(iopbend)
#ifdef USE_NVSHMEM_CUDA
         ipe     =     (i-1)/nangle_pe
         ind     = mod((i-1),nangle_pe) +1
         ia      = d_iang(ipe)%pel(1,ind)
         ib      = d_iang(ipe)%pel(2,ind)
         ic      = d_iang(ipe)%pel(3,ind)
         id      = d_iang(ipe)%pel(4,ind)
#else
         ia      = iang(1,i)
         ib      = iang(2,i)
         ic      = iang(3,i)
         id      = iang(4,i)
#endif
         force   = opbk(iopbend)
!
!     decide whether to compute the current interaction
!
         if (use_group.and.IAND(tfea,__use_groups__).NE.0)&
            &call groups4_inl (fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
         proceed = (use(ia).or. use(ib).or. use(ic).or. use(id))
!
!     get the coordinates of the atoms at trigonal center
!
         if (proceed) then
            call ker_opbend(ia,ib,ic,id,opbtypInt,loc&
               &,use_group,use_polymer&
               &,opbunit,fgrp,force,copb,qopb,popb,sopb,x,y,z&
               &,eopb,e,deopb,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz&
               &,tver,tfea)
         end if
      end do
      call timer_exit( timer_eopbend1 )
   end
end module

subroutine eopbend1gpu
   use deriv
   use eopbend1gpu_inl
   implicit none
   call eopbend1gpu_(deopb)
end subroutine
