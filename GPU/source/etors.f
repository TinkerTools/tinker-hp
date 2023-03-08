c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine etors  --  torsional potential energy  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "etors" calculates the torsional potential energy
c
c
#include "tinker_macro.h"
      module etors_inl
#include "atomicOp.h.f"
        contains
#include "ker_tors.inc.f"
      end module

      subroutine etors
      implicit none
      call etors0a
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine etors0a  --  standard torsional energy  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "etors0a" calculates the torsional potential energy
c     using a standard sum of Fourier terms
c
c
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
      real(t_p),intent(in):: tors1(:,:),tors2(:,:),tors3(:,:),tors4(:,:)
     &         ,tors5(:,:),tors6(:,:)
      real(r_p),intent(inout):: det(1)
      integer itor,i,ia,ib,ic,id,tver,tfea,grp
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,fgrp
      logical proceed
      parameter(
     &      grp=__use_groups__,
     &     tver=__use_ene__,
     &     tfea=__use_groups__+__use_polymer__+__use_gamd__
     &         )
c
      if(deb_Path) write(*,*) 'etors'
c
c     zero out the torsional potential energy
c
      et = 0.0_ti_p
c
c     calculate the torsional angle energy term
c
!$acc parallel loop default(present) present(et) async
!$acc&         reduction(+:et)
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
         if (use_group.and.IAND(tfea,grp).NE.0)
     &      call groups4_inl (fgrp,ia,ib,ic,id,ngrp,grplist,wgrp)
c
c     decide whether to compute the current interaction
c
         proceed = (use(ia).or. use(ib).or. use(ic).or. use(id))
c
c     compute the value of the torsional angle
c
         if (proceed) then
            call ker_tors(i,ia,ib,ic,id,loc
     &              ,use_group,use_polymer
     &              ,torsunit,fgrp,tors1,tors2,tors3,tors4,tors5,tors6
     &              ,x,y,z
     &              ,et,e,det,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &              ,tver,tfea)
         end if
      end do
      end

      subroutine etors0a
      use deriv   ,only: det
      use tors
      implicit none
      interface
      subroutine etors0a_(tors1,tors2,tors3,tors4,tors5,tors6,det)
      real(t_p),intent(in):: tors1(:,:),tors2(:,:),tors3(:,:),tors4(:,:)
     &         ,tors5(:,:),tors6(:,:)
      real(r_p),intent(inout)::  det(1)
      end subroutine
      end interface
      call etors0a_(tors1,tors2,tors3,tors4,tors5,tors6,det)
      end subroutine
