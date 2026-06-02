!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##########################################################
!     ##                                                      ##
!     ##  subroutine etors3  --  torsional energy & analysis  ##
!     ##                                                      ##
!     ##########################################################
!
!
!     "etors3" calculates the torsional potential energy; also
!     partitions the energy among the atoms
!
!
#include "tinker_macro.h"
module etors3_inl
#include "atomicOp.h.f90"
contains
#include "ker_tors.inc.f90"
!
!
!     ###########################################################
!     ##                                                       ##
!     ##  subroutine etors3a  --  standard torsional analysis  ##
!     ##                                                       ##
!     ###########################################################
!
!
!     "etors3a" calculates the torsional potential energy using
!     a standard sum of Fourier terms and partitions the energy
!     among the atoms
!
!
   subroutine etors3a_(tors1,tors2,tors3,tors4,tors5,tors6,det)
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use domdec
      use energi
      use group
      use inform
      use iounit
      use math
      use tinheader ,only:ti_p,re_p
      use torpot
      use tors      ,only: ntorsloc,itors
      use usage
      use virial
      implicit none
      real(t_p),intent(in):: tors1(:,:),tors2(:,:),tors3(:,:),tors4(:,:)&
         &,tors5(:,:),tors6(:,:)
      real(r_p),intent(inout):: det(:,:)
      integer itor,i,ia,ib,ic,id,tver,tfea,grp
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,fgrp
      logical proceed,header,huge
      parameter(&
         &grp=__use_groups__,&
         &tver=__use_ene__+__use_act__,&
         &tfea=__use_groups__+__use_polymer__+__use_gamd__&
         &)

      if(deb_Path) write(*,*) 'etors3'
!
!     zero out the torsional energy and partitioning terms
!
      net = 0
      et  = 0
#ifndef _OPENACC
      aet = 0
#endif
      header = (rank.eq.0)
!
!     calculate the torsional angle energy term
!
!$acc parallel loop default(present) present(et) async &
!$acc         reduction(+:et,net)
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
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
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
            if(e.ne.0.0) net = net + 1
#if 0
            ib = loc(ib)
            ic = loc(ic)
            aet(ib) = aet(ib) + 0.5*e
            aet(ic) = aet(ic) + 0.5*e
!
!     print a message if the energy of this interaction is large
!
            huge = (e .gt. 3.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
10                format (/,' Individual Torsional Angle',&
                     &' Interactions :',&
                     &//,' Type',25x,'Atom Names',21x,'Angle',&
                     &6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ic,&
                  &name(ic),id,name(id),angle,e
20             format (' Torsion',3x,4(i7,'-',a3),f11.4,f12.4)
            end if
#endif
         end if
      end do
   end
end module

subroutine etors3
   use deriv   ,only: det
   use etors3_inl
   use tors
   use utilgpu ,only: lam_buff
   implicit  none
   real(r_p),pointer:: buff(:,:)
   buff(1:3,1:1) => lam_buff(1:3)
   call etors3a_(tors1,tors2,tors3,tors4,tors5,tors6,buff)
end subroutine
