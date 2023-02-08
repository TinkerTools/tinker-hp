c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine etors1  --  torsional energy & derivatives  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "etors1" calculates the torsional potential energy and first
c     derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      module etors1gpu_inl
#include "atomicOp.h.f"
        contains
#include "ker_tors.inc.f"
      end module

      subroutine etors1gpu
      implicit none
      call etors1agpu
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine etors1a  --  standard torsional energy & derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "etors1a" calculates the torsional potential energy and first
c     derivatives with respect to Cartesian coordinates using a
c     standard sum of Fourier terms
c
c
      subroutine etors1agpu_(tors1,tors2,tors3,tors4,tors5,tors6,det)
      use atmlst    ,only: torsglob
      use atoms
      use bound
      use deriv     ,only: deamdD
      use domdec
      use energi    ,only: et,eDaMD
      use etors1gpu_inl
      use group
      use inform    ,only: deb_Path,minmaxone
      use mamd
      use potent    ,only: use_amd_dih
      use tinheader ,only: ti_p,re_p
      use torpot
      use tors      ,only: ntorsloc,itors
      use usage
      use utilgpu   ,only: mem_move,rec_stream
      use tinMemory ,only: mipk
      use virial
      implicit none
      real(t_p),intent(in):: tors1(:,:),tors2(:,:),tors3(:,:),tors4(:,:)
     &         ,tors5(:,:),tors6(:,:)
      real(r_p),intent(inout):: det(:,:)
      integer itor,i,ia,ib,ic,id,tver,tfea,grp
      integer(mipk) siz
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) e,fgrp
      logical proceed
      parameter(
     &      grp=__use_groups__,
     &     tver=__use_grd__+__use_ene__+__use_vir__,
     &     tfea=__use_groups__+__use_polymer__+__use_gamd__
     &         )

      if (deb_Path) write(*,*) 'etors1gpu'
c
c     calculate the torsional angle energy and first derivatives
c
!$acc parallel loop async
!$acc&  present(det,et,torsglob,loc,x,y,z,use,
#ifndef USE_NVSHMEM_CUDA
!$acc&     itors,tors1,tors2,tors3,tors4,tors5,tors6,
#endif
!$acc&     g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&  reduction(+:et,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do itor = 1, ntorsloc
         i   = torsglob(itor)
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
c
c     increment the aMD/GaMD derivs table
c
      if (use_amd_dih) then
         siz = 3*nloc
         call mem_move(deamdD,det,siz,rec_stream)
!$acc serial async present(eDaMD,et)
         eDaMD = et
!$acc end serial
      end if
      end

      subroutine etors1agpu
      use deriv   ,only: det
      use tors
      implicit none
      interface
      subroutine etors1agpu_(tors1,tors2,tors3,tors4,tors5,tors6,det)
      real(t_p),intent(in):: tors1(:,:),tors2(:,:),tors3(:,:),tors4(:,:)
     &         ,tors5(:,:),tors6(:,:)
      real(r_p),intent(inout):: det(:,:)
      end subroutine
      end interface
      call etors1agpu_(tors1,tors2,tors3,tors4,tors5,tors6,det)
      end subroutine
