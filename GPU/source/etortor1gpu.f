c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine etortor1  --  torsion-torsion energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "etortor1" calculates the torsion-torsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
#if defined(_CUDA) && !defined(_OPENACC)
#define _CUDA_ONLY
#endif
      module etortor1gpu_inl
      use tinheader ,only: ti_p
#include "atomicOp.h.f"
        contains
#include "ker_tortor.inc.f"

      subroutine etortor1gpu_(i12,n12,typAtm,atomic,dett)
      use atmlst
      use atmtyp
      use atoms    ,only: n,x,y,z
      use bitor
      use bound
      use deriv    ,only: deamdD
      use domdec
      use energi
      use group
      use inform   ,only: deb_Path
      use ktrtor
      use mamd
      use math
      use potent   ,only: use_amd_dih
      use sizes
      use tinheader,only: ti_p,re_p
      use torpot
      use tortor   ,only: itt,ntortorloc
      use usage
      use virial
      implicit none
      integer  ,intent(in):: i12(:,:),n12(:)
     &         ,atomic(:),typAtm(:)
      real(r_p),intent(inout):: dett(:,:)

      integer   i,k,itortor,iitortor,ia,ib,ic,id,ie,ver,fea
      logical   proceed
      real(t_p) ftt(4),ft12(4),ft1(4),ft2(4),cw(4,4)
      real(t_p) fgrp
      parameter(
     &    ver=__use_grd__+__use_ene__+__use_vir__
     &   ,fea=__use_mpi__+__use_polymer__+__use_groups__
     &         )

      if (deb_Path) write(*,*) 'etortor1gpu'
c
c     calculate the torsion-torsion interaction energy term
c

!$acc parallel loop private(ft12,ftt,ft1,ft2,cw) async
!$acc&         present(tortorglob,ibitor,x,y,z,typAtm,atomic,i12,n12
!$acc&   ,grplist,wgrp,use,loc,itt,ttx,tty,tnx,tny,tbf,tbx,tby,tbxy)
!$acc&         present(ett,dett
!$acc&   ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         reduction(+:ett,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
      do itortor = 1, ntortorloc
         iitortor = tortorglob(itortor)
         i = itt(1,iitortor)
         k = itt(2,iitortor)
         if (itt(3,iitortor) .eq. 1) then
            ia = ibitor(1,i)
            ib = ibitor(2,i)
            ic = ibitor(3,i)
            id = ibitor(4,i)
            ie = ibitor(5,i)
         else
            ia = ibitor(5,i)
            ib = ibitor(4,i)
            ic = ibitor(3,i)
            id = ibitor(2,i)
            ie = ibitor(1,i)
         end if
c
c     decide whether to compute the current interaction
c
         if (use_group) 
     &      call groups5_inl(fgrp,ia,ib,ic,id,ie,ngrp,grplist,wgrp)
         proceed = merge(.true.,use(ia).or.use(ib).or.use(ic).or.
     &                          use(id).or.use(ie), useAll)
c
c     compute the values of the torsional angles
c
         if (proceed) then
            call ker_tortor
     &           (k,ia,ib,ic,id,ie,iitortor,ntortorloc,n,loc ! Input
     &           ,i12,n12,atomic,typAtm
     &           ,x,y,z,tnx,tny,ttx,tty,tbf,tbx,tby,tbxy,ttorunit,fgrp
     &           ,use_polymer,use_group,use_amd_dih,use_virial
     &           ,ett,dett            ! Output
     &           ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz
     &           ,ver,fea)
         end if
      end do
c
c     increment the aMD/GaMD derivs table
c
      if (use_amd_dih) then
!$acc parallel loop collapse(2) async
!$acc&         present(deamdD,dett)
         do k = 1,nloc; do i = 1,3
            deamdD(i,k) = deamdD(i,k) + dett(i,k)
         end do; end do
!$acc serial async present(eDaMD,ett)
         eDaMD = eDaMD + ett
!$acc end serial
      end if
      end
      end module

      subroutine etortor1gpu
      use atoms
      use atmtyp
      use couple
      use deriv
      use etortor1gpu_inl
      call etortor1gpu_(i12,n12,type,atomic,dett)
      end subroutine
