c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine etortor3  --  torsion-torsion energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "etortor3" calculates the torsion-torsion potential energy;
c     also partitions the energy terms among the atoms
c
c
#include "tinker_macro.h"
      module etortor3_inl
#include "atomicOp.h.f"
        contains
#include "ker_tortor.inc.f"

      subroutine etortor3_(i12,n12,typAtm,atomic,dett)
      use action   ,only: nett
      use analyz   ,only: aett
      use atoms    ,only: n,x,y,z
      use atmlst   ,only: tortorglob
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
      real(r_p),intent(inout):: dett(3,*)

      integer   i,k,itortor,iitortor,ia,ib,ic,id,ie,ver,fea
      logical   proceed,header
      real(t_p) ftt(4),ft12(4),ft1(4),ft2(4),cw(4,4)
      real(t_p) fgrp
      parameter(
     &    ver=__use_ene__+__use_act__
     &   ,fea=__use_mpi__+__use_polymer__+__use_groups__
     &         )
      nett   = 0
      header = (rank.eq.0)
c
c     calculate the torsion-torsion interaction energy term
c
!$acc parallel loop private(ft12,ftt,ft1,ft2,cw) async
!$acc&         present(tortorglob,ibitor,x,y,z,typAtm,atomic,i12,n12
!$acc&   ,grplist,wgrp,use,loc,itt,ttx,tty,tnx,tny,tbf,tbx,tby,tbxy)
!$acc&         present(ett,dett
!$acc&   ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz) copy(nett)
!$acc&         reduction(+:ett,nett)
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
            nett = nett + 1
#if 0
!$acc atomic
            aett(ibloc) = aett(ibloc) + e/3.0_ti_p
!$acc atomic
            aett(icloc) = aett(icloc) + e/3.0_ti_p
!$acc atomic
            aett(idloc) = aett(idloc) + e/3.0_ti_p
c
c     print a message if the energy of this interaction is large
c
            huge = (e .gt. 3.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Torsion-Torsion',
     &                      ' Interactions :',
     &                   //,' Type',17x,'Atom Numbers',16x,'Angle1',
     &                       4x,'Angle2',6x,'Energy',/)
               end if
               write (iout,20)  ia,ib,ic,id,ie,angle1,angle2,e
   20          format (' TorTor',2x,5i7,2x,2f10.4,f12.4)
            end if
#endif
         end if
      end do
      end subroutine
      end module

      subroutine etortor3
      use atoms
      use atmtyp
      use couple
      use domdec
      use deriv
      use etortor3_inl
      use inform

      if (deb_Path) write(*,*) 'etortor3'
      call etortor3_(i12,n12,type,atomic,dett)
      end
