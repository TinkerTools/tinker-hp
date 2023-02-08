c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine etortor  --  torsion-torsion cross term energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "etortor" calculates the torsion-torsion potential energy
c
c
#include "tinker_macro.h"
      module etortor_inl
#include "atomicOp.h.f"
        contains
#include "ker_tortor.inc.f"
      subroutine etortor_(i12,n12,typAtm,atomic,dett)
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
      real(r_p),intent(inout):: dett(3,*)

      integer   i,k,itortor,iitortor,ia,ib,ic,id,ie,ver,fea
      logical   proceed
      real(t_p) ftt(4),ft12(4),ft1(4),ft2(4),cw(4,4)
      real(t_p) fgrp
      parameter(
     &    ver=__use_ene__
     &   ,fea=__use_mpi__+__use_polymer__+__use_groups__
     &         )
c
c     calculate the torsion-torsion interaction energy term
c
!$acc parallel loop private(ft12,ftt,ft1,ft2,cw) async
!$acc&         present(tortorglob,ibitor,x,y,z,typAtm,atomic,i12,n12
!$acc&   ,grplist,wgrp,use,loc,itt,ttx,tty,tnx,tny,tbf,tbx,tby,tbxy)
!$acc&         present(ett,dett
!$acc&   ,g_vxx,g_vxy,g_vxz,g_vyy,g_vyz,g_vzz)
!$acc&         reduction(+:ett)
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
      end subroutine
      end module

      subroutine etortor
      use atoms
      use atmtyp
      use couple
      use deriv
      use etortor_inl
      call etortor_(i12,n12,type,atomic,dett)
      end subroutine
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine chkttor  --  check torsion-torsion chirality  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "chkttor" tests the attached atoms at a torsion-torsion central
c     site and inverts the angle values if the site is chiral
c
c     note that the sign convention used in this version is correct
c     for phi-psi torsion-torsion interactions as defined in the
c     AMOEBA protein force field; the code may need to be altered
c     for other uses of the torsion-torsion potential, and will not
c     correctly handle enantiomeric sugar rings in nucleic acids
c
c
      subroutine chkttor 
     &           (ib,ic,id,sign,value1,value2,x,y,z,type_,atomic)
!$acc routine seq
      use sizes
      use atoms  ,only:n
      use couple
      use tinheader
      implicit none
      integer i,j,k,m
      integer ia,ib,ic,id
      real(t_p),intent(in)::x(n),y(n),z(n)
      integer  ,intent(in)::type_(n),atomic(n)
      real(t_p),intent(out):: sign
      real(t_p),intent(inout):: value1
      real(t_p),intent(inout):: value2
      real(t_p) xac,yac,zac
      real(t_p) xbc,ybc,zbc
      real(t_p) xdc,ydc,zdc
      real(t_p) c1,c2,c3,vol
c
c
c
c     test for chirality at the central torsion-torsion site
c
      sign = 1.0_ti_p
      if (n12(ic) .eq. 4) then
         j = 0
         do i = 1, 4
            m = i12(i,ic)
            if (m.ne.ib .and. m.ne.id) then
               if (j .eq. 0) then
                  j = m
               else
                  k = m
               end if
            end if
         end do
         ia = 0
         if (type_(j) .gt. type_(k))  ia = j
         if (type_(k) .gt. type_(j))  ia = k
         if (atomic(j) .gt. atomic(k))  ia = j
         if (atomic(k) .gt. atomic(j))  ia = k
c
c     compute the signed parallelpiped volume at central site
c
         if (ia .ne. 0) then
            xac = x(ia) - x(ic)
            yac = y(ia) - y(ic)
            zac = z(ia) - z(ic)
            xbc = x(ib) - x(ic)
            ybc = y(ib) - y(ic)
            zbc = z(ib) - z(ic)
            xdc = x(id) - x(ic)
            ydc = y(id) - y(ic)
            zdc = z(id) - z(ic)
            c1 = ybc*zdc - zbc*ydc
            c2 = ydc*zac - zdc*yac
            c3 = yac*zbc - zac*ybc
            vol = xac*c1 + xbc*c2 + xdc*c3
c
c     invert the angle values if chirality has an inverted sign
c
            if (vol .lt. 0.0_ti_p) then
               sign = -1.0_ti_p
               value1 = -value1
               value2 = -value2
            end if
         end if
      end if
      end
