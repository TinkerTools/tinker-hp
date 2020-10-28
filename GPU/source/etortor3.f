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
#include "tinker_precision.h"
      module etortor3_inl
        contains
#include "image.f.inc"
      end module

      subroutine etortor3
      use sizes
      use action
      use analyz
      use atmlst
      use atmtyp ,only: atomic
      use atoms
      use bitor
      use bound
      use domdec
      use energi
      use etortor3_inl
      use group
      use inform
      use iounit
      use ktrtor
      use math
      use tinheader ,only:ti_p,re_p
      use torpot
      use tortor
      use usage
      implicit none
      integer i,k,itortor,iitortor
      integer pos1,pos2
      integer ia,ib,ic,id,ie
      integer ibloc,icloc,idloc
      integer nlo,nhi,nt
      integer xlo,ylo
      real(t_p) e,sign
      real(t_p) angle1,angle2
      real(t_p) value1,value2
      real(t_p) cosine1,cosine2
      real(t_p) xt,yt,zt,rt2
      real(t_p) xu,yu,zu,ru2
      real(t_p) xv,yv,zv,rv2
      real(t_p) rtru
      real(t_p) rurv
      real(t_p) x1l,x1u
      real(t_p) y1l,y1u
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xie,yie,zie
      real(t_p) xba,yba,zba
      real(t_p) xdc,ydc,zdc
      real(t_p) xcb,ycb,zcb
      real(t_p) xed,yed,zed
      real(t_p) ftt(4),ft12(4)
      real(t_p) ft1(4),ft2(4)
      logical proceed
      logical header,huge
!$acc routine(chkttor) seq
!$acc routine(bcuint) seq
c
c
c     zero out the torsion-torsion energy and partitioning terms
c
      if (rank.eq.0.and.tinkerdebug) write(*,*) 'etortor3'
      nett   = 0
      ett    = 0.0_re_p
      aett   = 0.0_ti_p
      if (rank.eq.0) then
         header = .true.
      else
         header=.false.
      end if
c
c     calculate the torsion-torsion interaction energy term
c
!$acc parallel loop default(present) present(ett) async
!$acc&         private(ftt,ft12,ft1,ft2)
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
         icloc = loc(ic)
         ibloc = loc(ib)
         idloc = loc(id)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                       .or. use(id) .or. use(ie))
c
c     compute the value of the torsional angles
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xed = xie - xid
            yed = yie - yid
            zed = zie - zid
            if (use_polymer) then
               call image_inl (xba,yba,zba)
               call image_inl (xcb,ycb,zcb)
               call image_inl (xdc,ydc,zdc)
               call image_inl (xed,yed,zed)
            end if
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            xv = ydc*zed - yed*zdc
            yv = zdc*xed - zed*xdc
            zv = xdc*yed - xed*ydc
            rv2 = xv*xv + yv*yv + zv*zv
            rurv = sqrt(ru2 * rv2)
            if (rtru.ne.0.0_ti_p .and. rurv.ne.0.0_ti_p) then
               cosine1 = (xt*xu + yt*yu + zt*zu) / rtru
               cosine1 = min(1.0_ti_p,max(-1.0_ti_p,cosine1))
               angle1 = radian * acos(cosine1)
               sign = xba*xu + yba*yu + zba*zu
               if (sign .lt. 0.0_ti_p)  angle1 = -angle1
               value1 = angle1
               cosine2 = (xu*xv + yu*yv + zu*zv) / rurv
               cosine2 = min(1.0_ti_p,max(-1.0_ti_p,cosine2))
               angle2 = radian * acos(cosine2)
               sign = xcb*xv + ycb*yv + zcb*zv
               if (sign .lt. 0.0_ti_p)  angle2 = -angle2
               value2 = angle2
c
c     check for inverted chirality at the central atom
c
               call chkttor (ib,ic,id,sign,value1,value2,x,y,z,
     &                       type,atomic)
c
c     use bicubic interpolation to compute spline values
c
               nlo = 1
               nhi = tnx(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi+nlo) / 2
                  if (ttx(nt,k) .gt. value1) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               xlo = nlo
               nlo = 1
               nhi = tny(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi + nlo)/2
                  if (tty(nt,k) .gt. value2) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               ylo = nlo
               x1l = ttx(xlo,k)
               x1u = ttx(xlo+1,k)
               y1l = tty(ylo,k)
               y1u = tty(ylo+1,k)
               pos2 = ylo*tnx(k) + xlo
               pos1 = pos2 - tnx(k)
               ftt(1) = tbf(pos1,k)
               ftt(2) = tbf(pos1+1,k)
               ftt(3) = tbf(pos2+1,k)
               ftt(4) = tbf(pos2,k)
               ft1(1) = tbx(pos1,k)
               ft1(2) = tbx(pos1+1,k)
               ft1(3) = tbx(pos2+1,k)
               ft1(4) = tbx(pos2,k)
               ft2(1) = tby(pos1,k)
               ft2(2) = tby(pos1+1,k)
               ft2(3) = tby(pos2+1,k)
               ft2(4) = tby(pos2,k)
               ft12(1) = tbxy(pos1,k)
               ft12(2) = tbxy(pos1+1,k)
               ft12(3) = tbxy(pos2+1,k)
               ft12(4) = tbxy(pos2,k)
               call bcuint (ftt,ft1,ft2,ft12,x1l,x1u,
     &                      y1l,y1u,value1,value2,e)
               e = ttorunit * e
c
c     increment the total torsion-torsion energy
c
               nett = nett + 1
               ett  = ett + e
!$acc atomic
               aett(ibloc) = aett(ibloc) + e/3.0_ti_p
!$acc atomic
               aett(icloc) = aett(icloc) + e/3.0_ti_p
!$acc atomic
               aett(idloc) = aett(idloc) + e/3.0_ti_p
#ifndef _OPENACC
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 3.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Torsion-Torsion',
     &                         ' Interactions :',
     &                      //,' Type',17x,'Atom Numbers',16x,'Angle1',
     &                          4x,'Angle2',6x,'Energy',/)
                  end if
                  write (iout,20)  ia,ib,ic,id,ie,angle1,angle2,e
   20             format (' TorTor',2x,5i7,2x,2f10.4,f12.4)
               end if
#endif
            end if
         end if
      end do
!$acc update self(aett) async
      end
