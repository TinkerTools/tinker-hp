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
#include "tinker_precision.h"
      module etortor_inl
        contains
#include "image.f.inc"
      end module

      subroutine etortor
      use atmlst
      use atmtyp ,only: atomic
      use atoms
      use bitor
      use bound
      use energi
      use etortor_inl
      use group
      use ktrtor
      use math
      use tinheader
      use torpot
      use tortor
      use usage
      implicit none
      integer i,k,itortor,iitortor
      integer pos1,pos2
      integer ia,ib,ic,id,ie
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
!$acc routine(chkttor) seq
!$acc routine(bcuint) seq
c
c     zero out the torsion-torsion energy
c
      ett = 0.0_re_p
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
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                       .or. use(id) .or. use(ie))
c
c     compute the values of the torsional angles
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
               ett = ett + e
            end if
         end if
      end do
      end
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
