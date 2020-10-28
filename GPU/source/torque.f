c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_precision.h"
      subroutine torque (i,trq,frcx,frcy,frcz,de)
      use atoms
      use deriv
      use domdec
      use mpole
      use sizes
      use tinheader ,only: ti_p,re_p
      implicit none
      integer i,j
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      real(t_p) du,dv,dw,dot
      real(t_p) usiz,vsiz,wsiz
      real(t_p) psiz,rsiz,ssiz
      real(t_p) t1siz,t2siz
      real(t_p) uvsiz,uwsiz,vwsiz
      real(t_p) ursiz,ussiz
      real(t_p) vssiz,wssiz
      real(t_p) delsiz,dphiddel
      real(t_p) uvcos,uwcos,vwcos
      real(t_p) urcos,uscos
      real(t_p) vscos,wscos
      real(t_p) rwcos,wpcos
      real(t_p) upcos,vpcos
      real(t_p) rucos,rvcos
      real(t_p) ut1cos,ut2cos
      real(t_p) uvsin,uwsin,vwsin
      real(t_p) ursin,ussin
      real(t_p) vssin,wssin
      real(t_p) rwsin,wpsin
      real(t_p) upsin,vpsin
      real(t_p) rusin,rvsin
      real(t_p) ut1sin,ut2sin
      real(t_p) dphidu,dphidv,dphidw
      real(t_p) dphidr,dphids
      real(t_p) trq(3),frcx(3)
      real(t_p) frcy(3),frcz(3)
      real(t_p) u(3),v(3),w(3)
      real(t_p) p(3),r(3),s(3)
      real(t_p) t1(3),t2(3)
      real(t_p) uv(3),uw(3),vw(3)
      real(t_p) ur(3),us(3)
      real(t_p) vs(3),ws(3)
      real(t_p) del(3),eps(3)
      real(t_p) de(3,*)
      character*8 axetyp

c        if(rank.eq.0.and.tinkerdebug) write(*,*) 'torque'

c
c
c     zero out force components on local frame-defining atoms
c
      do j = 1, 3
         frcz(j) = 0.0_ti_p
         frcx(j) = 0.0_ti_p
         frcy(j) = 0.0_ti_p
      end do
c
c     get the local frame type and the frame-defining atoms
c
      ia = zaxis(i)
      if (ia.gt.0) ialoc = loc(ia)
      ib = ipole(i)
      ibloc = loc(ib)
      ic = xaxis(i)
      if (ic.gt.0) icloc = loc(ic)
      id = yaxis(i)
      if (id.gt.0) idloc = loc(id)
      axetyp = polaxe(i)
      if (axetyp .eq. 'None')  return
c
c     construct the three rotation axes for the local frame
c
      u(1) = x(ia) - x(ib)
      u(2) = y(ia) - y(ib)
      u(3) = z(ia) - z(ib)
      usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
      if (axetyp .ne. 'Z-Only') then
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
         vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      else
         v(1) = 1.0_ti_p
         v(2) = 0.0_ti_p
         v(3) = 0.0_ti_p
         vsiz = 1.0_ti_p
         dot = u(1) / usiz
         if (abs(dot) .gt. 0.866_ti_p) then
            v(1) = 0.0_ti_p
            v(2) = 1.0_ti_p
         end if
      end if
      if (axetyp.eq.'Z-Bisect' .or. axetyp.eq.'3-Fold') then
         w(1) = x(id) - x(ib)
         w(2) = y(id) - y(ib)
         w(3) = z(id) - z(ib)
      else
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
      end if
      wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
      do j = 1, 3
         u(j) = u(j) / usiz
         v(j) = v(j) / vsiz
         w(j) = w(j) / wsiz
      end do
c
c     build some additional axes needed for the Z-Bisect method
c
      if (axetyp .eq. 'Z-Bisect') then
         r(1) = v(1) + w(1)
         r(2) = v(2) + w(2)
         r(3) = v(3) + w(3)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         s(1) = u(2)*r(3) - u(3)*r(2)
         s(2) = u(3)*r(1) - u(1)*r(3)
         s(3) = u(1)*r(2) - u(2)*r(1)
         ssiz = sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
            s(j) = s(j) / ssiz
         end do
      end if
c
c     find the perpendicular and angle for each pair of axes
c
      uv(1) = v(2)*u(3) - v(3)*u(2)
      uv(2) = v(3)*u(1) - v(1)*u(3)
      uv(3) = v(1)*u(2) - v(2)*u(1)
      uvsiz = sqrt(uv(1)*uv(1) + uv(2)*uv(2) + uv(3)*uv(3))
      uw(1) = w(2)*u(3) - w(3)*u(2)
      uw(2) = w(3)*u(1) - w(1)*u(3)
      uw(3) = w(1)*u(2) - w(2)*u(1)
      uwsiz = sqrt(uw(1)*uw(1) + uw(2)*uw(2) + uw(3)*uw(3))
      vw(1) = w(2)*v(3) - w(3)*v(2)
      vw(2) = w(3)*v(1) - w(1)*v(3)
      vw(3) = w(1)*v(2) - w(2)*v(1)
      vwsiz = sqrt(vw(1)*vw(1) + vw(2)*vw(2) + vw(3)*vw(3))
      do j = 1, 3
         uv(j) = uv(j) / uvsiz
         uw(j) = uw(j) / uwsiz
         vw(j) = vw(j) / vwsiz
      end do
      if (axetyp .eq. 'Z-Bisect') then
         ur(1) = r(2)*u(3) - r(3)*u(2)
         ur(2) = r(3)*u(1) - r(1)*u(3)
         ur(3) = r(1)*u(2) - r(2)*u(1)
         ursiz = sqrt(ur(1)*ur(1) + ur(2)*ur(2) + ur(3)*ur(3))
         us(1) = s(2)*u(3) - s(3)*u(2)
         us(2) = s(3)*u(1) - s(1)*u(3)
         us(3) = s(1)*u(2) - s(2)*u(1)
         ussiz = sqrt(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
         vs(1) = s(2)*v(3) - s(3)*v(2)
         vs(2) = s(3)*v(1) - s(1)*v(3)
         vs(3) = s(1)*v(2) - s(2)*v(1)
         vssiz = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
         ws(1) = s(2)*w(3) - s(3)*w(2)
         ws(2) = s(3)*w(1) - s(1)*w(3)
         ws(3) = s(1)*w(2) - s(2)*w(1)
         wssiz = sqrt(ws(1)*ws(1) + ws(2)*ws(2) + ws(3)*ws(3))
         do j = 1, 3
            ur(j) = ur(j) / ursiz
            us(j) = us(j) / ussiz
            vs(j) = vs(j) / vssiz
            ws(j) = ws(j) / wssiz
         end do
      end if
c
c     get sine and cosine of angles between the rotation axes
c
      uvcos = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      uvsin = sqrt(1.0_ti_p - uvcos*uvcos)
      uwcos = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
      uwsin = sqrt(1.0_ti_p - uwcos*uwcos)
      vwcos = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
      vwsin = sqrt(1.0_ti_p - vwcos*vwcos)
      if (axetyp .eq. 'Z-Bisect') then
         urcos = u(1)*r(1) + u(2)*r(2) + u(3)*r(3)
         ursin = sqrt(1.0_ti_p - urcos*urcos)
         uscos = u(1)*s(1) + u(2)*s(2) + u(3)*s(3)
         ussin = sqrt(1.0_ti_p - uscos*uscos)
         vscos = v(1)*s(1) + v(2)*s(2) + v(3)*s(3)
         vssin = sqrt(1.0_ti_p - vscos*vscos)
         wscos = w(1)*s(1) + w(2)*s(2) + w(3)*s(3)
         wssin = sqrt(1.0_ti_p - wscos*wscos)
      end if
c
c     compute the projection of v and w onto the ru-plane
c
      if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
            t1(j) = v(j) - s(j)*vscos
            t2(j) = w(j) - s(j)*wscos
         end do
         t1siz = sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
         t2siz = sqrt(t2(1)*t2(1)+t2(2)*t2(2)+t2(3)*t2(3))
         do j = 1, 3
            t1(j) = t1(j) / t1siz
            t2(j) = t2(j) / t2siz
         end do
         ut1cos = u(1)*t1(1) + u(2)*t1(2) + u(3)*t1(3)
         ut1sin = sqrt(1.0_ti_p - ut1cos*ut1cos)
         ut2cos = u(1)*t2(1) + u(2)*t2(2) + u(3)*t2(3)
         ut2sin = sqrt(1.0_ti_p - ut2cos*ut2cos)
      end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      dphidu = -trq(1)*u(1) - trq(2)*u(2) - trq(3)*u(3)
      dphidv = -trq(1)*v(1) - trq(2)*v(2) - trq(3)*v(3)
      dphidw = -trq(1)*w(1) - trq(2)*w(2) - trq(3)*w(3)
      if (axetyp .eq. 'Z-Bisect') then
         dphidr = -trq(1)*r(1) - trq(2)*r(2) - trq(3)*r(3)
         dphids = -trq(1)*s(1) - trq(2)*s(2) - trq(3)*s(3)
      end if
c
c     force distribution for the Z-Only local coordinate method
c
      if (axetyp .eq. 'Z-Only') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            de(j,ialoc) = de(j,ialoc) + du
            de(j,ibloc) = de(j,ibloc) - du
            frcz(j) = frcz(j) + du
         end do
c
c     force distribution for the Z-then-X local coordinate method
c
      else if (axetyp .eq. 'Z-then-X') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin)
            de(j,ialoc) = de(j,ialoc) + du
            de(j,icloc) = de(j,icloc) + dv
            de(j,ibloc) = de(j,ibloc) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Bisector local coordinate method
c
      else if (axetyp .eq. 'Bisector') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + 0.5_ti_p*uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5_ti_p*vw(j)*dphidw/vsiz
            de(j,ialoc) = de(j,ialoc) + du
            de(j,icloc) = de(j,icloc) + dv
            de(j,ibloc) = de(j,ibloc) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Z-Bisect local coordinate method
c
      else if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
            du = ur(j)*dphidr/(usiz*ursin) + us(j)*dphids/usiz
            dv = (vssin*s(j)-vscos*t1(j))*dphidu
     &              / (vsiz*(ut1sin+ut2sin))
            dw = (wssin*s(j)-wscos*t2(j))*dphidu
     &              / (wsiz*(ut1sin+ut2sin))
            de(j,ialoc) = de(j,ialoc) + du
            de(j,icloc) = de(j,icloc) + dv
            de(j,idloc) = de(j,idloc) + dw
            de(j,ibloc) = de(j,ibloc) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
c
c     force distribution for the 3-Fold local coordinate method
c
      else if (axetyp .eq. '3-Fold') then
         p(1) = u(1) + v(1) + w(1)
         p(2) = u(2) + v(2) + w(2)
         p(3) = u(3) + v(3) + w(3)
         psiz = sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3)) 
         do j = 1, 3 
            p(j) = p(j) / psiz
         end do
         wpcos = w(1)*p(1) + w(2)*p(2) + w(3)*p(3)
         upcos = u(1)*p(1) + u(2)*p(2) + u(3)*p(3)
         vpcos = v(1)*p(1) + v(2)*p(2) + v(3)*p(3)
         wpsin = sqrt(1.0_ti_p - wpcos*wpcos)
         upsin = sqrt(1.0_ti_p - upcos*upcos)
         vpsin = sqrt(1.0_ti_p - vpcos*vpcos)
         r(1) = u(1) + v(1)
         r(2) = u(2) + v(2)
         r(3) = u(3) + v(3)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
         end do
         rwcos = r(1)*w(1) + r(2)*w(2) + r(3)*w(3)
         rwsin = sqrt(1.0_ti_p - rwcos*rwcos)
         dphidr = -trq(1)*r(1) - trq(2)*r(2) - trq(3)*r(3)
         del(1) = r(2)*w(3) - r(3)*w(2)
         del(2) = r(3)*w(1) - r(1)*w(3) 
         del(3) = r(1)*w(2) - r(2)*w(1)
         delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))
         do j = 1, 3
            del(j) = del(j) / delsiz
         end do
         dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
         eps(1) = del(2)*w(3) - del(3)*w(2)
         eps(2) = del(3)*w(1) - del(1)*w(3)
         eps(3) = del(1)*w(2) - del(2)*w(1)
         do j = 1, 3
            dw = del(j)*dphidr/(wsiz*rwsin)
     &              + eps(j)*dphiddel*wpcos/(wsiz*psiz) 
            de(j,idloc) = de(j,idloc) + dw
            de(j,ibloc) = de(j,ibloc) - dw
            frcy(j) = frcy(j) + dw
         end do
         r(1) = v(1) + w(1)
         r(2) = v(2) + w(2)
         r(3) = v(3) + w(3)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
         end do
         rucos = r(1)*u(1) + r(2)*u(2) + r(3)*u(3)
         rusin = sqrt(1.0_ti_p - rucos*rucos) 
         dphidr = -trq(1)*r(1) - trq(2)*r(2) - trq(3)*r(3)
         del(1) = r(2)*u(3) - r(3)*u(2)
         del(2) = r(3)*u(1) - r(1)*u(3)
         del(3) = r(1)*u(2) - r(2)*u(1)
         delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))
         do j = 1, 3
            del(j) = del(j) / delsiz
         end do
         dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
         eps(1) = del(2)*u(3) - del(3)*u(2)
         eps(2) = del(3)*u(1) - del(1)*u(3)
         eps(3) = del(1)*u(2) - del(2)*u(1)
         do j = 1, 3
            du = del(j)*dphidr/(usiz*rusin)
     &              + eps(j)*dphiddel*upcos/(usiz*psiz)
            de(j,ialoc) = de(j,ialoc) + du
            de(j,ibloc) = de(j,ibloc) - du
            frcz(j) = frcz(j) + du
         end do
         r(1) = u(1) + w(1)
         r(2) = u(2) + w(2)
         r(3) = u(3) + w(3)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
         end do
         rvcos = r(1)*v(1) + r(2)*v(2) + r(3)*v(3) 
         rvsin = sqrt(1.0_ti_p - rvcos*rvcos)
         dphidr = -trq(1)*r(1) - trq(2)*r(2) - trq(3)*r(3)
         del(1) = r(2)*v(3) - r(3)*v(2)
         del(2) = r(3)*v(1) - r(1)*v(3)
         del(3) = r(1)*v(2) - r(2)*v(1)
         delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))
         do j = 1, 3
            del(j) = del(j) / delsiz
         end do
         dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
         eps(1) = del(2)*v(3) - del(3)*v(2)
         eps(2) = del(3)*v(1) - del(1)*v(3)
         eps(3) = del(1)*v(2) - del(2)*v(1)
         do j = 1, 3
            dv = del(j)*dphidr/(vsiz*rvsin)
     &              + eps(j)*dphiddel*vpcos/(vsiz*psiz)
            de(j,icloc) = de(j,icloc) + dv
            de(j,ibloc) = de(j,ibloc) - dv
            frcx(j) = frcx(j) + dv
         end do
      end if
      return
      end
c
      subroutine torque_rec (i,trq,frcx,frcy,frcz,de)
      use atoms
      use deriv
      use domdec
      use mpole
      use tinheader ,only: ti_p,re_p
      implicit none
      integer i,j
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      real(t_p) du,dv,dw,dot
      real(t_p) usiz,vsiz,wsiz
      real(t_p) psiz,rsiz,ssiz
      real(t_p) t1siz,t2siz
      real(t_p) uvsiz,uwsiz,vwsiz
      real(t_p) ursiz,ussiz
      real(t_p) vssiz,wssiz
      real(t_p) delsiz,dphiddel
      real(t_p) uvcos,uwcos,vwcos
      real(t_p) urcos,uscos
      real(t_p) vscos,wscos
      real(t_p) rwcos,wpcos
      real(t_p) upcos,vpcos
      real(t_p) rucos,rvcos
      real(t_p) ut1cos,ut2cos
      real(t_p) uvsin,uwsin,vwsin
      real(t_p) ursin,ussin
      real(t_p) vssin,wssin
      real(t_p) rwsin,wpsin
      real(t_p) upsin,vpsin
      real(t_p) rusin,rvsin
      real(t_p) ut1sin,ut2sin
      real(t_p) dphidu,dphidv,dphidw
      real(t_p) dphidr,dphids
      real(t_p) trq(3),frcx(3)
      real(t_p) frcy(3),frcz(3)
      real(t_p) u(3),v(3),w(3)
      real(t_p) p(3),r(3),s(3)
      real(t_p) t1(3),t2(3)
      real(t_p) uv(3),uw(3),vw(3)
      real(t_p) ur(3),us(3)
      real(t_p) vs(3),ws(3)
      real(t_p) del(3),eps(3)
      real(t_p) de(3,*)
      character*8 axetyp
c
c
c     zero out force components on local frame-defining atoms
c
      do j = 1, 3
         frcz(j) = 0.0_ti_p
         frcx(j) = 0.0_ti_p
         frcy(j) = 0.0_ti_p
      end do
c
c     get the local frame type and the frame-defining atoms
c
      ia = zaxis(i)
      if (ia.gt.0) ialoc = locrec1(ia)
      ib = ipole(i)
      ibloc = locrec1(ib)
      ic = xaxis(i)
      if (ic.gt.0) icloc = locrec1(ic)
      id = yaxis(i)
      if (id.gt.0) idloc = locrec1(id)
      axetyp = polaxe(i)
      if (axetyp .eq. 'None')  return
c
c     construct the three rotation axes for the local frame
c
      u(1) = x(ia) - x(ib)
      u(2) = y(ia) - y(ib)
      u(3) = z(ia) - z(ib)
      usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
      if (axetyp .ne. 'Z-Only') then
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
         vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      else
         v(1) = 1.0_ti_p
         v(2) = 0.0_ti_p
         v(3) = 0.0_ti_p
         vsiz = 1.0_ti_p
         dot = u(1) / usiz
         if (abs(dot) .gt. 0.866_ti_p) then
            v(1) = 0.0_ti_p
            v(2) = 1.0_ti_p
         end if
      end if
      if (axetyp.eq.'Z-Bisect' .or. axetyp.eq.'3-Fold') then
         w(1) = x(id) - x(ib)
         w(2) = y(id) - y(ib)
         w(3) = z(id) - z(ib)
      else
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
      end if
      wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
      do j = 1, 3
         u(j) = u(j) / usiz
         v(j) = v(j) / vsiz
         w(j) = w(j) / wsiz
      end do
c
c     build some additional axes needed for the Z-Bisect method
c
      if (axetyp .eq. 'Z-Bisect') then
         r(1) = v(1) + w(1)
         r(2) = v(2) + w(2)
         r(3) = v(3) + w(3)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         s(1) = u(2)*r(3) - u(3)*r(2)
         s(2) = u(3)*r(1) - u(1)*r(3)
         s(3) = u(1)*r(2) - u(2)*r(1)
         ssiz = sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
            s(j) = s(j) / ssiz
         end do
      end if
c
c     find the perpendicular and angle for each pair of axes
c
      uv(1) = v(2)*u(3) - v(3)*u(2)
      uv(2) = v(3)*u(1) - v(1)*u(3)
      uv(3) = v(1)*u(2) - v(2)*u(1)
      uvsiz = sqrt(uv(1)*uv(1) + uv(2)*uv(2) + uv(3)*uv(3))
      uw(1) = w(2)*u(3) - w(3)*u(2)
      uw(2) = w(3)*u(1) - w(1)*u(3)
      uw(3) = w(1)*u(2) - w(2)*u(1)
      uwsiz = sqrt(uw(1)*uw(1) + uw(2)*uw(2) + uw(3)*uw(3))
      vw(1) = w(2)*v(3) - w(3)*v(2)
      vw(2) = w(3)*v(1) - w(1)*v(3)
      vw(3) = w(1)*v(2) - w(2)*v(1)
      vwsiz = sqrt(vw(1)*vw(1) + vw(2)*vw(2) + vw(3)*vw(3))
      do j = 1, 3
         uv(j) = uv(j) / uvsiz
         uw(j) = uw(j) / uwsiz
         vw(j) = vw(j) / vwsiz
      end do
      if (axetyp .eq. 'Z-Bisect') then
         ur(1) = r(2)*u(3) - r(3)*u(2)
         ur(2) = r(3)*u(1) - r(1)*u(3)
         ur(3) = r(1)*u(2) - r(2)*u(1)
         ursiz = sqrt(ur(1)*ur(1) + ur(2)*ur(2) + ur(3)*ur(3))
         us(1) = s(2)*u(3) - s(3)*u(2)
         us(2) = s(3)*u(1) - s(1)*u(3)
         us(3) = s(1)*u(2) - s(2)*u(1)
         ussiz = sqrt(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
         vs(1) = s(2)*v(3) - s(3)*v(2)
         vs(2) = s(3)*v(1) - s(1)*v(3)
         vs(3) = s(1)*v(2) - s(2)*v(1)
         vssiz = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
         ws(1) = s(2)*w(3) - s(3)*w(2)
         ws(2) = s(3)*w(1) - s(1)*w(3)
         ws(3) = s(1)*w(2) - s(2)*w(1)
         wssiz = sqrt(ws(1)*ws(1) + ws(2)*ws(2) + ws(3)*ws(3))
         do j = 1, 3
            ur(j) = ur(j) / ursiz
            us(j) = us(j) / ussiz
            vs(j) = vs(j) / vssiz
            ws(j) = ws(j) / wssiz
         end do
      end if
c
c     get sine and cosine of angles between the rotation axes
c
      uvcos = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      uvsin = sqrt(1.0_ti_p - uvcos*uvcos)
      uwcos = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
      uwsin = sqrt(1.0_ti_p - uwcos*uwcos)
      vwcos = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
      vwsin = sqrt(1.0_ti_p - vwcos*vwcos)
      if (axetyp .eq. 'Z-Bisect') then
         urcos = u(1)*r(1) + u(2)*r(2) + u(3)*r(3)
         ursin = sqrt(1.0_ti_p - urcos*urcos)
         uscos = u(1)*s(1) + u(2)*s(2) + u(3)*s(3)
         ussin = sqrt(1.0_ti_p - uscos*uscos)
         vscos = v(1)*s(1) + v(2)*s(2) + v(3)*s(3)
         vssin = sqrt(1.0_ti_p - vscos*vscos)
         wscos = w(1)*s(1) + w(2)*s(2) + w(3)*s(3)
         wssin = sqrt(1.0_ti_p - wscos*wscos)
      end if
c
c     compute the projection of v and w onto the ru-plane
c
      if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
            t1(j) = v(j) - s(j)*vscos
            t2(j) = w(j) - s(j)*wscos
         end do
         t1siz = sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
         t2siz = sqrt(t2(1)*t2(1)+t2(2)*t2(2)+t2(3)*t2(3))
         do j = 1, 3
            t1(j) = t1(j) / t1siz
            t2(j) = t2(j) / t2siz
         end do
         ut1cos = u(1)*t1(1) + u(2)*t1(2) + u(3)*t1(3)
         ut1sin = sqrt(1.0_ti_p - ut1cos*ut1cos)
         ut2cos = u(1)*t2(1) + u(2)*t2(2) + u(3)*t2(3)
         ut2sin = sqrt(1.0_ti_p - ut2cos*ut2cos)
      end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      dphidu = -trq(1)*u(1) - trq(2)*u(2) - trq(3)*u(3)
      dphidv = -trq(1)*v(1) - trq(2)*v(2) - trq(3)*v(3)
      dphidw = -trq(1)*w(1) - trq(2)*w(2) - trq(3)*w(3)
      if (axetyp .eq. 'Z-Bisect') then
         dphidr = -trq(1)*r(1) - trq(2)*r(2) - trq(3)*r(3)
         dphids = -trq(1)*s(1) - trq(2)*s(2) - trq(3)*s(3)
      end if
c
c     force distribution for the Z-Only local coordinate method
c
      if (axetyp .eq. 'Z-Only') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            de(j,ialoc) = de(j,ialoc) + du
            de(j,ibloc) = de(j,ibloc) - du
            frcz(j) = frcz(j) + du
         end do
c
c     force distribution for the Z-then-X local coordinate method
c
      else if (axetyp .eq. 'Z-then-X') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin)
            de(j,ialoc) = de(j,ialoc) + du
            de(j,icloc) = de(j,icloc) + dv
            de(j,ibloc) = de(j,ibloc) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Bisector local coordinate method
c
      else if (axetyp .eq. 'Bisector') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + 0.5_ti_p*uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5_ti_p*vw(j)*dphidw/vsiz
            de(j,ialoc) = de(j,ialoc) + du
            de(j,icloc) = de(j,icloc) + dv
            de(j,ibloc) = de(j,ibloc) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Z-Bisect local coordinate method
c
      else if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
            du = ur(j)*dphidr/(usiz*ursin) + us(j)*dphids/usiz
            dv = (vssin*s(j)-vscos*t1(j))*dphidu
     &              / (vsiz*(ut1sin+ut2sin))
            dw = (wssin*s(j)-wscos*t2(j))*dphidu
     &              / (wsiz*(ut1sin+ut2sin))
            de(j,ialoc) = de(j,ialoc) + du
            de(j,icloc) = de(j,icloc) + dv
            de(j,idloc) = de(j,idloc) + dw
            de(j,ibloc) = de(j,ibloc) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
c
c     force distribution for the 3-Fold local coordinate method
c
      else if (polaxe(i) .eq. '3-Fold') then
         p(1) = u(1) + v(1) + w(1)
         p(2) = u(2) + v(2) + w(2)
         p(3) = u(3) + v(3) + w(3)
         psiz = sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3)) 
         do j = 1, 3 
            p(j) = p(j) / psiz
         end do
         wpcos = w(1)*p(1) + w(2)*p(2) + w(3)*p(3)
         upcos = u(1)*p(1) + u(2)*p(2) + u(3)*p(3)
         vpcos = v(1)*p(1) + v(2)*p(2) + v(3)*p(3)
         wpsin = sqrt(1.0_ti_p - wpcos*wpcos)
         upsin = sqrt(1.0_ti_p - upcos*upcos)
         vpsin = sqrt(1.0_ti_p - vpcos*vpcos)
         r(1) = u(1) + v(1)
         r(2) = u(2) + v(2)
         r(3) = u(3) + v(3)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
         end do
         rwcos = r(1)*w(1) + r(2)*w(2) + r(3)*w(3)
         rwsin = sqrt(1.0_ti_p - rwcos*rwcos)
         dphidr = -trq(1)*r(1) - trq(2)*r(2) - trq(3)*r(3)
         del(1) = r(2)*w(3) - r(3)*w(2)
         del(2) = r(3)*w(1) - r(1)*w(3) 
         del(3) = r(1)*w(2) - r(2)*w(1)
         delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))
         do j = 1, 3
            del(j) = del(j) / delsiz
         end do
         dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
         eps(1) = del(2)*w(3) - del(3)*w(2)
         eps(2) = del(3)*w(1) - del(1)*w(3)
         eps(3) = del(1)*w(2) - del(2)*w(1)
         do j = 1, 3
            dw = del(j)*dphidr/(wsiz*rwsin)
     &              + eps(j)*dphiddel*wpcos/(wsiz*psiz) 
            de(j,idloc) = de(j,idloc) + dw
            de(j,ibloc) = de(j,ibloc) - dw
            frcy(j) = frcy(j) + dw
         end do
         r(1) = v(1) + w(1)
         r(2) = v(2) + w(2)
         r(3) = v(3) + w(3)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
         end do
         rucos = r(1)*u(1) + r(2)*u(2) + r(3)*u(3)
         rusin = sqrt(1.0_ti_p - rucos*rucos) 
         dphidr = -trq(1)*r(1) - trq(2)*r(2) - trq(3)*r(3)
         del(1) = r(2)*u(3) - r(3)*u(2)
         del(2) = r(3)*u(1) - r(1)*u(3)
         del(3) = r(1)*u(2) - r(2)*u(1)
         delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))
         do j = 1, 3
            del(j) = del(j) / delsiz
         end do
         dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
         eps(1) = del(2)*u(3) - del(3)*u(2)
         eps(2) = del(3)*u(1) - del(1)*u(3)
         eps(3) = del(1)*u(2) - del(2)*u(1)
         do j = 1, 3
            du = del(j)*dphidr/(usiz*rusin)
     &              + eps(j)*dphiddel*upcos/(usiz*psiz)
            de(j,ialoc) = de(j,ialoc) + du
            de(j,ibloc) = de(j,ibloc) - du
            frcz(j) = frcz(j) + du
         end do
         r(1) = u(1) + w(1)
         r(2) = u(2) + w(2)
         r(3) = u(3) + w(3)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
         end do
         rvcos = r(1)*v(1) + r(2)*v(2) + r(3)*v(3) 
         rvsin = sqrt(1.0_ti_p - rvcos*rvcos)
         dphidr = -trq(1)*r(1) - trq(2)*r(2) - trq(3)*r(3)
         del(1) = r(2)*v(3) - r(3)*v(2)
         del(2) = r(3)*v(1) - r(1)*v(3)
         del(3) = r(1)*v(2) - r(2)*v(1)
         delsiz = sqrt(del(1)*del(1) + del(2)*del(2) + del(3)*del(3))
         do j = 1, 3
            del(j) = del(j) / delsiz
         end do
         dphiddel = -trq(1)*del(1) - trq(2)*del(2) - trq(3)*del(3)
         eps(1) = del(2)*v(3) - del(3)*v(2)
         eps(2) = del(3)*v(1) - del(1)*v(3)
         eps(3) = del(1)*v(2) - del(2)*v(1)
         do j = 1, 3
            dv = del(j)*dphidr/(vsiz*rvsin)
     &              + eps(j)*dphiddel*vpcos/(vsiz*psiz)
            de(j,icloc) = de(j,icloc) + dv
            de(j,ibloc) = de(j,ibloc) - dv
            frcx(j) = frcx(j) + dv
         end do
      end if
      return
      end
