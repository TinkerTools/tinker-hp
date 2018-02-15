c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torque  --  convert single site torque to force  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torque" takes the torque values on a single site defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame
c
c     literature reference:
c
c     P. L. Popelier and A. J. Stone, "Formulae for the First and
c     Second Derivatives of Anisotropic Potentials with Respect to
c     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)
c
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torque2  --  convert all site torques to forces  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torque2" takes the torque values on all sites defined by
c     local coordinate frames and finds the total Cartesian force
c     components on original sites and sites specifying local frames
c
c     literature reference:
c
c     P. L. Popelier and A. J. Stone, "Formulae for the First and
c     Second Derivatives of Anisotropic Potentials with Respect to
c     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)
c
c
      subroutine torque2 (trq,torquerec)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'mpole.i'
      include 'openmp.i'
      integer i,j,iipole,iglob
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      real*8 du,dv,dw,random
      real*8 usiz,vsiz,wsiz
      real*8 rsiz,ssiz
      real*8 t1siz,t2siz
      real*8 uvsiz,uwsiz,vwsiz
      real*8 ursiz,ussiz
      real*8 vssiz,wssiz
      real*8 uvcos,uwcos,vwcos
      real*8 urcos,uscos
      real*8 vscos,wscos
      real*8 ut1cos,ut2cos
      real*8 uvsin,uwsin,vwsin
      real*8 ursin,ussin
      real*8 vssin,wssin
      real*8 ut1sin,ut2sin
      real*8 dphidu,dphidv,dphidw
      real*8 dphidr,dphids
      real*8 u(3),v(3),w(3)
      real*8 r(3),s(3)
      real*8 t1(3),t2(3)
      real*8 uv(3),uw(3),vw(3)
      real*8 ur(3),us(3)
      real*8 vs(3),ws(3)
      real*8 trq(3,*),torquerec(3,*)
      character*8 axetyp
c
c
c     get the local frame type and the frame-defining atoms
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         axetyp = polaxe(iipole)
         if (axetyp .eq. 'None')  goto 10
         ia = zaxis(iipole)
c         if (ia.eq.0) write(*,*) 'beware zaxis of atom ',iglob,
c     $    ' not defined'
         if (ia.ne.0) ialoc = locrec1(ia)
         ib = iglob
         ibloc = locrec1(ib)
         ic = xaxis(iipole)
         if (ic.ne.0) icloc = locrec1(ic)
         id = yaxis(iipole)
         if (id.ne.0) idloc = locrec1(id)
c
c     construct the three rotation axes for the local frame
c
         u(1) = x(ia) - x(ib)
         u(2) = y(ia) - y(ib)
         u(3) = z(ia) - z(ib)
         if (axetyp .ne. 'Z-Only') then
            v(1) = x(ic) - x(ib)
            v(2) = y(ic) - y(ib)
            v(3) = z(ic) - z(ib)
         else
            v(1) = random ()
            v(2) = random ()
            v(3) = random ()
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
         usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
         vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
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
            s(1) = u(2)*r(3) - u(3)*r(2)
            s(2) = u(3)*r(1) - u(1)*r(3)
            s(3) = u(1)*r(2) - u(2)*r(1)
            rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
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
         uw(1) = w(2)*u(3) - w(3)*u(2)
         uw(2) = w(3)*u(1) - w(1)*u(3)
         uw(3) = w(1)*u(2) - w(2)*u(1)
         vw(1) = w(2)*v(3) - w(3)*v(2)
         vw(2) = w(3)*v(1) - w(1)*v(3)
         vw(3) = w(1)*v(2) - w(2)*v(1)
         uvsiz = sqrt(uv(1)*uv(1) + uv(2)*uv(2) + uv(3)*uv(3))
         uwsiz = sqrt(uw(1)*uw(1) + uw(2)*uw(2) + uw(3)*uw(3))
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
            us(1) = s(2)*u(3) - s(3)*u(2)
            us(2) = s(3)*u(1) - s(1)*u(3)
            us(3) = s(1)*u(2) - s(2)*u(1)
            vs(1) = s(2)*v(3) - s(3)*v(2)
            vs(2) = s(3)*v(1) - s(1)*v(3)
            vs(3) = s(1)*v(2) - s(2)*v(1)
            ws(1) = s(2)*w(3) - s(3)*w(2)
            ws(2) = s(3)*w(1) - s(1)*w(3)
            ws(3) = s(1)*w(2) - s(2)*w(1)
            ursiz = sqrt(ur(1)*ur(1) + ur(2)*ur(2) + ur(3)*ur(3))
            ussiz = sqrt(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
            vssiz = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
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
         uvsin = sqrt(1.0d0 - uvcos*uvcos)
         uwcos = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
         uwsin = sqrt(1.0d0 - uwcos*uwcos)
         vwcos = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
         vwsin = sqrt(1.0d0 - vwcos*vwcos)
         if (axetyp .eq. 'Z-Bisect') then
            urcos = u(1)*r(1) + u(2)*r(2) + u(3)*r(3)
            ursin = sqrt(1.0d0 - urcos*urcos)
            uscos = u(1)*s(1) + u(2)*s(2) + u(3)*s(3)
            ussin = sqrt(1.0d0 - uscos*uscos)
            vscos = v(1)*s(1) + v(2)*s(2) + v(3)*s(3)
            vssin = sqrt(1.0d0 - vscos*vscos)
            wscos = w(1)*s(1) + w(2)*s(2) + w(3)*s(3)
            wssin = sqrt(1.0d0 - wscos*wscos)
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
            ut1sin = sqrt(1.0d0 - ut1cos*ut1cos)
            ut2cos = u(1)*t2(1) + u(2)*t2(2) + u(3)*t2(3)
            ut2sin = sqrt(1.0d0 - ut2cos*ut2cos)
         end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
         dphidu = -trq(1,i)*u(1) - trq(2,i)*u(2) - trq(3,i)*u(3)
         dphidv = -trq(1,i)*v(1) - trq(2,i)*v(2) - trq(3,i)*v(3)
         dphidw = -trq(1,i)*w(1) - trq(2,i)*w(2) - trq(3,i)*w(3)
         if (axetyp .eq. 'Z-Bisect') then
            dphidr = -trq(1,i)*r(1) - trq(2,i)*r(2) - trq(3,i)*r(3)
            dphids = -trq(1,i)*s(1) - trq(2,i)*s(2) - trq(3,i)*s(3)
         end if
c
c     force distribution for the Z-Only local coordinate method
c
         if (axetyp .eq. 'Z-Only') then
            do j = 1, 3
               du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
               torquerec(j,ialoc) = torquerec(j,ialoc) + du
               torquerec(j,ibloc) = torquerec(j,ibloc) - du
            end do
c
c     force distribution for the Z-then-X local coordinate method
c
         else if (axetyp .eq. 'Z-then-X') then
            do j = 1, 3
               du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
               dv = -uv(j)*dphidu/(vsiz*uvsin)
               torquerec(j,ialoc) = torquerec(j,ialoc) + du
               torquerec(j,icloc) = torquerec(j,icloc) + dv
               torquerec(j,ibloc) = torquerec(j,ibloc) - du - dv
            end do
c
c     force distribution for the Bisector local coordinate method
c
         else if (axetyp .eq. 'Bisector') then
            do j = 1, 3
               du = uv(j)*dphidv/(usiz*uvsin) + 0.5d0*uw(j)*dphidw/usiz
               dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5d0*vw(j)*dphidw/vsiz
               torquerec(j,ialoc) = torquerec(j,ialoc) + du
               torquerec(j,icloc) = torquerec(j,icloc) + dv
               torquerec(j,ibloc) = torquerec(j,ibloc) - du - dv
            end do
c
c     force distribution for the Z-Bisect local coordinate method
c
         else if (axetyp .eq. 'Z-Bisect') then
            do j = 1, 3
               du = ur(j)*dphidr/(usiz*ursin) + us(j)*dphids/usiz
               dv = (vssin*s(j)-vscos*t1(j))*dphidu
     &                 / (vsiz*(ut1sin+ut2sin))
               dw = (wssin*s(j)-wscos*t2(j))*dphidu
     &                 / (wsiz*(ut1sin+ut2sin))
               torquerec(j,ialoc) = torquerec(j,ialoc) + du
               torquerec(j,icloc) = torquerec(j,icloc) + dv
               torquerec(j,idloc) = torquerec(j,idloc) + dw
               torquerec(j,ibloc) = torquerec(j,ibloc) - du - dv - dw
            end do
c
c     force distribution for the 3-Fold local coordinate method
c        (correct for uv, uw and vw angles all equal to 90)
c
         else if (axetyp .eq. '3-Fold') then
            do j = 1, 3
               du = uw(j)*dphidw/(usiz*uwsin)
     &                 + uv(j)*dphidv/(usiz*uvsin)
     &                 - uw(j)*dphidu/(usiz*uwsin)
     &                 - uv(j)*dphidu/(usiz*uvsin)
               dv = vw(j)*dphidw/(vsiz*vwsin)
     &                 - uv(j)*dphidu/(vsiz*uvsin)
     &                 - vw(j)*dphidv/(vsiz*vwsin)
     &                 + uv(j)*dphidv/(vsiz*uvsin)
               dw = -uw(j)*dphidu/(wsiz*uwsin)
     &                 - vw(j)*dphidv/(wsiz*vwsin)
     &                 + uw(j)*dphidw/(wsiz*uwsin)
     &                 + vw(j)*dphidw/(wsiz*vwsin)
               du = du / 3.0d0
               dv = dv / 3.0d0
               dw = dw / 3.0d0
               torquerec(j,ialoc) = torquerec(j,ialoc) + du
               torquerec(j,icloc) = torquerec(j,icloc) + dv
               torquerec(j,idloc) = torquerec(j,idloc) + dw
               torquerec(j,ibloc) = torquerec(j,ibloc) - du - dv - dw
            end do
         end if
   10    continue
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torque3  --  convert torque to force components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torque3" takes the torque values on a single site defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame; this
c     version also returns the individual atomic components
c
c     literature reference:
c
c     P. L. Popelier and A. J. Stone, "Formulae for the First and
c     Second Derivatives of Anisotropic Potentials with Respect to
c     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)
c
c
      subroutine torque3 (i,trq1,trq2,frcx,frcy,frcz,torquedirt,
     $   torquedirpt)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'openmp.i'
      integer i,j
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      real*8 du,dv,dw,random
      real*8 usiz,vsiz,wsiz
      real*8 rsiz,ssiz
      real*8 t1siz,t2siz
      real*8 uvsiz,uwsiz,vwsiz
      real*8 ursiz,ussiz
      real*8 vssiz,wssiz
      real*8 uvcos,uwcos,vwcos
      real*8 urcos,uscos
      real*8 vscos,wscos
      real*8 ut1cos,ut2cos
      real*8 uvsin,uwsin,vwsin
      real*8 ursin,ussin
      real*8 vssin,wssin
      real*8 ut1sin,ut2sin
      real*8 dphidu,dphidv,dphidw
      real*8 dphidr,dphids
      real*8 trq1(3),trq2(3)
      real*8 frcx(3),frcy(3),frcz(3)
      real*8 u(3),v(3),w(3)
      real*8 r(3),s(3)
      real*8 t1(3),t2(3)
      real*8 uv(3),uw(3),vw(3)
      real*8 ur(3),us(3)
      real*8 vs(3),ws(3)
      real*8 torquedirt(3,*),torquedirpt(3,*)
      character*8 axetyp
c
c     zero out force components on local frame-defining atoms
c
      do j = 1, 3
         frcz(j) = 0.0d0
         frcx(j) = 0.0d0
         frcy(j) = 0.0d0
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
      if (axetyp .ne. 'Z-Only') then
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
      else
         v(1) = random ()
         v(2) = random ()
         v(3) = random ()
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
      usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
      vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
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
         s(1) = u(2)*r(3) - u(3)*r(2)
         s(2) = u(3)*r(1) - u(1)*r(3)
         s(3) = u(1)*r(2) - u(2)*r(1)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
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
      uw(1) = w(2)*u(3) - w(3)*u(2)
      uw(2) = w(3)*u(1) - w(1)*u(3)
      uw(3) = w(1)*u(2) - w(2)*u(1)
      vw(1) = w(2)*v(3) - w(3)*v(2)
      vw(2) = w(3)*v(1) - w(1)*v(3)
      vw(3) = w(1)*v(2) - w(2)*v(1)
      uvsiz = sqrt(uv(1)*uv(1) + uv(2)*uv(2) + uv(3)*uv(3))
      uwsiz = sqrt(uw(1)*uw(1) + uw(2)*uw(2) + uw(3)*uw(3))
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
         us(1) = s(2)*u(3) - s(3)*u(2)
         us(2) = s(3)*u(1) - s(1)*u(3)
         us(3) = s(1)*u(2) - s(2)*u(1)
         vs(1) = s(2)*v(3) - s(3)*v(2)
         vs(2) = s(3)*v(1) - s(1)*v(3)
         vs(3) = s(1)*v(2) - s(2)*v(1)
         ws(1) = s(2)*w(3) - s(3)*w(2)
         ws(2) = s(3)*w(1) - s(1)*w(3)
         ws(3) = s(1)*w(2) - s(2)*w(1)
         ursiz = sqrt(ur(1)*ur(1) + ur(2)*ur(2) + ur(3)*ur(3))
         ussiz = sqrt(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
         vssiz = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
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
      uvsin = sqrt(1.0d0 - uvcos*uvcos)
      uwcos = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
      uwsin = sqrt(1.0d0 - uwcos*uwcos)
      vwcos = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
      vwsin = sqrt(1.0d0 - vwcos*vwcos)
      if (axetyp .eq. 'Z-Bisect') then
         urcos = u(1)*r(1) + u(2)*r(2) + u(3)*r(3)
         ursin = sqrt(1.0d0 - urcos*urcos)
         uscos = u(1)*s(1) + u(2)*s(2) + u(3)*s(3)
         ussin = sqrt(1.0d0 - uscos*uscos)
         vscos = v(1)*s(1) + v(2)*s(2) + v(3)*s(3)
         vssin = sqrt(1.0d0 - vscos*vscos)
         wscos = w(1)*s(1) + w(2)*s(2) + w(3)*s(3)
         wssin = sqrt(1.0d0 - wscos*wscos)
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
         ut1sin = sqrt(1.0d0 - ut1cos*ut1cos)
         ut2cos = u(1)*t2(1) + u(2)*t2(2) + u(3)*t2(3)
         ut2sin = sqrt(1.0d0 - ut2cos*ut2cos)
      end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      dphidu = -trq1(1)*u(1) - trq1(2)*u(2) - trq1(3)*u(3)
      dphidv = -trq1(1)*v(1) - trq1(2)*v(2) - trq1(3)*v(3)
      dphidw = -trq1(1)*w(1) - trq1(2)*w(2) - trq1(3)*w(3)
      if (axetyp .eq. 'Z-Bisect') then
         dphidr = -trq1(1)*r(1) - trq1(2)*r(2) - trq1(3)*r(3)
         dphids = -trq1(1)*s(1) - trq1(2)*s(2) - trq1(3)*s(3)
      end if
c
c     force distribution for the Z-Only local coordinate method
c
      if (axetyp .eq. 'Z-Only') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            torquedirt(j,ialoc) = torquedirt(j,ialoc) + du
            torquedirt(j,ibloc) = torquedirt(j,ibloc) - du
            frcz(j) = frcz(j) + du
         end do
c
c     force distribution for the Z-then-X local coordinate method
c
      else if (axetyp .eq. 'Z-then-X') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin)
            torquedirt(j,ialoc) = torquedirt(j,ialoc) + du
            torquedirt(j,icloc) = torquedirt(j,icloc) + dv
            torquedirt(j,ibloc) = torquedirt(j,ibloc) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Bisector local coordinate method
c
      else if (axetyp .eq. 'Bisector') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + 0.5d0*uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5d0*vw(j)*dphidw/vsiz
            torquedirt(j,ialoc) = torquedirt(j,ialoc) + du
            torquedirt(j,icloc) = torquedirt(j,icloc) + dv
            torquedirt(j,ibloc) = torquedirt(j,ibloc) - du - dv
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
            torquedirt(j,ialoc) = torquedirt(j,ialoc) + du
            torquedirt(j,icloc) = torquedirt(j,icloc) + dv
            torquedirt(j,idloc) = torquedirt(j,idloc) + dw
            torquedirt(j,ibloc) = torquedirt(j,ibloc) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
c
c     force distribution for the 3-Fold local coordinate method
c        (correct for uv, uw and vw angles all equal to 90)
c
      else if (axetyp .eq. '3-Fold') then
         do j = 1, 3
            du = uw(j)*dphidw/(usiz*uwsin)
     &              + uv(j)*dphidv/(usiz*uvsin)
     &              - uw(j)*dphidu/(usiz*uwsin)
     &              - uv(j)*dphidu/(usiz*uvsin)
            dv = vw(j)*dphidw/(vsiz*vwsin)
     &              - uv(j)*dphidu/(vsiz*uvsin)
     &              - vw(j)*dphidv/(vsiz*vwsin)
     &              + uv(j)*dphidv/(vsiz*uvsin)
            dw = -uw(j)*dphidu/(wsiz*uwsin)
     &              - vw(j)*dphidv/(wsiz*vwsin)
     &              + uw(j)*dphidw/(wsiz*uwsin)
     &              + vw(j)*dphidw/(wsiz*vwsin)
            du = du / 3.0d0
            dv = dv / 3.0d0
            dw = dw / 3.0d0
            torquedirt(j,ialoc) = torquedirt(j,ialoc) + du
            torquedirt(j,icloc) = torquedirt(j,icloc) + dv
            torquedirt(j,idloc) = torquedirt(j,idloc) + dw
            torquedirt(j,ibloc) = torquedirt(j,ibloc) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
      end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      dphidu = -trq2(1)*u(1) - trq2(2)*u(2) - trq2(3)*u(3)
      dphidv = -trq2(1)*v(1) - trq2(2)*v(2) - trq2(3)*v(3)
      dphidw = -trq2(1)*w(1) - trq2(2)*w(2) - trq2(3)*w(3)
      if (axetyp .eq. 'Z-Bisect') then
         dphidr = -trq2(1)*r(1) - trq2(2)*r(2) - trq2(3)*r(3)
         dphids = -trq2(1)*s(1) - trq2(2)*s(2) - trq2(3)*s(3)
      end if
c
c     force distribution for the Z-Only local coordinate method
c
      if (axetyp .eq. 'Z-Only') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            torquedirpt(j,ialoc) = torquedirpt(j,ialoc) + du
            torquedirpt(j,ibloc) = torquedirpt(j,ibloc) - du
            frcz(j) = frcz(j) + du
         end do
c
c     force distribution for the Z-then-X local coordinate method
c
      else if (axetyp .eq. 'Z-then-X') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin)
            torquedirpt(j,ialoc) = torquedirpt(j,ialoc) + du
            torquedirpt(j,icloc) = torquedirpt(j,icloc) + dv
            torquedirpt(j,ibloc) = torquedirpt(j,ibloc) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Bisector local coordinate method
c
      else if (axetyp .eq. 'Bisector') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + 0.5d0*uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5d0*vw(j)*dphidw/vsiz
            torquedirpt(j,ialoc) = torquedirpt(j,ialoc) + du
            torquedirpt(j,icloc) = torquedirpt(j,icloc) + dv
            torquedirpt(j,ibloc) = torquedirpt(j,ibloc) - du - dv
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
            torquedirpt(j,ialoc) = torquedirpt(j,ialoc) + du
            torquedirpt(j,icloc) = torquedirpt(j,icloc) + dv
            torquedirpt(j,idloc) = torquedirpt(j,idloc) + dw
            torquedirpt(j,ibloc) = torquedirpt(j,ibloc) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
c
c     force distribution for the 3-Fold local coordinate method
c        (correct for uv, uw and vw angles all equal to 90)
c
      else if (axetyp .eq. '3-Fold') then
         do j = 1, 3
            du = uw(j)*dphidw/(usiz*uwsin)
     &              + uv(j)*dphidv/(usiz*uvsin)
     &              - uw(j)*dphidu/(usiz*uwsin)
     &              - uv(j)*dphidu/(usiz*uvsin)
            dv = vw(j)*dphidw/(vsiz*vwsin)
     &              - uv(j)*dphidu/(vsiz*uvsin)
     &              - vw(j)*dphidv/(vsiz*vwsin)
     &              + uv(j)*dphidv/(vsiz*uvsin)
            dw = -uw(j)*dphidu/(wsiz*uwsin)
     &              - vw(j)*dphidv/(wsiz*vwsin)
     &              + uw(j)*dphidw/(wsiz*uwsin)
     &              + vw(j)*dphidw/(wsiz*vwsin)
            du = du / 3.0d0
            dv = dv / 3.0d0
            dw = dw / 3.0d0
            torquedirpt(j,ialoc) = torquedirpt(j,ialoc) + du
            torquedirpt(j,icloc) = torquedirpt(j,icloc) + dv
            torquedirpt(j,idloc) = torquedirpt(j,idloc) + dw
            torquedirpt(j,ibloc) = torquedirpt(j,ibloc) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
      end if
      return
c
      end
