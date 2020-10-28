c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrbnd1   --  stretch-bend energy and derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrbnd1" calculates the stretch-bend potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_precision.h"
      subroutine estrbnd1
      use angle
      use angpot
      use atmlst
      use atoms
      use bond
      use bound
      use deriv
      use domdec
      use energi
      use group
      use math
      use mamd
      use potent,only:use_amd_wat1
      use strbnd
      use tinheader ,only:ti_p,re_p
      use usage
      use virial
      implicit none
      integer i,j,k,istrbnd,istrbndloc
      integer ia,ib,ic
      integer ialoc,ibloc,icloc
      real(t_p) e,dr1,dr2,dt
      real(t_p) angle1
      real(t_p) force1,force2
      real(t_p) dot,cosine
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xab,yab,zab
      real(t_p) xcb,ycb,zcb
      real(t_p) rab,rab2
      real(t_p) rcb,rcb2
      real(t_p) xp,yp,zp,rp
      real(t_p) term1,term2
      real(t_p) termr,term1t,term2t
      real(t_p) ddtdxia,ddtdyia,ddtdzia
      real(t_p) ddtdxic,ddtdyic,ddtdzic
      real(t_p) ddrdxia,ddrdyia,ddrdzia
      real(t_p) ddrdxic,ddrdyic,ddrdzic
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      logical proceed
c
!$acc update host(deba,vir)
c
c     zero out the energy and first derivative components
c
      eba = 0.0_ti_p
c
c     calculate the stretch-bend energy and first derivatives
c
      do istrbndloc = 1, nstrbndloc
         istrbnd = strbndglob(istrbndloc)
         i = isb(1,istrbnd)
         ia = iang(1,i)
         ialoc = loc(ia)
         ib = iang(2,i)
         ibloc = loc(ib)
         ic = iang(3,i)
         icloc = loc(ic)
         force1 = sbk(1,istrbnd)
         force2 = sbk(2,istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
c
c     get the coordinates of the atoms in the angle
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
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
            end if
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            if (rab2.ne.0.0_ti_p .and. rcb2.ne.0.0_ti_p) then
               rab = sqrt(rab2)
               rcb = sqrt(rcb2)
               xp = ycb*zab - zcb*yab
               yp = zcb*xab - xcb*zab
               zp = xcb*yab - ycb*xab
               rp = sqrt(xp*xp + yp*yp + zp*zp)
               rp = max(rp,0.001_ti_p)
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle1 = radian * acos(cosine)
c
c     find chain rule terms for the bond angle deviation
c
               dt = angle1 - anat(i)
               term1 = -radian / (rab2*rp)
               term2 = radian / (rcb2*rp)
               ddtdxia = term1 * (yab*zp-zab*yp)
               ddtdyia = term1 * (zab*xp-xab*zp)
               ddtdzia = term1 * (xab*yp-yab*xp)
               ddtdxic = term2 * (ycb*zp-zcb*yp)
               ddtdyic = term2 * (zcb*xp-xcb*zp)
               ddtdzic = term2 * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the bond length deviations
c
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               dr1 = rab - bl(j)
               term1 = 1.0_ti_p / rab
               dr2 = rcb - bl(k)
               term2 = 1.0_ti_p / rcb
               ddrdxia = term1 * xab
               ddrdyia = term1 * yab
               ddrdzia = term1 * zab
               ddrdxic = term2 * xcb
               ddrdyic = term2 * ycb
               ddrdzic = term2 * zcb
c
c     abbreviations used in defining chain rule terms
c
               term1 = stbnunit * force1
               term2 = stbnunit * force2
               termr = term1*dr1 + term2*dr2
               term1t = term1 * dt
               term2t = term2 * dt
c
c     get the energy and master chain rule terms for derivatives
c
               e = termr * dt
               dedxia = term1t*ddrdxia + termr*ddtdxia
               dedyia = term1t*ddrdyia + termr*ddtdyia
               dedzia = term1t*ddrdzia + termr*ddtdzia
               dedxic = term2t*ddrdxic + termr*ddtdxic
               dedyic = term2t*ddrdyic + termr*ddtdyic
               dedzic = term2t*ddrdzic + termr*ddtdzic
               dedxib = -dedxia - dedxic
               dedyib = -dedyia - dedyic
               dedzib = -dedzia - dedzic
c
c     increment the total stretch-bend energy and derivatives
c
               eba = eba + e
               deba(1,ibloc) = deba(1,ibloc) + dedxib
               deba(2,ibloc) = deba(2,ibloc) + dedyib
               deba(3,ibloc) = deba(3,ibloc) + dedzib
c
               deba(1,ialoc) = deba(1,ialoc) + dedxia
               deba(2,ialoc) = deba(2,ialoc) + dedyia
               deba(3,ialoc) = deba(3,ialoc) + dedzia
c
               deba(1,icloc) = deba(1,icloc) + dedxic
               deba(2,icloc) = deba(2,icloc) + dedyic
               deba(3,icloc) = deba(3,icloc) + dedzic
c
c     aMD storage if waters are considered
c
               if (use_amd_wat1) then
               if (type(ia) == aMDwattype(1) .or. type(ib)
     $         == aMDwattype(1) .or. type(ic) == aMDwattype(1)) then
                  eW1aMD = eW1aMD + e
                  deW1aMD(1,ialoc) = deW1aMD(1,ialoc) + dedxia
                  deW1aMD(2,ialoc) = deW1aMD(2,ialoc) + dedyia
                  deW1aMD(3,ialoc) = deW1aMD(3,ialoc) + dedzia
                  deW1aMD(1,ibloc) = deW1aMD(1,ibloc) + dedxib
                  deW1aMD(2,ibloc) = deW1aMD(2,ibloc) + dedyib
                  deW1aMD(3,ibloc) = deW1aMD(3,ibloc) + dedzib
                  deW1aMD(1,icloc) = deW1aMD(1,icloc) + dedxic
                  deW1aMD(2,icloc) = deW1aMD(2,icloc) + dedyic
                  deW1aMD(3,icloc) = deW1aMD(3,icloc) + dedzic
               end if
               end if
c
c     increment the internal virial tensor components
c
               vxx = xab*dedxia + xcb*dedxic
               vyx = yab*dedxia + ycb*dedxic
               vzx = zab*dedxia + zcb*dedxic
               vyy = yab*dedyia + ycb*dedyic
               vzy = zab*dedyia + zcb*dedyic
               vzz = zab*dedzia + zcb*dedzic
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
            end if
         end if
      end do
!$acc update device(deba,vir)
      return
      end
