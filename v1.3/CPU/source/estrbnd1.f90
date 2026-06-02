!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine estrbnd1   --  stretch-bend energy and derivs  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "estrbnd1" calculates the stretch-bend potential energy and
!     first derivatives with respect to Cartesian coordinates
!
!
!> @brief 
!> calculates the stretch-bend potential energy and
!> first derivatives with respect to Cartesian coordinates
!> @param no params
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
   use inform
   use iounit
   use math
   use strbnd
   use usage
   use virial
   implicit none
   integer i,j,k,istrbnd,istrbndloc
   integer ia,ib,ic
   integer ialoc,ibloc,icloc
   real*8 e,dr1,dr2,dt
   real*8 angle1
   real*8 force1,force2
   real*8 dot,cosine
   real*8 xia,yia,zia
   real*8 xib,yib,zib
   real*8 xic,yic,zic
   real*8 xab,yab,zab
   real*8 xcb,ycb,zcb
   real*8 rab,rab2
   real*8 rcb,rcb2
   real*8 xp,yp,zp,rp
   real*8 term1,term2
   real*8 termr,term1t,term2t
   real*8 ddtdxia,ddtdyia,ddtdzia
   real*8 ddtdxic,ddtdyic,ddtdzic
   real*8 ddrdxia,ddrdyia,ddrdzia
   real*8 ddrdxic,ddrdyic,ddrdzic
   real*8 dedxia,dedyia,dedzia
   real*8 dedxib,dedyib,dedzib
   real*8 dedxic,dedyic,dedzic
   real*8 vxx,vyy,vzz
   real*8 vyx,vzx,vzy
   real*8 fgrp
   logical proceed
!
   if (deb_Path) write(iout,*), 'estrbnd1 '
!
!
!     zero out the energy and first derivative components
!
   eba = 0.0d0
!
!     calculate the stretch-bend energy and first derivatives
!
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
!
!     decide whether to compute the current interaction
!
      if (use_group)  call groups (fgrp,ia,ib,ic,0,0,0)
      proceed = (use(ia) .or. use(ib) .or. use(ic))
!
!     get the coordinates of the atoms in the angle
!
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
!
!     compute the value of the bond angle
!
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
         if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
            rab = sqrt(rab2)
            rcb = sqrt(rcb2)
            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            rp = sqrt(xp*xp + yp*yp + zp*zp)
            rp = max(rp,0.001d0)
            dot = xab*xcb + yab*ycb + zab*zcb
            cosine = dot / (rab*rcb)
            cosine = min(1.0d0,max(-1.0d0,cosine))
            angle1 = radian * acos(cosine)
!
!     find chain rule terms for the bond angle deviation
!
            dt = angle1 - anat(i)
            term1 = -radian / (rab2*rp)
            term2 = radian / (rcb2*rp)
            ddtdxia = term1 * (yab*zp-zab*yp)
            ddtdyia = term1 * (zab*xp-xab*zp)
            ddtdzia = term1 * (xab*yp-yab*xp)
            ddtdxic = term2 * (ycb*zp-zcb*yp)
            ddtdyic = term2 * (zcb*xp-xcb*zp)
            ddtdzic = term2 * (xcb*yp-ycb*xp)
!
!     find chain rule terms for the bond length deviations
!
            j = isb(2,istrbnd)
            k = isb(3,istrbnd)
            dr1 = rab - bl(j)
            term1 = 1.0d0 / rab
            dr2 = rcb - bl(k)
            term2 = 1.0d0 / rcb
            ddrdxia = term1 * xab
            ddrdyia = term1 * yab
            ddrdzia = term1 * zab
            ddrdxic = term2 * xcb
            ddrdyic = term2 * ycb
            ddrdzic = term2 * zcb
!
!     abbreviations used in defining chain rule terms
!
            term1 = stbnunit * force1
            term2 = stbnunit * force2
            termr = term1*dr1 + term2*dr2
            term1t = term1 * dt
            term2t = term2 * dt
!
!     scale the interaction based on its group membership
!
            if (use_group) then
               termr = termr * fgrp
               term1t = term1t * fgrp
               term2t = term2t * fgrp
            end if
!
!     get the energy and master chain rule terms for derivatives
!
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
!
!     increment the total stretch-bend energy and derivatives
!
            eba = eba + e
            deba(1,ibloc) = deba(1,ibloc) + dedxib
            deba(2,ibloc) = deba(2,ibloc) + dedyib
            deba(3,ibloc) = deba(3,ibloc) + dedzib
!
            deba(1,ialoc) = deba(1,ialoc) + dedxia
            deba(2,ialoc) = deba(2,ialoc) + dedyia
            deba(3,ialoc) = deba(3,ialoc) + dedzia
!
            deba(1,icloc) = deba(1,icloc) + dedxic
            deba(2,icloc) = deba(2,icloc) + dedyic
            deba(3,icloc) = deba(3,icloc) + dedzic
!
!     increment the internal virial tensor components
!
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
   return
end
