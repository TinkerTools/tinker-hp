!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine ebond1  --  bond stretch energy & derivatives  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "ebond1" calculates the bond stretching energy and
!     first derivatives with respect to Cartesian coordinates
!
!
#include "tinker_macro.h"
subroutine ebond1
   use atmlst
   use atoms
   use bndpot
   use bond
   use bound
   use deriv
   use domdec
   use energi
   use group
   use mamd
   use potent ,only: use_amd_wat1
   use tinheader
   use usage
   use virial
   implicit none
   integer i,ia,ib,ialoc,ibloc
   integer ibond
   real(t_p) e,ideal,force
   real(t_p) expterm,bde
   real(t_p) dt,dt2,deddt
   real(t_p) de,dedx,dedy,dedz
   real(t_p) xab,yab,zab,rab
   real(t_p) vxx,vyy,vzz
   real(t_p) vyx,vzx,vzy
   real(t_p) fgrp
   logical proceed
!
!
!     zero out the bond energy and first derivatives
!
   eb = 0.0_ti_p
!
!     calculate the bond stretch energy and first derivatives
!
   do ibond = 1, nbondloc
      i  = bndglob(ibond)
      ia = ibnd(1,i)
      ib = ibnd(2,i)
      ialoc = loc(ia)
      ibloc = loc(ib)
      ideal = bl(i)
      force = bk(i)
!
!     decide whether to compute the current interaction
!
      proceed = (use(ia) .or. use(ib))
      if (use_group)  call groups (fgrp,ia,ib,0,0,0,0)
!
!     compute the value of the bond length deviation
!
      if (proceed) then
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         if (use_polymer)  call image (xab,yab,zab)
         rab = sqrt(xab*xab + yab*yab + zab*zab)
         dt = rab - ideal
!
!     harmonic potential uses Taylor expansion of Morse potential
!     through the fourth power of the bond length deviation
!
         if (bndtyp .eq. 'HARMONIC') then
            dt2 = dt * dt
            e = bndunit * force * dt2 * (1.0_ti_p+cbnd*dt+qbnd*dt2)
            deddt = 2.0_ti_p * bndunit * force * dt&
               &* (1.0_ti_p+1.5_ti_p*cbnd*dt+2.0_ti_p*qbnd*dt2)
!
!     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
!     with the approximations alpha = sqrt(ForceConst/BDE) = -2
!     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
!
         else if (bndtyp .eq. 'MORSE') then
            expterm = exp(-2.0_ti_p*dt)
            bde = 0.25_ti_p * bndunit * force
            e = bde * (1.0_ti_p-expterm)**2
            deddt = 4.0_ti_p * bde * (1.0_ti_p-expterm) * expterm
         end if
!
!     scale the interaction based on its group membership
!
         if (use_group) then
            e = e * fgrp
            deddt = deddt * fgrp
         end if
!
!     compute chain rule terms needed for derivatives
!
         if (rab .eq. 0.0_ti_p) then
            de = 0.0_ti_p
         else
            de = deddt / rab
         end if
         dedx = de * xab
         dedy = de * yab
         dedz = de * zab
!
!     increment the total bond energy and first derivatives
!
         eb = eb + e
         deb(1,ialoc) = deb(1,ialoc) + dedx
         deb(2,ialoc) = deb(2,ialoc) + dedy
         deb(3,ialoc) = deb(3,ialoc) + dedz
!
         deb(1,ibloc) = deb(1,ibloc) - dedx
         deb(2,ibloc) = deb(2,ibloc) - dedy
         deb(3,ibloc) = deb(3,ibloc) - dedz
!
!     aMD storage if waters are considered
!
         if (use_amd_wat1) then
            if (type(ia) == aMDwattype(1) .or. type(ib)&
               &== aMDwattype(1)) then
               eW1aMD = eW1aMD + e
               deW1aMD(1,ialoc) = deW1aMD(1,ialoc) + dedx
               deW1aMD(2,ialoc) = deW1aMD(2,ialoc) + dedy
               deW1aMD(3,ialoc) = deW1aMD(3,ialoc) + dedz
               deW1aMD(1,ibloc) = deW1aMD(1,ibloc) - dedx
               deW1aMD(2,ibloc) = deW1aMD(2,ibloc) - dedy
               deW1aMD(3,ibloc) = deW1aMD(3,ibloc) - dedz
            end if
         end if
!
!     increment the internal virial tensor components
!
         vxx = xab * dedx
         vyx = yab * dedx
         vzx = zab * dedx
         vyy = yab * dedy
         vzy = zab * dedy
         vzz = zab * dedz
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
   end do
end
