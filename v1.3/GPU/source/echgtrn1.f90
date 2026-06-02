!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine echgtrn1  --  charge transfer energy & derivs  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "echgtrn1" calculates the charge transfer energy and first
!     derivatives with respect to Cartesian coordinates
!
!
subroutine echgtrn1
   implicit none
!
!
!     choose method for summing over charge transfer interactions
!
   call echgtrn1c
   return
end
!
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine echgtrn1c  --  charge transfer derivs via list  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "echgtrn1c" calculates the charge transfer energy and first
!     derivatives using a pairwise neighbor list
!
!
subroutine echgtrn1c
   use atoms
   use atmlst
   use bound
   use chgpot
   use chgtrn
   use cell
   use couple
   use ctrpot
   use cutoff
   use deriv
   use domdec
   use energi
   use group
   use inform
   use iounit
   use mplpot
   use mpole
   use mutant
   use neigh
   use potent
   use shunt
   use usage
   use virial
   implicit none
   integer i,j,iipole,iglob,kglob,nnelst
   integer ii,kkk,kbis,kkpole
   real*8 e,de,f,fgrp
   real*8 rr1,r,r2
   real*8 r3,r4,r5
   real*8 xi,yi,zi
   real*8 xr,yr,zr
   real*8 chgi,chgk
   real*8 chgik
   real*8 alphai,alphak
   real*8 alphaik
   real*8 expi,expk
   real*8 expik
   real*8 frcx,frcy,frcz
   real*8 vxx,vyy,vzz
   real*8 vxy,vxz,vyz
   real*8 taper,dtaper
   real*8 s,ds,ctrnshortcut2,facts,factds
   real*8, allocatable :: mscale(:)
   real*8 delambdaetemp
   logical proceed,usei
   logical muti,mutk
   logical testcut,shortrange,longrange,fullrange
   character*11 mode
   character*80 :: RoutineName
!
   if (deb_Path) write(iout,*), 'echgtrn1c '
!
   shortrange = use_chgtrnshort
   longrange  = use_chgtrnlong
   fullrange  = .not.(shortrange.or.longrange)

   if (shortrange) then
      RoutineName = 'echgtrnshort1c'
      mode        = 'SHORTCHGTRN'
   else if (longrange) then
      RoutineName = 'echgtrnlong1c'
      mode        = 'CHGTRN'
   else
      RoutineName = 'echgtrnl1c'
      mode        = 'CHGTRN'
   endif
!
!
!     zero out the charge transfer energy and first derivatives
!
   ect = 0.0d0
   dect = 0d0
   if (nct .eq. 0)  return
!
!     perform dynamic allocation of some local arrays
!
   allocate (mscale(n))
!
!     initialize connected atom exclusion coefficients
!
   mscale = 1.0d0
!
!     set the coefficients for the switching function
!
   call switch (mode)
   ctrnshortcut2 = (ctrnshortcut-shortheal)**2
!
!     set conversion factor
!
   f = electric / dielec
!
!     compute the charge transfer energy and derivatives
!
   do ii = 1, npolelocnl
      iipole = poleglobnl(ii)
      iglob = ipole(iipole)
      i = loc(iglob)
      xi = x(iglob)
      yi = y(iglob)
      zi = z(iglob)
      chgi = chgct(iipole)
      alphai = dmpct(iipole)
      if (alphai .eq. 0.0d0)  alphai = 100.0d0
      usei = use(iglob)
      muti = mut(iglob)
!
!     set exclusion coefficients for connected atoms
!
      do j = 1, n12(iglob)
         mscale(i12(j,iglob)) = m2scale
      end do
      do j = 1, n13(iglob)
         mscale(i13(j,iglob)) = m3scale
      end do
      do j = 1, n14(iglob)
         mscale(i14(j,iglob)) = m4scale
      end do
      do j = 1, n15(iglob)
         mscale(i15(j,iglob)) = m5scale
      end do
!
!     evaluate all sites within the cutoff distance
!
      nnelst = merge(nshortelst(ii),&
         &nelst     (ii),&
         &shortrange&
         &)
      do kkk = 1, nnelst
         kkpole = merge(shortelst(kkk,ii),&
            &elst     (kkk,ii),&
            &shortrange&
            &)
         kglob = ipole(kkpole)
         kbis = loc(kglob)
         mutk = mut(kglob)
         proceed = (usei .or. use(kglob))
         if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
         if (proceed) then
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            testcut = merge(r2 .le. off2.and.r2.ge.ctrnshortcut2,&
            &r2 .le. off2,&
            &longrange&
            &)
            if (testcut) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               chgk = chgct(kkpole)
               alphak = dmpct(kkpole)
               if (alphak .eq. 0.0d0)  alphak = 100.0d0
               if (ctrntyp .eq. 'SEPARATE') then
                  expi = exp(-alphai*r)
                  expk = exp(-alphak*r)
                  e = -chgi*expk - chgk*expi
                  de = chgi*expk*alphak + chgk*expi*alphai
               else
                  chgik = sqrt(abs(chgi*chgk))
                  alphaik = 0.5d0 * (alphai+alphak)
                  expik = exp(-alphaik*r)
                  e = -chgik * expik
                  de = -e * alphaik
               end if
               e = f * e * mscale(kglob)
               de = f * de * mscale(kglob)
!
!     use energy switching if near the cutoff distance
!
               if(longrange.or.fullrange) then
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3&
                     &+ c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3&
                     &+ 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
               end if
!
!     scale the interaction based on its group membership
!
               if (use_group) then
                  e = e * fgrp
                  de = de * fgrp
               end if
!
!     apply lambda scaling for interaction annihilation 
!
               if (muti .or. mutk) then
                  if (use_lambdadyn) then
                    delambdaetemp = e
                  end if
                  e = e * elambda
                  de = de * elambda
               end if
!
!     use energy switching if close the cutoff distance (at short range)
!
               if(shortrange .or. longrange)&
               &call switch_respa(r,ctrnshortcut,shortheal,s,ds)

               if(shortrange) then
                  facts =          s
                  factds =        ds
               else if(longrange) then
                  facts  = 1.0d0 - s
                  factds =       -ds
               else
                  facts  = 1.0d0
                  factds = 0.0d0
               endif
               de = de * facts + e * factds
               e  = e  * facts

               if (use_lambdadyn) then
                  delambdaetemp = facts*delambdaetemp
                  if (shortrange) then
                     delambdaesave = delambdaesave + delambdaetemp
                  else
                     delambdae = delambdae + delambdaetemp
                  end if
               end if
!
!     compute the force components for this interaction
!
               frcx = de * xr * rr1
               frcy = de * yr * rr1
               frcz = de * zr * rr1
!
!     increment the total charge transfer energy and derivatives
!
               ect = ect + e
               dect(1,i) = dect(1,i) - frcx
               dect(2,i) = dect(2,i) - frcy
               dect(3,i) = dect(3,i) - frcz
               dect(1,kbis) = dect(1,kbis) + frcx
               dect(2,kbis) = dect(2,kbis) + frcy
               dect(3,kbis) = dect(3,kbis) + frcz
!
!     increment the internal virial tensor components
!
               vxx = xr * frcx
               vxy = yr * frcx
               vxz = zr * frcx
               vyy = yr * frcy
               vyz = zr * frcy
               vzz = zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
         end if
      end do
!
!     reset exclusion coefficients for connected atoms
!
      do j = 1, n12(iglob)
         mscale(i12(j,iglob)) = 1.0d0
      end do
      do j = 1, n13(iglob)
         mscale(i13(j,iglob)) = 1.0d0
      end do
      do j = 1, n14(iglob)
         mscale(i14(j,iglob)) = 1.0d0
      end do
      do j = 1, n15(iglob)
         mscale(i15(j,iglob)) = 1.0d0
      end do
   end do
!
!     perform deallocation of some local arrays
!
   deallocate (mscale)
   return
end
