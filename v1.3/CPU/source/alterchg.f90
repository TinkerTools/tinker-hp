!
!
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine alterchg  --  modification of partial charges  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "alterchg" calculates the change in atomic partial charge or
!     monopole values due to bond and angle charge flux coupling
!
!     literature reference:
!
!     C. Liu, J.-P. Piquemal and P. Ren, "Implementation of Geometry-
!     Dependent Charge Flux into the Polarizable AMOEBA+ Potential",
!     Journal of Physical Chemistry Letters, 11, 419-426 (2020)
!
!
!> @brief 
!> calculates the change in atomic partial charge or monopole values due to bond and angle charge flux coupling 
!> @param no params
subroutine alterchg
   use atoms
   use atmlst
   use charge
   use chgpen
   use domdec
   use inform
   use iounit
   use mplpot
   use mpole
   use potent
   implicit none
   integer i,k,ii
   integer iloc,ierr
   real*8, allocatable :: pdelta(:)
   logical header
!
   if (deb_Path) write(iout,*), 'alterchg '
!
!     perform dynamic allocation of some local arrays
!
   allocate (pdelta(n))
!
!     zero out the change in charge value at each site
!
   do i = 1, n
      pdelta(i) = 0.0d0
   end do
!
!     find charge modifications due to charge flux
!
   call bndchg (pdelta)
   call angchg (pdelta)
!
!     communicate neighboring values of delta to modify charges
!
   call commchgflx(pdelta)
   call MPI_BARRIER(hostcomm,ierr)
!
!     alter atomic partial charge values for charge flux
!
   header = .true.
   do iloc = 1, nionbloc
      i = chgglob(iloc)
      k = iion(i)
      pchg(i) = pchg0(i) + pdelta(k)
      if (debug .and. pdelta(k).ne.0.0d0) then
         if (header) then
            header = .false.
            write (iout,10)
10          format (/,' Charge Flux Modification of Partial',&
            &' Charges :',&
            &//,4x,'Atom',14x,'Base Value',7x,'Actual',/)
         end if
         write (iout,20)  k,pchg0(i),pchg(i)
20       format (i8,9x,2f14.5)
      end if
   end do
!
!     alter monopoles and charge penetration for charge flux
!
   header = .true.
   do ii = 1, npolebloc
      i = poleglob(ii)
      k = ipole(i)
      pole(1,i) = mono0(i) + pdelta(k)
      if (use_chgpen)  pval(i) = pval0(i) + pdelta(k)
      if (debug .and. pdelta(k).ne.0.0d0) then
         if (header) then
            header = .false.
            write (iout,30)
30          format (/,' Charge Flux Modification of Atomic',&
            &' Monopoles :',&
            &//,4x,'Atom',14x,'Base Value',7x,'Actual',/)
         end if
         write (iout,40)  k,mono0(i),pole(1,i)
40       format (i8,9x,2f14.5)
      end if
   end do
!
!     perform deallocation of some local arrays
!
   deallocate (pdelta)
   return
end
!
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine bndchg  --  charge flux bond stretch coupling  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "bndchg" computes modifications to atomic partial charges or
!     monopoles due to bond stretch using a charge flux formulation
!
!
!> @brief 
!> "bndchg" computes modifications to atomic partial charges or
!> monopoles due to bond stretch using a charge flux formulation
!> @param no params
subroutine bndchg (pdelta)
   use sizes
   use atmlst
   use atoms
   use bond
   use bound
   use cflux
   implicit none
   integer i,ia,ib
   integer ibond
   real*8 xab,yab,zab
   real*8 rab,rab0
   real*8 pb,dq
   real*8 pdelta(*)
!
!
!     loop over all the bond distances in the system
!
   do ibond = 1, nbondloc
      i = bndglob(ibond)
      ia = ibnd(1,i)
      ib = ibnd(2,i)
      pb = bflx(i)
!
!     compute the bond length value for the current bond
!
      xab = x(ia) - x(ib)
      yab = y(ia) - y(ib)
      zab = z(ia) - z(ib)
      if (use_polymer)  call image (xab,yab,zab)
      rab = sqrt(xab*xab + yab*yab + zab*zab)
!
!     find the charge flux increment for the current bond
!
      rab0 = bl(i)
      dq = pb * (rab-rab0)
      pdelta(ia) = pdelta(ia) - dq
      pdelta(ib) = pdelta(ib) + dq
   end do
   return
end
!
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine angchg  --  charge flux angle bend coupling  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "angchg" computes modifications to atomic partial charges or
!     monopoles due to angle bending using a charge flux formulation
!
!
!> @brief 
!> "angchg" computes modifications to atomic partial charges or
!> monopoles due to angle bending using a charge flux formulation
!> @param no params
subroutine angchg (pdelta)
   use sizes
   use angle
   use atmlst
   use atoms
   use bond
   use bound
   use cflux
   use math
   implicit none
   integer i,ia,ib,ic
   integer iangle
   real*8 angle1,eps
   real*8 rab,rcb
   real*8 xia,yia,zia
   real*8 xib,yib,zib
   real*8 xic,yic,zic
   real*8 xab,yab,zab
   real*8 xcb,ycb,zcb
   real*8 dot,cosine
   real*8 pa1,pa2
   real*8 pb1,pb2
   real*8 theta0
   real*8 rab0,rcb0
   real*8 dq1,dq2
   real*8 pdelta(*)
!
!
!     loop over all the bond angles in the system
!
   eps = 0.0001d0
   do iangle = 1, nangleloc
      i = angleglob(iangle)
      ia = iang(1,i)
      ib = iang(2,i)
      ic = iang(3,i)
      pa1 = aflx(1,i)
      pa2 = aflx(2,i)
      pb1 = abflx(1,i)
      pb2 = abflx(2,i)
!
!     calculate the angle values and included bond lengths
!
      xia = x(ia)
      yia = y(ia)
      zia = z(ia)
      xib = x(ib)
      yib = y(ib)
      zib = z(ib)
      xic = x(ic)
      yic = y(ic)
      zic = z(ic)
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
      rab = sqrt(max(xab*xab+yab*yab+zab*zab,eps))
      rcb = sqrt(max(xcb*xcb+ycb*ycb+zcb*zcb,eps))
      dot = xab*xcb + yab*ycb + zab*zcb
      cosine = dot / (rab*rcb)
      cosine = min(1.0d0,max(-1.0d0,cosine))
      angle1 = radian * acos(cosine)
!
!     find the charge flux increment for the current angle
!
      theta0 = anat(i)
      rab0 = bl(balist(1,i))
      rcb0 = bl(balist(2,i))
      dq1 = pb1*(rcb-rcb0) + pa1*(angle1-theta0)/radian
      dq2 = pb2*(rab-rab0) + pa2*(angle1-theta0)/radian
      pdelta(ia) = pdelta(ia) + dq1
      pdelta(ib) = pdelta(ib) - dq1 - dq2
      pdelta(ic) = pdelta(ic) + dq2
   end do
   return
end
