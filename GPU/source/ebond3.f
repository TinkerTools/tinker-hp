c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ebond3  --  bond stretch energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ebond3" calculates the bond stretching energy; also
c     partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      subroutine ebond3
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bndpot
      use bond
      use bound
      use domdec
      use energi
      use group
      use inform
      use iounit
      use tinheader
      use usage
      implicit none
      integer i,ia,ib,ialoc,ibloc,ibond
      real(t_p) e,ideal,force
      real(t_p) expterm,bde
      real(t_p) dt,dt2
      real(t_p) xab,yab,zab,rab
      logical proceed
      logical header,huge
c
c
c     zero out the bond energy and partitioning terms
c
      neb = 0
      eb = 0.0_ti_p
      aeb = 0.0_ti_p
      header = .true.
c
c     calculate the bond stretching energy term
c
      do ibond = 1, nbondloc
         i = bndglob(ibond)
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ialoc = loc(ia)
         ibloc = loc(ib)
         ideal = bl(i)
         force = bk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            if (use_polymer)  call image (xab,yab,zab)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
c
c     harmonic potential uses Taylor expansion of Morse potential
c     through the fourth power of the bond length deviation
c
            if (bndtyp .eq. 'HARMONIC') then
               dt2 = dt * dt
               e = bndunit * force * dt2 * (1.0_ti_p+cbnd*dt+qbnd*dt2)
c
c     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            else if (bndtyp .eq. 'MORSE') then
               expterm = exp(-2.0_ti_p*dt)
               bde = 0.25_ti_p * bndunit * force
               e = bde * (1.0_ti_p-expterm)**2
            end if
c
c     increment the total bond energy and partition between atoms
c
            neb = neb + 1
            eb = eb + e
            aeb(ialoc) = aeb(ialoc) + 0.5_ti_p*e
            aeb(ibloc) = aeb(ibloc) + 0.5_ti_p*e
c
c     print a message if the energy of this interaction is large
c
            huge = (e .gt. 5.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Bond Stretching',
     &                       ' Interactions :',
     &                    //,' Type',14x,'Atom Names',22x,'Ideal',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ideal,rab,e
   20          format (' Bond',6x,2(i7,'-',a3),13x,2f10.4,f12.4)
            end if
         end if
      end do
      return
      end
