c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine ebond  --  bond stretch potential energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "ebond" calculates the bond stretching energy
c
c
#include "tinker_precision.h"
      module ebond_inl
        contains
#include "image.f.inc"
      end module

      subroutine ebond
      use atmlst
      use atoms
      use bndpot
      use bond
      use bound
      use ebond_inl
      use energi
      use group
      use nvshmem
      use tinheader
      use timestat ,only: timer_enter,timer_exit,timer_ebond
      use usage
      implicit none
      integer i,ia,ib,ibond
      integer ind,ipe
      real(t_p) e,ideal,force
      real(t_p) expterm,bde
      real(t_p) dt,dt2
      real(t_p) xab,yab,zab,rab
      logical proceed
c
c
c     zero out the bond stretching energy
c
      call timer_enter( timer_ebond )
      eb = 0.0_re_p
c
c     calculate the bond stretching energy term
c
!$acc parallel loop default(present) present(eb) async
      do ibond = 1, nbondloc
         i     = bndglob(ibond)
#ifdef USE_NVSHMEM_CUDA
         ipe   = (i-1)/nbond_pe
         ind   = mod((i-1),nbond_pe) +1
         ia    = d_ibnd(ipe)%pel(1,ind)
         ib    = d_ibnd(ipe)%pel(2,ind)
         ideal = d_bl  (ipe)%pel(ind)
         force = d_bk  (ipe)%pel(ind)
#else
         ia    = ibnd(1,i)
         ib    = ibnd(2,i)
         ideal = bl(i)
         force = bk(i)
#endif
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
            if (use_polymer)  call image_inl (xab,yab,zab)
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
c     increment the total bond stretching energy
c
            eb = eb + e
         end if
      end do
      call timer_exit( timer_ebond )
      end
