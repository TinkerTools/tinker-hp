c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine beeman  --  Beeman molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "beeman" performs a single molecular dynamics time step
c     via the Beeman multistep recursion formula; uses original
c     coefficients or Bernie Brooks' "Better Beeman" values
c
c     literature references:
c
c     D. Beeman, "Some Multistep Methods for Use in Molecular
c     Dynamics Calculations", Journal of Computational Physics,
c     20, 130-139 (1976)
c
c     B. R. Brooks, "Algorithms for Molecular Dynamics at Constant
c     Temperature and Pressure", DCRT Report, NIH, April 1988
c
c
#include "tinker_precision.h"
      subroutine beeman (istep,dt)
      use atmtyp
      use atoms
      use cutoff
      use domdec
      use energi
      use freeze
      use mdstuf
      use moldyn
      use timestat
      use units
      use usage
      use mpi
      implicit none
      integer i,j,istep,iglob
      real(t_p) dt,dt_x,factor
      real(r_p) etot,eksum,epot
      real(t_p) temp,pres
      real(t_p) part1,part2
      real(r_p) ekin(3,3)
      real(t_p) stress(3,3)
      real(r_p) time0,time1
      real(r_p), allocatable :: derivs(:,:)
c
c     set time values and coefficients for Beeman integration
c
      factor = real(bmnmix,t_p)
      dt_x = dt / factor
      part1 = 0.5_ti_p*factor + 1.0_ti_p
      part2 = part1 - 2.0_ti_p
c
c     store the current atom positions, then find half-step
c     velocities and full-step positions via Beeman recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = v(j,iglob) + (part1*a(j,iglob)-
     $           aalt(j,iglob))*dt_x
            end do
            xold(iglob) = x(iglob)
            yold(iglob) = y(iglob)
            zold(iglob) = z(iglob)
            x(iglob) = x(iglob) + v(1,iglob)*dt
            y(iglob) = y(iglob) + v(2,iglob)*dt
            z(iglob) = z(iglob) + v(3,iglob)*dt
            end if
      end do
!$acc update device(x(:),y(:),z(:))
      factor = real(bmnmix,t_p)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
      call reassign
c
c     -> reciprocal space
      call reassignpme(.false.)
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt)
c
c     communicate positions
c
      call commpos
      call commposrec
c
      allocate (derivs(3,nbloc))
      derivs = 0_ti_p
c
      call reinitnl(istep)
c
      call mechanicstep(istep)
c
      call allocstep
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     MPI : get total energy
c
      call reduceen(epot)
c
c     make half-step temperature and pressure corrections
c
      call temper2 (temp)
      call pressure2 (epot,temp)
c
c     communicate forces
c
      call commforces(derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Beeman recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               aalt(j,iglob) = a(j,iglob)
               a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
               v(j,iglob) = v(j,iglob) +
     $          (part2*a(j,iglob)+aalt(j,iglob))*dt_x
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     make full-step temperature and pressure corrections
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + esum
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
      call mdrest (istep)
!$acc update device(v,a)
      return
      end
