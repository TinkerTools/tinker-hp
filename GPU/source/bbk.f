c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine bbk  --  BBK Langevin molecular dynamics step  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "bbk" performs a single molecular dynamics time step
c     via the BBK recursion formula
c
c     literature reference:
c
c     A. Brünger, C. L. Brooks III, M. Karplus, Stochastic boundary 
c     conditions fro molecular dynamics simulations of ST2 water. 
c     Chem. Phys. Letters, 1984, 105 (5) 495-500
c
c
#include "tinker_precision.h"
      subroutine bbk (istep,dt)
      use atmtyp
      use atoms
      use bath
      use cutoff
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use random_mod
      use timestat
      use units
      use usage
      use mpi
      implicit none
      integer i,j,istep,iglob
      real(t_p) dt,dt_2
      real(t_p) etot,eksum,epot
      real(t_p) temp,pres
      real(t_p) a1,a2
      real(t_p) ekin(3,3)
      real(t_p) stress(3,3)
      real(t_p) time0,time1
      real(t_p), allocatable :: derivs(:,:)
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
c
c     set time values and coefficients for BBK integration
c
      a1 = 1-dt_2*gamma
      a2 = 0.5*sqrt(2*boltzmann*kelvin*gamma*dt)
c
c     find half step velocities and full step positions via BBK recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = a1*v(j,iglob) 
     $       +  dt_2*a(j,iglob)
     $       +  a2*Rn(j,i)/sqrt(mass(iglob))
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
c     aMD/GaMD contributions
c
      call aMD (derivs,epot)
c
c     compute random part
c
      deallocate (Rn)
      allocate (Rn(3,nloc))
      do i = 1, nloc
        do j = 1, 3
          Rn(j,i) = normal()
        end do
      end do
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BBK recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
               v(j,iglob) = (v(j,iglob) +
     $          0.5*dt*a(j,iglob) + a2*Rn(j,i)/sqrt(mass(iglob)))/
     $          (1+dt*gamma/2)
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
