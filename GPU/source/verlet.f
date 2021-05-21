c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine verlet  --  Verlet molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "verlet" performs a single molecular dynamics time step
c     via the velocity Verlet multistep recursion formula
c
c
#include "tinker_precision.h"

      subroutine verlet (istep,dt)
      use atmtyp
      use atomsMirror
      use cutoff
      use domdec
      use deriv  ,only: info_forces,prtEForces,cDef
      use energi ,only: info_energy
      use freeze
      use inform
      use moldyn
      use timestat
      use units
      use usage
      use utils  ,only:set_to_zero1m
      use utilgpu,only:prmem_requestm,rec_queue
      use mpi
      use sizes
      use integrate_ws
      implicit none
      integer i,j,istep
      integer iglob
      real(r_p) dt,dt_2
      real(r_p) time0,time1
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5_re_p * dt

      if (istep.eq.1) then
!$acc enter data create(etot,epot,eksum,temp,pres,ekin,stress)
      end if

!$acc data present(etot,epot,eksum,temp,pres,ekin,stress)
!$acc&     present(x,y,z,xold,yold,zold,v,a,mass,glob,use)
c
c     store the current atom positions, then find half-step
c     velocities and full-step positions via Verlet recursion
c
      if (use_rattle) call save_atoms_pos
!$acc parallel loop async default(present)
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
!$acc loop seq
            do j = 1, 3
               v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
            end do
            x(iglob) = x(iglob) + v(1,iglob)*dt
            y(iglob) = y(iglob) + v(2,iglob)*dt
            z(iglob) = z(iglob) + v(3,iglob)*dt
         end if
      end do
c
c     Reassign the particules that have changed of domain
c
c     -> real space
      call reassign
c
c     -> reciprocal space
      call reassignpme(.false.)
c
c     communicate positions
c
      call commpos
      call commposrec
      call reCast_position
c
      call prmem_requestm(derivs,3,nbloc,async=.true.)
      call set_to_zero1m(derivs,3*nbloc,rec_queue)
c
      call reinitnl(istep)
c
      call mechanicstep(istep)

      call allocstep
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt)
c
c     rebuild the neighbor list
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
c     Debug print information
c
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  call info_forces(cDef)
      if (deb_Atom)   call info_minmax_pva
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Verlet recursion
c
!$acc parallel loop async collapse(2) present(derivs)
      do i = 1, nloc
         do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               a(j,iglob) = -convert * derivs(j,i) / mass(iglob)
               v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
            end if
         end do
      end do
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     make full-step temperature and pressure corrections
c
      call temper   (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
!$acc serial async
      etot = eksum + epot
!$acc end serial
c
c     compute statistics and save trajectory for this step
c
      if (nproc.eq.1.and.tinkerdebug.eq.64) call prtEForces(derivs,epot)
      call mdsave (istep,dt,epot)
      call mdrestgpu (istep)
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
!$acc end data
      end
