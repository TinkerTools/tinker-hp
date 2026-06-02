!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine bbk  --  BBK Langevin molecular dynamics step  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "bbk" performs a single molecular dynamics time step
!     via the BBK recursion formula
!
!     literature reference:
!
!     A. Brünger, C. L. Brooks III, M. Karplus, Stochastic boundary
!     conditions fro molecular dynamics simulations of ST2 water.
!     Chem. Phys. Letters, 1984, 105 (5) 495-500
!
!
#include "tinker_macro.h"
module bbk_inl
contains
#include "convert.inc.f90"
end module

subroutine bbk (istep,dt)
   use atmtyp
   use atomsMirror
   use bath
   use bbk_inl
   use cutoff
   use domdec
   use deriv    ,only:info_forces,cBond,cNBond,comm_forces,tdes_l&
      &,dr_stride,de_tot,derivx
   use energi
   use freeze
   use inform
   use langevin
   use mdstuf
   use mdstuf1
   use moldyn
   use mpi
   use random_mod
   use timestat
   use tinheader
   use units
   use usage
   use utilgpu
   implicit none
   integer i,j,istep,iglob
   real(r_p) dt,dt_2
   real(r_p) a1,a2
   mdyn_rtyp dx
!
!     set some time values for the dynamics integration
!
   dt_2 = 0.5_re_p * dt
!
!     set time values and coefficients for BBK integration
!
   a1 = 1-dt_2*gamma
   a2 = 0.5_re_p*sqrt(2*boltzmann*kelvin*gamma*dt)
!
!     find half step velocities and full step positions via BBK recursion
!
   if (use_rattle) call save_atoms_pos
!$acc parallel loop collapse(2) default(present) async
   do i = 1,nloc; do j = 1,3
         iglob = glob(i)
         if (use(iglob)) then
            v(j,iglob) = a1*v(j,iglob) +  dt_2*a(j,iglob)&
               &+ a2*Rn(j,i)/sqrt(mass(iglob))
         end if
      end do; end do
   call integrate_pos(dt)
!
!     Reassign the particules that have changed of domain
!
!     -> real space
   call reassign
!
!     -> reciprocal space
   call reassignpme(.false.)
!
!     get constraint-corrected positions and half-step velocities
!
   if (use_rattle)  call rattle (dt)
!
!     communicate positions
!
   call commpos
   call commposrec
!
   call reinitnl(istep)
!
   call mechanicstep(istep)
!
   call allocstep
!
!     rebuild the neighbor lists
!
   if (use_list) call nblist(istep)
!
!     get the potential energy and atomic forces
!
   call gradient (epot,derivs)
!
!     MPI : get total energy
!
   call reduceen(epot)
!
!     make half-step temperature and pressure corrections
!
   call temper2 (temp)
   call pressure2 (epot,temp)
!
!     communicate forces
!
   call comm_forces(derivs)
!
!     aMD/GaMD contributions
!
   call aMD (derivs,epot)
!
!     Debug print information
!
   if (deb_Energy) call info_energy(rank)
   if (deb_Force)  call info_forces(cNBond)
   if (deb_Atom)   call info_minmax_pva
   if (abort)      __TINKER_FATAL__
!
!     compute random part
!
   call prmem_request(Rn,3,nloc+1,async=.false.)
#ifdef _OPENACC
   call normalgpu(Rn,3*nloc)
#endif
   if (host_rand_platform) then
      call normalvec(Rn,3*nloc)
!$acc update device(Rn) async
   end if
!
!     use Newton's second law to get the next accelerations;
!     find the full-step velocities using the BBK recursion
!
   if (tdes_l) then
!$acc parallel loop collapse(2) async default(present)
      do i = 1,nloc; do j = 1,3
            iglob = glob(i)
            if (i.le.nloc) then
               if (useAll.or.use(iglob)) then
                  dx         = derivx(i + (j-1)*dr_stride )
                  a(j,iglob) = -convert * mdr2md(dx)/mass(iglob)
                  v(j,iglob) = (v(j,iglob) + dt_2*a(j,iglob)&
                     &+ a2*Rn(j,i)/sqrt(mass(iglob)))/(1.0+dt*gamma/2)
                  derivx(i + (j-1)*dr_stride ) = 0
               end if
            else
               derivx(i + (j-1)*dr_stride ) = 0
            end if
         end do; end do
   else
!$acc parallel loop collapse(2) async default(present)
      do i = 1,nloc; do j = 1,3
            iglob = glob(i)
            if (useAll.or.use(iglob)) then
               a(j,iglob) = -convert * mdr2md(de_tot(j,i))/mass(iglob)
               v(j,iglob) = (v(j,iglob) + dt_2*a(j,iglob)&
                  &+ a2*Rn(j,i)/sqrt(mass(iglob)))/(1.0+dt*gamma/2)
               de_tot(j,i)= 0
            end if
         end do; end do
   end if
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle)  call rattle2 (dt)
!
!     make full-step temperature and pressure corrections
!
   call temper (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
!
!     total energy is sum of kinetic and potential energies
!
   if (calc_e) then
!$acc serial present(etot,eksum,esum) async
      etot = eksum + esum
!$acc end serial
   end if
!
!     compute statistics and save trajectory for this step
!
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)
   call mdsave (istep,dt,epot)
   call mdrestgpu (istep)
end
