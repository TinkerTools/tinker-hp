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
#include "tinker_macro.h"
      module bbk_inl
      contains
#include "convert.f.inc"
      end module

      subroutine bbk (istep,dt)
      use atmtyp
      use atomsMirror
      use bath
      use bbk_inl
      use cutoff
      use domdec
      use deriv    ,only:info_forces,cBond,cNBond,comm_forces,tdes_l
     &             ,dr_stride,de_tot,derivx
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
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5_re_p * dt
c
c     set time values and coefficients for BBK integration
c
      a1 = 1-dt_2*gamma
      a2 = 0.5_re_p*sqrt(2*boltzmann*kelvin*gamma*dt)
c
c     find half step velocities and full step positions via BBK recursion
c
      if (use_rattle) call save_atoms_pos
!$acc parallel loop collapse(2) default(present) async
      do i = 1,nloc; do j = 1,3
         iglob = glob(i)
         if (use(iglob)) then
            v(j,iglob) = a1*v(j,iglob) +  dt_2*a(j,iglob)
     &                 + a2*Rn(j,i)/sqrt(mass(iglob))
         end if
      end do; end do
      call integrate_pos(dt)
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
      call comm_forces(derivs)
c
c     aMD/GaMD contributions
c
      call aMD (derivs,epot)
c
c     Debug print information
c
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  call info_forces(cNBond)
      if (deb_Atom)   call info_minmax_pva
      if (abort)      __TINKER_FATAL__
c
c     compute random part
c
      call prmem_request(Rn,3,nloc+1,async=.false.)
#ifdef _OPENACC
      call normalgpu(Rn,3*nloc)
#endif
      if (host_rand_platform) then
         call normalvec(Rn,3*nloc)
!$acc update device(Rn) async
      end if
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BBK recursion
c
      if (tdes_l) then
!$acc parallel loop collapse(2) async default(present)
      do i = 1,nloc; do j = 1,3
         iglob = glob(i)
         if (i.le.nloc) then
         if (useAll.or.use(iglob)) then
            dx         = derivx(i + (j-1)*dr_stride )
            a(j,iglob) = -convert * mdr2md(dx)/mass(iglob)
            v(j,iglob) = (v(j,iglob) + dt_2*a(j,iglob)
     &                 + a2*Rn(j,i)/sqrt(mass(iglob)))/(1.0+dt*gamma/2)
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
            v(j,iglob) = (v(j,iglob) + dt_2*a(j,iglob)
     &                 + a2*Rn(j,i)/sqrt(mass(iglob)))/(1.0+dt*gamma/2)
            de_tot(j,i)= 0
         end if
      end do; end do
      end if
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
      if (calc_e) then
!$acc serial present(etot,eksum,esum) async
         etot = eksum + esum
!$acc end serial
      end if
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
      call mdrestgpu (istep)
      end
