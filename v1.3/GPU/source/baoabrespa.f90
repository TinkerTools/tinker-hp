!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  subroutine baoab  --  baoab Langevin molecular dynamics step  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     "baoab" performs a single molecular dynamics time step
!     via the respa-baoab recursion formula
!
!     literature references:
!
!     Efficient molecular dynamics using geodesic integration
!     and solvent-solute splitting B. Leimkuhler and C. Matthews,
!     Proceedings of the Royal Society A, 472: 20160138, 2016
!
!     D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-
!     Time-Step Molecular Dynamics Algorithm for Macromolecules",
!     Journal of Physical Chemistry, 98, 6885-6892 (1994)
!
!     X. Qian and T. Schlick, "Efficient Multiple-Time-Step Integrators
!     with Distance-Based Force Splitting for Particle-Mesh-Ewald
!     Molecular Dynamics Simulations", Journal of Chemical Physics,
!     115, 4019-4029 (2001)
!
!
#include "tinker_macro.h"
subroutine baoabrespa (istep,dt)
   use mlpot
   use atmtyp
   use atomsMirror
   use bath
   use boxes,only: volbox
   use cutoff
   use domdec
   use deriv,only:info_forces,cBond,cNBond,ftot_l,comm_forces&
      &, zero_forces_rec
   use energi
   use freeze
   use group
   use inform
   use mdstuf
   use mdstuf1
   use moldyn
   use mpi
   use potent
   use random_mod
   use timestat
   use tortor
   use units
   use uprior ,only: use_pred
   use usage
   use utils,only: set_to_zero1m
   use utilgpu
   use utilbaoab
   use virial
   use spectra
   use utilbaoab
   implicit none
   real(r_p), intent(in) :: dt
   integer, intent(in) :: istep
   integer i,j,k,iglob
   real(r_p) dt_x,factor
   real(r_p) dta,dta_2,dt_2
   real(r_p) part1,part2
   real(r_p), allocatable, save :: dip(:,:),dipind(:,:)
   real(8) time0,time1
   logical save_mlpot
   logical save_pred
!
!     set some time values for the dynamics integration
!
   dt_2  = 0.5_re_p * dt
   dta   = dt / dshort
   dta_2 = 0.5_re_p * dta

   if (istep.eq.1) then
      allocate(dip(3,nalt),dipind(3,nalt))
!$acc enter data create(dip,dipind)
      if (use_piston) pres = atmsph
      call set_langevin_thermostat_coeff(dta)
   end if
   if(use_ml_embedding) use_mlpot=.FALSE.
!
!     find quarter step velocities and half step positions via baoab recursion
!
   call integrate_vel( a,dt_2 )
!
   if (use_rattle) call rattle2(dt_2)
   if (use_rattle) call save_atoms_pos
   !initialize virial from fast-evolving potential energy terms
   if (use_virial) call zero_virial(viralt)

   if (use_piston) then
      call apply_a_piston(dt_2,-1,.false.)
      call apply_o_piston(dt)
      if (use_rattle) call rattle (dta_2)
      if (use_rattle) call rattle2(dta_2)
   end if
!
!     respa inner loop
!
   do stepfast = 1, nalt
!
!     find fast-evolving velocities and positions via BAOAB recursion
!
      call integrate_vel( aalt,dta_2 )
!
      if (use_rattle) call rattle2 (dta_2)
      if (use_rattle) call save_atoms_pos
!
      call integrate_pos( dta_2 )
!
      if (use_rattle) call rattle (dta_2)
      if (use_rattle) call rattle2(dta_2)

      call apply_langevin_thermostat(dta)
!
      if (use_rattle) call rattle2(dta_2)
      if (use_rattle) call save_atoms_pos
!
      call integrate_pos( dta_2 )
!
      if (use_rattle) call rattle (dta_2)
      if (use_rattle) call rattle2(dta_2)
!
!       Reassign the particules that have changed of domain
!
!       -> real space
      call reassignrespa(stepfast,nalt)
!
!       communicate positions
!
      call commposrespa(stepfast.ne.nalt)
!
      if (.not.ftot_l) then
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m(derivs,3*nbloc,rec_queue)
      end if
!
      call mechanicsteprespa(istep,.true.)
      call allocsteprespa(.true.)
!
!     get the fast-evolving potential energy and atomic forces
!
      call gradfast ( ealt,derivs )
!
!       communicate forces
!
      call comm_forces( derivs,cBond )
!
!       MPI : get total energy
!
      call reduceen(ealt)
!
!     Debug print information
!
      if (deb_Energy)call info_energy(rank)
      if (deb_Force) call info_forces(cBond)
      if (deb_Atom)  call info_minmax_pva
!
!     use Newton's second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the Verlet recursion
!
      call integrate_vel( derivs,aalt,dta_2 )
!
      if (use_rattle)  call rattle2 (dta_2)
!
!     increment average virial from fast-evolving potential terms
!
      if (use_virial) then
!$acc parallel loop collapse(2) async present(vir,viralt)
         do i = 1,3; do j = 1,3
               viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
            end do; end do
      end if

      if (ir .and. stepfast<nalt) then
!          call rotpolegpu
         call compute_dipole(dip(:,stepfast)&
            &,dipind(:,stepfast),.FALSE.)
      end if
   end do

   if (use_piston) then
      if (use_rattle)  call save_atoms_pos
      call apply_a_piston(dt_2,-1,.false.)
      if (use_rattle) call rattle (dta_2)
      if (use_rattle) call rattle2(dta_2)
   end if
!
!     Reassign the particules that have changed of domain
!
!     -> real space
   call reassignrespa(nalt,nalt)
!
!     -> reciprocal space
   call reassignpme(.false.)
!
!     communicate positions
!
   call commposrespa(.false.)
   call commposrec
!
   call reinitnl(istep)
!
   call mechanicsteprespa(istep,.false.)

   call allocsteprespa(.false.)
!
!     rebuild the neighbor lists
!
   if (use_list) call nblist(istep)
!
   if (.not.ftot_l) then
      call prmem_requestm(derivs,3,nbloc,async=.true.)
      call set_to_zero1m(derivs,3*nbloc,rec_queue)
   end if
!
!     get the slow-evolving potential energy and atomic forces
!
   call gradslow (epot,derivs)
!
!     if necessary, communicate some forces
!
   call comm_forces( derivs,cNBond )
!
!     MPI : get total energy
!
   call reduceen(epot)

   if(use_ml_embedding) then
!$acc parallel loop collapse(2) async
      do i = 1, 3; do j = 1, 3
            viralt(j,i) = vir(j,i) + viralt(j,i)
         end do; end do
!     COMPUTE ML DELTA CONTRIBUTION (ml_embedding_mode)
      use_mlpot=  .TRUE.
      call zero_forces_rec
      save_pred = use_pred
      use_pred  = .FALSE.
      if(use_embd_potoff) then
         call gradembedding2 (eml,derivs)
      else
         call gradient (eml,derivs)
      endif
      use_pred  = save_pred
      call reduceen(eml)
      call comm_forces(derivs)
!$acc serial async
      epot = epot+eml
!$acc end serial
   endif
!
!     use Newton's second law to get the next accelerations;
!     find the full-step velocities using the BAOAB recursion
!
   call integrate_vel( derivs,a,dt_2 )
!
!     Debug print information
!
   if (deb_Energy) call info_energy(rank)
   if (deb_Force)  call info_forces(cNBond)
   if (deb_Atom)   call info_minmax_pva
   if (abort)   __TINKER_FATAL__
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle) call rattle2 (dt_2)
!
!     total potential and virial from sum of fast and slow parts
!
   if (calc_e.or.use_virial) then
!$acc serial async present(epot,ealt,vir,viralt)
      epot = epot + ealt
      do i = 1,3; do j = 1,3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do; end do
!$acc end serial
   end if
!
   call temper   (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
   call pressure2 (epot,temp)
   if(ir) then
      call compute_dipole(dip(:,nalt)&
         &,dipind(:,nalt),full_dipole)
!$acc update host(dip,dipind) async
!$acc wait
      call save_dipole_respa(dip,dipind)
   endif
!
!     total energy is sum of kinetic and potential energies
!
   if (calc_e) then
!$acc serial async present(etot,epot,eksum)
      etot = eksum + epot
!$acc end serial
   end if
!
!     compute statistics and save trajectory for this step
!
   call mdsave (istep,dt,epot)
   call mdrestgpu (istep)
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)

   if (use_piston) then
!$acc update host(pres,stress) async
!$acc wait
      call apply_b_piston(dt,pres,stress)
      call ddpme3dnpt(1.0_re_p,istep)
   end if

end
