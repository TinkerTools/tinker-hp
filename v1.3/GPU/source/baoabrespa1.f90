
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################################
!     ##                                                                               ##
!     ##  subroutine baoabrespa1  --  baoab r-RESPA1 Langevin molecular dynamics step  ##
!     ##                                                                               ##
!     ###################################################################################
!
!
!     "baoabrespa1" performs a single multiple time step molecular dynamics
!     step using the reversible reference system propagation algorithm
!     (r-RESPA) via a BAOAB core with the potential split into fast-
!     intermediate and slow-evolving portions
!
!     literature references:
!
!     Pushing the Limits of Multiple-Time-Step Strategies
!     wfor Polarizable Point Dipole Molecular Dynamics
!     L Lagardère, F Aviat, JP Piquemal
!     The journal of physical chemistry letters 10, 2593-2599
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
!     Ruhong Zhou, Edward Harder, Huafeng Xu and B. J. Berne,
!     "Efficient multiple time step method for use with Ewald and
!     particle mesh Ewald for large biomolecular systems",
!     J. Chem. Phys. 115, 2348-2358 (2001)
!
!
#include "tinker_macro.h"
subroutine baoabrespa1(istep,dt)
   use mlpot
   use atmtyp
   use atomsMirror
   use cutoff
   use domdec
   use deriv
   use energi
   use freeze
   use group
   use inform
   use mdstuf1
   use moldyn
   use mpi
   use potent
   use timestat
   use units
   use uprior  ,only: use_pred
   use usage
   use utilgpu
   use utilbaoab
   use utils
   use virial
   implicit none
   integer i,j,iglob
   integer istep
   real(r_p) dt,dt_2
   real(r_p) dta,dta_2,dts
   logical save_pred
!
!     set some time values for the dynamics integration
!
   dt_2  = 0.5_re_p * dt
   dta   = dt / dinter
   dta_2 = 0.5_re_p * dta
   dts   = dta / dshort
   if (istep.eq.1) call set_langevin_thermostat_coeff(dts)
!
   if (use_ml_embedding) use_mlpot=.FALSE.
!
!     store the current atom positions, then find half-step
!     velocities via BAOAB recursion
!
   call integrate_vel(a,dt_2)
!
   if (use_rattle) call rattle2(dt_2)
!
!     find intermediate-evolving velocities and positions via BAOAB recursion
!
   call baoabrespaint1(ealt,viralt,dta,dts)
!
!     Reassign the particules that have changed of domain
!
!     -> real space
!      call reassignrespa(.false.,nalt,nalt)
!
!     -> reciprocal space
   call reassignpme(.false.)
!
!     communicate positions
!
!     call commposrespa(.false.)
   call commposrec
!
!
   call reinitnl(istep)
!
   call mechanicsteprespa1(istep,2)

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
!     communicate some forces
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
!     COMPUTE ML DELTA CONTRIBUTION (ml_embedding_mode=2)
      use_mlpot=  .TRUE.
      call zero_forces_rec
      save_pred = use_pred
      use_pred = .FALSE.
      call gradient (eml,derivs)
      use_pred = save_pred
      call reduceen(eml)
      call comm_forces(derivs)
!$acc serial async
      epot = epot+eml
!$acc end serial
   endif
!
!     make half-step temperature and pressure corrections
!
!     call temper2 (temp)
!     call pressure2 (epot,temp)
!
!     Debug print information
!
   if (deb_Energy) call info_energy(rank)
   if (deb_Force) then
      call info_forces(cNBond)
      !if (ftot_l) call minmaxone(de_tot,3*nloc,'de_tot')
   end if
   if (deb_Atom)   call info_minmax_pva
!
!     use Newton's second law to get the slow accelerations;
!     find full-step velocities using BAOAB recursion
!
   call integrate_vel(derivs,a,dt_2)
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle)  call rattle2 (dt)
!
!     total potential and virial from sum of fast and slow parts
!
   if (calc_e.or.use_virial) then
!$acc serial default(present) present(epot,ealt) async
      epot = epot + ealt
      do i = 1,3; do j = 1,3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do; end do
!$acc end serial
      if(rank.eq.0) call chk_energy_fluct(epot,ealt,abort)
   end if
!
!     make full-step temperature and pressure corrections
!
   call temper   (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
   call pressure2 (epot,temp)
!
!     total energy is sum of kinetic and potential energies
!
   if (calc_e) then
!$acc serial present(etot,eksum,epot) async
      etot = eksum + epot
!$acc end serial
   end if

   ! Fatal Instructions
   if(abort)       call emergency_save
   if(abort)       call fatal
!
!     compute statistics and save trajectory for this step
!
   call mdsave (istep,dt,epot)
   call mdrestgpu (istep)
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)
end
!
!     subroutine baoabrespaint1 :
!     find intermediate-evolving velocities and positions via BAOAB recursion
!
subroutine baoabrespaint1(ealt,viralt,dta,dts)
   use mlpot
   use atmtyp
   use atomsMirror
   use cutoff
   use deriv
   use domdec
   use deriv
   use energi
   use freeze
   use inform
   use mdstuf1 ,only: derivs,ealt2,viralt2
   use moldyn
   use mutant
   use potent
   use timestat
   use units
   use uprior  ,only: use_pred
   use usage
   use utils
   use utilgpu
   use virial
   use mpi
   implicit none
   integer i,j,k,iglob
   real(r_p) dta,dta_2,dts
   real(r_p) ealt,ealtml
   real(r_p) viralt(3,3)
   logical save_pred

   dta_2 = 0.5_re_p * dta
!
!     initialize virial from fast-evolving potential energy terms
!
   if (use_virial) call zero_virial(viralt)

   do stepint = 1, nalt
      call integrate_vel(aalt,dta_2)
!
      if (use_rattle)  call rattle2 (dta_2)
!
!     find fast-evolving velocities and positions via BAOAB recursion
!
      call baoabrespafast1(ealt2,viralt2,dts)
!
!       Reassign the particules that have changed of domain
!
!       -> real space
!       call reassignrespa(.false.,nalt,nalt)
!
!      communicate positions
!
!        call commposshort(.false.)
!
      call mechanicsteprespa1(-1,1)
      call allocsteprespa(.false.)
!
      if (.not.ftot_l) then
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m (derivs,3*nbloc,rec_queue)
      end if
!
!     get the fast-evolving potential energy and atomic forces
!
      call gradint(ealt,derivs)
!
!     communicate forces
!
      call comm_forces( derivs,cSNBond )
!
!     MPI : get total energy
!
      call reduceen(ealt)
!
!     Debug print information
!
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  then
         call info_forces(cSNBond)
         !if (ftot_l) call minmaxone(de_tot,3*nloc,'de_tot')
      end if
      if (deb_Atom)   call info_minmax_pva
      if (abort)      call emergency_save
      if (abort)      call fatal

!      if(use_ml_embedding) then
!!$acc parallel loop collapse(2) async
!        do i = 1, 3; do j = 1, 3
!          viralt2(j,i) = vir(j,i) + viralt2(j,i)
!        end do; end do
!c     COMPUTE ML DELTA CONTRIBUTION (ml_embedding_mode=2)
!        use_mlpot= .TRUE.
!!$acc data create(ealtml)
!        call resetForcesRec
!        save_pred = use_pred
!        use_pred  = .FALSE.
!        call gradient (ealtml,derivs)
!        use_pred  = save_pred
!        call reduceen(ealtml)
!        call comm_forces(derivs)
!!$acc serial async present(ealt,ealtml)
!         ealt = ealt + ealtml
!!$acc end serial
!!$acc end data
!        use_mlpot= .FALSE.
!      endif
!
!     use Newton's second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the BAOAB recursion
!
      call integrate_vel(derivs,aalt,dta_2)

      if (use_rattle)  call rattle2 (dta)
!
!     increment average virial from fast-evolving potential terms
!
      if (calc_e.or.use_virial) then
!$acc serial async default(present) present(ealt,ealt2)
         ealt = ealt + ealt2
         do i = 1,3; do j = 1,3
               viralt(j,i) = viralt(j,i)+ ( viralt2(j,i)+vir(j,i) )/dinter
            end do; end do
!$acc end serial
      end if
   end do
end
!
!     subroutine baoabrespafast1 :
!     find fast-evolving velocities and positions via BAOAB recursion
!
subroutine baoabrespafast1(ealt,viralt,dta)
   use atmtyp
   use atomsMirror
   use bath
   use cutoff
   use domdec
   use deriv
   use energi
   use freeze
   use inform
   use langevin
   use mdstuf1  ,only: derivs
   use moldyn
   use random_mod
   use timestat
   use units
   use utils
   use utilgpu
   use usage
   use virial
   use mpi
   use utilbaoab
   implicit none
   integer i,j,k,iglob
   real(r_p) dta,dta_2
   real(r_p) a1,a2
   real(r_p) ealt,viralt(3,3)
!
!     set time values and coefficients for BAOAB integration
!
   dta_2 = 0.5_re_p * dta
!
!     initialize virial from fast-evolving potential energy terms
!
   if (use_virial) call zero_virial(viralt)

   do stepfast = 1, nalt2
      call integrate_vel(aalt2,dta_2)
!
      if (use_rattle) then
         call rattle2 (dta_2)
         call save_atoms_pos
      end if
!
      call integrate_pos(dta_2)
!
      if (use_rattle) then
         call rattle(dta_2)
         call rattle2(dta_2) !TODO Ask L about this call
      end if

      call apply_langevin_thermostat(dta)
!
      if (use_rattle) then
         call rattle2(dta)
         call save_atoms_pos
      end if
!
      call integrate_pos(dta_2)
!
!        Reassign the particules that have changed of domain
!        -> real space
      call reassignrespa(stepfast,nalt2)
!
!        communicate positions
!
      call commposrespa(stepfast.ne.nalt2)
!
      if (.not.ftot_l) then
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m (derivs,3*nbloc,rec_queue)
      end if
!
      if (stepfast.eq.nalt2) call mechanicsteprespa1(-1,0)
      if (stepfast.eq.nalt2) call allocsteprespa(.true.)
!
!     get the fast-evolving potential energy and atomic forces
!
      call gradfast(ealt,derivs)
!
!       communicate forces
!
      call comm_forces( derivs,cBond )
!
!       MPI : get total energy
!
      call reduceen(ealt)
!
!     aMD/GaMD contributions
!
      call aMD (derivs,ealt)
!
!     Debug print information
!
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  call info_forces(cBond)
      if (deb_Atom)   call info_minmax_pva
      if (abort)      call emergency_save
      if( abort)      call fatal
!
!     use Newton's second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the BAOAB recursion
!
      call integrate_vel(derivs,aalt2,dta_2)
!
      if (use_rattle)  call rattle2 (dta_2)
!
!     increment average virial from fast-evolving potential terms
!
      if (use_virial) then
!$acc parallel loop collapse(2) async default(present)
         do i = 1,3; do j = 1,3
               viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
            end do; end do
      end if
   end do
end
