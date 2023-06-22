c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine baoab  --  baoab Langevin molecular dynamics step  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "baoab" performs a single molecular dynamics time step
c     via the respa-baoab recursion formula
c
c     literature references:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c     D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-
c     Time-Step Molecular Dynamics Algorithm for Macromolecules",
c     Journal of Physical Chemistry, 98, 6885-6892 (1994)
c
c     X. Qian and T. Schlick, "Efficient Multiple-Time-Step Integrators
c     with Distance-Based Force Splitting for Particle-Mesh-Ewald
c     Molecular Dynamics Simulations", Journal of Chemical Physics,
c     115, 4019-4029 (2001)
c
c
#include "tinker_macro.h"
      subroutine baoabrespa (istep,dt)
      use ani
      use atmtyp
      use atomsMirror
      use bath
      use cutoff
      use domdec
      use deriv,only:info_forces,cBond,cNBond,ftot_l,comm_forces
     &         , zero_forces_rec
      use energi
      use freeze
      use group
      use inform
      use langevin
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
      implicit none
      integer i,j,k,istep,iglob
      real(r_p) dt,dt_x,factor
      real(r_p) dta,dta_2,dt_2
      real(r_p) part1,part2
      real(8) time0,time1
      logical save_mlpot
      logical save_pred 
c
c     set some time values for the dynamics integration
c
      dt_2  = 0.5_re_p * dt
      dta   = dt / dshort
      dta_2 = 0.5_re_p * dta

      if (istep.eq.1) then
         if (use_piston) pres = atmsph
         call set_langevin_thermostat_coeff(dta)
      end if
      if(use_ml_embedding) use_mlpot=.FALSE.
c
c     find quarter step velocities and half step positions via baoab recursion
c
      if (use_piston) call apply_b_piston(dt_2,pres)

      call integrate_vel( a,dt_2 )
c
      if (use_rattle) call rattle2(dt_2)
      if (use_rattle) call save_atoms_pos
c
c     respa inner loop
c
      !initialize virial from fast-evolving potential energy terms
      if (use_virial) call zero_virial(viralt)

      if (use_piston) then
         call apply_a_piston(dt_2,istep,.false.)
         call apply_o_piston(dt)
         if (use_rattle) call rattle (dta_2)
         if (use_rattle) call rattle2(dta_2)

         call reassignrespa(nalt,nalt)
         call commposrespa(.true.)
         if (.not.ftot_l) then
            call prmem_requestm(derivs,3,nbloc,async=.true.)
            call set_to_zero1m(derivs,3*nbloc,rec_queue)
         end if
         call mechanicsteprespa(istep,.true.)
         call allocsteprespa(.true.)
         call gradfast (ealt,derivs)
         call comm_forces( derivs,cBond )
         if (deb_Energy)call info_energy(rank)
         if (deb_Force) call info_forces(cBond)
      end if
c
c     find fast-evolving velocities and positions via BAOAB recursion
c
      do stepfast = 1, nalt
c
        if (use_piston.and.stepfast.eq.1) then
           call integrate_vel( derivs,aalt,dta_2 )
        else
           call integrate_vel( aalt,dta_2 )
        end if
c
        if (use_rattle) call rattle2 (dta_2)
        if (use_rattle) call save_atoms_pos
c
        call integrate_pos( dta_2 )
c
        if (use_rattle) call rattle (dta_2)
        if (use_rattle) call rattle2(dta_2)

        call apply_langevin_thermostat
c
        if (use_rattle) call rattle2(dta_2)
        if (use_rattle) call save_atoms_pos
c
        call integrate_pos( dta_2 )
c
        if (use_rattle) call rattle (dta_2)
        if (use_rattle) call rattle2(dta_2)
c
c       Reassign the particules that have changed of domain
c
c       -> real space
        call reassignrespa(stepfast,nalt)
c
c       communicate positions
c
        call commposrespa(stepfast.ne.nalt)
c
        if (.not.ftot_l) then
           call prmem_requestm(derivs,3,nbloc,async=.true.)
           call set_to_zero1m(derivs,3*nbloc,rec_queue)
        end if
c
        call mechanicsteprespa(istep,.true.)
        call allocsteprespa(.true.)
c
c     get the fast-evolving potential energy and atomic forces
c
        call gradfast ( ealt,derivs )
c
c       communicate forces
c
        call comm_forces( derivs,cBond )
c
c       MPI : get total energy
c
        call reduceen(ealt)
c
c     Debug print information
c
        if (deb_Energy)call info_energy(rank)
        if (deb_Force) call info_forces(cBond)
        if (deb_Atom)  call info_minmax_pva
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
        call integrate_vel( derivs,aalt,dta_2 )
c
        if (use_rattle)  call rattle2 (dta_2)
c
c     increment average virial from fast-evolving potential terms
c
        if (use_virial) then
!$acc parallel loop collapse(2) async present(vir,viralt)
           do i = 1,3; do j = 1,3
              viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
           end do; end do
        end if
      end do

      if (use_piston) then
         if (use_rattle)  call save_atoms_pos
         call apply_a_piston(dt_2,istep,.false.)
         if (use_rattle) call rattle (dta_2)
         if (use_rattle) call rattle2(dta_2)
      end if

c
c     Reassign the particules that have changed of domain
c
c     -> real space
      call reassignrespa(nalt,nalt)
c
c     -> reciprocal space
      call reassignpme(.false.)
c
c     communicate positions
c
      call commposrespa(.false.)
      call commposrec
c
      call reinitnl(istep)
c
      call mechanicsteprespa(istep,.false.)

      call allocsteprespa(.false.)
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)
c
      if (.not.ftot_l) then
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m(derivs,3*nbloc,rec_queue)
      end if
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslow (epot,derivs)
c
c     if necessary, communicate some forces
c
      call comm_forces( derivs,cNBond )
c
c     MPI : get total energy
c
      call reduceen(epot)

      if(use_ml_embedding) then
!$acc parallel loop collapse(2) async
        do i = 1, 3; do j = 1, 3
          viralt(j,i) = vir(j,i) + viralt(j,i)
        end do; end do
c     COMPUTE ML DELTA CONTRIBUTION (ml_embedding_mode)
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
        call commforces(derivs)
!$acc serial async
         epot = epot+eml
!$acc end serial
      endif
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
      call integrate_vel( derivs,a,dt_2 )
c
c     Debug print information
c
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  call info_forces(cNBond)
      if (deb_Atom)   call info_minmax_pva
      if (abort)   __TINKER_FATAL__
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle) call rattle2 (dt_2)
c
c     total potential and virial from sum of fast and slow parts
c
      if (calc_e.or.use_virial) then
!$acc serial async present(epot,ealt,vir,viralt)
         epot = epot + ealt
         do i = 1,3; do j = 1,3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do; end do
!$acc end serial
      end if
c
      call temper   (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
      call pressure2 (epot,temp)
c
c     total energy is sum of kinetic and potential energies
c
      if (calc_e) then
!$acc serial async present(etot,epot,eksum)
         etot = eksum + epot
!$acc end serial
      end if
c
c     compute statistics and save trajectory for this step
c
      call mdsave (istep,dt,epot)
      call mdrestgpu (istep)
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)

      if (use_piston) then
!$acc update host(pres) async
!$acc wait
         call apply_b_piston(dt_2,pres)
      end if

      end
