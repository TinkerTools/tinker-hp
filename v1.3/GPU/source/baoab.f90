!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  subroutine baoab  --  BAOAB Langevin molecular dynamics step    ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     "baoab" performs a single molecular dynamics time step
!     via the BAOAB recursion formula
!
!     literature reference:
!
!     Efficient molecular dynamics using geodesic integration
!     and solvent-solute splitting B. Leimkuhler and C. Matthews,
!     Proceedings of the Royal Society A, 472: 20160138, 2016
!
!
#include "tinker_macro.h"
subroutine baoab (istep,dt)
   use atmtyp
   use atomsMirror
   use bath
   use cutoff
   use domdec
   use deriv    ,only:info_forces,cDef,ftot_l,comm_forces
   use energi   ,only: info_energy,calc_e,chk_energy_fluct
   use freeze
   use inform
   use langevin
   use mdstuf
   use mdstuf1
   use moldyn
   use mpi
   use random_mod
   use spectra
   use timestat
   use tinMemory,only: prmem_requestm
   use tinheader,only: re_p
   use units
   use usage
   use utilbaoab
   use utils    ,only: set_to_zero1m
   use utilgpu  ,only: rec_queue,def_queue
   use virial
   implicit none
   integer  ,intent(in):: istep
   real(r_p),intent(in):: dt
   real(r_p),save :: dip(3),dipind(3)
   integer i,j,iglob
   real(r_p) dt_2

   dt_2 = 0.5_re_p*dt

   if (istep.eq.1) then
!$acc enter data create(dip, dipind)
      if (use_piston) pres = atmsph
      call set_langevin_thermostat_coeff(dt)
   end if

   !if (use_piston) call apply_b_piston(dt_2,pres,stress)
!
!     find quarter step velocities and half step positions via BAOAB recursion
!
   call integrate_vel( a,dt_2 )
!
   if (use_rattle) then
      if (rank == 0) then
         write(0,*) "Error: RATTLE not compatible with ",&
            &"BAOAB integrator yet."
         call fatal
      endif
      call rattle2(dt)
      call save_atoms_pos
   end if
!
   if (use_piston) then
      call apply_a_piston(dt_2,-1,.true.)
      call apply_o_piston(dt)
   else
      call integrate_pos( dt_2 )
   end if
!
   if (use_rattle) call rattle (dt_2)
   if (use_rattle) call rattle2(dt_2)
!
   call apply_langevin_thermostat(dt)
!
   if (use_rattle) call rattle2(dt)
   if (use_rattle) call save_atoms_pos
!
!     find full step positions via BAOAB recursion
!
   if(use_piston) then
      call apply_a_piston(dt_2,-1,.TRUE.)
   else
      call integrate_pos( dt_2 )
   end if
!
   if (use_rattle) call rattle (dt_2)
   if (use_rattle) call rattle2(dt_2)
!
!     Reassign the particules that have changed of domain
!
!     -> real space
   call reassign
!
!     -> reciprocal space
   call reassignpme(.false.)
!
!     communicate positions
!
   call commpos
   call commposrec
   call reCast_position
!
   if (.not.ftot_l) then
      call prmem_requestm(derivs,3,nbloc,async=.true.)
      call set_to_zero1m(derivs,3*nbloc,rec_queue)
   end if
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
!     communicate forces
!
   call comm_forces( derivs )
!
!     aMD/GaMD contributions
!
   call aMD (derivs,epot)
!
!     Debug print information
!
   if (deb_Energy) call info_energy(rank)
   if (deb_Force)  call info_forces(cDef)
   if (deb_Atom)   call info_minmax_pva
   if (abort)      call emergency_save
   if (abort)      __TINKER_FATAL__

   call temper   (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
   if(ir) then
      call compute_dipole(dip,dipind,full_dipole)
!$acc update host(dip,dipind) async
!$acc wait
      call save_dipole_traj(dip,dipind)
   endif
!
!     use Newton's second law to get the next accelerations;
!     find the full-step velocities using the BAOAB recursion
!
   call integrate_vel(derivs,a,dt_2)
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle)  call rattle2 (dt)
!
!
!     make half-step temperature and pressure corrections
!
   call pressure2 (epot,temp)
!
!     total energy is sum of kinetic and potential energies
!
   if (calc_e) then
!$acc serial present(epot,eksum,etot) async
      etot = eksum + epot
!$acc end serial
      if(rank.eq.0) call chk_energy_fluct(epot,eksum,abort)
   end if
!
!     compute statistics and save trajectory for this step
!
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)
   call mdsave (istep,dt,epot)
   call mdrestgpu (istep)

   if (use_piston) then
!$acc update host(pres,stress) async
!$acc wait
      call apply_b_piston(dt,pres,stress)
      call ddpme3dnpt(1.0_re_p,istep)
   endif

end
