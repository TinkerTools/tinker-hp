!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine verlet  --  Verlet molecular dynamics step  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "verlet" performs a single molecular dynamics time step
!     via the velocity Verlet multistep recursion formula
!
!
#include "tinker_macro.h"
subroutine verlet (istep,dt)
   use mlpot
   use atmtyp
   use atomsMirror
   use cutoff
   use domdec
   use deriv  ,only: info_forces,prtEForces,cDef,zero_forces_rec&
      &,ftot_l,comm_forces
   use energi ,only: info_energy,calc_e,chk_energy_fluct
   use freeze
   use inform
   use mdstuf1
   use moldyn
   use potent
   use timestat
   use units
   use usage
   use utils  ,only:set_to_zero1m
   use utilgpu,only:prmem_requestm,rec_queue
   use mpi
   use sizes
   use spectra
   implicit none
   integer i,j,istep
   integer iglob
   real(r_p),save :: dip(3),dipind(3)
   real(r_p) dt,dt_2
   real(r_p) time0,time1
!
!     set some time values for the dynamics integration
!
   dt_2 = 0.5_re_p * dt
!
!     store the current atom positions, then find half-step
!     velocities and full-step positions via Verlet recursion
!
   if (use_rattle) call save_atoms_pos
   call integrate_vel( a,dt_2 )
   call integrate_pos( dt )
!
!     get constraint-corrected positions and half-step velocities
!
   if (use_rattle)  call rattle (dt)
!
!     Reassign the particules that have changed of domain
!
!     -> real space
   call reassign
!     -> reciprocal space
   call reassignpme(.false.)
!
!     communicate positions
!
   call commpos
   call commposrec
!
   if (.not.ftot_l) then
      call prmem_requestm(derivs,3,nbloc,async=.true.)
      call set_to_zero1m (derivs,3*nbloc,rec_queue)
   end if
!
   call reinitnl(istep)
!
   call mechanicstep(istep)

   call allocstep
!
!     rebuild the neighbor list
!
   if (use_list) call nblist(istep)
!
!     get the potential energy and atomic forces
!
   call gradient (epot,derivs)
!
!     communicate forces
!
   call comm_forces( derivs )
!
!     MPI : get total energy
!
   call reduceen(epot)
!
!     compute hybrid machine learning potential embedding
!
   if(use_embd_potoff) then
      ml_embedding_mode=3
      call zero_forces_rec
      call gradembedding (eml,derivs)
      call reduceen(eml)
      call comm_forces( derivs )
!$acc serial async
      epot = epot+eml
!$acc end serial
      ml_embedding_mode=1
   endif
!
!     make half-step temperature and pressure corrections
!
   call temper2 (temp)
   call pressure2 (epot,temp)
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
!
!     use Newton's second law to get the next accelerations;
!     find the full-step velocities using the Verlet recursion
!
   call integrate_vel(derivs,a,dt_2)
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle)  call rattle2 (dt)
!
!     make full-step temperature and pressure corrections
!
   call temper   (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
   if(ir) then
      call compute_dipole(dip,dipind,full_dipole)
!$acc update host(dip,dipind) async
!$acc wait
      call save_dipole_traj(dip,dipind)
   endif
!
!     total energy is sum of kinetic and potential energies
!
   if (calc_e) then
!$acc serial present(etot,eksum,epot) async
      etot = eksum + epot
!$acc end serial
      if(rank.eq.0) call chk_energy_fluct(epot,eksum,abort)
   end if
   if (nproc.eq.1.and.tinkerdebug.eq.64) call prtEForces(derivs,epot)
!
!     compute statistics and save trajectory for this step
!
   if (abort)      call emergency_save
   if (abort)      __TINKER_FATAL__
   call mdsave (istep,dt,epot)
   call mdrestgpu (istep)
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)
end
!
!
!     ############################################################################
!     ##                                                                        ##
!     ##  subroutine gradembedding  --  embedding energy & gradient components  ##
!     ##                                                                        ##
!     ############################################################################
!
!
!     "gradembedding" calculates the potential energy and first derivatives
!     for the chosen embedding potential energy terms
!
!
subroutine gradembedding (energy,derivs)
   use potent
   implicit none
   real(r_p) energy
   real(r_p) derivs(3,*)
   logical save_embd_bond,save_embd_angle
   logical save_embd_strbnd,save_embd_urey
   logical save_embd_angang,save_embd_opbend
   logical save_embd_opdist,save_embd_improp
   logical save_embd_imptor,save_embd_tors
   logical save_embd_pitors,save_embd_strtor
   logical save_embd_tortor
!
!
!     save the original state of intra potential energy terms
!
   save_embd_bond   = use_bond
   save_embd_angle  = use_angle
   save_embd_strbnd = use_strbnd
   save_embd_urey   = use_urey
   save_embd_angang = use_angang
   save_embd_opbend = use_opbend
   save_embd_opdist = use_opdist
   save_embd_improp = use_improp
   save_embd_imptor = use_imptor
   save_embd_tors   = use_tors
   save_embd_pitors = use_pitors
   save_embd_strtor = use_strtor
   save_embd_tortor = use_tortor

!
!     turn off chosen intra potential energy terms
!
   if (.not. use_embd_bond)  use_bond = .false.
   if (.not. use_embd_angle)  use_angle = .false.
   if (.not. use_embd_strbnd)  use_strbnd = .false.
   if (.not. use_embd_urey)  use_urey = .false.
   if (.not. use_embd_angang)  use_angang = .false.
   if (.not. use_embd_opbend)  use_opbend = .false.
   if (.not. use_embd_opdist)  use_opdist = .false.
   if (.not. use_embd_improp)  use_improp = .false.
   if (.not. use_embd_imptor)  use_imptor = .false.
   if (.not. use_embd_tors)  use_tors = .false.
   if (.not. use_embd_pitors)  use_pitors = .false.
   if (.not. use_embd_strtor)  use_strtor = .false.
   if (.not. use_embd_tortor)  use_tortor = .false.
!
!     get energy and gradient for potential terms
!
   call gradient (energy,derivs)
!
!     restore the original state of fintra potential energy term
!
   use_bond   = save_embd_bond
   use_angle  = save_embd_angle
   use_strbnd = save_embd_strbnd
   use_urey   = save_embd_urey
   use_angang = save_embd_angang
   use_opbend = save_embd_opbend
   use_opdist = save_embd_opdist
   use_improp = save_embd_improp
   use_imptor = save_embd_imptor
   use_tors   = save_embd_tors
   use_pitors = save_embd_pitors
   use_strtor = save_embd_strtor
   use_tortor = save_embd_tortor

end
