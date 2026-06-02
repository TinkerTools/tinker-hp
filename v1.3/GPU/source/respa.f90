!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine respa  --  r-RESPA molecular dynamics step  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "respa" performs a single multiple time step molecular dynamics
!     step using the reversible reference system propagation algorithm
!     (r-RESPA) via a Verlet core with the potential split into fast-
!     and slow-evolving portions
!
!     literature references:
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
subroutine respa(istep,dt)
   use mlpot
   use atmtyp
   use atomsMirror
   use bath     ,only: barostat
   use cutoff
   use domdec
   use deriv    ,only: info_forces,cNBond,cBond,ftot_l,comm_forces&
      &,zero_forces_rec
   use energi   ,only: info_energy,calc_e,chk_energy_fluct
   use freeze
   use group
   use inform
   use mdstuf1
   use moldyn
   use timestat
   use tinMemory,only: prmem_requestm
   use tors
   use units
   use usage
   use potent
   use uprior ,only: use_pred
   use utils
   use utilgpu  ,only: rec_queue,openacc_abort
   use virial
   use improp
   use mpi
   use spectra
   implicit none
   integer i,j,k,iglob
   integer istep
   logical save_pred
   real(r_p) dt,dt_2,dt_in
   real(r_p) dta,dta_2
   real(r_p), allocatable, save :: dip(:,:),dipind(:,:)
   logical  ,save:: f_in=.true.
!
!     set some time values for the dynamics integration
!
   dt_2  = 0.5_re_p * dt
   dta   = dt / dshort
   dta_2 = 0.5_re_p * dta

   if (f_in) then
      if (ir) then
         allocate(dip(3,nalt),dipind(3,nalt))
!$acc enter data create(dip,dipind)
      end if
      f_in = .false.
   end if

   if (use_ml_embedding) use_mlpot=.FALSE.
!
!     make half-step temperature and pressure corrections
!
!     call temper (dt)
!
!     store the current atom positions, then find half-step
!     velocities via velocity Verlet recursion
!
   call integrate_vel(a,dt_2,aalt,dta_2)
!
!     initialize virial from fast-evolving potential energy terms
!
   if (use_virial) call zero_virial(viralt)
!
!     find fast-evolving velocities and positions via Verlet recursion
!
   do stepfast = 1, nalt
      if (use_rattle) call save_atoms_pos
      call integrate_pos(dta)

      if (use_rattle) then
         ! call openacc_abort("Rattle not tested with openacc")
         call rattle (dta)
      endif
!
!     Reassign the particules that have changed of domain
!
!        -> real space
      call reassignrespa(stepfast,nalt)
!
!        communicate positions
!
      call commposrespa(stepfast.ne.nalt)
!
      if ( .not.ftot_l ) then
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m( derivs,3*nbloc,rec_queue)
      end if
!
      call mechanicsteprespa(istep,.true.)
      call allocsteprespa(.true.)
!
!     get the fast-evolving potential energy and atomic forces
!
      call gradfast (ealt,derivs)
!
!        communicate forces
!
      call comm_forces( derivs,cBond )
!
!        aMD/GaMD contributions
!
      call aMD (derivs,ealt)
!
!        MPI : get total energy
!
      call reduceen(ealt)
!
!        Debug information
!
      if(deb_Energy) call info_energy(rank)
      if(deb_Force)  call info_forces(cBond)
      if(abort)      call emergency_save
      if(abort)      call fatal

      dt_in = merge(dta_2,dta,use_rattle)  ! level 1
      dt_in = merge(dt_in,dta_2,stepfast.ne.nalt) ! level 0
      call integrate_vel(derivs,aalt,dt_in)
!
!     use Newtons second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the Verlet recursion
!
      if (use_rattle) then
         ! call openacc_abort("Rattle2 not tested with openacc")
         call integrate_vel(aalt,dta_2)
         call rattle2 (dta)
      end if
      if(deb_Atom)   call info_minmax_pva(1)
!
!     increment average virial from fast-evolving potential terms
!
      if (use_virial) then
!$acc parallel loop collapse(2) async default(present)
         do i = 1, 3; do j = 1, 3
               viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
            end do; end do
      end if

      if(ir .and. stepfast<nalt) then
!           call rotpolegpu
         call compute_dipole(dip(:,stepfast)&
            &,dipind(:,stepfast),.FALSE.)
      endif
   end do
!
!     Reassign the particules that have changed of domain
!
!     -> real space
!
   !call reassignrespa(nalt,nalt)
!
!     -> reciprocal space
!
   call reassignpme(.false.)
!
!     communicate positions
!
   call commposrespa(.false.)
   call commposrec
!
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
      call set_to_zero1m( derivs,3*nbloc,rec_queue)
   end if
!
!     get the slow-evolving potential energy and atomic forces
!
   call gradslow (epot,derivs)
!
!     if necessary, communicate some forces
!
   call comm_forces( derivs )
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
      use_mlpot =  .TRUE.
      call zero_forces_rec
      save_pred = use_pred
      use_pred = .FALSE.
      if(use_embd_potoff) then
         call gradembedding2 (eml,derivs)
      else
         call gradient (eml,derivs)
      endif
      use_pred  = save_pred
      call reduceen(eml)
      call comm_forces( derivs )
!$acc serial async present(epot,eml)
      epot = epot+eml
!$acc end serial
   endif
!
!     use Newton's second law to get the slow accelerations;
!     find full-step velocities using velocity Verlet recursion
!
   call integrate_vel(derivs,a,dt_2)
!
!     Debug print information
!
   if(deb_Energy) call info_energy(rank)
   if(deb_Force)  call info_forces(cNBond)
   if(deb_Atom)   call info_minmax_pva
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle) then
      !  call openacc_abort("Rattle2 not tested with openacc")
      call rattle2 (dt)
   end if
!
!     total potential and virial from sum of fast and slow parts
!
   if (calc_e) then
!$acc serial async present(epot,ealt)
      epot = epot + ealt
!$acc end serial
      if(rank.eq.0) call chk_energy_fluct(epot,ealt,abort)
   end if
!
   if (use_virial) then
!$acc parallel loop collapse(2) async default(present)
      do i = 1, 3; do j = 1, 3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do; end do
   end if
!
!     make full-step temperature and pressure corrections
!
   call temper (dt,eksum,ekin,temp)
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
!$acc serial async present(etot,eksum,epot)
      etot = eksum + epot
!$acc end serial
   end if

   ! Fatal Instructions
   if(abort)      call emergency_save
   if(abort)      call fatal
!
!     compute statistics and save trajectory for this step
!
   call mdsave     (istep,dt,epot)
   call mdrestgpu  (istep)
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)
end
!
!
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine gradfast  --  fast energy & gradient components  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "gradfast" calculates the potential energy and first derivatives
!     for the fast-evolving local valence potential energy terms
!
!
subroutine gradfast (energy,derivs)
   use cutoff
   use potent
   use colvars
   use plumed
   implicit none
   real(r_p) energy
   real(r_p) derivs(3,*)
   logical save_vdw,save_charge
   logical save_dipole
   logical save_mpole,save_polar
   logical save_rxnfld,save_solv
   logical save_repuls,save_disp,save_chgtrn
   logical save_list
   logical save_smdvel, save_smdfor
   logical save_colvars
   logical save_plumed
!
!     save the original state of slow-evolving potentials
!
   save_vdw    = use_vdw
   save_charge = use_charge
   save_mpole  = use_mpole
   save_polar  = use_polar
   save_repuls = use_repuls
   save_disp   = use_disp
   save_chgtrn = use_chgtrn
   save_solv   = use_solv
   save_list   = use_list
   save_smdvel = use_smd_velconst
   save_smdfor = use_smd_forconst
   save_colvars = use_colvars
   save_plumed = lplumed
!
!     turn off slow-evolving nonbonded potential energy terms
!
   use_vdw    = .false.
   use_charge = .false.
   use_mpole  = .false.
   use_polar  = .false.
   use_solv   = .false.
   use_repuls = .false.
   use_disp   = .false.
   use_chgtrn = .false.
   use_list   = .false.
   use_smd_velconst = .false.
   use_smd_forconst = .false.
   nonbonded_l      = .false.
   use_colvars = .false.
   lplumed = .false.
!
!     get energy and gradient for fast-evolving potential terms
!
   call gradient (energy,derivs)
!
!     restore the original state of slow-evolving potentials
!
   use_vdw    = save_vdw
   use_charge = save_charge
   use_mpole  = save_mpole
   use_polar  = save_polar
   use_repuls = save_repuls
   use_disp   = save_disp
   use_chgtrn = save_chgtrn
   use_solv   = save_solv
   use_list   = save_list
   use_smd_velconst = save_smdvel
   use_smd_forconst = save_smdfor
   use_colvars = save_colvars
   lplumed = save_plumed
   nonbonded_l      = .true.
end
!
!
!     ##################################################################
!     ##                                                              ##
!     ##  subroutine gradslow  --  slow energy & gradient components  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "gradslow" calculates the potential energy and first derivatives
!     for the slow-evolving nonbonded potential energy terms
!
!
subroutine gradslow (energy,derivs)
   use deriv  ,only: remove_desave
   use domdec
   use energi ,only: esum,esave,calc_e
   use mdstuf
   use moldyn ,only: dinter
   use potent
   use virial
   implicit none
   integer i,j
   real(r_p) energy
   real(r_p) derivs(3,nbloc)
   logical save_bond,save_angle
   logical save_strbnd,save_urey
   logical save_angang,save_opbend
   logical save_opdist,save_improp
   logical save_imptor,save_tors
   logical save_pitors,save_angtor,save_strtor
   logical save_tortor,save_geom
   logical save_metal,save_extra
   logical respa1_l
   logical save_mlpot
!
   respa1_l = (index(integrate,'RESPA1').gt.0)
!
!     save the original state of fast-evolving potentials
!
   save_bond   = use_bond
   save_angle  = use_angle
   save_strbnd = use_strbnd
   save_urey   = use_urey
   save_angang = use_angang
   save_opbend = use_opbend
   save_opdist = use_opdist
   save_improp = use_improp
   save_imptor = use_imptor
   save_tors   = use_tors
   save_pitors = use_pitors
   save_strtor = use_strtor
   save_angtor = use_angtor
   save_tortor = use_tortor
   save_geom   = use_geom
   save_extra  = use_extra
   save_mlpot  = use_mlpot

!
!     turn off fast-evolving valence potential energy terms
!
   use_bond   = .false.
   use_angle  = .false.
   use_strbnd = .false.
   use_urey   = .false.
   use_angang = .false.
   use_opbend = .false.
   use_opdist = .false.
   use_improp = .false.
   use_imptor = .false.
   use_tors   = .false.
   use_pitors = .false.
   use_strtor = .false.
   use_angtor = .false.
   use_tortor = .false.
   use_geom   = .false.
   use_mlpot  = .false.
   use_extra  = .false.
   bonded_l   = .false.
   if (respa1_l) then
      use_vdwlong   = .true.
      use_clong     = .true.
      use_repulslong= .true.
      use_displong  = .true.
      use_mpolelong = .true.
      use_chgtrnlong= .true.
   end if
!
!     get energy and gradient for slow-evolving potential terms
!
   call gradient (energy,derivs)
!
!     restore the original state of fast-evolving potentials
!
   use_bond   = save_bond
   use_angle  = save_angle
   use_strbnd = save_strbnd
   use_urey   = save_urey
   use_angang = save_angang
   use_opbend = save_opbend
   use_opdist = save_opdist
   use_improp = save_improp
   use_imptor = save_imptor
   use_tors   = save_tors
   use_pitors = save_pitors
   use_strtor = save_strtor
   use_angtor = save_angtor
   use_tortor = save_tortor
   use_geom   = save_geom
   use_mlpot  = save_mlpot
   use_extra  = save_extra
   bonded_l   = .true.
   if (respa1_l) then
      use_vdwlong   = .false.
      use_clong     = .false.
      use_mpolelong = .false.
      use_repulslong= .false.
      use_displong  = .false.
      use_chgtrnlong= .false.
!
!        substract the previously stored short range energy and forces
!
      call remove_desave(derivs)

      if (calc_e.or.use_virial) then
!$acc serial present(esum,energy,esave,vir,virsave) async
         esum   = esum - esave
         energy = energy - esave
         do i = 1,3; do j = 1,3
               vir    (j,i) = vir(j,i) - virsave(j,i)/dinter
               virsave(j,i) = 0.0
            end do; end do
!$acc end serial
      end if
   end if

end
!
!
!     #############################################################################
!     ##                                                                         ##
!     ##  subroutine gradembedding2  --  embedding energy & gradient components  ##
!     ##                                                                         ##
!     #############################################################################
!
!
!     "gradembedding2" calculates the potential energy and first derivatives
!     for the chosen embedding potential energy terms to substract for respa
!
!
subroutine gradembedding2 (energy,derivs)
   use cutoff
   use potent
   implicit none
   real(r_p) energy
   real(r_p) derivs(3,*)
   logical save_embd_bond,save_embd_angle
   logical save_embd_strbnd,save_embd_urey
   logical save_embd_angang,save_embd_opbend
   logical save_embd_opdist,save_embd_improp
   logical save_embd_imptor,save_embd_tors
   logical save_embd_pitors,save_embd_angtor,save_embd_strtor
   logical save_embd_tortor,save_embd_geom
   logical save_embd_metal,save_embd_extra
   logical save_embd_mlpot
   logical save_embd_vdw,save_embd_charge
   logical save_embd_dipole
   logical save_embd_mpole,save_embd_polar
   logical save_embd_rxnfld,save_embd_solv
   logical save_embd_list
   logical save_embd_smdvel, save_embd_smdfor
!
!
!     save the original state of potential energy terms
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
   save_embd_angtor = use_angtor
   save_embd_tortor = use_tortor
   save_embd_geom   = use_geom
   save_embd_extra  = use_extra
   save_embd_mlpot  = use_mlpot
   save_embd_vdw    = use_vdw
   save_embd_charge = use_charge
   save_embd_mpole  = use_mpole
   save_embd_polar  = use_polar
   save_embd_solv   = use_solv
   save_embd_list   = use_list
   save_embd_smdvel = use_smd_velconst
   save_embd_smdfor = use_smd_forconst
!
!     turn on only chosen intra potential energy terms
!     that will be substract for respa
!
   call potoff

   use_geom         = save_embd_geom
   use_extra        = save_embd_extra
   use_solv         = save_embd_solv
   use_smd_velconst = save_embd_smdvel
   use_smd_forconst = save_embd_smdfor

   if (.not. use_embd_bond)  use_bond = .true.
   if (.not. use_embd_angle)  use_angle = .true.
   if (.not. use_embd_strbnd)  use_strbnd = .true.
   if (.not. use_embd_urey)  use_urey  = .true.
   if (.not. use_embd_angang)  use_angang = .true.
   if (.not. use_embd_opbend)  use_opbend = .true.
   if (.not. use_embd_opdist)  use_opdist = .true.
   if (.not. use_embd_improp)  use_improp = .true.
   if (.not. use_embd_imptor)  use_imptor = .true.
   if (.not. use_embd_tors)  use_tors = .true.
   if (.not. use_embd_pitors)  use_pitors = .true.
   if (.not. use_embd_strtor)  use_strtor = .true.
   if (.not. use_embd_tortor)  use_tortor = .true.
!
!     get energy and gradient for potential terms
!
   call gradient (energy,derivs)
!
!     restore the original state of potential energy term
!
   use_bond         = save_embd_bond
   use_angle        = save_embd_angle
   use_strbnd       = save_embd_strbnd
   use_urey         = save_embd_urey
   use_angang       = save_embd_angang
   use_opbend       = save_embd_opbend
   use_opdist       = save_embd_opdist
   use_improp       = save_embd_improp
   use_imptor       = save_embd_imptor
   use_tors         = save_embd_tors
   use_pitors       = save_embd_pitors
   use_strtor       = save_embd_strtor
   use_angtor       = save_embd_angtor
   use_tortor       = save_embd_tortor
   use_mlpot        = save_embd_mlpot
   use_vdw          = save_embd_vdw
   use_charge       = save_embd_charge
   use_mpole        = save_embd_mpole
   use_polar        = save_embd_polar
   use_list         = save_embd_list

end
