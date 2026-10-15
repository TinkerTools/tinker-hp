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
   use atomsMirror, only: integrate_pos,integrate_vel,save_atoms_pos
   use cutoff, only: use_list
   use domdec, only: nbloc,rank
   use deriv    ,only: info_forces,cNBond,cBond,ftot_l,comm_forces&
      &,save_total_forces,restore_total_forces
   use energi   ,only: info_energy,calc_e,chk_energy_fluct,reduceen
   use freeze, only: use_rattle
   use inform, only: abort,deb_Atom,deb_Energy,deb_Force,info_minmax_pva
   use mdstuf1, only: derivs,ealt,ekin,eksum,eml,epot,etot,pres,stress,temp,viralt
   use moldyn, only: a,aalt,dshort,nalt,stepfast
   use parallel, only: commpos, commposrec, commposrespa, reassign_respa
   use tinMemory,only: prmem_requestm
   use potent
   use utils, only: set_to_zero1m
   use mdevice  ,only: reciprocal_work
   use virial, only: use_virial,vir,virsave,zero_virial
   use spectra, only: compute_dipole,full_dipole,ir,save_dipole_respa
   use tinheader, only: re_p
   implicit none
   integer i,j,k,iglob
   integer istep
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
      call reassign_respa(stepfast,nalt)
!
!        communicate positions
!
      call commposrespa(stepfast.ne.nalt)
!
      if ( .not.ftot_l ) then
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m( derivs,3*nbloc,reciprocal_work%queue)
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
      if (calc_e.or.use_virial) call reduceen(ealt)
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
!$acc parallel loop collapse(2) async(reciprocal_work%queue) default(present)
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
      call set_to_zero1m( derivs,3*nbloc,reciprocal_work%queue)
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
   if (calc_e.or.use_virial) call reduceen(epot)

   call grad_embd_slow(istep,eml,epot,derivs)
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
   if(deb_Force)  call info_forces(cBond)
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
   if (calc_e.or.use_virial) then
!$acc serial async(reciprocal_work%queue) present(epot,ealt)
      epot = epot + ealt
!$acc end serial
      if (calc_e.and.rank.eq.0) call chk_energy_fluct(epot,ealt,abort)
   end if
!
   if (use_virial) then
!$acc parallel loop collapse(2) async(reciprocal_work%queue) default(present)
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
   use cutoff, only: use_list
   use energi , only: calc_e,esave,emlpot
   use mlpot  , only: use_dmts,model_file_cur,model_file2
   use potent
   use colvars, only: use_colvars
   use plumed, only: lplumed
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
   if (use_dmts) then
     model_file_cur = model_file2
   endif
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

   if (use_dmts.and.calc_e) then
!$acc serial present(esave,emlpot) async
      esave = emlpot
!$acc end serial
   end if
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
   use domdec, only: nbloc
   use energi ,only: esum,esave,calc_e
   use mdstuf, only: integrate
   use moldyn ,only: dinter
   use inform ,only: minmaxone
   use mlpot
   use potent
   use virial, only: use_virial,vir,virsave
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
   model_file_cur = model_file
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
   use_mlpot  = use_dmts
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
   if (use_dmts) then
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
!     #########################################################################
!     ##                                                                     ##
!     ##  subroutine grad_embd_slow  -- shared slow embedding correction  ##
!     ##                                                                     ##
!     #########################################################################
!
!
subroutine grad_embd_slow(istep,energy,epot,derivs)
   use deriv  , only: comm_forces,restore_total_forces,save_total_forces
   use energi , only: calc_e,reduceen
   use mdstuf1, only: viralt
   use mdevice, only: reciprocal_work
   use mlpot  , only: ml_embedding_mode,mlpotscale
   use potent , only: use_embd_potoff,use_ml_embedding,use_mlpot
   use uprior , only: use_pred
   use virial , only: use_virial,vir
   implicit none
   integer i,j
   integer, intent(in) :: istep
   real(r_p), intent(inout) :: energy,epot
   real(r_p) derivs(3,*)
   logical save_mlpot,save_pred
   integer save_ml_mode
   real(t_p) save_mlpotscale

   if (.not.use_ml_embedding) return

   if (use_virial) then
!$acc parallel loop collapse(2) async(reciprocal_work%queue) present(vir,viralt)
      do i = 1, 3; do j = 1, 3
         viralt(j,i) = vir(j,i) + viralt(j,i)
      end do; end do
   end if
   call save_total_forces

   save_ml_mode = ml_embedding_mode
   save_mlpotscale = mlpotscale
   save_mlpot = use_mlpot
   use_mlpot = .true.
   save_pred = use_pred
   use_pred = .false.
   if (use_embd_potoff) then
      call gradembedding2(energy,derivs)
   else
      call gradient_helper(energy,derivs,0,.true.)
   endif
   use_pred = save_pred
   use_mlpot = save_mlpot
   ml_embedding_mode = save_ml_mode
   mlpotscale = save_mlpotscale

   if (calc_e.or.use_virial) call reduceen(energy)
   call comm_forces(derivs)
   if (istep.eq.1 .and. .not.use_embd_potoff) then
      call audit_respa_embedding_forces(derivs)
   else
      call restore_total_forces
   endif
   if (calc_e.or.use_virial) then
!$acc serial async(reciprocal_work%queue) present(epot,energy)
      epot = epot + energy
!$acc end serial
   endif
end
!
!
!     #########################################################################
!     ##                                                                     ##
!     ##  subroutine gradient_helper  -- selected embedding gradient        ##
!     ##                                                                     ##
!     #########################################################################
!
!
!     Evaluate the selected potential terms. Mode 0 selects both fast and slow
!     terms, modes 1 and 2 select one part, and mode 3 selects intermediate
!     short-range terms.
!     When eval_ml is true, the ML contribution is evaluated with embedding
!     weights; otherwise this is an AMOEBA-only helper call.
!
subroutine gradient_helper (energy,derivs,mode,eval_ml)
   use energi , only: calc_e,ep,esave
   use colvars, only: use_colvars
   use cutoff , only: use_list
   use mdevice, only: reciprocal_work
   use polpot , only: polalg,polalgshort,tcgorder,tcgordershort,tcgprec,tcgprecshort &
              , tcgguess,tcgguessshort,tcgpeek,tcgpeekshort,tcgomega,tcgomegashort
   use plumed , only: lplumed
   use potent
   implicit none
   real(r_p) energy
   real(r_p) derivs(3,*)
   integer, intent(in) :: mode
   logical, intent(in) :: eval_ml
   integer, parameter  :: mode_default=0,mode_fast=1,mode_slow=2,mode_inter=3
   logical fast_l,slow_l,inter_l
   logical save_bond,save_angle
   logical save_strbnd,save_urey
   logical save_angang,save_opbend
   logical save_opdist,save_improp
   logical save_imptor,save_tors
   logical save_pitors,save_angtor,save_strtor
   logical save_tortor,save_geom,save_extra
   logical save_mlpot,save_bonded_l,save_nonbonded_l
   logical save_vdw,save_charge,save_mpole,save_polar
   logical save_repuls,save_disp,save_chgtrn,save_solv,save_list
   logical save_crec,save_mrec,save_prec,save_disprec
   logical save_cself,save_mself,save_pself,save_dispself
   logical save_cshortreal,save_mpshort,save_vdwshort
   logical save_pshort,save_repulsshort,save_dispshort
   logical save_dsreal,save_chgtrnshort
   integer save_polalg,save_tcgorder
   logical save_tcgprec,save_tcgguess,save_tcgpeek
   logical save_shortnb
   logical save_smd_vel,save_smd_for
   logical save_colvars,save_plumed
   real(t_p) save_tcgomega

   save_bond        = use_bond
   save_angle       = use_angle
   save_strbnd      = use_strbnd
   save_urey        = use_urey
   save_angang      = use_angang
   save_opbend      = use_opbend
   save_opdist      = use_opdist
   save_improp      = use_improp
   save_imptor      = use_imptor
   save_tors        = use_tors
   save_pitors      = use_pitors
   save_strtor      = use_strtor
   save_angtor      = use_angtor
   save_tortor      = use_tortor
   save_geom        = use_geom
   save_extra       = use_extra
   save_mlpot       = use_mlpot
   save_bonded_l    = bonded_l
   save_nonbonded_l = nonbonded_l
   save_vdw         = use_vdw
   save_charge      = use_charge
   save_mpole       = use_mpole
   save_polar       = use_polar
   save_repuls      = use_repuls
   save_disp        = use_disp
   save_chgtrn      = use_chgtrn
   save_solv        = use_solv
   save_list        = use_list
   save_smd_vel     = use_smd_velconst
   save_smd_for     = use_smd_forconst
   save_colvars     = use_colvars
   save_plumed      = lplumed

   select case (mode)
   case (mode_default)
      fast_l  = .true.
      slow_l  = .true.
      inter_l = .false.
   case (mode_fast)
      fast_l  = .true.
      slow_l  = .false.
      inter_l = .false.
   case (mode_slow)
      fast_l  = .false.
      slow_l  = .true.
      inter_l = .false.
   case (mode_inter)
      fast_l  = .false.
      slow_l  = .false.
      inter_l = .true.
   case default
      call fatal
   end select

   if (inter_l) then
      save_crec        = use_crec
      save_mrec        = use_mrec
      save_prec        = use_prec
      save_disprec     = use_disprec
      save_cself       = use_cself
      save_mself       = use_mself
      save_pself       = use_pself
      save_dispself    = use_dispself
      save_cshortreal  = use_cshortreal
      save_mpshort     = use_mpoleshortreal
      save_vdwshort    = use_vdwshort
      save_pshort      = use_polarshortreal
      save_repulsshort = use_repulsshort
      save_dispshort   = use_dispshort
      save_dsreal      = use_dispshortreal
      save_chgtrnshort = use_chgtrnshort
      save_polalg      = polalg
      save_tcgorder    = tcgorder
      save_tcgprec     = tcgprec
      save_tcgguess    = tcgguess
      save_tcgpeek     = tcgpeek
      save_tcgomega    = tcgomega
      save_shortnb     = shortnonbonded_l

      use_bond     = .false.
      use_angle    = .false.
      use_strbnd   = .false.
      use_urey     = .false.
      use_angang   = .false.
      use_opbend   = .false.
      use_opdist   = .false.
      use_improp   = .false.
      use_imptor   = .false.
      use_tors     = .false.
      use_pitors   = .false.
      use_strtor   = .false.
      use_angtor   = .false.
      use_tortor   = .false.
      use_geom     = .false.
      use_extra    = .false.
      use_crec     = .false.
      use_mrec     = .false.
      use_prec     = .false.
      use_disprec  = .false.
      use_cself    = .false.
      use_mself    = .false.
      use_pself    = .false.
      use_dispself = .false.
      use_cshortreal     = .true.
      use_mpoleshortreal = .true.
      use_vdwshort       = .true.
      use_polarshortreal = .true.
      use_repulsshort    = .true.
      use_dispshort      = .true.
      use_dispshortreal  = .true.
      use_chgtrnshort    = .true.
      polalg             = polalgshort
      tcgorder           = tcgordershort
      tcgprec            = tcgprecshort
      tcgguess           = tcgguessshort
      tcgpeek            = tcgpeekshort
      tcgomega           = tcgomegashort
      use_smd_velconst   = .false.
      use_smd_forconst   = .false.
      use_colvars        = .false.
      lplumed            = .false.
      shortnonbonded_l   = .true.
   else
      if (.not.slow_l) then
         use_vdw          = .false.
         use_charge       = .false.
         use_mpole        = .false.
         use_polar        = .false.
         use_repuls       = .false.
         use_disp         = .false.
         use_chgtrn       = .false.
         use_solv         = .false.
         use_list         = .false.
         use_smd_velconst = .false.
         use_smd_forconst = .false.
         lplumed          = .false.
      endif
      if (.not.fast_l) then
         use_bond    = .false.
         use_angle   = .false.
         use_strbnd  = .false.
         use_urey    = .false.
         use_angang  = .false.
         use_opbend  = .false.
         use_opdist  = .false.
         use_improp  = .false.
         use_imptor  = .false.
         use_tors    = .false.
         use_pitors  = .false.
         use_strtor  = .false.
         use_angtor  = .false.
         use_tortor  = .false.
         use_geom    = .false.
         use_extra   = .false.
         use_colvars = .false.
         lplumed     = .false.
      endif
   endif
   use_mlpot   = eval_ml
   if (inter_l) then
      bonded_l    = .false.
      nonbonded_l = save_nonbonded_l
   else
      bonded_l    = fast_l
      nonbonded_l = slow_l
   endif

   call gradient (energy,derivs)

   use_bond    = save_bond
   use_angle   = save_angle
   use_strbnd  = save_strbnd
   use_urey    = save_urey
   use_angang  = save_angang
   use_opbend  = save_opbend
   use_opdist  = save_opdist
   use_improp  = save_improp
   use_imptor  = save_imptor
   use_tors    = save_tors
   use_pitors  = save_pitors
   use_strtor  = save_strtor
   use_angtor  = save_angtor
   use_tortor  = save_tortor
   use_geom    = save_geom
   use_extra   = save_extra
   use_mlpot   = save_mlpot
   bonded_l    = save_bonded_l
   nonbonded_l = save_nonbonded_l
   use_vdw     = save_vdw
   use_charge  = save_charge
   use_mpole   = save_mpole
   use_polar   = save_polar
   use_repuls  = save_repuls
   use_disp    = save_disp
   use_chgtrn  = save_chgtrn
   use_solv    = save_solv
   use_list    = save_list
   use_smd_velconst = save_smd_vel
   use_smd_forconst = save_smd_for
   use_colvars = save_colvars
   lplumed     = save_plumed

   if (inter_l) then
      use_crec           = save_crec
      use_mrec           = save_mrec
      use_prec           = save_prec
      use_disprec        = save_disprec
      use_cself          = save_cself
      use_mself          = save_mself
      use_pself          = save_pself
      use_dispself       = save_dispself
      use_cshortreal     = save_cshortreal
      use_mpoleshortreal = save_mpshort
      use_vdwshort       = save_vdwshort
      use_polarshortreal = save_pshort
      use_repulsshort    = save_repulsshort
      use_dispshort      = save_dispshort
      use_dispshortreal  = save_dsreal
      use_chgtrnshort    = save_chgtrnshort
      polalg             = save_polalg
      tcgorder           = save_tcgorder
      tcgprec            = save_tcgprec
      tcgguess           = save_tcgguess
      tcgpeek            = save_tcgpeek
      tcgomega           = save_tcgomega
      shortnonbonded_l   = save_shortnb
   endif

   if (inter_l .and. calc_e) then
!$acc serial present(esave,ep) async(reciprocal_work%queue)
      esave = ep
!$acc end serial
   endif
end

!
!     #########################################################################
!     ##                                                                     ##
!     ##  subroutine audit_respa_embedding_forces  -- MTS force identity     ##
!     ##                                                                     ##
!     #########################################################################
!
!
!     Check the standard RESPA slow-stage force after restoration without
!     reusing the non-MTS Verlet audit.  The production force is held in
!     desave while the independent slow decomposition is rebuilt.
!
subroutine audit_respa_embedding_forces(derivs)
   use deriv    , only: comm_forces,ftot_l,get_ftot,restore_total_forces,zero_forces
   use domdec   , only: glob,nbloc,nproc,rank
   use group    , only: grplist,wgrp
   use mlpot    , only: ml_embedding_mode,mlpotscale
   use mdevice  , only: reciprocal_work
   use potent   , only: shortnonbonded_l,use_embd_potoff,use_ml_embedding,use_mlpot
   use sizes    , only: tinkerdebug
   use tinheader, only: re_p
   use uprior   , only: use_pred
   implicit none
   integer i,j,n,saved_mode
   real(r_p) derivs(3,*)
   real(r_p) audit_energy
   real(r_p) max_expected,max_abs,max_abs_restored
   real(r_p) ml_w_norm,tolerance
   real(r_p), allocatable :: amoeba_slow(:,:),amoeba_p(:,:),ml_p(:,:)
   real(r_p), allocatable :: delta_ml(:,:),expected(:,:),observed(:,:)
   real(r_p), allocatable :: restored(:,:)
   real(r_p), allocatable :: snapshot(:,:)
   real(t_p), allocatable :: wgrp_saved(:,:)
   real(t_p) saved_scale
   logical saved_use_mlpot,saved_pred,saved_shortnb,is_p

   if (.not.use_ml_embedding .or. .not.ftot_l .or.&
      &rank.ne.0 .or. nproc.ne.1 .or. tinkerdebug.eq.0) then
      call restore_total_forces
      return
   endif

   n = nbloc
   allocate(amoeba_slow(3,n),amoeba_p(3,n),ml_p(3,n),delta_ml(3,n))
   allocate(expected(3,n),observed(3,n),restored(3,n),snapshot(3,n))
   allocate(wgrp_saved(size(wgrp,1),size(wgrp,2)))
   wgrp_saved = wgrp
   saved_mode = ml_embedding_mode
   saved_scale = mlpotscale
   saved_use_mlpot = use_mlpot
   saved_pred = use_pred
   saved_shortnb = shortnonbonded_l
   shortnonbonded_l = .false.
!$acc data create(audit_energy,snapshot,observed,restored)
   ! Independently rebuild the slow AMOEBA reference A_slow(P union W).
   use_mlpot = .false.
   use_pred = .false.
   ml_embedding_mode = 0
   call zero_forces
   call gradient_helper(audit_energy,snapshot,2,.false.)
   call comm_forces(snapshot)
   call get_ftot(snapshot,nbloc)
!$acc update host(snapshot) async(reciprocal_work%queue)
!$acc wait(reciprocal_work%queue)
   amoeba_slow = snapshot

   ! Rebuild the full isolated AMOEBA reference A(P) with embedding weights.
   ml_embedding_mode = 2
   call set_embedding_weights
   call zero_forces
   call gradient_helper(audit_energy,snapshot,0,.false.)
   call comm_forces(snapshot)
   call get_ftot(snapshot,nbloc)
!$acc update host(snapshot) async(reciprocal_work%queue)
!$acc wait(reciprocal_work%queue)
   amoeba_p = -snapshot

   ! Rebuild ML(P) independently; its mode-2 total contains ML(P)-A(P).
   wgrp = wgrp_saved
!$acc update device(wgrp) async(reciprocal_work%queue)
   ml_embedding_mode = 2
   use_mlpot = .true.
   call zero_forces
   call gradient_helper(audit_energy,snapshot,0,.true.)
   call comm_forces(snapshot)
   call get_ftot(snapshot,nbloc)
!$acc update host(snapshot) async(reciprocal_work%queue)
!$acc wait(reciprocal_work%queue)
   ml_p = snapshot + amoeba_p
   delta_ml = ml_p - amoeba_p

   expected = amoeba_slow + delta_ml

   ! Rebuild the actual production correction after the reference calls so
   ! the caller receives the normal polarization and virial state.
   wgrp = wgrp_saved
!$acc update device(wgrp) async(reciprocal_work%queue)
   ml_embedding_mode = saved_mode
   mlpotscale = saved_scale
   use_mlpot = .true.
   use_pred = .false.
   call zero_forces
   if (use_embd_potoff) then
      call gradembedding2(audit_energy,derivs)
   else
      call gradient_helper(audit_energy,derivs,0,.true.)
   endif
   call comm_forces(derivs)
   call get_ftot(observed,nbloc)
!$acc update host(observed) async(reciprocal_work%queue)
!$acc wait(reciprocal_work%queue)
   call restore_total_forces
   call get_ftot(restored,nbloc)
!$acc update host(restored) async(reciprocal_work%queue)
!$acc wait(reciprocal_work%queue)

   max_expected = maxval(abs(expected))
   max_abs = 0.0_re_p
   max_abs_restored = 0.0_re_p
   ml_w_norm = 0.0_re_p
   do i = 1,n
      do j = 1,3
         max_abs = max(max_abs,abs(observed(j,i)-delta_ml(j,i)))
         max_abs_restored = max(max_abs_restored,&
            &abs(restored(j,i)-expected(j,i)))
      enddo
      is_p = .false.
      if (glob(i).le.size(grplist)) is_p = grplist(glob(i)).eq.1
      if (.not.is_p) then
         do j = 1,3
            ml_w_norm = ml_w_norm + ml_p(j,i)*ml_p(j,i)
         enddo
      endif
   enddo
   ml_w_norm = sqrt(ml_w_norm)
   tolerance = 1.0e-5_re_p * (1.0_re_p + max_expected)
   write(*,'(A)') 'ML RESPA force decomposition audit:'
   write(*,'(A,2ES16.8)') '  delta ML max abs/tolerance: ',max_abs,tolerance
   write(*,'(A,ES16.8)') '  restored total max abs: ',max_abs_restored
   write(*,'(A,ES16.8)') '  correction norm on W: ',ml_w_norm

   wgrp = wgrp_saved
!$acc update device(wgrp) async(reciprocal_work%queue)
   ml_embedding_mode = saved_mode
   mlpotscale = saved_scale
   use_mlpot = saved_use_mlpot
   use_pred = saved_pred
   shortnonbonded_l = saved_shortnb
!$acc end data
   deallocate(amoeba_slow,amoeba_p,ml_p,delta_ml,expected,observed,&
      &restored,snapshot)
   deallocate(wgrp_saved)

   if (max_abs.gt.tolerance .or. max_abs_restored.gt.tolerance .or.&
      &ml_w_norm.gt.tolerance) then
      write(0,*) 'ML RESPA force decomposition audit failed'
      __TINKER_FATAL__
   endif
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
   use cutoff, only: use_list
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
