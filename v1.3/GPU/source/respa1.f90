!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine respa1  --  r-RESPA1 molecular dynamics step  ##
!     ##                                                           ##
!     ###############################################################
!
!
!     "respa1" performs a single multiple time step molecular dynamics
!     step using the reversible reference system propagation algorithm
!     (r-RESPA) via a Verlet core with the potential split into fast-
!     intermediate and slow-evolving portions
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
!     Ruhong Zhou, Edward Harder, Huafeng Xu and B. J. Berne,
!     "Efficient multiple time step method for use with Ewald and
!     particle mesh Ewald for large biomolecular systems",
!     J. Chem. Phys. 115, 2348-2358 (2001)
!
!
#include "tinker_macro.h"
subroutine respa1(istep,dt)
   use mlpot
   use atmtyp
   use atomsMirror
   use cutoff
   use domdec
   use deriv
   use energi  ,only: info_energy,calc_e
   use freeze
   use group
   use inform
   use mdstuf1
   use moldyn
   use timestat
   use potent
   use uprior  ,only: use_pred
   use utils   ,only: set_to_zero1m
   use utilgpu ,only: prmem_requestm,rec_queue
   use units
   use usage
   use virial
   use mpi
   implicit none
   integer i,j,iglob
   integer istep
   real(r_p) dt,dt_2
   real(r_p) dta,dta_2,dta2
   real(8) time0,time1
   logical save_pred
!
!     set some time values for the dynamics integration
!
   dt_2  = 0.5_re_p * dt
   dta   = dt / dinter
   dta_2 = 0.5_re_p * dta
   dta2  = dta / dshort
   if(use_ml_embedding) use_mlpot=.FALSE.
!
!     store the current atom positions, then find half-step
!     velocities via velocity Verlet recursion
!
   call integrate_vel( a,dt_2 )
!
!     find intermediate-evolving velocities and positions via velocity Verlet recursion
!
   call respaint1(ealt,viralt,dta,dta2)
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
      use_mlpot= .TRUE.
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
!     use Newton's second law to get the slow accelerations;
!     find full-step velocities using velocity Verlet recursion
!
   call integrate_vel( derivs,a,dt_2 )
!
!     Debug print information
!
   if(deb_Energy) call info_energy(rank)
   if(deb_Force)  call info_forces(cNBond)
   if(deb_Atom)   call info_minmax_pva
   if(abort)      call fatal
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle)  call rattle2 (dt)
!
!     total potential and virial from sum of fast and slow parts
!
   if (calc_e.or.use_virial) then
!$acc serial async default(present) present(epot,ealt)
      epot = epot + ealt
      do i = 1,3; do j = 1,3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do; end do
!$acc end serial
   end if
!
!     make full-step temperature and pressure corrections
!
   call temper   (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
!
!     total energy is sum of kinetic and potential energies
!
   if (calc_e) then
!$acc serial present(etot,eksum,epot) async
      etot = eksum + epot
!$acc end serial
   end if
!
!     compute statistics and save trajectory for this step
!
   call mdsave (istep,dt,epot)
   call mdrestgpu (istep)
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)
end
!
!     ##########################################################################
!     ##                                                                      ##
!     ##  subroutine gradint  --  intermediate energy & gradient components  ##
!     ##                                                                      ##
!     ##########################################################################
!
!
!     for the fast-evolving local valence potential energy terms
!
!
subroutine gradint (energy,derivs)
   use cutoff
   use deriv
   use domdec
   use energi
   use polpot
   use potent
   use colvars
   use plumed
   implicit none
   real(r_p) energy
   real(r_p) derivs(3,*)
   real(t_p) save_tcgomega
   integer i,j
   logical save_bond,save_angle
   logical save_strbnd,save_urey
   logical save_angang,save_opbend
   logical save_opdist,save_improp
   logical save_imptor,save_tors
   logical save_pitors,save_strtor
   logical save_angtor,save_tortor,save_geom
   logical save_extra
   logical save_mrec,save_disprec
   logical save_prec,save_crec
   integer save_polalg,save_tcgorder
   logical save_tcgprec,save_tcgguess,save_tcgpeek
   logical save_smdvel, save_smdfor
   logical save_colvars
   logical save_plumed
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
   save_crec   = use_crec
   save_mrec   = use_mrec
   save_prec   = use_prec
   save_disprec  = use_disprec
   save_polalg = polalg
   save_tcgorder = tcgorder
   save_tcgprec  = tcgprec
   save_tcgguess = tcgguess
   save_tcgpeek  = tcgpeek
   save_tcgomega = tcgomega
   save_smdvel   = use_smd_velconst
   save_smdfor   = use_smd_forconst
   save_colvars = use_colvars
   save_plumed = lplumed
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
   use_extra  = .false.
   use_crec   = .false.
   use_mrec   = .false.
   use_prec   = .false.
   use_disprec= .false.
   use_cself  = .false.
   use_mself  = .false.
   use_pself  = .false.
   use_dispself       = .false.
   use_cshortreal     = .true.
   use_mpoleshortreal = .true.
   use_vdwshort       = .true.
   use_polarshortreal = .true.
   use_repulsshort    = .true.
   use_dispshort      = .true.
   use_dispshortreal  = .true.
   use_chgtrnshort    = .true.
   polalg     = polalgshort
   tcgorder   = tcgordershort
   tcgprec    = tcgprecshort
   tcgguess   = tcgguessshort
   tcgpeek    = tcgpeekshort
   tcgomega   = tcgomegashort
   use_smd_velconst   = .false.
   use_smd_forconst   = .false.
   bonded_l    = .false.
   shortnonbonded_l   = .true.
   use_colvars = .false.
   lplumed = .false.
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
   use_extra  = save_extra
   use_crec   = save_crec
   use_mrec   = save_mrec
   use_prec   = save_prec
   use_disprec= save_disprec
   use_cself  = .true.
   use_mself  = .true.
   use_pself  = .true.
   use_dispself       = .true.
   use_cshortreal     = .false.
   use_mpoleshortreal = .false.
   use_vdwshort       = .false.
   use_polarshortreal = .false.
   use_repulsshort    = .false.
   use_dispshort      = .false.
   use_dispshortreal  = .false.
   use_chgtrnshort    = .false.
   polalg     = save_polalg
   tcgorder   = save_tcgorder
   tcgprec    = save_tcgprec
   tcgguess   = save_tcgguess
   tcgpeek    = save_tcgpeek
   use_smd_velconst   = save_smdvel
   use_smd_forconst   = save_smdfor
   bonded_l   = .true.
   shortnonbonded_l   = .false.
   use_colvars = save_colvars
   lplumed = save_plumed

   if (calc_e) then
!$acc serial present(esave,ep) async
      esave = ep
!$acc end serial
   end if
end
!
!     subroutine respaint1 :
!     find intermediate-evolving velocities and positions via velocity Verlet recursion
!
subroutine respaint1(ealt,viralt,dta,dta2)
   use mlpot
   use atomsMirror ,only: integrate_vel
   use atmtyp
   use cutoff
   use deriv
   use domdec
   use deriv   ,only: zero_forces_rec
   use energi  ,only: info_energy,calc_e
   use freeze
   use inform
   use mdstuf1 ,only: ealt2,viralt2,derivs
   use moldyn
   use mutant
   use potent
   use timestat
   use uprior  ,only: use_pred
   use utils   ,only: set_to_zero1m
   use utilgpu ,only:prmem_requestm,rec_queue
   use units
   use usage
   use virial
   use mpi
   implicit none
   integer i,j,k,iglob
   logical save_pred
   real(r_p) dta,dta_2,dta2
   real(r_p) ealt,viralt(3,3)

   dta_2 = 0.5_re_p * dta
!
!     initialize virial from fast-evolving potential energy terms
!
   if (use_virial) call zero_virial(viralt)

   do stepint = 1, nalt
      call integrate_vel( aalt,dta_2 )
!
!     find fast-evolving velocities and positions via velocity Verlet recursion
!
      call respafast1(ealt2,viralt2,dta2)
!
      call mechanicsteprespa1(-1,1)
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
!      communicate forces
!
      call comm_forces( derivs,cSNBond )
!
!      MPI : get total energy
!
      call reduceen(ealt)
!
!     Debug print information
!
      if(deb_Energy) call info_energy(rank)
      if(deb_Force ) call info_forces(cSNBond)
      if(deb_Atom  ) call info_minmax_pva
      if(abort     ) __TINKER_FATAL__
!
!     use Newton's second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the Verlet recursion
!
      call integrate_vel( derivs,aalt,dta_2 )
      if (use_rattle)  call rattle2 (dta)
!
!     increment average virial from fast-evolving potential terms
!
      if (calc_e.or.use_virial) then
!$acc serial async default(present) present(ealt,ealt2)
         ealt = ealt + ealt2
         do i = 1,3; do j = 1,3
               viralt(j,i) = viralt(j,i) + (viralt2(j,i) + vir(j,i))/dinter
            end do; end do
!$acc end serial
      end if
   end do
end
!
!     subroutine respafast1 :
!     find fast-evolving velocities and positions via velocity Verlet recursion
!
subroutine respafast1(ealt,viralt,dta)
   use atmtyp
   use atomsMirror
   use cutoff
   use domdec
   use deriv
   use energi  ,only: info_energy
   use freeze
   use inform
   use mdstuf1 ,only: derivs
   use moldyn
   use mpi
   use timestat
   use utils   ,only: set_to_zero1m
   use utilgpu ,only: prmem_requestm,rec_queue
   use units
   use usage
   use virial
   implicit none
   integer i,j,k,iglob
   real(r_p) dta,dta_2
   real(r_p) ealt,viralt(3,3)

   dta_2 = 0.5_re_p * dta
!
!     initialize virial from fast-evolving potential energy terms
   if (use_virial) call zero_virial(viralt)

   do stepfast = 1, nalt2
      call integrate_vel( aalt2,dta_2 )
      if (use_rattle) call save_atoms_pos
      call integrate_pos( dta )
      if (use_rattle)  call rattle (dta)
!
!       Reassign the particules that have changed of domain
!
!       -> real space
      call reassignrespa(stepfast,nalt2)
!
!       communicate positions
!
      call commposrespa(stepfast.ne.nalt2)
!
      if (.not.ftot_l) then
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m(derivs,3*nbloc,rec_queue)
      end if
!
      if (stepfast.eq.nalt2) call mechanicsteprespa1(-1,0)
      if (stepfast.eq.nalt2) call allocsteprespa(.true.)
!
!       get the fast-evolving potential energy and atomic forces
!
      call gradfast(ealt,derivs)
!
!       communicate forces
!
      call comm_forces( derivs,cBond )
!
!       aMD/GaMD contributions
!
      call aMD (derivs,ealt)
!
!       MPI : get total energy
!
      call reduceen(ealt)
!
!     Debug print information
!
      if(deb_Energy) call info_energy(rank)
      if(deb_Force)  call info_forces(cBond)
      if(deb_Atom)   call info_minmax_pva
!
!     use Newton's second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the Verlet recursion
!
      call integrate_vel( derivs,aalt2,dta_2 )
      if (use_rattle)  call rattle2 (dta)
!
!       increment average virial from fast-evolving potential terms
!
      if (use_virial) then
!$acc parallel loop collapse(2) default(present) async
         do i = 1,3; do j = 1,3
               viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
            end do; end do
      end if
   end do
end
