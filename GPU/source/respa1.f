c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine respa1  --  r-RESPA1 molecular dynamics step  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "respa1" performs a single multiple time step molecular dynamics
c     step using the reversible reference system propagation algorithm
c     (r-RESPA) via a Verlet core with the potential split into fast-
c     intermediate and slow-evolving portions
c
c     literature references:
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
c     Ruhong Zhou, Edward Harder, Huafeng Xu and B. J. Berne,
c     "Efficient multiple time step method for use with Ewald and
c     particle mesh Ewald for large biomolecular systems",
c     J. Chem. Phys. 115, 2348-2358 (2001)
c
c
#include "tinker_macro.h"
      subroutine respa1(istep,dt)
      use ani
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
c
c     set some time values for the dynamics integration
c
      dt_2  = 0.5_re_p * dt
      dta   = dt / dinter
      dta_2 = 0.5_re_p * dta
      dta2  = dta / dshort
      if(use_ml_embedding) use_mlpot=.FALSE.
c
c     store the current atom positions, then find half-step
c     velocities via velocity Verlet recursion
c
      call integrate_vel( a,dt_2 )
c
c     find intermediate-evolving velocities and positions via velocity Verlet recursion
c
      call respaint1(ealt,viralt,dta,dta2)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c      call reassignrespa(.false.,nalt,nalt)
c
c     -> reciprocal space
      call reassignpme(.false.)
c
c     communicate positions
c
c     call commposrespa(.false.)
      call commposrec
c
      call reinitnl(istep)
c
      call mechanicsteprespa1(istep,2)

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
c     communicate some forces
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
c     COMPUTE ML DELTA CONTRIBUTION (ml_embedding_mode=2)
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
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using velocity Verlet recursion
c
      call integrate_vel( derivs,a,dt_2 )
c
c     Debug print information
c
      if(deb_Energy) call info_energy(rank)
      if(deb_Force)  call info_forces(cNBond)
      if(deb_Atom)   call info_minmax_pva
      if(abort)      call fatal 
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     total potential and virial from sum of fast and slow parts
c
      if (calc_e.or.use_virial) then
!$acc serial async default(present) present(epot,ealt)
      epot = epot + ealt
      do i = 1,3; do j = 1,3
         vir(j,i) = vir(j,i) + viralt(j,i)
      end do; end do
!$acc end serial
      end if
c
c     make full-step temperature and pressure corrections
c
      call temper   (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
      if (calc_e) then
!$acc serial present(etot,eksum,epot) async
      etot = eksum + epot
!$acc end serial
      end if
c
c     compute statistics and save trajectory for this step
c
      call mdsave (istep,dt,epot)
      call mdrestgpu (istep)
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      end
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine gradint  --  intermediate energy & gradient components  ##
c     ##                                                                      ##
c     ##########################################################################
c
c
c     for the fast-evolving local valence potential energy terms
c
c
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
c
c     save the original state of fast-evolving potentials
c
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
c
c     turn off fast-evolving valence potential energy terms
c
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
c
c     get energy and gradient for slow-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of fast-evolving potentials
c
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
c
c     subroutine respaint1 :
c     find intermediate-evolving velocities and positions via velocity Verlet recursion
c
      subroutine respaint1(ealt,viralt,dta,dta2)
      use ani
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
c
c     initialize virial from fast-evolving potential energy terms
c
      if (use_virial) call zero_virial(viralt)

      do stepint = 1, nalt
         call integrate_vel( aalt,dta_2 )
c
c     find fast-evolving velocities and positions via velocity Verlet recursion
c
         call respafast1(ealt2,viralt2,dta2)
c
         call mechanicsteprespa1(-1,1)
c
         if (.not.ftot_l) then
            call prmem_requestm(derivs,3,nbloc,async=.true.)
            call set_to_zero1m (derivs,3*nbloc,rec_queue)
         end if
c
c     get the fast-evolving potential energy and atomic forces
c
         call gradint(ealt,derivs)
c
c      communicate forces
c
         call comm_forces( derivs,cSNBond )
c
c      MPI : get total energy
c
         call reduceen(ealt)
c
c     Debug print information
c
         if(deb_Energy) call info_energy(rank)
         if(deb_Force ) call info_forces(cSNBond)
         if(deb_Atom  ) call info_minmax_pva
         if(abort     ) __TINKER_FATAL__
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
         call integrate_vel( derivs,aalt,dta_2 )
         if (use_rattle)  call rattle2 (dta)
c
c     increment average virial from fast-evolving potential terms
c
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
c
c     subroutine respafast1 :
c     find fast-evolving velocities and positions via velocity Verlet recursion
c
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
c
c     initialize virial from fast-evolving potential energy terms
      if (use_virial) call zero_virial(viralt)

      do stepfast = 1, nalt2
         call integrate_vel( aalt2,dta_2 )
         if (use_rattle) call save_atoms_pos
         call integrate_pos( dta )
         if (use_rattle)  call rattle (dta)
c
c       Reassign the particules that have changed of domain
c
c       -> real space
        call reassignrespa(stepfast,nalt2)
c
c       communicate positions
c
        call commposrespa(stepfast.ne.nalt2)
c
        if (.not.ftot_l) then
           call prmem_requestm(derivs,3,nbloc,async=.true.)
           call set_to_zero1m(derivs,3*nbloc,rec_queue)
        end if
c
        if (stepfast.eq.nalt2) call mechanicsteprespa1(-1,0)
        if (stepfast.eq.nalt2) call allocsteprespa(.true.)
c
c       get the fast-evolving potential energy and atomic forces
c
        call gradfast(ealt,derivs)
c
c       communicate forces
c
        call comm_forces( derivs,cBond )
c
c       aMD/GaMD contributions
c
        call aMD (derivs,ealt)
c
c       MPI : get total energy
c
        call reduceen(ealt)
c
c     Debug print information
c
        if(deb_Energy) call info_energy(rank)
        if(deb_Force)  call info_forces(cBond)
        if(deb_Atom)   call info_minmax_pva
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
        call integrate_vel( derivs,aalt2,dta_2 )
        if (use_rattle)  call rattle2 (dta)
c
c       increment average virial from fast-evolving potential terms
c
        if (use_virial) then
!$acc parallel loop collapse(2) default(present) async
           do i = 1,3; do j = 1,3
              viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
           end do; end do
        end if
      end do
      end
