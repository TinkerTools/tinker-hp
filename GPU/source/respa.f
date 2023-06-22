c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine respa  --  r-RESPA molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "respa" performs a single multiple time step molecular dynamics
c     step using the reversible reference system propagation algorithm
c     (r-RESPA) via a Verlet core with the potential split into fast-
c     and slow-evolving portions
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
c
#include "tinker_macro.h"
      subroutine respa(istep,dt)
      use ani
      use atmtyp
      use atomsMirror
      use bath     ,only: barostat
      use cutoff
      use domdec
      use deriv    ,only: info_forces,cNBond,cBond,ftot_l,comm_forces
     &             ,zero_forces_rec
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
      implicit none
      integer i,j,k,iglob
      integer istep
      logical save_pred
      real(r_p) dt,dt_2,dt_in
      real(r_p) dta,dta_2
c
c     set some time values for the dynamics integration
c
      dt_2  = 0.5_re_p * dt
      dta   = dt / dshort
      dta_2 = 0.5_re_p * dta

      if (use_ml_embedding) use_mlpot=.FALSE.
c
c     make half-step temperature and pressure corrections
c
c     call temper (dt)
c
c     store the current atom positions, then find half-step
c     velocities via velocity Verlet recursion
c
      call integrate_vel(a,dt_2,aalt,dta_2)
c
c     initialize virial from fast-evolving potential energy terms
c
      if (use_virial) call zero_virial(viralt)
c
c     find fast-evolving velocities and positions via Verlet recursion
c
      do stepfast = 1, nalt
         if (use_rattle) call save_atoms_pos
         call integrate_pos(dta)

         if (use_rattle) then
            call openacc_abort("Rattle not tested with openacc")
            call rattle (dta)
         endif
c
c     Reassign the particules that have changed of domain
c
c        -> real space
         call reassignrespa(stepfast,nalt)
c
c        communicate positions
c
         call commposrespa(stepfast.ne.nalt)
c
         if ( .not.ftot_l ) then
            call prmem_requestm(derivs,3,nbloc,async=.true.)
            call set_to_zero1m( derivs,3*nbloc,rec_queue)
         end if
c
         call mechanicsteprespa(istep,.true.)
         call allocsteprespa(.true.)
c
c     get the fast-evolving potential energy and atomic forces
c
         call gradfast (ealt,derivs)
c
c        communicate forces
c
         call comm_forces( derivs,cBond )
c
c        aMD/GaMD contributions
c
         call aMD (derivs,ealt)
c
c        MPI : get total energy
c
         call reduceen(ealt)
c
c        Debug information
c
         if(deb_Energy) call info_energy(rank)
         if(deb_Force)  call info_forces(cBond)
         if(abort)      call emergency_save
         if(abort)      call fatal

         dt_in = merge(dta_2,dta,use_rattle)  ! level 1
         dt_in = merge(dt_in,dta_2,stepfast.ne.nalt) ! level 0
         call integrate_vel(derivs,aalt,dt_in)
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
         if (use_rattle) then
            call openacc_abort("Rattle2 not tested with openacc")
            call rattle2 (dta)
            call integrate_vel(aalt,dta_2)
         end if
         if(deb_Atom)   call info_minmax_pva(1)
c
c     increment average virial from fast-evolving potential terms
c
         if (use_virial) then
!$acc parallel loop collapse(2) async default(present)
            do i = 1, 3; do j = 1, 3
               viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
            end do; end do
         end if
      end do
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      !call reassignrespa(nalt,nalt)
c
c     -> reciprocal space
c
      call reassignpme(.false.)
c
c     communicate positions
c
      call commposrespa(.false.)
      call commposrec
c
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
         call set_to_zero1m( derivs,3*nbloc,rec_queue)
      end if
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslow (epot,derivs)
c
c     if necessary, communicate some forces
c
      call comm_forces( derivs )
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
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using velocity Verlet recursion
c
      call integrate_vel(derivs,a,dt_2)
c
c     Debug print information
c
      if(deb_Energy) call info_energy(rank)
      if(deb_Force)  call info_forces(cNBond)
      if(deb_Atom)   call info_minmax_pva
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle) then
         call openacc_abort("Rattle2 not tested with openacc")
         call rattle2 (dt)
      end if
c
c     total potential and virial from sum of fast and slow parts
c
      if (calc_e) then
!$acc serial async present(epot,ealt)
         epot = epot + ealt
!$acc end serial
         if(rank.eq.0) call chk_energy_fluct(epot,ealt,abort)
      end if
c
      if (use_virial) then
!$acc parallel loop collapse(2) async default(present)
         do i = 1, 3; do j = 1, 3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do; end do
      end if
c
c     make full-step temperature and pressure corrections
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
      call pressure2 (epot,temp)
c
c     total energy is sum of kinetic and potential energies
c
      if (calc_e) then
!$acc serial async present(etot,eksum,epot)
         etot = eksum + epot
!$acc end serial
      end if

      ! Fatal Instructions
      if(abort)      call emergency_save
      if(abort)      call fatal
c
c     compute statistics and save trajectory for this step
c
      call mdsave     (istep,dt,epot)
      call mdrestgpu  (istep)
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      end
c
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradfast  --  fast energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradfast" calculates the potential energy and first derivatives
c     for the fast-evolving local valence potential energy terms
c
c
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
c
c     save the original state of slow-evolving potentials
c
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
c
c     turn off slow-evolving nonbonded potential energy terms
c
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
c
c     get energy and gradient for fast-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of slow-evolving potentials
c
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
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradslow  --  slow energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradslow" calculates the potential energy and first derivatives
c     for the slow-evolving nonbonded potential energy terms
c
c
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
c
      respa1_l = (index(integrate,'RESPA1').gt.0)
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
      save_mlpot  = use_mlpot

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
c
c        substract the previously stored short range energy and forces
c
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
c
c
c     #############################################################################
c     ##                                                                         ##
c     ##  subroutine gradembedding2  --  embedding energy & gradient components  ##
c     ##                                                                         ##
c     #############################################################################
c
c
c     "gradembedding2" calculates the potential energy and first derivatives
c     for the chosen embedding potential energy terms to substract for respa
c
c
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
c
c
c     save the original state of potential energy terms
c
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
c
c     turn on only chosen intra potential energy terms
c     that will be substract for respa
c
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
c
c     get energy and gradient for potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of potential energy term
c
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
