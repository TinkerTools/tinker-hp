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
   use atmtyp, only: mass
   use atomsMirror, only: integrate_pos,integrate_vel,save_atoms_pos
   use cutoff, only: use_list
   use domdec, only: nbloc,nproc,rank
   use deriv   ,only: info_forces,prtEForces,cDef,ftot_l,comm_forces
   use energi  ,only: info_energy,calc_e,chk_energy_fluct,reduceen
   use freeze, only: use_rattle
   use inform, only: abort,deb_Atom,deb_Energy,deb_Force,info_minmax_pva
   use mdevice, only: reciprocal_work
   use mdstuf1, only: derivs,ekin,eksum,epot,etot,pres,stress,temp
   use moldyn, only: a
   use parallel,only: commpos,commposrec,migrate_atoms
   use potent, only: use_embd_potoff,use_ml_embedding
   use utils    ,only:set_to_zero1m
   use tinMemory, only: prmem_requestm
   use sizes, only: tinkerdebug
   use virial, only: use_virial
   use spectra, only: compute_dipole,full_dipole,ir,save_dipole_traj
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
   call migrate_atoms
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
      call set_to_zero1m (derivs,3*nbloc,reciprocal_work%queue)
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

   if (istep.eq.1 .and. nproc.eq.1 .and. .not.use_embd_potoff) then
      if (use_ml_embedding .and. tinkerdebug.ne.0)&
         &call audit_embedding_forces
   endif
!
!     get the potential energy and atomic forces
!
   if (use_ml_embedding) then
      call gradient_embd(epot,derivs)
   else
      call gradient (epot,derivs)
   endif
!
!     communicate forces
!
   call comm_forces( derivs )
!
!     MPI : get total energy
!
   if (calc_e.or.use_virial) call reduceen(epot)
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
!     #########################################################################
!     ##                                                                     ##
!     ##  subroutine gradient_embedding_verlet  -- explicit hybrid assembly  ##
!     ##                                                                     ##
!     #########################################################################
!
!
!     Build A(P union W) and the standalone Delta ML(P) correction separately.
!     This is required because mutual polarization makes group-weighted AMOEBA
!     forces non-linear in the group weights.
!
subroutine gradient_embd(energy,derivs)
   use deriv, only: ftot_l,restore_total_forces,save_total_forces,zero_forces
   use energi, only: calc_e,esave
   use mlpot, only: ml_embedding_mode,mlpotscale
   use potent, only: use_embd_potoff,use_mlpot
   use mdstuf1, only: viralt
   use mdevice, only: reciprocal_work
   use uprior, only: use_pred
   use virial, only: use_virial,vir
   implicit none
   integer i,j
   real(r_p) energy
   real(r_p) derivs(3,*)
   real(t_p) saved_scale
   integer saved_mode
   logical saved_use_mlpot,saved_pred

   if (.not.ftot_l) then
      write(0,*) 'ML embedding audit requires the de_tot force path'
      __TINKER_FATAL__
   endif

   saved_mode = ml_embedding_mode
   saved_scale = mlpotscale
   saved_use_mlpot = use_mlpot
   saved_pred = use_pred

   call zero_forces

   ! Full AMOEBA(P union W), without the ML component.
   use_mlpot = .false.
   ml_embedding_mode = 0
   call gradient(energy,derivs)
   if (calc_e.or.use_virial) then
!$acc serial present(energy,esave) async(reciprocal_work%queue)
      esave = energy
!$acc end serial
   endif
   if (use_virial) then
!$acc parallel loop collapse(2) async(reciprocal_work%queue) present(vir,viralt)
      do i = 1,3; do j = 1,3
         viralt(j,i) = vir(j,i)
      end do; end do
   end if
   call save_total_forces

   ! Delta ML(P) = ML(P) - AMOEBA(P), evaluated by the existing mode-2 path.
   ! The potential-off mode uses the same mode-2 correction semantics.
   use_mlpot = .true.
   use_pred = .false.
   if (use_embd_potoff) then
      ml_embedding_mode = 2
      call gradembedding2(energy,derivs)
   else
      ml_embedding_mode = 2
      call gradient(energy,derivs)
   endif
   if (calc_e.or.use_virial) then
!$acc serial present(energy,esave) async(reciprocal_work%queue)
      energy = energy + esave
!$acc end serial
   endif
   if (use_virial) then
!$acc parallel loop collapse(2) async(reciprocal_work%queue) present(vir,viralt)
      do i = 1,3; do j = 1,3
         viralt(j,i) = viralt(j,i) + vir(j,i)
         vir(j,i) = viralt(j,i)
      end do; end do
   end if
   call restore_total_forces

   ml_embedding_mode = saved_mode
   mlpotscale = saved_scale
   use_mlpot = saved_use_mlpot
   use_pred = saved_pred
end

!
!     #########################################################################
!     ##                                                                     ##
!     ##  subroutine audit_embedding_forces  -- verify the hybrid derivative ##
!     ##                                                                     ##
!     #########################################################################
!
!
!     Compare the production Verlet derivative with the embedding identity at
!     one frozen, one-rank snapshot.  The final gradient call deliberately
!     rebuilds the production buffers before returning to the integrator.
!
subroutine audit_embedding_forces
   use deriv, only: ftot_l,get_ftot,zero_forces
   use domdec, only: glob,nbloc,nproc,rank
   use mdstuf1, only: epot
   use mdevice, only: reciprocal_work
   use mlpot, only: ml_embedding_mode,mlpotscale
   use potent, only: use_ml_embedding,use_mlpot
   use group, only: grplist,wgrp
   use tinheader, only: re_p
   implicit none
   integer i,j,n
   integer saved_mode
   real(r_p) max_expected,max_abs,max_abs_p,max_abs_w
   real(r_p) max_rel,max_rel_p,max_rel_w
   real(r_p) ml_w_norm,tolerance
   real(r_p) sum_expected(3),sum_production(3)
   real(r_p), allocatable :: full(:,:),isolated(:,:),ml(:,:)
   real(r_p), allocatable :: expected(:,:),production(:,:)
   real(r_p), allocatable :: snapshot(:,:)
   real(t_p), allocatable :: wgrp_saved(:,:)
   real(t_p) saved_scale
   logical saved_use_mlpot,is_p

   if (rank.ne.0 .or. nproc.ne.1 .or. .not.use_ml_embedding) return
   if (.not.ftot_l) then
      write(0,*) 'ML embedding audit requires the de_tot force path'
      return
   endif
   n = nbloc
   allocate(full(3,n),isolated(3,n),ml(3,n))
   allocate(expected(3,n),production(3,n))
   allocate(snapshot(3,n))
   allocate(wgrp_saved(size(wgrp,1),size(wgrp,2)))
   wgrp_saved = wgrp
   saved_mode = ml_embedding_mode
   saved_scale = mlpotscale
   saved_use_mlpot = use_mlpot

!$acc data create(snapshot)

   ! A(P union W).
   use_mlpot = .false.
   call zero_forces
   call gradient(epot,snapshot)
   call get_ftot(snapshot,nbloc)
   !$acc update host(snapshot) async(reciprocal_work%queue)
   !$acc wait(reciprocal_work%queue)
   full = snapshot

   ! A(P), recovered from the AMOEBA part of the existing mode-2 correction.
   ml_embedding_mode = 2
   call set_embedding_weights
   call zero_forces
   call gradient(epot,snapshot)
   call get_ftot(snapshot,nbloc)
   !$acc update host(snapshot) async(reciprocal_work%queue)
   !$acc wait(reciprocal_work%queue)
   isolated = -snapshot

   ! The mode-2 call supplies both ML(P) in dmlpot and the correction path.
   wgrp = wgrp_saved
   !$acc update device(wgrp) async(reciprocal_work%queue)
   use_mlpot = .true.
   ml_embedding_mode = 2
   call zero_forces
   call gradient(epot,snapshot)
   call get_ftot(snapshot,nbloc)
   !$acc update host(snapshot) async(reciprocal_work%queue)
   !$acc wait(reciprocal_work%queue)
   ml = snapshot + isolated

   ! Refresh the full reference after the correction evaluation changed the
   ! iterative polarization state, then rebuild the non-MTS force explicitly.
   wgrp = wgrp_saved
   !$acc update device(wgrp) async(reciprocal_work%queue)
   use_mlpot = .false.
   ml_embedding_mode = 0
   call zero_forces
   call gradient(epot,snapshot)
   call get_ftot(snapshot,nbloc)
   !$acc update host(snapshot) async(reciprocal_work%queue)
   !$acc wait(reciprocal_work%queue)
   full = snapshot
   use_mlpot = .true.
   call zero_forces
   call gradient_embd(epot,snapshot)
   call get_ftot(snapshot,nbloc)
   !$acc update host(snapshot) async(reciprocal_work%queue)
   !$acc wait(reciprocal_work%queue)
   production = snapshot

   expected = full + ml - isolated
   max_expected = maxval(abs(expected))
   max_abs = 0.0_re_p
   max_abs_p = 0.0_re_p
   max_abs_w = 0.0_re_p
   max_rel = 0.0_re_p
   max_rel_p = 0.0_re_p
   max_rel_w = 0.0_re_p
   ml_w_norm = 0.0_re_p
   sum_expected = 0.0_re_p
   sum_production = 0.0_re_p
   do i = 1,n
      is_p = .false.
      if (glob(i).le.size(grplist)) is_p = grplist(glob(i)).eq.1
      do j = 1,3
         max_abs = max(max_abs,abs(production(j,i)-expected(j,i)))
         max_rel = max(max_rel,abs(production(j,i)-expected(j,i))&
            &/max(1.0e-12_re_p,abs(expected(j,i))))
         if (is_p) then
            max_abs_p = max(max_abs_p,abs(production(j,i)-expected(j,i)))
            max_rel_p = max(max_rel_p,abs(production(j,i)-expected(j,i))&
               &/max(1.0e-12_re_p,abs(expected(j,i))))
         else
            max_abs_w = max(max_abs_w,abs(production(j,i)-expected(j,i)))
            max_rel_w = max(max_rel_w,abs(production(j,i)-expected(j,i))&
               &/max(1.0e-12_re_p,abs(expected(j,i))))
            ml_w_norm = ml_w_norm + ml(j,i)*ml(j,i)
         endif
         sum_expected(j) = sum_expected(j) + expected(j,i)
         sum_production(j) = sum_production(j) + production(j,i)
      enddo
   enddo
   ml_w_norm = sqrt(ml_w_norm)
   tolerance = 1.0e-5_re_p * (1.0_re_p + max_expected)
   write(*,'(A)') 'ML embedding derivative audit (dV/dx, before acceleration):'
   write(*,'(A,2ES16.8)') '  all max abs/rel: ',max_abs,max_rel
   write(*,'(A,2ES16.8)') '  P   max abs/rel: ',max_abs_p,max_rel_p
   write(*,'(A,2ES16.8)') '  W   max abs/rel: ',max_abs_w,max_rel_w
   write(*,'(A,ES16.8)') '  ML derivative norm on W: ',ml_w_norm
   write(*,'(A,3ES16.8)') '  expected derivative sum: ',sum_expected
   write(*,'(A,3ES16.8)') '  de_tot derivative sum:   ',sum_production
   write(*,'(A,ES16.8)') '  tolerance:               ',tolerance

   wgrp = wgrp_saved
   !$acc update device(wgrp) async(reciprocal_work%queue)
   ml_embedding_mode = saved_mode
   mlpotscale = saved_scale
   use_mlpot = saved_use_mlpot
!$acc end data
   deallocate(full,isolated,ml,expected,production,snapshot,wgrp_saved)
   call zero_forces

   if (max_abs.gt.tolerance .or. ml_w_norm.gt.tolerance) then
      write(0,*) 'ML embedding derivative audit failed'
      __TINKER_FATAL__
   endif
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
