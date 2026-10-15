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
   use atomsMirror, only: integrate_vel
   use cutoff, only: use_list
   use domdec, only: nbloc,rank
   use deriv, only: cNBond,comm_forces,ftot_l,info_forces
   use energi, only: calc_e,info_energy,reduceen
   use freeze, only: use_rattle
   use inform, only: abort,deb_Atom,deb_Energy,deb_Force,info_minmax_pva
   use mdstuf1, only: ealt,ekin,eksum,epot,etot,derivs,pres,stress,temp,viralt
   use parallel, only: commposrec
   use potent
   use utils, only: set_to_zero1m
   use tinMemory, only: prmem_requestm
   use mdevice, only: reciprocal_work
   use moldyn, only: a,dinter,dshort
   use tinheader, only: re_p
   use virial, only: use_virial,vir
   implicit none
   integer i,j,iglob
   integer istep
   real(r_p) dt,dt_2
   real(r_p) dta,dta_2,dta2
   real(8) time0,time1
!
!     set some time values for the dynamics integration
!
   dt_2  = 0.5_re_p * dt
   dta   = dt / dinter
   dta_2 = 0.5_re_p * dta
   dta2  = dta / dshort
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
      call set_to_zero1m(derivs,3*nbloc,reciprocal_work%queue)
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
!
!     subroutine respaint1 :
!     find intermediate-evolving velocities and positions via velocity Verlet recursion
!
subroutine respaint1(ealt,viralt,dta,dta2)
   use atomsMirror, only: integrate_vel
   use domdec, only: nbloc,rank
   use deriv, only: cSNBond,comm_forces,ftot_l,info_forces
   use energi, only: calc_e,info_energy,reduceen
   use freeze, only: use_rattle
   use inform, only: abort,deb_Atom,deb_Energy,deb_Force,info_minmax_pva
   use mdstuf1 ,only: ealt2,viralt2,derivs
   use moldyn, only: aalt,dinter,nalt,stepint
   use potent
   use utils, only: set_to_zero1m
   use tinMemory, only: prmem_requestm
   use mdevice, only: reciprocal_work
   use tinheader, only: re_p
   use virial, only: use_virial,vir,zero_virial
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
         call set_to_zero1m (derivs,3*nbloc,reciprocal_work%queue)
      end if
!
!     get the fast-evolving potential energy and atomic forces
!
      call gradient_helper(ealt,derivs,3,.false.)
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
   use atomsMirror, only: integrate_pos,integrate_vel,save_atoms_pos
   use domdec, only: nbloc,rank
   use deriv, only: cBond,comm_forces,ftot_l,info_forces
   use energi, only: info_energy,reduceen
   use freeze, only: use_rattle
   use inform, only: deb_Atom,deb_Energy,deb_Force,info_minmax_pva
   use mdstuf1 ,only: derivs
   use moldyn, only: aalt2,dshort,nalt2,stepfast
   use utils, only: set_to_zero1m
   use tinMemory, only: prmem_requestm
   use mdevice, only: reciprocal_work
   use parallel, only: commposrespa,reassign_respa
   use tinheader, only: re_p
   use virial, only: use_virial,vir,zero_virial
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
      call reassign_respa(stepfast,nalt2)
!
!       communicate positions
!
      call commposrespa(stepfast.ne.nalt2)
!
      if (.not.ftot_l) then
         call prmem_requestm(derivs,3,nbloc,async=.true.)
         call set_to_zero1m(derivs,3*nbloc,reciprocal_work%queue)
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
