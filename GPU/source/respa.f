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
#include "tinker_precision.h"
      module respa_mod
         real(r_p),allocatable::derivs(:,:)
      end module

      subroutine respa(istep,dt)
      use atmtyp
      use atomsMirror
      use cutoff
      use domdec
      use deriv  ,only: info_forces,cNBond,cBond
      use energi ,only: info_energy
      use freeze
      use inform
      use moldyn
      use respa_mod
      use timestat
      use units
      use usage
      use utils
      use utilgpu,only:prmem_requestm
     &           ,openacc_abort
      use virial
      use mpi
      implicit none
      integer i,j,k,iglob
      integer istep
      real(r_p) dt,dt_2
      real(r_p) dta,dta_2
      real(r_p),save:: epot,etot
      real(r_p),save:: eksum
      real(r_p),save:: temp,pres
      real(r_p),save:: ealt
      real(r_p),save:: ekin(3,3)
      real(r_p),save:: stress(3,3)
      real(r_p),save:: viralt(3,3)
      logical  ,save:: f_in=.true.
c
c     set some time values for the dynamics integration
c
      dt_2  = 0.5_re_p * dt
      dta   = dt / dshort
      dta_2 = 0.5_re_p * dta

      if (f_in) then
!$acc enter data create(ekin,ealt,stress,viralt)
!$acc&           create(etot,epot,eksum,temp,pres)
         f_in = .false.
      end if

!$acc data present(ekin,ealt,stress,viralt,etot,epot,eksum,
!$acc&         temp,pres)
!$acc&     present(x,y,z,v,a,glob,use,aalt,mass,
!$acc&         xold,yold,zold,vir)
c
c     make half-step temperature and pressure corrections
c
c     call temper (dt)
c
c     store the current atom positions, then find half-step
c     velocities via velocity Verlet recursion
c
!$acc parallel loop collapse(2) async
      do i = 1, nloc
         do j = 1, 3
            iglob = glob(i)
            if (use(iglob))
     &         v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
         end do
      end do
c
c     initialize virial from fast-evolving potential energy terms
c
!$acc parallel loop collapse(2) async
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0_re_p
         end do
      end do
c
c     find fast-evolving velocities and positions via Verlet recursion
c
      do k = 1, nalt
!$acc parallel loop async
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
!$acc loop seq
               do j = 1, 3
                  v(j,iglob) = v(j,iglob) + aalt(j,iglob)*dta_2
               end do
               xold(iglob) = x(iglob)
               yold(iglob) = y(iglob)
               zold(iglob) = z(iglob)
               x(iglob) = x(iglob) + v(1,iglob)*dta
               y(iglob) = y(iglob) + v(2,iglob)*dta
               z(iglob) = z(iglob) + v(3,iglob)*dta
            end if
         end do
         if (use_rattle) then
            call openacc_abort("Rattle not tested with openacc")
            call rattle (dta)
         endif
c
c     Reassign the particules that have changed of domain
c
c        -> real space
         call reassignrespa(k,nalt)
c
c        communicate positions
c
         call commposrespa(k.ne.nalt)
         call reCast_position
c
         call prmem_requestm(derivs,3,nbloc,async=.true.)
c
!$acc parallel loop collapse(2) present(derivs) async
         do i = 1, nbloc
            do j = 1, 3
               derivs(j,i) = 0.0_re_p
            end do
         end do
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
         call commforcesrespa(derivs,.true.)
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
         if(deb_Atom)   call info_minmax_pva
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
!$acc parallel loop collapse(2) present(derivs) async
         do i = 1, nloc
            do j = 1, 3
               iglob = glob(i)
               if (use(iglob)) then
                  aalt(j,iglob) = -convert * derivs(j,i) / mass(iglob)
                  v(j,iglob)    = v(j,iglob) + aalt(j,iglob)*dta_2
               end if
            end do
         end do
         if (use_rattle) then
            call openacc_abort("Rattle2 not tested with openacc")
            call rattle2 (dta)
         end if
c
c     increment average virial from fast-evolving potential terms
c
!$acc parallel loop collapse(2) async
         do i = 1, 3
            do j = 1, 3
               viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
            end do
         end do
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
      call reCast_position
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
      call prmem_requestm(derivs,3,nbloc,async=.true.)
c
!$acc parallel loop collapse(2) present(derivs) async
      do i = 1,nbloc
         do j = 1,3
            derivs(j,i) = 0.0_re_p
         end do
      end do
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslow (epot,derivs)
c
c     if necessary, communicate some forces
c
      call commforcesrespa(derivs,.false.)
c
c     MPI : get total energy
c
      call reduceen(epot)
c
c     Debug print information
c
      if(deb_Energy) call info_energy(rank)
      if(deb_Force)  call info_forces(cNBond)
      if(deb_Atom)   call info_minmax_pva
c
c     make half-step temperature and pressure corrections
c
      call temper2 (temp)
      call pressure2 (epot,temp)
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using velocity Verlet recursion
c
!$acc parallel loop collapse(2) present(derivs) async
      do i = 1, nloc
         do j = 1, 3
            iglob = glob(i)
            if (use(iglob)) then
               a(j,iglob) = -convert * derivs(j,i) / mass(iglob)
               v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
            end if
         end do
      end do
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
!$acc serial async
      epot = epot + ealt
!$acc end serial
c
!$acc parallel loop collapse(2) async
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do
      end do
c
c     make full-step temperature and pressure corrections
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
!$acc serial async
      etot = eksum + epot
!$acc end serial
c
c     compute statistics and save trajectory for this step
c
      call mdsave     (istep,dt,epot)
      call mdrestgpu  (istep)
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
!$acc end data
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
      implicit none
      real(r_p) energy
      real(r_p) derivs(3,*)
      logical save_vdw,save_charge
      logical save_dipole
      logical save_mpole,save_polar
      logical save_rxnfld,save_solv
      logical save_list
      logical save_smdvel, save_smdfor
c
c
c     save the original state of slow-evolving potentials
c
      save_vdw    = use_vdw
      save_charge = use_charge
      save_mpole  = use_mpole
      save_polar  = use_polar
      save_solv   = use_solv
      save_list   = use_list
      save_smdvel = use_smd_velconst
      save_smdfor = use_smd_forconst
c
c     turn off slow-evolving nonbonded potential energy terms
c
      use_vdw    = .false.
      use_charge = .false.
      use_mpole  = .false.
      use_polar  = .false.
      use_solv   = .false.
      use_list   = .false.
      use_smd_velconst = .false.
      use_smd_forconst = .false.
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
      use_solv   = save_solv
      use_list   = save_list
      use_smd_velconst = save_smdvel
      use_smd_forconst = save_smdfor
      return
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
      use potent
      implicit none
      real(r_p) energy
      real(r_p) derivs(3,*)
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_tortor,save_geom
      logical save_metal,save_extra
c
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
      save_tortor = use_tortor
      save_geom   = use_geom
      save_extra  = use_extra
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
      use_tortor = .false.
      use_geom   = .false.
      use_extra  = .false.
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
      use_tortor = save_tortor
      use_geom   = save_geom
      use_extra  = save_extra
      return
      end
c
