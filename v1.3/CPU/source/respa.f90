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
!> @brief 
!> "respa" performs a single multiple time step molecular dynamics
!> step using the reversible reference system propagation algorithm
!> (r-RESPA) via a Verlet core with the potential split into fast-
!> and slow-evolving portions
!> @param[in] istep: index of timestep
!> @param[in] dt: timestep
subroutine respa(istep,dt)
   use atmtyp
   use atoms
   use cutoff
   use deriv
   use domdec
   use freeze
   use inform
   use iounit
   use moldyn
   use timestat
   use units
   use usage
   use virial
   use mpi
   use spectra
   implicit none
   integer i,j,iglob
   integer istep
   real*8 dt,dt_2
   real*8 dta,dta_2
   real*8 epot,etot
   real*8 eksum
   real*8 temp,pres
   real*8 ealt
   real*8 ekin(3,3)
   real*8 stress(3,3)
   real*8 viralt(3,3)
   real*8 time0,time1
   real*8, allocatable :: derivs(:,:)
   real*8 dip(3,nalt),dipind(3,nalt)
   time0 = mpi_wtime()
!
   if (deb_Path) write(iout,*), 'respa '
!
!
!     set some time values for the dynamics integration
!
   dt_2 = 0.5d0 * dt
   dta = dt / dshort
   dta_2 = 0.5d0 * dta
!
!     store the current atom positions, then find half-step
!     velocities via velocity Verlet recursion
!
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         do j = 1, 3
            v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
         end do
      end if
   end do
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
!     find fast-evolving velocities and positions via velocity Verlet recursion
!
   call respafast(ealt,viralt,dta,istep,dip,dipind)
!
!     Reassign the particules that have changed of domain
!
!     -> real space
!
   time0 = mpi_wtime()
   call reassignrespa(nalt,nalt)
!
!     -> reciprocal space
!
   call reassignpme(.false.)
   time1 = mpi_wtime()
   timereneig = timereneig + time1 - time0
!
!     communicate positions
!
   time0 = mpi_wtime()
   call commposrespa(.false.)
   call commposrec
   time1 = mpi_wtime()
   timecommpos = timecommpos + time1 - time0
!
   time0 = mpi_wtime()
   call reinitnl(istep)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
   time0 = mpi_wtime()
   call mechanic_up_para_respa(istep,.false.)
   time1 = mpi_wtime()
   timeparam = timeparam + time1 - time0

   time0 = mpi_wtime()
   call allocsteprespa(.false.)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
!     rebuild the neighbor lists
!
   time0 = mpi_wtime()
   if (use_list) call nblist(istep)
   time1 = mpi_wtime()
   timenl = timenl + time1 - time0
!
   time0 = mpi_wtime()
   allocate (derivs(3,nbloc))
   derivs = 0d0
   time1 = mpi_wtime()
   timeinte = timeinte + time1 - time0
!
!     get the slow-evolving potential energy and atomic forces
!
   time0 = mpi_wtime()
   call gradslow (epot,derivs)
   time1 = mpi_wtime()
   timegrad = timegrad + time1 - time0
!
!     if necessary, communicate some forces
!
   time0 = mpi_wtime()
   call commforcesrespa(derivs,.false.)
   time1 = mpi_wtime()
   timecommforces = timecommforces + time1-time0
!
!     MPI : get total energy
!
   time0 = mpi_wtime()
   call reduceen(epot)
   time1 = mpi_wtime()
   timered = timered + time1 - time0
!
!     Debug print information
!
   if(deb_Energy) call info_energy(rank)
   if(deb_Force)  call info_forces(cNBond)
   if(deb_Atom)   call info_minmax_pva
!
!     make half-step temperature and pressure corrections
!
   time0 = mpi_wtime()
   call temper2 (temp)
   call pressure2 (epot,temp)
   time1 = mpi_wtime()
   timetp = timetp + time1-time0
!
!     use Newton's second law to get the slow accelerations;
!     find full-step velocities using velocity Verlet recursion
!
   time0 = mpi_wtime()
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         do j = 1, 3
            a(j,iglob) = -convert * derivs(j,i) / mass(iglob)
            v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
         end do
      end if
   end do
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle)  call rattle2 (dt)
!
!     total potential and virial from sum of fast and slow parts
!
   epot = epot + ealt
   do i = 1, 3
      do j = 1, 3
         vir(j,i) = vir(j,i) + viralt(j,i)
      end do
   end do
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
!     make full-step temperature and pressure corrections
!
   time0 = mpi_wtime()
   call temper (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
   if(ir) then
      call compute_dipole(dip(:,nalt)&
      &,dipind(:,nalt),full_dipole)
      call save_dipole_respa(dip,dipind)
   endif
   time1 = mpi_wtime()
   timetp = timetp + time1-time0
!
!     total energy is sum of kinetic and potential energies
!
   time0 = mpi_wtime()
   etot = eksum + epot
!
!     compute statistics and save trajectory for this step
!
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)
   call mdsave (istep,dt,epot,derivs)
   call mdrest (istep)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
!     perform deallocation of some local arrays
!
   deallocate (derivs)
   return
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
!> @brief 
!> calculates the potential energy and first derivatives
!> for the fast-evolving local valence potential energy terms
!> @param[in] energy: potential energy associated to fast potential terms
!> @param[in] energy: derivatives of potential energy associated to fast potential terms
subroutine gradfast (energy,derivs)
   use cutoff
   use inform
   use iounit
   use potent
#ifdef COLVARS
   use colvars
#endif
#ifdef PLUMED
   use plumed
#endif
   implicit none
   real*8 energy
   real*8 derivs(3,*)
   logical save_vdw,save_charge
   logical save_mpole,save_polar
   logical save_repuls,save_disp,save_chgtrn
   logical save_list
#ifdef COLVARS
   logical save_colvars
#endif
#ifdef PLUMED
   logical save_plumed
#endif
!
   if (deb_Path) write(iout,*), 'gradfast '
!
!
!     save the original state of slow-evolving potentials
!
   save_vdw = use_vdw
   save_charge = use_charge
   save_mpole = use_mpole
   save_polar = use_polar
   save_repuls = use_repuls
   save_disp = use_disp
   save_chgtrn = use_chgtrn
   save_list = use_list
#ifdef COLVARS
   save_colvars = use_colvars
#endif
#ifdef PLUMED
   save_plumed = lplumed
#endif
!
!     turn off slow-evolving nonbonded potential energy terms
!
   use_vdw = .false.
   use_charge = .false.
   use_mpole = .false.
   use_polar = .false.
   use_repuls = .false.
   use_disp = .false.
   use_chgtrn = .false.
   use_list = .false.
#ifdef COLVARS
   use_colvars = .false.
#endif
#ifdef PLUMED
   lplumed = .false.
#endif
!
!     get energy and gradient for fast-evolving potential terms
!
   call gradient (energy,derivs)
!
!     restore the original state of slow-evolving potentials
!
   use_vdw = save_vdw
   use_charge = save_charge
   use_mpole = save_mpole
   use_polar = save_polar
   use_repuls = save_repuls
   use_disp = save_disp
   use_chgtrn = save_chgtrn
   use_list = save_list
#ifdef COLVARS
   use_colvars = save_colvars
#endif
#ifdef PLUMED
   lplumed = save_plumed
#endif
   return
end
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
!> @brief 
!> calculates the potential energy and first derivatives
!> for the slow-evolving nonbonded potential energy terms
!> @param[in] energy: potential energy associated to slow potential terms
!> @param[in] energy: derivatives of potential energy associated to slow potential terms
subroutine gradslow (energy,derivs)
   use inform
   use iounit
   use potent
   implicit none
   real*8 energy
   real*8 derivs(3,*)
   logical save_bond,save_angle
   logical save_strbnd,save_urey
   logical save_angang,save_opbend
   logical save_opdist,save_improp
   logical save_imptor,save_tors
   logical save_pitors,save_strtor
   logical save_angtor
   logical save_tortor,save_geom
   logical save_extra
!
   if (deb_Path) write(iout,*), 'gradslow '
!
!
!
!     save the original state of fast-evolving potentials
!
   save_bond = use_bond
   save_angle = use_angle
   save_strbnd = use_strbnd
   save_urey = use_urey
   save_angang = use_angang
   save_opbend = use_opbend
   save_opdist = use_opdist
   save_improp = use_improp
   save_imptor = use_imptor
   save_tors = use_tors
   save_pitors = use_pitors
   save_strtor = use_strtor
   save_angtor = use_angtor
   save_tortor = use_tortor
   save_geom = use_geom
   save_extra = use_extra
!
!     get energy and gradient for slow-evolving potential terms
!
!     turn off fast-evolving valence potential energy terms
!
   use_bond = .false.
   use_angle = .false.
   use_strbnd = .false.
   use_urey = .false.
   use_angang = .false.
   use_opbend = .false.
   use_opdist = .false.
   use_improp = .false.
   use_imptor = .false.
   use_tors = .false.
   use_pitors = .false.
   use_strtor = .false.
   use_angtor = .false.
   use_tortor = .false.
   use_geom = .false.
   use_extra = .false.
!
!     get energy and gradient for slow-evolving potential terms
!
   call gradient (energy,derivs)
!
!     restore the original state of fast-evolving potentials
!
   use_bond = save_bond
   use_angle = save_angle
   use_strbnd = save_strbnd
   use_urey = save_urey
   use_angang = save_angang
   use_opbend = save_opbend
   use_opdist = save_opdist
   use_improp = save_improp
   use_imptor = save_imptor
   use_tors = save_tors
   use_pitors = save_pitors
   use_strtor = save_strtor
   use_angtor = save_angtor
   use_tortor = save_tortor
   use_geom = save_geom
   use_extra = save_extra
!
!     get energy and gradient for slow-evolving potential terms
   return
end
!
!     subroutine respafast :
!     find fast-evolving velocities and positions via velocity Verlet recursion
!
!> @brief 
!> find fast-evolving velocities and positions via velocity Verlet recursion
!> computes and displays the total potential energy
!> @param[in] ealt: local potential energy
!> @param[in] viralt: local virial
!> @param[in] dta: intermediate timestep
!> @param[in] istep: index of timestep
!> @param[in] dip: dipole moment
!> @param[in] dipind: induced dipole moment
subroutine respafast(ealt,viralt,dta,istep,dip,dipind)
   use atmtyp
   use atoms
   use cutoff
   use deriv
   use domdec
   use freeze
   use inform
   use iounit
   use moldyn
   use timestat
   use units
   use usage
   use virial
   use mpi
   use spectra
   implicit none
   integer i,j,k,iglob
   integer istep
   real*8 dip(3,nalt),dipind(3,nalt)
   real*8 dta,dta_2
   real*8 ealt
   real*8 time0,time1
   real*8, allocatable :: derivs(:,:)
   real*8 viralt(3,3)
   time0 = mpi_wtime()
!
   if (deb_Path) write(iout,*), 'respafast '
!

   dta_2 = 0.5d0 * dta
!
!     initialize virial from fast-evolving potential energy terms
!
   do i = 1, 3
      do j = 1, 3
         viralt(j,i) = 0.0d0
      end do
   end do
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0

   do k = 1, nalt
      time0 = mpi_wtime()
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
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
      if (use_rattle)  call rattle (dta)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
      if(deb_Atom)   call info_minmax_pva(1)
!
!       Reassign the particules that have changed of domain
!
!       -> real space
!
      time0 = mpi_wtime()
!
      call reassignrespa(k,nalt)
      time1 = mpi_wtime()
      timereneig = timereneig + time1 - time0
!
!
!       communicate positions
!
      time0 = mpi_wtime()
      call commposrespa(k.ne.nalt)
      time1 = mpi_wtime()
      timecommpos = timecommpos + time1 - time0
!
      time0 = mpi_wtime()
      allocate (derivs(3,nbloc))
      derivs = 0d0
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
!
      time0 = mpi_wtime()
      call mechanic_up_para_respa(istep,.true.)
      call allocsteprespa(.true.)
      time1 = mpi_wtime()
      timeparam = timeparam + time1 - time0
!
!     get the fast-evolving potential energy and atomic forces
!
      time0 = mpi_wtime()
      call gradfast (ealt,derivs)
      time1 = mpi_wtime()
      timegrad = timegrad + time1 - time0
!
!       communicate forces
!
      time0 = mpi_wtime()
      call commforcesrespa(derivs,.true.)
      time1 = mpi_wtime()
      timecommforces = timecommforces + time1-time0
!
!       MPI : get total energy
!
      time0 = mpi_wtime()
      call reduceen(ealt)
      time1 = mpi_wtime()
      timered = timered + time1 - time0
!
!        Debug information
!
      if(deb_Energy) call info_energy(rank)
      if(deb_Force)  call info_forces(cBond)
!
!     use Newton's second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the Verlet recursion
!
      time0 = mpi_wtime()
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               aalt(j,iglob) = -convert *&
               &derivs(j,i) / mass(iglob)
               v(j,iglob) = v(j,iglob) + aalt(j,iglob)*dta_2
            end do
         end if
      end do
      deallocate (derivs)
      if (use_rattle)  call rattle2 (dta)
!
!     increment average virial from fast-evolving potential terms
!
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = viralt(j,i) + vir(j,i)/dshort
         end do
      end do

      if(ir .and. k<nalt) then
         !call rotpole
         call compute_dipole(dip(:,k)&
         &,dipind(:,k),.FALSE.)
      endif
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
   end do
   return
end
