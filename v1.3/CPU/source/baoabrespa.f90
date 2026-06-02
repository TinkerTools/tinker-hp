!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  subroutine baoab  --  baoab Langevin molecular dynamics step  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     "baoab" performs a single molecular dynamics time step
!     via the respa-baoab recursion formula
!
!     literature references:
!
!     Efficient molecular dynamics using geodesic integration
!     and solvent-solute splitting B. Leimkuhler and C. Matthews,
!     Proceedings of the Royal Society A, 472: 20160138, 2016
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
!> "baoabrespa" performs a single PIMD time step
!> via the BAOAB-RESPA recursion formula
!> @param[in] istep: number of timestep
!> @param[in] dt: length of timestep
subroutine baoabrespa (istep,dt)
   use atmtyp
   use atoms
   use bath
   use cutoff
   use domdec
   use deriv
   use energi
   use freeze
   use inform
   use iounit
   use mdstuf
   use moldyn
   use mpi
   use timestat
   use tortor
   use units
   use usage
   use virial
   use spectra
   use utilbaoab
   use mdstuf1
   implicit none
   real*8, intent(in) :: dt
   integer, intent(in) :: istep
   integer i,iglob,stepfast
   real*8 dta,dta_2,dt_2
   real*8 time0,time1
   real*8 dip(3,nalt),dipind(3,nalt)
!
   if (deb_Path) write(iout,*), 'baoabrespa '
!

   if (istep.eq.1) then
      pres=atmsph
   end if
   time0 = mpi_wtime()
!
!     set some time values for the dynamics integration
!
   dt_2  = 0.5d0 * dt
   dta   = dt / dshort
   dta_2 = 0.5d0 * dta
!
!     find quarter step velocities and half step positions via baoab recursion
!
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         v(:,iglob) = v(:,iglob) + dt_2*a(:,iglob)
      end if
   end do
!
   if (use_rattle) call rattle2(dt_2)
!     initialize virial from fast-evolving potential energy terms
   viralt(:,:)=0.0d0

   if(use_piston) then
      call apply_a_piston(dt_2,-1,.FALSE.)
      call apply_o_piston(dt)
   endif
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
!     respa inner loop
!
   do stepfast = 1, nalt
!
!     find fast-evolving velocities and positions via BAOAB recursion
!
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            v(:,iglob) = v(:,iglob) + dta_2*aalt(:,iglob)
         end if
      end do
!
      if (use_rattle)  call rattle2 (dta_2)
!
      call apply_a_block(dta_2)
      call apply_o_block(dta)
      call apply_a_block(dta_2)
!
!       Reassign the particules that have changed of domain
!
!       -> real space
      call reassignrespa(stepfast,nalt)
!
!       communicate positions
!
      call commposrespa(stepfast.ne.nalt)
      !call commposrespa(.TRUE.)
!
      if(allocated(derivs)) deallocate(derivs)
      allocate(derivs(3,nbloc))
      derivs(:,:)=0.d0
!
      call mechanic_up_para_respa(istep,.true.)
      call allocsteprespa(.true.)
!
!     get the fast-evolving potential energy and atomic forces
!
      call gradfast (ealt,derivs)
!
!       communicate forces
!
      call commforcesrespa(derivs,.true.)
!
!       MPI : get total energy
!
      call reduceen(ealt)
!
!     use Newton's second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the Verlet recursion
!
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            aalt(:,iglob) = -convert *&
            &derivs(:,i) / mass(iglob)
            v(:,iglob) = v(:,iglob) + aalt(:,iglob)*dta_2
         end if
      end do
      deallocate(derivs)
!
      if (use_rattle)  call rattle2 (dta_2)
!
!     Debug print information
!
      if (deb_Energy)call info_energy(rank)
      if (deb_Force) call info_forces(cBond)
      if (deb_Atom)  call info_minmax_pva
!
!     increment average virial from fast-evolving potential terms
      viralt(:,:)=viralt(:,:)+vir(:,:)/dshort

      if(ir .and. stepfast<nalt) then
         !call rotpole
         call compute_dipole(dip(:,stepfast)&
         &,dipind(:,stepfast),.FALSE.)
      endif
   end do

   if(use_piston) call apply_a_piston(dt_2,-1,.FALSE.)
!
!     Reassign the particules that have changed of domain
!
!     -> real space
!
   time0 = mpi_wtime()
   call reassignrespa(nalt,nalt)
!
!     -> reciprocal space
   call reassignpme(.false.)
!
!     communicate positions
!
   call commposrespa(.false.)
   call commposrec
   time1 = mpi_wtime()
   timecommpos = timecommpos + time1 - time0
!
!
   time0 = mpi_wtime()
   call reinitnl(istep)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
   time0 = mpi_wtime()
!
   call mechanic_up_para_respa(istep,.false.)

   call allocsteprespa(.false.)
   time1 = mpi_wtime()
   timeinte = timeinte + time1 - time0
!
!     rebuild the neighbor lists
!
   time0 = mpi_wtime()
   if (use_list) call nblist(istep)
   time1 = mpi_wtime()
   timenl = timenl + time1 - time0
!
   time0 = mpi_wtime()
   if(allocated(derivs)) deallocate(derivs)
   allocate(derivs(3,nbloc))
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
!     communicate forces
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
!     use Newton's second law to get the next accelerations;
!     find the full-step velocities using the BAOAB recursion
!
   time0 = mpi_wtime()
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         a(:,iglob) = -convert * derivs(:,i)/mass(iglob)
         v(:,iglob) = v(:,iglob) + dt_2*a(:,iglob)
      end if
   end do
!
!     find the constraint-corrected full-step velocities
!
   if (use_rattle)  call rattle2 (dt)
!
!     Debug print information
!
   if (deb_Energy) call info_energy(rank)
   if (deb_Force)  call info_forces(cNBond)
   if (deb_Atom)   call info_minmax_pva
!
!     total potential and virial from sum of fast and slow parts
!
   epot = epot + ealt
   vir = vir + viralt

   time1 = mpi_wtime()
   timeinte = timeinte + time1 - time0
!
   time0 = mpi_wtime()
   call temper (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
   call pressure2 (epot,temp)
   if(ir) then
      call compute_dipole(dip(:,nalt)&
      &,dipind(:,nalt),full_dipole)
      call save_dipole_respa(dip,dipind)
   endif
   if(use_piston) then
      call apply_b_piston(dt,pres,stress)
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
   call mdsave (istep,dt,epot,derivs)
   call mdrest (istep)
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
   return
end
