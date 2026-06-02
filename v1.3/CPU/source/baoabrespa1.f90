!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################################
!     ##                                                                               ##
!     ##  subroutine baoabrespa1  --  baoab r-RESPA1 Langevin molecular dynamics step  ##
!     ##                                                                               ##
!     ###################################################################################
!
!
!     "baoabrespa1" performs a single multiple time step molecular dynamics
!     step using the reversible reference system propagation algorithm
!     (r-RESPA) via a BAOAB core with the potential split into fast-
!     intermediate and slow-evolving portions
!
!     literature references:
!
!     Pushing the Limits of Multiple-Time-Step Strategies
!     for Polarizable Point Dipole Molecular Dynamics
!     L Lagardere, F Aviat, JP Piquemal
!     The journal of physical chemistry letters 10, 2593-2599
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
!     Ruhong Zhou, Edward Harder, Huafeng Xu and B. J. Berne,
!     "Efficient multiple time step method for use with Ewald and
!     particle mesh Ewald for large biomolecular systems",
!     J. Chem. Phys. 115, 2348-2358 (2001)
!
!
!> @brief 
!> "baoabrespa1" performs a single PIMD time step
!> via the BAOAB-RESPA1 recursion formula
!> @param[in] istep: number of timestep
!> @param[in] dt: length of timestep
subroutine baoabrespa1(istep,dt)
   use atmtyp
   use atoms
   use cutoff
   use domdec
   use deriv
   use energi
   use freeze
   use inform
   use iounit
   use moldyn
   use mpi
   use timestat
   use units
   use usage
   use virial
   use utilbaoab
   use bath, only: use_piston
   use mdstuf1
   implicit none
   integer i,j,iglob
   integer istep
   real*8 dt,dt_2
   real*8 dta,dta_2,dta2
   real*8    time0,time1
!
   if (deb_Path) write(iout,*), 'baoabrespa1 '
!

   if(istep==1) then
      if(use_piston .and. rank==0) then
         write(0,*) "BAOABRESPA1 is not compatible with "&
         &//"LANGEVIN PISTON barostat"
         call fatal
      endif
   endif

   time0 = mpi_wtime()
!
!     set some time values for the dynamics integration
!
   dt_2  = 0.5d0 * dt
   dta   = dt / dinter
   dta_2 = 0.5d0 * dta
!
   dta2 = dta / dshort
!
!     store the current atom positions, then find half-step
!     velocities via BAOAB recursion
!
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         v(:,iglob) = v(:,iglob) + a(:,iglob)*dt_2
      end if
   end do
!
   if (use_rattle) call rattle2(dt_2)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
!     find intermediate-evolving velocities and positions via BAOAB recursion
!
   call baoabrespaint1(dta,dta2)
!
!     Reassign the particules that have changed of domain
!
!     -> reciprocal space
!
   time0 = mpi_wtime()
   call reassignpme(.false.)
!
!     communicate positions
!
   time0 = mpi_wtime()
   call commposrec
   time1 = mpi_wtime()
   timecommpos = timecommpos + time1 - time0
!
   time0 = mpi_wtime()
   call reinitnl(istep)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
   call mechanic_up_para_respa1(istep,2)

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
   if (allocated(derivs)) deallocate(derivs)
   allocate(derivs(3,nbloc))
   derivs = 0d0
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
!     get the slow-evolving potential energy and atomic forces
!
   time0 = mpi_wtime()
   call gradslow1 (epot,derivs)
   time1 = mpi_wtime()
   timegrad = timegrad + time1 - time0
!
!     communicate some forces
!
   time0 = mpi_wtime()
   call commforcesrespa1(derivs,2)
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
   if (deb_Energy) call info_energy(rank)
   if (deb_Force) call info_forces(cNBond)
   if (deb_Atom)   call info_minmax_pva
!
!     make half-step temperature and pressure corrections
!
   time0 = mpi_wtime()
   call temper2 (temp)
   time1 = mpi_wtime()
   timetp = timetp + time1-time0
!
!     use Newton's second law to get the slow accelerations;
!     find full-step velocities using BAOAB recursion
!
   time0 = mpi_wtime()
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         a(:,iglob) = -convert * derivs(:,i) / mass(iglob)
         v(:,iglob) = v(:,iglob) + a(:,iglob)*dt_2
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
   call pressure2 (epot,temp)
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

end
!
!     subroutine baoabrespaint1 :
!     find intermediate-evolving velocities and positions via BAOAB recursion
!
!> @brief 
!> "baoabrespaint1" performs a single time step of the intermediate forces
!> via the BAOAB-RESPA recursion formula
!> @param[in] dta: intermediate timestep
!> @param[in] dta2: short timestep
subroutine baoabrespaint1(dta,dta2)
   use atmtyp
   use cutoff
   use deriv
   use domdec
   use deriv
   use energi
   use freeze
   use inform
   use iounit
   use moldyn
   use potent
   use timestat
   use units
   use usage
   use virial
   use mpi
   use mdstuf1
   implicit none
   integer i,j,iglob,stepint
   real*8 dta,dta_2,dta2
   real*8 time0,time1
   time0 = mpi_wtime()
   dta_2 = 0.5d0 * dta
!
   if (deb_Path) write(iout,*), 'baoabrespaint1 '
!
!
!     initialize virial from fast-evolving potential energy terms
!
   viralt(:,:) = 0.d0
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0

   do stepint = 1, nalt
      time0=mpi_wtime()
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            v(:,iglob) = v(:,iglob) + aalt(:,iglob)*dta_2
         end if
      end do
!
      if (use_rattle)  call rattle2 (dta_2)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
!
!     find fast-evolving velocities and positions via BAOAB recursion
!
      call baoabrespafast1(dta2)
!
      time0 = mpi_wtime()
      call mechanic_up_para_respa1(-1,1)
      time1 = mpi_wtime()
      timeparam = timeparam + time1 - time0

      time0 = mpi_wtime()
      call allocsteprespa(.false.)
!
      if(allocated(derivs)) deallocate(derivs)
      allocate(derivs(3,nbloc))
      derivs(:,:)=0.d0

      if (allocated(desave)) deallocate(desave)
      allocate(desave(3,nbloc))
      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
!
!     get the fast-evolving potential energy and atomic forces
!
      time0 = mpi_wtime()
      call gradint1(ealt,derivs)
      time1 = mpi_wtime()
      timegrad = timegrad + time1 - time0
!
!     communicate forces
!
      time0 = mpi_wtime()
      call commforcesrespa1(derivs,1)
      time1 = mpi_wtime()
      timecommforces = timecommforces + time1-time0
!
!     MPI : get total energy
!
      time0 = mpi_wtime()
      call reduceen(ealt)
      time1 = mpi_wtime()
      timered = timered + time1 - time0
!
!     Debug print information
!
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  call info_forces(cSNBond)
      if (deb_Atom)   call info_minmax_pva
!
!     use Newton's second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the BAOAB recursion
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
      ealt = ealt + ealt2
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = viralt(j,i) + (viralt2(j,i)&
            &+ vir(j,i))/dinter
         end do
      end do
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
      if (use_lambdadyn.and.stepint.ne.nalt) then
         delambdaesave = 0d0
         delambdavsave = 0d0
      end if
   end do
end
!
!     subroutine baoabrespafast1 :
!     find fast-evolving velocities and positions via BAOAB recursion
!
!> @brief 
!> "baoabrespafast1" performs a single time step of the bonded forces
!> via the BAOAB-RESPA recursion formula
!> @param[in] dta: short timestep
subroutine baoabrespafast1(dta)
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
   use moldyn
   use timestat
   use units
   use usage
   use virial
   use mpi
   use utilbaoab
   use mdstuf1
   implicit none
   integer i,j,iglob,stepfast
   real*8 dta,dta_2
   real*8 time0,time1
!
   if (deb_Path) write(iout,*), 'baoabrespafast1 '
!

   time0 = mpi_wtime()
!
!     set time values and coefficients for BAOAB integration
!
   dta_2 = 0.5d0 * dta
!
!     initialize virial from fast-evolving potential energy terms
   viralt2(:,:)=0.d0
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0

   do stepfast = 1, nalt2
      time0=mpi_wtime()
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            v(:,iglob) = v(:,iglob) + aalt2(:,iglob)*dta_2
         end if
      end do
!
      if (use_rattle) then
         call rattle2 (dta_2)
      end if

      call apply_a_block(dta_2)
      call apply_o_block(dta)
      call apply_a_block(dta_2)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
!
!       Reassign the particules that have changed of domain
!
!       -> real space
      time0 = mpi_wtime()
      call reassignrespa(stepfast,nalt2)
      time1 = mpi_wtime()
      timereneig = timereneig + time1 - time0
!
!       communicate positions
!
      time0 = mpi_wtime()
      call commposrespa(stepfast.ne.nalt2)
      time1 = mpi_wtime()
      timecommpos = timecommpos + time1 - time0
!
      time0 = mpi_wtime()
      if (allocated(derivs)) deallocate(derivs)
      allocate(derivs(3,nbloc))
      derivs(:,:)=0.d0
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
!
      if (stepfast.eq.nalt2) then
         time0 = mpi_wtime()
         call mechanic_up_para_respa1(-1,0)
         timeparam = timeparam + time1 - time0
         time0 = mpi_wtime()
         call allocsteprespa(.true.)
         time1 = mpi_wtime()
         timeinte = timeinte + time1 - time0
      endif
!
!     get the fast-evolving potential energy and atomic forces
!
      time0 = mpi_wtime()
      call gradfast1(ealt2,derivs)
      time1 = mpi_wtime()
      timegrad = timegrad + time1 - time0
!
!       communicate forces
!
      time0 = mpi_wtime()
      call commforcesrespa1(derivs,0)
      time1 = mpi_wtime()
      timecommforces = timecommforces + time1-time0
!
!       MPI : get total energy
!
      time0 = mpi_wtime()
      call reduceen(ealt2)
      time1 = mpi_wtime()
      timered = timered + time1 - time0
!
!     Debug print information
!
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  call info_forces(cBond)
      if (deb_Atom)   call info_minmax_pva
!
!     use Newton's second law to get fast-evolving accelerations;
!     update fast-evolving velocities using the BAOAB recursion
!
      time0 = mpi_wtime()
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            aalt2(:,iglob) = -convert *&
            &derivs(:,i) / mass(iglob)
            v(:,iglob) = v(:,iglob) + aalt2(:,iglob)*dta_2
         end if
      end do
      deallocate(derivs)
!
      if (use_rattle)  call rattle2 (dta_2)
!
!     increment average virial from fast-evolving potential terms
!
      do i = 1, 3
         do j = 1, 3
            viralt2(j,i) = viralt2(j,i) + vir(j,i)/dshort
         end do
      end do
      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
   end do
end
