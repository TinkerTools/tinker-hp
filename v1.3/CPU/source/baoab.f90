!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  subroutine baoab  --  BAOAB Langevin molecular dynamics step    ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     "baoab" performs a single molecular dynamics time step
!     via the BAOAB recursion formula
!
!     literature reference:
!
!     Efficient molecular dynamics using geodesic integration
!     and solvent-solute splitting B. Leimkuhler and C. Matthews,
!     Proceedings of the Royal Society A, 472: 20160138, 2016
!
!
!> @brief 
!> "baoab" performs a single molecular dynamics time step
!> via the BAOAB recursion formula
!> @param[in] istep: number of timestep
!> @param[in] dt: length of timestep
subroutine baoab (istep,dt)
   use atmtyp
   use atoms
   use bath
   use cutoff
   use domdec
   use energi
   use freeze
   use deriv
   use inform
   use iounit
   use mdstuf
   use moldyn
   use timestat
   use units
   use usage
   use mpi
   use sizes
   use spectra
   use utilbaoab
   use mdstuf1
   implicit none
   real*8, intent(in) :: dt
   integer, intent(in) :: istep
   integer i,iglob
   real*8 :: dip(3),dipind(3)
   real*8 :: dt_2
   real*8 :: time0,time1
   time0 = mpi_wtime()
!
   if (deb_Path) write(iout,*), 'baoab '
!
!
!     set time values and coefficients for BAOAB integration
!
   dt_2 = 0.5d0 * dt

   if (istep.eq.1) then
      if(use_piston) pres=atmsph
   end if

   if(use_piston) call apply_b_piston(dt_2,pres,stress)
!     find quarter step velocities and half step positions via BAOAB recursion
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) v(:,iglob) = v(:,iglob) + dt_2*a(:,iglob)
   end do
!
   if (use_rattle) then
      if (rank == 0) then
         write(0,*) "Error: RATTLE not compatible with ",&
         &"BAOAB integrator yet."
         call fatal
      endif
      call rattle2(dt)
   endif

   if(use_piston) then
      call apply_a_piston(dt_2,-1,.TRUE.)
      call apply_o_piston(dt)
   else
      call apply_a_block(dt_2)
   endif
   call apply_o_block(dt)
   if(use_piston) then
      call apply_a_piston(dt_2,-1,.TRUE.)
   else
      call apply_a_block(dt_2)
   endif
!
!
!     make half-step temperature and pressure corrections
!
   time0 = mpi_wtime()
   call pressure2 (epot,temp)
   time1 = mpi_wtime()
   timetp = timetp + time1 - time0
   time0=time1
!
!     Reassign the particules that have changed of domain
!
!     -> real space
   call reassign
!
!     -> reciprocal space
   call reassignpme(.false.)
!
!     communicate positions
!
   call commpos
   call commposrec
   time1 = mpi_wtime()
   timecommpos = timecommpos + time1 - time0
!
   time0 = mpi_wtime()
   if (allocated(derivs)) deallocate (derivs)
   allocate (derivs(3,nbloc))
   derivs(:,:) = 0.d0
!
   call reinitnl(istep)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
   time0 = mpi_wtime()
   call mechanic_up_para(istep)
   time1 = mpi_wtime()
   timeparam = timeparam + time1 - time0
!
   time0 = mpi_wtime()
   call allocstep
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
!     get the potential energy and atomic forces
!
   time0 = mpi_wtime()
   call gradient (epot,derivs)
   time1 = mpi_wtime()
   timegrad = timegrad + time1 - time0
!
!     MPI : get total energy
!
   time0 = mpi_wtime()
   call reduceen(epot)
   time1 = mpi_wtime()
   timered = timered + time1 - time0
!
!     communicate forces
!
   time0 = mpi_wtime()
   call commforces(derivs)
   time1 = mpi_wtime()
   timecommforces = timecommforces + time1-time0
!
!     Debug print information
!
   if (deb_Energy) call info_energy(rank)
   if (deb_Force)  call info_forces(cDef)
   if (deb_Atom)   call info_minmax_pva
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
   time1 = mpi_wtime()
   timeinte = timeinte + time1 - time0
!
   time0 = mpi_wtime()
   call temper (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
   if(ir) then
      call compute_dipole(dip,dipind,full_dipole)
      call save_dipole_traj(dip,dipind)
   endif

   if(use_piston) then
      call apply_b_piston(dt_2,pres,stress)
      call ddpme3dnpt(1.d0,istep)
   endif
!
!     total energy is sum of kinetic and potential energies
!
   etot = eksum + epot
   time1 = mpi_wtime()
   timetp = timetp + time1-time0
   time0=time1
!
!     compute statistics and save trajectory for this step
!
   call mdsave (istep,dt,epot,derivs)
   call mdrest (istep)
   call mdstat (istep,dt,etot,epot,eksum,temp,pres)

   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0

end
