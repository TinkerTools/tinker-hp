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
!> @brief 
!> performs a single molecular dynamics time step
!> via the velocity Verlet multistep recursion formula
!> @param[in] istep: index of timestep
!> @param[in] dt: length of timestep
subroutine verlet (istep,dt)
   use atmtyp
   use atoms
   use cutoff
   use domdec
   use deriv
   use freeze
   use inform
   use iounit
   use moldyn
   use timestat
   use units
   use usage
   use mpi
   use spectra
   implicit none
   integer i,j,istep
   integer iglob
   real*8 dt,dt_2
   real*8 etot,epot
   real*8 eksum
   real*8 temp,pres
   real*8 ekin(3,3)
   real*8 stress(3,3)
   real*8 time0,time1
   real*8, allocatable :: derivs(:,:)
   real*8 dip(3),dipind(3)
   time0 = mpi_wtime()
!
   if (deb_Path) write(iout,*), 'verlet '
!
!
!     set some time values for the dynamics integration
!
   dt_2 = 0.5d0 * dt
!
!     store the current atom positions, then find half-step
!     velocities and full-step positions via Verlet recursion
!
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         do j = 1, 3
            v(j,iglob) = v(j,iglob) + a(j,iglob)*dt_2
         end do
         xold(iglob) = x(iglob)
         yold(iglob) = y(iglob)
         zold(iglob) = z(iglob)
         x(iglob) = x(iglob) + v(1,iglob)*dt
         y(iglob) = y(iglob) + v(2,iglob)*dt
         z(iglob) = z(iglob) + v(3,iglob)*dt
      end if
   end do
!
!     get constraint-corrected positions and half-step velocities
!
   if (use_rattle)  call rattle (dt)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
!     Reassign the particules that have changed of domain
!
!     -> real space
!
   time0 = mpi_wtime()
!
   call reassign
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
   call commpos
   call commposrec
   time1 = mpi_wtime()
   timecommpos = timecommpos + time1 - time0
!
   time0 = mpi_wtime()
   allocate (derivs(3,nbloc))
   derivs = 0d0
!
   call reinitnl(istep)
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
   time0 = mpi_wtime()
   call mechanic_up_para(istep)
   time1 = mpi_wtime()
   timeparam = timeparam + time1 - time0

   time0 = mpi_wtime()
   call allocstep
   time1 = mpi_wtime()
   timeinte = timeinte + time1 - time0
!
!     rebuild the neighbor list
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
!     make half-step temperature and pressure corrections
!
   time0 = mpi_wtime()
   call temper2 (temp)
   call pressure2 (epot,temp)
   time1 = mpi_wtime()
   timetp = timetp + time1-time0
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
!     find the full-step velocities using the Verlet recursion
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
   time1 = mpi_wtime()
   timeinte = timeinte + time1-time0
!
!     make full-step temperature and pressure corrections
!
   time0 = mpi_wtime()
   call temper (dt,eksum,ekin,temp)
   call pressure (dt,ekin,pres,stress,istep)
   if(ir) then
      call compute_dipole(dip,dipind,full_dipole)
      call save_dipole_traj(dip,dipind)
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
