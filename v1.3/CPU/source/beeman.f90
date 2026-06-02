!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine beeman  --  Beeman molecular dynamics step  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "beeman" performs a single molecular dynamics time step
!     via the Beeman multistep recursion formula; uses original
!     coefficients or Bernie Brooks' "Better Beeman" values
!
!     literature references:
!
!     D. Beeman, "Some Multistep Methods for Use in Molecular
!     Dynamics Calculations", Journal of Computational Physics,
!     20, 130-139 (1976)
!
!     B. R. Brooks, "Algorithms for Molecular Dynamics at Constant
!     Temperature and Pressure", DCRT Report, NIH, April 1988
!
!
!> @brief 
!> "beeman" performs a single molecular dynamics time step
!> via the Beeman multistep recursion formula; uses original
!> coefficients or Bernie Brooks' "Better Beeman" values
!> @param no params
subroutine beeman (istep,dt)
   use atmtyp
   use atoms
   use cutoff
   use deriv
   use domdec
   use energi
   use freeze
   use inform
   use iounit
   use mdstuf
   use moldyn
   use timestat
   use units
   use usage
   use mpi
   implicit none
   integer i,j,istep,iglob
   real*8 dt,dt_x,factor
   real*8 etot,eksum,epot
   real*8 temp,pres
   real*8 part1,part2
   real*8 ekin(3,3)
   real*8 stress(3,3)
   real*8 time0,time1
   real*8, allocatable :: derivs(:,:)
   time0 = mpi_wtime()
!
   if (deb_Path) write(iout,*), 'beeman '
!
!
!     set time values and coefficients for Beeman integration
!
   factor = dble(bmnmix)
   dt_x = dt / factor
   part1 = 0.5d0*factor + 1.0d0
   part2 = part1 - 2.0d0
!
!     store the current atom positions, then find half-step
!     velocities and full-step positions via Beeman recursion
!
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         do j = 1, 3
            v(j,iglob) = v(j,iglob) + (part1*a(j,iglob)-&
            &aalt(j,iglob))*dt_x
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
   call reassign
!
!     -> reciprocal space
!
   time0 = mpi_wtime()
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
!     Debug print information
!
   if (deb_Energy) call info_energy(rank)
   if (deb_Force)  call info_forces(cDef)
   if (deb_Atom)   call info_minmax_pva
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
!     use Newton's second law to get the next accelerations;
!     find the full-step velocities using the Beeman recursion
!
   time0 = mpi_wtime()
   do i = 1, nloc
      iglob = glob(i)
      if (use(iglob)) then
         do j = 1, 3
            aalt(j,iglob) = a(j,iglob)
            a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
            v(j,iglob) = v(j,iglob) +&
            &(part2*a(j,iglob)+aalt(j,iglob))*dt_x
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
   time1 = mpi_wtime()
   timetp = timetp + time1-time0
!
!     total energy is sum of kinetic and potential energies
!
   time0 = mpi_wtime()
   etot = eksum + esum
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
