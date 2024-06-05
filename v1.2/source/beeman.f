c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine beeman  --  Beeman molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "beeman" performs a single molecular dynamics time step
c     via the Beeman multistep recursion formula; uses original
c     coefficients or Bernie Brooks' "Better Beeman" values
c
c     literature references:
c
c     D. Beeman, "Some Multistep Methods for Use in Molecular
c     Dynamics Calculations", Journal of Computational Physics,
c     20, 130-139 (1976)
c
c     B. R. Brooks, "Algorithms for Molecular Dynamics at Constant
c     Temperature and Pressure", DCRT Report, NIH, April 1988
c
c
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
c
      if (deb_Path) write(iout,*), 'beeman '
c
c
c     set time values and coefficients for Beeman integration
c
      factor = dble(bmnmix)
      dt_x = dt / factor
      part1 = 0.5d0*factor + 1.0d0
      part2 = part1 - 2.0d0
c
c     store the current atom positions, then find half-step
c     velocities and full-step positions via Beeman recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = v(j,iglob) + (part1*a(j,iglob)-
     $           aalt(j,iglob))*dt_x
            end do
            xold(iglob) = x(iglob)
            yold(iglob) = y(iglob)
            zold(iglob) = z(iglob)
            x(iglob) = x(iglob) + v(1,iglob)*dt
            y(iglob) = y(iglob) + v(2,iglob)*dt
            z(iglob) = z(iglob) + v(3,iglob)*dt
            end if
      end do
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0 
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      call reassign
c
c     -> reciprocal space
c
      time0 = mpi_wtime()
      call reassignpme(.false.)
      time1 = mpi_wtime()
      timereneig = timereneig + time1 - time0
c
c     communicate positions
c
      time0 = mpi_wtime()
      call commpos
      call commposrec
      time1 = mpi_wtime()
      timecommpos = timecommpos + time1 - time0
c
      time0 = mpi_wtime()
      allocate (derivs(3,nbloc))
      derivs = 0d0
c
      call reinitnl(istep)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
c
      time0 = mpi_wtime()
      call mechanic_up_para(istep)
      time1 = mpi_wtime()
      timeparam = timeparam + time1 - time0
c
      time0 = mpi_wtime()
      call allocstep
      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
c
c     rebuild the neighbor lists
c
      time0 = mpi_wtime()
      if (use_list) call nblist(istep)
      time1 = mpi_wtime()
      timenl = timenl + time1 - time0
c
c     get the potential energy and atomic forces
c
      time0 = mpi_wtime()
      call gradient (epot,derivs)
      time1 = mpi_wtime()
      timegrad = timegrad + time1 - time0
c
c     MPI : get total energy
c
      time0 = mpi_wtime()
      call reduceen(epot)
      time1 = mpi_wtime()
      timered = timered + time1 - time0
c
c     Debug print information
c
      if (deb_Energy) call info_energy(rank)
      if (deb_Force)  call info_forces(cDef)
      if (deb_Atom)   call info_minmax_pva
c
c     make half-step temperature and pressure corrections
c
      time0 = mpi_wtime()
      call temper2 (temp)
      call pressure2 (epot,temp)
      time1 = mpi_wtime()
      timetp = timetp + time1-time0
c
c     communicate forces
c
      time0 = mpi_wtime()
      call commforces(derivs)
      time1 = mpi_wtime()
      timecommforces = timecommforces + time1-time0
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Beeman recursion
c
      time0 = mpi_wtime()
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               aalt(j,iglob) = a(j,iglob)
               a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
               v(j,iglob) = v(j,iglob) +
     $          (part2*a(j,iglob)+aalt(j,iglob))*dt_x
            end do
         end if
      end do
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
c
c     make full-step temperature and pressure corrections
c
      time0 = mpi_wtime()
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress,istep)
      time1 = mpi_wtime()
      timetp = timetp + time1-time0
c
c     total energy is sum of kinetic and potential energies
c
      time0 = mpi_wtime()
      etot = eksum + esum
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,derivs)
      call mdrest (istep)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
