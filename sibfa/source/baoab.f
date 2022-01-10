c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine baoab  --  BAOAB Langevin molecular dynamics step    ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "baoab" performs a single molecular dynamics time step
c     via the BAOAB recursion formula
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c
      subroutine baoab (istep,dt)
      use atmtyp
      use atoms
      use bath
      use cutoff
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use timestat
      use units
      use usage
      use mpi
      implicit none
      integer i,j,istep,iglob
      real*8 dt,dt_2,factor
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 part1,part2
      real*8 a1,a2,normal
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
c
c     store the current atom positions, then find half-step
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
c
c     set time values and coefficients for BAOAB integration
c
      a1 = exp(-gamma*dt)
      a2 = sqrt((1-a1**2)*boltzmann*kelvin)
      
      dt_2=0.5d0*dt
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = v(j,iglob) + dt_2*a(j,iglob)
            end do
         end if
      end do
c
      if (use_rattle) call rattle2(dt)
c
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
          xold(iglob) = x(iglob)
          yold(iglob) = y(iglob)
          zold(iglob) = z(iglob)
          x(iglob) = x(iglob) + v(1,iglob)*dt_2
          y(iglob) = y(iglob) + v(2,iglob)*dt_2
          z(iglob) = z(iglob) + v(3,iglob)*dt_2
        end if
      end do
c
      if (use_rattle) call rattle(dt_2,xold,yold,zold)
      if (use_rattle) call rattle2(dt_2)
c
      if (use_rattle) call rattle2(dt)
c
c     compute random part
c
      deallocate (Rn)
      allocate (Rn(3,nloc))
      do i = 1, nloc
        do j = 1, 3
          Rn(j,i) = normal()
        end do
      end do
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = a1*v(j,iglob) + 
     $            a2*Rn(j,i)/sqrt(mass(iglob))
            end do
         end if
      end do
c
      if (use_rattle) call rattle2(dt)
c
c     find full step positions via BAOAB recursion
c     
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            xold(iglob) = x(iglob)
            yold(iglob) = y(iglob)
            zold(iglob) = z(iglob)
            x(iglob) = x(iglob) + v(1,iglob)*dt_2
            y(iglob) = y(iglob) + v(2,iglob)*dt_2
            z(iglob) = z(iglob) + v(3,iglob)*dt_2
         end if
      end do
c
      if (use_rattle) call rattle(dt_2,xold,yold,zold)
      if (use_rattle) call rattle2(dt_2)
c
c
c     make half-step temperature and pressure corrections
c
c      call temper (dt,eksum,ekin,temp)
      call pressure2 (epot,temp)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      time0 = mpi_wtime()
c
      call reassign
c
c     -> reciprocal space
c
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
      timecommstep = timecommstep + time1 - time0
c
      allocate (derivs(3,nbloc))
      derivs = 0d0
c
      call reinitnl(istep)
c
      time0 = mpi_wtime()
      call mechanicstep(istep)
      time1 = mpi_wtime()
c
      timeparam = timeparam + time1 - time0
c
      time0 = mpi_wtime()
      call allocstep
      time1 = mpi_wtime()
      timeclear = timeclear  + time1 - time0
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     MPI : get total energy
c
      call reduceen(epot)
c
c     communicate forces
c
      call commforces(derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
               v(j,iglob) = v(j,iglob) + dt_2*a(j,iglob)
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + esum
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
!      call mdrest (istep)
      return
      end
