c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine bbk  --  BBK Langevin molecular dynamics step  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "bbk" performs a single molecular dynamics time step
c     via the BBK recursion formula
c
c     literature reference:
c
c     A. Brünger, C. L. Brooks III, M. Karplus, Stochastic boundary 
c     conditions fro molecular dynamics simulations of ST2 water. 
c     Chem. Phys. Letters, 1984, 105 (5) 495-500
c
c
      subroutine bbk (istep,dt)
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
      real*8 dt,dt_2
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 a1,a2,normal
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
      time0 = mpi_wtime()

      dt_2 = 0.5d0*dt
c
c     set time values and coefficients for BBK integration
c
      a1 = 1-dt_2*gamma
      a2 = 0.5d0*sqrt(2*boltzmann*kelvin*gamma*dt)
c
c     find half step velocities and full step positions via BBK recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = a1*v(j,iglob) 
     $       +  dt_2*a(j,iglob)
     $       +  a2*Rn(j,i)/sqrt(mass(iglob))
            end do
            xold(iglob) = x(iglob)
            yold(iglob) = y(iglob)
            zold(iglob) = z(iglob)
            x(iglob) = x(iglob) + v(1,iglob)*dt
            y(iglob) = y(iglob) + v(2,iglob)*dt
            z(iglob) = z(iglob) + v(3,iglob)*dt
            end if
      end do
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0 
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
c     get constraint-corrected positions and half-step velocities
c
      time0 = mpi_wtime()
      if (use_rattle)  call rattle (dt)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0 
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
      call mechanicstep(istep)
      time1 = mpi_wtime()
      timeparam = timeparam + time1 - time0
c
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
c     make half-step temperature and pressure corrections
c
      time0 = mpi_wtime()
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
c     compute random part
c
      time0 = mpi_wtime()
      deallocate (Rn)
      allocate (Rn(3,nloc))
      do i = 1, nloc
        do j = 1, 3
          Rn(j,i) = normal()
        end do
      end do
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BBK recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
               v(j,iglob) = (v(j,iglob) +
     $          + dt_2*a(j,iglob) + a2*Rn(j,i)/sqrt(mass(iglob)))/
     $          (1+dt_2*gamma)
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
      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
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
      call mdsave (istep,dt,epot)
      call mdrest (istep)
      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0
      return
      end
