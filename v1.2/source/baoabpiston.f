c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################################################
c     ##                                                                                           ##
c     ##  subroutine baoabpiston  --  BAOAB/Langein Piston NPT Langevin molecular dynamics step    ##
c     ##                                                                                           ##
c     ###############################################################################################
c
c
c     "baoabpiston" performs a single molecular dynamics time step
c     via the BAOAB recursion formula and using a Langevin Piston barostat
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c     
c      Constant pressure molecular dynamics simulation: The Langevin piston method
c      J. Chem. Phys. 103, 4613 (1995)
c      Scott E. Feller, Yuhong Zhang,Richard W. Pastor, Bernie Brooks
c
      subroutine baoabpiston (istep,dt)
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
      integer i,j,istep,iglob,ierr
      real*8 dt,dt_2
      real*8 etot,eksum,epot
      real*8 temp,pres

      real*8 a1,a2,normal
      real*8 a1piston,a2piston,R
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
c
c     set time values and coefficients for BAOAB integration
c
      a1 = exp(-gamma*dt)
      a2 = sqrt((1-a1**2)*boltzmann*kelvin)
      a1piston = exp(-gammapiston*dt)
      a2piston = sqrt((1-a1piston**2)*boltzmann*kelvin)
      dt_2 = 0.5d0*dt
      
c
c     find quarter step velocities and half step positions via BAOAB recursion
c
      vextvol = vextvol + dt_2*aextvol
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
      extvolold = extvol
      extvol = extvol + vextvol*dt_2
      call rescale(istep)
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
      if (use_rattle) call rattle(dt_2)
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
      if (rank.eq.0) then
        R = normal()
        vextvol = a1piston*vextvol + a2piston*R/sqrt(masspiston)
      end if
      call MPI_BCAST(vextvol,1,MPI_REAL8,0,COMM_TINKER,ierr)
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
      if (use_rattle) call rattle(dt_2)
      if (use_rattle) call rattle2(dt_2)

      extvolold = extvol
      extvol = extvol + vextvol*dt_2
      call rescale(istep)
c
c
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

      call temper (dt,eksum,ekin,temp)
c      call pressure2 (epot,temp)

      call pressure (dt,ekin,pres,stress,istep)
c
      aextvol = convert*(pres-atmsph)/(prescon*masspiston)
      vextvol = vextvol + dt_2*aextvol
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + esum
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
      call mdrest (istep)
      return
      end
