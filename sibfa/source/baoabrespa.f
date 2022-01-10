c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine baoab  --  baoab Langevin molecular dynamics step  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "baoab" performs a single molecular dynamics time step
c     via the respa-baoab recursion formula
c
c     literature references:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c     D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-
c     Time-Step Molecular Dynamics Algorithm for Macromolecules",
c     Journal of Physical Chemistry, 98, 6885-6892 (1994)
c
c     X. Qian and T. Schlick, "Efficient Multiple-Time-Step Integrators
c     with Distance-Based Force Splitting for Particle-Mesh-Ewald
c     Molecular Dynamics Simulations", Journal of Chemical Physics,
c     115, 4019-4029 (2001)
c
c
      subroutine baoabrespa (istep,dt)
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
      use tortor
      use units
      use usage
      use virial
      use mpi
      implicit none
      integer i,j,k,istep,iglob
      integer nalt
      real*8 dt,factor
      real*8 dta,dta_2,dt_2,eps
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 part1,part2
      real*8 a1,a2,normal
      real*8 ealt,dalt
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 viralt(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
c
c     set some time values for the dynamics integration
c
      eps =  0.00000001d0
c      dalt = 0.00025d0
      dalt = 0.001d0
      nalt = int(dt/(dalt+eps)) + 1
      dalt = dble(nalt)
      dt_2 = 0.5d0 * dt
      dta = dt / dalt
      dta_2 = 0.5d0 * dta
c
c     set time values and coefficients for BAOAB integration
c
      a1 = exp(-gamma*dta)
      a2 = sqrt((1-a1**2)*boltzmann*kelvin)
c
c
c     find quarter step velocities and half step positions via baoab recursion
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
      if (use_rattle) call rattle2(dt_2)
c
c     respa inner loop
c
c     initialize virial from fast-evolving potential energy terms
c
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0d0
         end do
      end do
c
c     find fast-evolving velocities and positions via BAOAB recursion
c
      do k = 1, nalt
        do i = 1, nloc
           iglob = glob(i)
           if (use(iglob)) then
              do j = 1, 3
                 v(j,iglob) = v(j,iglob) + dta_2*aalt(j,iglob)
              end do
           end if
        end do
c
        if (use_rattle)  call rattle2 (dta_2)
c
        do i = 1, nloc
          iglob = glob(i)
          if (use(iglob)) then
            xold(iglob) = x(iglob)
            yold(iglob) = y(iglob)
            zold(iglob) = z(iglob)
            x(iglob) = x(iglob) + v(1,iglob)*dta_2
            y(iglob) = y(iglob) + v(2,iglob)*dta_2
            z(iglob) = z(iglob) + v(3,iglob)*dta_2
          end if
        end do
c
        if (use_rattle) call rattle(dta_2,xold,yold,zold)
        if (use_rattle) call rattle2(dta_2)
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
     $              a2*Rn(j,i)/sqrt(mass(iglob))
              end do
           end if
        end do
c
        if (use_rattle) call rattle2(dta)
c
        do i = 1, nloc
          iglob = glob(i)
          if (use(iglob)) then
            xold(iglob) = x(iglob)
            yold(iglob) = y(iglob)
            zold(iglob) = z(iglob)
            x(iglob) = x(iglob) + v(1,iglob)*dta_2
            y(iglob) = y(iglob) + v(2,iglob)*dta_2
            z(iglob) = z(iglob) + v(3,iglob)*dta_2
          end if
        end do
c
        if (use_rattle) call rattle(dta_2,xold,yold,zold)
        if (use_rattle) call rattle2(dta_2)
c
c       Reassign the particules that have changed of domain
c
c       -> real space
c
        time0 = mpi_wtime()
c
        call reassignrespa(.true.,k,nalt)
c
        time1 = mpi_wtime()
        timereneig = timereneig + time1 - time0
c
c       communicate positions
c
        time0 = mpi_wtime()
        call commposrespa(k.ne.nalt)
        time1 = mpi_wtime()
        timecommstep = timecommstep + time1 - time0
c
        allocate (derivs(3,nbloc))
        derivs = 0d0
c
        time0 = mpi_wtime()
        call mechanicsteprespa(istep,.true.)
        time1 = mpi_wtime()
        timeparam = timeparam + time1 - time0
        time0 = mpi_wtime()
        call allocsteprespa(.true.)
        time1 = mpi_wtime()
        timeclear = timeclear + time1 - time0
c
c     get the fast-evolving potential energy and atomic forces
c
        call gradfast (ealt,derivs)
c
c       communicate forces
c
        call commforcesrespa(derivs,.true.)
c
c       MPI : get total energy
c
        call reduceen(ealt)
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
          do i = 1, nloc
             iglob = glob(i)
             if (use(iglob)) then
                do j = 1, 3
                   aalt(j,iglob) = -convert *
     $                derivs(j,i) / mass(iglob)
                   v(j,iglob) = v(j,iglob) + aalt(j,iglob)*dta_2
                end do
             end if
          end do
        deallocate (derivs)
c
        if (use_rattle)  call rattle2 (dta_2)
c
c     increment average virial from fast-evolving potential terms
c
        do i = 1, 3
           do j = 1, 3
              viralt(j,i) = viralt(j,i) + vir(j,i)/dalt
           end do
        end do
      end do
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      time0 = mpi_wtime()
c
      call reassignrespa(.false.,nalt,nalt)
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
      call commposrespa(.false.)
      call commposrec
      time1 = mpi_wtime()
      timecommstep = timecommstep + time1 - time0
c
c
      call reinitnl(istep)
c
      time0 = mpi_wtime()
      call mechanicsteprespa(istep,.false.)
      time1 = mpi_wtime()
      timeparam = timeparam + time1 - time0

      time0 = mpi_wtime()
      call allocsteprespa(.false.)
      time1 = mpi_wtime()
      timeclear = timeclear + time1 - time0
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)
c
      allocate (derivs(3,nbloc))
      derivs = 0d0
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslow (epot,derivs)
c
c     if necessary, communicate some forces
c
      call commforcesrespa(derivs,.false.)
c
c     MPI : get total energy
c
      call reduceen(epot)
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
c     total potential and virial from sum of fast and slow parts
c
      epot = epot + ealt
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do
      end do
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
!      call mdrest (istep)
      return
      end
