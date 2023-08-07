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
      use deriv
      use inform
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
      integer i,j,iglob
      real*8 :: dip(3),dipind(3)
      real*8 :: dt_2,factor
      real*8 :: part1,part2
      real*8 :: a1,a2
      real*8 :: time0,time1
      time0 = mpi_wtime()
c
c     set time values and coefficients for BAOAB integration
c
      dt_2 = 0.5d0 * dt

      if (istep.eq.1) then
         if(use_piston) pres=atmsph
      end if

      if(use_piston) call apply_b_piston(dt_2,pres,stress)
c     find quarter step velocities and half step positions via BAOAB recursion
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) v(:,iglob) = v(:,iglob) + dt_2*a(:,iglob)
      end do
c
      if (use_rattle) then
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
c
c
c     make half-step temperature and pressure corrections
c
      time0 = mpi_wtime()
      call pressure2 (epot,temp)
      time1 = mpi_wtime()
      timetp = timetp + time1 - time0
      time0=time1
c
c     Reassign the particules that have changed of domain
c
c     -> real space
      call reassign
c
c     -> reciprocal space
      call reassignpme(.false.)
c
c     communicate positions
c    
      call commpos
      call commposrec
      time1 = mpi_wtime()
      timecommpos = timecommpos + time1 - time0
c
      time0 = mpi_wtime()
      if (allocated(derivs)) deallocate (derivs)
      allocate (derivs(3,nbloc))
      derivs(:,:) = 0.d0
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
c     communicate forces
c
      time0 = mpi_wtime()
      call commforces(derivs)
      time1 = mpi_wtime()
      timecommforces = timecommforces + time1-time0
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
      time0 = mpi_wtime()
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then            
          a(:,iglob) = -convert * derivs(:,i)/mass(iglob)
          v(:,iglob) = v(:,iglob) + dt_2*a(:,iglob)
        end if
      end do
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
      time1 = mpi_wtime()
      timeinte = timeinte + time1 - time0
c
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
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
      time1 = mpi_wtime()
      timetp = timetp + time1-time0
      time0=time1
c
c     compute statistics and save trajectory for this step
c
      call mdsave (istep,dt,epot,derivs)      
      call mdrest (istep)
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)

      time1 = mpi_wtime()
      timeinte = timeinte + time1-time0

      end
