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
      subroutine baoabpistonqtb (istep,dt)
      use adqtb
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
      use qtb
      implicit none
      integer i,j,istep,iglob,ierr
      real*8 dt,dt_x,factor
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 part1,part2
      real*8 a1,a2,normal
      real*8 a1piston,a2piston,R
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
      real*8 dt_2,scale,velfact
c
c     set time values and coefficients for BAOAB integration
c
      dt_2=0.5d0*dt
c      a1 = exp(-gamma*dt)
c      a2 = sqrt((1.d0-a1**2)*boltzmann*kelvin)
      a1piston = exp(-gammapiston*dt)
      a2piston = sqrt((1.d0-a1piston**2)*boltzmann*kelvin/masspiston)
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
      scale=exp(dt_2*vextvol) 
      velfact=sinh(dt_2*vextvol)/vextvol
      extvol = extvol*(scale**3) ! + vextvol*0.5*dt
      call rescale(istep)
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
          xold(iglob) = x(iglob)
          yold(iglob) = y(iglob)
          zold(iglob) = z(iglob)          
          x(iglob) = x(iglob)*scale + v(1,iglob)*velfact
          y(iglob) = y(iglob)*scale + v(2,iglob)*velfact
          z(iglob) = z(iglob)*scale + v(3,iglob)*velfact
          v(:,iglob)=v(:,iglob)/scale          
        end if
      end do
c
      if (use_rattle) call rattle(dt_2,xold,yold,zold)
      if (use_rattle) call rattle2(dt_2)
c
c      if (use_rattle) call rattle2(dt)
c
c     compute random part
c
!      deallocate (Rn)
!      allocate (Rn(3,nloc))
!      do i = 1, nloc
!        do j = 1, 3
!          Rn(j,i) = normal()
!        end do
!      end do
!      do i = 1, nloc
!         iglob = glob(i)
!         if (use(iglob)) then
!            do j = 1, 3
!               v(j,iglob) = a1*v(j,iglob) + 
!     $            a2*Rn(j,i)/sqrt(mass(iglob))
!            end do
!         end if
!      end do

      call qtbrandom(dt,istep)

      if (rank.eq.0) then
        R = normal()
        vextvol = a1piston*vextvol + a2piston*R
      end if
      call MPI_BCAST(vextvol,1,MPI_REAL8,0,COMM_BEAD,ierr)
c
      if (use_rattle) call rattle2(dt)
c
c     find full step positions via BAOAB recursion
c
      extvolold = extvol
      scale=exp(dt_2*vextvol) 
      velfact=sinh(dt_2*vextvol)/vextvol
      extvol = extvol*(scale**3) ! + vextvol*0.5*dt
      call rescale(istep)
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
          xold(iglob) = x(iglob)
          yold(iglob) = y(iglob)
          zold(iglob) = z(iglob)          
          x(iglob) = x(iglob)*scale + v(1,iglob)*velfact
          y(iglob) = y(iglob)*scale + v(2,iglob)*velfact
          z(iglob) = z(iglob)*scale + v(3,iglob)*velfact
          v(:,iglob)=v(:,iglob)/scale          
        end if
      end do
c
      if (use_rattle) call rattle(0.5*dt,xold,yold,zold)
      if (use_rattle) call rattle2(0.5*dt)
c
c
c
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      time0 = mpi_wtime()
c
      call reassignqtb(istep)
      
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
      call pressure (dt,epot,ekin,temp,pres,stress,istep)

      aextvol = 3.d0*convert*(extvol*(pres-atmsph)/prescon 
     &  +gasconst*kelvin)/masspiston
      vextvol = vextvol + dt_2*aextvol

c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c      call temper (dt,eksum,ekin,temp)
c      call pressure (dt,epot,ekin,temp,pres,stress,istep)
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
      if ((mod(istep,nseg).eq.0).and.adaptive) then
                  compteur=compteur+1
            if(compteur .ge. skipseg) then
            call fftw(istep,dt)
            call adHnoise(dt)
            endif
      endif
      if ((mod(istep,nseg).eq.0)) call convolseg
      return
      end
