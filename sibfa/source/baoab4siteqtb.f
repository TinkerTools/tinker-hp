c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #########################################################################################
c     ##                                                                                     ##
c     ##  subroutine baoab4siteqtb  --  Baoab 4-site water MD step with Quantum Thermal Bath ##
c     ##                                                                                     ##
c     #########################################################################################
c
c     literature references:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c     this is a special version for rigid 4-site water models such as
c     TIP4P; the 4th site is expressed as a linear combination of the
c     other three to avoid a degeneracy in the RATTLE constraints
c
c     literature reference:
c
c     G. Ciccotti, M. Ferrario and J. P. Ryckaert, "Molecular Dynamics
c     of Rigid Systems in Cartesian Coordinates. A General Formulation",
c     Molecular Physics, 47, 1253-1264 (1982)
c
c
      subroutine baoab4siteqtb(istep,dt)
      use adqtb
      use atmtyp
      use atoms
      use bath
      use bound
      use cutoff
      use deriv
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use qtb
      use timestat
      use units
      use usage
      use neigh
      use mpi
      implicit none
      integer i,j,istep,nm,iglob,iloc
      integer iloc1,iloc2,iloc3,ilocrec1,ilocrec2,ilocrec3
      real*8 a1,a2,normal
      real*8 dt,dt_2,factor
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 part1,part2
      real*8 oterm,hterm
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
c
      dt_2 = 0.5d0*dt
c
c     set linear combination factors for TIP4P or Dang-Chang water
c
      oterm = 0.73612d0
      hterm = 0.13194d0
c      oterm = 0.7439758702395579d0
c      hterm = 0.1280120648475456d0
c      do i = 1, n
c         if (story(1)(1:10) .eq. 'DANG-CHANG') then
c            oterm = 0.6330320807091532d0
c            hterm = 0.1834839596454234d0
c         end if
c      end do
cc
cc     set time values and coefficients for BAOAB integration
cc
c      a1 = exp(-gamma*dt)
c      a2 = sqrt((1-a1**2)*boltzmann*kelvin)
c
c
      do i = 1, nloc
         iglob = glob(i)
         if (atomic(iglob).eq.0) then
           v(:,iglob) = 0d0
         else if (use(iglob)) then
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
        if (atomic(iglob).eq.0) then
           v(:,iglob) = 0d0
           xold(iglob) = x(iglob)
           yold(iglob) = y(iglob)
           zold(iglob) = z(iglob)
        else if (use(iglob)) then
          xold(iglob) = x(iglob)
          yold(iglob) = y(iglob)
          zold(iglob) = z(iglob)
          x(iglob) = x(iglob) + v(1,iglob)*dt_2
          y(iglob) = y(iglob) + v(2,iglob)*dt_2
          z(iglob) = z(iglob) + v(3,iglob)*dt_2
        end if
      end do

      if (use_rattle) call rattle(dt_2,xold,yold,zold)
      if (use_rattle) call rattle2(dt_2)
c
      if (use_rattle) call rattle2(dt)
c
c     compute random part
c
      call qtbrandom(dt,istep)

      if (use_rattle) call rattle2(dt)

      do i = 1, nloc
         iglob = glob(i)
         if (atomic(iglob).eq.0) then
           v(:,iglob) = 0d0
           xold(iglob) = x(iglob)
           yold(iglob) = y(iglob)
           zold(iglob) = z(iglob)
         else if (use(iglob)) then
            xold(iglob) = x(iglob)
            yold(iglob) = y(iglob)
            zold(iglob) = z(iglob)
            x(iglob) = x(iglob) + v(1,iglob)*dt_2
            y(iglob) = y(iglob) + v(2,iglob)*dt_2
            z(iglob) = z(iglob) + v(3,iglob)*dt_2
         end if
      end do

      if (use_rattle) call rattle(dt_2,xold,yold,zold)
      if (use_rattle) call rattle2(dt_2)
c
c     make half-step pressure corrections
c
      call pressure2 (epot,temp)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      time0 = mpi_wtime()
c
c      call reassign4site
c      call reassign
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
c
c     move M center to linear combination of H-O-H coordinates
c
      do i = 1, nlocnl
         iglob = ineignl(i)
c      do i = 1, nbloc
c         iglob = glob(i)
         if (atomic(iglob) .eq. 0) then
            x(iglob) = oterm*x(iglob-3) + hterm*x(iglob-2) + 
     $           hterm*x(iglob-1)
            y(iglob) = oterm*y(iglob-3) + hterm*y(iglob-2) +
     $           hterm*y(iglob-1)
            z(iglob) = oterm*z(iglob-3) + hterm*z(iglob-2) + 
     $           hterm*z(iglob-1)
         end if
      end do
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
c     transfer the forces on the M center to the H-O-H sites
c
      do i = 1, nlocnl
         iglob = ineignl(i)
         iloc = loc(iglob)
         if (atomic(iglob) .eq. 0) then
            iloc1 = loc(iglob-1)
            iloc2 = loc(iglob-2)
            iloc3 = loc(iglob-3)
            if ((iloc.gt.0).and.(iloc.le.nbloc)) then
            if ((iloc3.gt.0).and.(iloc3.le.nbloc)) then
            do j = 1, 3
               derivs(j,iloc3) = derivs(j,iloc3) + oterm*derivs(j,iloc)
            end do
            end if
            if ((iloc2.gt.0).and.(iloc2.le.nbloc)) then
            do j = 1, 3
               derivs(j,iloc2) = derivs(j,iloc2) + hterm*derivs(j,iloc)
            end do
            end if
            if ((iloc1.gt.0).and.(iloc1.le.nbloc)) then
            do j = 1, 3
               derivs(j,iloc1) = derivs(j,iloc1) + hterm*derivs(j,iloc)
            end do
            end if
            end if
         end if
      end do
      do i = 1, nlocrec
         iglob = globrec(i)
         iloc = locrec1(iglob)
         if (atomic(iglob) .eq. 0) then
            iloc1 = locrec1(iglob-1)
            iloc2 = locrec1(iglob-2)
            iloc3 = locrec1(iglob-3)
            if (iloc.ne.0) then
            if (iloc3.ne.0) then
            do j = 1, 3
               decrec(j,iloc3) = decrec(j,iloc3) + oterm*decrec(j,iloc)
            end do
            end if
            if (iloc2.ne.0) then
            do j = 1, 3
               decrec(j,iloc2) = decrec(j,iloc2) + hterm*decrec(j,iloc)
            end do
            end if
            if (iloc1.ne.0) then
            do j = 1, 3
               decrec(j,iloc1) = decrec(j,iloc1) + hterm*decrec(j,iloc)
            end do
            end if
            end if
         end if
      end do
c
c     MPI : get total energy
c
      call reduceen(epot)

c      call temper2 (dt,temp)
c      call pressure2 (epot,temp)
c
c     communicate forces
c
      call commforces(derivs)
cc
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Beeman recursion
c
      do i = 1, nloc
         iglob = glob(i)
c         write(*,*) 'i = ',iglob,'deriv = ',derivs(1,i)
         if (atomic(iglob) .eq. 0) then
            do j = 1, 3
               v(j,iglob) = 0.0d0
               a(j,iglob) = 0.0d0
            end do
         else if (use(iglob)) then
c           write(*,*) 'i = ',iglob,'deriv = ',derivs(1,i)
            do j = 1, 3
               aalt(j,iglob) = a(j,iglob)
               a(j,iglob) = -convert * derivs(j,i) / mass(iglob)
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
c     make full-step temperature and pressure corrections
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
c      call mdrest (istep)
      if ((mod(istep,nseg).eq.0).and.(adaptive)) then
            if (ir) then
              call irspectra
            endif
            compteur=compteur+1
            write(*,*) compteur
            if(compteur .ge. skipseg) then
            call fftw(istep,dt)
            call adHnoise(dt)
            endif
      endif
      if ((mod(istep,nseg).eq.0)) call convolseg
      return
      end
