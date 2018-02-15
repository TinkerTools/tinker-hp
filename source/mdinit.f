c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mdinit  --  initialize a dynamics trajectory  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mdinit" initializes the velocities and accelerations
c     for a molecular dynamics trajectory, including restarts
c
c
      subroutine mdinit
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'mdstuf.i'
      include 'molcul.i'
      include 'moldyn.i'
      include 'mpole.i'
      include 'couple.i'
      include 'rgddyn.i'
      include 'units.i'
      include 'uprior.i'
      include 'usage.i'
      include 'openmp.i'
      integer i,j,k,idyn,iglob,nh
      integer next
      integer lext,freeunit
      integer ierr
      real*8 e
      real*8 maxwell,speed
      real*8 hmax,hmass
      real*8 sum,dmass
      real*8 vec(3)
      real*8, allocatable :: derivs(:,:)
      logical exist,heavy
      character*7 ext
      character*20 keyword
      character*120 dynfile
      character*120 record
      character*120 string
c
c     set default parameters for the dynamics trajectory
c
      integrate = 'BEEMAN'
      bmnmix = 8
      nfree = 0
      irest = 1
      heavy = .false.
      velsave = .false.
      frcsave = .false.
      uindsave = .false.
      use_pred = .true.
      polpred = 'ASPC'
      iprint = 100
c
c     set default values for temperature and pressure control
c
      thermostat = 'BUSSI'
      tautemp = 0.2d0
      collide = 0.1d0
      do i = 1, maxnose
         vnh(i) = 0.0d0
         qnh(i) = 0.0d0
         gnh(i) = 0.0d0
      end do
      barostat = 'BERENDSEN'
      anisotrop = .false.
      taupres = 2.0d0
      compress = 0.000046d0
      vbar = 0.0d0
      qbar = 0.0d0
      gbar = 0.0d0
      eta = 0.0d0
      voltrial = 20
      volmove = 100.0d0
      volscale = 'MOLECULAR'
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         else if (keyword(1:14) .eq. 'BEEMAN-MIXING ') then
            read (string,*,err=10,end=10)  bmnmix
         else if (keyword(1:16) .eq. 'DEGREES-FREEDOM ') then
            read (string,*,err=10,end=10)  nfree
         else if (keyword(1:15) .eq. 'REMOVE-INERTIA ') then
            read (string,*,err=10,end=10)  irest
         else if (keyword(1:15) .eq. 'HEAVY-HYDROGEN ') then
            heavy = .true.
         else if (keyword(1:14) .eq. 'SAVE-VELOCITY ') then
            velsave = .true.
         else if (keyword(1:11) .eq. 'SAVE-FORCE ') then
            frcsave = .true.
         else if (keyword(1:13) .eq. 'SAVE-INDUCED ') then
            uindsave = .true.
         else if (keyword(1:14) .eq. 'POLAR-PREDICT ') then
            use_pred = .true.
            call getword (record,polpred,next)
            call upcase (polpred)
         else if (keyword(1:11) .eq. 'THERMOSTAT ') then
            call getword (record,thermostat,next)
            call upcase (thermostat)
         else if (keyword(1:16) .eq. 'TAU-TEMPERATURE ') then
            read (string,*,err=10,end=10)  tautemp
         else if (keyword(1:10) .eq. 'COLLISION ') then
            read (string,*,err=10,end=10)  collide
         else if (keyword(1:9) .eq. 'BAROSTAT ') then
            call getword (record,barostat,next)
            call upcase (barostat)
         else if (keyword(1:15) .eq. 'ANISO-PRESSURE ') then
            anisotrop = .true.
         else if (keyword(1:13) .eq. 'TAU-PRESSURE ') then
            read (string,*,err=10,end=10)  taupres
         else if (keyword(1:9) .eq. 'COMPRESS ') then
            read (string,*,err=10,end=10)  compress
         else if (keyword(1:13) .eq. 'VOLUME-TRIAL ') then
            read (string,*,err=10,end=10)  voltrial
         else if (keyword(1:12) .eq. 'VOLUME-MOVE ') then
            read (string,*,err=10,end=10)  volmove
         else if (keyword(1:13) .eq. 'VOLUME-SCALE ') then
            call getword (record,volscale,next)
            call upcase (volscale)
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=10,end=10)  iprint
         end if
   10    continue
      end do
c
c     repartition hydrogen masses to allow long time steps
c
      if (heavy) then
         if (hostrank.ne.0) goto 11
         hmax = 4.0d0
         do i = 1, n
            nh = 0
            sum = mass(i)
            do j = 1, n12(i)
               k = i12(j,i)
               if (atomic(k) .eq. 1) then
                  nh = nh + 1
                  sum = sum + mass(k)
               end if
            end do
            hmass = max(hmax,sum/dble(nh+1))
            do j = 1, n12(i)
               k = i12(j,i)
               if (atomic(k) .eq. 1) then
                  dmass = hmass - mass(k)
                  mass(k) = mass(k) + dmass
                  mass(i) = mass(i) - dmass
               end if
            end do
         end do
 11      call MPI_BARRIER(hostcomm,ierr)
      end if
c
c     make sure all atoms or groups have a nonzero mass
c
      do i = 1, n
         if (use(i) .and. mass(i).le.0.0d0) then
            mass(i) = 1.0d0
            totmass = totmass + 1.0d0
            write (iout,30)  i
   30       format (/,' MDINIT  --  Warning, Mass of Atom',i6,
     &                 ' Set to 1.0 for Dynamics')
         end if
      end do
c
c     perform dynamic allocation of some pointer arrays
c
      if (use_pred) then
         if (associated(udalt))  deallocate (udalt)
         if (associated(upalt))  deallocate (upalt)
         allocate (udalt(maxualt,3,n))
         allocate (upalt(maxualt,3,n))
c
c     set the Gear predictor binomial coefficients
c
         gear(1) = 6.0d0
         gear(2) = -15.0d0
         gear(3) = 20.0d0
         gear(4) = -15.0d0
         gear(5) = 6.0d0
         gear(6) = -1.0d0
         gear(7) = 0.0d0
c
c     set always stable predictor-corrector (ASPC) coefficients
c
         aspc(1) = 22.0d0 / 7.0d0
         aspc(2) = -55.0d0 / 14.0d0
         aspc(3) = 55.0d0 / 21.0d0
         aspc(4) = -22.0d0 / 21.0d0
         aspc(5) = 5.0d0 / 21.0d0
         aspc(6) = -1.0d0 / 42.0d0
         aspc(7) = 0.0d0
c
c    initialize prior values of induced dipole moments
c
         nualt = 0
         do i = 1, npole
            do j = 1, 3
               do k = 1, maxualt
                  udalt(k,j,i) = 0.0d0
                  upalt(k,j,i) = 0.0d0
               end do
            end do
         end do
      end if
c
c     enforce use of velocity Verlet with Andersen thermostat
c
      if (thermostat .eq. 'ANDERSEN') then
         if (integrate .eq. 'BEEMAN')  integrate = 'VERLET'
      end if
c
c     enforce use of Bussi thermostat/barostat with integrator
c
      if (integrate .eq. 'BUSSI') then
         thermostat = 'BUSSI'
         barostat = 'BUSSI'
      else if (thermostat.eq.'BUSSI' .and. barostat.eq.'BUSSI') then
         integrate = 'BUSSI'
      end if
c
c     set the number of degrees of freedom for the system
c
      if (nfree .eq. 0) then
         nfree = 3 * nuse
         if (isothermal .and. thermostat.ne.'ANDERSEN'.and.
     &       integrate.ne.'STOCHASTIC' .and. integrate.ne.'GHMC') then
            if (use_bounds) then
               nfree = nfree - 3
            else
               nfree = nfree - 6
            end if
         end if
         if (barostat .eq. 'BUSSI')  nfree = nfree + 1
      end if
c
c     check for a nonzero number of degrees of freedom
c
      if (nfree .lt. 0)  nfree = 0
      if (debug) then
         write (iout,50)  nfree
   50    format (/,' Number of Degrees of Freedom for Dynamics :',i10)
      end if
      if (nfree .eq. 0) then
         write (iout,60)
   60    format (/,' MDINIT  --  No Degrees of Freedom for Dynamics')
         call fatal
      end if
c
c     decide whether to remove center of mass motion
c
      dorest = .true.
      if (irest .eq. 0)  dorest = .false.
      if (nuse. ne. n)  dorest = .false.
      if (isothermal .and. thermostat.eq.'ANDERSEN')  dorest = .false.
c
c     try to restart using prior velocities and accelerations
c
      dynfile = filename(1:leng)//'.dyn'
      call version (dynfile,'old')
      inquire (file=dynfile,exist=exist)
      if (exist) then
         idyn = freeunit ()
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
         call readdyn (idyn)
         close (unit=idyn)
c
c     Do the domain decomposition
c
         call ddpme3d
         call reinitnl(0)
         call mechanicstep(0)
         call nblist(0)
c
c     set velocities and fast/slow accelerations for RESPA method
c
      else if (integrate .eq. 'RESPA') then
         allocate (derivs(3,nbloc))
         derivs = 0d0
         call allocsteprespa(.false.)
         call gradslow (e,derivs)
         call commforcesrespa(derivs,.false.)
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
               speed = maxwell (mass(iglob),kelvin)
               call ranvec (vec)
               do j = 1, 3
                  v(j,iglob) = speed * vec(j)
                  a(j,iglob) = -convert * derivs(j,i) / mass(iglob)
               end do
            else
               do j = 1, 3
                  v(j,i) = 0.0d0
                  a(j,i) = 0.0d0
                  aalt(j,i) = 0.0d0
               end do
            end if
         end do
         derivs = 0d0
         call allocsteprespa(.true.)
         call gradfast (e,derivs)
         call commforcesrespa(derivs,.true.)
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
               do j = 1, 3
                  aalt(j,iglob) = -convert *
     $              derivs(j,i) / mass(iglob)
               end do
            else
               do j = 1, 3
                  v(j,iglob) = 0.0d0
                  a(j,iglob) = 0.0d0
                  aalt(j,iglob) = 0.0d0
               end do
            end if
         end do
         deallocate (derivs)
         if (nuse .eq. n)  call mdrest (0)
c
c     set velocities and accelerations for Cartesian dynamics
c
      else
         allocate (derivs(3,nbloc))
         derivs = 0d0
         call allocstep
         call gradient (e,derivs)
         call commforces(derivs)
c
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
               speed = maxwell (mass(iglob),kelvin)
               call ranvec (vec)
               do j = 1, 3
                  v(j,iglob) = speed * vec(j)
                  a(j,iglob) = -convert * derivs(j,i) / mass(iglob)
                  aalt(j,iglob) = a(j,iglob)
               end do
            else
               do j = 1, 3
                  v(j,iglob) = 0.0d0
                  a(j,iglob) = 0.0d0
                  aalt(j,iglob) = 0.0d0
               end do
            end if
         end do
         deallocate (derivs)
         if (nuse .eq. n)  call mdrest (0)
      end if
c
c     check for any prior dynamics coordinate sets
c
      i = 0
      exist = .true.
      do while (exist)
         i = i + 1
         lext = 3
         call numeral (i,ext,lext)
         dynfile = filename(1:leng)//'.'//ext(1:lext)
         inquire (file=dynfile,exist=exist)
         if (.not.exist .and. i.lt.100) then
            lext = 2
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
         if (.not.exist .and. i.lt.10) then
            lext = 1
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
      end do
      nprior = i - 1
      return
      end
c
