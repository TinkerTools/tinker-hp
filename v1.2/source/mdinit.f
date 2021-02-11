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
      subroutine mdinit(dt)
      use atmtyp
      use atoms
      use bath
      use bound
      use couple
      use cutoff
      use deriv
      use domdec
      use files
      use group
      use keys
      use freeze
      use inform
      use iounit
      use langevin
      use mdstuf
      use molcul
      use moldyn
      use mpole
      use neigh
      use polpot
      use units 
      use uprior
      use usage
      use virial
#ifdef PLUMED
      use plumed
      use mpi
#endif
      implicit none
      integer i,j,k,idyn,iglob
      integer next
      integer lext,freeunit
      integer ierr
      real*8 e
      real*8 maxwell,speed
      real*8 normal
      real*8 vec(3)
      real*8 dt
      real*8, allocatable :: derivs(:,:)
      real*8 eps
      logical exist
      character*7 ext
      character*20 keyword
      character*240 dynfile
      character*240 record
      character*240 string
c
c     set default parameters for the dynamics trajectory
c
      integrate = 'BEEMAN'
      bmnmix = 8
      nfree = 0
      irest = 1
      velsave = .false.
      frcsave = .false.
      uindsave = .false.
      use_pred = .true.
      polpred = 'ASPC'
      iprint = 100
c
c      set default for neighbor list update, reference = 20 for a 2fs time step
c
      ineigup = int(20*0.002d0/dt)
c
c     set default values for temperature and pressure control
c
      thermostat = 'BUSSI'
      tautemp = 0.2d0
      collide = 0.1d0
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
c      volscale = 'ATOMIC'
      gamma = 1.0d0
      gammapiston = 20.0d0
      masspiston = 0.000200d0
      virsave = 0d0
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         else if (keyword(1:14) .eq. 'BEEMAN-MIXING ') then
            read (string,*,err=10,end=10)  bmnmix
         else if (keyword(1:16) .eq. 'DEGREES-FREEDOM ') then
            read (string,*,err=10,end=10)  nfree
         else if (keyword(1:15) .eq. 'REMOVE-INERTIA ') then
            read (string,*,err=10,end=10)  irest
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
            if (polpred .eq. 'NONE') use_pred = .false.
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
         else if (keyword(1:14) .eq. 'NLUPDATE ') then
            read (string,*,err=10,end=10) ineigup
         else if (keyword(1:9) .eq. 'FRICTION ') then
            read (string,*,err=10,end=10) gamma
         else if (keyword(1:15) .eq. 'FRICTIONPISTON ') then
            read (string,*,err=10,end=10) gammapiston
         else if (keyword(1:12) .eq. 'MASSPISTON ') then
            read (string,*,err=10,end=10) masspiston
#ifdef PLUMED
c
c initialise PLUMED
c
         else if (keyword(1:7) .eq. 'PLUMED ') then
            lplumed = .true.
            allocate(pl_force(3,n))
            allocate(pl_pos  (3,n))
            allocate(pl_mass (  n))
            allocate(pl_glob (  n))
            call getword (record,pl_input,next)
            call getword (record,pl_output,next)
c            if (rank.eq.0) then
c              write(0,*)"  CREATING PLUMED FROM THE PROGRAM"
c            endif
c
c conversion factors to go from TINKER's units to PLUMED
c
c kcal/mol -> kJ/mol
            energyUnits = 4.184d0
c angstrom -> nm
            lengthUnits = 0.1d0
c fs -> ps
            timeUnits = 0.001d0

            call plumed_f_gcreate()
            call plumed_f_gcmd("setRealPrecision"//char(0),8)
            call plumed_f_gcmd("setMPIFComm"//char(0),COMM_TINKER)
            call plumed_f_gcmd("setMDEnergyUnits"//char(0),energyUnits)
            call plumed_f_gcmd("setMDLengthUnits"//char(0),lengthUnits)
            call plumed_f_gcmd("setMDTimeUnits"//char(0),timeUnits)
            call plumed_f_gcmd("setPlumedDat"//char(0),trim(pl_input)//
     $          char(0))
            call plumed_f_gcmd("setLogFile"//char(0),trim(pl_output)//
     $          char(0))
            call plumed_f_gcmd("setLogFile"//char(0),pl_output)
            call plumed_f_gcmd("setNatoms"//char(0),n)
            call plumed_f_gcmd("setMDEngine"//char(0),"TinkerHP");
            call plumed_f_gcmd("setTimestep"//char(0),dt);
            call plumed_f_gcmd("init"//char(0),0);
c
c end PLUMED initialisation
c
#endif
         end if
   10    continue
      end do
c
c     enforce the use of monte-carlo barostat with the TCG family of solvers
c
      if ((isobaric).and.(polalg.eq.3)) then
        barostat = 'MONTECARLO' 
        if (rank.eq.0) write(iout,*)'enforcing the use of Monte Carlo ',
     $    'Barostat with TCG'
      end if
c
c     enforce the use of baoabpiston integrator with Langevin Piston barostat
c
      if (barostat.eq.'LANGEVIN') then
        integrate = 'BAOABPISTON'
        call plangevin()
      end if
c
c     default time steps for respa and respa1 integrators
c
      if ((integrate.eq.'RESPA').or.(integrate.eq.'BAOABRESPA')) then
        eps =  0.00000001d0
        dshort = 0.00025
      else if ((integrate.eq.'RESPA1').or.(integrate.eq.'BAOABRESPA1'))
     $  then
        dinter = 0.002
        dshort = 0.00025
      end if
c
c     keywords for respa and respa1 integrators
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:7) .eq. 'DSHORT ') then
            read (string,*,err=20,end=20) dshort
         else if (keyword(1:7) .eq. 'DINTER ') then
            read (string,*,err=20,end=20) dinter
         end if
   20    continue
      end do
c
      if ((integrate.eq.'RESPA').or.(integrate.eq.'BAOABRESPA')) then
        eps =  0.00000001d0
        nalt = int(dt/(dshort+eps)) + 1
        dshort = dble(nalt)
      else if ((integrate.eq.'RESPA1').or.(integrate.eq.'BAOABRESPA1'))
     $  then
        eps =  0.00000001d0
        nalt = int(dt/(dinter+eps)) + 1
        nalt2 = int(dinter/(dshort+eps)) + 1
        dinter = dble(nalt)
        dshort = dble(nalt2)
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
c     enforce use of velocity Verlet with Andersen thermostat
c
      if (thermostat .eq. 'ANDERSEN') then
         if (integrate .eq. 'BEEMAN')  integrate = 'VERLET'
      end if
c
c     check for use of Monte Carlo barostat with constraints
c
      if (barostat.eq.'MONTECARLO' .and. volscale.eq.'ATOMIC') then
         if (use_rattle) then
            write (iout,40)
   40       format (/,' MDINIT  --  Atom-based Monte Carlo',
     &                 ' Barostat Incompatible with RATTLE')
            call fatal
         end if
      end if

c
c     nblist initialization for respa-n integrators
c
      if ((integrate.eq.'RESPA1').or.(integrate.eq.'BAOABRESPA1')) 
     $  then
        use_shortmlist = use_mlist
        use_shortclist = use_clist
        use_shortvlist = use_vlist
      end if

c    
c     initialization for Langevin dynamics
c
      if ((integrate .eq. 'BBK').or.(integrate.eq.'BAOAB').or.
     $ (integrate.eq.'BAOABRESPA').or.(integrate.eq.'BAOABRESPA1').or.
     $ (integrate.eq.'BAOABPISTON'))
     $   then
         if (.not.(isothermal)) then
           if (rank.eq.0) then
             write(*,*) 'Langevin integrators only available with NVT',
     $ ' or NPT ensembles'
            end if
           call fatal
         end if
         thermostat = 'none'
         if (allocated(Rn)) deallocate (Rn)
         allocate (Rn(3,nloc))
         do i = 1, nloc
           do j = 1, 3
             Rn(j,i) = normal()
           end do
         end do
         dorest = .false.
      end if
c
c     set the number of degrees of freedom for the system
c
      if (nfree .eq. 0) then
         nfree = 3 * nuse
         if (isothermal .and. thermostat.ne.'ANDERSEN') then
            if (use_bounds) then
               nfree = nfree - 3
            else
               nfree = nfree - 6
            end if
         end if
         if (use_rattle) then
            nfree = nfree - nrat
            do i = 1, nratx
               nfree = nfree - kratx(i)
            end do
         end if
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
         call reassignpme(.false.)
         call mechanicstep(0)
         call allocstep
         call nblist(0)

         if (barostat.eq.'LANGEVIN') then
           integrate = 'BAOABPISTON'
           call plangevin()
         end if
c    
c     set velocities and fast/slow accelerations for RESPA method
c
      else if ((integrate.eq.'RESPA').or.
     $   (integrate.eq.'BAOABRESPA'))
     $  then
c
c     Do the domain decomposition
c
         call ddpme3d
         call reinitnl(0)
         call reassignpme(.false.)
         call mechanicstep(0)
         call nblist(0)
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
                  v(j,iglob) = 0.0d0
                  a(j,iglob) = 0.0d0
                  aalt(j,iglob) = 0.0d0
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
c     set velocities and fast/inter/slow accelerations for RESPA-n method
c
      else if ((integrate.eq.'RESPA1').or.(integrate.eq.'BAOABRESPA1'))
     $ then
c
c     Do the domain decomposition
c
         call ddpme3d
         call reinitnl(0)
         call reassignpme(.false.)
         call mechanicstep(0)
         call nblist(0)
         allocate (derivs(3,nbloc))
         derivs = 0d0
         if (allocated(desave)) deallocate (desave)
         allocate (desave(3,nbloc))
         desave = 0d0
         call allocsteprespa(.false.)

         call gradfast1(e,derivs)
         call commforcesrespa1(derivs,0)
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
               do j = 1, 3
                  aalt2(j,iglob) = -convert *
     $              derivs(j,i) / mass(iglob)
               end do

            else
               do j = 1, 3
                  aalt2(j,iglob) = 0.0d0
               end do
            end if
         end do

         derivs = 0d0
         call gradint1(e,derivs)
         call commforcesrespa1(derivs,1)
         do i = 1, nloc
            iglob = glob(i)
            if (use(iglob)) then
               do j = 1, 3
                  aalt(j,iglob) = -convert *
     $              derivs(j,i) / mass(iglob)
               end do
            else
               do j = 1, 3
                  aalt(j,iglob) = 0.0d0
               end do
            end if
         end do

         derivs = 0d0
         call gradslow1(e,derivs)
         call commforcesrespa1(derivs,2)
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
                  v(j,iglob) = 0.0d0
                  a(j,iglob) = 0.0d0
               end do
            end if
         end do

         deallocate (derivs)
         if (nuse .eq. n)  call mdrest (0)
c
      else
c
c     Do the domain decomposition
c
        call ddpme3d
        call reinitnl(0)
        call reassignpme(.false.)
        call mechanicstep(0)
        call nblist(0)
c    
c     set velocities and accelerations for cartesian dynamics
c
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
