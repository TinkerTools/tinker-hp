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
      use beads, only: path_integral_md
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
      use qtb
      use potent
      use units 
      use uprior
      use usage
      use virial
      use replicas
#ifdef PLUMED
      use plumed
      use mpi
#endif
      use math, only:pi
      use spectra
      implicit none
      integer i,j,idyn,iglob,ierr
      integer next
      integer lext,freeunit,ios
      real*8 e
      real*8 speed
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
      interface
       function maxwell (mass,temper)
       real*8 maxwell
       real*8 mass
       real*8 temper
       end function
      end interface
c
c     set default parameters for the dynamics trajectory
c
      integrate = 'BEEMAN'
      mts = .false.
      bmnmix = 8
      nfree = 0
      dorest = .true.
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
      if(path_integral_md) then
         barostat = 'LANGEVIN'
         thermostat = 'LANGEVIN'
      else
        barostat   = 'BERENDSEN'
        thermostat = 'BUSSI'
      endif
      tautemp = 0.2d0
      collide = 0.1d0
      anisotrop = .false.
      taupres = 2.0d0
      compress = 0.000046d0
      vbar = 0.0d0
      qbar = 0.0d0
      gbar = 0.0d0
      eta = 0.0d0
      voltrial = 20
      volmove = 100.0d0
c      volscale = 'MOLECULAR'
      volscale = 'ATOMIC'
      gamma = 1.0d0
      gammapiston = 20.0d0
      masspiston=100000.d0 
      virsave = 0d0

      ! spectra parameters
      omegacut = 15000.d0/cm1 !0.45*pi/dt   
      Tseg = 1
      ir=.false.
      full_dipole=.TRUE.
      deconvolute_IR = .FALSE.
      niter_deconvolution = 20
      startsavespec=25

      qtb_thermostat=.false.
      adaptive_qtb=.false.
c
c     set default values for lambda dynamic
c
      use_lambdadyn = .false.
      use_OSRW = .false.
c
c     set default for pbc unwrapping
c
      pbcunwrap = .false.
c
c     allocate array of positions to be printed
c
      if (allocated(xwrite)) deallocate (xwrite)
      allocate (xwrite(n))
      xwrite = 0d0
      if (allocated(ywrite)) deallocate (ywrite)
      allocate (ywrite(n))
      ywrite = 0d0
      if (allocated(zwrite)) deallocate (zwrite)
      allocate (zwrite(n))
      zwrite = 0d0
c
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
        ios = 0
        next = 1
        record = keyline(i)
        call gettext (record,keyword,next)
        call upcase (keyword)
        string = record(next:240)
        
        select case(trim(keyword))
        case ('INTEGRATOR')
          call getword (record,integrate,next)
          call upcase (integrate)
        case('BEEMAN-MIXING')
          read (string,*,iostat=ios)  bmnmix
        case ('REMOVE-INERTIA')
          read (string,*,iostat=ios)  irest
        case ('SAVE-VELOCITY')
          velsave = .true.
        case ('SAVE-FORCE')
          frcsave = .true.
        case ('SAVE-INDUCED')
          uindsave = .true.
        case ('POLAR-PREDICT')
          use_pred = .true.
          call getword (record,polpred,next)
          call upcase (polpred)
          if (polpred .eq. 'NONE') use_pred = .false.
        case ('THERMOSTAT')
          call getword (record,thermostat,next)
          call upcase (thermostat)
        case ('TAU-TEMPERATURE')
          read (string,*,iostat=ios)  tautemp
        case ('COLLISION')
          read (string,*,iostat=ios)  collide
        case ('BAROSTAT')
          call getword (record,barostat,next)
          call upcase (barostat)
        case ('ANISO-PRESSURE')
          anisotrop = .true.
          write(*,*) "ANISO-PRESSURE"
        case ('TAU-PRESSURE')
          read (string,*,iostat=ios)  taupres
        case ('COMPRESS')
          read (string,*,iostat=ios)  compress
        case ('VOLUME-TRIAL')
          read (string,*,iostat=ios)  voltrial
        case ('VOLUME-MOVE')
          read (string,*,iostat=ios)  volmove
        case ('PRINTOUT')
          read (string,*,iostat=ios)  iprint
        case ('PBCUNWRAP')
          pbcunwrap = .true.
        case ('NLUPDATE')
          read (string,*,iostat=ios) ineigup
        case ('FRICTION')
          read (string,*,iostat=ios) gamma
        case ('OMEGACUT')
          read (string,*,iostat=ios) omegacut
          omegacut = omegacut/cm1
        case ('TSEG')
          read (string,*,iostat=ios) TSEG
        case ('FRICTIONPISTON')
          read (string,*,iostat=ios) gammapiston
        case ('MASSPISTON')
          read (string,*,iostat=ios) masspiston
        case ('VIRNUM')
          virnum=.true.
          if(ranktot.eq.0) then
            write(*,*) 'WARNING - virial pressure',
     &        ' computed using finite differences.'
          endif
        case ('KIN_PRES_AVG')
          kin_instant=.false.
        case ('STARTSAVESPEC')
          read (string,*,iostat=ios) startsavespec
        case ('IR_SPECTRA')
          ir = .true.
        case ('IR_LINEAR_DIPOLE')
          full_dipole = .false.
        case ('IR_DECONVOLUTION')
          deconvolute_IR = .true.
        case ('NITER_DECONV')
          read (string,*,iostat=ios) niter_deconvolution
        case ('LAMBDADYN')
          use_lambdadyn = .true.
        case ('OSRW')
          use_lambdadyn = .true.
          use_OSRW = .true.
#ifdef PLUMED
c
c initialise PLUMED
c
        case ('PLUMED')
          lplumed = .true.
          allocate(pl_force(3,n))
          allocate(pl_pos  (3,n))
          allocate(pl_mass (  n))
          allocate(pl_glob (  n))
          call getword (record,pl_input,next)
          call getword (record,pl_output,next)
c          if (rank.eq.0) then
c            write(0,*)"  CREATING PLUMED FROM THE PROGRAM"
c          endif
c
c          conversion factors to go from TINKER units to PLUMED
          energyUnits = 4.184d0 ! kcal/mol -> kJ/mol
          lengthUnits = 0.1d0 ! angstrom -> nm
          timeUnits = 0.001d0 ! fs      
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
          call plumed_f_gcmd("setNatoms"//char(0),n)
          call plumed_f_gcmd("setMDEngine"//char(0),"TinkerHP");
          call plumed_f_gcmd("setTimestep"//char(0),dt);
          call plumed_f_gcmd("init"//char(0),0);
c
c end PLUMED initialisation
c
#endif
        end select

        if (ios /= 0) then
            write(*,*) "Warning: keyword ",trim(keyword)
     &         ," not correctly read!"
          endif
      end do
c
c     pbcwrapindex initialization
c
      if (allocated (pbcwrapindex)) deallocate (pbcwrapindex)
      allocate (pbcwrapindex(3,n))
      pbcwrapindex = 0
c
c     lambda dynamic initialization
c
      if ((use_lambdadyn).or.(use_osrw)) then
c
c     lambda-dynamics not compatible with pmecore
c
        if (use_pmecore) then
          if (rank.eq.0) then
            write(iout,*) 'Lambda-dynamics not compatible with',
     $      'separate cores for PME'
          end if
          call MPI_BARRIER(COMM_TINKER,ierr)
          call fatal
        end if
        call def_lambdadyn_init
      end if
c
c     enforce the use of monte-carlo barostat with the TCG family of solvers
c
      if ((isobaric).and.(polalg.eq.3)) then
        barostat = 'MONTECARLO' 
        if (ranktot.eq.0) write(iout,*)'enforcing the use of Monte ',
     $    'Carlo Barostat with TCG'
      end if
c
c     default time steps for respa and respa1 integrators
c
      if ((integrate.eq.'RESPA').or.(integrate.eq.'BAOABRESPA')) then
        eps =  0.00000001d0
        dshort = 0.00025
        mts = .true.
      else if ((integrate.eq.'RESPA1').or.(integrate.eq.'BAOABRESPA1'))
     $  then
        dinter = 0.002
        dshort = 0.00025
        mts = .true.
      end if
c
c     keywords fr respa and respa1 integrators
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
        nalt2 = 1
        dshort = dble(nalt)
        dinter=1.0d0
      else if ((integrate.eq.'RESPA1').or.(integrate.eq.'BAOABRESPA1'))
     $  then
        eps =  0.00000001d0
        nalt = int(dt/(dinter+eps)) + 1
        nalt2 = int(dinter/(dshort+eps)) + 1
        dinter = dble(nalt)
        dshort = dble(nalt2)
      else
        dshort = 1.0d0
        dinter = 1.0d0
        nalt   = 1
        nalt2  = 1
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
      if (thermostat .eq. 'ANDERSEN') then
c     enforce use of velocity Verlet with Andersen thermostat
         if (integrate .eq. 'BEEMAN')  integrate = 'VERLET'
      elseif(thermostat.eq.'QTB') then
        qtb_thermostat=.TRUE.
        adaptive_qtb=.FALSE.
      elseif(thermostat.eq.'ADQTB') then
        qtb_thermostat=.TRUE.
        adaptive_qtb=.TRUE.
      elseif(path_integral_md) then
        thermostat='LANGEVIN'
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
        use_shortdlist = use_dlist
      end if

c    
c     initialization for Langevin dynamics
c
      if ((integrate.eq.'BAOAB').or.
     $ (integrate.eq.'BAOABRESPA').or.(integrate.eq.'BAOABRESPA1'))
     $   then
         if (.not.(isothermal)) then
           if (ranktot.eq.0) then
             write(*,*) 'Langevin integrators only available with NVT',
     $ ' or NPT ensembles'
            end if
           call fatal
         end if
         if(.not.qtb_thermostat) then
           thermostat = 'LANGEVIN'
           if (allocated(Rn)) deallocate (Rn)
           allocate (Rn(3,nloc))
           do i = 1, nloc
             do j = 1, 3
               Rn(j,i) = normal()
             end do
           end do
         endif
      end if

c
c     decide whether to remove center of mass motion
c
      if (irest .eq. 0)  dorest = .false.
      if (nuse. ne. n)  dorest = .false.
      if (isothermal .and. thermostat.eq.'ANDERSEN')  dorest = .false.
      if (qtb_thermostat) dorest = .false. 
c
c     set the number of degrees of freedom for the system
c
      nfree = 3 * nuse
      if (dorest) then
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
c     try to restart using prior velocities and accelerations
c
      dynfile = filename(1:leng)//'.dyn'
      call version (dynfile,'old')
      inquire (file=dynfile,exist=exist)
      if(.not. exist) then
 33   format (/," Tinker-HP find No Restart file for ",A," system",/,
     $       6x,"--- Initiate dynamic from begining ---",/)
         if (verbose.and.ranktot.eq.0) write(*,33) filename(1:leng)
      endif

      if (exist) then
 32   format(/," --- Tinker-HP: loading restart file for ",A,
     &         " system --- ")
         if (verbose.and.rank.eq.0) then
            write(*,32) filename(1:leng)
         end if
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
         call reassignpme(.true.)
         call mechanic_init_para
         call allocstep
         call nblist(0)
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
         call mechanic_init_para
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
         call mechanic_init_para
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
        call mechanic_init_para
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

      ! initialize langevin piston 
      if (isobaric .AND. barostat.eq.'LANGEVIN') then
        use_piston = .TRUE.
        call initialize_langevin_piston()
      end if

      !dshort=1. if not respa(1); !dinter=1. if not respa1
      if(ir) then
       if (use_reps) then
          write(0,*) "Error: IR not compatible with replicas yet"
          call fatal
        endif
       call initialize_spectra(dt/dshort/dinter)
      endif
      
c     INITIALIZE QTB
      if(qtb_thermostat) then
        call check_qtb_compatibility()
        call qtbinit(dt/dinter/dshort)
      else
         corr_fact_qtb=1.d0
         use_corr_fact_qtb=.FALSE.
      endif

      if (use_reps) call mdinitreps

      return
      end
c

      subroutine check_qtb_compatibility()
      use domdec
      use replicas, only : use_reps
      use freeze, only: use_rattle
      use bath, only: barostat, isobaric
      use mdstuf, only: integrate
      implicit none

      if(ranktot/=0) return

      !if(path_integral_md) then
      !  write(0,*) 'PIQTB not available yet!'
      !  call fatal
      !endif

      if (use_reps) then
        write(0,*) "Error: QTB not compatible with replicas yet"
        call fatal
      endif

      if(use_rattle) then
        write(0,*) "Error: QTB not compatible with RATTLE"
        call fatal
      endif

      if(isobaric .AND. barostat /= 'LANGEVIN') then
        write(0,*) "Error: QTB only compatible",
     &      " with LANGEVIN barostat"
        call fatal
      endif

      if(.NOT. (integrate.eq."BAOAB"
     &     .OR. integrate.eq."BAOABRESPA"
     &     .OR. integrate.eq."BAOABRESPA1" )) then

        write(0,*) "Error: QTB only compatible with ",
     &         "BAOAB, BAOABRESPA and BAOABRESPA1 integrators"
        call fatal
      endif

      end subroutine check_qtb_compatibility

