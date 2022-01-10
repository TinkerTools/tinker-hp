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
      use domdec
      use files
      use keys
      use freeze
      use inform
      use iounit
      use langevin
      use math
      use mdstuf
      use molcul
      use moldyn
      use mpole
      use neigh
      use units 
      use uprior
      use adqtb
      use qtb
      use usage
      use virial
      implicit none
      integer i,j,k,idyn,iglob,nh,nm
      integer next
      integer lext,freeunit
      integer ierr
      real*8 e
      real*8 maxwell,speed
      real*8 hmax,hmass
      real*8 sum,dmass
      real*8 normal
      real*8 vec(3)
      real*8 dt
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
      voltrial = 25
      volmove = 100.0d0
      volscale = 'MOLECULAR'
      virnum=.false.
      kin_instant=.true.
c      volscale = 'ATOMIC'
      gamma = 1.0d0
      omegamax = pi/dt
      omegacut = 0.45*omegamax
      domega = 10d0
      omegasmear = omegamax/100.0
      Tseg = 1
      A_gamma=0.1
      compteur=0
      skipseg=3
      gammapiston=20.0
      masspiston=100000
      noQTB=.FALSE.
      corr_pot=.false.
      adaptive=.false.
      startsavespec=25
      corr_fact_qtb=1.0d0
      register_spectra=.false.
      A_gamma_piston=0.1
      ir=.false.
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
         else if (keyword(1:14) .eq. 'NLUPDATE ') then
            read (string,*,err=10,end=10) ineigup
         else if (keyword(1:9) .eq. 'FRICTION ') then
            read (string,*,err=10,end=10) gamma
         else if (keyword(1:9) .eq. 'OMEGACUT ') then
            read (string,*,err=10,end=10) omegacut
         else if (keyword(1:11) .eq. 'OMEGASMEAR ') then
            read (string,*,err=10,end=10) omegasmear
         else if (keyword(1:9) .eq. 'OMEGAMAX ') then
            read (string,*,err=10,end=10) omegamax
         else if (keyword(1:5) .eq. 'TSEG ') then
            read (string,*,err=10,end=10) TSEG
         else if (keyword(1:8) .eq. 'A_GAMMA ') then
            read (string,*,err=10,end=10) a_gamma
         else if (keyword(1:8) .eq. 'SKIPSEG ') then
            read (string,*,err=10,end=10) skipseg
         else if (keyword(1:15) .eq. 'FRICTIONPISTON ') then
            read (string,*,err=10,end=10) gammapiston
         else if (keyword(1:12) .eq. 'MASSPISTON ') then
            read (string,*,err=10,end=10) masspiston
         else if (keyword(1:6) .eq. 'ADQTB ') then
            adaptive = .true.
         else if (keyword(1:6) .eq. 'NOQTB ') then
            noQTB = .true.
         else if (keyword(1:7) .eq. 'VIRNUM ') then
            virnum=.true.
            if(ranktot.eq.0) then
              write(*,*) 'WARNING - virial pressure',
     &               ' computed using finite differences.'
            endif
         else if (keyword(1:13) .eq. 'KIN_PRES_AVG ') then
            kin_instant=.false.
         else if (keyword(1:14) .eq. 'STARTSAVESPEC ') then
            read (string,*,err=10,end=10) startsavespec
         else if (keyword(1:14) .eq. 'CORR_FACT_QTB ') then
            read (string,*,err=10,end=10) corr_fact_qtb
            write(*,*) 'WARNING _ USed of Kinetic correction'
            write(*,*) 'corr_fact_qtb=', corr_fact_qtb
         else if (keyword(1:9) .eq. 'CORR_POT ') then
            corr_pot = .true.
            write(*,*) 'WARNING -  used of potential',
     &              ' correction.'
         else if (keyword(1:17) .eq. 'REGISTER_SPECTRA ') then
            write(*,*) 'Register the different spectras while',
     &            ' using the QTB method'
            register_spectra = .true.
         else if (keyword(1:8) .eq. 'A_GAMMA_piston ') then
            read (string,*,err=10,end=10) a_gamma_piston
         else if (keyword(1:11) .eq. 'IR_SPECTRA ') then
            ir = .true.
         end if
   10    continue
      end do
c

      if ((register_spectra).and.(adaptive.eq. .false.)) then
        adaptive=.true.
        a_gamma=0
        a_gamma_piston=0
      endif
          
      if(ir) then
            nseg=nint(Tseg/dt)
            Tseg=nseg*dt
            allocate(vad(3,n,nseg))
      endif
        
      if ((rank.eq.0).and.adaptive) then
        write(*,*) 'Using adQTB with adaptive gamma coefficients'
        write(*,*) 'a_gamma=', a_gamma
        write(*,*) 'a_gamma_piston=', a_gamma_piston
      endif

c     enforce the use of baoabpiston integrator with Langevin Piston
c     thermostat
c      
       if (isobaric) then
           
         if (integrate.eq.'QTB4SITE') then
           integrate = 'PQTB4SITE'
            if(barostat/='LANGEVIN') then
              barostat='LANGEVIN'
              write(*,*) 'WARNING - enforced langevin piston for QTB'
            endif
         elseif (integrate.eq.'QTB') then
            integrate = 'PISTONQTB'
            if(barostat/='LANGEVIN') then
              barostat='LANGEVIN'
              write(*,*) 'WARNING - enforced langevin piston for QTB'
            endif
         endif    
       endif
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
            mass(i) = 1.d0           
            if(atomic(i)/=0) then
              totmass = totmass + 1.0d0
              write (iout,30)  i
   30         format (/,' MDINIT  --  Warning, Mass of Atom',i6,
     &                 ' Set to 1.0 for Dynamics')
            endif
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
c     initialization for Langevin dynamics
c
      if ((integrate .eq. 'BBK').or.(integrate.eq.'BAOAB').or.
     $ (integrate.eq.'BAOABRESPA').or.(integrate.eq.'QTB').or.
     $  (integrate.eq.'BAOAB4SITE').or.
     $  (integrate.eq.'QTB4SITE').or.(integrate.eq.'PQTB4SITE'))
     $  then
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
         
          nm = 0
          do i = 1, n
             if (atomic(i) .eq. 0)  nm = nm + 1
          end do
          nfree = 3*(n - nm)         
         
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
      dorest=.true. 
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
         call drivermpi
         call reinitnl(0)
         call reassignpme(.false.)
         call mechanicstep(0)
         call nblist(0)

c    
c     initialization for Langevin dynamics
c
         if (integrate .eq. 'BBK') then
            thermostat = 'none'
            if (allocated(Rn)) deallocate (Rn)
            allocate (Rn(3,nloc))
            do i = 1, nloc
              do j = 1, 3
                Rn(j,i) = normal()
              end do
            end do
         end if
            
c
c     set velocities and fast/slow accelerations for RESPA method
c
      else if ((integrate.eq.'RESPA').or.
     $   (integrate.eq.'BAOABRESPA')) then
c
c     Do the domain decomposition
c
         call drivermpi
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
      else
c
c     Do the domain decomposition
c
        call drivermpi
        call reinitnl(0)
        call reassignpme(.false.)
        call mechanicstep(0)
        call nblist(0)
c    
c       initialization for Langevin dynamics
c
        if (integrate .eq. 'BBK') then
           thermostat = 'none'
           if (allocated(Rn)) deallocate (Rn)
           allocate (Rn(3,nloc))
           do i = 1, nloc
             do j = 1, 3
               Rn(j,i) = normal()
             end do
           end do
        end if
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
c     initialization for QTB
c
      if ((integrate.eq.'QTB')
     &  .or.(integrate.eq.'PISTONQTB').or.(integrate.eq.'QTB4SITE').or.
     $   (integrate.eq.'PQTB4SITE'))
     $    then
            if(.not.kin_instant) then
              if(ranktot.eq.0) then
                write(*,*) "WARNING - enforcing the use of "
     &             ,"velocities for pressure estimation"
              endif
              kin_instant=.TRUE.
            endif
            nseg=nint(Tseg/dt)
            Tseg=nseg*dt
            call qtbinit(dt)
      endif   
      
      if (isobaric) then
        if (integrate.eq.'PQTB4SITE') call init_langevin_piston()
        if (integrate.eq.'PISTONQTB') call init_langevin_piston()
      endif
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
