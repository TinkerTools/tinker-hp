c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  program pimd  --  run pimd molecular or stochastic dynamics  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "pimd" computes a molecular dynamics trajectory
c     in one of the standard statistical mechanical ensembles and using
c     any of several possible integration methods
c
c
      program pimd
      use mpi
      implicit none
      integer ierr,nthreadsupport
c      call MPI_INIT(ierr)
      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
      call pimd_bis
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      end
c
      subroutine pimd_bis
      use atoms
      use bath
      use beads
      use bound
      use boxes, only: volbox
      use domdec
      use keys
      use inform
      use iounit
      use mdstuf
      use moldyn
      use neigh
      use timestat
      use mpi
      implicit none
      integer i,istep,nstep,ierr,k,ibead,ibeadglob!,nstep_therm
      integer mode,next,nseg
      real*8 dt,dtdump,time0,time1
      real*8 timesteppimd
      real*8 ekvir,ekprim,presvir
      logical exist,query,restart
      character*20 keyword
      character*120 record
      real*8 pres_tmp
      character*120 string
      real*8 maxwell
      logical isobaric_save
c
c
 1000 Format(' Time for 100 Steps: ',f15.4,/,
     $  ' Ave. Time per step: ',f15.4)
 1010 Format(' ns per day: ',f15.4)
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call getxyz
c
c     check for keywords containing any altered parameters
c
      nbeads = 1 
      nbeadsloc = 1 
      nproc = 1
!      nstep_therm=0
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:7) .eq. 'NBEADS ') then
            read (string,*,err=5,end=5) nbeads
!         else if (keyword(1:12) .eq. 'NSTEP_THERM ') then
!            read (string,*,err=5,end=5) nstep_therm
         end if
   5  continue
      end do

      if (nproctot.lt.nbeads) then
        nproc = 1
      else if (mod(nproctot,nbeads).ne.0) then
        if (ranktot.eq.0) then
          write(iout,*) 'inconsistent number of process for parallelism'
          write(iout,*) 'the total number of processors should be lower
     $     or a multiple of nbeads'
          call fatal
        end if
      else
        nproc = nproctot/nbeads
      end if

c
      call initmpi
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call cutoffs
      call unitcell
      call lattice
c
c     setup for MPI
c
c
c     get nprocforce to do correct allocation of MPI arrays
c
      call drivermpi
      call reinitnl(0)
c
c     allocate some arrays
c
      if (allocated(v)) deallocate (v)
      allocate (v(3,n))
      if (allocated(a)) deallocate (a)
      allocate (a(3,n))
      if (allocated(aalt)) deallocate (aalt)
      allocate (aalt(3,n))
c
      a = 0d0
      v = 0d0
      aalt = 0d0
c
      call mechanic
      call nblist(0)
c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
c
c     check for keywords containing any altered parameters
c
      integrate = 'BAOAB'
c      nbeads = 1 
c      nprocbeads = 1
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
c         else if (keyword(1:7) .eq. 'NBEADS ') then
c            read (string,*,err=5,end=5) nbeads
c         else if (keyword(1:11) .eq. 'NPROCBEADS ') then
c            read (string,*,err=5,end=5) nprocbeads
         end if
      end do
c      write(*,*) 'nbeads = ',nbeads
c      write(*,*) 'nbeadsloc = ',nbeadsloc
c      if (nprocbeads.gt.nbeads) then
c        write(iout,*) 'too many processes for beads parallelism'
c      else
c         nbeadsloctemp = nbeads/nprocbeads
c      end if
c
c
c     initialize the simulation length as number of time steps
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  nstep
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30)  nstep
   30    format (i10)
      end if

c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      do while (dt .lt. 0.0d0)
         write (iout,50)
   50    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,60,err=70)  dt
   60    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   70    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  dtdump
   80 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,90)
   90    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,100,err=110)  dtdump
  100    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  110    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=120,end=120)  mode
  120    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,130)
  130       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,140,err=150)  mode
  140       format (i10)
            if (mode .le. 0)  mode = 1
  150       continue
         end do
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=170,end=170)  kelvin
  170       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,180)
  180          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,190,err=200)  kelvin
  190          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  200          continue
            end do
         end if
        if (mode.eq.1 .or. mode.eq.3) then
            if(ranktot.eq.0) then
               write(iout,*) 'NVE, NPH and NPT not implemented yet!'
               call fatal
            endif
        endif
        if(mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=210,end=210)  atmsph
  210       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,220)
  220          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,230,err=240)  atmsph
  230          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  240          continue
            end do
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  mode
  250    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,260)
  260       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,270,err=280)  mode
  270       format (i10)
            if (mode .le. 0)  mode = 1
  280       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=290,end=290)  kelvin
  290       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,300)
  300          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,310,err=320)  kelvin
  310          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  320          continue
            end do
         end if
      end if
c
c     setup dynamics
c
      call mdinit(dt)      
      call allocpi(.true.,0)      

      do ibead = 1, nbeadsloc
        call mdinitbead(ibead,dt,restart)       
        if(.not.restart) then
          call dedvcalc()
          call savebeadnl(ibead,0)          
        endif
        call initbead(ibead,0)
      end do


      if(isobaric) then
        extvol = volbox
        extvolold = volbox
        masspiston=masspiston*nbeads
        if(ranktot.eq.0) then
          vextvol = maxwell(masspiston,kelvin)
          aextvol = 0d0
        endif
      endif      
      
      if (restart) then
        call allocpi(.false.,0)
        do ibead=1,nbeadsloc
          call pushbead(ibead,0)
          call nblist(0)
          call dedvcalc()
          call savebeadnl(ibead,0)          
          call initbead(ibead,0)
        enddo
      endif
c
c     print out a header line for the dynamics computation
c

       if (ranktot.eq.0) then
        write (iout,'(A,A,A,i3,A)') ' Path Integral Molecular Dynamics',
     &       ' Trajectory via BAOAB Algorithm',
     &       ' with ',nbeads,' beads'
        endif

c
c     integrate equations of motion to take a time step
c

      
!      if(nstep_therm>0) then
!        isobaric_save=isobaric
!        isobaric=.false.
!        if(ranktot.eq.0) then
!          write(*,*) nstep_therm, "STEPS OF THERMALIZATION."
!        endif
!        do istep=1,nstep_therm
!          call integrate_pimd(istep,dt)
!        enddo
!        if(ranktot.eq.0) write(*,*) "THERMALIZATION DONE."
!        isobaric=isobaric_save
!      endif

      time0 = mpi_wtime()
      do istep = 1, nstep
        call integrate_pimd(istep,dt,ekvir,ekprim,presvir)
        if (mod(istep,iprint).eq.0) then
        time1 = mpi_wtime()
        timestep = time1-time0
        if (ranktot.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,timestep,1,MPI_REAL8,MPI_SUM,
     $      0,MPI_COMM_WORLD,ierr)
        else
         call MPI_REDUCE(timestep,timestep,1,MPI_REAL8,MPI_SUM,0,
     $      MPI_COMM_WORLD,ierr)
        end if
        if (ranktot.eq.0) then
         write(6,1000) (timestep)/nproctot,(timestep)/
     $       dble(nproctot*iprint)
         write(6,1010) 86400*dt*dble(iprint*nproctot)/(1000*timestep)
          end if
        time0 = time1
        end if
      end do

c     perform any final tasks before program exit
c
      call final

      end subroutine

      subroutine integrate_pimd(istep,dt,ekvir,ekprim,presvir)
        use atoms
        use bath
        use beads
        use bound
        use boxes, only: volbox
        use domdec
        use keys
        use inform
        use iounit
        use mdstuf
        use moldyn
        use timestat
        use mpi
        implicit none
        integer, intent(in) :: istep
        real*8, intent(in) :: dt
        real*8 ekvir,ekprim,presvir
        integer :: ibead
        real*8 time0,time1,time00,time01
        real*8 timesteppimd, timebeads,timepush,timereass,timeinit
        real*8 timededvcalc,timebaoab2,timeaoa
        real*8 pres_tmp

        timesteppimd=0d0
        timepush=0d0
        timereass=0d0
        timeinit=0d0
c
c     loop on local beads
c
        time0 = mpi_wtime()
        do ibead = 1, nbeadsloc
c
c     put bead variable in global ones used for MD
c
c
          call pushbead(ibead,1)
          if (integrate.eq.'BAOAB4SITE') then
            call baoabpi4site1(istep,dt)
          else
            call baoabpi1(istep,dt)
          end if
c          timestep = timestep+time1-time0
          call initbead(ibead,istep)
        end do
        time1 = mpi_wtime()

       timesteppimd=timesteppimd+time1-time0
         
c        PROPAGATE AOA IN NORMAL MODES
         if (integrate.eq.'BAOAB4SITE') then
           call aoapi4site(istep,dt)
         else
           call aoapi(istep,dt,ekvir,ekprim,presvir)
         end if
         time01=mpi_wtime()
         timeaoa=time01-time00

        
c        reassign all positions in beads before 
c        reconstruction of neighborlist
      time0=mpi_wtime()
        do ibead = 1, nbeadsloc
          time00=mpi_wtime()
          call pushbead(ibead,istep)
          time01=mpi_wtime()
          timepush=timepush+time01-time00
          call reassignpi(istep)
          time00=mpi_wtime()
          timereass=timereass+time00-time01
          call initbead(ibead,istep)
          time01=mpi_wtime()
          timeinit=timeinit+time01-time00
        enddo 
        time1 = mpi_wtime()   
        timebeads=time1-time0

        call allocpi(.false.,istep)

         
        do ibead = 1, nbeadsloc
            ibeadsloc=ibead
            time00=mpi_wtime()
          call pushbead(ibead,istep)
            time01=mpi_wtime()
            timepush=timepush+time01-time00

          time00=mpi_wtime()
          if (integrate.eq.'BAOAB4SITE') then
            call baoabpi4site2(istep,dt,ibeadsloc)
          else
            call baoabpi2(istep,dt)
          end if
          time01=mpi_wtime()
          timebaoab2=timebaoab2+time01-time00
          time00=mpi_wtime()
          call dedvcalc()
          time01=mpi_wtime()
          timededvcalc=timededvcalc+time01-time00
          call savebeadnl(ibead,istep)
         ! call initbead(ibead,1)
          time00=mpi_wtime()
         call initbead(ibead,istep)
         time01=mpi_wtime()
         timeinit=timeinit+time01-time00

        end do
c        if(ranktot.eq.0) then
c            write(*,*) 'Time baoab 1', timesteppimd
c            write(*,*) 'Time pushbeads etc ', timebeads
c            write(*,*) 'Time pushbead ', timepush
c            write(*,*) 'Time reassigpi etc ', timereass
c            write(*,*) 'Time initbead etc ', timeinit
c            write(*,*) 'Time baoab 2 (gradient)  ', timebaoab2
c            write(*,*) 'Time dedv ', timededvcalc
c            write(*,*) 'Time aoa ', timeaoa
c      endif

      end subroutine
