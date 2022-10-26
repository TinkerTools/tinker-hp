c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  program dynamic  --  run molecular or stochastic dynamics  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dynamic" computes a molecular dynamics trajectory
c     in one of the standard statistical mechanical ensembles and using
c     any of several possible integration methods
c
c
      program dynamic
      use mpi
      implicit none
      integer ierr
      call MPI_INIT(ierr)
      call dynamic_bis
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      end
c
      subroutine dynamic_bis
      use atoms
      use bath
      use bound
      use domdec
      use keys
      use inform
      use iounit
      use mdstuf
      use moldyn
      use mpi
      use potent
      use timestat

#ifdef COLVARS
      use colvars
#endif
      implicit none
      integer i,istep,nstep,ierr
      integer mode,next
      real*8 dt,dtdump,time0,time1
      logical exist,query
      character*20 keyword
      character*240 record
      character*240 string
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
      nproc = nproctot
      call initmpi
      call unitcell
      call cutoffs
      call lattice
c
c     setup for MPI
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
      if (allocated(aalt2)) deallocate (aalt2)
      allocate (aalt2(3,n))
c
      a = 0d0
      v = 0d0
      aalt = 0d0
      aalt2 = 0d0
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
      integrate = 'BEEMAN'
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         end if
      end do
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
         if (mode.eq.3 .or. mode.eq.4) then
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
      call mdinit(dt)

#ifdef COLVARS
      dt_sim = dt*1000d0
c
c     only the master does colvars computations, but other ranks need to allocate coord arrays
c
      if (rank.eq.0) then
         call allocate_colvars
      end if
      call MPI_BCAST(use_colvars,1,MPI_LOGICAL,0,COMM_TINKER,ierr)
      if (use_colvars) then
        call MPI_BCAST(ncvatoms,1,MPI_INT,0,COMM_TINKER,ierr)
        if (rank.gt.0) then
          allocate (cvatoms_ids(ncvatoms))
        end if
        call MPI_BCAST(cvatoms_ids,ncvatoms,MPI_INT,0,COMM_TINKER,
     $       ierr)
        if (rank.gt.0) then
           allocate (cv_pos(3,ncvatoms))
           allocate (decv(3,ncvatoms),decv_tot(3,ncvatoms))
           cv_pos = 0d0
           decv = 0d0
           decv_tot = 0d0
        end if
      end if
#endif
      if (use_lambdadyn .and. .not.(use_colvars)) then
        if (rank.eq.0) then
          write(iout,*) 'cannot run lambda dynamics without colvars'
        end if
        call fatal
      end if
c
c     print out a header line for the dynamics computation
c
      if (integrate .eq. 'VERLET') then
         if (rank.eq.0) write (iout,330)
  330    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'RESPA') then
         if (rank.eq.0) write (iout,390)
  390    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA MTS Algorithm')
      else if (integrate .eq. 'BBK') then
         if (rank.eq.0) write (iout,400)
  400    format (/,' Langevin Molecular Dynamics Trajectory via',
     &              ' BBK Algorithm')
      else if (integrate .eq. 'BAOAB') then
         if (rank.eq.0) write (iout,410)
  410    format (/,' Langevin Molecular Dynamics Trajectory via',
     &              ' BAOAB Algorithm')
      else if (integrate .eq. 'BAOABRESPA') then
         if (rank.eq.0) write (iout,420)
  420    format (/,' Langevin Molecular Dynamics Trajectory via',
     &              ' BAOAB-RESPA Algorithm')
      else if (integrate .eq. 'BAOABRESPA1') then
         if (rank.eq.0) write (iout,430)
  430    format (/,' Langevin Molecular Dynamics Trajectory via',
     &              ' BAOAB-RESPA-1 Algorithm')
      else if (integrate .eq. 'RESPA1') then
         if (rank.eq.0) write (iout,440)
  440    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA-1 MTS Algorithm')
      else if (integrate .eq. 'BAOABPISTON') then
         if (rank.eq.0) write (iout,450)
  450    format (/,' Molecular Dynamics Trajectory via',
     &              ' BAOAB-PISTON Algorithm')
      else
         if (rank.eq.0) write (iout,460)
  460    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if
c
c     integrate equations of motion to take a time step
c
      do istep = 1, nstep
         time0 = mpi_wtime()
         if (integrate .eq. 'VERLET') then
            call verlet (istep,dt)
         else if (integrate .eq. 'RESPA') then
            call respa (istep,dt)
         else if (integrate .eq. 'BBK') then
           call bbk(istep,dt)
         else if (integrate .eq. 'BAOAB') then
           call baoab(istep,dt)
         else if (integrate .eq. 'BAOABPISTON') then
           call baoabpiston(istep,dt)
         else if (integrate .eq. 'BAOABRESPA') then
           call baoabrespa(istep,dt)
         else if (integrate .eq. 'BAOABRESPA1') then
           call baoabrespa1(istep,dt)
         else if (integrate .eq. 'RESPA1') then
           call respa1(istep,dt)
         else
            call beeman (istep,dt)
        end if
         time1 = mpi_wtime()
         timestep = timestep + time1-time0
         call mdstattime(istep,dt)
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
