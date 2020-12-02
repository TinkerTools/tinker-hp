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
#include "tinker_precision.h"
      program dynamic
      use mpi
#ifdef _OPENACC
      use utilgpu,only: bind_gpu
#endif
      implicit none
      integer ierr,nthreadsupport
#ifdef _OPENACC
      call bind_gpu
#endif
c      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
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
      use utils
      use utilgpu ,only: rec_queue
      use timestat
      use tinMemory
      implicit none
      integer i,istep,nstep,ierr
      integer mode,next
      real(r_p) dt,dtdump
      real(8) time0,time1,timestep
      logical exist,query
      character*20 keyword
      character*120 record
      character*120 string
c
c
 1000 Format(' Time for ',I4,' Steps :',f12.4)
 1010 Format(' ns per day',10x,':',f12.4)
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call initmpi
      call getxyz
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
      call prmem_requestm(    v,3,n)
      call prmem_requestm(    a,3,n)
      call prmem_requestm( aalt,3,n)
      call prmem_requestm(aalt2,3,n)
c
      call set_to_zero2m(   a,    v,3*n,rec_queue)
      call set_to_zero2m(aalt,aalt2,3*n,rec_queue)
c
      call mechanic
      call nblist(0)
c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0_re_p
      atmsph = 0.0_re_p
      isothermal = .false.
      isobaric = .false.
c
c     check for keywords containing any altered parameters
c
      integrate = 'VERLET'
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
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
         write (iout,'(a)',advance='no')
     &         ' Enter the Number of Dynamics Steps to be',
     &         ' Taken :  '
         read (input,30)  nstep
   30    format (i10)
      end if
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0_re_p
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      do while (dt .lt. 0.0_re_p)
         write (iout,'(a)',advance='no')
     &         ' Enter the Time Step Length in Femtoseconds',
     &         ' [1.0] :  '
         read (input,60,err=70)  dt
   60    format (f20.0)
         if (dt .le. 0.0_re_p)  dt = 1.0_re_p
   70    continue
      end do
      dt = 0.001_re_p * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0_re_p
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  dtdump
   80 continue
      do while (dtdump .lt. 0.0_re_p)
         write (iout,'(a)',advance='no')
     &         ' Enter Time between Dumps in Picoseconds',
     &         ' [0.1] :  '
         read (input,100,err=110)  dtdump
  100    format (f20.0)
         if (dtdump .le. 0.0_re_p)  dtdump = 0.1_re_p
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
            write (iout,130,advance='no')
  130       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ')
            read (input,140,err=150)  mode
  140       format (i10)
            if (mode .le. 0)  mode = 1
  150       continue
         end do
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0_re_p
            call nextarg (string,exist)
            if (exist)  read (string,*,err=170,end=170)  kelvin
  170       continue
            do while (kelvin .lt. 0.0_re_p)
               write (iout,180,advance='no')
  180          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ')
               read (input,190,err=200)  kelvin
  190          format (f20.0)
               if (kelvin .le. 0.0_re_p)  kelvin = 298.0_re_p
  200          continue
            end do
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0_re_p
            call nextarg (string,exist)
            if (exist)  read (string,*,err=210,end=210)  atmsph
  210       continue
            do while (atmsph .lt. 0.0_re_p)
               write (iout,'(a)',advance='no')
     &               ' Enter the Desired Pressure in Atm',
     &               ' [1.0] :  '
               read (input,230,err=240)  atmsph
  230          format (f20.0)
               if (atmsph .le. 0.0_re_p)  atmsph = 1.0_re_p
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
            write (iout,260,advance='no')
  260       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &               /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ')
            read (input,270,err=280)  mode
  270       format (i10)
            if (mode .le. 0)  mode = 1
  280       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0_re_p
            call nextarg (string,exist)
            if (exist)  read (string,*,err=290,end=290)  kelvin
  290       continue
            do while (kelvin .lt. 0.0_re_p)
               write (iout,300,advance='no')
  300          format (/,' Enter the Desired Temperature in Degrees',
     &                   ' K [298] :  ')
               read (input,310,err=320)  kelvin
  310          format (f20.0)
               if (kelvin .le. 0.0_re_p)  kelvin = 298.0_re_p
  320          continue
            end do
         end if
      end if
c
c     setup dynamics
c
      call mdinit(dt)
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
      time0 = 0
      do istep = 1, nstep
         call timer_enter( timer_timestep )
         if      (integrate .eq. 'VERLET') then
            call verlet (istep,dt)
         else if (integrate .eq. 'RESPA') then
            call respa (istep,dt)
         else if (integrate .eq. 'BBK') then
            call bbk(istep,dt)
         else if (integrate .eq. 'BAOAB') then
            call baoab(istep,dt)
         else if (integrate .eq. 'BAOABPISTON') then
            print*,"BAOABPISTON integrator unavailable for now"
            call fatal
            !call baoabpiston(istep,dt)
         else if (integrate .eq. 'BAOABRESPA') then
            call baoabrespa(istep,dt)
         else if (integrate .eq. 'BAOABRESPA1') then
            call baoabrespa1(istep,dt)
         else if (integrate .eq. 'RESPA1') then
            call respa1(istep,dt)
         else
            call beeman (istep,dt)
         end if
         call timer_exit( timer_timestep )

         ! Print mean simulation/execution time over last iprint timesteps
         if (mod(istep,iprint).eq.0) then
           time1 = timer_get_total( timer_timestep )
           timestep = time1 - time0
           time0 = time1
           if (rank.eq.0) then
              call MPI_REDUCE(MPI_IN_PLACE,timestep,1,MPI_REAL8,MPI_SUM,
     $          0,MPI_COMM_WORLD,ierr)
           else
              call MPI_REDUCE(timestep,timestep,1,MPI_REAL8,MPI_SUM,0,
     $          MPI_COMM_WORLD,ierr)
           end if
           if (verbose.and.rank.eq.0) then
              write(6,1000) iprint, (timestep)/nproc
              write(6,1010) 86400*dt*real(iprint*nproc,t_p)/
     $                                   (1000*timestep)
           end if
         end if

         ! Abort if any problem detected
         if (abort) call fatal
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
