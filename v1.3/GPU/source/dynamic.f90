!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  program dynamic  --  run molecular or stochastic dynamics  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "dynamic" computes a molecular dynamics trajectory
!     in one of the standard statistical mechanical ensembles and using
!     any of several possible integration methods
!
!
#include "tinker_macro.h"
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
!      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
   call MPI_INIT(ierr)
   call dynamic_bis
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)
end
!
subroutine dynamic_bis
   use mlpot
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
   use potent  ,only: use_ml_embedding
   use utils
   use utilgpu ,only: rec_queue,ti_p,re_p
   use timestat
   use tinMemory
   use qtb, only: qtb_thermostat,adaptive_qtb
   implicit none
   integer i,istep,nstep,ierr
   integer mode,next
   real(r_p) dt,dtdump
   real(8) time0,time1,timestep
   logical exist,query
   character*20 keyword
   character*240 record
   character*240 string
!
!
1000 Format(' Time for ',I4,' Steps :',f12.4)
1010 Format(' ns per day',10x,':',f12.4)
   ! Sign running program
   app_id = dynamic_a
!
!     set up the structure and molecular mechanics calculation
!
   call initial
   call init_keys
   call initmpi
   call getxyz
   call unitcell
   call cutoffs
   call lattice
!
!     setup for MPI
!
   call drivermpi
   !call kewald_2  ! FIXME  Is it necessary
   call reinitnl(0)
!
!     allocate some arrays
!
   call prmem_requestm(    v,3,n)
   call prmem_requestm(    a,3,n)
   call prmem_requestm( aalt,3,n)
   call prmem_requestm(aalt2,3,n)
!
   call set_to_zero2m(   a,    v,3*n,rec_queue)
   call set_to_zero2m(aalt,aalt2,3*n,rec_queue)
!
   call mechanic
   call nblist(0)
!
!     initialize the temperature, pressure and coupling baths
!
   kelvin = 0.0_re_p
   atmsph = 0.0_re_p
   isothermal = .false.
   isobaric = .false.
!
!     check for keywords containing any altered parameters
!
   integrate = 'VERLET'
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
!
!     initialize the simulation length as number of time steps
!
   query = .true.
   call nextarg (string,exist)
   if (.not.exist)  then
     if (ranktot.eq.0) write (iout,*) 'You Need To Enter the Number of Steps'
     call MPI_BARRIER(COMM_TINKER,ierr)
     __TINKER_FATAL__
   else
!   if (exist) then
      read (string,*,err=10,end=10)  nstep
      query = .false.
   end if
10 continue
!
!     get the length of the dynamics time step in picoseconds
!
   dt = -1.0_re_p
   call nextarg (string,exist)
   if (.not.exist)  then
     if (ranktot.eq.0) write (iout,*) 'You Need To Enter the Time Step Length in Femtoseconds'
     call MPI_BARRIER(COMM_TINKER,ierr)
     __TINKER_FATAL__
   else
     read (string,*,err=40,end=40)  dt
   end if
40 continue
   dt = 0.001_re_p * dt
!
!     enforce bounds on thermostat and barostat coupling times
!
   tautemp = max(tautemp,dt)
   taupres = max(taupres,dt)
!
!     set the time between trajectory snapshot coordinate dumps
!
   dtdump = -1.0_re_p
   call nextarg (string,exist)
   if (.not.exist) then
     if (ranktot.eq.0) write (iout,*) 'You Need To Enter the Time Between Dumps in Picoseconds'
     call MPI_BARRIER(COMM_TINKER,ierr)
     __TINKER_FATAL__
   else
     read (string,*,err=80,end=80)  dtdump
   end if
80 continue
   iwrite = nint(dtdump/dt)
!
!     get choice of statistical ensemble for periodic system
!
   if (use_bounds) then
      mode = -1
      call nextarg (string,exist)
      if (.not.exist) then
        if (ranktot.eq.0) write (iout,*) 'You Need To Enter the Statistical Ensemble:'
        if (ranktot.eq.0) write (iout,130)
130      format (/,' Available Statistical Mechanical Ensembles :',&
         &//,4x,'(1) Microcanonical (NVE)',&
         &/ ,4x,'(2) Canonical (NVT)',&
         &/ ,4x,'(4) Isothermal-Isobaric (NPT)')
        call MPI_BARRIER(COMM_TINKER,ierr)
        __TINKER_FATAL__
      else
        read (string,*,err=120,end=120)  mode
      end if
120   continue
      if (mode.le.0)  mode= 1
      if (mode.eq.2 .or. mode.eq.4) then
         isothermal = .true.
         kelvin = -1.0d0
         call nextarg (string,exist)
         if (.not.exist) then
            if (ranktot.eq.0) write (iout,*) 'You Need To Enter the Temperature in Kelvin'
            call MPI_BARRIER(COMM_TINKER,ierr)
            __TINKER_FATAL__
         else
            read (string,*,err=170,end=170)  kelvin
         end if
170      continue
         if (kelvin .le. 0.0d0)  kelvin = 298.0d0
      end if
      if (mode.eq.4) then
         isobaric = .true.
         atmsph = -1.0_re_p
         call nextarg (string,exist)
         if (.not.exist) then
            if (ranktot.eq.0) write (iout,*) 'You Need To Enter the Pressure in Atm'
            call MPI_BARRIER(COMM_TINKER,ierr)
            __TINKER_FATAL__
         else 
            read (string,*,err=210,end=210)  atmsph
         end if
210      continue
         if (atmsph .le. 0.0_re_p)  atmsph = 1.0_re_p
      end if
   end if
!
!     setup dynamics
!
   call mdinit(dt)
!
!     print out a header line for the dynamics computation
!
   if (integrate .eq. 'VERLET') then
      if (rank.eq.0) write (iout,330)
330   format (/,' Molecular Dynamics Trajectory via',&
         &' Velocity Verlet Algorithm')
   else if (integrate .eq. 'RESPA') then
      if (rank.eq.0) write (iout,390)
390   format (/,' Molecular Dynamics Trajectory via',&
         &' r-RESPA MTS Algorithm')
   else if (integrate .eq. 'BBK') then
      if (rank.eq.0) write (iout,400)
400   format (/,' Langevin Molecular Dynamics Trajectory via',&
         &' BBK Algorithm')
   else if (integrate .eq. 'BAOAB') then
      if (rank.eq.0) write (iout,410)
410   format (/,' Langevin Molecular Dynamics Trajectory via',&
         &' BAOAB Algorithm')
   else if (integrate .eq. 'BAOABRESPA') then
      if (rank.eq.0) write (iout,420)
420   format (/,' Langevin Molecular Dynamics Trajectory via',&
         &' BAOAB-RESPA Algorithm')
   else if (integrate .eq. 'BAOABRESPA1') then
      if (rank.eq.0) write (iout,430)
430   format (/,' Langevin Molecular Dynamics Trajectory via',&
         &' BAOAB-RESPA-1 Algorithm')
   else if (integrate .eq. 'RESPA1') then
      if (rank.eq.0) write (iout,440)
440   format (/,' Molecular Dynamics Trajectory via',&
         &' r-RESPA-1 MTS Algorithm')
   else
      if (rank.eq.0) write (iout,470)
470   format (/,' Molecular Dynamics Trajectory via',&
         &' Modified Beeman Algorithm')
   end if

   if(qtb_thermostat .and. rank==0) then
      if(adaptive_qtb) then
         write(iout,*) " using adaptive QTB thermostat"
      else
         write(iout,*) " using standard QTB thermostat"
      endif
   endif
!
!     Inform about ML Embedding
!
   if (use_ml_embedding) then
      if(ranktot==0) write(*,*) 'USING ML EMBEDDING'
      ml_embedding_mode=2
   endif
!
!     integrate equations of motion to take a time step
!
   time0 = 0
   do istep = 1, nstep
      step_c = istep
      call timer_enter( timer_timestep )
      if      (integrate .eq. 'VERLET') then
         call verlet (istep,dt)
      else if (integrate .eq. 'RESPA') then
         call respa (istep,dt)
      else if (integrate .eq. 'BBK') then
         call bbk(istep,dt)
      else if (integrate .eq. 'BAOAB') then
         call baoab(istep,dt)
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
            call MPI_REDUCE(MPI_IN_PLACE,timestep,1,MPI_REAL8,MPI_SUM,&
               &0,COMM_TINKER,ierr)
         else
            call MPI_REDUCE(timestep,timestep,1,MPI_REAL8,MPI_SUM,0,&
               &COMM_TINKER,ierr)
         end if
         if (verbose.and.rank.eq.0) then
            write(6,1000) iprint, (timestep)/nproc
            write(6,1010) 86400*dt*real(iprint*nproc,t_p)/&
               &(1000*timestep)
         end if
      end if

      ! Abort if any problem detected
      if (abort) __TINKER_FATAL__
   end do
!
!     perform any final tasks before program exit
!
   call final
end
