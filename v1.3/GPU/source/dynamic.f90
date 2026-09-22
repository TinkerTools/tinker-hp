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
   use config, only: raw_cli_t, capture_raw_cli
   use domdec, only: rank
   use iounit, only: iout
   use mdstate, only: track_mds
   use status, only: status_t, failed
   use mdevice, only: bind_device
   implicit none
   integer ierr,nthreadsupport,world_rank
   type(raw_cli_t) raw_cli
   type(status_t) status_value
   call bind_device(status_value)
   if (failed(status_value)) stop 1
   call capture_raw_cli(raw_cli)
!      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,world_rank,ierr)
   iout = 6
   rank = world_rank
   track_mds = .false.
   call dynamic_bis(raw_cli,status_value)
   if (failed(status_value)) then
      if (world_rank.eq.0) write(iout,*) trim(status_value%message)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      stop 1
   end if
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)
end
!
subroutine dynamic_bis(raw_cli,status_value)
   use mlpot
   use atoms
   use bath
   use bound
   use config, only: raw_cli_t, config_t, parse_dynamic_cli
   use domdec
   use keys
   use inform
   use iounit
   use mdstuf
   use moldyn
   use mpi
   use potent  ,only: use_ml_embedding
   use parallel, only: parallel_t, release
   use tinheader, only: re_p
   use utils
   use mdevice, only: channel_t, finalize_device
   use initialization, only: dynamic_application, initialize
   use timestat
   use tinMemory
   use qtb, only: qtb_thermostat,adaptive_qtb
   use structure, only: structure_t
   use status, only: status_t, failed, enrich_status
   implicit none
   integer i,istep,nstep,ierr
   integer mode,next
   real(r_p) dt,dtdump
   real(r_p) timestep_fs,dump_interval_ps,temperature,pressure
   real(8) time0,time1,timestep
   character*20 keyword
   character*240 record
   character*240 string
   type(raw_cli_t), intent(in) :: raw_cli
   type(config_t) :: configuration
   type(structure_t) :: molecular_structure
   type(parallel_t) :: platform
   type(channel_t) :: reciprocal_work
   type(status_t), intent(out) :: status_value
   type(status_t) :: cleanup_status
!
!
1000 Format(' Time for ',I4,' Steps :',f12.4)
1010 Format(' ns per day',10x,':',f12.4)
   call initialize(dynamic_application, raw_cli, configuration, molecular_structure, &
      &platform, reciprocal_work, status_value)
   if (failed(status_value)) then
      return
   end if
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
!     parse Dynamic-specific command-line values without mutating raw argv
!
   call parse_dynamic_cli(raw_cli,use_bounds,nstep,timestep_fs,dump_interval_ps,&
      &mode,temperature,pressure,status_value)
   if (failed(status_value)) then
      call enrich_status(status_value,'Dynamic Command Line')
      final_message_suppressed = .true.
      call final
      call finalize_device(cleanup_status)
      call release(platform, cleanup_status)
      return
   end if
   dt = 0.001_re_p * real(timestep_fs,re_p)
   dtdump = real(dump_interval_ps,re_p)
   tautemp = max(tautemp,dt)
   taupres = max(taupres,dt)
   iwrite = nint(dtdump/dt)
   kelvin = real(temperature,re_p)
   atmsph = real(pressure,re_p)
   if (use_bounds) then
      isothermal = mode.eq.2 .or. mode.eq.4
      isobaric = mode.eq.4
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
   call final(status_value)
   call finalize_device(cleanup_status)
   call release(platform, cleanup_status)
   if (failed(cleanup_status)) status_value = cleanup_status
end
