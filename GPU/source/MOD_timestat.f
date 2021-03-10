c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c   module timestat : Timers to measure time spent in regions of the code
c
c     tinkertime     logical flag to enable timers verbosity

#include "tinker_precision.h"
      module timestat
        use domdec, only: rank,nproc,nproctot,COMM_TINKER
        use mpi   , only: MPI_Wtime, MPI_COMM_WORLD, MPI_REAL8
     &            , MPI_SUM, MPI_MAX, MPI_MIN, MPI_ALLREDUCE
     &            , MPI_REDUCE, MPI_IN_PLACE
        implicit none
        ! Timers Parameters
        integer tinkertime
        integer, parameter, private:: numslots= 2
        integer, parameter, private:: name_len= 25
        logical, parameter :: quiet_timers = .true.
        logical, parameter :: print_timers = .false.
        enum,bind(c)
          enumerator normal_disp
          enumerator   stat_disp
          enumerator    ave_disp
        end enum
        enum,bind(C)
          enumerator sumy_time
          enumerator path_time
          enumerator sync_time
          enumerator all_time=8
        end enum
        ! Timers ids used in code
        enum, bind(C)
        enumerator ::
     &           timer_ebond=1,
     &           timer_eangle,
     &           timer_eurey1,timer_eurey3,
     &           timer_eopbend1,timer_eopbend3,
     &           timer_geom,
     &           timer_nl,
     &           timer_clistcell,timer_vlistcell,timer_mlistcell,
     &           timer_vdw,
     &           timer_ehal3,timer_elj3,
     &           timer_elec,
     &           timer_echarge,timer_ecreal,timer_ecrecip,
     &           timer_empole3,timer_emreal,timer_emrecip,
     &           timer_efld0_direct,timer_efld0_recip,
     &           timer_tmatxb_pmevec,timer_tmatxb_pmegpu,
     &           timer_tmatxrecipgpu,
     &           timer_polar,timer_polarsolve,
     &           timer_epreal,timer_eprecip,
     &           timer_real,timer_realdip,timer_rec,timer_recdip,
     &           timer_grid1,timer_grid2,timer_ffts,timer_scalar,
     &           timer_torque,
     &           timer_bonded,timer_nonbonded,
     &           timer_fcomm,timer_dirreccomm,timer_rectorcomm,
     &           timer_dirbondfcomm,timer_recreccomm,
     &           timer_polsolvcomm,timer_polfieldcomm,
     &           timer_timestep,timer_eneig,timer_commstep,
     &           timer_param,timer_other,timer_clear,
     &           timer_debug,timer_scalfactor,
     &           timer_ulspred,timer_fmanage,timer_tinker,
     &           timer_io,timer_prog,timer_plumed,timer_reduceen,
     &           timer_b1,timer_b2,timer_b3,
     &           last_timer
        end enum
        ! number Timers
        integer, parameter:: max_timer_count=last_timer

        integer :: stat_timers(30) =
     &          [timer_param,timer_nl,
     &           timer_bonded,timer_vdw,
     &           timer_real,timer_realdip,timer_rec,timer_recdip,
     &           timer_grid1,timer_grid2,timer_scalar,timer_ffts,
     &           timer_ulspred,timer_fmanage,timer_plumed,
     &           timer_nonbonded,timer_torque,timer_other,timer_clear,
     &           timer_fcomm,timer_dirbondfcomm,timer_dirreccomm,
     &           timer_recreccomm,timer_commstep,timer_polsolvcomm,
     &           timer_polfieldcomm,timer_eneig,timer_reduceen,
     &           timer_timestep,timer_b1
     &          ]
        integer :: comm_timers(8) =
     &          [timer_param,timer_fcomm,timer_recreccomm,
     &           timer_polsolvcomm,timer_polfieldcomm,
     &           timer_eneig,timer_commstep,timer_reduceen
     &          ]
        integer :: comp_timers(14) =
     &          [timer_fmanage,timer_bonded,timer_vdw,timer_nl,
     &           timer_real,timer_realdip,timer_grid1,timer_grid2,
     &           timer_scalar,timer_ffts,timer_clear,
     &           timer_other,timer_ulspred,timer_plumed
     &          ]

        ! Cumulated times on regions
        real*8 timeclear,timereneig,timecommstep,timeparam
        real*8 timeforcescomm,timedirreccomm
        real*8 timereciptorquescomm
        real*8 timedirbondfcomm,timebonded,timenonbonded
        real*8 timenl
        real*8 timereal,timerec,timerecreccomm
        real*8 timerealdip,timerecdip,timepolarsolve
        real*8 timegrid1,timeffts,timescalar,timegrid2

        ! Timer names
        character(name_len),protected :: timer_name(max_timer_count)
        ! Time measurement
        real*8, protected:: timer_time     (max_timer_count)
        real*8, protected:: timer_starttime(max_timer_count)
        real*8, protected:: timer_savetime (max_timer_count,numslots)
        real*8, protected:: timer_lasttime (max_timer_count)

        ! Timer subroutines
        public :: initiate_timers,
     &      timer_enter, timer_exit, timer_save,
     &      timer_get_total, timer_get_last, timer_get_save,
     &      display_timers
      contains

      ! Init module timestat
      subroutine init_timestat()
        timer_time(:)       = -1
        timer_starttime(:)  = -2
        timer_lasttime(:)   = 0
        timer_savetime(:,:) = 0
      end subroutine

      ! Create new timer identifier
      subroutine timer_init( name, id )
        character(*), intent(in) :: name
        integer     , intent(in) ::  id

        timer_name(id) = name
        timer_time(id) = 0
      end subroutine

      subroutine initiate_timers
        implicit none
        integer ierr
        character*32 timeval
        integer,parameter::sucess=0

         ! Fetch if possible "TINKER_TIME" From Environment
        call get_environment_variable('TINKER_TIME',timeval,
     &       status=ierr)

        if (ierr.eq.sucess) then
           read(timeval,*) tinkertime
        else
           tinkertime = 0
           if (nproctot.gt.1) tinkertime=1
        end if

        call init_timestat()
        call timer_init( "bonded"          , timer_bonded)
        call timer_init( "ebond"           , timer_ebond)
        call timer_init( "eangle"          , timer_eangle)
        call timer_init( "eopbend1"        , timer_eopbend1)
        call timer_init( "eopbend3"        , timer_eopbend3)
        call timer_init( "eurey1"          , timer_eurey1)
        call timer_init( "eurey3"          , timer_eurey3)
        call timer_init( "vdw"             , timer_vdw)
        call timer_init( "elj"             , timer_elj3)
        call timer_init( "ehal3"           , timer_ehal3)
        call timer_init( "elec"            , timer_elec)
        call timer_init( "echarge"         , timer_echarge)
        call timer_init( "ecreal"          , timer_ecreal)
        call timer_init( "ecrecip"         , timer_ecrecip)
        call timer_init( "empole3"         , timer_empole3)
        call timer_init( "emreal1cgpu"     , timer_emreal)
        call timer_init( "emrecip1gpu"     , timer_emrecip)
        call timer_init( "geom"            , timer_geom)
        call timer_init( "neighbor list"   , timer_nl)
        call timer_init( "clistcell"       , timer_clistcell)
        call timer_init( "vlistcell"       , timer_vlistcell)
        call timer_init( "mlistcell"       , timer_mlistcell)
        call timer_init( "nonbonded"       , timer_nonbonded)
        call timer_init( "polar"           , timer_polar)
        call timer_init( "efld0_direct"    , timer_efld0_direct)
        call timer_init( "efld0_recip"     , timer_efld0_recip)
        call timer_init( "tmatxb_pmevec"   , timer_tmatxb_pmevec)
        call timer_init( "tmatxb_pmegpu"   , timer_tmatxb_pmegpu)
        call timer_init( "tmatxrecipgpu"   , timer_tmatxrecipgpu)
        call timer_init( "epreal1c"        , timer_epreal)
        call timer_init( "eprecip1gpu"     , timer_eprecip)
        call timer_init( "polar solv  comm", timer_polsolvcomm)
        call timer_init( "polar field comm", timer_polfieldcomm)
        call timer_init( "polar solver"    , timer_polarsolve)
        call timer_init( "real space"      , timer_real)
        call timer_init( "rec space"       , timer_rec)
        call timer_init( "realdip space"   , timer_realdip)
        call timer_init( "recdip space"    , timer_recdip)
        call timer_init( "ulspred"         , timer_ulspred)
        call timer_init( "fill grid"       , timer_grid1)
        call timer_init( "extract grid"    , timer_grid2)
        call timer_init( "ffts"            , timer_ffts)
        call timer_init( "scalar product"  , timer_scalar)
        call timer_init( "torquegpu"       , timer_torque)
        call timer_init( "timestep"        , timer_timestep)
        call timer_init( "reassign neig"   , timer_eneig)
        call timer_init( "scalfactor"      , timer_scalfactor)
        call timer_init( "position    comm", timer_commstep)
        call timer_init( "param"           , timer_param)
        call timer_init( "clear force"     , timer_clear)
        call timer_init( "other"           , timer_other)
        call timer_init( "forces set/add"  , timer_fmanage)
        call timer_init( "forces      comm", timer_fcomm)
        call timer_init( "dirrec      comm", timer_dirreccomm)
        call timer_init( "rectorque   comm", timer_rectorcomm)
        call timer_init( "dirbondf    comm", timer_dirbondfcomm)
        call timer_init( "recrec      comm", timer_recreccomm)
        call timer_init( "tinker"          , timer_tinker)
        call timer_init( "io"              , timer_io)
        call timer_init( "program"         , timer_prog)
        call timer_init( "reduceen"        , timer_reduceen)
        call timer_init( "plumed"          , timer_plumed)
        call timer_init( "buff1"           , timer_b1)
        call timer_init( "buff2"           , timer_b2)
        call timer_init( "buff3"           , timer_b3)
        call timer_init( "last"            , last_timer)

        timeclear      = 0.0
        timereneig     = 0.0
        timecommstep   = 0.0
        timeparam      = 0.0
        timeforcescomm = 0.0
        timedirreccomm = 0.0
        timedirbondfcomm  = 0.0
        timebonded     = 0.0
        timenonbonded  = 0.0
        timepolarsolve = 0.0
        timereal       = 0.0
        timerealdip    = 0.0
        timegrid1      = 0.0
        timeffts       = 0.0
        timescalar     = 0.0
        timegrid2      = 0.0
        timerecreccomm = 0.0
        timerec        = 0.0
        timerecdip     = 0.0
      end subroutine

      subroutine verify_timer_id( id )
        integer, intent(in) ::  id
        if( id < 1 .or. id > max_timer_count ) then
          print *, "WARNING : timer id out of range : ", id
          return
        end if
      end subroutine

      ! Start timer, enter region
      subroutine timer_enter( id )
        use nvtx
        integer, intent(in) ::  id

        call verify_timer_id( id )
        if( timer_starttime(id) > 0 )
     &    print *, "WARNING : timer started but not stopped : ",
     &             timer_name(id), id
        call nvtxStartRange( timer_name(id) )
        if (btest(tinkertime,sync_time)) then
!$acc wait
        end if
        timer_starttime(id) = MPI_Wtime()
      end subroutine

      ! End timer, exit region, print timer if TINKER_TIME is on
      subroutine timer_exit( id, quiet_ )
        use nvtx
        integer, intent(in) ::  id
        logical, optional :: quiet_ ! Never print time
        real*8 :: time_diff
        logical:: quiet

        quiet = .not.btest(tinkertime,path_time)
        if(present(quiet_)) quiet = quiet_

        call verify_timer_id( id )
        if( timer_starttime(id) < 0 )
     &     print *, "WARNING : timer stopped but not started : ",
     &              timer_name(id), id
        call nvtxEndRange()
        if (btest(tinkertime,sync_time)) then
!$acc wait
        end if
        time_diff = MPI_Wtime() - timer_starttime(id)
        timer_starttime(id) = -1 ! for timer started but not stopped verification
        timer_lasttime(id) = time_diff
        timer_time(id) = timer_time(id) + time_diff
        ! Print on process 0 if TINKER_TIME env variable is 1
        if (rank.eq.0 .and. .not. quiet) then
          write(*,'(A,A18,A,F10.6,A,F13.6)') 'timer: ', timer_name(id),
     &                   'last:', time_diff ,' total:', timer_time(id)
        end if
      end subroutine

      ! Save timer
      subroutine timer_save( id,slot )
        implicit none
        integer,intent(in),optional:: id(:)
        integer,intent(in),optional:: slot
        integer sz, slot_

        if (present(slot)) then
           slot_ = slot
        else
           slot_ = 1
        end if

        if (present(id)) then
           timer_savetime(id,slot_) = timer_time(id)
        else
           timer_savetime(: ,slot_) = timer_time(: )
        end if
      end subroutine

      ! Reset timers
      subroutine reset_timer( id )
        integer, intent(in) ::  id
        call verify_timer_id( id )
        timer_time(id) = -1
        timer_starttime(id) = -2
        timer_lasttime(id) = 0
      end subroutine

      ! Get cumulated time
      real(8) function timer_get_total( id ) result( time )
        integer, intent(in) ::  id
        call verify_timer_id( id )
        time = timer_time( id )
      end function

      ! Get last time
      real(8) function timer_get_last( id ) result( time )
        integer, intent(in) ::  id
        call verify_timer_id( id )
        time = timer_lasttime( id )
      end function

      ! Get save time
      real(8) function timer_get_save( id, slot ) result( time )
        integer, intent(in) :: id
        integer, optional   :: slot

        call verify_timer_id( id )
        if (present(slot)) then
           time = timer_time( id ) - timer_savetime( id,slot )
        else
           time = timer_time( id ) - timer_savetime( id,1    )
        end if
      end function

      !display subroutine time of one rank
      subroutine display_timers( ids,config,slot,iter )
         implicit none
         integer,intent(in),optional :: ids(:)
         integer,intent(in),optional :: config,slot,iter
         real(8) buffer(max_timer_count,3)
         real(8) t0,t1,t_timestep,t_comm,t_comp,t_comm1,t_comm2
         real(8) tcomp1,tcomp2
         integer i,j,ierr,cfg,slt,it

   20    format(59('-'))
   30    format(/,25('-'),A,25('-'))

         cfg = stat_disp
         slt = 1
         it  = 1

         if (present(config)) cfg = config
         if (present(slot  )) slt = slot
         if (present(iter  )) it  = iter

         buffer(:,3) = timer_time(:) - timer_savetime(:,slt)
         if (nproc.Eq.1) then
            buffer(:,1) = buffer(:,3)
            buffer(:,2) = buffer(:,3)
         else
            call MPI_REDUCE(buffer(1,3),buffer(1,1),max_timer_count,
     &           MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(buffer(1,3),buffer(1,2),max_timer_count,
     &           MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
         end if

         if (cfg.eq.stat_disp) then
         else if (cfg.eq.ave_disp) then
            buffer = buffer / iter
         end if
         t_timestep = buffer(timer_timestep,1)

         if (rank.eq.0) write(*,30) ' Timers '

         if (present(ids)) then
            do j = 1,size(ids)
               i = ids(j)
               t0 = buffer(i,1)
               t1 = buffer(i,2)
               if (cfg.eq.ave_disp) then
                  t1 = 100.0*t0/(t_timestep)
               end if
               if ( rank.eq.0) then
               call timer_print( timer_name(i),t0,t1,cfg )
               end if
            end do
         else
            do i = 1,max_timer_count
               t0 = buffer(i,1)
               t1 = buffer(i,2)
               if (cfg.eq.ave_disp) then
                  t1 = 100.0*t0/(t_timestep)
               end if
               if ( rank.eq.0 ) then
               call timer_print( timer_name(i),t0,t1,cfg )
               end if
            end do
         end if

         t_comm = 0
         t_comp = 0
         do i = 1,size(comm_timers)
            t_comm = t_comm + buffer(comm_timers(i),3)
c           write(*,'(A20,2F14.4)') timer_name(comm_timers(i)),
c    &            buffer(comm_timers(i),3),t_comm
         end do

         call MPI_REDUCE(t_comm,t_comm1,1,MPI_REAL8,MPI_MAX,0
     &                  ,COMM_TINKER,ierr)
         call MPI_REDUCE(t_comm,t_comm2,1,MPI_REAL8,MPI_MIN,0
     &                  ,COMM_TINKER,ierr)

         do i = 1,size(comp_timers)
            t_comp = t_comp + buffer(comp_timers(i),3)
c           if (rank.eq.0)
c    &      write(*,'(A20,2F14.4)') timer_name(comp_timers(i)),
c    &            buffer(comp_timers(i),3),t_comp
         end do
         call MPI_REDUCE(t_comp,tcomp2,1,MPI_REAL8,MPI_MAX,0
     &                  ,COMM_TINKER,ierr)
         call MPI_REDUCE(t_comp,tcomp1,1,MPI_REAL8,MPI_MIN,0
     &                  ,COMM_TINKER,ierr)

         if (cfg.eq.ave_disp.and.rank.eq.0) then
            print*
            if (nproc.ne.1) then
            call timer_print( " MPI Slow-rk",t_comm1,
     &           100.0*t_comm1/t_timestep,ave_disp)
            call timer_print( " MPI Fast-rk",t_comm2,
     &           100.0*t_comm2/t_timestep,ave_disp)
            end if
            
            if (nproc.ne.1) then
            call timer_print( "Compute Slow-rk",tcomp1,
     &           100.0*tcomp1/t_timestep,ave_disp)
            call timer_print( "Compute Fast-rk",tcomp2,
     &           100.0*tcomp2/t_timestep,ave_disp)
            end if
         end if

         if (rank.eq.0) write(*,20)

      end subroutine

      subroutine timer_print( name,time0,time1,config )
      implicit none
      character(*) name
      real(8),intent(in):: time0
      real(8),intent(in),optional:: time1
      integer,intent(in),optional:: config
      integer cfg
      real(8) t1
      real(8),parameter:: min_duration=1e-5

  10  format(2x,a16,f12.6)
  11  format(2x,a16,' min',f12.6,' max',f12.6)
  12  format('Ave Time for ',a16,1x,f12.6,' sec',f7.2,' %')
  13  format('timer: ',a16,'last:',F10.6,' total',F13.6)

      if (time0.lt.1e-5) return

      t1  = 0
      cfg = normal_disp

      if (present(config)) cfg = config
      if (present(time1))  t1  = time1

      if (cfg.eq.normal_disp) then
         write(*,13) name,time0,t1
      else if (cfg.eq.stat_disp.and.rank.eq.0) then
         if (nproc.eq.1) then
            write(*,10) name,time0
         else
            write(*,11) name,time0,t1
         end if
      else if (cfg.eq.ave_disp) then
         write(*,12) name,time0,t1
      end if

      end subroutine

      end module
