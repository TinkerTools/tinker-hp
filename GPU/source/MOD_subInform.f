c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module subInform  --  submodule for Inform                  ##
c     ##                                                              ##
c     ##################################################################
c
#include "tinker_precision.h"

      submodule(inform) subInform
      use domdec ,only: rank,MasterRank,COMM_TINKER
     &           ,glob,nloc
      use atomsMirror ,only: x,y,z,n
      use cutoff
      use moldyn ,only: v,a
      use mpi
      use neigh  ,only: ineigup, lbuffer
      use potent
      use usage
      use sizes  ,only: tinkerdebug

      contains

      module subroutine initDebugEnv()
      implicit none
      integer ierr,length
      character*32 dvalue
      integer,parameter::success=0

      ! Default debug switch value
      deb_Path    = .false.
      deb_Force   = .false.
      deb_Energy  = .false.
      deb_Polar   = .false.
      tinkerdebug = 0

      ! Fetch is possible TINKER_DEBUG From environment
      call get_environment_variable("TINKER_DEBUG",dvalue,length,
     &         status=ierr)

      if (ierr.eq.success) read(dvalue,*) tinkerdebug

      ! Set debug switches
      if (tinkerdebug>0) then
         !if (rank.eq.0) print*, 'debug enable',tinkerdebug
         if (btest(tinkerdebug,tindPath)) then
            if (rank.eq.MasterRank) deb_Path = .true.
         end if
         if (btest(tinkerdebug,tindForce )) deb_Force  = .true.
         if (btest(tinkerdebug,tindEnergy)) deb_Energy = .true.
         if (btest(tinkerdebug,tindAtom  )) deb_Atom   = .true.
         if (btest(tinkerdebug,tindPolar )) deb_Polar  = .true.
      end if

      end subroutine

      module subroutine info_dyn()
      implicit none
      integer i

 13   format(A20,I10)
 14   format(A20,F15.6)

      print 13, 'nlupdate'     , ineigup
      print 14, 'list buffer'  , lbuffer
      if (use_vdw) then
         print 14, 'vdw cutoff'   , vdwcut
         print 14, 'vdw short cutoff'   , vdwshortcut
         print 14, 'vdw taper'    , vdwtaper
         print 14, 'shortheal'    , shortheal
      end if
      if (use_mpole.or.use_polar) then
         print 14, 'mpole cutoff' , mpolecut
         print 14, 'mpole short cutoff' , mpoleshortcut
         print 14, 'ewaldcut'     , ewaldcut
      end if
      if (use_charge) then
         print 14, 'charge cutoff', mpolecut
         print 14, 'chg taper'    , chgtaper
      end if
      end subroutine
c
c     Print information on position, velocities and aceleration
c
      module subroutine info_minmax_pva()
      implicit none
      integer i,iglob
      integer,parameter::nel=5
      real(r_p) minmax(2*nel)
      logical tinker_isnan

!$acc wait
!$acc update host(x,y,z,v,a)
 20   format (80('_'))

      do i = 1, nloc
         iglob      = glob(i)
         if (use(iglob)) then
         minmax(01) = min(minmax(01),x(iglob))
         minmax(02) = min(minmax(02),y(iglob))
         minmax(03) = min(minmax(03),z(iglob))
         minmax(04) = min(minmax(04),v(1,iglob))
         minmax(04) = min(minmax(04),v(2,iglob))
         minmax(04) = min(minmax(04),v(3,iglob))
         minmax(05) = min(minmax(05),a(1,iglob))
         minmax(05) = min(minmax(05),a(2,iglob))
         minmax(05) = min(minmax(05),a(3,iglob))
         minmax(06) = max(minmax(06),x(iglob))
         minmax(07) = max(minmax(07),y(iglob))
         minmax(08) = max(minmax(08),z(iglob))
         minmax(09) = max(minmax(09),v(1,iglob))
         minmax(09) = max(minmax(09),v(2,iglob))
         minmax(09) = max(minmax(09),v(3,iglob))
         minmax(10) = max(minmax(10),a(1,iglob))
         minmax(10) = max(minmax(10),a(2,iglob))
         minmax(10) = max(minmax(10),a(3,iglob))
         end if
      end do

      if (rank.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,minmax,nel,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,minmax(6),nel,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
      else
         call MPI_REDUCE(minmax,minmax,nel,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(minmax(6),minmax(6),nel,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
      end if

      if (rank.eq.0) then
         print 20
         print*,"min max_x ", minmax(01),minmax(06)
         print*,"min max_y ", minmax(02),minmax(07)
         print*,"min max_z ", minmax(03),minmax(08)
         print*,"min max_v ", minmax(04),minmax(09)
         print*,"min max_a ", minmax(05),minmax(10)
      end if
      end subroutine

      end submodule
