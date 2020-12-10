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
#include "tinker_types.h"

      submodule(inform) subInform
      use domdec ,only: rank,MasterRank,COMM_TINKER
     &           ,glob,nloc
      use atomsMirror ,only: x,y,z,n
      use bath
      use cutoff
      use domdec ,only: nbloc
      use moldyn ,only: v,a
      use mpi
      use neigh  ,only: ineigup, lbuffer
      use potent
      use usage
      use vdw    ,only: skipvdw12
      use sizes  ,only: tinkerdebug

      contains

#include "convert.f.inc"

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
 15   format(A20,A15)
 16   format(A20,5x,L4)

      if (n.ne.nbloc) print 13, 'nbloc', nbloc
      print 13, 'nlupdate'     , ineigup
      print 14, 'list buffer'  , lbuffer
      if (use_vdw) then
         print 14, 'vdw cutoff'   , vdwcut
         print 14, 'vdw short cutoff'   , vdwshortcut
         print 14, 'vdw taper'    , vdwtaper
         print 14, 'shortheal'    , shortheal
         print 16, 'skipvdw12'    , skipvdw12
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
      print 15, 'thermostat', thermostat
      print 15, 'barostat', barostat
      end subroutine

#ifdef USE_DETERMINISTIC_REDUCTION
      module subroutine minmaxonef( vector,sz,name,mi_,ma_,on_ )
      implicit none
      integer sz
      mdyn_rtyp vector(*)
      character(*),optional,intent(in)::name
      real(8)     ,optional,intent(inout):: mi_,ma_,on_
      real(8) mi,ma,on
      integer i
      real(md_p) val
      mi=0;ma=0;on=0
!$acc parallel loop async present(vector(1:sz))
      do i = 1, sz
         val = mdr2md(vector(i))
         mi  = min( mi,val )
         ma  = max( ma,val )
         on  = on + abs( val )
      end do
!$acc wait
12    format(A6,3F16.6)
      write(*,12) name,mi,ma,on
      if (present(mi_)) mi_=mi
      if (present(ma_)) ma_=ma
      if (present(on_)) on_=on
      end subroutine
#endif

      module subroutine minmaxonet( vector,sz,name )
      implicit none
      integer sz
      real(t_p) vector(*)
      character(*),optional,intent(in)::name
      real(8) mi,ma,on
      integer i
      real(8) val
      mi=0;ma=0;on=0
!$acc parallel loop async present(vector(1:sz))
      do i = 1, sz
         val = (vector(i))
         mi  = min( mi,val )
         ma  = max( ma,val )
         on  = on + abs( val )
      end do
!$acc wait
12    format(A6,3F16.6)
      write(*,12) name,mi,ma,on
      end subroutine

#if TINKER_MIXED_PREC
      module subroutine minmaxoner( vector,sz,name )
      implicit none
      integer sz
      real(r_p) vector(*)
      character(*),optional,intent(in)::name
      real(8) mi,ma,on
      integer i
      real(r_p) val
      mi=0;ma=0;on=0
!$acc parallel loop async present(vector(1:sz))
      do i = 1, sz
         val = vector(i)
         mi  = min( mi,val )
         ma  = max( ma,val )
         on  = on + abs( val )
      end do
!$acc wait
12    format(A6,3F16.6)
      write(*,12) name,mi,ma,on
      end subroutine
#endif
c
c     Print information on position, velocities and aceleration
c
      module subroutine info_minmax_pva()
      implicit none
      integer i,iglob
      integer,parameter::nel=5
      real(r_p) minmax(3*nel)
      logical tinker_isnan

      minmax=0
!$acc wait
 20   format (80('_'))
!$acc data copy(minmax)

!$acc parallel loop default(present)
      do i = 1, nloc
         iglob      = glob(i)
         if (use(iglob)) then
!$acc atomic update
         minmax(01) = min(minmax(01),x(iglob))
!$acc atomic update
         minmax(02) = min(minmax(02),y(iglob))
!$acc atomic update
         minmax(03) = min(minmax(03),z(iglob))
!$acc atomic update
         minmax(04) = min(minmax(04),v(1,iglob))
!$acc atomic update
         minmax(04) = min(minmax(04),v(2,iglob))
!$acc atomic update
         minmax(04) = min(minmax(04),v(3,iglob))
!$acc atomic update
         minmax(05) = min(minmax(05),a(1,iglob))
!$acc atomic update
         minmax(05) = min(minmax(05),a(2,iglob))
!$acc atomic update
         minmax(05) = min(minmax(05),a(3,iglob))
!$acc atomic update
         minmax(06) = max(minmax(06),x(iglob))
!$acc atomic update
         minmax(07) = max(minmax(07),y(iglob))
!$acc atomic update
         minmax(08) = max(minmax(08),z(iglob))
!$acc atomic update
         minmax(09) = max(minmax(09),v(1,iglob))
!$acc atomic update
         minmax(09) = max(minmax(09),v(2,iglob))
!$acc atomic update
         minmax(09) = max(minmax(09),v(3,iglob))
!$acc atomic update
         minmax(10) = max(minmax(10),a(1,iglob))
!$acc atomic update
         minmax(10) = max(minmax(10),a(2,iglob))
!$acc atomic update
         minmax(10) = max(minmax(10),a(3,iglob))
!$acc atomic update
         minmax(11) = minmax(11)+v(1,iglob)
!$acc atomic update
         minmax(11) = minmax(11)+v(2,iglob)
!$acc atomic update
         minmax(11) = minmax(11)+v(3,iglob)
!$acc atomic update
         minmax(12) = minmax(12)+a(1,iglob)
!$acc atomic update
         minmax(12) = minmax(12)+a(2,iglob)
!$acc atomic update
         minmax(12) = minmax(12)+a(3,iglob)
         end if
      end do

!$acc end data
      if (rank.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,minmax,nel,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,minmax(6),nel,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,minmax(11),2,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      else
         call MPI_REDUCE(minmax,minmax,nel,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(minmax(6),minmax(6),nel,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(minmax(11),minmax(11),2,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      end if

      if (rank.eq.0) then
30    format(A10,2F18.10)
32    format(A10,3F18.10)
         print 20
         print 30,"min max_x ",minmax(01),minmax(06)
         print 30,"min max_y ",minmax(02),minmax(07)
         print 30,"min max_z ",minmax(03),minmax(08)
         print 32,"min max_v ",minmax(04),minmax(09),minmax(11)
         print 32,"min max_a ",minmax(05),minmax(10),minmax(12)
      end if
      end subroutine

      end submodule
