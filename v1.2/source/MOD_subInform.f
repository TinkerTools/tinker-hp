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
      submodule(inform) subInform
      use atoms ,only: x,y,z,n
      use bath
      use cutoff
      use domdec ,only: rank,ranktot,COMM_TINKER
     &           ,glob,nloc,nlocnl,nbloc,nproc
      use moldyn ,only: v,a,aalt,aalt2
      use mpi
      use neigh  ,only: ineigup, lbuffer
      use potent
      use polpot
      use usage  ,only: muse=>use
      use sizes  ,only: tinkerdebug

      contains

      module subroutine initDebugEnv()
      implicit none
      integer ierr,length,sav
      character*32 dvalue
      integer,parameter::success=0

      ! Default debug switch value
      deb_Path    = .false.
      deb_Force   = .false.
      deb_Energy  = .false.
      deb_Polar   = .false.
      tinkerdebug = 0
      sav         = 0
      dd_verbose  = .true.

      ! Fetch if possible TINKER_DEBUG From environment
      call get_environment_variable("TINKER_DEBUG",dvalue,length,
     &         status=ierr)

      if (ierr.eq.success) read(dvalue,'(I16)') tinkerdebug

      sav = abs(tinkerdebug)
      if (tinkerdebug.lt.0) tinkerdebug=ior(ishft(1,31),sav)


      ! Set debug switches
      if (sav.gt.0) then
         if (btest(sav,tindPath)) then
            if (ranktot.eq.0) deb_Path = .true.
         end if
         if (btest(sav,tindForce )) deb_Force  = .true.
         if (btest(sav,tindEnergy)) deb_Energy = .true.
         if (btest(sav,tindAtom  )) deb_Atom   = .true.
         if (btest(sav,tindPolar )) deb_Polar  = .true.
      end if

      end subroutine

      module subroutine info_dyn()
      implicit none
      integer i

 13   format(A20,2I10)
 14   format(A20,F15.6)
 15   format(A20,A15)
 16   format(A20,5x,L4)
 17   format(A20,G15.4E1)

      print 13, 'natoms', n
      if (n.ne.nbloc) print 13, 'nbloc', nbloc
      print 13, 'natoms loc/nl', nloc,nlocnl
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
         print 17, 'polar solver tol', poleps
      end if
      if (use_charge) then
         print 14, 'charge cutoff', mpolecut
         print 14, 'chg taper'    , chgtaper
      end if
      print 15, 'thermostat', thermostat
      print 15, 'barostat', barostat
      end subroutine

      subroutine lookfor(val,n,array,find,ind)
      implicit none
      integer val,ind,n
      integer array(*)
      logical find
      integer i
      find=.false.
      ind=0
      do i = 1,n
         if(array(i).eq.val) then
            ind=i; find=.true.
            exit
         end if
      end do
      end subroutine

      module subroutine minmaxonei( vector,sz,name )
      implicit none
      integer sz
      integer vector(*)
      character(*),optional,intent(in)::name
      integer mi,ma
      integer(8) on,on1
      integer i
      integer val
      mi=huge(mi);ma=-mi;on=0
      do i = 1, sz
         val = (vector(i))
         mi  = min( mi,val )
         ma  = max( ma,val )
         on  = on + abs(val)
      end do
12    format(A6,2I8,1x,I12,I5)
13    format(A6,2I8,1x,I12,I5,I0)
      if (nproc.gt.1) then
         call MPI_ALLREDUCE(on,on1,1,MPI_INTEGER
     &       ,MPI_SUM,COMM_TINKER,i)
      end if
      !on1 = sqrt(on1)
      !on  = sqrt(on )
      if (rank.eq.0.and.nproc.gt.1) then
      write(*,13) name,mi,ma,on,rank,on1
      else
      write(*,12) name,mi,ma,on,rank
      end if
      end subroutine

      module subroutine minmaxonet( vector,sz,name )
      implicit none
      integer sz
      real*8 vector(*)
      character(*),optional,intent(in)::name
      real(8) mi,mi1,ma,ma1,on,on1
      integer i
      real(8) val
      mi=huge(mi);ma=tiny(ma);on=0
      do i = 1, sz
         val = (vector(i))
         mi  = min( mi,val )
         ma  = max( ma,val )
         on  = on + abs(val)
      end do
12    format(A6,2F20.8,1x,F21.5,I5)
13    format(A6,2F20.8,1x,F21.5,I5,F19.5)
      if (nproc.gt.1) then
         call MPI_ALLREDUCE(on,on1,1,MPI_REAL8
     &       ,MPI_SUM,COMM_TINKER,i)
         call MPI_ALLREDUCE(mi,mi1,1,MPI_REAL8
     &       ,MPI_MIN,COMM_TINKER,i)
         call MPI_ALLREDUCE(ma,ma1,1,MPI_REAL8
     &       ,MPI_MAX,COMM_TINKER,i)
      end if
      !on1 = sqrt(on1)
      !on  = sqrt(on )
      if (rank.eq.0.and.nproc.gt.1) then
      write(*,13) name,mi1,ma1,on,rank,on1
      else
      write(*,12) name,mi,ma,on,rank
      end if
      end subroutine

      module subroutine normt( array,n,val,p )
      implicit none
      integer n
      integer,optional::p
      real*8 array(*)
      real*8 val
      integer i,p_

      val=0
      if (present(p)) then
         do i = 1,n
            val = val + (abs(array(i)))**p
         end do
      else
         do i = 1,n
            val = val + (abs(array(i)))**2
         end do
      end if
      end subroutine

      subroutine endiam(array,n)
      implicit none
      real*8 array(*)
      integer n
      integer i
      real(8),save:: are,mi
      logical,save:: f_in=.true.

      if (f_in) then
          f_in=.false.
      end if

      are = 0
      mi  = 0
      do i = 1,n
         are = are +     array(i)
         mi  =  mi + abs(array(i))
      end do
      print*, are,mi
      end subroutine

      subroutine endiam1(array,n)
      implicit none
      real*8 array(*)
      integer n
      integer i
      real(8),save:: are,mi
      logical,save:: f_in=.true.

      if (f_in) then
          f_in=.false.
      end if

      are = 0; mi = 0
      do i = 1,n
         are = max(are , abs(array(i)))
         mi  = mi + abs(array(i))
      end do
      print*, are,mi
      end subroutine

c
c     Print information on position, velocities and aceleration
c
      module subroutine info_minmax_pva(opt)
      implicit none
      integer,optional:: opt
      integer i,iglob,opt_
      integer,parameter::nel=5
      real*8 minmax(3*nel)
      real*8 tmp
      logical tinker_isnan

      opt_ = 0
      if (present(opt)) opt_=opt
      minmax=0
 20   format (80('_'))

      do i = 1, nloc
         iglob      = glob(i)
         if (muse(iglob)) then
         minmax(01) = min(minmax(01),x(iglob))
         minmax(02) = min(minmax(02),y(iglob))
         minmax(03) = min(minmax(03),z(iglob))
         tmp = min(min(v(1,iglob),v(2,iglob)),v(3,iglob))
         minmax(04) = min(minmax(04),tmp)
         if (opt_.eq.0) then
            tmp = min(min(a(1,iglob),a(2,iglob)),a(3,iglob))
         else if(opt_.eq.1) then
            tmp = min(min(aalt(1,iglob),aalt(2,iglob)),aalt(3,iglob))
         end if
         minmax(05) = min(minmax(05),tmp)
         minmax(06) = max(minmax(06),x(iglob))
         minmax(07) = max(minmax(07),y(iglob))
         minmax(08) = max(minmax(08),z(iglob))

         tmp = max(max(v(1,iglob),v(2,iglob)),v(3,iglob))
         minmax(09) = max(minmax(09),v(1,iglob))
         if (opt_.eq.0) then
            tmp = max(max(a(1,iglob),a(2,iglob)),a(3,iglob))
         else if(opt_.eq.1) then
            tmp = max(max(aalt(1,iglob),aalt(2,iglob)),aalt(3,iglob))
         end if
         minmax(10) = max(minmax(10),tmp)

         tmp = v(1,iglob)+v(2,iglob)+v(3,iglob)
         minmax(11) = minmax(11)+tmp
         if (opt_.eq.0) then
            tmp = a(1,iglob)+a(2,iglob)+a(3,iglob)
         else if(opt_.eq.1) then
            tmp = aalt(1,iglob)+aalt(2,iglob)+aalt(3,iglob)
         end if
         minmax(12) = minmax(12)+tmp
         end if
      end do

      if (rank.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,minmax,nel,MPI_REAL8,
     &                   MPI_MIN,0,COMM_TINKER,i)
         call MPI_REDUCE(MPI_IN_PLACE,minmax(6),nel,MPI_REAL8,
     &                   MPI_MAX,0,COMM_TINKER,i)
         call MPI_REDUCE(MPI_IN_PLACE,minmax(11),2,MPI_REAL8,
     &                   MPI_SUM,0,COMM_TINKER,i)
      else
         call MPI_REDUCE(minmax,minmax,nel,MPI_REAL8,
     &                   MPI_MIN,0,COMM_TINKER,i)
         call MPI_REDUCE(minmax(6),minmax(6),nel,MPI_REAL8,
     &                   MPI_MAX,0,COMM_TINKER,i)
         call MPI_REDUCE(minmax(11),minmax(11),2,MPI_REAL8,
     &                   MPI_SUM,0,COMM_TINKER,i)
      end if

      if (rank.eq.0) then
30    format(A10,2F18.8)
32    format(A10,3F18.8)
         print 20
         print 30,"min max_x ",minmax(01),minmax(06)
         print 30,"min max_y ",minmax(02),minmax(07)
         print 30,"min max_z ",minmax(03),minmax(08)
         print 32,"min max_v ",minmax(04),minmax(09),minmax(11)
         if (opt_.eq.0) then
         print 32,"min max_a ",minmax(05),minmax(10),minmax(12)
         elseif (opt_.eq.1) then
         print 32,"min max_a1",minmax(05),minmax(10),minmax(12)
         elseif (opt_.eq.2) then
         print 32,"min max_a2",minmax(05),minmax(10),minmax(12)
         end if
      end if
      end subroutine
      end submodule
