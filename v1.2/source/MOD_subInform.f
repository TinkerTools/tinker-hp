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
      use freeze,only: use_rattle
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
      end if
      end subroutine

      module subroutine info_dyn()
      use iounit
      implicit none
      integer i

 13   format(A20,2I10)
 14   format(A20,F15.6)
 15   format(A20,A15)
 16   format(A20,5x,L4)
 17   format(A20,G15.4E1)

      write(iout,13) 'natoms', n
      if (n.ne.nbloc) write(iout,13), 'nbloc', nbloc
      write(iout,13), 'natoms loc/nl', nloc,nlocnl
      write(iout,13), 'nlupdate'     , ineigup
      write(iout,14), 'list buffer'  , lbuffer
      if (use_vdw) then
         write(iout,14), 'vdw cutoff'   , vdwcut
         write(iout,14), 'vdw short cutoff'   , vdwshortcut
         write(iout,14), 'vdw taper'    , vdwtaper
         write(iout,14), 'shortheal'    , shortheal
      end if
      if (use_mpole.or.use_polar) then
         write(iout,14), 'mpole cutoff' , mpolecut
         write(iout,14), 'mpole short cutoff' , mpoleshortcut
         write(iout,14), 'ewaldcut'     , ewaldcut
         write(iout,17), 'polar solver tol', poleps
      end if
      if (use_charge) then
         write(iout,14), 'charge cutoff', mpolecut
         write(iout,14), 'chg taper'    , chgtaper
      end if
      write(iout,15), 'thermostat', thermostat
      write(iout,15), 'barostat', barostat
      write(iout,16), 'rattle', use_rattle
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

      ! Debug routine on Forces
      subroutine info_forces(rule)
      use atoms
#ifdef COLVARS
      use colvars
#endif
      use deriv
      use domdec
      use inform
      use iounit
      use mpi
      use potent
      implicit none
      integer,intent(in) :: rule
      integer,parameter:: nf=28  !number of forces
      integer i,j,sze,ids
      real(8) mmx(3*nf),maxi
      logical tinker_isnan_m, save_arc,save_dcdio,focus_nbond,abortf
      real(8) dt
      integer,save:: doin = 1

      enum,bind(C)
      enumerator ::commBonded
      enumerator ::commShortNonBonded
      enumerator ::commNonBonded
      end enum

      mmx    = 0
      sze    = 3*nloc

      if (deb_Path) write(iout,*) 'info_forces', rule

      if (btest(rule,idBond)) then

         if(use_bond)   call comm_forces_dd(deb  ,cBond)
         if(use_angle)  call comm_forces_dd(dea  ,cBond)
         if(use_strbnd) call comm_forces_dd(deba ,cBond)
         if(use_urey)   call comm_forces_dd(deub ,cBond)
         if(use_angtor) call comm_forces_dd(deat ,cBond)
         if(use_improp) call comm_forces_dd(deid ,cBond)
         if(use_imptor) call comm_forces_dd(deit ,cBond)
         if(use_tors)   call comm_forces_dd(det  ,cBond)
         if(use_pitors) call comm_forces_dd(dept ,cBond)
         if(use_tortor) call comm_forces_dd(dett ,cBond)
         if(use_opbend) call comm_forces_dd(deopb,cBond)
         if(use_strtor) call comm_forces_dd(debt ,cBond)
         if(use_opdist) call comm_forces_dd(deopd,cBond)
         if(use_angang) call comm_forces_dd(deaa ,cBond)
         if(use_geom)   call comm_forces_dd(deg  ,cBond)
         if(use_extra)  call comm_forces_dd(dex  ,cBond)
      end if


      if (btest(rule,idSNBond)) then

      if (use_vdw) then
         call comm_forces_dd(dev ,cSNbond)
         call minmaxone1(mmx(16),mmx(nf+16),mmx(2*nf+16),dev
     &                  ,sze,'dev');
      end if
      if (use_charge) then
         call comm_forces_dd(dec,cSNbond)
         call minmaxone1(mmx(17),mmx(nf+17),mmx(2*nf+17),dec
     &                  ,sze,'dec');
      end if
      if (use_mpole) then
         call comm_forces_dd(dem,cSNbond)
         call minmaxone1(mmx(19),mmx(nf+19),mmx(2*nf+19),dem
     &                  ,sze,'dem');
      end if
      if (use_polar) then
         call comm_forces_dd(dep,cSNbond)
         call minmaxone1(mmx(21),mmx(nf+21),mmx(2*nf+21),dep
     &                  ,sze,'dep');
      end if

      end if

      if (btest(rule,idNBond)) then

      if (use_vdw) then
         call comm_forces_dd(dev,cNBond)
         call minmaxone1(mmx(16),mmx(nf+16),mmx(2*nf+16),dev
     &                  ,sze,'dev')
      end if

      if(use_charge) then
         call comm_forces_dd(dec,cNBond)
         call minmaxone1(mmx(17),mmx(nf+17),mmx(2*nf+17),dec
     &                  ,sze,'dec');
         call commforcesrec1(3)
         call minmaxone1(mmx(18),mmx(nf+18),mmx(2*nf+18),dec
     &                  ,sze,'decsum');
      end if

      if(use_mpole) then
         call comm_forces_dd(dem,cNBond)
         call minmaxone1(mmx(19),mmx(nf+19),mmx(2*nf+19),dem
     &                  ,sze,'dem');
         call commforcesrec1(1)
         call minmaxone1(mmx(20),mmx(nf+20),mmx(2*nf+20),dem
     &                  ,sze,'demsum');
      end if

      if(use_polar) then
         call comm_forces_dd(dep,cNBond)
         call minmaxone1(mmx(21),mmx(nf+21),mmx(2*nf+21),dep
     &                  ,sze,'dep');
         call commforcesrec1(2)
         call minmaxone1(mmx(22),mmx(nf+22),mmx(2*nf+22),dep
     &                  ,sze,'depsum');
      end if

      if (use_chgtrn) then
         call comm_forces_dd(dect,cNBond)
         call minmaxone1(mmx(28),mmx(nf+28),mmx(2*nf+28),dect
     &                  ,sze,'dect');
      end if

      end if

      if (btest(rule,idBond)) then

      if(use_bond)
     &call minmaxone1(mmx(01),mmx(nf+01),mmx(2*nf+01),deb
     &     ,sze,'deb');
      if(use_angle)
     &call minmaxone1(mmx(02),mmx(nf+02),mmx(2*nf+02),dea
     &     ,sze,'dea');
      if(use_strbnd)
     &call minmaxone1(mmx(03),mmx(nf+03),mmx(2*nf+03),deba
     &     ,sze,'deba');
      if(use_urey)
     &call minmaxone1(mmx(04),mmx(nf+04),mmx(2*nf+04),deub
     &     ,sze,'deub');
      if(use_angang)
     &call minmaxone1(mmx(05),mmx(nf+05),mmx(2*nf+05),deaa
     &     ,sze,'deaa');
      if(use_improp)
     &call minmaxone1(mmx(06),mmx(nf+06),mmx(2*nf+06),deid
     &     ,sze,'deid');
      if(use_imptor)
     &call minmaxone1(mmx(07),mmx(nf+07),mmx(2*nf+07),deit
     &     ,sze,'deit');
      if(use_tors)
     &call minmaxone1(mmx(08),mmx(nf+08),mmx(2*nf+08),det
     &     ,sze,'det');
      if(use_pitors)
     &call minmaxone1(mmx(09),mmx(nf+09),mmx(2*nf+09),dept
     &     ,sze,'dept');
      if(use_strtor)
     &call minmaxone1(mmx(10),mmx(nf+10),mmx(2*nf+10),debt
     &     ,sze,'debt');
      if(use_tortor)
     &call minmaxone1(mmx(11),mmx(nf+11),mmx(2*nf+11),dett
     &     ,sze,'dett');
      if(use_angtor)
     &call minmaxone1(mmx(25),mmx(nf+25),mmx(2*nf+25),deat
     &     ,sze,'deat');
      if(use_opbend)
     &call minmaxone1(mmx(12),mmx(nf+12),mmx(2*nf+12),deopb
     &     ,sze,'deopb');
      if(use_opdist)
     &call minmaxone1(mmx(13),mmx(nf+13),mmx(2*nf+13),deopd
     &     ,sze,'deopd');
      if(use_geom)
     &call minmaxone1(mmx(14),mmx(nf+14),mmx(2*nf+14),deg
     &     ,sze);
      if(use_extra)
     &call minmaxone1(mmx(15),mmx(nf+15),mmx(2*nf+15),dex
     &     ,sze,'dex');

      end if


c      if (btest(rule,4)) then
c
cc      if (abort) then
cc         dint1 = minval(inte); dint2=maxval(inte)
cc         call searchpair(nshortvlst,shortvlst,maxvlst
cc     &                  ,dint1,dint2)
cc      end if
c
c      else if (btest(rule,5)) then
c
c      end if

#ifdef COLVARS
      if (use_colvars.and.ncvatoms.gt.0) then
      call minmaxone1(mmx(26),mmx(nf+26),mmx(2*nf+26),decv_tot
     &              ,3*ncvatoms,'decolv');
      call minmaxone1(mmx(27),mmx(nf+27),mmx(2*nf+27),decv
     &              ,3*ncvatoms,'decolv');
      end if
#endif

      if (nproc.gt.1) then
      if (rank.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,mmx,nf,MPI_REAL8,
     &                   MPI_MIN,0,COMM_TINKER,i)
         call MPI_REDUCE(MPI_IN_PLACE,mmx(nf+1),nf,MPI_REAL8,
     &                   MPI_MAX,0,COMM_TINKER,i)
         call MPI_REDUCE(MPI_IN_PLACE,mmx(2*nf+1),nf,MPI_REAL8,
     &                   MPI_SUM,0,COMM_TINKER,i)
      else
         call MPI_REDUCE(mmx,mmx,nf,MPI_REAL8,
     &                   MPI_MIN,0,COMM_TINKER,i)
         call MPI_REDUCE(mmx(nf+1),mmx(nf+1),nf,MPI_REAL8,
     &                   MPI_MAX,0,COMM_TINKER,i)
         call MPI_REDUCE(mmx(2*nf+1),mmx(2*nf+1),nf,MPI_REAL8,
     &                   MPI_SUM,0,COMM_TINKER,i)
      end if
      end if

      maxi=0
c      ids = merge(2,1,fdebs_l)
      do i = 1,26; maxi = max(mmx(nf+i),maxi); end do
c      if (focus_nbond) maxi = maxval(mmx(nf+16:nf+21))

c      if ((maxi.gt.f_ulim).and.f_ulim.gt.0)
c     &   then
c         n_adjust = n_adjust + 1
c 11      format('info_forces: above frc uplimit',F8.2,' cnt',I10)
c         if(rank.eq.0) write(*,11) f_ulim,n_adjust
c
c         if (maxi.gt.f_ulim+7) then
c         read(arg(3),*) dt
c         dt = 1d-3*dt
c         new_restart = .true.
c         f_mdsave    = .true.
c         save_dcdio  = dcdio
c         save_arc    = archive
c         dcdio       = .false.
c         archive     = .false.
c         n_fwriten   = n_fwriten + 1
c         call mdsave(step_c,dt,epot)
c         new_restart = .false.
c         f_mdsave    = .false.
c         dcdio       = save_dcdio
c         archive     = save_arc
c         end if
c      end if

      abortf = merge(.true.,.false.,maxi.gt.10000)
      if (abortf) then
 24   format(' info_forces ! Abnormal forces detected ',/)
         write(0,24)
         abort=.true.
      end if

      if (rank.eq.0) then
 30      format(a10,3F20.8)
 40      format(80('='))
         print 40
c         if (fdebs_l) then
c         if(mmx(2*nf+01)/=0.0) print 30,'de_tot =>',
c     &      mmx(01),mmx(01+nf),mmx(01+nf*2)
c         else
         if(mmx(2*nf+01)/=0.0) print 30,'deb    =>',
     &      mmx(01),mmx(01+nf),mmx(01+nf*2)
         if(mmx(2*nf+02)/=0.0) print 30,'dea    =>',
     &      mmx(02),mmx(02+nf),mmx(02+nf*2)
         if(mmx(2*nf+03)/=0.0) print 30,'deba   =>',
     &      mmx(03),mmx(03+nf),mmx(03+nf*2)
         if(mmx(2*nf+04)/=0.0) print 30,'deub   =>',
     &      mmx(04),mmx(04+nf),mmx(04+nf*2)
         if(mmx(2*nf+05)/=0.0) print 30,'deaa   =>',
     &      mmx(05),mmx(05+nf),mmx(05+nf*2)
         if(mmx(2*nf+06)/=0.0) print 30,'deid   =>',
     &      mmx(06),mmx(06+nf),mmx(06+nf*2)
         if(mmx(2*nf+07)/=0.0) print 30,'deit   =>',
     &      mmx(07),mmx(07+nf),mmx(07+nf*2)
         if(mmx(2*nf+08)/=0.0) print 30,'det    =>',
     &      mmx(08),mmx(08+nf),mmx(08+nf*2)
         if(mmx(2*nf+09)/=0.0) print 30,'dept   =>',
     &      mmx(09),mmx(09+nf),mmx(09+nf*2)
         if(mmx(2*nf+10)/=0.0) print 30,'debt   =>',
     &      mmx(10),mmx(10+nf),mmx(10+nf*2)
         if(mmx(2*nf+11)/=0.0) print 30,'dett   =>',
     &      mmx(11),mmx(11+nf),mmx(11+nf*2)
         if(mmx(2*nf+25)/=0.0) print 30,'deat   =>',
     &      mmx(25),mmx(25+nf),mmx(25+nf*2)
         if(mmx(2*nf+12)/=0.0) print 30,'deopb  =>',
     &      mmx(12),mmx(12+nf),mmx(12+nf*2)
         if(mmx(2*nf+13)/=0.0) print 30,'deopd  =>',
     &      mmx(13),mmx(13+nf),mmx(13+nf*2)
         if(mmx(2*nf+14)/=0.0) print 30,'deg    =>',
     &      mmx(14),mmx(14+nf),mmx(14+nf*2)
         if(mmx(2*nf+15)/=0.0) print 30,'dex    =>',
     &      mmx(15),mmx(15+nf),mmx(15+nf*2)
         if(mmx(2*nf+16)/=0.0) print 30,'dev    =>',
     &      mmx(16),mmx(16+nf),mmx(16+nf*2)
         if(mmx(2*nf+17)/=0.0) print 30,'dec    =>',
     &      mmx(17),mmx(17+nf),mmx(17+nf*2)
         if(mmx(2*nf+18)/=0.0) print 30,'decsum =>',
     &      mmx(18),mmx(18+nf),mmx(18+nf*2)
         if(mmx(2*nf+19)/=0.0) print 30,'dem    =>',
     &      mmx(19),mmx(19+nf),mmx(19+nf*2)
         if(mmx(2*nf+20)/=0.0) print 30,'demsum =>',
     &      mmx(20),mmx(20+nf),mmx(20+nf*2)
         if(mmx(2*nf+21)/=0.0) print 30,'dep    =>',
     &      mmx(21),mmx(21+nf),mmx(21+nf*2)
         if(mmx(2*nf+22)/=0.0) print 30,'depsum =>',
     &      mmx(22),mmx(22+nf),mmx(22+nf*2)
         if(mmx(2*nf+28)/=0.0) print 30,'dect   =>',
     &      mmx(28),mmx(28+nf),mmx(28+nf*2)
         if(mmx(2*nf+23)/=0.0) print 30,'deamdD =>',
     &      mmx(23),mmx(23+nf),mmx(23+nf*2)
         if(mmx(2*nf+24)/=0.0) print 30,'deW1aMD=>',
     &      mmx(24),mmx(24+nf),mmx(24+nf*2)
         if(mmx(2*nf+26)/=0.0) print 30,'declvi =>',
     &      mmx(26),mmx(26+nf),mmx(26+nf*2)
         if(mmx(2*nf+27)/=0.0) print 30,'declvo =>',
     &      mmx(27),mmx(27+nf),mmx(27+nf*2)
         if(mmx(2*nf+28)/=0.0) print 30,'dmlpot =>',
     &      mmx(28),mmx(28+nf),mmx(28+nf*2)
      end if

c      call reset_forces_d
c      if (btest(rule,idNBond)) then
c         call zero_forces_rec
c      end if
c
c      doin = doin +1
      end subroutine
