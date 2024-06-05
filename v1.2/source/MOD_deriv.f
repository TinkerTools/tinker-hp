c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module deriv  --  Cartesian coordinate derivative components  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     desum   total energy Cartesian coordinate derivatives
c     deb     bond stretch Cartesian coordinate derivatives
c     dea     angle bend Cartesian coordinate derivatives
c     deba    stretch-bend Cartesian coordinate derivatives
c     deub    Urey-Bradley Cartesian coordinate derivatives
c     deaa    angle-angle Cartesian coordinate derivatives
c     deopb   out-of-plane bend Cartesian coordinate derivatives
c     deopd   out-of-plane distance Cartesian coordinate derivatives
c     deid    improper dihedral Cartesian coordinate derivatives
c     deit    improper torsion Cartesian coordinate derivatives
c     det     torsional Cartesian coordinate derivatives
c     dept    pi-orbital torsion Cartesian coordinate derivatives
c     debt    stretch-torsion Cartesian coordinate derivatives
c     deat    angle-torsion Cartesian coordinate derivatives
c     dett    torsion-torsion Cartesian coordinate derivatives
c     dev     van der Waals Cartesian coordinate derivatives
c     der     repulsion Cartesian coordinate derivatives
c     dedsp   dispersion Cartesian coordinate derivatives
c     dedsprec   reciprocal dispersion Cartesian coordinate derivatives
c     dect    charge transfer Cartesian coordinate derivatives
c     dec     charge-charge Cartesian coordinate derivatives
c     decrec  reciprocal charge-charge Cartesian coordinate derivatives
c     dem     multipole Cartesian coordinate derivatives
c     demrec  reciprocal multipole Cartesian coordinate derivatives
c     dep     polarization Cartesian coordinate derivatives
c     dep     reciprocal polarization Cartesian coordinate derivatives
c     deg     geometric restraint Cartesian coordinate derivatives
c     dex     extra energy term Cartesian coordinate derivatives
c     desave  stored Cartesian coordinate derivatives
c     desmd   extra smd energy term Cartesian coordinate derivatives
c
c     Lambda-dynamics derivatives
c
c     delambda           hamiltonian derivative with respect to lambda (to be sent to colvar)
c     delambdae          hamiltonian derivative with respect to elambda
c     delambdav          hamiltonian derivative with respect to vlambda
c     delambdaesave      stored hamiltonian derivative with respect to elambda
c     delambdavsave      stored hamiltonian derivative with respect to vlambda
c     dlambdaelambda     derivative of elambda with respect to lambda
c     dlambdavlambda     derivative of vlambda with respect to lambda     
c
c     Orthogonal Space Random Walk - note x stands for Cartesian coordinates
c     dxdelambda         hamiltonian double derivative with respect to x and lambda (to be sent to colvar)
c     dxdelambdae        hamiltonian double derivative with respect to x and elambda (electrostatic interactions)
c     dxdelambdav        hamiltonian double derivative with respect to x and vlambda (vdw interactions)
c     d2edlambda2         hamiltonian double derivative with respect to lambda (to be sent to colvar)
c     d2edlambdae2        hamiltonian double derivative with respect to elambda (electrostatic interactions)
c     d2edlambdav2        hamiltonian double derivative with respect to vlambda (vdw interactions)
c
c     dotstgrad : flag when the main program is testgrad (communication
c      of the forces one by one)
c
c
      module deriv
      implicit none
      real*8, allocatable :: desum(:,:),deb(:,:),dea(:,:),deba(:,:)
      real*8, allocatable :: deub(:,:),deaa(:,:),deopb(:,:),deopd(:,:)
      real*8, allocatable :: deid(:,:),det(:,:),dept(:,:),deit(:,:)
      real*8, allocatable :: deat(:,:)
      real*8, allocatable :: debt(:,:),dett(:,:),dev(:,:),dec(:,:)
      real*8, allocatable :: der(:,:),dedsp(:,:),dect(:,:)
      real*8, allocatable :: dedsprec(:,:)
      real*8, allocatable :: dem(:,:),dep(:,:)
      real*8, allocatable :: deg(:,:),dex(:,:)
      real*8, allocatable :: decrec(:,:),demrec(:,:),deprec(:,:)
      real*8, allocatable :: debond(:,:),desave(:,:)
      real*8, allocatable :: desmd(:,:)
      real*8 :: delambda,delambdae,delambdav
      real*8 :: delambdaesave,delambdavsave
      real*8 :: d2edlambda2,d2edlambdae2,d2edlambdav2
      real*8 :: dlambdaelambda, dlambdavlambda
      real*8, allocatable  :: dxdelambda(:,:)
      real*8, allocatable :: dxdelambdae(:,:), dxdelambdav(:,:)
      logical dotstgrad
      logical abortall
      integer inte(2)
      integer cBond,cNBond,cSNBond,cDef
      enum,bind(C)
      enumerator idBond,idSNBond,idNBond
      end enum
      parameter( cBond=1,cSNBond=2,cNBond=4,cDef=5 )
      save

      contains

      subroutine resetForcesRec
      implicit none
        if(allocated(demrec)) demrec = 0.0d0
        if(allocated(decrec)) decrec = 0.0d0
        if(allocated(deprec)) deprec = 0.0d0
        if(allocated(dedsprec)) dedsprec = 0.0d0
      end subroutine resetForcesRec

      ! Debug routine on Forces
      subroutine info_forces(rule)
      use atoms
#ifdef COLVARS
      use colvars
#endif
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

      subroutine minmaxone1( mi,ma,on,vector,sz,name )
      use atoms
      use domdec
      use inform
      use mpi
      implicit none
      integer sz
      integer,parameter::lgli=10
      real(8) mi,ma,on
      real*8 vector(*)
      character(*),optional,intent(in)::name
      integer i,j,i1,i2,iglob,cap,cap1,cap2
      integer gli(lgli,nproc)
      real(8) val

      abortall = .false.
      if (present(name)) then
         cap = 1
         gli = 0
         if (name.eq.'devi') then
            do i = 1, sz/3
               cap2  = 0
               iglob = glob(i)
               do j  = 1,3
               val   = vector((i-1)*3+j)
               if (abs(val).gt.90.0) then
                  print*,j,iglob,rank,x(iglob),y(iglob)
     &                  ,z(iglob),val
                  abort=.true.
                  cap2 = cap2 + 1
               end if
               end do
               if (cap2.gt.0) then
               cap1 = cap
               cap  = cap + 1
               if (cap1.le.lgli) gli(cap1,rank+1) = iglob
               end if
            end do
            do i = 1, nproc
               if (rank.eq.i-1) abortall=abort
               call MPI_BCAST(abortall,1,MPI_LOGICAL,i-1,COMM_TINKER,i1)
               if (abortall) then
                  abort = .true.
                  call MPI_AllGather(MPI_IN_PLACE,lgli,MPI_DATATYPE_NULL
     $                              ,gli,lgli,MPI_INT,COMM_TINKER,i1)
                  exit
               end if
            end do

            if (abort) then
            
            cap  = 0
            cap1 = 0
            do i1 = 1, nproc; do i2 = 1,lgli
               if ( gli(i2,i1).ne.0 ) then
                  if (cap.eq.0) then
                     cap = 1
                     inte(cap) = gli(i2,i1)
                  else if (cap.gt.0.and.abs(gli(i2,i1)-inte(1)).lt.5)
     &                 then
                     cap = cap + 1
                     inte(2) = gli(i2,i1)
                  else
                     cap = cap + 1
                  end if
               end if
            end do; end do
            
            if (cap.ne.2.and.rank.eq.0) then
               print*,' more than one interactions found '
     &               ,rank,cap
            do i = 1,nproc
               do j = 1,lgli
                  if (gli(j,i).ne.0) write(*,'(I10,$)') gli(j,i)
               end do
               print*
            end do
            end if
            
            end if
         end if

      end if

      do i = 1, sz
         val = vector(i)
         mi = min( mi,val )
         ma = max( ma,val )
         on = on + abs(val)
      end do
      end subroutine

      end
