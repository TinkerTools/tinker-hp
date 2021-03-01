c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module subderiv                                               ##
c     ##  --  Routines associated to derivatives component managment    ##
c     ##                                                                ##
c     ####################################################################
c
#include "tinker_precision.h"
#include "tinker_types.h"

      submodule(deriv) subderiv
      use atoms    ,only: x,y,z,n
      use domdec
      use inform   ,only:deb_Path,abort,dint1,dint2
      use mpi
      use neigh
      use potent
      use sizes
      use tinheader,only:ti_p,re_p
      implicit none
      logical abortall
      integer inte(2)

      contains

#include "convert.f.inc"

      module subroutine ConfigPotential
      implicit none
      integer watchPot
      watchPot = 0

      PotentialAmoeba = use_bond.and.use_angle.and.use_strbnd
     &      .and.use_urey.and..not.use_angang.and.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and.use_tors.and.use_pitors.and..not.use_strtor
     &      .and.use_tortor.and.use_vdw.and..not.use_charge
     &      .and.use_mpole.and.use_polar.and..not.use_extra
      PotentialAmoeba18 = use_bond.and.use_angle.and.use_strbnd
     &      .and.use_urey.and..not.use_angang.and.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and.use_tors.and.use_pitors.and.use_strtor
     &      .and.use_tortor.and.use_vdw.and..not.use_charge
     &      .and.use_mpole.and.use_polar.and..not.use_extra
      PotentialWaterAmoeba = use_bond.and.use_angle.and..not.use_strbnd
     &      .and.use_urey.and..not.use_angang.and..not.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and..not.use_tors.and..not.use_pitors.and..not.use_strtor
     &      .and..not.use_tortor.and.use_vdw.and..not.use_charge
     &      .and.use_mpole.and.use_polar.and..not.use_extra
      PotentialWaterCharmm = use_bond.and.use_angle.and..not.use_strbnd
     &      .and.use_urey.and..not.use_angang.and..not.use_opbend
     &      .and..not.use_opdist.and..not.use_improp.and..not.use_imptor
     &      .and..not.use_tors.and..not.use_pitors.and..not.use_strtor
     &      .and..not.use_tortor.and.use_vdw.and.use_charge
     &      .and..not.use_mpole.and..not.use_polar.and..not.use_extra
      PotentialCharmm = use_bond.and.use_angle.and..not.use_strbnd
     &      .and.use_urey.and..not.use_angang.and..not.use_opbend
     &      .and..not.use_opdist.and.use_improp.and..not.use_imptor
     &      .and.use_tors.and..not.use_pitors.and..not.use_strtor
     &      .and..not.use_tortor.and.use_vdw.and.use_charge
     &      .and..not.use_mpole.and..not.use_polar.and..not.use_extra

      if (PotentialAmoeba) then
         if (deb_Path) print*,'Using Amoeba Potential'
         PotentialAll=.false.
         watchPot = watchPot + 1
         resetForces_p => resetForcesAmoeba
         addForces_p   =>   addForcesAmoeba
      end if
      if (PotentialAmoeba18) then
         if (deb_Path) print*,'Using Amoeba18 Potential'
         PotentialAll=.false.
         watchPot = watchPot + 1
         resetForces_p => resetForcesAmoeba18
         addForces_p   =>   addForcesAmoeba18
      end if
      if (PotentialWaterAmoeba) then
         if (deb_Path) print*,'Using Water Amoeba Potential'
         PotentialAll=.false.
         watchPot = watchPot + 1
         resetForces_p => resetForcesWaterAm
         addForces_p   =>   addForcesWaterAm
      end if
      if (PotentialWaterCharmm) then
         PotentialAll=.true.
         watchPot = watchPot + 1
         resetForces_p => resetForcesDefault
         addForces_p   =>   addForcesDefault
      end if
      if (PotentialCharmm) then
         PotentialAll=.false.
         if (deb_Path) print*,'Using Charmm Potential'
         watchPot = watchPot + 1
         resetForces_p => resetForcesCharmm
         addForces_p   =>   addForcesCharmm
      end if
      if (PotentialAll) then
         if (deb_Path) print*,'Using Default Potential'
         watchPot = watchPot + 1
         resetForces_p => resetForcesDefault
         addForces_p   =>   addForcesDefault
      end if

      if (watchPot.gt.1) then
12       format ( " WARNING !!! ConfigPotential Routine ",/,
     &   " Found two of more differents configurations on running",
     &   " potentials ",/,
     &   " Switched to default potential by reserving space for all",
     &   " forces ")
         print 12

         PotentialAll         = .true.
         PotentialAmoeba      = .false.
         PotentialAmoeba18    = .false.
         PotentialCharmm      = .false.
         PotentialWaterAmoeba = .false.
         PotentialWaterCharmm = .false.
         resetForces_p => resetForcesDefault
         addForces_p   =>   addForcesDefault
      end if

      end subroutine

      subroutine updateForcesHost
      implicit none

      if (PotentialAll) then
!$acc update host(dea,deb,deba,deub,deaa,deopb,deopd,
!$acc&       deg,deid,deit,det,dept,debt,dett,dex,
!$acc&       deamdD,dev,
!$acc&       dem,demrec,dep,deprec,dec,decrec,desum)
      else if (PotentialAmoeba) then
!$acc update host(dea,deb,deba,deub,deopb,
!$acc&       det,dept,dett,
!$acc&       dev,dem,demrec,dep,deprec,desum)
         if (allocated(deg)) then
!$acc update host(deg)
         end if
      else if (PotentialAmoeba18) then
!$acc update host(dea,deb,deba,deub,deopb,
!$acc&       det,dept,debt,dett,
!$acc&       dev,dem,demrec,dep,deprec,desum)
      else if (PotentialCharmm) then
!$acc update host(dea,deb,deub,
!$acc&       deid,
!$acc&       dev,dec,decrec,desum)
      else if (PotentialWaterAmoeba) then
!$acc update host(dea,deb,deub,
!$acc&       dev,dem,demrec,dep,deprec,desum)
      else
         print*,'Unkown potential for force Updates'
      end if

      end subroutine

      module subroutine resetForcesAmoeba
      implicit none
      integer i,j
      if (deb_Path) print*, 'resetForcesAmoeba'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            deb  (j,i) = 0.0_re_p ! ebond
            dea  (j,i) = 0.0_re_p ! eangle
            deba (j,i) = 0.0_re_p ! estrbnd
            deub (j,i) = 0.0_re_p ! eurey
            deopb(j,i) = 0.0_re_p ! eopbend
            det  (j,i) = 0.0_re_p ! etors
            dept (j,i) = 0.0_re_p ! epitors
            dett (j,i) = 0.0_re_p ! etortor
            dev  (j,i) = 0.0_re_p ! ehal1
            dem  (j,i) = 0        ! empole
            dep  (j,i) = 0        ! epolar
         end do
      end do
      if (use_geom) call resetForcesRestrain
      if (use_smd_velconst .or. use_smd_forconst) call resetForcesSMD
      end subroutine
      module subroutine addForcesAmoeba
      implicit none
      integer i,j
      if (deb_Path) print*, 'addForcesAmoeba'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            desum(j,i) = deb  (j,i) + dea(j,i) + deba(j,i) + deub(j,i)
     &                 + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                 + dev  (j,i) + mdr2md(dem(j,i) + dep (j,i))

c          debond(j,i) = deb  (j,i) + dea(j,i) + deba(j,i) + deub(j,i)
c    &                 + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
         end do
      end do
      if (use_geom) call addForcesRestrain
      if (use_smd_velconst .or. use_smd_forconst) call addForcesSMD
      end subroutine

      module subroutine resetForcesAmoeba18
      implicit none
      integer i,j
      if (deb_Path) print*, 'resetForcesAmoeba18'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            deb  (j,i) = 0.0_re_p ! ebond
            dea  (j,i) = 0.0_re_p ! eangle
            deba (j,i) = 0.0_re_p ! estrbnd
            deub (j,i) = 0.0_re_p ! eurey
            deopb(j,i) = 0.0_re_p ! eopbend
            det  (j,i) = 0.0_re_p ! etors
            dept (j,i) = 0.0_re_p ! epitors
            debt (j,i) = 0.0_re_p ! estrtor
            dett (j,i) = 0.0_re_p ! etortor
            dev  (j,i) = 0.0_re_p ! ehal1
            dem  (j,i) = 0        ! empole
            dep  (j,i) = 0        ! epolar
         end do
      end do
      if (use_geom) call resetForcesRestrain
      if (use_smd_velconst .or. use_smd_forconst) call resetForcesSMD
      end subroutine
      module subroutine addForcesAmoeba18
      implicit none
      integer i,j

      if (deb_Path) print*, 'addForcesAmoeba18'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            desum(j,i) = deb  (j,i) + dea(j,i) + deba(j,i) + deub(j,i)
     &                 + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
     &                 + debt (j,i)
     &                 + dev  (j,i) + mdr2md(dem(j,i) + dep (j,i))

c          debond(j,i) = deb  (j,i) + dea(j,i) + deba(j,i) + deub(j,i)
c    &                 + deopb(j,i) + det(j,i) + dept(j,i) + dett(j,i)
c    &                 + debt (j,i)
         end do
      end do
      if (use_geom) call addForcesRestrain
      if (use_smd_velconst .or. use_smd_forconst) call addForcesSMD
      end subroutine

      module subroutine resetForcesWaterAm
      implicit none
      integer i,j
      if (deb_Path) print*, 'resetForcesWaterAm'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            deb  (j,i) = 0.0_re_p ! ebond
            dea  (j,i) = 0.0_re_p ! eangle
            deub (j,i) = 0.0_re_p ! eurey
            dev  (j,i) = 0.0_re_p ! ehal1
            dem  (j,i) = 0        ! empole
            dep  (j,i) = 0        ! epolar
         end do
      end do
      if (use_geom) call resetForcesRestrain
      if (use_smd_velconst .or. use_smd_forconst) call resetForcesSMD
      end subroutine
      module subroutine addForcesWaterAm
      implicit none
      integer i,j

      if (deb_Path) print*, 'addForcesWaterAm'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deub(j,i)
     &                 + dev(j,i) + mdr2md(dem(j,i) + dep (j,i))
c          debond(j,i) = deb(j,i) + dea(j,i) + deub(j,i)
         end do
      end do
      if (use_geom) call addForcesRestrain
      if (use_smd_velconst .or. use_smd_forconst) call addForcesSMD
      end subroutine

      module subroutine resetForcesCharmm
      implicit none
      integer i,j
      if (deb_Path) print*,'resetForcesCharmm'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            deb  (j,i) = 0.0_re_p ! ebond
            dea  (j,i) = 0.0_re_p ! eangle
            deub (j,i) = 0.0_re_p ! eurey
            deid (j,i) = 0.0_re_p ! eimprop
            det  (j,i) = 0.0_re_p ! etors
            dec  (j,i) = 0        ! echarge
            dev  (j,i) = 0.0_re_p ! ehal1
         end do
      end do
      if (use_geom) call addForcesRestrain
      if (use_smd_velconst .or. use_smd_forconst) call addForcesSMD
      end subroutine

      module subroutine addForcesCharmm
      implicit none
      integer i,j

      if (deb_Path) print*, 'addForcesCharmm'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            desum(j,i) = dea(j,i) + deb(j,i) + deid(j,i) + deub(j,i)
     &                 + det(j,i) + mdr2md(dec(j,i)) + dev(j,i)
         end do
      end do
      if (use_geom) call addForcesRestrain
      if (use_smd_velconst .or. use_smd_forconst) call addForcesSMD
      end subroutine

      module subroutine resetForcesSMD
      implicit none
      integer i,j
      if (deb_Path) print*, 'resetForcesSMD'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
           desmd(j,i) = 0.0_re_p
         end do
      end do
      end subroutine
      module subroutine addForcesSMD
      implicit none
      integer i,j
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            desum(j,i) = desum(j,i) + desmd(j,i)
         end do
      end do
      end subroutine

      subroutine resetForcesRestrain
      implicit none
      integer i,j
      if (deb_Path) print*,'resetForcesRestrain'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            deg(j,i) = 0.0_re_p
         end do
      end do
      end subroutine
      subroutine addForcesRestrain
      implicit none
      integer i,j
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
c           debond(j,i) = debond(j,i) + deg(j,i)
            desum (j,i) = desum (j,i) + deg(j,i)
         end do
      end do
      end subroutine

      module subroutine resetForcesAMD
      implicit none
      integer i,j
      if (deb_Path) print*,'resetForcesAMD'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
           deamdD(j,i) = 0.0_re_p
           deamdP(j,i) = 0.0_re_p
          deW1aMD(j,i) = 0.0_re_p
         end do
      end do
      end subroutine

      module subroutine resetForcesDefault
      implicit none
      integer i,j
      if (deb_Path) print*,'resetForcesDefault'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            deb  (j,i) = 0.0_re_p ! ebond
            dea  (j,i) = 0.0_re_p ! eangle
            deba (j,i) = 0.0_re_p ! estrbnd
            deub (j,i) = 0.0_re_p ! eurey
            deaa (j,i) = 0.0_re_p ! eangang
            deopb(j,i) = 0.0_re_p ! eopbend
            deopd(j,i) = 0.0_re_p ! eopdist
            deid (j,i) = 0.0_re_p ! eimprop
            deit (j,i) = 0.0_re_p ! eimptor
            det  (j,i) = 0.0_re_p ! etors
            dept (j,i) = 0.0_re_p ! epitors
            debt (j,i) = 0.0_re_p ! estrtor
            dett (j,i) = 0.0_re_p ! etortor
            dev  (j,i) = 0.0_re_p ! ehal1
            dec  (j,i) = 0        ! echarge
            dem  (j,i) = 0        ! empole
            dep  (j,i) = 0        ! epolar
            deg  (j,i) = 0.0_re_p ! egeom
            desmd(j,i) = 0.0_re_p ! esmd
            dex  (j,i) = 0.0_re_p ! extra
         end do
      end do
      end subroutine

      module subroutine addForcesDefault
      implicit none
      integer i,j

      if (deb_Path) print*, 'addForcesDefault'
!$acc parallel loop collapse(2) default(present) async
      do i = 1, nbloc
         do j = 1, 3
            desum(j,i) = deaa(j,i) + deba(j,i) + dea  (j,i) + deb  (j,i)
     &         + mdr2md(dec (j,i) + dem (j,i) + dep (j,i)) + deopb(j,i)
     &                 + deid(j,i) + dev (j,i) + deub (j,i) + deopd(j,i)
     &                 + deit(j,i) + det (j,i) + dept (j,i) + debt (j,i)
     &                 + dett(j,i) + deg (j,i) + dex  (j,i) + desmd(j,i)

c          debond(j,i) = deaa (j,i) + deba (j,i) + dea (j,i) + deb (j,i)
c    &                 + deopb(j,i) + deopd(j,i) + deid(j,i) + deub(j,i)
c    &                 + deit (j,i) + det  (j,i) + dept(j,i) + debt(j,i)
c    &                 + dett (j,i) + deg  (j,i)
         end do
      end do
      end subroutine

      subroutine info_Forces_one( array,name )
      implicit none
      character(*),intent(in) :: name
      real(r_p),intent(inout):: array(:,:)
      integer i,j
      real(r_p) mini,maxi
      real(8) norm_l1,temp

      mini = huge(mini)
      maxi = tiny(maxi)
      norm_l1 = 0.0d0
!$acc data present(array)

      call commforce_one(array)

!$acc wait
!$acc parallel loop
      do i = 1,nloc
         do j = 1,3
            mini = min( mini,array(j,i) )
            maxi = max( maxi,array(j,i) )
            norm_l1 = norm_l1 + abs(array(j,i))
         end do
      end do

!$acc end data
      if (rank.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,mini,1,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,maxi,1,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,norm_l1,1,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      else
         call MPI_REDUCE(mini,mini,1,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(maxi,maxi,1,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(norm_l1,norm_l1,1,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      end if

 30   format(a7,a3,3F20.8)
      if (rank.eq.0) print 30,name,'=>',mini,maxi,norm_l1
      end subroutine

      module subroutine prtEForces(des,etot)
      use atoms  ,only: n
      use domdec
      use argue
      implicit none
      integer i, freeunit, iunit
      character(30) Ffile
      real(r_p) des(:,:)
      real(r_p) etot

      if (nproc.ne.1) then
         print '(10x,A)', "WARNING prtForces routine !!!"
         print '(10x,A)', "Not Operational in parallel execution " 
         return
      end if
      print*, 'prtForces'

!$acc wait
!$acc update host(des,etot)
      iunit = freeunit()
#if TINKER_SINGLE_PREC
      Ffile = trim(arg(1))//"_fs.txt"
#elif TINKER_MIXED_PREC
#  ifdef USE_DETERMINISTIC_REDUCTION
      Ffile = trim(arg(1))//"_ff.txt"
#  else
      Ffile = trim(arg(1))//"_fm.txt"
#  endif
#else
      Ffile = trim(arg(1))//"_fd.txt"
#endif
      call version(Ffile,'new')
      open( unit=iunit,file=Ffile,status='new' )
 12   format(3F18.10)
 13   format(F30.10)
      write (iunit,'(I)') n
      write (iunit,13 ) etot

      do i = 1,n
         write(iunit,12) des(1:3,i)
      end do

      close(iunit)
      call fatal
      end subroutine


      ! Debug routine on Forces
      module subroutine info_forces(rule)
      !use atoms
      implicit none
      integer,intent(in) :: rule
      integer,parameter:: nf=24  !number of forces
      integer i,j,sze
      real(r_p),dimension(3,nloc)::demsum,depsum,decsum
      real(8) mmx(3*nf),smm(nf)
      logical tinker_isnan_m
      real(t_p) dat,dat1
      integer,save:: doin = 1

      enum,bind(C)
      enumerator ::commBonded
      enumerator ::commShortNonBonded
      enumerator ::commNonBonded
      end enum

      mmx    = 0
      sze    = 3*nloc

      if (deb_Path) write(*,*) 'info_forces', rule

      !TODO Add desave and desmd to debug functions
      !TODO Reduce communication amount in this function
!$acc wait
      if (btest(rule,commBonded)) then

      if(use_bond)   call commforce_one(deb  ,commBonded)
      if(use_angle)  call commforce_one(dea  ,commBonded)
      if(use_strbnd) call commforce_one(deba ,commBonded)
      if(use_urey)   call commforce_one(deub ,commBonded)
      if(use_angang) call commforce_one(deaa ,commBonded)
      if(use_improp) call commforce_one(deid ,commBonded)
      if(use_imptor) call commforce_one(deit ,commBonded)
      if(use_tors)   call commforce_one(det  ,commBonded)
      if(use_pitors) call commforce_one(dept ,commBonded)
      if(use_strtor) call commforce_one(debt ,commBonded)
      if(use_tortor) call commforce_one(dett ,commBonded)
      if(use_opbend) call commforce_one(deopb,commBonded)
      if(use_opdist) call commforce_one(deopd,commBonded)
      if(use_geom)   call commforce_one(deg  ,commBonded)
      if(use_extra)  call commforce_one(dex  ,commBonded)

      end if

      if (btest(rule,commShortNonBonded)) then

      if(use_charge) call commforce_one(dec ,commShortNonBonded)
      if(use_vdw)    call commforce_one(dev ,commShortNonBonded)
      if(use_mpole)  call commforce_one(dem ,commShortNonBonded)
      if(use_polar)  call commforce_one(dep ,commShortNonBonded)

      end if

      if (btest(rule,commNonBonded)) then

!$acc enter data create(depsum,decsum,demsum)

!$acc parallel loop default(present)
      do i = 1, nloc
         do j = 1, 3
            demsum(j,i) = real(0.0,r_p)
            depsum(j,i) = real(0.0,r_p)
            decsum(j,i) = real(0.0,r_p)
         end do
      end do


      if(use_vdw)    call commforce_one(dev,commNonBonded)

#ifndef USE_DETERMINISTIC_REDUCTION
      if(use_charge) call commforce_one(dec,commNonBonded)
      if(use_charge) call commforcerecdir_one(dec,decrec,decsum)

      if(use_mpole)  call commforce_one(dem,commNonBonded)
      if(use_mpole)  call commforcerecdir_one(dem,demrec,demsum)

      if(use_polar)  call commforce_one(dep,commNonBonded)
      if(use_polar)  call commforcerecdir_one(dep,deprec,depsum)
#endif

      end if
!$acc wait

c      call updateForcesHost
c!$acc update host(glob)

      if (btest(rule,commBonded)) then

      if(use_bond)
     &call minmaxone(mmx(01),mmx(nf+01),mmx(2*nf+01),deb
     &     ,3*nloc,'deb');
      if(use_angle)
     &call minmaxone(mmx(02),mmx(nf+02),mmx(2*nf+02),dea
     &     ,3*nloc,'dea');
      if(use_strbnd)
     &call minmaxone(mmx(03),mmx(nf+03),mmx(2*nf+03),deba
     &     ,3*nloc,'deba');
      if(use_urey)
     &call minmaxone(mmx(04),mmx(nf+04),mmx(2*nf+04),deub
     &     ,3*nloc,'deub');
      if(use_angang)
     &call minmaxone(mmx(05),mmx(nf+05),mmx(2*nf+05),deaa
     &     ,3*nloc,'deaa');
      if(use_improp)
     &call minmaxone(mmx(06),mmx(nf+06),mmx(2*nf+06),deid
     &     ,sze,'deid');
      if(use_imptor)
     &call minmaxone(mmx(07),mmx(nf+07),mmx(2*nf+07),deit
     &     ,sze,'deit');
      if(use_tors)
     &call minmaxone(mmx(08),mmx(nf+08),mmx(2*nf+08),det
     &     ,sze,'det');
      if(use_pitors)
     &call minmaxone(mmx(09),mmx(nf+09),mmx(2*nf+09),dept
     &     ,sze,'dept');
      if(use_strtor)
     &call minmaxone(mmx(10),mmx(nf+10),mmx(2*nf+10),debt
     &     ,sze,'debt');
      if(use_tortor)
     &call minmaxone(mmx(11),mmx(nf+11),mmx(2*nf+11),dett
     &     ,sze,'dett');
      if(use_opbend)
     &call minmaxone(mmx(12),mmx(nf+12),mmx(2*nf+12),deopb
     &     ,sze,'deopb');
      if(use_opdist)
     &call minmaxone(mmx(13),mmx(nf+13),mmx(2*nf+13),deopd
     &     ,sze,'deopd');
      if(use_geom)
     &call minmaxone(mmx(14),mmx(nf+14),mmx(2*nf+14),deg
     &     ,sze);
      if(use_extra)
     &call minmaxone(mmx(15),mmx(nf+15),mmx(2*nf+15),dex
     &     ,sze,'dex');

      end if

      if (btest(rule,commShortNonBonded)) then

      if(use_charge)
     &call minmaxone1(mmx(16),mmx(nf+16),mmx(2*nf+16),dec
     &     ,sze,'dec');
      if(use_vdw)
     &call minmaxone(mmx(17),mmx(nf+17),mmx(2*nf+17),dev
     &     ,sze,'dev');
      if(use_mpole)
     &call minmaxone1(mmx(18),mmx(nf+18),mmx(2*nf+18),dem
     &     ,sze,'dem');
      if(use_polar)
     &call minmaxone1(mmx(20),mmx(nf+20),mmx(2*nf+20),dep
     &     ,sze,'dep');

c     if (abort.and.vlst_enable) then
c        dint1 = minval(inte); dint2=maxval(inte)
c        call searchpair(nshortvlst,shortvlst,maxvlst
c    &        ,dint1,dint2)
c     end if

      else if (btest(rule,commNonBonded)) then

      if(use_charge) then
      call minmaxone1(mmx(16),mmx(nf+16),mmx(2*nf+16),dec
     &     ,sze,'dec');
      call minmaxone(mmx(22),mmx(nf+22),mmx(2*nf+22),decsum
     &     ,sze,'decsum');
      end if
      if(use_vdw)
     &call minmaxone(mmx(17),mmx(nf+17),mmx(2*nf+17),dev
     &     ,sze,'dev');
      if(use_mpole) then
      call minmaxone1(mmx(18),mmx(nf+18),mmx(2*nf+18),dem
     &     ,sze,'dem');
      call minmaxone(mmx(19),mmx(nf+19),mmx(2*nf+19),demsum
     &     ,sze,'demsum');
      end if
      if(use_polar) then
      call minmaxone1(mmx(20),mmx(nf+20),mmx(2*nf+20),dep
     &     ,sze,'dep');
      call minmaxone(mmx(21),mmx(nf+21),mmx(2*nf+21),depsum
     &     ,sze,'depsum');
      end if

c     if (abort.and.vlst_enable) then
c        dint1 = minval(inte); dint2=maxval(inte)
c        call searchpair(nshortvlst,shortvlst,maxvlst,
c    &        dint1,dint2)
c     end if

      end if

c     if (allocated(deamdD)) then
c     call minmaxone(mmx(23),mmx(nf+23),mmx(2*nf+23),deamdD ,3*nloc);
c     call minmaxone(mmx(24),mmx(nf+24),mmx(2*nf+24),deW1aMD,3*nloc);
c     end if

      if (nproc.gt.1) then
      if (rank.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,mmx,nf,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,mmx(nf+1),nf,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(MPI_IN_PLACE,mmx(2*nf+1),nf,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      else
         call MPI_REDUCE(mmx,mmx,nf,MPI_RPREC,
     &                   MPI_MIN,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(mmx(nf+1),mmx(nf+1),nf,MPI_RPREC,
     &                   MPI_MAX,0,MPI_COMM_WORLD,i)
         call MPI_REDUCE(mmx(2*nf+1),mmx(2*nf+1),nf,MPI_RPREC,
     &                   MPI_SUM,0,MPI_COMM_WORLD,i)
      end if
      end if

c     if (doin>4) then
c35   format(I8,I8,3F16.6)
c        do i = 1, nloc
c           print 35, i, glob(i), dev(:,i)
c        end do
c     end if

      if (rank.eq.0) then
 30      format(a10,3F20.8)
 40      format(80('='))
         print 40
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
         if(mmx(2*nf+12)/=0.0) print 30,'deopb  =>',
     &      mmx(12),mmx(12+nf),mmx(12+nf*2)
         if(mmx(2*nf+13)/=0.0) print 30,'deopd  =>',
     &      mmx(13),mmx(13+nf),mmx(13+nf*2)
         if(mmx(2*nf+14)/=0.0) print 30,'deg    =>',
     &      mmx(14),mmx(14+nf),mmx(14+nf*2)
         if(mmx(2*nf+15)/=0.0) print 30,'dex    =>',
     &      mmx(15),mmx(15+nf),mmx(15+nf*2)
         if(mmx(2*nf+16)/=0.0) print 30,'dec    =>',
     &      mmx(16),mmx(16+nf),mmx(16+nf*2)
         if(mmx(2*nf+17)/=0.0) print 30,'dev    =>',
     &      mmx(17),mmx(17+nf),mmx(17+nf*2)
         if(mmx(2*nf+18)/=0.0) print 30,'dem    =>',
     &      mmx(18),mmx(18+nf),mmx(18+nf*2)
         if(mmx(2*nf+19)/=0.0) print 30,'demsum =>',
     &      mmx(19),mmx(19+nf),mmx(19+nf*2)
         if(mmx(2*nf+20)/=0.0) print 30,'dep    =>',
     &      mmx(20),mmx(20+nf),mmx(20+nf*2)
         if(mmx(2*nf+21)/=0.0) print 30,'depsum =>',
     &      mmx(21),mmx(21+nf),mmx(21+nf*2)
         if(mmx(2*nf+22)/=0.0) print 30,'decsum =>',
     &      mmx(22),mmx(22+nf),mmx(22+nf*2)
         if(mmx(2*nf+23)/=0.0) print 30,'deamdD =>',
     &      mmx(23),mmx(23+nf),mmx(23+nf*2)
         if(mmx(2*nf+24)/=0.0) print 30,'deW1aMD =>',
     &      mmx(24),mmx(24+nf),mmx(24+nf*2)
      end if

      if (btest(rule,commNonBonded)) then

!$acc exit data delete(decsum,depsum,demsum)

      if (use_mpole.and.use_polar) then
!$acc    parallel loop default(present) async
         do i = 1, size(demrec)
            demrec(i,1) = 0
            deprec(i,1) = 0
         end do
      end if
      if (use_charge) then
!$acc    parallel loop default(present) async
         do i = 1, size(deprec)
            decrec(i,1) = 0
         end do
      end if
      if (use_mpole.and..not.use_polar) then
!$acc    parallel loop default(present) async
         do i = 1, size(decrec)
            demrec(i,1) = 0
         end do
      end if

      end if

      doin = doin +1
      end subroutine

      subroutine minmaxone( mi,ma,on,vector,sz,name )
      implicit none
      integer sz
      real(8) mi,ma,on
      real(r_p) vector(*)
      character(*),optional,intent(in)::name
      integer i,j,i1,i2,iglob,cap,cap1,cap2
      integer gli(2,nproc)
      real(8) val

      abortall = .false.
      if (present(name)) then
         cap = 1
         gli = 0
         if (name.eq.'devi') then
!$acc parallel loop default(present) copy(abort,cap,gli)
            do i = 1, sz/3
               cap2  = 0
               iglob = glob(i)
!$acc loop seq
               do j  = 1,3
               val   = vector((i-1)*3+j)
               if (abs(val).gt.300.0_re_p) then
                  print*,i,iglob,rank,x(iglob),y(iglob),z(iglob),val
!$acc atomic write
                  abort=.true.
                  cap2 = cap2 + 1
               end if
               end do
               if (cap2.gt.0) then
!$acc atomic capture
               cap1 = cap
               cap  = cap + 1
!$acc end atomic
               if (cap1.le.2) gli(cap1,rank+1) = iglob
               end if
            end do
            do i = 1, nproc
               if (rank.eq.i-1) abortall=abort
               call MPI_BCAST(abortall,1,MPI_LOGICAL,i-1,COMM_TINKER,i1)
               if (abortall) then
                  abort = .true.
                  call MPI_AllGather(MPI_IN_PLACE,2,MPI_DATATYPE_NULL
     $                              ,gli,2,MPI_INT,COMM_TINKER,i1)
                  exit
               end if
            end do

            if (abort) then
            
            cap  = 0
            cap1 = 0
            do i1 = 1, nproc
               do i2 = 1,2
                  if ( gli(i2,i1).ne.0 ) then
                     cap = cap + 1
                     if (cap.lt.3) inte(cap) = gli(i2,i1)
                  else
                     cap1 = cap1 + 1
                  end if
               end do
            end do
            
            if (cap.ne.2) print*,' more than one interactions find '
     &                   ,rank,cap1
            !print*, 'interaction', inte,rank
            
            end if
         end if

      end if

!$acc parallel loop present(vector(1:sz))
      do i = 1, sz
         val = vector(i)
         mi = min( mi,val )
         ma = max( ma,val )
         on = on + abs(val)
      end do
      end subroutine

      subroutine minmaxone1( mi,ma,on,vector,sz,name )
      implicit none
      integer sz
      real(8) mi,ma,on
      mdyn_rtyp vector(*)
      character(*),optional,intent(in)::name
      integer i,j
      real(8) val

      abortall = .false.

!$acc parallel loop present(vector(1:sz))
      do i = 1, sz
         val = mdr2md(vector(i))
         mi = min( mi,val )
         ma = max( ma,val )
         on = on + abs(val)
      end do
      end subroutine
      end submodule
