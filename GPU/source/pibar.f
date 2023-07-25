c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  program pibar  --  thermodynamic values from simulation data  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "pibar" computes the free energy, enthalpy and entropy difference
c     between two states via Zwanzig free energy perturbation (FEP)
c     and Bennett acceptance ratio (BAR) methods for PI simulations
c
c     note current version takes as input the trajectory archives A
c     and B, originally generated for states 0 and 1, respectively;
c     then finds the total potential energy for all frames of each
c     trajectory under control of key files for both states 0 and 1;
c     finally the FEP and BAR algorithms are used to compute the
c     free energy for state 0 --> state 1, and similar estimators
c     provide the enthalpy and entropy for state 0 --> state 1
c
c     modifications for NPT simulations by Chengwen Liu, University
c     of Texas at Austin, October 2015; enthalpy and entropy methods
c     by Aaron Gordon, Washington University, December 2016
c
c     literature references:
c
c     C. H. Bennett, "Efficient Estimation of Free Energy Differences
c     from Monte Carlo Data", Journal of Computational Physics, 22,
c     245-268 (1976)
c
c     K. B. Daly, J. B. Benziger, P. G. Debenedetti and
c     A. Z. Panagiotopoulos, "Massively Parallel Chemical Potential
c     Calculation on Graphics Processing Units", Computer Physics
c     Communications, 183, 2054-2062 (2012)  [modification for NPT]
c
c     M. A. Wyczalkowski, A. Vitalis and R. V. Pappu, "New Estimators
c     for Calculating Solvation Entropy and Enthalpy and Comparative
c     Assessments of Their Accuracy and Precision, Journal of Physical
c     Chemistry, 114, 8166-8180 (2010)  [entropy and enthalpy]
c
c
#include "tinker_precision.h"
#include "tinker_types.h"

      program pibar
      use mpi
      implicit none
      integer ierr,nthreadsupport
c      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
      call MPI_INIT(ierr)
      call pibar_bis
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      end
c
      subroutine pibar_bis
      use iounit
      use inform
      use domdec
      use beads
      use files
      implicit none
      integer mode
      logical exist,query
      character*240 record
      character*240 string
      character*20 keyword
c
      ! Sign running program
      app_id = pibar_a

      call initial
      nproc=nproctot
      call initmpi
c
c     find thermodynamic perturbation procedure to perform
c
      mode  = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         query = .false.
      end if
   10 continue
      if (query) then
         if (ranktot.eq.0) write (iout,20)
   20    format (/,' The TINKER Thermodynamic Perturbation Utility'
     &            ,' Can :'
     &         ,//,4x,'(1) Create BAR File with Perturbed Potential'
     &            ,   ' Energies'
     &          ,/,4x,'(2) Compute Thermodynamic Values from TINKER'
     &            ,   ' BAR File'
     &          )
         do while (mode.lt.1 .or. mode.gt.2)
            mode = 0
            if (ranktot.eq.0) write (iout,'(a,$)')
     &         ' Enter the Number of the Desired Choice :  '
            read (input,40,err=50,end=50) mode
   40       format (i10)
   50       continue
         end do
      end if
c
c     create TINKER BAR file or compute thermodynamic values
c
      if (mode .eq. 1)  call makebarpi
      if (mode .eq. 2)  call barcalcpi
      call final
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine makebar  --  store trajectory potential energy  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine makebarpi
      use atoms
      use boxes
      use dcdmod
      use domdec
      use files
      use energi ,only: info_energy,calc_e
      use inform
      use iounit
      use keys
      use mutant
      use titles
      use tinheader
      use beads
      use virial
      use mpi
      use utilbaoabpi
      use commstuffpi
      use cutoff
      implicit none
      integer i,j,iarc,ibar,ierr,ibead
      integer lenga,lengb
      integer ltitlea,ltitleb
      integer nkey0,nkey1
      integer nfrma,nfrmb
      integer maxframe
      integer freeunit
      integer trimtext
      real(en_p) energy
      real(en_p) tempa,tempb
      real(en_p), allocatable :: ua0(:)
      real(en_p), allocatable :: ua1(:)
      real(en_p), allocatable :: ub0(:)
      real(en_p), allocatable :: ub1(:)
      real(en_p), allocatable :: ua0_ctr(:)
      real(en_p), allocatable :: ua1_ctr(:)
      real(en_p), allocatable :: ub0_ctr(:)
      real(en_p), allocatable :: ub1_ctr(:)
      real(r_p) , allocatable :: vola(:)
      real(r_p) , allocatable :: volb(:)
      real(en_p) epot,epot_ctr
      logical exist, only_long
      character*240 string
      character*240 exten
      character*240 fnamea
      character*240 fnameb
      character*240 titlea
      character*240 titleb
      character*240 arcfile
      character*240 dcdfile
      character*240 barfile
      character*240, allocatable :: keys0(:)
      character*240, allocatable :: keys1(:)
      logical :: contract_a, contract_b
      type(dcdinfo_t) :: dcdinfo,dcdinfoB
      character*3 numberbeads
      integer :: nbeads_old,nbeads_ctr_old

      calc_e = .true.
      use_virial = .false.
c
c
c     get trajectory A archive and setup mechanics calculation
c
c      call initial
      call init_keys
      call read_bead_config

      call getxyz
      call cutoffs
      call unitcell
      call lattice
c
c     setup for MPI
c
      call drivermpi

      call reset_observables(.true.)
      !call allocpi()  
c     call reinitnl(0)

c     call mechanic
c     call nblist(0)
c
c     find the original temperature value for trajectory A
c
      tempa = En_p(-1.0)
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  tempa
   10 continue
      do while (tempa .lt. En_p(0.0))
         if (ranktot.eq.0) write (iout,'(a)', advance='no')
     &      ' Enter Trajectory A Temperature in Degrees K [298] : '
         read (input,30,err=40) tempa
   30    format (f20.0)
         if (tempa .le. En_p(0.0))  tempa = En_p(298.0)
   40    continue
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys0(maxkey))
      allocate (keys1(maxkey))
c
c     store filename and keywords for trajectory A and state 0
c
      fnamea  = filename
      lenga   = leng
      titlea  = title
      ltitlea = ltitle
      nkey0   = nkey
      do i = 1, nkey0
         keys0(i) = keyline(i)
      end do
      contract_a = contract .or. centroid_longrange
c
c     get trajectory B archive and setup mechanics calculation
c
      keys_already_read=.false.
      call getxyz
      call read_bead_config

      !call ddpme3d
c     call AllDirAssign
c     call reassignpme(.true.)
c     call reinitnl(0)
c     call mechanic
      silent = .true.
c
c     find the original temperature value for trajectory B
c
      tempb = En_p(-1.0)
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  tempb
   50 continue
      do while (tempb .lt. En_p(0.0))
         if (ranktot.eq.0) write (iout,'(a)',advance='no')
     &      ' Enter Trajectory B Temperature in Degrees K [298] :  '
         read (input,70,err=80)  tempb
   70    format (f20.0)
         if (tempb .le. 0.0d0)  tempb = En_p(298.0)
   80    continue
      end do
c
c     store filename and keywords for trajectory B and state 1
c
      fnameb  = filename
      lengb   = leng
      titleb  = title
      ltitleb = ltitle
      nkey1   = nkey
      do i = 1, nkey1
         keys1(i) = keyline(i)
      end do
      contract_b = contract .or. centroid_longrange
c
c     perform dynamic allocation of some local arrays
c
      maxframe = 1000000
      allocate (ua0(maxframe))
      allocate (ua1(maxframe))
      allocate (ub0(maxframe))
      allocate (ub1(maxframe))
      if(contract_a .or. contract_b) then
        allocate (ua0_ctr(maxframe))
        allocate (ua1_ctr(maxframe))
        allocate (ub0_ctr(maxframe))
        allocate (ub1_ctr(maxframe))
      endif
      allocate (vola(maxframe))
      allocate (volb(maxframe))
c
c     reopen trajectory A using the parameters for state 0
c
      nkey = nkey0
      do i = 1, nkey
        keyline(i) = keys0(i)
      end do
      !set polymer for trajectory A
      call read_bead_config
      call allocpi
      
      if(.not. allocated(polymer%dcdinfo)) then
        allocate(polymer%dcdinfo(polymer%nbeads))
      endif
!$acc parallel loop async default(present) collapse(2)
      DO i=1,n; do j=1,3
        polymer%eigpos(j,i,1)=0.d0
      ENDDO; ENDDO
!$acc wait
      do ibead=1,nbeads
        write(numberbeads, '(i3.3)') ibead
        exten='_beads'//numberbeads
        if (dcdio) then
          dcdfile = fnamea(1:lenga)//trim(exten)//'.dcd'
          call dcdfile_open(polymer%dcdinfo(ibead),dcdfile)
          call dcdfile_read_header(polymer%dcdinfo(ibead),.false.
     &       ,verbose=.false.)
          call dcdfile_read_next(polymer%dcdinfo(ibead))
          call dcdfile_skip_next(polymer%dcdinfo(ibead),0)
        else
          iarc    = freeunit ()
          polymer%dcdinfo(ibead)%idcd = iarc
          arcfile = fnamea(1:lenga)//trim(exten)
          call suffix (arcfile,'arc','old')
          open (unit=iarc,file=arcfile,status ='old')
          rewind (unit=iarc)
          call readxyz (iarc)
        end if
!$acc wait
!$acc parallel loop default(present)
        DO i=1,n
          polymer%pos(1,i,ibead)=x(i)
          polymer%pos(2,i,ibead)=y(i)
          polymer%pos(3,i,ibead)=z(i)
          polymer%eigpos(1,i,1)=polymer%eigpos(1,i,1)+x(i)/nbeads
          polymer%eigpos(2,i,1)=polymer%eigpos(2,i,1)+y(i)/nbeads
          polymer%eigpos(3,i,1)=polymer%eigpos(3,i,1)+z(i)/nbeads
        ENDDO
      enddo
      
      
      !call ddpme3d
c
c     the box shape can change between frames
c
      call load_bead(polymer, 0,.true.)
      call lattice
      call AllDirAssign
      call reinitnl(0)
      call mechanic
      call nblist(0)
      if (abort) then
         if (ranktot.eq.0) write (iout,90)
   90    format (/,' BAR  --  No Coordinate Frames Available',
     &              ' from First Input File')
         call fatal
      end if
c
c     find potential energies for trajectory A in state 0
c
      if (ranktot.eq.0) write (iout,100)
  100 format (/,' Initial Processing for Trajectory A :',/)
      i = 0
      maina1: do while (.not. abort)
         i = i + 1
         call cutoffs
         call compute_energy_beads(epot,polymer
     &     ,.TRUE.,.TRUE.,.TRUE.,.TRUE.)
         if(ranktot==0) then
          call MPI_REDUCE(MPI_IN_PLACE,epot
     &       ,1,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         else
          call MPI_REDUCE(epot,epot
     &        ,1,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         endif

         ua0(i)  = epot
         vola(i) = volbox

          if(centroid_recip .or. contract) then
            only_long = centroid_longrange .and.
     &         (.not. (contract.and.nbeads_ctr==1) )
            if(only_long) then
              use_shortmlist = use_mlist
              use_shortclist = use_clist
              use_shortvlist = use_vlist
              call load_bead(polymer, 0,.true.)
              call reinitnl(0)
              call mechanicstep(0)
              call nblist(0)
            endif
            epot_ctr=0.d0
            if(centroid_recip) then
             call load_bead(polymer,0)
             call compute_eslow_centroid(epot
     &          ,only_long, polar_allbeads)
             epot_ctr = epot_ctr + epot
            endif
   
            call compute_energy_beads(epot,polymer
     &            ,.not. (centroid_longrange .or. contract)  !LONG 
     &            ,.not. contract ! INT
     &            , .TRUE.     !SHORT 
     &            ,polar_allbeads)
            epot_ctr = epot_ctr + epot

            if(contract .and. nbeads_ctr>1) then
              call contract_polymer(polymer,polymer_ctr)
              call comm_for_gradient(polymer_ctr,.FALSE.,.FALSE.)

              call compute_energy_beads(epot,polymer_ctr
     &       ,.not.centroid_longrange, .TRUE., .FALSE.  ![LONG, INT  , SHORT]
     &       ,polar_allbeads)
              epot_ctr = epot_ctr + epot
            endif

            if(ranktot==0) then
              call MPI_REDUCE(MPI_IN_PLACE,epot_ctr
     &          ,1,MPI_RPREC
     &          ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            else
              call MPI_REDUCE(epot_ctr,epot_ctr
     &           ,1,MPI_RPREC
     &          ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            endif

            ua0_ctr(i) = epot_ctr
            if (verbose.and.ranktot==0) then
            write(iout,'(i11,2x,2f16.4)')  i,ua0(i), ua0_ctr(i)
            endif
         else
            if (verbose.and.ranktot==0) then
            write(iout,'(i11,2x,1f16.4)') i,ua0(i)
            endif
         endif

!$acc parallel loop async default(present) collapse(2)
         DO i=1,n; do j=1,3
           polymer%eigpos(j,i,1)=0.d0
         ENDDO; ENDDO
!$acc wait
         do ibead=1,nbeads
            if (dcdio) then
              call dcdfile_read_next(polymer%dcdinfo(ibead))
              call dcdfile_skip_next(polymer%dcdinfo(ibead),0)
              if (abort) cycle maina1
            else
              call readxyz (polymer%dcdinfo(ibead)%idcd)
            end if
!$acc wait
!$acc parallel loop default(present)
            DO i=1,n
              polymer%pos(1,i,ibead)=x(i)
              polymer%pos(2,i,ibead)=y(i)
              polymer%pos(3,i,ibead)=z(i)
              polymer%eigpos(1,i,1)=polymer%eigpos(1,i,1)+x(i)/nbeads
              polymer%eigpos(2,i,1)=polymer%eigpos(2,i,1)+y(i)/nbeads
              polymer%eigpos(3,i,1)=polymer%eigpos(3,i,1)+z(i)/nbeads
            ENDDO
         enddo
         if (n.eq.0) cycle
c
c        the box shape can change between frames
c
         use_shortmlist = .false.
         use_shortclist = .false.
         use_shortvlist = .false.
         call load_bead(polymer, 0,.true.)
         call lattice
         !call ddpme3d
         call AllDirAssign
         call AllRecAssign
         call reinitnl(0)
         call mechanicstep(0)
         call nblist(0)
         if (i .ge. maxframe)  abort = .true.
         if (mod(i,100).eq.0 .or. abort) then
            if (ranktot.eq.0) write (iout,110)  i
  110       format (7x,'Completed',i8,' Coordinate Frames')
            flush (iout)
         end if
      end do maina1
c
c     reset trajectory A using the parameters for state 1
c
!$acc parallel loop async default(present) collapse(2)
      DO i=1,n; do j=1,3
        polymer%eigpos(j,i,1)=0.d0
      ENDDO; ENDDO
!$acc wait
      do ibead=1,nbeads
        write(numberbeads, '(i3.3)') ibead
        exten='_beads'//numberbeads
        if (dcdio) then
          call dcdfile_close(polymer%dcdinfo(ibead))
          dcdfile = fnamea(1:lenga)//trim(exten)//'.dcd'
          call dcdfile_open(polymer%dcdinfo(ibead),dcdfile)
          call dcdfile_read_header(polymer%dcdinfo(ibead),.false.
     &       ,verbose=.false.)
          call dcdfile_read_next(polymer%dcdinfo(ibead))
          call dcdfile_skip_next(polymer%dcdinfo(ibead),0)
          abort = .false.
        else
          rewind (unit=polymer%dcdinfo(ibead)%idcd)
          call readxyz (polymer%dcdinfo(ibead)%idcd)
        end if
!$acc wait
!$acc parallel loop default(present)
        DO i=1,n
          polymer%pos(1,i,ibead)=x(i)
          polymer%pos(2,i,ibead)=y(i)
          polymer%pos(3,i,ibead)=z(i)
          polymer%eigpos(1,i,1)=polymer%eigpos(1,i,1)+x(i)/nbeads
          polymer%eigpos(2,i,1)=polymer%eigpos(2,i,1)+y(i)/nbeads
          polymer%eigpos(3,i,1)=polymer%eigpos(3,i,1)+z(i)/nbeads
        ENDDO
      enddo
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
c
c     the box shape can change between frames
c
      call load_bead(polymer, 0,.true.)
      call lattice
      call AllDirAssign
      call reinitnl(0)
      call mechanic
      call nblist(0)
c
c     find potential energies for trajectory A in state 1
c
      if (verbose) then
         if (ranktot.eq.0) write (iout,120)
  120    format (/,' Potential Energy Values for Trajectory A :',
     &           //,7x,'Frame',9x,'State 0',9x,'State 1',12x,'Delta',/)
      end if
      i = 0
      maina2: do while (.not. abort)
         i = i + 1
         call cutoffs
         call compute_energy_beads(epot,polymer
     &     ,.TRUE.,.TRUE.,.TRUE.,.TRUE.)
         if(ranktot==0) then
          call MPI_REDUCE(MPI_IN_PLACE,epot
     &       ,1,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         else
          call MPI_REDUCE(epot,epot
     &        ,1,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         endif

         ua1(i) = epot

          if(centroid_recip .or. contract) then
            only_long = centroid_longrange .and.
     &         (.not. (contract.and.nbeads_ctr==1) )
            if(only_long) then
              use_shortmlist = use_mlist
              use_shortclist = use_clist
              use_shortvlist = use_vlist
              call load_bead(polymer, 0,.true.)
              call reinitnl(0)
              call mechanicstep(0)
              call nblist(0)
            endif
            epot_ctr=0.d0
            if(centroid_recip) then
             call load_bead(polymer,0)
             call compute_eslow_centroid(epot
     &          ,only_long, polar_allbeads)
             epot_ctr = epot_ctr + epot
            endif
   
            call compute_energy_beads(epot,polymer
     &            ,.not. (centroid_longrange .or. contract)  !LONG 
     &            ,.not. contract ! INT
     &            , .TRUE.     !SHORT 
     &            ,polar_allbeads)
            epot_ctr = epot_ctr + epot

            if(contract .and. nbeads_ctr>1) then
              call contract_polymer(polymer,polymer_ctr)
              call comm_for_gradient(polymer_ctr,.FALSE.,.FALSE.)

              call compute_energy_beads(epot,polymer_ctr
     &       ,.not.centroid_longrange, .TRUE., .FALSE.  ![LONG, INT  , SHORT]
     &       ,polar_allbeads)
              epot_ctr = epot_ctr + epot
            endif

            if(ranktot==0) then
              call MPI_REDUCE(MPI_IN_PLACE,epot_ctr
     &          ,1,MPI_RPREC
     &          ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            else
              call MPI_REDUCE(epot_ctr,epot_ctr
     &           ,1,MPI_RPREC
     &          ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            endif

            ua1_ctr(i) = epot_ctr
            if (verbose) then
               if (ranktot.eq.0) write (iout,'(i11,2x,6f16.4)')
     $              i,ua0(i),ua1(i),
     $              ua1(i)-ua0(i),ua0_ctr(i),ua1_ctr(i),
     $              ua1_ctr(i)-ua0_ctr(i)
            end if
         else
            if (verbose) then
               if (ranktot.eq.0) write (iout,'(i11,2x,3f16.4)')  
     &               i,ua0(i),ua1(i),
     $              ua1(i)-ua0(i)
            end if
         endif

!$acc parallel loop async default(present) collapse(2)
         DO i=1,n; do j=1,3
           polymer%eigpos(j,i,1)=0.d0
         ENDDO; ENDDO
!$acc wait
         do ibead=1,nbeads
            if (dcdio) then
              call dcdfile_read_next(polymer%dcdinfo(ibead))
              call dcdfile_skip_next(polymer%dcdinfo(ibead),0)
              if (abort) cycle maina2
            else
              call readxyz (polymer%dcdinfo(ibead)%idcd)
            end if
!$acc wait
!$acc parallel loop default(present)
            DO i=1,n
              polymer%pos(1,i,ibead)=x(i)
              polymer%pos(2,i,ibead)=y(i)
              polymer%pos(3,i,ibead)=z(i)
              polymer%eigpos(1,i,1)=polymer%eigpos(1,i,1)+x(i)/nbeads
              polymer%eigpos(2,i,1)=polymer%eigpos(2,i,1)+y(i)/nbeads
              polymer%eigpos(3,i,1)=polymer%eigpos(3,i,1)+z(i)/nbeads
            ENDDO
         enddo
         if (n.eq.0) cycle
c
c        the box shape can change between frames
c
         use_shortmlist = .false.
         use_shortclist = .false.
         use_shortvlist = .false.
         call load_bead(polymer, 0,.true.)
         call lattice
         !call ddpme3d
         call AllDirAssign
         call AllRecAssign
         call reinitnl(0)
         call mechanicstep(0)
         call nblist(0)
         if (i .ge. maxframe)  abort = .true.
      end do maina2
      nfrma = i
      do ibead=1,nbeads
        if (dcdio) then
          call dcdfile_close(polymer%dcdinfo(ibead))
        else
          close (unit=polymer%dcdinfo(ibead)%idcd)
        end if
      enddo
      !deallocate(polymer%dcdinfo)
c
c     reopen trajectory B using the parameters for state 0
c
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      !set polymer for trajectory B
      nbeads_old=nbeads
      nbeads_ctr_old=nbeads_ctr
      call read_bead_config
      if(nbeads_old.ne.nbeads) then
        write(iout,*) "ERROR: number of beads in trajectory"
     &      ," A and B must be the same"
        call fatal
      elseif(nbeads_ctr_old.ne.nbeads_ctr) then
        write(iout,*) "ERROR: number of contracted beads in trajectory"
     &      ," A and B must be the same"
        call fatal
      endif
      !call allocpi
!
      !allocate(polymer%dcdinfo(polymer%nbeads))

!$acc parallel loop async default(present) collapse(2)
      DO i=1,n; do j=1,3
        polymer%eigpos(j,i,1)=0.d0
      ENDDO; ENDDO
!$acc wait
      do ibead=1,nbeads
        write(numberbeads, '(i3.3)') ibead
        exten='_beads'//numberbeads
        if (dcdio) then
          dcdfile = fnameb(1:lengb)//trim(exten)//'.dcd'
          call dcdfile_open(polymer%dcdinfo(ibead),dcdfile)
          call dcdfile_read_header(polymer%dcdinfo(ibead),.false.
     &       ,verbose=.false.)
          call dcdfile_read_next(polymer%dcdinfo(ibead))
          call dcdfile_skip_next(polymer%dcdinfo(ibead),0)
          abort = .false.
        else
          iarc    = freeunit ()
          polymer%dcdinfo(ibead)%idcd = iarc
          arcfile = fnameb(1:lengb)//trim(exten)
          call suffix (arcfile,'arc','old')
          open (unit=iarc,file=arcfile,status ='old')
          rewind (unit=iarc)
          call readxyz (iarc)
        end if
!$acc wait
!$acc parallel loop default(present)
        DO i=1,n
          polymer%pos(1,i,ibead)=x(i)
          polymer%pos(2,i,ibead)=y(i)
          polymer%pos(3,i,ibead)=z(i)
          polymer%eigpos(1,i,1)=polymer%eigpos(1,i,1)+x(i)/nbeads
          polymer%eigpos(2,i,1)=polymer%eigpos(2,i,1)+y(i)/nbeads
          polymer%eigpos(3,i,1)=polymer%eigpos(3,i,1)+z(i)/nbeads
        ENDDO
      enddo

      if (abort) then
         if (ranktot.eq.0) write (iout,140)
  140    format (/,' BAR  --  No Coordinate Frames Available',
     &              ' from Second Input File')
         call fatal
      end if

      ! Reset State 0 parameters
      nkey = nkey0
      do i = 1, nkey
         keyline(i) = keys0(i)
      end do
      call load_bead(polymer, 0,.true.)
      call lattice
      call AllDirAssign
      call reinitnl(0)
      call mechanic
      call nblist(0)
      
c
c     find potential energies for trajectory B in state 0
c
      if (ranktot.eq.0) write (iout,150)
  150 format (/,' Initial Processing for Trajectory B :',/)
      i = 0
      mainb1: do while (.not. abort)
         i = i + 1
         call cutoffs
         call compute_energy_beads(epot,polymer
     &     ,.TRUE.,.TRUE.,.TRUE.,.TRUE.)
         if(ranktot==0) then
          call MPI_REDUCE(MPI_IN_PLACE,epot
     &       ,1,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         else
          call MPI_REDUCE(epot,epot
     &        ,1,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         endif
 
         ub0 (i) = epot
         volb(i) = volbox

          if(centroid_recip .or. contract) then
            only_long = centroid_longrange .and.
     &         (.not. (contract.and.nbeads_ctr==1) )
            if(only_long) then
              use_shortmlist = use_mlist
              use_shortclist = use_clist
              use_shortvlist = use_vlist
              call load_bead(polymer, 0,.true.)
              call reinitnl(0)
              call mechanicstep(0)
              call nblist(0)
            endif
            epot_ctr=0.d0
            if(centroid_recip) then
             call load_bead(polymer,0)
             call compute_eslow_centroid(epot
     &          ,only_long, polar_allbeads)
             epot_ctr = epot_ctr + epot
            endif
   
            call compute_energy_beads(epot,polymer
     &            ,.not. (centroid_longrange .or. contract)  !LONG 
     &            ,.not. contract ! INT
     &            , .TRUE.     !SHORT 
     &            ,polar_allbeads)
            epot_ctr = epot_ctr + epot

            if(contract .and. nbeads_ctr>1) then
              call contract_polymer(polymer,polymer_ctr)
              call comm_for_gradient(polymer_ctr,.FALSE.,.FALSE.)

              call compute_energy_beads(epot,polymer_ctr
     &       ,.not.centroid_longrange, .TRUE., .FALSE.  ![LONG, INT  , SHORT]
     &       ,polar_allbeads)
              epot_ctr = epot_ctr + epot
            endif

            if(ranktot==0) then
              call MPI_REDUCE(MPI_IN_PLACE,epot_ctr
     &          ,1,MPI_RPREC
     &          ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            else
              call MPI_REDUCE(epot_ctr,epot_ctr
     &           ,1,MPI_RPREC
     &          ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            endif

            ub0_ctr(i) = epot_ctr
            if (verbose.and.ranktot==0) then
            write(iout,'(i11,2x,2f16.4)') i,ub0(i), ub0_ctr(i)
            endif
         else
            if (verbose.and.ranktot==0) then
            write(iout,'(i11,2x,1f16.4)') i,ub0(i)
            endif
         endif

!$acc parallel loop async default(present) collapse(2)
         DO i=1,n; do j=1,3
           polymer%eigpos(j,i,1)=0.d0
         ENDDO; ENDDO
!$acc wait
         do ibead=1,nbeads
            if (dcdio) then
              call dcdfile_read_next(polymer%dcdinfo(ibead))
              call dcdfile_skip_next(polymer%dcdinfo(ibead),0)
              if (abort) cycle mainb1
            else
              call readxyz (polymer%dcdinfo(ibead)%idcd)
            end if
            if (n.eq.0) cycle mainb1
!$acc wait
!$acc parallel loop default(present)
            DO i=1,n
              polymer%pos(1,i,ibead)=x(i)
              polymer%pos(2,i,ibead)=y(i)
              polymer%pos(3,i,ibead)=z(i)
              polymer%eigpos(1,i,1)=polymer%eigpos(1,i,1)+x(i)/nbeads
              polymer%eigpos(2,i,1)=polymer%eigpos(2,i,1)+y(i)/nbeads
              polymer%eigpos(3,i,1)=polymer%eigpos(3,i,1)+z(i)/nbeads
            ENDDO
         enddo
c
c     the box shape can change between frames
c
         use_shortmlist = .false.
         use_shortclist = .false.
         use_shortvlist = .false.
         call load_bead(polymer, 0,.true.)
         call lattice
         call AllDirAssign
         call AllRecAssign
         call reinitnl(0)
         call mechanicstep(0)
         call nblist(0)
         if (i .ge. maxframe)  abort = .true.
         if (mod(i,100).eq.0 .or. abort) then
            if (ranktot.eq.0) write (iout,160)  i
  160       format (7x,'Completed',i8,' Coordinate Frames')
            flush (iout)
         end if

      end do mainb1
c
c     reset trajectory B
c
!$acc parallel loop async default(present) collapse(2)
      DO i=1,n; do j=1,3
        polymer%eigpos(j,i,1)=0.d0
      ENDDO; ENDDO
!$acc wait
      do ibead=1,nbeads
        write(numberbeads, '(i3.3)') ibead
        exten='_beads'//numberbeads
        if (dcdio) then
          call dcdfile_close(polymer%dcdinfo(ibead))
          dcdfile = fnameb(1:lengb)//trim(exten)//'.dcd'
          call dcdfile_open(polymer%dcdinfo(ibead),dcdfile)
          call dcdfile_read_header(polymer%dcdinfo(ibead),.false.
     &       ,verbose=.false.)
          call dcdfile_read_next(polymer%dcdinfo(ibead))
          call dcdfile_skip_next(polymer%dcdinfo(ibead),0)
          abort = .false.
        else
          rewind (unit=polymer%dcdinfo(ibead)%idcd)
          call readxyz (polymer%dcdinfo(ibead)%idcd)
        end if
!$acc wait
!$acc parallel loop default(present)
        DO i=1,n
          polymer%pos(1,i,ibead)=x(i)
          polymer%pos(2,i,ibead)=y(i)
          polymer%pos(3,i,ibead)=z(i)
          polymer%eigpos(1,i,1)=polymer%eigpos(1,i,1)+x(i)/nbeads
          polymer%eigpos(2,i,1)=polymer%eigpos(2,i,1)+y(i)/nbeads
          polymer%eigpos(3,i,1)=polymer%eigpos(3,i,1)+z(i)/nbeads
        ENDDO
      enddo
c
c     Reset State 1 parameters
c
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      call load_bead(polymer, 0,.true.)
      call cutoffs
      call lattice
      call AllDirAssign
      call reinitnl(0)
      call mechanic
      call nblist(0)
c
c     find potential energies for trajectory B in state 1
c
      if (verbose) then
         if (ranktot.eq.0) write (iout,170)
  170    format (/,' Potential Energy Values for Trajectory B :',
     &           //,7x,'Frame',9x,'State 0',9x,'State 1',12x,'Delta',/)
      end if
      i = 0
      mainb2: do while (.not. abort)
         i = i + 1

         call compute_energy_beads(epot,polymer
     &     ,.TRUE.,.TRUE.,.TRUE.,.TRUE.)
         if(ranktot==0) then
          call MPI_REDUCE(MPI_IN_PLACE,epot
     &       ,1,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         else
          call MPI_REDUCE(epot,epot
     &        ,1,MPI_RPREC
     &       ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         endif

         ub1(i) = epot 

         if(centroid_recip .or. contract) then
            only_long = centroid_longrange .and.
     &         (.not. (contract.and.nbeads_ctr==1) )
            if(only_long) then
              use_shortmlist = use_mlist
              use_shortclist = use_clist
              use_shortvlist = use_vlist
              call load_bead(polymer, 0,.true.)
              call reinitnl(0)
              call mechanicstep(0)
              call nblist(0)
            endif
            epot_ctr=0.d0
            if(centroid_recip) then
             call load_bead(polymer,0)
             call compute_eslow_centroid(epot
     &          ,only_long, polar_allbeads)
             epot_ctr = epot_ctr + epot
            endif
   
            call compute_energy_beads(epot,polymer
     &            ,.not. (centroid_longrange .or. contract)  !LONG 
     &            ,.not. contract ! INT
     &            , .TRUE.     !SHORT 
     &            ,polar_allbeads)
            epot_ctr = epot_ctr + epot

            if(contract .and. nbeads_ctr>1) then
              call contract_polymer(polymer,polymer_ctr)
              call comm_for_gradient(polymer_ctr,.FALSE.,.FALSE.)

              call compute_energy_beads(epot,polymer_ctr
     &       ,.not.centroid_longrange, .TRUE., .FALSE.  ![LONG, INT  , SHORT]
     &       ,polar_allbeads)
              epot_ctr = epot_ctr + epot
            endif

            if(ranktot==0) then
              call MPI_REDUCE(MPI_IN_PLACE,epot_ctr
     &          ,1,MPI_RPREC
     &          ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            else
              call MPI_REDUCE(epot_ctr,epot_ctr
     &           ,1,MPI_RPREC
     &          ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            endif

            ub1_ctr(i) = epot_ctr
            if (verbose) then
               if (ranktot.eq.0) write (iout,'(i11,2x,6f16.4)')
     $              i,ub0(i),ub1(i),
     $              ub1(i)-ub0(i),ub0_ctr(i),ub1_ctr(i),
     $              ub1_ctr(i)-ub0_ctr(i)
            end if
         else
            if (verbose) then
               if (ranktot.eq.0) write (iout,'(i11,2x,3f16.4)')  
     &               i,ub0(i),ub1(i),
     $              ub1(i)-ub0(i)
            end if
         endif

!$acc parallel loop async default(present) collapse(2)
         DO i=1,n; do j=1,3
           polymer%eigpos(j,i,1)=0.d0
         ENDDO; ENDDO
!$acc wait
         do ibead=1,nbeads
            if (dcdio) then
              call dcdfile_read_next(polymer%dcdinfo(ibead))
              call dcdfile_skip_next(polymer%dcdinfo(ibead),0)
              if (abort) cycle mainb2
            else
              call readxyz (polymer%dcdinfo(ibead)%idcd)
            end if
!$acc wait
!$acc parallel loop default(present)
            DO i=1,n
              polymer%pos(1,i,ibead)=x(i)
              polymer%pos(2,i,ibead)=y(i)
              polymer%pos(3,i,ibead)=z(i)
              polymer%eigpos(1,i,1)=polymer%eigpos(1,i,1)+x(i)/nbeads
              polymer%eigpos(2,i,1)=polymer%eigpos(2,i,1)+y(i)/nbeads
              polymer%eigpos(3,i,1)=polymer%eigpos(3,i,1)+z(i)/nbeads
            ENDDO
         enddo
         if (n.eq.0) cycle

         !The box shape can change between frames
         use_shortmlist = .false.
         use_shortclist = .false.
         use_shortvlist = .false.
         call load_bead(polymer, 0,.true.)
         call lattice
         call AllDirAssign
         call AllRecAssign
         call reinitnl(0)
         call mechanicstep(0)
         call nblist(0)

         if (i .ge. maxframe)  abort = .true.
      end do mainb2
      nfrmb = i
      do ibead=1,nbeads
        if (dcdio) then
          call dcdfile_close(polymer%dcdinfo(ibead))
        else
          close (unit=polymer%dcdinfo(ibead)%idcd)
        endif
      enddo
c
c     perform deallocation of some local arrays
c
      deallocate (keys0)
      deallocate (keys1)
c
c     save potential energies and system volumes to a file
c
      if(ranktot==0) then
        ibar = freeunit ()
        barfile = fnamea(1:lenga)//'.bar'
        call version (barfile,'new')
        open (unit=ibar,file=barfile,status ='new')
        write (ibar,'(i8,f10.2,2x,a)')  nfrma,tempa,titlea(1:ltitlea)
        do i = 1, nfrma
           if (vola(i) .eq. 0.0d0) then
              write (ibar,'(i8,2x,2f18.4)')  i,ua0(i),ua1(i)
           else
              write (ibar,'(i8,2x,3f18.4)')  i,ua0(i),ua1(i),vola(i)
           end if
        end do
        write (ibar,'(i8,f10.2,2x,a)')  nfrmb,tempb,titleb(1:ltitleb)
        do i = 1, nfrmb
           if (volb(i) .eq. 0.0d0) then
              write (ibar,'(i8,2x,2f18.4)')  i,ub0(i),ub1(i)
           else
              write (ibar,'(i8,2x,3f18.4)')  i,ub0(i),ub1(i),volb(i)
           end if
        end do
        write (iout,'(A,A)')  ' Potential Energy Values Written To :  ',
     &             barfile(1:trimtext(barfile))
        flush (iout)
        close (unit=ibar)

        if(contract_a .or. contract_b) then
          if(.not. contract_a) then
            ua0_ctr = ua0
            ua1_ctr = ua1
          endif
          if(.not. contract_b) then
            ub0_ctr = ub0
            ub1_ctr = ub1
          endif
          ibar = freeunit ()
          barfile=trim(barfile)//'.ctr'
          open (unit=ibar,file=barfile,status ='new')
          write (ibar,'(i8,f10.2,2x,a)')  nfrma,tempa,titlea(1:ltitlea)
          do i = 1, nfrma
             if (vola(i) .eq. 0.0d0) then
                write (ibar,'(i8,2x,2f18.4)')  i,ua0_ctr(i),ua1_ctr(i)
             else
                write (ibar,'(i8,2x,3f18.4)')  i,ua0_ctr(i),ua1_ctr(i)
     &             ,vola(i)
             end if
          end do
          write (ibar,'(i8,f10.2,2x,a)')  nfrmb,tempb,titleb(1:ltitleb)
          do i = 1, nfrmb
             if (volb(i) .eq. 0.0d0) then
                write (ibar,'(i8,2x,2f18.4)')  i,ub0_ctr(i),ub1_ctr(i)
             else
                write (ibar,'(i8,2x,3f18.4)')  i,ub0_ctr(i),ub1_ctr(i)
     &             ,volb(i)
             end if
          end do
          write (iout,'(A,A)') ' Potential Energy Values Written To :  '
     &               ,barfile(1:trimtext(barfile))
          flush (iout)
          close (unit=ibar)
        endif
      endif
c
c     perform deallocation of some local arrays
c
      deallocate (ua0)
      deallocate (ua1)
      deallocate (ub0)
      deallocate (ub1)
      deallocate (vola)
      deallocate (volb)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine barcalc  --  get thermodynamics via FEP & BAR  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine barcalcpi
      use ascii
      use files
      use inform
      use iounit
      use titles
      use random_mod
      use units, only: gasconst, prescon
      use domdec
      implicit none
      integer i,j,k,ibar
      integer size,next
      integer iter,maxiter
      integer ltitlea,ltitleb
      integer nfrma,nfrmb
      integer nfrm,nbst
      integer starta,startb
      integer stopa,stopb
      integer stepa,stepb
      integer maxframe
      integer freeunit
      integer trimtext
      integer, allocatable :: bsta(:)
      integer, allocatable :: bstb(:)
      real(en_p) rt,term
      real(en_p) rta,rtb
      real(en_p) delta,eps
      real(en_p) frma,frmb
      real(en_p) tempa,tempb
      real(en_p) cold,cnew
      real(en_p) top,top2
      real(en_p) bot,bot2
      real(en_p) fterm,rfrm
      real(en_p) sum,sum2,bst
      real(en_p) vavea,vaveb
      real(en_p) vstda,vstdb
      real(en_p) vasum,vbsum
      real(en_p) vasum2,vbsum2
      real(en_p) stdev,patm
      real(en_p) ratio
      real(en_p) cfore,cback
      real(en_p) cfsum,cbsum
      real(en_p) cfsum2,cbsum2
      real(en_p) stdcf,stdcb
      real(en_p) uave0,uave1
      real(en_p) u0sum,u1sum
      real(en_p) u0sum2,u1sum2
      real(en_p) stdev0,stdev1
      real(en_p) hfore,hback
      real(en_p) sfore,sback
      real(en_p) hfsum,hbsum
      real(en_p) hfsum2,hbsum2
      real(en_p) stdhf,stdhb
      real(en_p) fore,back
      real(en_p) epv,stdpv
      real(en_p) hdir,hbar
      real(en_p) hsum,hsum2
      real(en_p) sbar,tsbar
      real(en_p) fsum,bsum
      real(en_p) fvsum,bvsum
      real(en_p) fbvsum,vsum
      real(en_p) fbsum0,fbsum1
      real(en_p) alpha0,alpha1
      real(en_p), allocatable :: ua0(:)
      real(en_p), allocatable :: ua1(:)
      real(en_p), allocatable :: ub0(:)
      real(en_p), allocatable :: ub1(:)
      real(en_p), allocatable :: vola(:)
      real(en_p), allocatable :: volb(:)
      real(en_p), allocatable :: vloga(:)
      real(en_p), allocatable :: vlogb(:)
      logical exist,query,done
      character*240 record
      character*240 string
      character*240 titlea
      character*240 titleb
      character*240 barfile
      character*10 use_contract_str
      logical :: use_contract

c
c
c     ask the user for file with potential energies and volumes
c
      call nextarg (barfile,exist)
      if (exist) then
         call basefile (barfile)
         call suffix (barfile,'bar','old')
         inquire (file=barfile,exist=exist)
      end if
      do while (.not. exist)
         if (ranktot.eq.0) write (iout,10)
   10    format (/,' Enter Potential Energy BAR File Name :  ',$)
         read (input,20)  barfile
   20    format (a240)
         call basefile (barfile)
         call suffix (barfile,'bar','old')
         inquire (file=barfile,exist=exist)
      end do

      inquire(file=trim(barfile)//".ctr",exist=exist)
      if(exist) then
        call nextarg (use_contract_str,exist)
        exist = exist .and. str_to_bool(use_contract_str)/="U"
        if(.not. exist) then
          use_contract_str=""
          do while (str_to_bool(use_contract_str)=="U") 
            if(ranktot==0) write(iout,*) 
     &       "contraction file detected, do you want to use it ?"
            read (input,'(a10)') use_contract_str
          enddo
        endif
        use_contract = str_to_logical(use_contract_str)
        if(use_contract) then
          barfile=trim(barfile)//".ctr"
          if (ranktot==0) write(*,*) 'Using RP contractions'
        endif
      endif
c
c     perform dynamic allocation of some local arrays
c
      maxframe = 1000000
      allocate (ua0(maxframe))
      allocate (ua1(maxframe))
      allocate (ub0(maxframe))
      allocate (ub1(maxframe))
      allocate (vola(maxframe))
      allocate (volb(maxframe))
      allocate (vloga(maxframe))
      allocate (vlogb(maxframe))
c
c     set beginning and ending frame for trajectory A
c
      starta = 0
      stopa  = 0
      stepa  = 0
      query  = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=30,end=30)  starta
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  stopa
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  stepa
   30 continue
      if (query) then
         if (ranktot.eq.0) write (iout,40)
   40    format (/,' First & Last Frame and Step Increment',
     &              ' for Trajectory A :  ',$)
         read (input,50)  record
   50    format (a240)
         read (record,*,err=60,end=60)  starta,stopa,stepa
   60    continue
      end if
      if (starta .eq. 0)  starta = 1
      if (stopa .eq. 0)  stopa = maxframe
      if (stepa .eq. 0)  stepa = 1
c
c     set beginning and ending frame for trajectory B
c
      startb = 0
      stopb  = 0
      stepb  = 0
      query  = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=70,end=70)  startb
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  stopb
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  stepb
   70 continue
      if (query) then
         if (ranktot.eq.0) write (iout,80)
   80    format (/,' First & Last Frame and Step Increment',
     &              ' for Trajectory B :  ',$)
         read (input,90)  record
   90    format (a240)
         read (record,*,err=100,end=100)  startb,stopb,stepb
  100    continue
      end if
      if (startb .eq. 0) startb = 1
      if (stopb  .eq. 0) stopb  = maxframe
      if (stepb  .eq. 0) stepb  = 1
c
c     read potential energies and volumes for trajectory A
c
      ibar = freeunit ()
      open (unit=ibar,file=barfile,status='old')
      rewind (unit=ibar)
      nfrma = 0
      tempa = En_p(0.0)
      next = 1
      read (ibar,110,err=250,end=250)  record
  110 format (a240)
      size = trimtext (record)
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  nfrma
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  tempa
      titlea = record(next:trimtext(record))
      call trimhead (titlea)
      ltitlea = trimtext(titlea)
      do i = 1, starta-1
         read (ibar,120,err=160,end=160)
  120    format ()
      end do
      j = 0
      stopa = (min(stopa,nfrma)-starta+1)/stepa
      do i = 1, stopa
         read (ibar,130,err=160,end=160)  record
  130    format (a240)
         j = j + 1
         ua0(j)  = En_p(0.0)
         ua1(j)  = En_p(0.0)
         vola(j) = En_p(0.0)
         read (record,*,err=140,end=140)  k,ua0(j),ua1(j),vola(j)
  140    continue
         do k = 1, stepa-1
            read (ibar,150,err=160,end=160)
  150       format ()
         end do
      end do
  160 continue
      nfrma = j
c
c     reset the file position to the beginning of trajectory B
c
      rewind (unit=ibar)
      next = 1
      read (ibar,170,err=190,end=190)  record
  170 format (a240)
      size = trimtext (record)
      call gettext (record,string,next)
      read (string,*,err=190,end=190)  k
      do i = 1, k
         read (ibar,180,err=190,end=190)
  180    format ()
      end do
  190 continue
c
c     read potential energies and volumes for trajectory B
c
      nfrmb = 0
      tempb = En_p(0.0)
      next = 1
      read (ibar,200,err=250,end=250)  record
  200 format (a240)
      size = trimtext (record)
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  nfrmb
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  tempb
      titleb = record(next:trimtext(record))
      call trimhead (titleb)
      ltitleb = trimtext(titleb)
      do i = 1, startb-1
         read (ibar,210,err=250,end=250)
  210    format ()
      end do
      j = 0
      stopb = (min(stopb,nfrmb)-startb+1)/stepb
      do i = 1, stopb
         read (ibar,220,err=250,end=250)  record
  220    format (a240)
         j = j + 1
          ub0(j) = En_p(0.0)
          ub1(j) = En_p(0.0)
         volb(j) = En_p(0.0)
         read (record,*,err=230,end=230)  k,ub0(j),ub1(j),volb(j)
  230    continue
         do k = 1, stepb-1
            read (ibar,240,err=250,end=250)
  240       format ()
         end do
      end do
  250 continue
      nfrmb = j
      close (unit=ibar)
c
c     provide info about trajectories and number of frames
c
      if (ranktot.eq.0) write (iout,260)  titlea(1:ltitlea)
  260 format (/,' Simulation Trajectory A and Thermodynamic',
     &           ' State 0 :  ',//,' ',a)
      if (ranktot.eq.0) write (iout,270)  nfrma,tempa
  270 format (' Number of Frames :',4x,i8,/,' Temperature :',7x,f10.2)
      if (ranktot.eq.0) write (iout,280)  titleb(1:ltitleb)
  280 format (/,' Simulation Trajectory B and Thermodynamic',
     &           ' State 1 :  ',//,' ',a)
      if (ranktot.eq.0) write (iout,290)  nfrmb,tempb
  290 format (' Number of Frames :',4x,i8,/,' Temperature :',7x,f10.2)
c
c     set the frame ratio, temperature and Boltzmann factor
c
      term = random ()
      frma = dble(nfrma)
      frmb = dble(nfrmb)
      rfrm = frma / frmb
      rta = gasconst * tempa
      rtb = gasconst * tempb
      rt = 0.5d0 * (rta+rtb)
c
c     set the number of bootstrap trials to be generated
c
      nfrm = max(nfrma,nfrmb)
      nbst = min(100000,nint(1.0d8/dble(nfrm)))
      bst = dble(nbst)
      ratio = bst / (bst-1.0d0)
c
c     find average volumes and corrections for both trajectories
c
      vasum  = En_p(0.0)
      vasum2 = En_p(0.0)
      vbsum  = En_p(0.0)
      vbsum2 = En_p(0.0)
      do k = 1, nbst
         sum = En_p(0.0)
         do i = 1, nfrma
            j = max( int(frma*random()),1 )
            sum = sum + vola(j)
         end do
         vavea = sum / frma
         vasum = vasum + vavea
         vasum2 = vasum2 + vavea*vavea
         sum = En_p(0.0)
         do i = 1, nfrmb
            j = max( int(frmb*random()),1 )
            sum = sum + volb(j)
         end do
         vaveb = sum / frmb
         vbsum = vbsum + vaveb
         vbsum2 = vbsum2 + vaveb*vaveb
      end do
      vavea = vasum / bst
      vstda = sqrt(ratio*(vasum2/bst-vavea*vavea))
      vaveb = vbsum / bst
      vstdb = sqrt(ratio*(vbsum2/bst-vaveb*vaveb))
      if (vavea .ne. En_p(0.0)) then
         do i = 1, nfrma
            if (vola(i) .ne. En_p(0.0))
     &         vloga(i) = -rta * log(vola(i)/vavea)
         end do
      end if
      if (vaveb .ne. En_p(0.0)) then
         do i = 1, nfrmb
            if (volb(i) .ne. En_p(0.0))
     &         vlogb(i) = -rtb * log(volb(i)/vaveb)
         end do
      end if
c
c     get the free energy change via thermodynamic perturbation
c
      if (ranktot.eq.0) write (iout,300)
  300 format (/,' Free Energy Difference via FEP Method :',/)
      cfsum  = En_p(0.0)
      cfsum2 = En_p(0.0)
      cbsum  = En_p(0.0)
      cbsum2 = En_p(0.0)
      do k = 1, nbst
         sum = En_p(0.0)
         do i = 1, nfrma
            j   = max( int(frma*random()),1 )
            sum = sum + exp((ua0(j)-ua1(j)+vloga(j))/rta)
         end do
         cfore  = -rta * log(sum/frma)
         cfsum  = cfsum + cfore
         cfsum2 = cfsum2 + cfore*cfore
         sum = En_p(0.0)
         do i = 1, nfrmb
            j   = max( int(frmb*random()),1 )
            sum = sum + exp((ub1(j)-ub0(j)+vlogb(j))/rtb)
         end do
         cback  = rtb * log(sum/frmb)
         cbsum  = cbsum + cback
         cbsum2 = cbsum2 + cback*cback
      end do
      cfore = cfsum / bst
      stdcf = sqrt(ratio*(cfsum2/bst-cfore*cfore))
      cback = cbsum / bst
      stdcb = sqrt(ratio*(cbsum2/bst-cback*cback))
      if (ranktot.eq.0) write (iout,310)  cfore,stdcf
  310 format (' Free Energy via Forward FEP',9x,f12.4,
     &           ' +/-',f9.4,' Kcal/mol')
      if (ranktot.eq.0) write (iout,320)  cback,stdcb
  320 format (' Free Energy via Backward FEP',8x,f12.4,
     &           ' +/-',f9.4,' Kcal/mol')
c
c     determine the initial free energy via the BAR method
c
      if (ranktot.eq.0) write (iout,330)
  330 format (/,' Free Energy Difference via BAR Method :',/)
      maxiter = 100
      eps     = En_p(0.0001)
      done    = .false.
      iter    = 0
      cold    = En_p(0.0)
      top     = En_p(0.0)
      top2    = En_p(0.0)
      do i = 1, nfrmb
         fterm = 1.0d0 / (1.0d0+exp((ub0(i)-ub1(i)+vlogb(i)+cold)/rtb))
         top   = top + fterm
         top2  = top2 + fterm*fterm
      end do
      bot  = En_p(0.0)
      bot2 = En_p(0.0)
      do i = 1, nfrma
         fterm = 1.0d0 / (1.0d0+exp((ua1(i)-ua0(i)+vloga(i)-cold)/rta))
         bot = bot + fterm
         bot2 = bot2 + fterm*fterm
      end do
      cnew = rt*log(rfrm*top/bot) + cold
      stdev = sqrt((bot2-bot*bot/frma)/(bot*bot)
     &                + (top2-top*top/frmb)/(top*top))
      delta = abs(cnew-cold)
      if (ranktot.eq.0) write (iout,340)  iter,cnew
  340 format (' BAR Iteration',i4,19x,f12.4,' Kcal/mol')
      if (delta .lt. eps) then
         done = .true.
         if (ranktot.eq.0) write (iout,350)  cnew,stdev
  350    format (' BAR Free Energy Estimate',8x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
      end if
c
c     iterate the BAR equation to converge the free energy
c
      do while (.not. done)
         iter = iter + 1
         cold = cnew
         top  = En_p(0.0)
         top2 = En_p(0.0)
         do i = 1, nfrmb
            fterm = 1.0d0 / (1.0d0+exp((ub0(i)-ub1(i)+vlogb(i)
     &                                     +cold)/rtb))
            top = top + fterm
            top2 = top2 + fterm*fterm
         end do
         bot  = En_p(0.0)
         bot2 = En_p(0.0)
         do i = 1, nfrma
            fterm = 1.0d0 / (1.0d0+exp((ua1(i)-ua0(i)+vloga(i)
     &                                        -cold)/rta))
            bot  = bot + fterm
            bot2 = bot2 + fterm*fterm
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         stdev = sqrt((bot2-bot*bot/frma)/(bot*bot)
     &                   + (top2-top*top/frmb)/(top*top))
         delta = abs(cnew-cold)
         if (ranktot.eq.0) write (iout,360)  iter,cnew
  360    format (' BAR Iteration',i4,19x,f12.4,' Kcal/mol')
         if (delta .lt. eps) then
            done = .true.
            if (ranktot.eq.0) write (iout,370)  cnew,stdev
  370       format (/,' Free Energy via BAR Iteration',7x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         end if
         if (iter.ge.maxiter .and. .not.done) then
            done = .true.
            if (ranktot.eq.0) write (iout,380)  maxiter
  380       format (/,' BAR Free Energy Estimate not Converged',
     &                 ' after',i4,' Iterations')
            call fatal
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (bsta(nfrm))
      allocate (bstb(nfrm))
c
c     use bootstrap analysis to estimate statistical error
c
      sum  = En_p(0.0)
      sum2 = En_p(0.0)
      do k = 1, nbst
         done = .false.
         iter = 0
         cold = En_p(0.0)
         top  = En_p(0.0)
         do i = 1, nfrmb
            j = max( int(frmb*random()),1 )
            bstb(i) = j
            top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)
     &                                   +vlogb(i)+cold)/rtb))
         end do
         bot  = 0.0d0
         do i = 1, nfrma
            j = max( int(frma*random()),1 )
            bsta(i) = j
            bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)
     &                                   +vloga(i)-cold)/rta))
         end do
         cnew  = rt*log(rfrm*top/bot) + cold
         delta = abs(cnew-cold)
         do while (.not. done)
            iter = iter + 1
            cold = cnew
            top  = En_p(0.0)
            do i = 1, nfrmb
               j   = bstb(i)
               top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)
     &                                      +vlogb(i)+cold)/rtb))
            end do
            bot = En_p(0.0)
            do i = 1, nfrma
               j   = bsta(i)
               bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)
     &                                      +vloga(i)-cold)/rta))
            end do
            cnew  = rt*log(rfrm*top/bot) + cold
            delta = abs(cnew-cold)
            if (delta .lt. eps) then
               done = .true.
               sum  = sum + cnew
               sum2 = sum2 + cnew*cnew
            end if
         end do
      end do
      cnew  = sum / bst
      ratio = bst / (bst-En_p(1.0))
      stdev = sqrt(ratio*(sum2/bst-cnew*cnew))
      if (ranktot.eq.0) write (iout,390)  cnew,stdev
  390 format (' Free Energy via BAR Bootstrap',7x,f12.4,
     &           ' +/-',f9.4,' Kcal/mol')
c
c     perform deallocation of some local arrays
c
      deallocate (bsta)
      deallocate (bstb)
c
c     find the enthalpy directly via average potential energy
c
      if (ranktot.eq.0) write (iout,400)
  400 format (/,' Enthalpy from Potential Energy Averages :',/)
      patm   = En_p(1.0)
      epv    = (vaveb-vavea) * patm / prescon
      stdpv  = (vstda+vstdb) * patm / prescon
      u0sum  = En_p(0.0)
      u0sum2 = En_p(0.0)
      u1sum  = En_p(0.0)
      u1sum2 = En_p(0.0)
      hsum   = En_p(0.0)
      hsum2  = En_p(0.0)
      do k = 1, nbst
         uave0 = En_p(0.0)
         uave1 = En_p(0.0)
         do i = 1, nfrma
            j = max( int(frma*random()),1 )
            uave0 = uave0 + ua0(j)
         end do
         do i = 1, nfrmb
            j = max( int(frmb*random()),1 )
            uave1 = uave1 + ub1(j)
         end do
         uave0  = uave0 / frma
         uave1  = uave1 / frmb
         u0sum  = u0sum + uave0
         u0sum2 = u0sum2 + uave0*uave0
         u1sum  = u1sum + uave1
         u1sum2 = u1sum2 + uave1*uave1
         hdir   = uave1 - uave0 + epv
         hsum   = hsum + hdir
         hsum2  = hsum2 + hdir*hdir
      end do
      uave0  = u0sum / bst
      stdev0 = sqrt(ratio*(u0sum2/bst-uave0*uave0))
      uave1  = u1sum / bst
      stdev1 = sqrt(ratio*(u1sum2/bst-uave1*uave1))
      hdir   = hsum / bst
      stdev  = sqrt(ratio*(hsum2/bst-hdir*hdir))
      if (ranktot.eq.0) write (iout,410)  uave0,stdev0
  410 format (' Average Energy for State 0',10x,f12.4,
     &           ' +/-',f9.4,' Kcal/mol')
      if (ranktot.eq.0) write (iout,420)  uave1,stdev1
  420 format (' Average Energy for State 1',10x,f12.4,
     &           ' +/-',f9.4,' Kcal/mol')
      if (epv .ne. 0.0d0) then
         if (ranktot.eq.0) write (iout,430)  epv,stdpv
  430    format (' PdV Work Term for 1 Atm',13x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
      end if
      if (ranktot.eq.0) write (iout,440)  hdir,stdev
  440 format (' Enthalpy via Direct Estimate',8x,f12.4,
     &           ' +/-',f9.4,' Kcal/mol')
c
c     calculate the enthalpy via thermodynamic perturbation
c
      if (ranktot.eq.0) write (iout,450)
  450 format (/,' Enthalpy and Entropy via FEP Method :',/)
      hfsum  = En_p(0.0)
      hfsum2 = En_p(0.0)
      hbsum  = En_p(0.0)
      hbsum2 = En_p(0.0)
      do k = 1, nbst
         top = En_p(0.0)
         bot = En_p(0.0)
         do i = 1, nfrma
            j    = max( int(frma*random()),1 )
            term = exp((ua0(j)-ua1(j)+vloga(j))/rta)
            top  = top + ua1(j)*term
            bot  = bot + term
         end do
         hfore  = (top/bot) - uave0
         hfsum  = hfsum + hfore
         hfsum2 = hfsum2 + hfore*hfore
         top    = En_p(0.0)
         bot    = En_p(0.0)
         do i = 1, nfrmb
            j    = max( int(frmb*random()),1 )
            term = exp((ub1(j)-ub0(j)+vlogb(j))/rtb)
            top  = top + ub0(j)*term
            bot  = bot + term
         end do
         hback  = -(top/bot) + uave1
         hbsum  = hbsum + hback
         hbsum2 = hbsum2 + hback*hback
      end do
      hfore = hfsum / bst
      stdhf = sqrt(ratio*(hfsum2/bst-hfore*hfore))
      stdhf = stdhf + stdev0
      sfore = (hfore-cfore) / tempa
      hback = hbsum / bst
      stdhb = sqrt(ratio*(hbsum2/bst-hback*hback))
      stdhb = stdhb + stdev1
      sback = (hback-cback) / tempb
      if (ranktot.eq.0) write (iout,460)  hfore,stdhf
  460 format (' Enthalpy via Forward FEP',12x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
      if (ranktot.eq.0) write (iout,470)  sfore
  470 format (' Entropy via Forward FEP',13x,f12.6,' Kcal/mol/K')
      if (ranktot.eq.0) write (iout,480)  -tempa*sfore
  480 format (' Forward FEP -T*dS Value',13x,f12.4,' Kcal/mol')
      if (ranktot.eq.0) write (iout,490)  hback,stdhb
  490 format (/,' Enthalpy via Backward FEP',11x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
      if (ranktot.eq.0) write (iout,500)  sback
  500 format (' Entropy via Backward FEP',12x,f12.6,' Kcal/mol/K')
      if (ranktot.eq.0) write (iout,510)  -tempb*sback
  510 format (' Backward FEP -T*dS Value',12x,f12.4,' Kcal/mol')
c
c     determine the enthalpy and entropy via the BAR method
c
      if (ranktot.eq.0) write (iout,520)
  520 format (/,' Enthalpy and Entropy via BAR Method :',/)
      hsum  = En_p(0.0)
      hsum2 = En_p(0.0)
      do k = 1, nbst
         fsum   = En_p(0.0)
         fvsum  = En_p(0.0)
         fbvsum = En_p(0.0)
         vsum   = En_p(0.0)
         fbsum0 = En_p(0.0)
         do i = 1, nfrma
            j     = max( int(frma*random()),1 )
            fore  = 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)+vloga(j)-cnew)/rta))
            back  = 1.0d0/(1.0d0+exp((ua0(j)-ua1(j)+vloga(j)+cnew)/rta))
            fsum  = fsum + fore
            fvsum = fvsum + fore*ua0(j)
           fbvsum = fbvsum + fore*back*(ua1(j)-ua0(j)+vloga(j))
           vsum   = vsum + ua0(j)
           fbsum0 = fbsum0 + fore*back
         end do
         alpha0 = fvsum - fsum*(vsum/frma) + fbvsum
         bsum   = En_p(0.0)
         bvsum  = En_p(0.0)
         fbvsum = En_p(0.0)
         vsum   = En_p(0.0)
         fbsum1 = En_p(0.0)
         do i = 1, nfrmb
            j = max( int(frmb*random()),1 )
            fore  = 1.0d0/(1.0d0+exp((ub1(j)-ub0(j)+vlogb(j)-cnew)/rtb))
            back  = 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+vlogb(j)+cnew)/rtb))
            bsum  = bsum + back
            bvsum = bvsum + back*ub1(j)
           fbvsum = fbvsum + fore*back*(ub1(j)-ub0(j)+vlogb(j))
           vsum   = vsum + ub1(j)
           fbsum1 = fbsum1 + fore*back
         end do
         alpha1 = bvsum - bsum*(vsum/frmb) - fbvsum
         hbar   = (alpha0-alpha1) / (fbsum0+fbsum1)
         hsum   = hsum + hbar
         hsum2  = hsum2 + hbar*hbar
      end do
      hbar  = hsum / bst
      stdev = sqrt(ratio*(hsum2/bst-hbar*hbar))
      tsbar = hbar - cnew
      sbar  = tsbar / (En_p(0.5)*(tempa+tempb))
      if (ranktot.eq.0) write (iout,530)  hbar,stdev
  530 format (' Enthalpy via BAR Estimate',11x,f12.4
     &           ' +/-',f9.4,' Kcal/mol')
      if (ranktot.eq.0) write (iout,540)  sbar
  540 format (' Entropy via BAR Estimate',12x,f12.6,' Kcal/mol/K')
      if (ranktot.eq.0) write (iout,550)  -tsbar
  550 format (' BAR Estimate of -T*dS',15x,f12.4,' Kcal/mol')
c
c     perform deallocation of some local arrays
c
      deallocate (ua0)
      deallocate (ua1)
      deallocate (ub0)
      deallocate (ub1)
      deallocate (vola)
      deallocate (volb)
      deallocate (vloga)
      deallocate (vlogb)
      end
