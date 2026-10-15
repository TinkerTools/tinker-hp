!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  program analyze_beads  --  energy partitioning and analysis  ##
!     ##                               for PIMD                        ##
!     ###################################################################
!
!
!     "analyze_beads" computes and displays the total potential energy;
!
!
#include "tinker_macro.h"
program analyze_beads
   use mpi
   use config, only: raw_cli_t, capture_raw_cli
   use iounit, only: iout
   use mdevice, only: bind_device
   use status, only: failed, status_t
   implicit none
   integer ierr,world_rank
   type(raw_cli_t) :: raw_cli
   type(status_t) :: status_value
   call bind_device(status_value)
   if (failed(status_value)) stop 1
   call capture_raw_cli(raw_cli)
   call MPI_INIT(ierr)
!      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,nthreadsupport,ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,world_rank,ierr)
   iout = 6
   call analyze_beads_bis(raw_cli,status_value)
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
subroutine analyze_beads_bis(raw_cli,status_value)
   use ascii, only: int_to_str
   use atomsMirror
   use bath
   use commstuffpi
   use dcdmod
   use domdec
   use energi
   use files
   use inform
   use config, only: config_t, raw_cli_t
   use initialization, only: analyze_beads_application, initialize
   use mdevice, only: channel_t, finalize_device
   use parallel, only: parallel_t, release, assign_direct, assign_reciprocal, &
      &exchange_polymer_gradient_state
   use status, only: failed, status_t
   use structure, only: read_frame, structure_t
   use iounit
   use mpi
   use beads
   use utilbaoabpi
   use cutoff
   use keys
   use units
   use tinMemory, only: prmem_requestm
   use utils
   use virial
   use moldyn
   implicit none
   integer i,ierr,j,k,ibead,ibeadglob,ixyz
   integer frame,nbeads_para
   integer opened_beads
   integer mode,next,nseg
   integer freeunit
   integer trimtext
   real(r_p) :: epot,epot_full,epot_ctr
   real(r_p) energy
   logical doenergy,dodipoltot,dodipolmol
   logical exist
   character*1 letter
   character*3 numberbeads
   character*20 keyword
   character*240 record
   character*240 string
   character*240 xyzfile
   character*240 dcdfile,exten
   logical :: only_long
   type(dcdinfo_t) :: dcdinfo
   type(raw_cli_t), intent(in) :: raw_cli
   type(config_t) :: configuration
   type(structure_t) :: molecular_structure
   type(parallel_t) :: platform
   type(channel_t) :: reciprocal_work
   type(status_t) :: cleanup_status
   type(status_t), intent(out) :: status_value

!
   call initialize(analyze_beads_application, raw_cli, configuration, molecular_structure, &
      &platform, reciprocal_work, status_value)
   if (failed(status_value)) return
   opened_beads = 0
!
!     call nblist(0)

   call reset_observables(.true.)
   call allocpi(status_value)
   if (failed(status_value)) return

!
!     get the desired types of analysis to be performed
!
   call nextarg (string,exist)
   if (.not. exist) then
      write (iout,10)
10    format (/,' The TINKER Analysis Facility can Provide :',&
         &/,' Total Potential Energy and its Components [E]',&
         &/,' Total Dipolar Moment [D]',&
         &/,' Molecular Dipolar Moments [M]')
20    continue
      write (iout,30)
30    format (/,' Enter the Desired Analysis Types',&
         &' [E,D,M] :  ',$)
      read (input,40,err=20)  string
40    format (a240)
   end if

!
!     set option control flags based desired analysis types
!
   doenergy   = .false.
   dodipoltot = .false.
   dodipolmol = .false.
   call upcase (string)
   do i = 1, trimtext(string)
      letter = string(i:i)
      if (letter .eq. 'E')  doenergy = .true.
      if (letter .eq. 'D')  dodipoltot = .true.
      if (letter .eq. 'M')  dodipolmol = .true.
   end do

   if(dodipolmol .or. dodipoltot) then
      if(ranktot==0) then
         write(iout,*) 'Only energy analysis is implementd for beads'
      endif
      dodipolmol=.false.
      dodipoltot=.false.
   end if

   kelvin = -1.0d0
   call nextarg (string,exist)
   if (exist)  read (string,*,err=50,end=50)  kelvin
50 continue
   do while (kelvin .lt. 0.0d0)
      if (ranktot.eq.0) write (iout,'(a)', advance='no')&
         &' Enter Trajectory Temperature in Degrees K [298] : '
      read (input,'(f20.0)',err=60) kelvin
      if (kelvin .le. 0.d0)  kelvin = 298.0d0
60    continue
   end do
!
!     reopen the coordinates file and read the first structure
!
   frame = 0
   if(.not. allocated(polymer_main%dcdinfo)) then
      allocate(polymer_main%dcdinfo(polymer_main%nbeads))
   endif
!$acc parallel loop async default(present) collapse(2)
   DO i=1,n; do j=1,3
         polymer_main%eigpos(j,i,1)=0.d0
      ENDDO; ENDDO
!$acc wait
   do ibead=1,nbeads
      write(numberbeads, '(i3.3)') ibead
      exten='_beads'//numberbeads
      if (dcdio) then
         dcdfile = filename(1:leng)//trim(exten)//'.dcd'
         call dcdfile_open(polymer_main%dcdinfo(ibead),dcdfile)
         call dcdfile_read_header(polymer_main%dcdinfo(ibead),.false.&
            &,verbose=.false.)
         call dcdfile_read_next(polymer_main%dcdinfo(ibead))
         call dcdfile_skip_next(polymer_main%dcdinfo(ibead),0)
      else
         ixyz = freeunit ()
         polymer_main%dcdinfo(ibead)%idcd = ixyz
         xyzfile = filename(1:leng)//trim(exten)
         call suffix (xyzfile,'arc','old')
         open (unit=ixyz,file=xyzfile,status ='old')
         rewind (unit=ixyz)
         opened_beads = ibead
         call read_frame(ixyz,molecular_structure,status_value)
         if (failed(status_value)) goto 900
      end if
      if (dcdio) opened_beads = ibead
!$acc wait
!$acc parallel loop default(present)
      DO i=1,n
         polymer_main%pos(1,i,ibead)=x(i)
         polymer_main%pos(2,i,ibead)=y(i)
         polymer_main%pos(3,i,ibead)=z(i)
         polymer_main%eigpos(1,i,1)=polymer_main%eigpos(1,i,1)+x(i)/nbeads
         polymer_main%eigpos(2,i,1)=polymer_main%eigpos(2,i,1)+y(i)/nbeads
         polymer_main%eigpos(3,i,1)=polymer_main%eigpos(3,i,1)+z(i)/nbeads
      ENDDO
   enddo
!
!     perform analysis for each successive coordinate structure
!
   calc_e = .true.
   use_virial = .false.
   do while (.not. abort)
      frame = frame + 1
      if (frame .gt. 1) then
         if (ranktot.eq.0) write (iout,90)  frame
90       format (/,' Analysis for Archive Structure :',8x,i8)
      end if
!
!       setup for MPI
!
      use_shortmlist = .false.
      use_shortclist = .false.
      use_shortvlist = .false.
      call load_bead(polymer_main, 0,.true.)
      call lattice
      call assign_direct
      call assign_reciprocal
      call reinitnl(0)
      call mechanicstep(0)
      call nblist(0)

      epot = energy()
      if(ranktot==0) then
         call MPI_REDUCE(MPI_IN_PLACE,epot&
            &,1,MPI_RPREC&
            &,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      else
         call MPI_REDUCE(epot,epot&
            &,1,MPI_RPREC&
            &,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      endif
      if(ranktot==0) then
         write(iout,'(A,f20.8,A)') 'E_centroid = '&
            &,epot,' Kcal/mole'
      endif


      call compute_energy_beads(epot,polymer_main&
         &,.TRUE.,.TRUE.,.TRUE.,.TRUE.)
      if(ranktot==0) then
         call MPI_REDUCE(MPI_IN_PLACE,epot&
            &,1,MPI_RPREC&
            &,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      else
         call MPI_REDUCE(epot,epot&
            &,1,MPI_RPREC&
            &,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      endif

      if(ranktot==0) then
         epot_full = epot
         write(iout,'(A,f20.8,A)') 'E_full = '&
            &,epot_full,' Kcal/mole'
      endif


      if(centroid_recip .or. contract) then
         only_long = centroid_longrange .and.&
            &(.not. (contract.and.nbeads_ctr==1) )
         if(only_long) then
            use_shortmlist = use_mlist
            use_shortclist = use_clist
            use_shortvlist = use_vlist
            call load_bead(polymer_main, 0,.true.)
            call reinitnl(0)
            call mechanicstep(0)
            call allocstep
            call nblist(0)
            call prmem_requestm(derivs,3,nbloc,nbeadsloc,async=.true.)
         endif
         epot_ctr=0.d0
         if(centroid_recip) then
            call load_bead(polymer_main,0)
            call compute_eslow_centroid(epot&
               &,only_long, polar_allbeads)
            epot_ctr = epot_ctr + epot
         endif

         call compute_energy_beads(epot,polymer_main&
            &,.not. (centroid_longrange .or. contract)&  !LONG
            &,.not. contract& ! INT
            &, .TRUE.&     !SHORT
            &,polar_allbeads)
         epot_ctr = epot_ctr + epot

         if(contract .and. nbeads_ctr>1) then
            call contract_polymer(polymer_main,polymer_ctr)
            call exchange_polymer_gradient_state(polymer_ctr,.FALSE.,.FALSE.)

            call compute_energy_beads(epot,polymer_ctr&
               &,.not.centroid_longrange, .TRUE., .FALSE.&  ![LONG, INT  , SHORT]
               &,polar_allbeads)
            epot_ctr = epot_ctr + epot
         endif

         if(ranktot==0) then
            call MPI_REDUCE(MPI_IN_PLACE,epot_ctr&
               &,1,MPI_RPREC&
               &,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         else
            call MPI_REDUCE(epot_ctr,epot_ctr&
               &,1,MPI_RPREC&
               &,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         endif

         if(ranktot==0) then
            write(iout,'(A,f20.8,A)') 'E_ctr = '&
               &,epot_ctr,' Kcal/mole'
            write(iout,'(A,e20.8)') 'reweight = '&
               &,exp(-(epot_full-epot_ctr)/(gasconst*kelvin))
            write(iout,'(A,2f20.8)') 'dE/nat = '&
               &,(epot_full-epot_ctr)/n
         endif

      endif

!
!     attempt to read next structure from the coordinate file
!

!$acc parallel loop async default(present) collapse(2)
      DO i=1,n; do j=1,3
            polymer_main%eigpos(j,i,1)=0.d0
         ENDDO; ENDDO
!$acc wait
      do ibead=1,nbeads
         if (.not.dcdio) then
            call read_frame(polymer_main%dcdinfo(ibead)%idcd,molecular_structure,status_value)
            if (failed(status_value)) goto 900
         else
            call dcdfile_read_next(polymer_main%dcdinfo(ibead))
            call dcdfile_skip_next(polymer_main%dcdinfo(ibead),0)
         end if
!$acc wait
!$acc parallel loop default(present)
         DO i=1,n
            polymer_main%pos(1,i,ibead)=x(i)
            polymer_main%pos(2,i,ibead)=y(i)
            polymer_main%pos(3,i,ibead)=z(i)
            polymer_main%eigpos(1,i,1)=polymer_main%eigpos(1,i,1)+x(i)/nbeads
            polymer_main%eigpos(2,i,1)=polymer_main%eigpos(2,i,1)+y(i)/nbeads
            polymer_main%eigpos(3,i,1)=polymer_main%eigpos(3,i,1)+z(i)/nbeads
         ENDDO
      enddo
   end do
!
!     perform any final tasks before program exit
!
900 continue
   do ibead=1,opened_beads
      if (dcdio) then
         call dcdfile_close(polymer_main%dcdinfo(ibead))
      else
         close (unit=polymer_main%dcdinfo(ibead)%idcd)
      end if
   enddo
   call final(status_value)
   call finalize_device(cleanup_status)
   call release(platform, cleanup_status)
   if (failed(cleanup_status)) status_value = cleanup_status
end
