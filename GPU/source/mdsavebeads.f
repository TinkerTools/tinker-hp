c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine mdsavebeads  --  save trajectory and restart files   ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "mdsavebeads" writes molecular dynamics trajectory snapshots and
c     auxiliary files with velocity, force or induced dipole data;
c     also checks for user requested termination of a simulation while
c     using path pimd
c
c
#include "tinker_precision.h"
      subroutine mdsavebeads (istep,dt,polymer)
      use atmtyp
      use atomsMirror
      use beads
      use bound
      use boxes
      use domdec
      use files
      use group
      use inform
      use iounit
      use mdstuf
      use moldyn
      use mpole
      use output
      use polar
      use potent
      use titles
      use units
      use mpi
      use dcdmod
      use random_mod
      implicit none
      integer, intent(in) :: istep
      real(r_p), intent(in) :: dt
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer i,j,k,iloc
      integer ixyz,iind
      integer ivel,ifrc
      integer iend1,idump,lext
      integer freeunit,trimtext
      integer moddump,ierr
      real(r_p) pico,wt,r
      integer ibead,ibead_save
      logical exist,save_all
      character*7 ext
      character*240 endfile
      character*240 xyzfile
      character*240 velfile
      character*240 frcfile
      character*240 indfile
      character*240 exten
      character*3 numberbeads
      integer :: nproc_write,ibead_beg,ibead_end
c
      moddump = mod(istep,iwrite)
      if (moddump .ne. 0)  return

      call update_direct_space_pi(polymer)

!$acc update host(polymer%pos,polymer%vel,polymer%forces) async
        if(allocated(polymer%forces_slow)) then
!$acc update host(polymer%forces_slow) async
        endif
!$acc wait

      SELECT CASE(trim(save_beads))
      CASE('ALL')
        save_all=.true.
      CASE('RANDOM')
        save_all=.false.
        if(ranktot==0) then
          r=random()
          ibead_save=floor(r*nbeads)+1
        endif
        call MPI_BCAST(ibead_save,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
        write(*,*) "saving bead ",ibead_save
      CASE('ONE')
        save_all=.false.
        ibead_save=1
      CASE DEFAULT
        save_all=.true.
      END SELECT

      nproc_write=min(nproctot,polymer%nbeads)
      if(nproctot==1) then
        ibead_beg=1
        ibead_end=polymer%nbeads
      else
        call comm_for_write_pi(polymer,nproc_write
     &   ,ibead_beg,ibead_end)      
      endif
      if(ranktot+1>nproc_write) then
        call check_endfile()
        return
      endif   
     
      
c     get the sequence number of the current trajectory frame
c
      idump = nprior + istep/iwrite
      lext = 3
      call numeral (idump,ext,lext)
c
c     print header for the instantaneous values at current step
c
      if (ranktot==0) then
        pico = dble(istep) * dt
        write (iout,10)  istep
   10   format (/,' Instantaneous Values for Frame saved at',
     &             i10,' Dynamics Steps')
        write (iout,20)  pico
   20   format (/,' Current Time',8x,f15.4,' Picosecond')
        write (iout,30)  epotpi
   30   format (' Current Potential',3x,f15.4,' Kcal/mole')
        if (use_bounds) then
           write (iout,40)  xbox,ybox,zbox
   40      format (' Lattice Lengths',6x,3f14.6)
           write (iout,50)  alpha,beta,gamma
   50      format (' Lattice Angles',7x,3f14.6)
        end if
        write (iout,60)  idump
   60   format (' Frame Number',13x,i10)
      endif
c
c     update the information needed to restart the trajectory
c
      call prtdynbeads(polymer,ibead_beg,ibead_end)

      if(dcdio) then
        if(.not. allocated(polymer%dcdinfo)) then
          allocate(polymer%dcdinfo(polymer%nbeads))
        endif
      endif

      if(.NOT. save_all) then
        if(ibead_save < ibead_beg 
     &     .or. ibead_save > ibead_end) then
          call check_endfile()
          return
        endif
        ibead_beg=ibead_save
        ibead_end=ibead_save
        exten=''
      endif

      do ibead = ibead_beg,ibead_end 
!$acc parallel loop default(present) async
        do i=1,n
          x(i)=polymer%pos(1,i,ibead)
          y(i)=polymer%pos(2,i,ibead)
          z(i)=polymer%pos(3,i,ibead)
        enddo
c
c     save coordinates to an archive or numbered structure file
c
        if(save_all) then
          write(numberbeads, '(i3.3)') ibead
          exten='_beads'//numberbeads
        endif
        !write(*,*) 'globead=' ,globbead(ibeadsloc), 'ranktot=',ranktot
        if (dcdio) then
!$acc update host(x,y,z) async
!$acc wait
          if(save_all) then
            call dcdio_write(polymer%dcdinfo(ibead)
     &         ,istep,dt,trim(exten))
          else
            call dcdio_write(polymer%dcdinfo(1)
     &         ,istep,dt,trim(exten))
          endif
        else
          ixyz = freeunit ()
          if (archive) then
           xyzfile = filename(1:leng)//trim(exten)
           call suffix (xyzfile,'arc','old')
           inquire (file=xyzfile,exist=exist)
           if (exist) then
              call openend (ixyz,xyzfile)
           else
              open (unit=ixyz,file=xyzfile,status='new')
           end if
          else
           xyzfile = filename(1:leng)//trim(exten)
     &                          //'.'//ext(1:lext)
           call version (xyzfile,'new')
           open (unit=ixyz,file=xyzfile,status='new')
          end if
          call prtxyz (ixyz)
          close (unit=ixyz)
c         write (iout,70)  xyzfile(1:trimtext(xyzfile))
c   70   format (' Coordinate File',12x,a)
        endif
c     
c       save the velocity vector components at the current step
c
        if (velsave) then
           ivel = freeunit ()
           if (archive) then
              velfile = filename(1:leng)//trim(exten)
              call suffix (velfile,'vel','old')
              inquire (file=velfile,exist=exist)
              if (exist) then
                 call openend (ivel,velfile)
              else
                 open (unit=ivel,file=velfile,status='new')
              end if
           else
              velfile = filename(1:leng)//trim(exten)//
     &                 '.'//ext(1:lext)//'v'
              call version (velfile,'new')
              open (unit=ivel,file=velfile,status='new')
           end if
           write (ivel,110)  n,title(1:ltitle)
  110      format (i6,2x,a)
           do i = 1, n
              write (ivel,120)  i,name(i),(polymer%vel(j,i,ibead),j=1,3)
  120         format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
           end do
           close (unit=ivel)
           write (iout,130)  velfile(1:trimtext(velfile))
  130      format (' Velocity File',15x,a)
        end if
c
c       save the force vector components for the current step
c
        if (frcsave) then
           ifrc = freeunit ()
           if (archive) then
              frcfile = filename(1:leng)//trim(exten)
              call suffix (frcfile,'frc','old')
              inquire (file=frcfile,exist=exist)
              if (exist) then
                 call openend (ifrc,frcfile)
              else
                 open (unit=ifrc,file=frcfile,status='new')
              end if
           else
              frcfile = filename(1:leng)//trim(exten)//
     &                 '.'//ext(1:lext)//'f'
              call version (frcfile,'new')
              open (unit=ifrc,file=velfile,status='new')
           end if
           write (ifrc,140)  n,title(1:ltitle)
  140      format (i6,2x,a)
           do i = 1, n
              write (ifrc,150)  i,name(i),
     &            (polymer%forces(j,i,ibead),j=1,3)
  150         format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
           end do
           close (unit=ifrc)
           write (iout,160)  frcfile(1:trimtext(frcfile))
  160      format (' Force Vector File',11x,a)
        end if

      enddo
      
      call check_endfile()
      end

      subroutine check_endfile()
      use domdec
      use files
      use iounit
      use mpi
      implicit none
      integer :: ierr,iend1
      logical :: exist
      character*240 :: endfile
      integer freeunit
c
c     test for requested termination of the dynamics calculation
c
      endfile = 'tinker.end'
      inquire (file=endfile,exist=exist)
      if (.not. exist) then
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            iend1 = freeunit ()
            open (unit=iend1,file=endfile,status='old')
            close (unit=iend1,status='delete')
         end if
      end if
      if (exist) then
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(ranktot==0) then 
         write (iout,200)
  200    format (/,' MDSAVE  --  Dynamics Calculation Ending',
     &              ' due to User Request')
         call fatal
        endif
      end if

      end subroutine check_endfile

      subroutine comm_for_write_pi(polymer,nproc_write
     &   ,ibead_beg_r,ibead_end_r)
        use atoms
        use beads
        use mpi
        use domdec
        implicit none
        type(POLYMER_COMM_TYPE), intent(inout) :: polymer
        integer, intent(in) :: nproc_write
        integer, intent(out) :: ibead_beg_r,ibead_end_r
        integer :: nbeads_write,nbeads_write_max
        integer :: natpi_max,iproc, ibead
        integer :: ibead_beg_s,ibead_end_s
        integer :: ibegpi_r,iendpi_r
        integer :: cs,cr,nsend,nsendbis
        integer i,j,k,ii,jj,kk,ierr
        integer, allocatable :: nlocr(:)
        real(r_p), allocatable :: buffers(:,:,:),buffer(:,:,:)
        integer, allocatable :: reqr(:),reqs(:),reqsize(:)


        nbeads_write=polymer%nbeads/nproc_write
        ibead_beg_r=ranktot*nbeads_write+1
        ibead_end_r=(ranktot+1)*nbeads_write
        if(ranktot==nproc_write-1) ibead_end_r=polymer%nbeads
  
        nbeads_write_max=polymer%nbeads-(nproc_write-1)*nbeads_write
  
        !SEND DATA
        allocate(buffers(13,nlocpi*nbeads_write_max
     &    ,nproc_write))
        allocate(reqs(nproc_write))
        allocate(reqsize(nproc_write))
        reqs(:)=MPI_REQUEST_NULL
        cs=0
        do iproc=1,nproc_write
          if(iproc-1==ranktot) CYCLE
          cs=cs+1
          ibead_beg_s=(iproc-1)*nbeads_write+1
          ibead_end_s=iproc*nbeads_write
          if(iproc==nproc_write) ibead_end_s=polymer%nbeads
          do k=1,ibead_end_s-ibead_beg_s+1
            kk=ibead_beg_s+k-1
            do i=1,nlocpi; do j=1,3
              ii=glob(i-1+ilocpi_beg(rank_polymer+1))
              buffers(1,i+nlocpi*(k-1),iproc)=real(ii,r_p)
              buffers(1+j,i+nlocpi*(k-1),iproc)
     &             =polymer%pos(j,ii,kk)
              buffers(1+j+3,i+nlocpi*(k-1),iproc)
     &             =polymer%vel(j,ii,kk)
              buffers(1+j+6,i+nlocpi*(k-1),iproc)
     &             =polymer%forces(j,ii,kk)
              buffers(1+j+9,i+nlocpi*(k-1),iproc)
     &             =polymer%forces_slow(j,ii,kk)
            enddo; enddo
          enddo
          nsend=13*nlocpi*(ibead_end_s-ibead_beg_s+1)
c          write(0,*) ranktot,iproc-1,nbeads_write
c     &       ,(ibead_end_s-ibead_beg_s+1)
          
          call MPI_ISEND(nlocpi,1, MPI_INT
     &      ,iproc-1,1,MPI_COMM_WORLD,reqsize(cs),ierr)
          call MPI_ISEND(buffers(1,1,iproc),nsend, MPI_RPREC
     &      ,iproc-1,0,MPI_COMM_WORLD,reqs(cs),ierr)
        enddo

        !RECEIVE SIZES
        if(ranktot<nproc_write) then
          allocate(reqr(nproctot))
          allocate(nlocr(nproctot))
          nlocr=0
          reqr(:)=MPI_REQUEST_NULL
          cr=0
          do iproc=1,nproctot
            if(iproc-1==ranktot) CYCLE
            cr=cr+1          
            call MPI_IRECV(nlocr(iproc),1, MPI_INT
     &        ,iproc-1,1,MPI_COMM_WORLD,reqr(cr),ierr)
          enddo
        endif

        !WAIT FOR SIZES 
        call MPI_WAITALL(cs,reqsize
     &    ,MPI_STATUSES_IGNORE,ierr)
        deallocate(reqsize)
        if(ranktot<nproc_write) then
          call MPI_WAITALL(cr,reqr
     &     ,MPI_STATUSES_IGNORE,ierr)
        endif
  
        !RECEIVE DATA
        if(ranktot<nproc_write) then
          natpi_max=maxval(nlocr)
          allocate(buffer(13,natpi_max*(ibead_end_r-ibead_beg_r+1)
     &      ,nproctot))
          reqr(:)=MPI_REQUEST_NULL
          cr=0
          do iproc=1,nproctot
            if(iproc-1==ranktot) CYCLE
            cr=cr+1          
            nsend=13*nlocr(iproc)*(ibead_end_r-ibead_beg_r+1)
            call MPI_IRECV(buffer(1,1,iproc),nsend, MPI_RPREC
     &        ,iproc-1,0,MPI_COMM_WORLD,reqr(cr),ierr)
          enddo
        endif
  
        !WAIT FOR DATA
        call MPI_WAITALL(cs,reqs
     &    ,MPI_STATUSES_IGNORE,ierr)
        deallocate(buffers,reqs)
        if(ranktot<nproc_write) then
          call MPI_WAITALL(cr,reqr
     &     ,MPI_STATUSES_IGNORE,ierr)
          deallocate(reqr)
        endif
  
        if(ranktot<nproc_write) then
          do iproc=1,nproctot
            if(iproc-1==ranktot) CYCLE
            do k=1,ibead_end_r-ibead_beg_r+1
              kk=ibead_beg_r+k-1
              do i=1,nlocr(iproc); do j=1,3
                ii=nint(buffer(1,i+nlocr(iproc)*(k-1),iproc))
                polymer%pos(j,ii,kk)
     &             =buffer(1+j,i+nlocr(iproc)*(k-1),iproc)
                polymer%vel(j,ii,kk)
     &             =buffer(1+j+3,i+nlocr(iproc)*(k-1),iproc)
                polymer%forces(j,ii,kk)
     &             =buffer(1+j+6,i+nlocr(iproc)*(k-1),iproc)
                polymer%forces_slow(j,ii,kk)
     &             =buffer(1+j+9,i+nlocr(iproc)*(k-1),iproc)
              enddo; enddo
            enddo
          enddo
          deallocate(buffer)
        endif

        end subroutine comm_for_write_pi
