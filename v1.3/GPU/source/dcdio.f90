!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     "dcdio_write" writes out a set of Cartesian coordinates
!     to an external disk file in the dcd format
!     based on libdcdfort: https://github.com/wesbarnett/dcdfort
!
#include "tinker_macro.h"
submodule(dcdmod) subdcdmod
   implicit none

contains

   subroutine unwrap_dcd_pos(pos,wrapIdx,na,wdim,coordinates)
      use cell
      implicit none
      integer(1) wrapIdx(*)
      integer   ,intent(in):: na, wdim
      real(r_p), intent(in) :: coordinates(*)
      real(4)   ,intent(inout):: pos(*)
      integer    kk,i

66    format(" Unknown option for unwrap_dcd_pos ",I0)
      select case (wdim)
       case (1)
         do i=1,na
            pos(i)= real(coordinates(i)+ int(wrapIdx(4*(i-1)+wdim))*xcell,4)
         end do
       case (2)
         do i=1,na
            pos(i)= real(coordinates(i)+ int(wrapIdx(4*(i-1)+wdim))*ycell,4)
         end do
       case (3)
         do i=1,na
            pos(i)= real(coordinates(i)+ int(wrapIdx(4*(i-1)+wdim))*zcell,4)
         end do
       case default
         print 66, wdim
         __TINKER_FATAL__
      end select

   end subroutine

   module subroutine dcdio_write(dcdinfo,trajectory_file,istep,dt,append,&
      &x_frame,y_frame,z_frame,status_value)
      use atmtyp
      use atoms      ,only: pbcunwrap,pbcWrapIdx
      use boxes
      use inform
      use replicas
      use iso_c_binding, only: C_NULL_CHAR
      use status, only: clear_status, failed, set_status, status_invalid_argument, status_ok, status_t
      implicit none
      type(dcdinfo_t), intent(inout) :: dcdinfo
      character(*), intent(in) :: trajectory_file
      integer, intent(in) :: istep
      real(r_p), intent(in) :: dt
      logical, intent(in) :: append
      real(r_p), intent(in) :: x_frame(:),y_frame(:),z_frame(:)
      type(status_t), intent(out) :: status_value
      type(status_t) :: cleanup_status
      integer i,j,k,idcd
      integer freeunit
      integer(kind=4) :: coord_size,ios
      real(kind=4), allocatable :: posw(:)
      real(kind=8) :: box(6)
      logical init,exist
      character*240 dcdfile
      character (len=79) :: info1,info2
      character (len=8) :: date
      character (len=10) :: time
      logical append_positioned

      call clear_status(status_value)

      dcdinfo%natoms     = size(x_frame)
      coord_size = 4*dcdinfo%natoms
      dcdinfo%timestep   = real(dt,4)

      dcdinfo%filename = trim(trajectory_file)
      init       = (dcdinfo%idcd.eq.0)
      idcd = dcdinfo%idcd
      append_positioned = .false.
!
      if (init) then
!
!     check if dcd trajectory file exists
!
         inquire (file=trim(dcdinfo%filename),exist=exist)
         if (exist .and. append) then
            open (newunit=idcd,file=trim(dcdinfo%filename)&
               &,access="stream",status='old',position='append',iostat=ios)
            if (ios /= 0) then
               call set_status(status_value,status_invalid_argument,'Unable to open DCD trajectory')
               return
            end if
            dcdinfo%idcd=idcd
!
!     check header
!
            call dcdio_read_header(dcdinfo,status_value)
            if (failed(status_value)) then
               call dcdio_close(dcdinfo,cleanup_status)
               if (failed(cleanup_status)) status_value = cleanup_status
               return
            end if
            append_positioned = .true.
!
!       nevery: frequency output in timesteps
!
            dcdinfo%nevery = iwrite
!
         else
!
!     open dcd trajectory file
!
            open(newunit=idcd, file=trim(dcdinfo%filename)&
               &, form="unformatted",access="stream", status="replace",iostat=ios)
            if (ios /= 0) then
               call set_status(status_value,status_invalid_argument,'Unable to open DCD trajectory')
               return
            end if

            dcdinfo%idcd=idcd
!
!     create the header
!
            do i = 1, 79
               info1(i:i) = ' '
               info2(i:i) = ' '
            end do
            info1 = "Created by Tinker-HP development version"
            call date_and_time(date=date,time=time)
            info2 = "Created on "//date//" "//time

            dcdinfo%nframes = 0
!
!       istart: first timestep
!
            dcdinfo%istart = iwrite
            dcdinfo%iend = dcdinfo%istart
!
!       nevery: frequency output in timesteps
!
            dcdinfo%nevery = iwrite

            write(idcd,iostat=ios) 84
            if (ios /= 0) goto 900
            write(idcd,iostat=ios) "CORD"
            if (ios /= 0) goto 900

            inquire(unit=idcd, pos=dcdinfo%nframes_pos)

! Number of snapshots in file
            write(idcd,iostat=ios) dcdinfo%nframes
            if (ios /= 0) goto 900
!
! Timestep of first snapshot
            write(idcd,iostat=ios) dcdinfo%istart
            if (ios /= 0) goto 900
!
! Save snapshots every this many steps
            write(idcd,iostat=ios) dcdinfo%nevery
            if (ios /= 0) goto 900
!
            inquire(unit=idcd, pos=dcdinfo%iend_pos)
! Timestep ! of ! last ! snapshot
            write(idcd,iostat=ios) dcdinfo%iend
            if (ios /= 0) goto 900

            do i = 1, 5
               write(idcd,iostat=ios) 0
               if (ios /= 0) goto 900
            end do

! Simulation timestep
            write(idcd,iostat=ios) dcdinfo%timestep
            if (ios /= 0) goto 900

! Has unit cell
            write(idcd,iostat=ios) 1
            if (ios /= 0) goto 900

            do i = 1, 8
               write(idcd,iostat=ios) 0
               if (ios /= 0) goto 900
            end do
            ! Pretend to be CHARMM version 24
            write(idcd,iostat=ios) 24
            if (ios /= 0) goto 900
            write(idcd,iostat=ios) 84
            if (ios /= 0) goto 900
            write(idcd,iostat=ios) 164
            if (ios /= 0) goto 900

            write(idcd,iostat=ios) 2
            if (ios /= 0) goto 900
            write(idcd,iostat=ios) info1//C_NULL_CHAR
            if (ios /= 0) goto 900
            write(idcd,iostat=ios) info2//C_NULL_CHAR
            if (ios /= 0) goto 900

            write(idcd,iostat=ios) 164
            if (ios /= 0) goto 900
            write(idcd,iostat=ios) 4
            if (ios /= 0) goto 900

            ! Number of atoms in each snapshot
            write(idcd,iostat=ios) dcdinfo%natoms
            if (ios /= 0) goto 900
            write(idcd,iostat=ios) 4
            if (ios /= 0) goto 900

            flush(idcd,iostat=ios)
            if (ios /= 0) goto 900
         end if
      end if

      allocate (posw(dcdinfo%natoms))

      if (append_positioned) then
         write(idcd,pos=dcdinfo%curr_pos,iostat=ios) 48
      else
         write(idcd,iostat=ios) 48
      end if
      if (ios /= 0) goto 900

      box(1) = xbox
      write(idcd,iostat=ios) box(1) ! A
      if (ios /= 0) goto 900
      box(6) = gamma
      write(idcd,iostat=ios) box(6) ! gamma
      if (ios /= 0) goto 900
      box(2) = ybox
      write(idcd,iostat=ios) box(2) ! B
      if (ios /= 0) goto 900
      box(5) = beta
      write(idcd,iostat=ios) box(5) ! beta
      if (ios /= 0) goto 900
      box(4) = alpha
      write(idcd,iostat=ios) box(4) ! alpha
      if (ios /= 0) goto 900
      box(3) = zbox
      write(idcd,iostat=ios) box(3) ! C
      if (ios /= 0) goto 900


      write(idcd,iostat=ios) 48
      if (ios /= 0) goto 900
      write(idcd,iostat=ios) coord_size
      if (ios /= 0) goto 900

      if (pbcunwrap) then
         call unwrap_dcd_pos(posw,pbcWrapIdx,dcdinfo%natoms,1,x_frame)
      else
         do i = 1,dcdinfo%natoms; posw(i) = real(x_frame(i),4); end do
      end if
      write(idcd,iostat=ios) posw(:)
      if (ios /= 0) goto 900

      write(idcd,iostat=ios) coord_size
      if (ios /= 0) goto 900
      write(idcd,iostat=ios) coord_size
      if (ios /= 0) goto 900

      if (pbcunwrap) then
         call unwrap_dcd_pos(posw,pbcWrapIdx,dcdinfo%natoms,2,y_frame)
      else
         do i = 1,dcdinfo%natoms; posw(i) = real(y_frame(i),4); end do
      end if
      write(idcd,iostat=ios) posw(:)
      if (ios /= 0) goto 900

      write(idcd,iostat=ios) coord_size
      if (ios /= 0) goto 900
      write(idcd,iostat=ios) coord_size
      if (ios /= 0) goto 900

      if (pbcunwrap) then
         call unwrap_dcd_pos(posw,pbcWrapIdx,dcdinfo%natoms,3,z_frame)
      else
         do i = 1,dcdinfo%natoms; posw(i) = real(z_frame(i),4); end do
      end if
      write(idcd,iostat=ios) posw(:)
      if (ios /= 0) goto 900

      write(idcd,iostat=ios) coord_size
      if (ios /= 0) goto 900

      inquire(unit=idcd, pos=dcdinfo%curr_pos)

      dcdinfo%nframes = dcdinfo%nframes+1
      dcdinfo%iend = dcdinfo%iend + dcdinfo%nevery

      ! Go back and update header
      write(idcd, pos=9, iostat=ios) dcdinfo%nframes
      if (ios /= 0) goto 900
      write(idcd, pos=21, iostat=ios) dcdinfo%iend
      if (ios /= 0) goto 900
      write(idcd, pos=dcdinfo%curr_pos, iostat=ios)
      if (ios /= 0) goto 900

      flush(idcd,iostat=ios)
      if (ios /= 0) goto 900

      deallocate (posw)
      call set_status(status_value,status_ok,'DCD trajectory frame written')
      return
900   continue
      if (allocated(posw)) deallocate(posw)
      call set_status(status_value,status_invalid_argument,'Unable to write DCD trajectory')
      return
   end

   module subroutine dcdio_read_header(dcdinfo,status_value)
      use status, only: clear_status, set_status, status_invalid_argument, &
         &status_ok, status_t
      implicit none
      type(dcdinfo_t), intent(inout) :: dcdinfo
      type(status_t), intent(out) :: status_value
      integer(kind=4) :: magic_number, ntitle, title_marker, atom_marker
      integer(kind=4) :: atom_count, charmm_version, has_extra_block
      integer(kind=4) :: ios
      integer(kind=8) :: first_frame, expected_size, file_size
      character(len=4) :: magic_string

      call clear_status(status_value)
      read(dcdinfo%idcd,pos=1,iostat=ios) magic_number
      if (ios /= 0) goto 910
      read(dcdinfo%idcd,pos=5,iostat=ios) magic_string
      if (ios /= 0 .or. magic_number /= 84 .or. magic_string /= 'CORD') goto 910
      read(dcdinfo%idcd,pos=9,iostat=ios) dcdinfo%nframes,dcdinfo%istart,&
         &dcdinfo%nevery,dcdinfo%iend
      if (ios /= 0) goto 910
      read(dcdinfo%idcd,pos=45,iostat=ios) dcdinfo%timestep
      if (ios /= 0) goto 910
      read(dcdinfo%idcd,pos=49,iostat=ios) has_extra_block
      if (ios /= 0 .or. has_extra_block /= 1) goto 910
      read(dcdinfo%idcd,pos=85,iostat=ios) charmm_version
      if (ios /= 0 .or. charmm_version == 0) goto 910
      read(dcdinfo%idcd,pos=97,iostat=ios) ntitle
      if (ios /= 0 .or. ntitle < 0) goto 910

      first_frame = 101 + int(80,kind=8)*int(ntitle,kind=8) + 16
      read(dcdinfo%idcd,pos=first_frame-16,iostat=ios) title_marker
      if (ios /= 0 .or. title_marker /= 164) goto 910
      read(dcdinfo%idcd,pos=first_frame-12,iostat=ios) atom_marker
      if (ios /= 0 .or. atom_marker /= 4) goto 910
      read(dcdinfo%idcd,pos=first_frame-8,iostat=ios) atom_count
      if (ios /= 0 .or. atom_count <= 0 .or. atom_count /= dcdinfo%natoms) goto 910
      read(dcdinfo%idcd,pos=first_frame-4,iostat=ios) atom_marker
      if (ios /= 0 .or. atom_marker /= 4) goto 910
      if (dcdinfo%nframes > 0) then
         read(dcdinfo%idcd,pos=first_frame,iostat=ios) atom_marker
         if (ios /= 0 .or. atom_marker /= 48) goto 910
      end if

      dcdinfo%framesize = dcdinfo%natoms*12 + 80
      inquire(unit=dcdinfo%idcd,size=file_size,iostat=ios)
      if (ios /= 0) goto 910
      expected_size = first_frame - 1 + dcdinfo%framesize*dcdinfo%nframes
      if (dcdinfo%nframes < 0 .or. file_size /= expected_size) goto 910
      dcdinfo%nframes_pos = 9
      dcdinfo%iend_pos = 21
      dcdinfo%curr_pos = file_size + 1
      call set_status(status_value,status_ok,'DCD append header read')
      return
910   continue
      call set_status(status_value,status_invalid_argument,&
         &'Unable to read DCD append header')
   end subroutine dcdio_read_header

   module subroutine dcdio_close(dcdinfo,status_value)
      use status, only: clear_status, set_status, status_invalid_argument, status_ok, status_t
      implicit none
      type(dcdinfo_t), intent(inout) :: dcdinfo
      type(status_t), intent(out) :: status_value
      integer :: ios
      logical :: opened

      call clear_status(status_value)
      opened = .false.
      if (dcdinfo%idcd /= 0) inquire(unit=dcdinfo%idcd,opened=opened)
      if (opened) then
         close(dcdinfo%idcd,iostat=ios)
         if (ios /= 0) then
            call set_status(status_value,status_invalid_argument,&
               &'Unable to close DCD trajectory')
            return
         end if
      end if
      if (allocated(dcdinfo%titlesdcd)) deallocate(dcdinfo%titlesdcd)
      dcdinfo%idcd = 0
      call set_status(status_value,status_ok,'DCD trajectory closed')
   end subroutine dcdio_close
!
!     subroutine dcdfile_open: opens file to read trajectory from
!
   module subroutine dcdfile_open(dcdinfo,dcdfile)
      use domdec
      use files
      use inform
      use iounit
      implicit none
      type(dcdinfo_t), intent(inout) :: dcdinfo
      character(*), intent(in) :: dcdfile
      character(len=*), parameter :: magic_string = "CORD"
      integer(kind=4), parameter :: magic_number = 84
      integer(kind=4) :: line1, charmm_version, has_extra_block,&
         &four_dimensions
      integer :: freeunit,idcd
      character(len=4) :: line2
      logical :: ex
!
      ! Does file exist?
      dcdinfo%filename = trim(dcdfile)
      inquire(file=trim(dcdfile), exist=ex, size=dcdinfo%filesize)
      if (.not.ex) then
         if (rank.eq.0) write(iout,*) 'The specified DCD file ',&
            &trim(dcdfile), ' does not exist'
         call fatal
      end if

      ! Open file in native endinness
      idcd = freeunit()
      open(unit=idcd, file=trim(dcdfile), form="unformatted",&
         &access="stream", status="old")
      dcdinfo%idcd = idcd

      ! Read in magic number of magic string
      read(idcd,pos=1) line1
      read(idcd) line2

      ! Ensure the magic number and string are correct
      if (line1 .ne. magic_number .or. line2 .ne. magic_string) then
         if (rank.eq.0) write(iout,*) 'This DCD file format is not ',&
            &'supported, or the file header is corrupt.'
         call fatal
      end if

      ! Check if the file identifies as CHARMM (LAMMPS pretends to be CHARMM v. 24)
      read(idcd, pos=85) charmm_version
      if (charmm_version .eq. 0) then
         if (rank.eq.0) write(iout,*) 'DCD file indicated it is ',&
            &'not CHARMM. Only CHARMM-style DCD files are supported'
         call fatal
      end if

      ! We only support files with the extra unitcell block
      read(idcd, pos=49) has_extra_block
      if (has_extra_block .ne. 1) then
         if (rank.eq.0) write(iout,*) 'DCD file indicated it does not ',&
            &'have unit information. Only DCD files with unit cell ',&
            &'information are supported'
         call fatal
      end if

      ! We don't support files with four dimensions
      read(idcd) four_dimensions
      if (four_dimensions .eq. 1) then
         if (rank.eq.0) write(iout,*) 'DCD file indicates it has four ',&
            &'dimensions Only DCD files with three dimensions'
         call fatal
      end if
      return
   end
!
!     read header of the dcd file
!
   module subroutine dcdfile_read_header(dcdinfo,dowrite,verbose)
      use atoms,only :n
      use domdec
      use iounit
      use iso_c_binding, only: C_NULL_CHAR
      implicit none
      type(dcdinfo_t), intent(inout) :: dcdinfo
      logical, intent(in) :: dowrite
      logical, intent(in), optional :: verbose
      integer(kind=4) :: i, ntitle, m, dummy, pos
      integer(kind=8) :: nframes2,newpos
      character(len=80) :: title_string
      integer :: idcd
      logical :: verbose_

      verbose_=.true.
      if(present(verbose)) verbose_=verbose

      idcd = dcdinfo%idcd
      read(idcd, pos=9) dcdinfo%nframes, dcdinfo%istart&
         &, dcdinfo%nevery, dcdinfo%iend

      read(idcd, pos=45) dcdinfo%timestep

      read(idcd, pos=97) ntitle
      if (ntitle > 0 .and. verbose_) then
         if (rank.eq.0) write(iout,*) "The following titles were found:"
      end if
      if (allocated(dcdinfo%titlesdcd)) deallocate(dcdinfo%titlesdcd)
      allocate(dcdinfo%titlesdcd(ntitle))
      do i = 1, ntitle
         read(idcd) title_string
         m = 1
         do while (m .le. 80 .and. title_string(m:m) .ne. C_NULL_CHAR)
            m = m + 1
         end do
         dcdinfo%titlesdcd(i) = trim(title_string(1:m-1))
         if (rank.eq.0 .and. verbose_)&
            &write(iout,*) dcdinfo%titlesdcd(i)
      end do

      read(idcd) dummy, dummy
      if (dummy .ne. 4) then
         if (rank.eq.0) write(iout,*) 'This DCD file format is not ',&
            &'supported or the file header is corrupt.'
         call fatal
      end if

      ! Number of atoms in each snapshot
      read(idcd) dcdinfo%natoms, dummy
!
!    check if the number of atoms is the one of the "topology" xyz file
!
      if (dcdinfo%natoms.ne.n) then
         if (rank.eq.0) write(iout,*) 'number of atoms in the ',&
            &'trajectory is not the same as in the xyz file'
         if (rank.eq.0) write(iout,*) 'respectively ',dcdinfo%natoms,n
         call fatal
      end if

      if (dummy .ne. 4) then
         if (rank.eq.0) write(iout,*) 'This DCD file format is not ',&
            &'supported or the file header is corrupt.'
         call fatal
      end if

      inquire(unit=idcd, pos=pos)
      pos = pos - 1

      ! Each frame has natoms*3 (4 bytes each) = natoms*12
      ! plus 6 box dimensions (8 bytes each) = 48
      ! Additionally there are 32 bytes of file information in each frame
      dcdinfo%framesize = dcdinfo%natoms*12 + 80
      ! Header is typically 276 bytes, but inquire gives us exact size
      ! Only check size if we are reading and not if we are appending a
      ! file
      if (dowrite) then
         ! Where are we?
         inquire(unit=idcd, pos=pos)
         newpos = pos + dcdinfo%framesize*dcdinfo%nframes - 4
         read(idcd, pos=newpos) dummy
         return
      end if
      nframes2 = (dcdinfo%filesize-pos)/dcdinfo%framesize
      if ( nframes2 .ne. dcdinfo%nframes) then
         write(iout,'(a,i0,a,i0,a)') "WARNING:&
            &         Header indicates ", dcdinfo%nframes,&
            &" frames, but file size indicates ", nframes2, "."
         dcdinfo%nframes = int(nframes2)
      end if
   end subroutine dcdfile_read_header
!
!      Closes DCD file
!
   module subroutine dcdfile_close(dcdinfo)
      implicit none
      type(dcdinfo_t), intent(inout) :: dcdinfo
      logical :: opened

      if (allocated(dcdinfo%titlesdcd)) deallocate(dcdinfo%titlesdcd)
      if (allocated(dcdinfo%filename)) deallocate(dcdinfo%filename)
      inquire(unit=dcdinfo%idcd,opened=opened)
      if (opened) close(dcdinfo%idcd)
      dcdinfo%idcd = 0
   end subroutine dcdfile_close

   !> @Reads next frame into memory
   module subroutine dcdfile_read_next(dcdinfo)
!      use atoms, only : x,y,z
      use atomsMirror,only:xm=>x,ym=>y,zm=>z&
         &,reCast_position,atomsmirror_init&
         &,download_position
      use boxes, only : xbox,ybox,zbox,alpha,beta,gamma
      use domdec
      use inform,only: abort
      use iounit
      use tinheader ,only:ti_p,re_p
      implicit none
      type(dcdinfo_t), intent(inout) :: dcdinfo
      integer(kind=4) :: dummy(6), nbytes, ndims, i
      integer(kind=8) :: pos,newpos
      real(kind=4), allocatable :: xd(:),yd(:),zd(:)
      real(kind=8) :: box(6)
      character*120 :: linepos
      integer :: idcd
!
!     abort if we are past the end of the file
!
      ! Where are we?
      idcd = dcdinfo%idcd
      inquire(unit=idcd, pos=pos)

      ! We subtract 4 bytes so that the next read of the 4-byte integer will line things up properly for the next read
      newpos = pos - 4
      if (newpos.ge.(dcdinfo%framesize*dcdinfo%nframes)) then
         abort = .true.
         return
      end if

      allocate(xd(dcdinfo%natoms))
      allocate(yd(dcdinfo%natoms))
      allocate(zd(dcdinfo%natoms))
      nbytes = size(xd)*4
      ndims  = 3

      if (ndims /= 3) then
         if (rank.eq.0) write(iout,*) 'Number of dimensions of xyz ',&
            &'array is incorrect'
         call fatal
      end if

      read(idcd) dummy(1)
      if (dummy(1) /= 48) then
         if (rank.eq.0) write(iout,*) 'Problem reading in DCD snapshot'
         call fatal
      end if

      !            A       gamma   B       beta    alpha   C
      read(idcd) box(1), box(6), box(2), box(5), box(4), box(3)
      if (box(1) < 0 .or. box(2) < 0 .or. box(3) < 0) then
         if (rank.eq.0) write(iout,*) 'Problm reading in DCD snapshot',&
            &' box dimensions.'
         call fatal
      end if
!
!     put size of the box in global Tinker variables
!
      xbox = box(1)
      ybox = box(2)
      zbox = box(3)
      alpha = box(4)
      beta  = box(5)
      gamma = box(6)


      ! 48, then no. of bytes for x coordinates, x coordinates (repeat for y and z coordinates)
      read(idcd) dummy(1:2), xd(:), dummy(3:4), yd(:), dummy(5:6), zd(:)


      if (dummy(1) /= 48) then
         if (rank.eq.0) write(iout,*) 'Problem reading in DCD snapshot ',&
            &'coordinates'
         call fatal
      end if

      do i = 2, 6
         if (dummy(i) /= nbytes) then
            write(iout,*) 'Number of bytes in DCD snapshot is incorrect',&
               &' for size of xyz array passed'
            call fatal
         end if
      end do

      read(idcd) dummy(1)
      if (dummy(1) .ne. nbytes) then
         if (rank.eq.0) write(iout,*) 'Problem reading in DCD snapshot'
         call fatal
      end if
!
!     put coordinates in global Tinker arrays
!
!     for reproducibility: first write truncated value in a char and then read them
!
      do i = 1, dcdinfo%natoms
         write(linepos,'(3F16.6)') xd(i),yd(i),zd(i)
         read(linepos,'(3F16.6)') xm(i),ym(i),zm(i)
      end do

!$acc update device(xm,ym,zm)
      if (t_p.ne.r_p) then
         call reCast_position
!$acc wait
         call download_position
      end if

      deallocate(xd,yd,zd)
   end subroutine dcdfile_read_next
!
!     Skips reading n-1 frames into memory
!
   module subroutine dcdfile_skip_next(dcdinfo,n)
      implicit none
      type(dcdinfo_t), intent(inout) :: dcdinfo
      integer(kind=4) :: dummy
      integer(kind=8) :: pos, newpos
      integer(kind=4), intent(in), optional :: n

      ! Where are we?
      inquire(unit=dcdinfo%idcd, pos=pos)

      ! We subtract 4 bytes so that the next read of the 4-byte integer will line things up properly for the next read
      if (.not. present(n)) then
         newpos = pos + dcdinfo%framesize - 4
      else
         newpos = pos + dcdinfo%framesize*n - 4
      end if

      read(dcdinfo%idcd, pos=newpos) dummy

   end subroutine dcdfile_skip_next

end submodule
