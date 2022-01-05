c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "dcdio_write" writes out a set of Cartesian coordinates
c     to an external disk file in the dcd format
c     based on libdcdfort: https://github.com/wesbarnett/dcdfort
c
      subroutine dcdio_write(istep)
      use atoms
      use atmtyp
      use boxes
      use dcdmod
      use files
      use inform
      use, intrinsic :: iso_c_binding
      implicit none
      integer i,istep
      integer freeunit
      integer(kind=4) :: coord_size
      real(kind=4), allocatable :: posw(:,:)
      real(kind=8) :: box(6)
      logical init,exist
      character*240 dcdfile
      character (len=79) :: info1,info2
      character (len=8) :: date
      character (len=10) :: time
c
      natoms = n
      coord_size = 4*natoms
      timestep = istep

      dcdfile = filename(1:leng)//'.dcd'
      init = (istep.eq.iwrite)
c
      if (init) then
        idcd = freeunit ()
c
c     check if dcd trajectory file exists
c
        inquire (file=trim(dcdfile),exist=exist)
        if (exist) then
          open (unit=idcd,file=trim(dcdfile),access="stream",
     $      status='old',position='append')
c
c     check header
c
        call dcdfile_read_header(.true.)
c
c       nevery: frequency output in timesteps
c
        nevery = iwrite 
c  
        else
c
c     open dcd trajectory file
c
        open(newunit=idcd, file=trim(dcdfile), form="unformatted",
     $ access="stream", status="replace")
c  
c     create the header
c
        do i = 1, 79
          info1(i:i) = ' '
          info2(i:i) = ' '
        end do
        info1 = "Created by Tinker-HP version 1.2"
        call date_and_time(date=date,time=time)
        info2 = "Created on "//date//" "//time

        nframes = 0
c
c       istart: first timestep
c
        istart = iwrite
        iend = istart
c
c       nevery: frequency output in timesteps
c
        nevery = iwrite 

        write(idcd) 84
        write(idcd) "CORD"

        inquire(unit=idcd, pos=nframes_pos)

! Number of snapshots in file
        write(idcd) nframes
c
! Timestep of first snapshot
        write(idcd) istart
c
! Save snapshots every this many steps
        write(idcd) nevery
c
        inquire(unit=idcd, pos=iend_pos)
! Timestep ! of ! last ! snapshot
         write(idcd) iend

         do i = 1, 5
           write(idcd) 0
         end do

! Simulation timestep
         write(idcd) timestep

! Has unit cell
        write(idcd) 1

        do i = 1, 8
          write(idcd) 0
        end do
         ! Pretend to be CHARMM version 24
        write(idcd) 24
        write(idcd) 84
        write(idcd) 164
         
        write(idcd) 2
        write(idcd) info1//C_NULL_CHAR
        write(idcd) info2//C_NULL_CHAR

        write(idcd) 164
        write(idcd) 4

        ! Number of atoms in each snapshot
        write(idcd) natoms
        write(idcd) 4

        flush(idcd)
      end if
      end if

      allocate (posw(3,natoms))
      posw(1,1:natoms) = real(x(1:natoms),4)
      posw(2,1:natoms) = real(y(1:natoms),4)
      posw(3,1:natoms) = real(z(1:natoms),4)

      write(idcd) 48
    
      box(1) = xbox
      write(idcd) box(1) ! A
      box(6) = gamma
      write(idcd) box(6) ! gamma
      box(2) = ybox
      write(idcd) box(2) ! B
      box(5) = beta
      write(idcd) box(5) ! beta
      box(4) = alpha
      write(idcd) box(4) ! alpha
      box(3) = zbox
      write(idcd) box(3) ! C

      write(idcd) 48
      write(idcd) coord_size
      
      write(idcd) posw(1,:)

      write(idcd) coord_size
      write(idcd) coord_size

      write(idcd) posw(2,:)

      write(idcd) coord_size
      write(idcd) coord_size

      write(idcd) posw(3,:)

      write(idcd) coord_size

      inquire(unit=idcd, pos=curr_pos)

      nframes = nframes+1
      iend = iend + nevery

      ! Go back and update header
      write(idcd, pos=9) nframes
      write(idcd, pos=21) iend
      write(idcd, pos=curr_pos)

      flush(idcd)

      deallocate (posw)
      return
      end
c
c     subroutine dcdfile_open: opens file to read trajectory from
c
      subroutine dcdfile_open(dcdfile)
      use dcdmod
      use domdec
      use files
      use inform
      use iounit
      implicit none
      character*240 dcdfile
      character(len=*), parameter :: magic_string = "CORD"
      integer(kind=4), parameter :: magic_number = 84
      integer(kind=4) :: line1, charmm_version, has_extra_block, 
     $  four_dimensions
      integer :: freeunit
      character(len=4) :: line2
      logical :: ex
c
      ! Does file exist?
      inquire(file=trim(dcdfile), exist=ex, size=filesize)
      if (ex .eqv. .false.) then
        if (rank.eq.0) write(iout,*) 'The specified DCD file ',
     $    trim(dcdfile), ' does not exist'
        call fatal
      end if

      ! Open file in native endinness
      idcd = freeunit ()
      open(newunit=idcd, file=trim(dcdfile), form="unformatted",
     $      access="stream", status="old")

      ! Read in magic number of magic string
      read(idcd,pos=1) line1
      read(idcd) line2

      ! Ensure the magic number and string are correct
      if (line1 .ne. magic_number .or. line2 .ne. magic_string) then
        if (rank.eq.0) write(iout,*) 'This DCD file format is not ',
     $  'supported, or the file header is corrupt.'
        call fatal
      end if

      ! Check if the file identifies as CHARMM (LAMMPS pretends to be CHARMM v. 24)
      read(idcd, pos=85) charmm_version
      if (charmm_version .eq. 0) then
        if (rank.eq.0) write(iout,*) 'DCD file indicated it is ',
     $    'not CHARMM. Only CHARMM-style DCD files are supported'
        call fatal 
      end if

      ! We only support files with the extra unitcell block
      read(idcd, pos=49) has_extra_block
      if (has_extra_block .ne. 1) then
        if (rank.eq.0) write(iout,*) 'DCD file indicated it does not ',
     $    'have unit information. Only DCD files with unit cell ',
     $    'information are supported'
          call fatal
      end if

      ! We don't support files with four dimensions
      read(idcd) four_dimensions
      if (four_dimensions .eq. 1) then
        if (rank.eq.0) write(iout,*) 'DCD file indicates it has four ',
     $    'dimensions Only DCD files with three dimensions'
          call fatal 
      end if
      return
      end 
c
c     read header of the dcd file
c
      subroutine dcdfile_read_header(dowrite)
      use atoms
      use dcdmod
      use domdec
      use iounit
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=4) :: i, ntitle, m, dummy, pos
      integer(kind=8) :: nframes2,newpos
      character(len=80) :: title_string
      logical :: dowrite

      read(idcd, pos=9) nframes, istart, nevery, iend

      read(idcd, pos=45) timestep

      read(idcd, pos=97) ntitle
      if (ntitle > 0) then
        if (rank.eq.0) write(iout,*) "The following titles were found:"
      end if
      if (allocated(titlesdcd)) deallocate(titlesdcd)
      allocate(titlesdcd(ntitle))
      do i = 1, ntitle
          read(idcd) title_string
          m = 1
          do while (m .le. 80 .and. title_string(m:m) .ne. C_NULL_CHAR)
              m = m + 1
          end do
          titlesdcd(i) = trim(title_string(1:m-1))
          if (rank.eq.0) write(iout,*) titlesdcd(i)
      end do

      read(idcd) dummy, dummy
      if (dummy .ne. 4) then
        if (rank.eq.0) write(iout,*) 'This DCD file format is not ',
     $      'supported or the file header is corrupt.'
        call fatal
      end if

      ! Number of atoms in each snapshot
      read(idcd) natoms, dummy
c
c    check if the number of atoms is the one of the "topology" xyz file
c
      if (natoms.ne.n) then
        if (rank.eq.0) write(iout,*) 'number of atoms in the ',
     $   'trajectory is not the same as in the xyz file'
        if (rank.eq.0) write(iout,*) 'respectively ',natoms,n
        call fatal
      end if

      if (dummy .ne. 4) then
        if (rank.eq.0) write(iout,*) 'This DCD file format is not ',
     $   'supported or the file header is corrupt.'
        call fatal
      end if

      inquire(unit=idcd, pos=pos)
      pos = pos - 1

      ! Each frame has natoms*3 (4 bytes each) = natoms*12
      ! plus 6 box dimensions (8 bytes each) = 48
      ! Additionally there are 32 bytes of file information in each frame
      framesize = natoms*12 + 80
      ! Header is typically 276 bytes, but inquire gives us exact size
      ! Only check size if we are reading and not if we are appending a
      ! file 
      if (dowrite) then
      ! Where are we?
        inquire(unit=idcd, pos=pos)
        newpos = pos + framesize*nframes - 4
        read(idcd, pos=newpos) dummy
        return
      end if
      nframes2 = (filesize-pos)/framesize
      if ( nframes2 .ne. nframes) then
          write(iout,'(a,i0,a,i0,a)') "WARNING:
     $         Header indicates ", nframes, 
     $         " frames, but file size indicates ", nframes2, "." 
          nframes = int(nframes2)
      end if
      end subroutine dcdfile_read_header
c
c      Closes DCD file
c
      subroutine dcdfile_close
      use dcdmod
      implicit none
      deallocate(titlesdcd)
      close(idcd)
      end subroutine dcdfile_close

      !> @Reads next frame into memory
      subroutine dcdfile_read_next
      use atoms
      use boxes
      use dcdmod
      use domdec
      use inform
      use iounit
      implicit none
      real(kind=4), allocatable :: xyz(:,:)
      real(kind=8) :: box(6)
      integer(kind=8) :: pos,newpos
      integer(kind=4) :: dummy(6), nbytes, ndims, i
      character*120 :: linepos
c
c     abort if we are past the end of the file
c
      ! Where are we?
      inquire(unit=idcd, pos=pos)

      ! We subtract 4 bytes so that the next read of the 4-byte integer will line things up properly for the next read
      newpos = pos - 4
      if (newpos.ge.(framesize*nframes)) then
        abort = .true.
        return
      end if

      allocate(xyz(natoms,3))
      nbytes = size(xyz,1)*4
      ndims = size(xyz,2)

      if (ndims /= 3) then
        if (rank.eq.0) write(iout,*) 'Number of dimensions of xyz ',
     $  'array is incorrect'
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
        if (rank.eq.0) write(iout,*) 'Problm reading in DCD snapshot',
     $  ' box dimensions.'
        call fatal
      end if
c
c     put size of the box in global Tinker variables
c
      xbox = box(1)
      ybox = box(2)
      zbox = box(3)
      alpha = box(4)
      beta  = box(5)
      gamma = box(6)
      

      ! 48, then no. of bytes for x coordinates, x coordinates (repeat for y and z coordinates)
      read(idcd) dummy(1:2), xyz(:,1), dummy(3:4), xyz(:,2), dummy(5:6),
     $   xyz(:,3)
     

      if (dummy(1) /= 48) then
        if (rank.eq.0) write(iout,*) 'Problem reading in DCD snapshot ',
     $  'coordinates'
        call fatal
      end if

      do i = 2, 6
        if (dummy(i) /= nbytes) then
          write(iout,*) 'Number of bytes in DCD snapshot is incorrect',
     $    ' for size of xyz array passed'
          call fatal 
        end if
      end do

      read(idcd) dummy(1)
      if (dummy(1) .ne. nbytes) then
        if (rank.eq.0) write(iout,*) 'Problem reading in DCD snapshot'
        call fatal
      end if
c
c     put coordinates in global Tinker arrays
c
c     for reproducibility: first write truncated value in a char and then read them
c
      do i = 1, natoms
        write(linepos,'(3F16.6)') xyz(i,1),xyz(i,2),xyz(i,3)
        read(linepos,'(3F16.6)') x(i),y(i),z(i)
      end do

      deallocate(xyz)
      end subroutine dcdfile_read_next
c
c     Skips reading n-1 frames into memory
c
      subroutine dcdfile_skip_next(n)
      use dcdmod
      use inform
      implicit none
      integer(kind=4) :: dummy
      integer(kind=8) :: pos, newpos
      integer(kind=4), intent(in) :: n
   
      ! Where are we?
      inquire(unit=idcd, pos=pos)

      ! We subtract 4 bytes so that the next read of the 4-byte integer will line things up properly for the next read
c      if (.not. present(n)) then
c          newpos = pos + framesize - 4
c      else
          newpos = pos + framesize*n - 4
c      end if

      read(idcd, pos=newpos) dummy
      end subroutine dcdfile_skip_next

