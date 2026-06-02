!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module dcdmod  --  dcd input/output (I/O) global variables  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     idcd  !<dcd file unit
!     nframes_pos !<current number of frames
!     iend_pos !<position of the end of the file
!     curr_pos !<current position in the file
!     nframes !<number of frames in the file
!     istart !<timestep of last snapshot
!     iend !<timestep of last snapshot
!     nevery !<number of steps between outputs
!
!     reader:
!     filesize !<size of the dcd file to be read
!     framesize !<size of a frame
!     titlesdcd: !<titles of the trajectory
!     timestep: !<timestep of the trajectory
!     natoms: !<number of atoms in the trajectory
!
!     dcdio: write/read trajectory in the dcd format
!
module dcdmod
   implicit none

   type dcdinfo_t
      integer(kind=4) :: idcd !<dcd file unit
      integer(kind=4) :: nframes !<number of frames in the file
      integer(kind=4) :: istart !<timestep of last snapshot
      integer(kind=4) :: iend !<timestep of last snapshot
      integer(kind=4) :: nevery !<number of steps between outputs
      integer(kind=8) :: nframes_pos !<current number of frames
      integer(kind=8) :: iend_pos !<position of the end of the file
      integer(kind=8) :: curr_pos !<current position in the file
      integer(kind=8) :: filesize !<size of the dcd file to be read
      integer(kind=8) :: framesize !<size of a frame
      integer(kind=4) :: natoms !<number of atoms in the trajectory
      real(kind=4) :: timestep !<timestep of the trajectory
      character(len=80), allocatable :: titlesdcd(:) !<titles of the trajectory
      character(:), allocatable :: filename !<title of the file
   end type dcdinfo_t

   logical dcdio
   save

   interface
      module subroutine dcdio_write(dcdinfo,istep,suffix)
         type(dcdinfo_t), intent(inout) :: dcdinfo
         character(*), intent(in) :: suffix
         integer, intent(in) :: istep
      end subroutine
   end interface

   interface
      module subroutine dcdfile_open(dcdinfo,dcdfile)
         type(dcdinfo_t), intent(inout) :: dcdinfo
         character(*), intent(in) :: dcdfile
      end subroutine
   end interface

   interface
      module subroutine dcdfile_read_header(dcdinfo,dowrite,verbose)
         type(dcdinfo_t), intent(inout) :: dcdinfo
         logical, intent(in) :: dowrite
         logical, intent(in), optional :: verbose
      end subroutine
   end interface

   interface
      module subroutine dcdfile_close(dcdinfo)
         type(dcdinfo_t), intent(inout) :: dcdinfo
      end subroutine
   end interface

   interface
      module subroutine dcdfile_read_next(dcdinfo)
         type(dcdinfo_t), intent(inout) :: dcdinfo
      end subroutine
   end interface

   interface
      module subroutine dcdfile_skip_next(dcdinfo,n)
         type(dcdinfo_t), intent(inout) :: dcdinfo
         integer(kind=4), intent(in), optional :: n
      end subroutine
   end interface

end
