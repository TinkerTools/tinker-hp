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
!     idcd  dcd file unit
!     nframes_pos current number of frames
!     iend_pos position of the end of the file
!     curr_pos current position in the file
!     nframes number of frames in the file
!     istart timestep of last snapshot
!     iend timestep of last snapshot
!     nevery number of steps between outputs
!
!     reader:
!     filesize size of the dcd file to be read
!     framesize size of a frame
!     titlesdcd: titles of the trajectory
!     timestep: timestep of the trajectory
!     natoms: number of atoms in the trajectory
!
!     dcdio: write/read trajectory in the dcd format
!
#include "tinker_precision.h"
module dcdmod
   implicit none

   type dcdinfo_t
      integer(kind=4) :: idcd,nframes,istart,iend,nevery
      integer(kind=8) :: nframes_pos,iend_pos,curr_pos
      integer(kind=8) :: filesize, framesize
      integer(kind=4) :: natoms
      real(kind=4) :: timestep
      character(len=80), allocatable :: titlesdcd(:)
      character(:), allocatable :: filename
   end type dcdinfo_t

   logical dcdio
   save

   interface
      module subroutine dcdio_write(dcdinfo,istep,dt,suffix)
         type(dcdinfo_t), intent(inout) :: dcdinfo
         character(*), intent(in) :: suffix
         integer, intent(in) :: istep
         real(r_p), intent(in) :: dt
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
