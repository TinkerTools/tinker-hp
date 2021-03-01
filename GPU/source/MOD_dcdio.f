c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module dcdmod  --  dcd input/output (I/O) global variables  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     idcd  dcd file unit
c     nframes_pos current number of frames
c     iend_pos position of the end of the file
c     curr_pos current position in the file
c     nframes number of frames in the file
c     istart timestep of last snapshot
c     iend timestep of last snapshot
c     nevery number of steps between outputs
c
c     reader:
c     filesize size of the dcd file to be read
c     framesize size of a frame
c     titlesdcd: titles of the trajectory
c     timestep: timestep of the trajectory
c     natoms: number of atoms in the trajectory
c
c     dcdio: write/read trajectory in the dcd format
c
      module dcdmod
      implicit none
      integer(kind=4) :: idcd,nframes,istart,iend,nevery
      integer(kind=8) :: nframes_pos,iend_pos,curr_pos
      integer(kind=8) :: filesize, framesize
      integer(kind=4) :: natoms
      real(kind=4) :: timestep
      character(len=80), allocatable :: titlesdcd(:)
      logical dcdio
      save
      end
