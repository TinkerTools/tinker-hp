c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module atoms  --  number, position and type of current atoms  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     x       current x-coordinate for each atom in the system
c     y       current y-coordinate for each atom in the system
c     z       current z-coordinate for each atom in the system
c     xold    last x-coordinate for each atom in the system
c     yold    last y-coordinate for each atom in the system
c     zold    last z-coordinate for each atom in the system
c     n       total number of atoms in the current system
c     nloop   First multiple of 16 after n
c     type    atom type number for each atom in the system
c
c
      module atoms
      implicit none
      integer n,nloop
      integer, pointer :: type(:)
      real*8, allocatable ::  x(:),y(:),z(:)
      real*8, allocatable ::  xold(:),yold(:),zold(:)
      save
      end
