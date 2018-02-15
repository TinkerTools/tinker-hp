c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  atoms.i  --  number, position and type of current atoms  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     x       current x-coordinate for each atom in the system
c     y       current y-coordinate for each atom in the system
c     z       current z-coordinate for each atom in the system
c     n       total number of atoms in the current system
c     type    atom type number for each atom in the system
c
c
      integer n
      integer, pointer :: type(:)
      real*8, pointer ::  x(:),y(:),z(:)
      real*8, pointer ::  xold(:),yold(:),zold(:)
      common /atoms/ x,y,z,xold,yold,zold,n,type
