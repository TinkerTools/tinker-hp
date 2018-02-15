c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  opbend.i  --  out-of-plane bends in the current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     opbk      force constant values for out-of-plane bending
c     nopbend   total number of out-of-plane bends in the system
c     iopb      bond angle numbers used in out-of-plane bending
c
c
      integer nopbend,nopbendloc
      integer, pointer :: iopb(:), nbopbend(:)
      real*8, pointer ::  opbk(:)
      common /opbend/ opbk,nopbend,iopb,nopbendloc,nbopbend
