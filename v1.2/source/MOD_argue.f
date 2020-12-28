c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  module argue  --  command line arguments at program startup ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxarg    maximum number of command line arguments
c
c     narg      number of command line arguments to the program
c     listarg   flag to mark available command line arguments
c     arg       strings containing the command line arguments
c
c
      module argue
      implicit none
      integer maxarg
      parameter (maxarg=20)
      integer narg
      logical listarg(0:maxarg)
      character*240 arg(0:maxarg)
      save
      end
