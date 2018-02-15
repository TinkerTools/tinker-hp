c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  argue.i  --  command line arguments at program startup  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxarg    maximum number of command line arguments
c
c     narg      number of command line arguments to the program
c     listarg   flag to mark available command line arguments
c     arg       strings containing the command line arguments
c
c
      integer maxarg
      parameter (maxarg=20)
      integer narg
      logical listarg
      character*120 arg
      common /argue/ narg,listarg(0:maxarg),arg(0:maxarg)
