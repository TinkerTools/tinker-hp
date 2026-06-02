!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module argue  --  command line arguments at program startup ##
!     ##                                                              ##
!     ##################################################################
!
!
!     maxarg    maximum number of command line arguments
!
!     narg      number of command line arguments to the program
!     listarg   flag to mark available command line arguments
!     arg       strings containing the command line arguments
!
!
#include "tinker_macro.h"
module argue
   implicit none
   integer maxarg
   parameter (maxarg=20)
   integer narg
   logical listarg(0:maxarg)
   character*240 arg(0:maxarg)
end
