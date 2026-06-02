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
module argue
   implicit none
   integer :: maxarg!<maximum number of command line arguments
   parameter (maxarg=20)
   integer :: narg!<number of command line arguments to the program
   logical :: listarg(0:maxarg)!<flag to mark available command line arguments
   character*240 :: arg(0:maxarg)!<strings containing the command line arguments
   save
end
