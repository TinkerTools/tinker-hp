!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module keys  --  contents of current keyword parameter file  ##
!     ##                                                               ##
!     ###################################################################
!
!
module keys
   use sizes
   implicit none
   integer :: nkey !<number of nonblank lines in the keyword file
   character*240 :: keyline(maxkey) !<contents of each individual keyword file line
   save
end
