!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module iounit  --  Fortran input/output (I/O) unit numbers  ##
!     ##                                                              ##
!     ##################################################################
!
!
!
module iounit
   implicit none
   integer :: iout !<Fortran I/O unit for main output (default=6)
   integer :: input !<Fortran I/O unit for main input (default=5)
   save
end
