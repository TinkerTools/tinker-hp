!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #######################################################################
!     ##                                                                   ##
!     ##  module kvdws  --  forcefield parameters for van der Waals terms  ##
!     ##                                                                   ##
!     #######################################################################
!
!
!
module kvdws
   use sizes
   implicit none
   real*8 :: rad(maxtyp) !<van der Waals radius parameter for each atom type
   real*8 :: eps(maxtyp) !<van der Waals well depth parameter for each atom type
   real*8 :: rad4(maxtyp) !<van der Waals radius parameter in 1-4 interactions
   real*8 :: eps4(maxtyp) !<van der Waals well depth parameter in 1-4 interactions
   real*8 :: reduct(maxtyp) !<van der Waals reduction factor for each atom type
   save
end
