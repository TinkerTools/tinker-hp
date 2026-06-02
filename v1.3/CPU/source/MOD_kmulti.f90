!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kmulti  --  forcefield parameters for atomic multipoles  ##
!     ##                                                                  ##
!     ######################################################################
!
!
module kmulti
   use sizes
   implicit none
   integer :: maxnmp !<maximum number of atomic multipole parameter entries
   parameter (maxnmp=2000)
   real*8 :: multip(13,maxnmp) !<atomic monopole, dipole and quadrupole values
   real*8 :: sibfacp(3,maxtyp) !<sibfa charge penetration parameters
   character*8 :: mpaxis(maxnmp) !<type of local axis definition for atomic multipoles
   character*16 :: kmp(maxnmp) !<string of atom types for atomic multipoles
   save
end
