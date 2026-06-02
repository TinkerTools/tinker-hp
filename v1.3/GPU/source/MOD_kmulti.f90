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
!     maxnmp   maximum number of atomic multipole parameter entries
!
!     multip   atomic monopole, dipole and quadrupole values
!     mpaxis   type of local axis definition for atomic multipoles
!     kmp      string of atom types for atomic multipoles
!
!
#include "tinker_macro.h"
module kmulti
   use sizes
   implicit none
   integer maxnmp
   parameter (maxnmp=2000)
   real(t_p) multip(13,maxnmp),sibfacp(3,maxtyp)
   character*8 mpaxis(maxnmp)
   character*16 kmp(maxnmp)
   save
end
