!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module kbonds  --  forcefield parameters for bond stretching  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     maxnb   maximum number of bond stretch parameter entries
!     maxnb5  maximum number of 5-membered ring bond stretch entries
!     maxnb4  maximum number of 4-membered ring bond stretch entries
!     maxnb3  maximum number of 3-membered ring bond stretch entries
!     maxnel  maximum number of electronegativity bond corrections
!
!     bcon    force constant parameters for harmonic bond stretch
!     blen    bond length parameters for harmonic bond stretch
!     bcon5   force constant parameters for 5-ring bond stretch
!     blen5   bond length parameters for 5-ring bond stretch
!     bcon4   force constant parameters for 4-ring bond stretch
!     blen4   bond length parameters for 4-ring bond stretch
!     bcon3   force constant parameters for 3-ring bond stretch
!     blen3   bond length parameters for 3-ring bond stretch
!     dlen    electronegativity bond length correction parameters
!     kb      string of atom classes for harmonic bond stretch
!     kb5     string of atom classes for 5-ring bond stretch
!     kb4     string of atom classes for 4-ring bond stretch
!     kb3     string of atom classes for 3-ring bond stretch
!     kel     string of atom classes for electronegativity corrections
!
!
#include "tinker_macro.h"
module kbonds
   implicit none
   integer maxnb,maxnb5,maxnb4
   integer maxnbm,maxnbm4,maxnbflt
   integer maxnb3,maxnel
   parameter (maxnb=2000)
   parameter (maxnb5=500)
   parameter (maxnb4=500)
   parameter (maxnb3=500)
   parameter (maxnel=500)
   parameter (maxnbm=500)
   parameter (maxnbm4=500)
   parameter (maxnbflt=500)
   real(r_p) bcon(maxnb),blen(maxnb)
   real(r_p) bcon5(maxnb5),blen5(maxnb5)
   real(r_p) bcon4(maxnb4),blen4(maxnb4)
   real(r_p) bcon3(maxnb3),blen3(maxnb3)
   real(r_p) dlen(maxnel)
   real(r_p) bmor(3,maxnbm),bmor4(3,maxnbm4),bflat(3,maxnbflt)
   character*8 kb(maxnb),kb5(maxnb5)
   character*8 kb4(maxnb4),kb3(maxnb3)
   character*8 kbm(maxnbm),kbm4(maxnbm4),kbfl(maxnbflt)
   character*12 kel(maxnel)
   save
end
