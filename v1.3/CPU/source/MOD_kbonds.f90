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
module kbonds
   implicit none
   integer :: maxnb !<maximum number of bond stretch parameter entries
   integer :: maxnb5 !<maximum number of 5-membered ring bond stretch entries
   integer :: maxnb4 !<maximum number of 4-membered ring bond stretch entries
   integer :: maxnb3 !<maximum number of 3-membered ring bond stretch entries
   integer :: maxnel !<maximum number of electronegativity bond corrections
   integer :: maxnbm !<maximum number of morse paramter entries
   integer :: maxnbm4 !<maximum number of morse4 paramter entries
   parameter (maxnb=2000)
   parameter (maxnbm=2000)
   parameter (maxnbm4=2000)
   parameter (maxnb5=500)
   parameter (maxnb4=500)
   parameter (maxnb3=500)
   parameter (maxnel=500)
   real*8 :: bcon(maxnb) !<force constant parameters for harmonic bond stretch
   real*8 :: blen(maxnb) !<bond length parameters for harmonic bond stretch
   real*8 :: bmor(3,maxnbm) !<morse parameters for morse potential
   real*8 :: bmor4(3,maxnbm4) !<morse4 parameters for morse4 potential
   real*8 :: bcon5(maxnb5) !<force constant parameters for 5-ring bond stretch
   real*8 :: blen5(maxnb5) !<bond length parameters for 5-ring bond stretch
   real*8 :: bcon4(maxnb4) !<force constant parameters for 4-ring bond stretch
   real*8 :: blen4(maxnb4) !<bond length parameters for 4-ring bond stretch
   real*8 :: bcon3(maxnb3) !<force constant parameters for 3-ring bond stretch
   real*8 :: blen3(maxnb3) !<bond length parameters for 3-ring bond stretch
   real*8 :: dlen(maxnel) !<electronegativity bond length correction parameters
   character*8 :: kb(maxnb) !<string of atom classes for harmonic bond stretch
   character*8 :: kb5(maxnb5) !<string of atom classes for 5-ring bond stretch
   character*8 :: kb4(maxnb4) !<string of atom classes for 4-ring bond stretch
   character*8 :: kb3(maxnb3) !<string of atom classes for 3-ring bond stretch
   character*8 :: kbm(maxnbm) !<string of atom classes for morse
   character*8 :: kbm4(maxnbm4) !<string of atom classes for morse4
   character*12 :: kel(maxnel) !<string of atom classes for electronegativity corrections
   save
end
