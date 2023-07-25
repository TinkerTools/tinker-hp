c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module kbonds  --  forcefield parameters for bond stretching  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     maxnb   maximum number of bond stretch parameter entries
c     maxnb5  maximum number of 5-membered ring bond stretch entries
c     maxnb4  maximum number of 4-membered ring bond stretch entries
c     maxnb3  maximum number of 3-membered ring bond stretch entries
c     maxnel  maximum number of electronegativity bond corrections
c
c     bcon    force constant parameters for harmonic bond stretch
c     blen    bond length parameters for harmonic bond stretch
c     bcon5   force constant parameters for 5-ring bond stretch
c     blen5   bond length parameters for 5-ring bond stretch
c     bcon4   force constant parameters for 4-ring bond stretch
c     blen4   bond length parameters for 4-ring bond stretch
c     bcon3   force constant parameters for 3-ring bond stretch
c     blen3   bond length parameters for 3-ring bond stretch
c     dlen    electronegativity bond length correction parameters
c     kb      string of atom classes for harmonic bond stretch
c     kb5     string of atom classes for 5-ring bond stretch
c     kb4     string of atom classes for 4-ring bond stretch
c     kb3     string of atom classes for 3-ring bond stretch
c     kel     string of atom classes for electronegativity corrections
c
c
#include "tinker_macro.h"
      module kbonds
      implicit none
      integer maxnb,maxnb5,maxnb4
      integer maxnbm,maxnbm4
      integer maxnb3,maxnel
      parameter (maxnb=2000)
      parameter (maxnb5=500)
      parameter (maxnb4=500)
      parameter (maxnb3=500)
      parameter (maxnel=500)
      parameter (maxnbm=500)
      parameter (maxnbm4=500)
      real(r_p) bcon(maxnb),blen(maxnb)
      real(r_p) bcon5(maxnb5),blen5(maxnb5)
      real(r_p) bcon4(maxnb4),blen4(maxnb4)
      real(r_p) bcon3(maxnb3),blen3(maxnb3)
      real(r_p) dlen(maxnel)
      real(r_p) bmor(3,maxnbm),bmor4(3,maxnbm4)
      character*8 kb(maxnb),kb5(maxnb5)
      character*8 kb4(maxnb4),kb3(maxnb3)
      character*8 kbm(maxnbm),kbm4(maxnbm4)
      character*12 kel(maxnel)
      save
      end
