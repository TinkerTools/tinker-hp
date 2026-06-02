!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module imptor  --  improper torsions in the current structure  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     itors1   1-fold amplitude and phase for each improper torsion
!     winitors1 window corresponding to itors1
!     itors2   2-fold amplitude and phase for each improper torsion
!     winitors2 window corresponding to itors2
!     itors3   3-fold amplitude and phase for each improper torsion
!     winitors3 window corresponding to itors3
!     nitors   total number of improper torsional angles in the system
!     winbimptors window corresponding to nbimptors
!     iitors   numbers of the atoms in each improper torsional angle
!     winiitors window corresponding to iitors
!
!     nbimptor number of improper torsions before each atom
!
!
#include "tinker_macro.h"
module imptor
   implicit none
   integer nitors,nitorsloc
   integer, pointer :: iitors(:,:),nbimptor(:)
   real(t_p), pointer :: itors1(:,:),itors2(:,:),itors3(:,:)
   integer :: winiitors,winnbimptor,winitors1,winitors2,winitors3
   save
end
