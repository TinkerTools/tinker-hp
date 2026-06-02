!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module bndpot  --  specifics of bond stretch functional forms  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     cbnd      cubic coefficient in bond stretch potential
!     qbnd      quartic coefficient in bond stretch potential
!     bndunit   convert bond stretch energy to kcal/mole
!     bndtyp    type of bond stretch potential energy function
!
!
#include "tinker_macro.h"
module bndpot
   implicit none
   real(t_p) cbnd,qbnd
   real(t_p) bndunit
   character*8 bndtyp
   integer bndtyp_i
   integer , pointer :: bndtypI(:)
   integer :: winbndtypI
   real(t_p) :: flatbottom_delta = 0.3
   enum,bind(C)
      enumerator BND_HARMONIC
      enumerator BND_MORSE
      enumerator BND_MORSE4
      enumerator BND_NO_TYPE
      enumerator BND_FLATBOTTOM
   end enum
   save
end
