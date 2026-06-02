!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module urypot  --  specifics of Urey-Bradley functional form  ##
!     ##                                                                ##
!     ####################################################################
!
!
!     cury       cubic coefficient in Urey-Bradley potential
!     qury       quartic coefficient in Urey-Bradley potential
!     ureyunit   convert Urey-Bradley energy to kcal/mole
!
!
#include "tinker_macro.h"
module urypot
   implicit none
   real(t_p) cury,qury
   real(t_p) ureyunit
   integer  , pointer :: ureytypI(:)
   integer :: winureytypI
   enum,bind(C)
      enumerator UREY_BRAD, UREY_ANGREP, UREY_QUARTIC
   end enum
   save
end
