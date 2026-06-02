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
!     rad      van der Waals radius parameter for each atom type
!     eps      van der Waals well depth parameter for each atom type
!     rad4     van der Waals radius parameter in 1-4 interactions
!     eps4     van der Waals well depth parameter in 1-4 interactions
!     reduct   van der Waals reduction factor for each atom type
!     radv     van der Waals radius parameter for each atom
!     epsv     van der Waals well depth parameter for each atom
!     vadradrule  pairwise vdw rule type
!
!
#include "tinker_macro.h"
module kvdws
   use sizes
   implicit none
   real(t_p) rad(maxtyp),eps(maxtyp)
   real(t_p) rad4(maxtyp),eps4(maxtyp)
   real(t_p) reduct(maxtyp)
   real(t_p),pointer:: radv(:),epsv(:)
   save
end
