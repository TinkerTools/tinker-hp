c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  module kvdws  --  forcefield parameters for van der Waals terms  ##
c     ##                                                                   ##
c     #######################################################################
c
c
c     rad      van der Waals radius parameter for each atom type
c     eps      van der Waals well depth parameter for each atom type
c     rad4     van der Waals radius parameter in 1-4 interactions
c     eps4     van der Waals well depth parameter in 1-4 interactions
c     reduct   van der Waals reduction factor for each atom type
c     radv     van der Waals radius parameter for each atom
c     epsv     van der Waals well depth parameter for each atom
c     vadradrule  pairwise vdw rule type
c
c
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
