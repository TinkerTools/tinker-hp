!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module  sizes  --  parameter values to set array dimensions  ##
!     ##                                                               ##
!     ###################################################################
!
!
!     "sizes.f" sets values for critical array dimensions used
!     throughout the software
!
!     parameter:      maximum allowed number of:
!
!     maxvalue          atoms directly bonded to an atom
!     maxgrp          user-defined groups of atoms
!     maxtyp          force field atom type definitions
!     maxclass        force field atom class definitions
!     maxprm          lines in the parameter file
!     maxkey          lines in the keyword file
!     maxvlst         neighbors in van der Waals pair list
!     maxelst         neighbors in electrostatics pair list
!     maxfft          grid points in each FFT dimension
!     maxring         3-, 4-, or 5-membered rings
!     maxbio          biopolymer atom definitions
!     maxres          residues in the macromolecule
!     maxele          elements in periodic table
!     maxamino        amino acid residue types
!     maxnuc          nucleic acid residue types
!     maxtors         torsional angles in molecular system
!     maxbitor        bitorsions in molecular system
!
!
#include "tinker_macro.h"
module sizes
   implicit none
   integer maxvalue,maxgrp
   integer maxtyp,maxclass
   integer maxprm,maxkey
   integer maxopt
   integer maxvlst,maxelst
   integer maxfft
   integer maxcell,maxref
   integer maxring,maxbio,maxres
   integer maxele,maxamino,maxnuc
   integer maxbnd,maxang,maxtors
   integer maxbitor,maxlp
   parameter (maxvalue =     8)
   parameter (maxgrp   =  1000)
   parameter (maxref   =    10)
   parameter (maxtyp   =  5000)
   parameter (maxclass =  2000)
   parameter (maxprm   = 25000)
   parameter (maxkey   =  5000)
#if TINKERHP_REL_BUILD
   parameter (maxvlst  =  2500)
   parameter (maxelst  =  1250)
#else
   parameter (maxvlst  =  1500)
   parameter (maxelst  =   700)
#endif
   parameter (maxfft   =   864)
   parameter (maxring  = 10000)
   parameter (maxbio   = 10000)
   parameter (maxres   = 10000)
   parameter (maxele   =   112)
   parameter (maxamino =    38)
   parameter (maxnuc   =    12)
   integer tinkerdebug

contains
   function size_i8_to_i( i8 ) result(i)
      integer(kind=8) :: i8
      integer :: i
      if( i8 >= huge(i) ) then
         print *, "size_i8_to_i : overflow ", i8, ">=", huge(i)
         stop
      endif
      i = int(i8,kind(i))
   end function

end
