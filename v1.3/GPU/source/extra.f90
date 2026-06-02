!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###########################################################
!     ##                                                       ##
!     ##  subroutine extra  --  user defined extra potentials  ##
!     ##                                                       ##
!     ###########################################################
!
!
!     "extra" calculates any additional user defined potential
!     energy contribution
!
!
#include "tinker_macro.h"
subroutine extra
   use atoms
   use domdec
   use energi
   use sizes
   use tinheader ,only:ti_p,re_p
   implicit none
   integer i,iglob
!
!
!     zero out the energy due to extra potential terms
!
   ex = 0.0_ti_p
!
!     add any user-defined extra potentials below here
!
!     Note that in Tinker-HP two sets of indexes exist (due to the spatial decomposition
!     used to run in parallel): the local and the global one.
!     The global index is the one defined by the xyz file, the local involves the atoms
!     treated by the local process and the ones belonging to the neighboring ones (closer
!     than half the larger cutoff involving non bonded interactions, see midpoint method
!     for more explanations), it is therefore updated at each time step.
!     The number of local atoms (belonging to the local process) is nloc, the number of
!     local + neighboring atoms is nbloc.
!     It is possible to switch between the two indexes by using the arrays "loc" and
!     "glob" that are in the 'openmp.i' common bloc.
!     The forces arrays such as dex are in the local index.
!
!     Many of the global parameters arrays such as the multipoles, the positions, the atom
!     types are set in the global index
!
!     e = ......
!     ex = ex + e
!      do i = 1, nloc
!        iglob = glob(i)
!        ex = ex + ...
!      end do
!
   return
end
