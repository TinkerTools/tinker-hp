!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ############################################################
!     ##                                                        ##
!     ##  subroutine extra1  --  user defined extra potentials  ##
!     ##                                                        ##
!     ############################################################
!
!
!     "extra1" calculates any additional user defined potential
!     energy contribution and its first derivatives
!
!
#include "tinker_macro.h"
subroutine extra1
   use atoms
   use deriv
   use domdec
   use energi
   use tinheader ,only:ti_p,re_p
   use sizes
   implicit none
   integer i,iglob
!
!
!     zero out the extra energy term and first derivatives
!
   ex = 0.0_ti_p
   dex = 0.0_ti_p
!
!     add any user-defined extra potentials and derivatives;
!     also increment intermolecular energy and virial as needed
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
!     ex = ex + e
!
!     loop over the local atom sites
!
!       do i = 1, nloc
!          iglob = glob(i)
!          xi = x(iglob)
!          yi = y(iglob)
!          zi = z(iglob)
!          dex(1,i) = ......
!          dex(2,i) = ......
!          dex(3,i) = ......
!       end do
!
!!$acc data present(ex,dex)
!!$acc update device(ex)
!!$acc end data
   return
end
