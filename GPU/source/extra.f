c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine extra  --  user defined extra potentials  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "extra" calculates any additional user defined potential
c     energy contribution
c
c
#include "tinker_precision.h"
      subroutine extra
      use atoms
      use domdec
      use energi
      use sizes
      use tinheader ,only:ti_p,re_p
      implicit none
      integer i,iglob
c
c
c     zero out the energy due to extra potential terms
c
      ex = 0.0_ti_p
c
c     add any user-defined extra potentials below here
c
c     Note that in Tinker-HP two sets of indexes exist (due to the spatial decomposition
c     used to run in parallel): the local and the global one.
c     The global index is the one defined by the xyz file, the local involves the atoms
c     treated by the local process and the ones belonging to the neighboring ones (closer
c     than half the larger cutoff involving non bonded interactions, see midpoint method
c     for more explanations), it is therefore updated at each time step.
c     The number of local atoms (belonging to the local process) is nloc, the number of
c     local + neighboring atoms is nbloc. 
c     It is possible to switch between the two indexes by using the arrays "loc" and
c     "glob" that are in the 'openmp.i' common bloc.
c     The forces arrays such as dex are in the local index.
c      
c     Many of the global parameters arrays such as the multipoles, the positions, the atom
c     types are set in the global index
c
c     e = ......
c     ex = ex + e
c      do i = 1, nloc
c        iglob = glob(i)
c        ex = ex + ... 
c      end do
c
      return
      end
