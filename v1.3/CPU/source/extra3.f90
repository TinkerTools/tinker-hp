!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ############################################################
!     ##                                                        ##
!     ##  subroutine extra3  --  user defined extra potentials  ##
!     ##                                                        ##
!     ############################################################
!
!
!     "extra3" calculates any additional user defined potential
!     contribution and also partitions the energy among the atoms
!
!
!> @brief 
!> calculates any additional user defined potential
!> energy contribution
!> @param no params
subroutine extra3
   use sizes
   use action
   use analyz
   use atoms
   use domdec
   use energi
   use inform
   use iounit
   implicit none
!
   if (deb_Path) write(iout,*), 'extra3 '
!
!
!     zero out the energy due to extra potential terms
!
   nex = 0
   ex = 0.0d0
   aex = 0.0d0
!
!     add any user-defined extra potentials and partitioning
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
!     "glob" that are in the 'domdec' module.
!     The forces arrays such as dex are in the local index.
!
!     Many of the global parameters arrays such as the multipoles, the positions, the atom
!     types are set in the global index
!
!     e = ......
!      do i = 1, nloc
!        iglob = glob(i)
!        nex = nex + 1
!        ex = ex + ...
!        aex(i) = ...
!      end do
!
   return
end
