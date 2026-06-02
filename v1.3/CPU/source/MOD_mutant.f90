!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module mutant  --  hybrid atoms for free energy perturbation  ##
!     ##                                                                ##
!     ####################################################################
!
!
module mutant
   implicit none
   integer :: nmut !<number of atoms mutated from initial to final state
   integer :: vcouple !<van der Waals lambda type (0=decouple, 1=annihilate)
   integer :: winimut !<window object corresponding to imut
   integer :: wintype0 !<window object corresponding to type0
   integer :: winclass0 !<window object corresponding to class0
   integer :: wintype1 !<window object corresponding to type1
   integer :: winclass1 !<window object corresponding to class1
   integer :: winmut !<window object corresponding to mut
   integer, pointer :: type1(:) !<atom type of each atom in the final state system
   integer, pointer :: class1(:) !<atom class of each atom in the final state system
   integer, pointer :: imut(:) !<atomic sites differing in initial and final state
   integer, pointer :: type0(:) !<atom type of each atom in the initial state system
   integer, pointer :: class0(:) !<atom class of each atom in the initial state system
   logical, pointer :: mut(:) !<true if an atom is to be mutated, false otherwise
   real*8 :: lambda !<generic weighting between initial and final states
   real*8 :: vlambda !<state weighting value for electrostatic potentials
   real*8 :: elambda !<state weighting value for van der Waals potentials
   real*8 :: tlambda !<state weighting value for torsional potential
   real*8 :: scexp !<softcore vdw parameter for vdw Halgren potential
   real*8 :: scalpha !<softcore vdw parameter for vdw Halgren and LJ potential
   real*8 :: sck !<softcore main exponent
   real*8 :: sct !<softcore external lambda exponent
   real*8 :: scs !<softcore internal lambda exponent
   real*8 :: bvlambda !<intervall bound for state weighting value vlambda
   real*8 :: belambda !<intervall bound for state weighting value elambda
   real*8 :: flambdabias !<scalar bias to be applied to flambda (osrw)
   real*8, allocatable :: deflambda(:,:,:) !<derivative of direct permanent field wrt elambda
   save
end
