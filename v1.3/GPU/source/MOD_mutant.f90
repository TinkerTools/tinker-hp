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
!     vcouple    van der Waals lambda type (0=decouple, 1=annihilate)
!     lambda     generic weighting between initial and final states
!     vlambda    state weighting value for electrostatic potentials
!     elambda    state weighting value for van der Waals potentials
!     tlambda    state weighting value for torsional potential
!     bvlambda   intervall bound for state weighting value vlambda
!     belambda   intervall bound for state weighting value elambda
!     bplambda   value of elambda from which pol is activated
!     flambdabias scalar bias to be applied to flambda (osrw)
!     scexp  softcore vdw parameter for vdw Halgren potential
!     scalpha softcore vdw parameter for vdw Halgren and LJ potential
!     softcore parameters for vdw LJ potential:
!     scvdw = rv*(scalpha*2**(-sck/6)*(1-lambda)**scs+rho**sck)**(1/sck)
!     V_sc = lamda**sct*(V_vdw(scvdw))
!     sck       softcore main exponent
!     sct       softcore external lambda exponent
!     scs       softcore internal lambda exponent
!     nmut       number of atoms mutated from initial to final state
!     imut       atomic sites differing in initial and final state
!     winimut    window object corresponding to imut
!     type0      atom type of each atom in the initial state system
!     wintype0    window object corresponding to type0
!     class0     atom class of each atom in the initial state system
!     winclass0    window object corresponding to class0
!     type1      atom type of each atom in the final state system
!     wintype1    window object corresponding to type1
!     class1     atom class of each atom in the final state system
!     winclass1    window object corresponding to class1
!     mut        true if an atom is to be mutated, false otherwise
!     mutInt        1 if an atom is to be mutated,     0 otherwise
!     winmut    window object corresponding to mut
!     deflambda derivative of direct permanent field wrt elambda
!
!
#include "tinker_macro.h"
module mutant
   implicit none
   integer nmut
   integer vcouple
   integer, pointer :: imut(:),type0(:),class0(:)
   integer :: winimut,wintype0,winclass0
   integer, pointer :: type1(:),class1(:)
   integer :: wintype1,winclass1
   logical, pointer :: mut(:)
   integer(1),pointer:: mutInt(:)
   integer :: winmut,winmutInt
   real(t_p) lambda
   real(t_p) vlambda,elambda,tlambda
   real(t_p) scexp,scalpha
   real(t_p) sck,sct,scs
   real(t_p) bvlambda, belambda
   real(t_p) bplambda
   real(t_p) flambdabias
   real(t_p), allocatable :: deflambda(:,:,:)
end
