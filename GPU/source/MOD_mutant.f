c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module mutant  --  hybrid atoms for free energy perturbation  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     vcouple    van der Waals lambda type (0=decouple, 1=annihilate)
c     lambda     generic weighting between initial and final states
c     vlambda    state weighting value for electrostatic potentials
c     elambda    state weighting value for van der Waals potentials
c     tlambda    state weighting value for torsional potential
c     bvlambda   intervall bound for state weighting value vlambda
c     belambda   intervall bound for state weighting value elambda
c     bplambda   value of elambda from which pol is activated
c     flambdabias scalar bias to be applied to flambda (osrw)     
c     scexp  softcore vdw parameter for vdw Halgren potential
c     scalpha softcore vdw parameter for vdw Halgren and LJ potential
c     softcore parameters for vdw LJ potential:
c     scvdw = rv*(scalpha*2**(-sck/6)*(1-lambda)**scs+rho**sck)**(1/sck)
c     V_sc = lamda**sct*(V_vdw(scvdw))
c     sck       softcore main exponent
c     sct       softcore external lambda exponent
c     scs       softcore internal lambda exponent
c     nmut       number of atoms mutated from initial to final state
c     imut       atomic sites differing in initial and final state
c     winimut    window object corresponding to imut
c     type0      atom type of each atom in the initial state system
c     wintype0    window object corresponding to type0
c     class0     atom class of each atom in the initial state system
c     winclass0    window object corresponding to class0
c     type1      atom type of each atom in the final state system
c     wintype1    window object corresponding to type1
c     class1     atom class of each atom in the final state system
c     winclass1    window object corresponding to class1
c     mut        true if an atom is to be mutated, false otherwise
c     mutInt        1 if an atom is to be mutated,     0 otherwise
c     winmut    window object corresponding to mut
c
c
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
      end
