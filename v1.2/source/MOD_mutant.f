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
c     scexp
c     scalpha
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
c     winmut    window object corresponding to mut
c
c
      module mutant
      implicit none
      integer nmut
      integer vcouple
      integer, pointer :: imut(:),type0(:),class0(:)
      integer :: winimut,wintype0,winclass0
      integer, pointer :: type1(:),class1(:)
      integer :: wintype1,winclass1
      logical, pointer :: mut(:)
      integer :: winmut
      real*8 lambda
      real*8 vlambda,elambda,tlambda
      real*8 scexp,scalpha
      save
      end
