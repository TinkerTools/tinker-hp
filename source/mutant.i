c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  mutant.i  --  hybrid atoms for free energy perturbation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     lambda     generic weighting between initial and final states
c     vlambda    state weighting value for electrostatic potentials
c     elambda    state weighting value for van der Waals potentials
c     scexp
c     scalpha
c     nmut       number of atoms mutated from initial to final state
c     imut       atomic sites differing in initial and final state
c     type0      atom type of each atom in the initial state system
c     class0     atom class of each atom in the initial state system
c     type1      atom type of each atom in the final state system
c     class1     atom class of each atom in the final state system
c     mut        true if an atom is to be mutated, false otherwise
c
c
      integer nmut
      integer, pointer :: imut(:),type0(:),class0(:)
      integer, pointer :: type1(:),class1(:)
      logical, pointer :: mut(:)
      real*8 lambda
      real*8 vlambda,elambda
      real*8 scexp,scalpha
      common /mutant/ lambda,vlambda,elambda,scexp,scalpha,
     &                nmut,imut,type0,class0,
     &                type1,class1,mut
