c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module imptor  --  improper torsions in the current structure  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     itors1   1-fold amplitude and phase for each improper torsion
c     winitors1 window corresponding to itors1
c     itors2   2-fold amplitude and phase for each improper torsion
c     winitors2 window corresponding to itors2
c     itors3   3-fold amplitude and phase for each improper torsion
c     winitors3 window corresponding to itors3
c     nitors   total number of improper torsional angles in the system
c     winbimptors window corresponding to nbimptors
c     iitors   numbers of the atoms in each improper torsional angle
c     winiitors window corresponding to iitors
c
c     nbimptor number of improper torsions before each atom
c
c
#include "tinker_macro.h"
      module imptor
      implicit none
      integer nitors,nitorsloc
      integer, pointer :: iitors(:,:),nbimptor(:)
      real(t_p), pointer :: itors1(:,:),itors2(:,:),itors3(:,:)
      integer :: winiitors,winnbimptor,winitors1,winitors2,winitors3
      save
      end
