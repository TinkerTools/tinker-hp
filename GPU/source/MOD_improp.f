c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module improp  --  improper dihedrals in the current structure  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     kprop    force constant values for improper dihedral angles
c     winkprop    window object corresponding to kprop
c     vprop    ideal improper dihedral angle value in degrees
c     winvprop    window object corresponding to vprop
c     niprop   total number of improper dihedral angles in the system
c     iiprop   numbers of the atoms in each improper dihedral angle
c     winiiprop    window object corresponding to iiprop
c     nbimprop number of improper diehedral before each atom
c     winnbimprop    window object corresponding to nbimprop
c
c
#include "tinker_precision.h"
      module improp
      implicit none
      integer niprop,niproploc
      integer, pointer :: iiprop(:,:),nbimprop(:)
      real(t_p), pointer :: kprop(:),vprop(:)
      integer :: winiiprop,winnbimprop,winkprop,winvprop
      save
      end
