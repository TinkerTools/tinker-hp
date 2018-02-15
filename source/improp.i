c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  improp.i  --  improper dihedrals in the current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     kprop    force constant values for improper dihedral angles
c     vprop    ideal improper dihedral angle value in degrees
c     niprop   total number of improper dihedral angles in the system
c     iiprop   numbers of the atoms in each improper dihedral angle
c
c
      integer niprop,niproploc
      integer, pointer :: iiprop(:,:),nbimprop(:)
      real*8, pointer :: kprop(:),vprop(:)
      common /improp/ kprop,vprop,niprop,
     &                iiprop,nbimprop,niproploc
