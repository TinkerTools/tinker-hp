c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module boxes  --  parameters for periodic boundary conditions  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     xbox        length of a-axis of periodic box in Angstroms
c     ybox        length of b-axis of periodic box in Angstroms
c     zbox        length of c-axis of periodic box in Angstroms
c     alpha       angle between b- and c-axes of box in degrees
c     beta        angle between a- and c-axes of box in degrees
c     gamma       angle between a- and b-axes of box in degrees
c     xbox2       half of the a-axis length of periodic box
c     ybox2       half of the b-axis length of periodic box
c     zbox2       half of the c-axis length of periodic box
c     box34       three-fourths axis length of truncated octahedron
c     lvec        real space lattice vectors as matrix rows
c     recip       reciprocal lattice vectors as matrix columns
c     volbox      volume in Ang**3 of the periodic box
c     beta_sin    sine of the beta periodic box angle
c     beta_cos    cosine of the beta periodic box angle
c     gamma_sin   sine of the gamma periodic box angle
c     gamma_cos   cosine of the gamma periodic box angle
c     beta_term   term used in generating triclinic box
c     gamma_term  term used in generating triclinic box
c     orthogonal  flag to mark periodic box as orthogonal
c     monoclinic  flag to mark periodic box as monoclinic
c     triclinic   flag to mark periodic box as triclinic
c     octahedron  flag to mark box as truncated octahedron
c     spacegrp    space group symbol for the unitcell type
c
c
#include "tinker_precision.h"
      module boxes
      implicit none
      real(r_p) xbox,ybox,zbox
      real(t_p) alpha,beta,gamma
      real(r_p) xbox2,ybox2,zbox2
      real(t_p) box34
      real(r_p) volbox
      real(r_p),target:: lvec(3,3)
      real(t_p),target:: recip(3,3)
      real(t_p) beta_sin,beta_cos
      real(t_p) gamma_sin,gamma_cos
      real(t_p) beta_term,gamma_term
      logical orthogonal,monoclinic
      logical triclinic,octahedron
      character*10 spacegrp
!$acc declare create(orthogonal,octahedron)
!$acc declare create(lvec,recip,volbox,xbox,ybox,zbox,xbox2,
!$acc& ybox2,zbox2,alpha,beta,gamma,box34)
      save
      end
