!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module boxes  --  parameters for periodic boundary conditions  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!     xbox        length of a-axis of periodic box in Angstroms
!     ybox        length of b-axis of periodic box in Angstroms
!     zbox        length of c-axis of periodic box in Angstroms
!     alpha       angle between b- and c-axes of box in degrees
!     beta        angle between a- and c-axes of box in degrees
!     gamma       angle between a- and b-axes of box in degrees
!     xbox2       half of the a-axis length of periodic box
!     ybox2       half of the b-axis length of periodic box
!     zbox2       half of the c-axis length of periodic box
!     box34       three-fourths axis length of truncated octahedron
!     lvec        real space lattice vectors as matrix rows
!     recip       reciprocal lattice vectors as matrix columns
!     volbox      volume in Ang**3 of the periodic box
!     beta_sin    sine of the beta periodic box angle
!     beta_cos    cosine of the beta periodic box angle
!     gamma_sin   sine of the gamma periodic box angle
!     gamma_cos   cosine of the gamma periodic box angle
!     beta_term   term used in generating triclinic box
!     gamma_term  term used in generating triclinic box
!     orthogonal  flag to mark periodic box as orthogonal
!     monoclinic  flag to mark periodic box as monoclinic
!     triclinic   flag to mark periodic box as triclinic
!     octahedron  flag to mark box as truncated octahedron
!     spacegrp    space group symbol for the unitcell type
!     ftc         fractional to cartesian transformation matrix
!     ftcmat      3*3 fractional to cartesian transformation matrix
!     ctfmat      3*3 cartesian to fractional transformation matrix
!
!
#include "tinker_macro.h"
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
   real(t_p) ctfmat(3,3),ftcmat(3,3)
   logical orthogonal,monoclinic,triclinic,octahedron
   character*10 spacegrp

!$acc declare create(orthogonal,octahedron,monoclinic,triclinic)
!$acc declare create(lvec,recip,ctfmat,ftcmat,volbox,xbox,ybox,zbox,xbox2,ybox2,zbox2 &
!$acc               ,beta_sin,beta_cos,beta_term,gamma_sin,gamma_cos,gamma_term &
!$acc               ,alpha,beta,gamma,box34)

end
