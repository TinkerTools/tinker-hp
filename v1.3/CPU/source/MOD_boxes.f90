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
module boxes
   implicit none
   real*8 :: xbox !<length of a-axis of periodic box in Angstroms
   real*8 :: ybox !<length of b-axis of periodic box in Angstroms
   real*8 :: zbox !<length of c-axis of periodic box in Angstroms
   real*8 :: alpha !<angle between b- and c-axes of box in degrees
   real*8 :: beta !<angle between a- and c-axes of box in degrees
   real*8 :: gamma !<angle between a- and b-axes of box in degrees
   real*8 :: xbox2 !<half of the a-axis length of periodic box
   real*8 :: ybox2 !<half of the b-axis length of periodic box
   real*8 :: zbox2 !<half of the c-axis length of periodic box
   real*8 :: box34 !<three-fourths axis length of truncated octahedron
   real*8 :: volbox !<volume in Ang**3 of the periodic box
   real*8 :: lvec(3,3) !<real space lattice vectors as matrix rows
   real*8 :: recip(3,3) !<reciprocal lattice vectors as matrix columns
   real*8 :: beta_sin !<sine of the beta periodic box angle
   real*8 :: beta_cos !<cosine of the beta periodic box angle
   real*8 :: gamma_sin !<sine of the gamma periodic box angle
   real*8 :: gamma_cos !<cosine of the gamma periodic box angle
   real*8 :: beta_term !<term used in generating triclinic box
   real*8 :: gamma_term !<term used in generating triclinic box
   real*8 :: ftc(10,10) !<fractional to cartesian transformation matrix
   real*8 :: ctfmat(3,3) !<3*3 cartesian to fractional transformation matrix
   real*8 :: ftcmat(3,3) !<3*3 fractional to cartesian transformation matrix
   logical :: orthogonal !<flag to mark periodic box as orthogonal
   logical :: monoclinic !<flag to mark periodic box as monoclinic
   logical :: triclinic !<flag to mark periodic box as triclinic
   logical :: octahedron !<flag to mark box as truncated octahedron
   character*10 :: spacegrp !<space group symbol for the unitcell type
   save
end
