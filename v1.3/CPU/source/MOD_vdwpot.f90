!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #####################################################################
!     ##                                                                 ##
!     ##  module vdwpot  --  specifics of van der Waals functional form  ##
!     ##                                                                 ##
!     #####################################################################
!
!
!
module vdwpot
   implicit none
   integer :: ngauss !<number of Gaussians used in fit to vdw potential
   real*8 :: abuck !<value of "A" constant in Buckingham vdw potential
   real*8 :: bbuck !<value of "B" constant in Buckingham vdw potential
   real*8 :: cbuck !<value of "C" constant in Buckingham vdw potential
   real*8 :: ghal !<value of "gamma" in buffered 14-7 vdw potential
   real*8 :: dhal !<value of "delta" in buffered 14-7 vdw potential
   real*8 :: v2scale !<factor by which 1-2 vdw interactions are scaled
   real*8 :: v3scale !<factor by which 1-3 vdw interactions are scaled
   real*8 :: v4scale !<factor by which 1-4 vdw interactions are scaled
   real*8 :: v5scale !<factor by which 1-5 vdw interactions are scaled
!   real*8 :: igauss(2,maxgauss) !<coefficients of Gaussian fit to vdw potential
   logical :: use_vcorr !<flag to use long range vdw der Waals correction
   character*5 :: vdwindex !<indexing mode (atom type or class) for vdw parameters
   character*5 :: radtyp !<type of parameter (sigma or R-min) for atomic size
   character*8 :: radsiz !<atomic size provided as radius or diameter
   character*8 :: gausstyp !<type of Gaussian fit to van der Waals potential
   character*10 :: radrule !<combining rule for atomic size parameters
   character*10 :: epsrule !<combining rule for vdw well depth parameters
   character*13 :: vdwtyp !<type of van der Waals potential energy function
   save
end
