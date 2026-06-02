!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module math  --  mathematical and geometrical constants  ##
!     ##                                                           ##
!     ###############################################################
!
!
module math
   implicit none
   real*8 :: radian !<conversion factor from radians to degrees
   real*8 :: pi !<numerical value of the geometric constant
   real*8 :: sqrtpi !<numerical value of the square root of Pi
   real*8 :: logten !<numerical value of the natural log of ten
   real*8 :: sqrttwo !<numerical value of the square root of two
   real*8 :: twosix !<numerical value of the sixth root of two
   parameter (radian=57.29577951308232088d0)
   parameter (pi=3.141592653589793238d0)
   parameter (sqrtpi=1.772453850905516027d0)
   parameter (logten=2.302585092994045684d0)
   parameter (sqrttwo=1.414213562373095049d0)
   parameter (twosix=1.122462048309372981d0)
   save
end
