c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  module math  --  mathematical and geometrical constants  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     radian   conversion factor from radians to degrees
c     pi       numerical value of the geometric constant
c     sqrtpi   numerical value of the square root of Pi
c     logten   numerical value of the natural log of ten
c     sqrttwo  numerical value of the square root of two
c     twosix   numerical value of the sixth root of two
c
c
#include "tinker_precision.h"
      module math
      implicit none
      integer,private,parameter :: math_prec=t_p !Get precision from tinker_precision.h
      real(t_p) radian,pi
      real(t_p) sqrtpi,logten
      real(t_p) sqrttwo,twosix
      parameter (radian =57.29577951308232088_math_prec
     &          ,pi     =3.141592653589793238_math_prec
     &          ,sqrtpi =1.772453850905516027_math_prec
     &          ,logten =2.302585092994045684_math_prec
     &          ,sqrttwo=1.414213562373095049_math_prec
     &          ,twosix =1.122462048309372981_math_prec)
      ! Do not remove precision on those parameters
      ! Differences apprears on the energy
      ! Probably due to ti_p appendix on real variable along the code
      end
