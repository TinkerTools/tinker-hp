c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module kangs  --  forcefield parameters for bond angle bending  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     maxna    maximum number of harmonic angle bend parameter entries
c     maxna5   maximum number of 5-membered ring angle bend entries
c     maxna4   maximum number of 4-membered ring angle bend entries
c     maxna3   maximum number of 3-membered ring angle bend entries
c     maxnap   maximum number of in-plane angle bend parameter entries
c     maxnaf   maximum number of Fourier angle bend parameter entries
c
c     acon     force constant parameters for harmonic angle bends
c     acon5    force constant parameters for 5-ring angle bends
c     acon4    force constant parameters for 4-ring angle bends
c     acon3    force constant parameters for 3-ring angle bends
c     aconp    force constant parameters for in-plane angle bends
c     aconf    force constant parameters for Fourier angle bends
c     ang      bond angle parameters for harmonic angle bends
c     ang5     bond angle parameters for 5-ring angle bends
c     ang4     bond angle parameters for 4-ring angle bends
c     ang3     bond angle parameters for 3-ring angle bends
c     angp     bond angle parameters for in-plane angle bends
c     angf     phase shift angle and periodicity for Fourier bends
c     ka       string of atom classes for harmonic angle bends
c     ka5      string of atom classes for 5-ring angle bends
c     ka4      string of atom classes for 4-ring angle bends
c     ka3      string of atom classes for 3-ring angle bends
c     kap      string of atom classes for in-plane angle bends
c     kaf      string of atom classes for Fourier angle bends
c
c
#include "tinker_macro.h"
      module kangs
      implicit none
      integer maxna,maxna5,maxna4
      integer maxna3,maxnap,maxnaf
      parameter (maxna=2000)
      parameter (maxna5=500)
      parameter (maxna4=500)
      parameter (maxna3=500)
      parameter (maxnap=2000)
      parameter (maxnaf=500)
      real*8 acon(maxna),acon5(maxna5),acon4(maxna4)
      real*8 acon3(maxna3),aconf(maxnaf)
      real*8 aconp(maxnap)
      real*8 ang(3,maxna),ang5(3,maxna5),ang4(3,maxna4)
      real*8 ang3(3,maxna3),angf(2,maxnaf),angp(2,maxnap)
      character*12 ka(maxna),ka5(maxna5),ka4(maxna4)
      character*12 ka3(maxna3),kaf(maxnaf)
      character*12 kap(maxnap)
      save
      end
