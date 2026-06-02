!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################################
!     ##                                                                  ##
!     ##  module kangs  --  forcefield parameters for bond angle bending  ##
!     ##                                                                  ##
!     ######################################################################
!
!
!     maxna    !<maximum number of harmonic angle bend parameter entries
!     maxna5   !<maximum number of 5-membered ring angle bend entries
!     maxna4   !<maximum number of 4-membered ring angle bend entries
!     maxna3   !<maximum number of 3-membered ring angle bend entries
!     maxnap   !<maximum number of in-plane angle bend parameter entries
!     maxnaf   !<maximum number of Fourier angle bend parameter entries
!
!     acon     !<force constant parameters for harmonic angle bends
!     acon5    !<force constant parameters for 5-ring angle bends
!     acon4    !<force constant parameters for 4-ring angle bends
!     acon3    !<force constant parameters for 3-ring angle bends
!     aconp    !<force constant parameters for in-plane angle bends
!     aconf    !<force constant parameters for Fourier angle bends
!     ang      !<bond angle parameters for harmonic angle bends
!     ang5     !<bond angle parameters for 5-ring angle bends
!     ang4     !<bond angle parameters for 4-ring angle bends
!     ang3     !<bond angle parameters for 3-ring angle bends
!     angp     !<bond angle parameters for in-plane angle bends
!     angf     !<phase shift angle and periodicity for Fourier bends
!     ka       !<string of atom classes for harmonic angle bends
!     ka5      !<string of atom classes for 5-ring angle bends
!     ka4      !<string of atom classes for 4-ring angle bends
!     ka3      !<string of atom classes for 3-ring angle bends
!     kap      !<string of atom classes for in-plane angle bends
!     kaf      !<string of atom classes for Fourier angle bends
!
!
module kangs
   implicit none
   integer :: maxna !<maximum number of harmonic angle bend parameter entries
   integer :: maxna5 !<maximum number of 5-membered ring angle bend entries
   integer :: maxna4 !<maximum number of 4-membered ring angle bend entries
   integer :: maxna3 !<maximum number of 3-membered ring angle bend entries
   integer :: maxnap !<maximum number of in-plane angle bend parameter entries
   integer :: maxnaf !<maximum number of Fourier angle bend parameter entries
   integer :: maxnaps !<maximum number of in-plane angle bend parameter entries
   parameter (maxna=2000)
   parameter (maxna5=500)
   parameter (maxna4=500)
   parameter (maxna3=500)
   parameter (maxnap=2000)
   parameter (maxnaps=2000)
   parameter (maxnaf=500)
   real*8 :: acon(maxna) !<force constant parameters for harmonic angle bends
   real*8 :: acon5(maxna5) !<force constant parameters for 5-ring angle bends
   real*8 :: acon4(maxna4) !<force constant parameters for 4-ring angle bends
   real*8 :: acon3(maxna3) !<force constant parameters for 3-ring angle bends
   real*8 :: aconf(maxnaf) !<force constant parameters for Fourier angle bends
   real*8 :: aconp(maxnap) !<force constant parameters for in-plane angle bends
   real*8 :: ang(3,maxna) !<bond angle parameters for harmonic angle bends
   real*8 :: ang5(3,maxna5) !<bond angle parameters for 5-ring angle bends
   real*8 :: ang4(3,maxna4) !<bond angle parameters for 4-ring angle bends
   real*8 :: ang3(3,maxna3) !<bond angle parameters for 3-ring angle bends
   real*8 :: angf(2,maxnaf) !<phase shift angle and periodicity for Fourier bends
   real*8 :: angp(2,maxnap) !<bond angle parameters for in-plane angle bends
   real*8 :: angps(3,maxnaps) !<bond angle parameters for in-plane angle bends
   character*12 :: ka(maxna) !<string of atom classes for harmonic angle bends
   character*12 :: ka5(maxna5) !<string of atom classes for 5-ring angle bends
   character*12 :: ka4(maxna4) !<string of atom classes for 4-ring angle bends
   character*12 :: ka3(maxna3) !<string of atom classes for 3-ring angle bends
   character*12 :: kaf(maxnaf) !<string of atom classes for Fourier angle bends
   character*12 :: kap(maxnap) !<string of atom classes for in-plane angle bends
   character*12 :: kaps(maxnaps) !<string of atom classes for in-plane angle bends
   save
end
