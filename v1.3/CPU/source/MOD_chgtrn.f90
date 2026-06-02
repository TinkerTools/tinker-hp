!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  module chgtrn  --  charge transfer for current structure  ##
!     ##                                                            ##
!     ################################################################
!
module chgtrn
   implicit none
   integer :: nct !<total number of dispersion sites in the system
   integer :: winchgct !<window associated to chgct array
   integer :: windmpct !<window associated to dmpct array
   real*8, pointer :: chgct(:) !<charge for charge transfer at each multipole site
   real*8, pointer :: dmpct(:) !<charge transfer damping factor at each multipole site
   save
end
