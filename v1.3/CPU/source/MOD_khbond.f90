!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ## module  khbond  --  forcefield parameters for H-bonding terms  ##
!     ##                                                                ##
!     ####################################################################
!
module khbond
   implicit none
   integer :: maxnhb !<maximum number of hydrogen bonding pair entries
   parameter (maxnhb=500)
   real*8 :: radhb(maxnhb) !<radius parameter for hydrogen bonding pairs
   real*8 :: epshb(maxnhb) !<well depth parameter for hydrogen bonding pairs
   character*8 :: khb(maxnhb) !<string of atom types for hydrogen bonding pairs
   save
end
