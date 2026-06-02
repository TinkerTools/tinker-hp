!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###################################################################
!     ##                                                               ##
!     ##  module output  --  control of coordinate output file format  ##
!     ##                                                               ##
!     ###################################################################
!
!
module output
   implicit none
   logical :: archive !<logical flag to save structures in an archive
   logical :: noversion !<logical flag governing use of filename versions
   logical :: overwrite !<logical flag to overwrite intermediate files inplace
   logical :: cyclesave !<logical flag to mark use of numbered cycle files
   character*9 :: coordtype !<selects Cartesian, internal, rigid body or none
   save
end
