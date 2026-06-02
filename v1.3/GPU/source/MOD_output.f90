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
!     archive    logical flag to save structures in an archive
!     noversion  logical flag governing use of filename versions
!     overwrite  logical flag to overwrite intermediate files inplace
!     cyclesave  logical flag to mark use of numbered cycle files
!     coordtype  selects Cartesian, internal, rigid body or none
!     new_restart logical flag to write restart  in a new file
!     f_mdsave    logical flag to force md traj writing
!
!
module output
   implicit none
   logical archive,noversion
   logical overwrite,cyclesave
   logical new_restart,f_mdsave
   character*9 coordtype
end
