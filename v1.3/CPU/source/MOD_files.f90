!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ####################################################################
!     ##                                                                ##
!     ##  module files  --  name and number of current structure files  ##
!     ##                                                                ##
!     ####################################################################
!
!
!
module files
   implicit none
   integer :: nprior !<number of previously existing cycle files
   integer :: ldir !<length in characters of the directory name
   integer :: leng !<length in characters of the base filename
   character*240 :: filename !<base filename used by default for all files
   character*240 :: outfile !<output filename used for intermediate results
   logical :: keys_already_read=.FALSE. !<info about parsing of key file
   save
end
