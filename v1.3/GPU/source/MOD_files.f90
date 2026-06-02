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
!     nprior     number of previously existing cycle files
!     ldir       length in characters of the directory name
!     leng       length in characters of the base filename
!     filename   base filename used by default for all files
!     outfile    output filename used for intermediate results
!
!
module files
   implicit none
   integer nprior,ldir,leng
   character*240,target:: filename,outfile
   logical :: keys_already_read=.FALSE.

end
