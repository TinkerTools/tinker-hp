!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine suffix  --  test for default file extension  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "suffix" checks a filename for the presence of an extension,
!     and appends an extension and version if none is found
!
!
#include "tinker_macro.h"
subroutine suffix (filename,extension,status)
   use ascii
   implicit none
   integer i,leng,lext
   integer last,trimtext
   logical exist
   character*1 letter
   character*3 status
   character*(*) filename
   character*(*) extension
!
!
!     get the length of the current filename
!
   leng = trimtext (filename)
   lext = trimtext (extension)
!
!     check for an extension on the current filename
!
   last = leng
   do i = 1, leng
      letter = filename(i:i)
      if (letter .eq. '/')  last = leng
!        if (letter .eq. '\')  last = leng
      if (ichar(letter) .eq. backslash)  last = leng
      if (letter .eq. ']')  last = leng
      if (letter .eq. ':')  last = leng
      if (letter .eq. '~')  last = leng
      if (letter .eq. '.')  last = i - 1
   end do
!
!     append an extension or version as appropriate
!
   if (last .eq. leng) then
      exist = .false.
      if (leng .ne. 0) then
         inquire (file=filename(1:leng),exist=exist)
      end if
      if (.not. exist) then
         filename = filename(1:leng)//'.'//extension(1:lext)
         call version (filename,status)
      end if
   else if (status .eq. 'new') then
      call version (filename,status)
   end if
   return
end
