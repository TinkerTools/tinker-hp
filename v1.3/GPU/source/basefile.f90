!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine basefile  --  get base prefix from a filename  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "basefile" extracts from an input filename the portion
!     consisting of any directory name and the base filename
!
!
#include "tinker_precision.h"
subroutine basefile (string)
   use ascii
   use files
   implicit none
   integer i,k,trimtext
   character*1 letter
   character*240 string
!
!
!     store the input filename and find its full length
!
   filename = string
   leng = trimtext (string)
!
!     count the number of characters prior to any extension
!
   k = leng
   do i = 1, leng
      letter = string(i:i)
      if (letter .eq. '/')  k = leng
!        if (letter .eq. '\')  k = leng
      if (ichar(letter) .eq. backslash)  k = leng
      if (letter .eq. ']')  k = leng
      if (letter .eq. ':')  k = leng
      if (letter .eq. '~')  k = leng
      if (letter .eq. '.')  k = i - 1
   end do
   leng = min(leng,k)
!
!     find the length of any directory name prefix
!
   k = 0
   do i = leng, 1, -1
      letter = string(i:i)
      if (letter .eq. '/')  k = i
!        if (letter .eq. '\')  k = i
      if (ichar(letter) .eq. backslash)  k = i
      if (letter .eq. ']')  k = i
      if (letter .eq. ':')  k = i
      if (letter .eq. '~')  k = i
!        if (letter .eq. '.')  k = i
      if (k .ne. 0)  goto 10
   end do
10 continue
   ldir = k
!
!     read and store the keywords from the keyfile
!
   call getkey
!
!     get the information level and output style
!
   call control
   return
end
