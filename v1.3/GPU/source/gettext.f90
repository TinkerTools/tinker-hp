!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##########################################################
!     ##                                                      ##
!     ##  subroutine gettext  --  extract text from a string  ##
!     ##                                                      ##
!     ##########################################################
!
!
!     "gettext" searches an input string for the first string of
!     non-blank characters; the region from a non-blank character
!     to the first space or tab is returned as "text"; if the
!     actual text is too long, only the first part is returned
!
!     variables and parameters:
!
!     string    input character string to be searched
!     text      output with the first text string found
!     next      input with first position of search string;
!                 output with the position following text
!
!
#include "tinker_precision.h"
subroutine gettext (string,text,next)
   use ascii
   implicit none
   integer i,j
   integer len,length
   integer size,next
   integer first,last
   integer code,extent
   integer initial,final
   character*(*) string
   character*(*) text
!
!
!     get the length of input string and output text
!
   length = len(string(next:))
   size = len(text)
!
!     move through the string one character at a time,
!     searching for the first non-blank character
!
   first = next
   last = 0
   initial = next
   final = next + length - 1
   do i = initial, final
      code = ichar(string(i:i))
      if (code.ne.space .and. code.ne.tab) then
         first = i
         do j = i+1, final
            code = ichar(string(j:j))
            if (code.eq.space .or. code.eq.tab) then
               last = j - 1
               next = j
               goto 10
            end if
         end do
         last = final
         next = last + 1
      end if
   end do
10 continue
!
!     trim the actual text if it is too long to return
!
   extent = next - first
   final = first + size - 1
   if (extent .gt. size)  last = final
!
!     transfer the text into the return string
!
   j = 0
   do i = first, last
      j = j + 1
      text(j:j) = string(i:i)
   end do
   do i = next, final
      j = j + 1
      text(j:j) = ' '
   end do
   return
end
