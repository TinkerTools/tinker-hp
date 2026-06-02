!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##############################################################
!     ##                                                          ##
!     ##  subroutine getstring  --  extract single quoted string  ##
!     ##                                                          ##
!     ##############################################################
!
!
!     "getstring" searches for a quoted text string within an input
!     character string; the region between the first and second
!     quotes is returned as the "text"; if the actual text is too
!     long, only the first part is returned
!
!     variables and parameters:
!
!     string    input character string to be searched
!     text      the quoted text found in the input string
!     next      input with first position of search string;
!                 output with the position following text
!
!
!> @brief 
!> searches for a quoted text string within an input
!> character string; the region between the first and second
!> quotes is returned as the "text"; if the actual text is too
!> long, only the first part is returned
!> @param[in] string: input character string to be searched
!> @param[in] text: the quoted text found in the input string
!> @param[in] next:  input with first position of search string;
!                 output with the position following text
!
subroutine getstring (string,text,next)
   use ascii
   implicit none
   integer i,j
   integer len,length
   integer size,next
   integer code,extent
   integer first,last
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
!     searching for the quoted text string characters
!
   first = next
   last = 0
   initial = next
   final = next + length - 1
   do i = initial, final
      code = ichar(string(i:i))
      if (code .eq. quote) then
         first = i + 1
         do j = first, final
            code = ichar(string(j:j))
            if (code .eq. quote) then
               last = j - 1
               next = j + 1
               goto 10
            end if
         end do
      end if
   end do
10 continue
!
!     trim the actual word if it is too long to return
!
   extent = last - first + 1
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
   do i = last+1, final
      j = j + 1
      text(j:j) = ' '
   end do
   return
end
