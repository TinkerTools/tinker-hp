!
!
!     ###################################################
!     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
!     ##              All Rights Reserved              ##
!     ###################################################
!
!     ############################################################
!     ##                                                        ##
!     ##  function nexttext  --  find next non-blank character  ##
!     ##                                                        ##
!     ############################################################
!
!
!     "nexttext" finds and returns the location of the first
!     non-blank character within an input text string; zero
!     is returned if no such character is found
!
!
!> @brief 
!> "nexttext" finds and returns the location of the first
!> non-blank character within an input text string; zero
!> is returned if no such character is found
!> @param[in] string: command line
function nexttext (string)
   implicit none
   integer i,size
   integer len,nexttext
   character*(*) string
!
!
!     move forward through the string, one character
!     at a time, looking for first non-blank character
!
   nexttext = 0
   size = len(string)
   do i = 1, size
      if (string(i:i) .gt. ' ') then
         nexttext = i
         goto 10
      end if
   end do
10 continue
   return
end
