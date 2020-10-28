c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine command  --  get any command line arguments  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "command" uses the standard Unix-like iargc/getarg routines
c     to get the number and values of arguments specified on the
c     command line at program runtime
c
c
#include "tinker_precision.h"
      subroutine command
      use argue
      implicit none
      integer i,iargc
      character*1 letter
      character*20 blank
c
c
c     initialize command line arguments as blank strings
c
      narg = 0
      blank = '                    '
      do i = 0, maxarg
         arg(i) = blank//blank//blank
      end do
c
c     get the number of arguments and store each in a string
c
      narg = command_argument_count ()
      if (narg .gt. maxarg)  narg = maxarg
      do i = 0, narg
         call get_command_argument (i,arg(i))
      end do
c
c     mark the command line options as unuseable for input
c
      listarg(0) = .false.
      do i = 1, narg
         listarg(i) = .true.
      end do
      do i = 1, narg
         letter = arg(i)(1:1)
         if (letter .eq. '-') then
            letter = arg(i)(2:2)
            call upcase (letter)
            if (letter.ge.'A' .and. letter.le.'Z') then
               listarg(i) = .false.
               listarg(i+1) = .false.
            end if
         end if
      end do
      return
      end
