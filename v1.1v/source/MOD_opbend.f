c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module opbend  --  out-of-plane bends in the current structure  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     nopbend   total number of out-of-plane bends in the system
c     nopbendloc   local number of out-of-plane bends in the system
c     iopb      bond angle numbers used in out-of-plane bending
c     npopbend  number of angle used in out-of-plane bending before each atom
c     opbk      force constant values for out-of-plane bending
c
c
      module opbend
      implicit none
      integer nopbend,nopbendloc
      integer, pointer :: iopb(:), nbopbend(:)
      real*8, pointer ::  opbk(:)
      save
      end
