c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module strbnd  --  stretch-bends in the current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     sbk       force constants for stretch-bend terms
c     winsbk    window object corresponding to sbk
c     nstrbnd   total number of stretch-bend interactions
c     nstrbndloc   local number of stretch-bend interactions
c     nbstrbnd   number of stretch-bend interactions before each atom
c     winnbstrbnd    window object corresponding to nbstrbnd
c     isb       angle and bond numbers used in stretch-bend
c     winisb    window object corresponding to isb
c
c
      module strbnd
      implicit none
      integer nstrbnd,nstrbndloc
      integer, pointer :: isb(:,:),nbstrbnd(:)
      real*8, pointer ::  sbk(:,:)
      integer :: winisb,winnbstrbnd,winsbk
      save
      end
