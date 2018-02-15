c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  strbnd.i  --  stretch-bends in the current structure  ##
c     ##                                                        ##
c     ############################################################
c
c
c     sbk       force constants for stretch-bend terms
c     nstrbnd   total number of stretch-bend interactions
c     isb       angle and bond numbers used in stretch-bend
c
c
      integer nstrbnd,nstrbndloc
      integer, pointer :: isb(:,:),nbstrbnd(:)
      real*8, pointer ::  sbk(:,:)
      common /strbnd/ sbk,nstrbnd,nstrbndloc,isb,nbstrbnd
