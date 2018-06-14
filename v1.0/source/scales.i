c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  scales.i  --  parameter scale factors for optimization  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     scale      multiplicative factor for each optimization parameter
c     set_scale  logical flag to show if scale factors have been set
c
c
      real*8 scale
      logical set_scale
      common /scales/ scale(maxvar),set_scale
