c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  module fields  --  molecular mechanics force field description  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     biotyp       force field atom type of each biopolymer type
c     forcefield   string used to describe the current forcefield
c
c
      module fields
      use sizes
      implicit none
      integer biotyp(maxbio)
      character*20 forcefield
      save
      end
