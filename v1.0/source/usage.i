c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  usage.i  --  atoms active during energy computation  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     nuse   total number of active atoms in energy calculation
c     iuse   numbers of the atoms active in energy calculation
c     use    true if an atom is active, false if inactive
c
c
      integer nuse
      integer, pointer :: iuse(:)
      logical, pointer :: use(:)
      common /usage/ nuse,iuse,use
