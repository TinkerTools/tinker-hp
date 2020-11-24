c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  module usage  --  atoms active during energy computation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nuse   total number of active atoms in energy calculation
c     iuse   numbers of the atoms active in energy calculation
c     winiuse window object corresponding to iuse
c     use    true if an atom is active, false if inactive
c     winuse window object corresponding to use
c
c
      module usage
      implicit none
      integer nuse
      integer, pointer :: iuse(:)
      logical, pointer :: use(:)
      integer winuse,winiuse
      save
      end
