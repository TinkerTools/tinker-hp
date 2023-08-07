c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  module virial  --  components of internal virial tensor  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     vir    total internal virial Cartesian tensor components
c     virsave stored internal virial Cartesian tensor components
c
c
      module virial
      implicit none
      real*8 dedv
      real*8 vir(3,3)
      real*8 virsave(3,3)
      logical virnum
      logical kin_instant
      save
      end
