c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module mdstuf1  --  control of molecular dynamics trajectory  ##
c     ##                                                               ##
c     ###################################################################
c
c

c     derivs  stores forces computes by gradient routine
c     etot    holds the system total energy at current timestep
c     epot    holds the system potential energy at current timestep
c     eksum   holds the system kinetic energy at current timestep
c     temp    holds the system temperature at current timestep
c     pres    holds the system pressure at current timestep
c
      module mdstuf1
      implicit none
      real*8,allocatable::derivs(:,:)
      real*8 etot,epot,eksum,ealt,ealt2,eml,ealtml
      real*8 temp,pres
      real*8 ekin(3,3),stress(3,3),viralt(3,3),viralt2(3,3)

      save
      end module mdstuf1
