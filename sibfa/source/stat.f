c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c   module stat : average values during a dynamic
c
c
      module stat
      implicit none
      real*8 etot_sum,etot2_sum
      real*8 eint_sum,eint2_sum
      real*8 epot_sum,epot2_sum
      real*8 ekin_sum,ekin2_sum
      real*8 temp_sum,temp2_sum
      real*8 pres_sum,pres2_sum
      real*8 dens_sum,dens2_sum
      real*8 vol_sum,vol2_sum
      real*8 pistontemp_sum,pistontemp2_sum
      real*8 etotpi_sum,etot2pi_sum
      real*8 eintpi_sum,eint2pi_sum
      real*8 epotpi_sum,epot2pi_sum
      real*8 ekinpi_sum,ekin2pi_sum
      real*8 temppi_sum,temp2pi_sum
      real*8 prespi_sum,pres2pi_sum
      real*8 denspi_sum,dens2pi_sum
      real*8 volpi_sum,vol2pi_sum
      real*8 pistontemppi_sum,pistontemp2pi_sum
      save
      end
