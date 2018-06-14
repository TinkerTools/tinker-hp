c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1" : driver for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine empole1
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'potent.i'
      include 'openmp.i'
c
c
c     choose the method for summing over multipole interactions
c
      call empole1pme
c
c     zero out energy and derivative terms which are not in use
c
      if (.not. use_mpole) then
         em = 0.0d0
      end if
      if (.not. use_polar) then
         ep = 0.0d0
      end if
      return
      end
