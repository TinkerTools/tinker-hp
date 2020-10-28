c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  function maxwell  --  Maxwell-Boltzmann distribution value  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "maxwell" returns a speed in Angstroms/picosecond randomly
c     selected from a 3-D Maxwell-Boltzmann distribution for the
c     specified particle mass and system temperature
c
c     literature reference:
c
c     P. W. Atkins, "Physical Chemistry, 4th Edition", W. H. Freeman,
c     New York, 1990; see section 24.2 for general discussion
c
c
#include "tinker_precision.h"
      function maxwell (mass,temper)
      use erf_mod
      use random_mod
      use units
      use tinheader ,only: ti_p,prec_eps
      implicit none
      real(r_p) maxwell
      real(r_p) mass
      real(r_p) temper
      real(t_p) rho
      real(r_p) beta
      real(r_p) xspeed,yspeed
      real(r_p) zspeed
c     real(t_p) random  !!!work on pgi compiler
c
c     set normalization factor for cumulative velocity distribution
c
      beta = sqrt(mass / (2.0_ti_p*boltzmann*temper))
c
c     pick a randomly distributed velocity along each of three axes
c
      rho = random ()
      if (rho.eq.1.0_ti_p) rho=rho-prec_eps  !Correct if equals to one
      xspeed = real(erfinv(rho),r_p) / beta
      rho = random ()
      if (rho.eq.1.0_ti_p) rho=rho-prec_eps
      yspeed = real(erfinv(rho),r_p) / beta
      rho = random ()
      if (rho.eq.1.0_ti_p) rho=rho-prec_eps
      zspeed = real(erfinv(rho),r_p) / beta
c
c     set the final value of the particle speed in 3-dimensions
c
      maxwell = sqrt(xspeed**2 + yspeed**2 + zspeed**2)
      return
      end
c
c     Device version of maxwell up define
c
#ifdef _OPENACC
      subroutine maxwellgpu (mass,temper,nloc,max_result)
      use atoms,only:n
      use domdec,only:glob
      use erf_mod
      use langevin
      use tinheader ,only: prec1_eps,ti_p
      use random_mod
      use units
      implicit none
      integer i,j,iglob
      integer  ,intent(in ):: nloc
      real(r_p),intent(in ):: mass(n)
      real(r_p),intent(in ):: temper
      real(r_p),intent(out):: max_result(nloc)
      real(r_p) beta
      real(r_p) xspeed,yspeed,zspeed
      logical :: check=.false.
c     real(t_p) random  !!!work on pgi compiler
      
      if (.not.allocated(Rn)) then
         allocate(Rn(3,nloc))
         check = .true.
!$acc enter data create(Rn) async
      end if
c
c     pick a randomly distributed velocity along each of three axes
c
!$acc data present(Rn,glob,mass,max_result) async
c
      call randomgpu(Rn(1,1),3*nloc)
      ! Make correction in Rn if 1.0 value appears.
      ! log(0) cannot be computed in erfinv
      ! -- Tend to happen with hugh sample
!$acc parallel loop collapse(2) async
      do i = 1, nloc
         do j = 1,3
            if (Rn(j,i).eq.1.0_ti_p)
     &         Rn(j,i) = Rn(j,i) - 2*prec1_eps
         end do
      end do
!$acc parallel loop async
      do i = 1, nloc
         iglob = glob(i)
c
c        set normalization factor for cumulative velocity distribution
c
         beta   = sqrt(mass(iglob) / (2.0_re_p*boltzmann*temper))
         xspeed = real(erfinv(Rn(1,i)),r_p) / beta
         yspeed = real(erfinv(Rn(2,i)),r_p) / beta
         zspeed = real(erfinv(Rn(3,i)),r_p) / beta
c
c        set the final value of the particle speed in 3-dimensions
c
         max_result(i) =  sqrt(xspeed**2 + yspeed**2 + zspeed**2)
      end do
c
!$acc end data

      if (check) then
!$acc exit data delete(Rn) async
         deallocate(Rn)
      end if

      end
#endif
