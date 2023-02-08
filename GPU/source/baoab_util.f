#include "tinker_macro.h"
      submodule(utilbaoab) baoabutil
      implicit none
      real(r_p):: a1=0.0,a2=0.0

      contains

      module subroutine apply_b_piston(dt,pres)
      use bath
      use units
      implicit none
      real(r_p), intent(in) :: dt,pres
      aextvol = 3.0*convert*(extvol*(pres-atmsph)/prescon 
     &          +gasconst*kelvin)/masspiston
      vextvol = vextvol + dt*aextvol

      end subroutine apply_b_piston

      module subroutine apply_o_piston(dt)
      use bath
      use units
      use random_mod
      use mpi
      use domdec
      implicit none
      real(r_p), intent(in) :: dt
      real(r_p) :: a1piston,a2piston
      integer :: ierr

      if (rank.eq.0) then
        a1piston = exp(-gammapiston*dt)
        a2piston = sqrt((1.d0-a1piston**2)*
     &         boltzmann*kelvin/masspiston)
        vextvol  = a1piston*vextvol + a2piston*normal()
      end if
      call MPI_BCAST(vextvol,1,MPI_RPREC,0,COMM_TINKER,ierr)
      temppiston = masspiston*vextvol**2/boltzmann

      end subroutine apply_o_piston

      module subroutine apply_a_piston(dt,istep,A_full)
      use atomsMirror
      use atmlst
      use moldyn
      use usage
      use domdec
      use freeze
      use bath
      implicit none
      real(r_p),intent(in):: dt
      integer  ,intent(in):: istep
      logical  ,intent(in):: A_full
      integer   i,j,iglob
      real(r_p) tau,scale

      extvolold = extvol
      scale     = exp(dt*vextvol) 
      tau       = 0.0
      if (A_full) tau=sinh(dt*vextvol)/vextvol  
      extvol    = extvol*(scale**3)
      call rescale_box(istep,scale)

!$acc parallel loop async default(present)
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
          x(iglob) = x(iglob)*scale + v(1,iglob)*tau
          y(iglob) = y(iglob)*scale + v(2,iglob)*tau
          z(iglob) = z(iglob)*scale + v(3,iglob)*tau
          v(1,iglob)= v(1,iglob)/scale
          v(2,iglob)= v(2,iglob)/scale
          v(3,iglob)= v(3,iglob)/scale
        end if
      end do
      call reCast_position
      end subroutine apply_a_piston

      module subroutine set_langevin_thermostat_coeff(dta)
      use langevin   ,only: gamma
      use bath       ,only: kelvin
      use units      ,only: boltzmann
      real(r_p),intent(in):: dta
c
c     set coefficients for BAOAB integration
c
      a1 = exp(-gamma*dta)
      a2 = sqrt((1-a1**2)*boltzmann*kelvin)
      end subroutine

      module subroutine apply_langevin_thermostat
      use atmtyp     ,only: mass
      use bath       ,only: kelvin
      use domdec     ,only: rank,glob,nloc
      use langevin   ,only: Rn,gamma
      use moldyn     ,only: v,a
      use random_mod ,only: normalvec,host_rand_platform
#ifdef _OPENACC
     &               ,normalgpu
#endif
      use tinMemory  ,only: prmem_request
      use units      ,only: boltzmann
      use usage
      implicit none
      integer i,j,iglob

      if (a1.eq.0.0) then
         write(0,*) ' langevin thermostat coeff uninitialized'
         __TINKER_FATAL__
      end if

      call prmem_request(Rn,3,nloc+1,async=.false.)
#ifdef _OPENACC
      call normalgpu(Rn(1,1),3*nloc)
#endif
      if (host_rand_platform) then
         call normalvec(Rn,3*nloc)
!$acc update device(Rn) async
      end if

!$acc host_data use_device(glob,use,v,Rn,mass)
!$acc parallel loop collapse(2) async
!$acc&         deviceptr(glob,use,v,Rn,mass)
      do i = 1,nloc; do j = 1,3
         iglob = glob(i)
         if (use(iglob)) then
            v(j,iglob) = a1*v(j,iglob) +
     $      a2*real(Rn(j,i),r_p)/sqrt(mass(iglob))
         end if
      end do; end do
!$acc end host_data
      end subroutine

      end submodule baoabutil
