#include "tinker_macro.h"
      submodule(utilbaoab) baoabutil
      implicit none
      real(r_p) :: a1=0.0,a2=0.0
      real(r_p) :: acc_nose, an

      contains

      module subroutine apply_b_piston(dt,pres,stress)
      use bath
      use units
      implicit none
      real(r_p), intent(in) :: dt,pres,stress(3,3)
      integer :: i

      if (anisotrop) then
        do i=1,3
          aextbox(i) = convert*(extvol*(stress(i,i)-atmsph)/prescon 
     &          +gasconst*kelvin)/masspiston
          vextbox(i) = vextbox(i) + dt*aextbox(i)
        enddo

      else
        aextvol = 3.0*convert*(extvol*(pres-atmsph)/prescon 
     &          +gasconst*kelvin)/masspiston
        vextvol = vextvol + dt*aextvol
      endif

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
      integer :: ierr,i

      if (rank.eq.0) then
        a1piston = exp(-gammapiston*dt)
        a2piston = sqrt((1.d0-a1piston**2)*
     &         boltzmann*kelvin/masspiston)
        if(anisotrop) then
          do i=1,3
            vextbox(i)=a1piston*vextbox(i)+a2piston*normal()
          enddo
        else
          vextvol=a1piston*vextvol + a2piston*normal()
        endif
      end if
      if(nproc>1) then
        if (anisotrop) then
         call MPI_BCAST(vextbox,3,MPI_RPREC,0,COMM_TINKER,ierr)
        else
         call MPI_BCAST(vextvol,1,MPI_RPREC,0,COMM_TINKER,ierr)
        endif
      endif
      if (anisotrop) then
        temppiston = masspiston*sum(vextbox**2)/(3.*boltzmann)      
      else
        temppiston = masspiston*vextvol**2/boltzmann
      endif

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
      real(r_p) scale(3),scalex,scaley,scalez
      real(r_p) tau(3),taux,tauy,tauz

      extvolold = extvol
      tau       = 0.0
      if(anisotrop) then
       scale(:) = exp(dt*vextbox(:))
       extbox(:)=extbox(:)*scale(:)
       if (A_full) tau=sinh(dt*vextbox)/vextbox  
      else       
       scale(:) = exp(dt*vextvol) 
       if (A_full) tau=sinh(dt*vextvol)/vextvol  
      endif
      extvol    = extvol*scale(1)*scale(2)*scale(3)

      call rescale_box(istep,scale)
      scalex=scale(1); scaley=scale(2); scalez=scale(3)
      taux=tau(1); tauy=tau(2); tauz=tau(3)

!$acc parallel loop async default(present)
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
          x(iglob) = x(iglob)*scalex + v(1,iglob)*taux
          y(iglob) = y(iglob)*scaley + v(2,iglob)*tauy
          z(iglob) = z(iglob)*scalez + v(3,iglob)*tauz
          v(1,iglob)= v(1,iglob)/scalex
          v(2,iglob)= v(2,iglob)/scaley
          v(3,iglob)= v(3,iglob)/scalez
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

      module subroutine apply_langevin_thermostat(dt)
      use atomsmirror,only: n
      use atmtyp     ,only: mass
      use bath       ,only: kelvin
      use domdec     ,only: rank,glob,nloc
      use langevin
      use moldyn     ,only: v,a
      use qtb
      use random_mod ,only: normalvec,host_rand_platform,normal
#ifdef _OPENACC
     &               ,normalgpu
#endif
      use tinMemory  ,only: prmem_request
      use units      ,only: boltzmann
      use usage
      implicit none
      real(r_p) ,intent(in):: dt
      integer i,j,iglob
      real(r_p) :: rnose,dt2
      real(r_p) :: kT,nkT
      logical, save :: f_in=.true.

      if (a1.eq.0.0) then
         write(0,*) ' langevin thermostat coeff uninitialized'
         __TINKER_FATAL__
      end if

      kT = boltzmann*kelvin
      dt2= 0.5*dt

      if(qtb_thermostat) then        
        call qtbrandom()  
      elseif(use_noselangevin) then
        if (f_in) then
!$acc enter data create(acc_nose,an) async
          f_in = .false.
        endif
!$acc serial async present(acc_nose,an)
        acc_nose = 0.0d0
        an = exp(-nose*dt2)
!$acc end serial
!$acc parallel loop collapse(2) async default(present)
!$acc& present(acc_nose)
        do i = 1, nloc
          do j = 1, 3
              iglob = glob(i)
              if (use(iglob)) then
                v(j,iglob) = an*v(j,iglob)
                acc_nose = acc_nose + mass(iglob)*v(j,iglob)**2
              end if
          end do
        end do
        rnose=a2*normal()/sqrt(nose_mass)
        nkT=3*n*kT
!$acc serial async
        acc_nose = dt2*(acc_nose - nkT)/nose_mass
        !nose = nose  + acc_nose
        !nose = a1*nose + rnose
        !nose = nose  + acc_nose
        nose = a1*(nose + acc_nose) + rnose + acc_nose
        an = exp(-nose*dt2)
!$acc end serial
!$acc parallel loop collapse(2) async default(present)
        do i = 1, nloc
          do j = 1, 3
              iglob = glob(i)
              if (use(iglob)) then
                v(j,iglob) = an*v(j,iglob)
              end if
          end do
        end do
      elseif(use_noselangevin_massive) then
        call prmem_request(Rn,3,nloc+1,async=.true.)
#ifdef _OPENACC
        if(.not. host_rand_platform) then
          call normalgpu(Rn,3*nloc)
        endif
#endif
        if (host_rand_platform) then
          do i = 1, nloc; do j = 1, 3
            Rn(j,i) = normal()
          enddo; enddo
!$acc update device(Rn) async
        end if
!$acc parallel loop collapse(2) async default(present)
        do i = 1, nloc; do j = 1, 3
          iglob = glob(i)
          if (use(iglob)) then
            noses(j,iglob) = noses(j,iglob) + dt2/nose_mass
     &           *(mass(iglob)*v(j,iglob)**2- kT)
            v(j,iglob) = exp(-noses(j,iglob)*dt2)*v(j,iglob)
            noses(j,iglob) = a1*noses(j,iglob) 
     &             + a2*Rn(j,i)/sqrt(nose_mass)
            v(j,iglob) = exp(-noses(j,iglob)*dt2)*v(j,iglob)
            noses(j,iglob) = noses(j,iglob) + dt2/nose_mass
     &           *(mass(iglob)*v(j,iglob)**2- kT)
          end if
        end do; end do 
               
      else
         call prmem_request(Rn,3,nloc+1,async=.false.)
#ifdef _OPENACC
         call normalgpu(Rn,3*nloc)
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
     &         a2*real(Rn(j,i),r_p)/sqrt(mass(iglob))
            end if
         end do; end do
!$acc end host_data
      endif
      end subroutine

      end submodule baoabutil
