
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Tex_sizeas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine adHnoise  --  modify the random noise kernel  ##
!     ##                                                           ##
!     ###############################################################
!
!     Adaptation of our random values and Ht arrays to create
!     colored noise 
!
!     Literature reference: 
!     J-L Barrat & D. Rodney JStatPhys (2011)
!     and H Dammak PRL 103, 190601 (2009

      subroutine adHnoise(dt)
      use adqtb
      use atmtyp
      use qtb
      use atoms
      use bath
      use cutoff
      use domdec
      use energi
      use freeze
      use langevin
      use math
      use mdstuf
      use moldyn
      use timestat
      use usage
      use units
      use mpi

      implicit none
      integer i,j,k,l,iglob
      real*8  C_omega   !sinus cardinal
      real*8 f_omega   !fermi dira! distribution
      real*8  theta_tilde
      real*8 dt
      real*8 h
      real*8 k_b                          !bolzmann in good units
      real*8 omega
      real*8 t
      real*8 normal
      real*8 g_ratio

!     constant in Kcal/mol and ps
!
      !hbar=(planck*1.d11*avogadro)/(2*pi)
      k_b=boltzmann
     
!     Arguments
!
!      a1=exp(-gamma*dt)
!      a2=sqrt(1-a1**2)
!
c      open(1,file='H_omegaad.out')


      Htilde=0.0

      do i=1,typemax
        do k=0,(3*nseg)/2
          omega=k*domega
            if  (k .eq.0) then
              Htilde(k,i)=sqrt(k_b*kelvin*(1.0d0+corr_pot_ratio(k)))
     &                    *gamma_type(k,i)/gamma
            else if (noQTB) then
              Htilde(k,i)=sqrt(k_b*kelvin)
              Htilde(3*nseg-k,i)=sqrt(k_b*kelvin)
            else
              C_omega=(1-2*exp(-gamma*dt)*cos(omega*dt)+exp(-2*gamma*dt)
     &                )/((gamma**2+omega**2)*dt**2)
              theta_tilde=hbar*abs(omega)*(1.0/2.0+1.0
     &                   /(exp(hbar*abs(omega)/(k_b*kelvin))-1))
              f_omega=1.0/(1+exp((abs(omega)-omegacut)/omegasmear))
              g_ratio=1.
              if(k.lt.nad-1) then
                  g_ratio=gamma_type(k,i)/gamma
     &            *(1.0d0+corr_pot_ratio(k))
              endif
              Htilde(k,i)=sqrt(theta_tilde*f_omega*C_omega*g_ratio)
              Htilde(3*nseg-k,i)=sqrt(theta_tilde
     &                      *f_omega*C_omega*g_ratio)
c                              write(1,'(2e20.5)') omega,
c    &                        Htilde(k,i)
              endif
c                  write(1,'(6e20.5)') omega, C_omega, theta
c     &                                _tilde, f_omega,Htilde(k)
c     &                  ,ns*dt*C_omega*theta_tilde*f_omega

        enddo
      enddo

      end
