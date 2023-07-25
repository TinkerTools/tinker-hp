#include "tinker_macro.h"

MODULE deconvolution
  IMPLICIT NONE
  
  ABSTRACT INTERFACE
    FUNCTION convolution_kernel(omega,omega0,gamma,kT)
      IMPLICIT NONE
      REAL(r_p), intent(in) :: omega,omega0,gamma,kT
      REAL(r_p) :: convolution_kernel
    END FUNCTION convolution_kernel
  END INTERFACE
  
  PRIVATE
  PUBLIC :: deconvolute_spectrum,kernel_lorentz, kernel_lorentz_pot &
            ,deconvolute_spectrum_symmetrize

CONTAINS

SUBROUTINE deconvolute_spectrum_symmetrize(s_in,domega,nIterations, &
                trans,thr,gamma,kT,kernel,s_out,s_rec,verbose)
  REAL(r_p),INTENT(in) :: s_in(:),domega,kT,gamma,thr
  INTEGER, INTENT(in) :: nIterations
  REAL(r_p), INTENT(inout) :: s_out(:)
  LOGICAL, INTENT(in) :: trans
  REAL(r_p), INTENT(inout), optional :: s_rec(:)
  LOGICAL, INTENT(in), optional :: verbose
  PROCEDURE(convolution_kernel):: kernel
  integer :: nom,k
  REAL(r_p), allocatable :: s_in_sym(:),omsym(:),s_out_sym(:)
  REAL(r_p), allocatable :: s_rec_sym(:)

  nom = size(s_in)
  allocate(s_in_sym(2*nom-1))
  allocate(omsym(2*nom-1))
  allocate(s_out_sym(2*nom-1))
  allocate(s_rec_sym(2*nom-1))
  omsym(nom)=0.d0
  s_in_sym(nom)=s_in(1)
  do k=1,nom-1
    omsym(nom+k)=k*domega
    omsym(nom-k)=-k*domega
    s_in_sym(nom+k)=s_in(k+1)
    s_in_sym(nom-k)=s_in(k+1)
  enddo

  call deconvolute_spectrum(s_in_sym,omsym,nIterations,trans   &
          ,thr,gamma, kT, kernel, s_out_sym,s_rec_sym,verbose)
  
  s_out(:)=s_out_sym(nom:2*nom-1)
  if(present(s_rec)) then
    s_rec(:)=s_rec_sym(nom:2*nom-1)
  endif

END SUBROUTINE deconvolute_spectrum_symmetrize

!! deconvolute the spectrum 's_in' using the method from ref:
!! J. Chem. Phys. 148, 102301 (2018); 
!! https://doi.org/10.1063/1.4990536@jcp.2018.NQE2018.issue-1
SUBROUTINE deconvolute_spectrum(s_in,omega,nIterations, &
                trans,thr,gamma,kT,kernel,s_out,s_rec,  &
                verbose)
  IMPLICIT NONE
  REAL(r_p),INTENT(in) :: s_in(:),omega(:),kT,gamma,thr
  INTEGER, INTENT(in) :: nIterations
  REAL(r_p), INTENT(inout) :: s_out(:)
  REAL(r_p), INTENT(inout), optional :: s_rec(:)
  LOGICAL, INTENT(in) :: trans
  LOGICAL, INTENT(in), optional :: verbose
  PROCEDURE(convolution_kernel):: kernel
  INTEGER :: nom,i,j,it,info
  REAL(r_p), ALLOCATABLE, SAVE :: D(:,:),K(:,:)
  REAL(r_p), ALLOCATABLE :: h(:),s_next(:),Ktr(:,:)
  REAL(r_p), ALLOCATABLE, SAVE :: K_norm(:,:),omega_norm(:)  
  REAL(r_p) :: domega,om,om0,den,rn,ln,diff
  real(r_p), save :: gamma_save=-1.d0
  LOGICAL :: verb,compute_kernel
  LOGICAL, save :: trans_save=.FALSE.
  
  verb = merge(verbose, .FALSE., mask=present(verbose))
  
  nom=size(omega)
  if(size(s_in)/=nom) STOP "[spectrum_deconvolution.f90] Error: size(s_in)/=nom"
  if(size(s_out)/=nom) STOP "[spectrum_deconvolution.f90] Error: size(s_out)/=nom"
  if(present(s_rec)) then
    if(size(s_rec)/=nom) STOP "[spectrum_deconvolution.f90] Error: size(s_rec)/=nom"
  endif
  
  domega=omega(2)-omega(1)
  
  allocate(h(nom),s_next(nom))
  compute_kernel=.FALSE.
  if(.not. allocated(D)) then
    ALLOCATE(D(nom,nom),K(nom,nom))
    compute_kernel=.TRUE.
  elseif(size(D,1)/=nom) then
    DEALLOCATE(D,K)
    ALLOCATE(D(nom,nom),K(nom,nom))
    compute_kernel=.TRUE.
  endif
  if(trans /= trans_save) compute_kernel=.TRUE.
  trans_save=trans
  if(gamma /= gamma_save) compute_kernel=.TRUE.
  gamma_save=gamma
  if (trans .AND. compute_kernel) then 
    if(.not. allocated(K_norm)) then
      ALLOCATE(K_norm(nom,2*nom-1),omega_norm(2*nom-1))
    elseif(size(K_norm)/=2*nom-1) then
      DEALLOCATE(K_norm,omega_norm)
      ALLOCATE(K_norm(nom,2*nom-1),omega_norm(2*nom-1))
    endif
      omega_norm(nom)=0.d0
      K_norm=0.d0
      do j=1,nom-1
        omega_norm(nom+j)=j*domega
        omega_norm(nom-j)=-j*domega
      enddo
      DO j=1,2*nom-1
        om0=omega_norm(j)
        DO i=1,nom
          om=omega(i)  
          K_norm(i,j)=kernel(om,om0,gamma,kT)
        ENDDO
      ENDDO
  endif

  if(verb .or. compute_kernel) then
    write(*,*) "Initializing deconvolution..."
  endif
  if(compute_kernel) then
    write(*,*) "    Computing kernel matrix..."
    K=0.d0
    DO j=1,nom
      om0=omega(j)
      DO i=1,nom
        om=omega(i)  
        K(i,j)=kernel(om,om0,gamma,kT)
      ENDDO
    ENDDO
    DO j=1,nom
      if (trans) then
        !write(*,*) omega(j)*au%THz,sum(K(j,:))*domega
  !      if (abs(omega(j))<maxval(omega)/5.0d0) K(j,:)=K(j,:)/sum(K(j,:))/domega 
        K(j,:)=K(j,:)/sum(K_norm(j,:))/domega 
      else
        K(:,j)=K(:,j)/sum(K(:,j))/domega 
      endif
    enddo

    !if(verb) then
      write(*,*) "    done!"
      write(*,*) "    Computing double convolution matrix..."
    !endif
    !allocate(Ktr(nom,nom))
    !Ktr=transpose(K)
    !D=matmul(Ktr,K)*domega
    call DGEMM('T','N',nom,nom,nom,domega,K,nom,K,nom,0.d0,D,nom,info)
    !deallocate(Ktr)
    write(*,*) "    done!"  
  endif
  call convolution(K,s_in,domega,h,'T')
  if(verb) write(*,*) "done!"
  
  if(verb) write(*,*) "Starting iterative process"
  s_out(:)=s_in(:)
  s_next(:)=0.d0
  DO it=1,nIterations
    if(verb) write(*,*) "Iteration ",it
    DO i=1,nom
      den=SUM(s_out*D(i,:))*domega
      s_next(i)=s_out(i)*h(i)/den
    ENDDO   
    diff=SUM((s_next-s_out)**2)/sum(s_out**2)
    if(verb) then
      write(*,*) "    relative squared difference: ",diff
      call compute_rn_ln(s_in,s_next,K,domega,rn,ln)
      write(*,*) "    rn=",rn,"ln=",ln
    endif
    s_out= s_next
    if(diff<=thr) exit  
  ENDDO  
  
  if(present(s_rec)) then
    call convolution(K,s_out,domega,s_rec,'N')
  endif
  if(verb) write(*,*) "done!"
  
END SUBROUTINE deconvolute_spectrum

SUBROUTINE convolution(kernel,s_in,domega,s_out,op)
  IMPLICIT NONE
  REAL(r_p), INTENT(in) :: kernel(:,:),s_in(:),domega
  REAL(r_p), INTENT(inout) :: s_out(:)
  character(1), intent(in) :: op
  integer :: nom,info
  nom=size(s_in)
  call DGEMV(op,nom,nom,domega,kernel,nom,s_in,1,0.d0,s_out,1,info)
  !s_out=matmul(kernel,s_in)*domega

END SUBROUTINE convolution

SUBROUTINE compute_rn_ln(s0,sn,K,domega,rn,ln)
  IMPLICIT NONE
  REAL(r_p), INTENT(in) :: s0(:),sn(:),K(:,:),domega
  REAL(r_p), INTENT(out) :: rn,ln
  REAL(r_p), ALLOCATABLE :: fc(:)
  INTEGER :: nom,i
  REAL(r_p) :: domega3
  
  domega3=domega**3
  
  nom=size(s0)
  ALLOCATE(fc(nom))
  call convolution(K,sn,domega,fc,"N")
  
  rn=SUM((fc-s0)**2)/SUM(s0**2)
  
  DEALLOCATE(fc)
  
  ln=0
  DO i=2,nom-1
    ! centered second derivative
    ln=ln+(sn(i+1)-2.d0*sn(i)+sn(i-1))**2/domega3
  ENDDO

END SUBROUTINE compute_rn_ln

FUNCTION kernel_lorentz(omega,omega0,gamma,kT) result(kern)
  IMPLICIT NONE
  REAL(r_p), intent(in) :: omega,omega0,gamma,kT
  REAL(r_p) :: kern
  REAL(r_p), parameter :: pi=4.d0*atan(1.d0)
  REAL(r_p) :: omegamin,omega2
  
  omegamin=TINY(omega)
  
  if(abs(omega)<=omegamin .AND. abs(omega0)<=omegamin) then
    kern = 1.d0/(pi*gamma)
  else  
    omega2=omega*omega
    kern = gamma*omega2/(pi*(omega2*gamma**2 + (omega2-omega0**2)**2))
  endif
  
END FUNCTION kernel_lorentz


FUNCTION kernel_lorentz_pot(omega,omega0,gamma,kT) result(kern)
  IMPLICIT NONE
  REAL(r_p), intent(in) :: omega,omega0,gamma,kT
  REAL(r_p) :: kern
  REAL(r_p), parameter :: pi=4.d0*atan(1.d0)
  REAL(r_p) :: omegamin,omega2
  
  omegamin=TINY(omega)
  
  if(abs(omega)<=omegamin .AND. abs(omega0)<=omegamin) then
    kern = 1.d0/(pi*gamma)
  else  
    omega2=omega*omega
    kern = gamma*omega2/(pi*(omega0**2*gamma**2 + (omega2-omega0**2)**2))
  endif
  
END FUNCTION kernel_lorentz_pot

END MODULE deconvolution



