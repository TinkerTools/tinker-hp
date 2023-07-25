c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Tex_sizeas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  module spectra-- store variables related to spectra      ##
c     ##     and compute IR spectra                                ##
c     ##                                                           ##
c     ###############################################################
c
#include "tinker_precision.h"

      module spectra
      implicit none
      integer :: nseg,nad
      integer ::  startsavespec
      real(r_p) :: omegacut
      real(r_p) :: domega
      real(r_p) :: Tseg
      real(r_p) :: dt_spec
      real(r_p) :: u_fact
      integer :: plan_ir
      logical :: ir=.FALSE.
      logical :: spectra_initialized = .FALSE.
      logical :: kubo=.FALSE.
      logical :: full_dipole=.TRUE.
      logical :: deconvolute_IR = .FALSE.
      integer :: niter_deconvolution = 20
      real(r_p), allocatable :: kubo_fact(:)

      real(r_p), allocatable :: mudot(:,:)
      real(r_p), allocatable :: Cmumu_average(:,:)

      logical :: piqtb_classical_centroid=.FALSE.
      logical :: brieuc_piqtb=.FALSE.


      TYPE compteur_type
          integer :: i_avg
          integer :: i
          integer :: nseg
          integer :: startsave
          integer :: nsample
          logical :: reset_step
      contains
          procedure :: initialize => initialize_compteur
          procedure :: update => update_compteur
      END TYPE compteur_type

      type(compteur_type), private :: compteur

      save

      contains 

      subroutine initialize_spectra(dt)
      use atmtyp
      use atoms
      use bath
      use domdec
      use langevin
      use math
      use moldyn
      use mdstuf, only: integrate
      use sizes
      use usage
      use units
      use mpi
#ifdef _OPENACC
      use cufft
#endif
      use utilgpu, only: prmem_requestm,prmem_request,rec_queue
      use utils, only: set_to_zero1m
      implicit none
      real(r_p), intent(in) ::  dt
      integer, pointer :: pnull
      integer :: ierr    
      real(r_p) :: u  

      if(spectra_initialized) return
      spectra_initialized=.TRUE.

      dt_spec = dt
      nseg=nint(Tseg/dt)
      if(mod(nseg,2)/=0) nseg=nseg+1
      Tseg=nseg*dt
      domega = (2.*pi)/(3.*nseg*dt)
      nad=int(omegacut/domega)
      u_fact = 0.5d0*hbar_planck/(kelvin*boltzmann)
      
      if(ranktot==0)
     &  write(*,*) "initialize spectra with dt=",dt
     &            ,"nseg=",nseg
      call compteur%initialize(nseg,startsavespec)

      if(ir) then
        call prmem_requestm(mudot,3,nseg, async=.TRUE.)
        call set_to_zero1m(mudot,3*nseg,rec_queue)
        if(ranktot==0) then
          if (allocated(Cmumu_average)) deallocate (Cmumu_average)
          allocate(Cmumu_average(3*nseg,3))
          Cmumu_average(:,:)=0.0d0
#ifdef _OPENACC
          pnull=>null()
          ierr = cufftPlanMany(plan_ir,1,(/3*nseg/)
     &          ,pnull,0,0
     &          ,pnull,0,0
     &          ,CUFFT_C2C, 3 )
#endif
        endif
      endif

      end subroutine initialize_spectra

      subroutine set_kubo(piqtb)
      use atmtyp
      use atoms
      use bath
      use domdec
      use langevin
      use math
      use moldyn
      use mdstuf, only: integrate
      use sizes
      use usage
      use units
      use mpi
#ifdef _OPENACC
      use cufft
#endif
      implicit none
      integer :: i,k
      logical, intent(in), optional :: piqtb
      real(r_p) :: u,beta,w2

      if((.not.spectra_initialized) .and. ranktot==0) then
        write(0,*) "Spectra must be initialized to set Kubo"
        call fatal
      endif
      kubo=.TRUE.
      allocate(kubo_fact(nad))
      kubo_fact(1)=1.
      do i=2,nad
        kubo_fact(i)=1.d0/QTB_kernel((i-1)*domega,piqtb,1)
      enddo

      end subroutine set_kubo

      function QTB_kernel(omega,piqtb,ibead)
      use bath
c      use beads, only: omkpi,nbeads
      use units
      implicit none
      real(r_p), intent(in) :: omega
      logical, intent(in) :: piqtb
      integer, intent(in) :: ibead
      real(r_p) :: QTB_kernel
      real(r_p) :: u,w2
      integer :: k

      if (piqtb) then
        QTB_kernel=PIQTB_kernel(omega,ibead)
        return
      endif

      u=u_fact*omega
      QTB_kernel=u/tanh(u)

      end function QTB_kernel

      function PIQTB_kernel(omega,ibead)
      use bath
      use units
      use beads
      use tinheader
      implicit none
      real(r_p), intent(in) :: omega
      integer, intent(in) :: ibead
      real(r_p) :: PIQTB_kernel
      real(r_p) :: u,w2,omegak2,omegak,hu
      real(r_p) :: bead_factor
      integer :: jbead,kk
      real(r_p) :: ufact,a,b,alpha,x0,kd,beta

      if(brieuc_piqtb) then
        PIQTB_kernel=PIQTB_kernel_brieuc(omega,ibead)
        return
      endif

      
      if (omega <= omkpi(ibead)) then
        !if(piqtb_classical_centroid .and. ibead>1) then
        !  bead_factor=0.0d0
        !  do jbead=2,nbeads
        !    bead_factor=bead_factor + 1.d0/(u_fact*omkpi(jbead))**2
        !  enddo
        !  alpha=1.d0/(3.d0*bead_factor)
        !  !kd=PIQTB_kernel(omkpi(ibead)+domega,ibead)
        !  !beta=(kd-alpha)/domega*omkpi(ibead)
        !  !a=3*(alpha-1.d0)-beta
        !  !b=-2*(alpha-1.d0)+beta
        !  a=alpha-1.d0
        !  x0=omega/omkpi(ibead)
        !  PIQTB_kernel=1+a*x0**2 !+b*x0**3
        !  return
        !else
          PIQTB_kernel=1.0_re_p
          return
        !endif
      endif

      omegak2=omega**2-omkpi(ibead)**2
      u=u_fact*sqrt(omegak2)
      hu=u/tanh(u)

      bead_factor=0.0d0
      do jbead=2,nbeads
        bead_factor=bead_factor + omegak2/(omegak2+omkpi(jbead)**2)     
      enddo

      if (piqtb_classical_centroid) then
        if(ibead==1) then
          PIQTB_kernel = 1.0_re_p
        else
          PIQTB_kernel=(hu-1.d0)/bead_factor
        endif
      else
        PIQTB_kernel=hu/(1.d0+bead_factor)
      endif

      end function PIQTB_kernel

      function PIQTB_kernel_brieuc(omega,ibead)
      use bath
      use units
      use beads
      use domdec
      use tinheader
      implicit none
      real(r_p), intent(in) :: omega
      integer, intent(in) :: ibead
      real(r_p) :: PIQTB_kernel_brieuc
      real(r_p) :: u,w2,omegak2,omegak,hu
      real(r_p) :: bead_factor
      integer :: jbead,kk,iom,iter,jbead_start
      real(r_p) :: ufact,a,b,alpha,x0,kd,beta
      logical, save :: first_in=.TRUE.
      integer, save :: nom
      real(r_p), allocatable, save :: F(:),us(:)
      real(r_p), allocatable :: targ(:),Vn(:,:),Fp(:)
      real(r_p), allocatable :: uns(:,:)
      real(r_p) :: umax,un,tmax,nu,b1,err,un1,u12,un2,u2,un12
      real(r_p) :: tmp,u1
      logical :: converged

      if (piqtb_classical_centroid .and. ibead==1) then
        PIQTB_kernel_brieuc=1.0_re_p
        return
      endif

      nu=real(nbeads,r_p)

      if (first_in) then
        first_in=.FALSE.
        if (ranktot==0) 
     &     write(*,*) "--- initializing Brieuc-type PIQTB kernel ---"

        nom=(3*nseg)/2
        do while (.TRUE.)
          umax=(nom-1)*u_fact*domega
          tmax=nu*tanh(umax/nu)/tanh(umax)
          do jbead=2,nbeads
            un=sqrt(umax**2 + (omkpi(jbead)*u_fact)**2)
            tmax = tmax - (umax*tanh(umax/nu))/(un*tanh(un/nu))
          enddo
          if(abs(tmax-1.d0)<1.d-2) exit
          nom=nom*2
        enddo
        if(ranktot==0)
     &     write(*,*) "   PIQTB kernel: nom=",nom," (nseg=",nseg,")"

        b1=1./nu
        allocate(F(nom),targ(nom),Vn(nom,nbeads),us(nom))
        allocate(uns(nom,nbeads),Fp(nom))
        F(:)=1.d0
        Vn(1,:)=0.d0
        us(1)=0.d0
        if (piqtb_classical_centroid) then
          u12=(omkpi(2)*u_fact)**2
          targ(1)=u12/(3*(1-1.d0/nu))
          do iom=2,nom
            u=real(iom-1,r_p)*u_fact*domega
            us(iom)=u
            u2=u*u
            un2=u2 + u12
            tmp=un2*nu/(u/tanh(u/nu)-1.d0)
            targ(iom)=tmp*(u/tanh(u)-1.d0)/u2
            do jbead=3,nbeads
              un2=u2 + (omkpi(jbead)*u_fact)**2
              un=sqrt(un2)
              un1=sqrt(un2-u12)
              uns(iom,jbead)=un1
              Vn(iom,jbead)=tmp*(un1/tanh(un1/nu)-1.d0)/(nu*un2)
            enddo
          enddo
          Vn(1,:)=Vn(2,:)
          targ(:)=targ(:)*(nu-1.d0)/nu
        else
          targ(1)=1.d0
          do iom=2,nom
            u=real(iom-1,r_p)*u_fact*domega
            us(iom)=u
            targ(iom)=nu*tanh(u/nu)/tanh(u)
            do jbead=2,nbeads
              un=sqrt(u**2 + (omkpi(jbead)*u_fact)**2)
              uns(iom,jbead)=un
              Vn(iom,jbead)=(u*tanh(u/nu))/(un*tanh(un/nu))
            enddo
          enddo
        endif
        converged=.FALSE.
        jbead_start=2
        if (piqtb_classical_centroid) jbead_start=3
        do iter=1,200
          do iom=1,nom
            Fp(iom)=0.d0
            do jbead=jbead_start,nbeads
              Fp(iom)=Fp(iom)+Vn(iom,jbead)
     &     *linear_interpolation(uns(iom,jbead),us,F,F(1),F(nom))
            enddo
          enddo
          err = maxval(abs(F+Fp-targ)/targ)
          if (err<1.d-3) then
            if (ranktot==0) then
              write(*,*) "   PIQTB kernel: converged in "
     &                       ,iter," iterations"
            endif
            converged=.TRUE.
            exit
          endif
          F = b1*(targ-Fp) + (1.d0-b1)*F
        enddo
        if (.not.converged .and. ranktot==0) then
          write(*,*) "   PIQTB kernel: not converged !"
          call fatal
        endif

        deallocate(targ,Vn,uns,Fp)
      endif

      if(omega*u_fact<us(2)) then
        PIQTB_kernel_brieuc=1.0_re_p
        return
      endif

      u = omega*u_fact
      if(piqtb_classical_centroid) then
        u12=u*u-(omkpi(2)*u_fact)**2
        if(u12<=0.d0) then
          PIQTB_kernel_brieuc=1.0_re_p
          return
        endif
        u1=sqrt(u12)
        PIQTB_kernel_brieuc = (u1/tanh(u1/nu)-1.d0)/(nu-1.d0)
     &     *linear_interpolation(u1,us,F,1._re_p,F(nom))
      else
        PIQTB_kernel_brieuc = u/(nu*tanh(u/nu))
     &    *linear_interpolation(u,us,F,1._re_p,F(nom))
      endif

      end function PIQTB_kernel_brieuc

      function linear_interpolation(x,x_ref,y_ref,left,right)
      ! linearly interpolates y(x) from array of y_ref(x_ref)
      implicit none
      real(r_p), intent(in) :: x,x_ref(:),y_ref(:)
      real(r_p), intent(in) :: right,left
      real(r_p) :: linear_interpolation
      integer :: i
      real(r_p) :: alpha

      if(x<x_ref(1)) then
        linear_interpolation=left
        return
      elseif(x>x_ref(size(x_ref))) then
        linear_interpolation=right
        return
      endif

      do i=2,size(x_ref)
        if(x<x_ref(i)) then
          alpha = (x-x_ref(i-1))/(x_ref(i)-x_ref(i-1))
          linear_interpolation=(1.d0-alpha)*y_ref(i-1)+alpha*y_ref(i)
          return
        endif
      enddo

      end function linear_interpolation

      subroutine save_dipole_respa(dip,dipind)
        use moldyn, only: dshort,dinter,nalt
        use uprior
        implicit none
        real(r_p), intent(in) :: dip(3,nalt),dipind(3,nalt)
        real(r_p) :: dipindprev(3),dipindnew(3),dipindint(3,nalt)
        real(r_p), save :: dipsave(3,maxualt)=0.0d0
        real(t_p) :: x(maxualt),y(maxualt)
        real(t_p) :: b(maxualt),c(maxualt),d(maxualt)
        integer, save :: nsave = 0
        integer :: ialt,j,k
        real(r_p) :: alpha

        if(.not.full_dipole) then
          dipindint(:,1) = 0.0d0
          do ialt=1,nalt
            call save_dipole_traj(dip(:,ialt),dipindint(:,1))
          enddo
          return
        endif
        dipindint(:,nalt)=dipind(:,nalt) !+dip(:,nalt)

        if(nsave==0) then
          dipsave(:,2)=dipind(:,nalt) !+dip(:,nalt)
          nsave=1
          return
        elseif(nsave==1) then
          dipsave(:,1)=dipind(:,nalt) !+dip(:,nalt)
          nsave=2
          !LINEAR INTERPOLATION
          do ialt=1,nalt-1
            alpha = real(ialt,r_p)/real(nalt,r_p)
            dipindint(:,ialt)=dipsave(:,2)
     &       + alpha*(dipsave(:,1)-dipsave(:,2))
          enddo
        else
          !move up previously saved values
          nsave = min(nsave+1,maxualt)
          do k=nsave,2,-1
            dipsave(:,k)=dipsave(:,k-1)
          enddo
          dipsave(:,1)=dipind(:,nalt) !+dip(:,nalt)
          !SPLINE INTERPOLATION
          do ialt=1,nsave; x(ialt)=real(ialt-1,t_p); enddo
          do j=1,3
            y=dipsave(j,1:nsave)
            call cubic_spline(nsave,x,y,b,c,d)
            do ialt=1,nalt-1
              alpha = 1.d0 - real(ialt,r_p)/real(nalt,r_p)
              dipindint(j,ialt) = y(1) + b(1)*alpha 
     &           + c(1)*alpha**2 + d(1)*alpha**3
            enddo
          enddo
        endif

        do ialt=1,nalt
          call save_dipole_traj(dip(:,ialt),dipindint(:,ialt))
        enddo

      end subroutine save_dipole_respa

      subroutine save_dipole_traj(dip,dipind)
      use mpi
      use domdec
      use moldyn, only: dshort,dinter,nalt
      use random_mod
      use uprior
      use potent, only: use_polar
      implicit none
      integer :: iseg,i,j,iipole,iglob,iichg,ierr,k
      real(r_p), intent(in) :: dip(3),dipind(3)
      real(r_p) :: dipx,dipy,dipz
      real(r_p), save :: dipsave(3,2)=0.0d0
      real(r_p) :: diptot(3),dipdot(3),alpha
      logical :: first_in=.TRUE.

      diptot(:) = dip(:)
      if(full_dipole) diptot(:) = diptot(:) + dipind(:)

      if(first_in) then
        dipsave(:,2)=diptot
        dipsave(:,1)=diptot
        first_in=.FALSE.
        return
      endif

      dipdot(:)=(diptot(:)-dipsave(:,2))/(2*dt_spec)

      !move up previously saved values
      dipsave(:,2)=dipsave(:,1)
      dipsave(:,1)=diptot

      dipx=dipdot(1)
      dipy=dipdot(2)
      dipz=dipdot(3)
      iseg = compteur%i
!$acc serial async present(mudot)
        mudot(1,iseg)=dipx
        mudot(2,iseg)=dipy
        mudot(3,iseg)=dipz
!$acc end serial

      call compteur%update()
      if(compteur%reset_step) call irspectra()

      end subroutine save_dipole_traj

      subroutine irspectra()
      use atmtyp
      use atoms
      use bath
      use domdec
      use langevin
      use deconvolution
      use math
      use moldyn
      use mdstuf, only: integrate
      use sizes
      use usage
      use units
      use boxes, only: volbox
      use tinheader
#ifdef _OPENACC
      use cufft
#else
      use glassman
#endif
      use utilgpu, only: rec_queue
      use utils, only: set_to_zero1m
      use mpi
      implicit none
      integer :: i,j,k,l,ierr,idof
      real(r_p) :: mufact
      complex(4),  allocatable :: smu(:,:)
      real(r_p),  allocatable :: Cmumu(:,:)
      real(r_p),  allocatable :: Cmumu_sum(:)
      real(r_p),  allocatable :: Cmumu_dec(:)
      real(r_p),  allocatable :: Cmumu_rec(:)
#ifndef _OPENACC
      complex(t_p), allocatable :: smu1d(:),workfft(:)
#endif

!$acc wait
!$acc host_data use_device(mudot)
      if (ranktot.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,mudot,3*nseg,MPI_RPREC
     &      ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      else
        call MPI_REDUCE(mudot,mudot,3*nseg,MPI_RPREC
     &      ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      endif
!$acc end host_data

      if(ranktot/=0) then
        call set_to_zero1m(mudot,3*nseg,rec_queue)
        return
      endif

      allocate(smu(3*nseg,3))
      allocate(Cmumu(3*nseg,3))
!$acc data create(smu) copyout(Cmumu)
!$acc&     present( mudot)

!$acc parallel loop collapse(2)
      do k=1,3*nseg; do j=1,3
        if(k<=nseg) then
          smu(k,j)=cmplx(mudot(j,k),0.)
          mudot(j,k) = 0._re_p
        else
          smu(k,j)=cmplx(0.,0.)
        endif
      enddo; enddo

#ifdef _OPENACC
!$acc wait
!$acc host_data use_device(smu)
      ierr = cufftExecC2C(plan_ir,smu,smu,CUFFT_FORWARD)
!$acc end host_data
      if(ierr /= 0 ) then
        write(0,*) "Error in FFT of the dipole moment"
        call fatal
      endif
#else
      allocate(smu1d(3*nseg),workfft(3*nseg))
      do j=1,3
        smu1d(:)=smu(:,j)
        call SPCFFT(smu1d,3*nseg,-1,workfft)
        smu(:,j)=smu1d(:)
      enddo
      deallocate(smu1d,workfft)
#endif


!$acc parallel loop collapse(2)
      do j=1,3; do k=1,3*nseg
        Cmumu(k,j) = real(abs(smu(k,j)),r_p)**2*dt_spec/(3*nseg)
      enddo; enddo

!$acc end data

      if(kubo) then
        do j=1,3; do k=1,nad
          Cmumu(k,j)=Cmumu(k,j)*kubo_fact(k)
        enddo; enddo
      endif

      if(compteur%i_avg.le.compteur%startsave) then
        Cmumu_average(:,:)=0.d0
      endif
      Cmumu_average(:,:)=Cmumu_average(:,:)+Cmumu(:,:)


      ! constants+conversion factors => IR spectrum in cm^-1
      mufact = convert*coulomb *(2.d0*pi)
     &         *pi/(boltzmann*kelvin*3.d0*lightspd*volbox)
      Cmumu(:,:)=mufact*Cmumu_average(:,:)/real(compteur%nsample)
      do j=1,3
        call correct_window_spectrum(Cmumu(:,j))
      enddo

      allocate(Cmumu_sum(nad),Cmumu_dec(nad),Cmumu_rec(nad))
      Cmumu_sum(:)=Cmumu(1:nad,1)+Cmumu(1:nad,2)+Cmumu(1:nad,3)
      if(compteur%i_avg > compteur%startsave
     &    .AND. deconvolute_IR  ) then
        call deconvolute_spectrum_symmetrize(Cmumu_sum,domega
     &        ,niter_deconvolution, .FALSE., real(1.0e-10,r_p)
     &        ,gamma ,boltzmann*kelvin,kernel_lorentz
     &        ,Cmumu_dec,Cmumu_rec,verbose=.FALSE.)
      else
        Cmumu_dec = Cmumu_sum
        Cmumu_rec = Cmumu_sum
      endif
        
      open(13,file='IR_spectra.out')
      do k=1,nad
        write(13,'(4e16.8)') domega*(k-1)*cm1
     &       ,Cmumu_sum(k), Cmumu_dec(k), Cmumu_rec(k)
      enddo   
      close(13)     

      end subroutine irspectra

      subroutine update_compteur(this)
        implicit none
        class(compteur_type), intent(inout) :: this

        this%i=this%i+1
        if(this%i > this%nseg) then    
          this%i=1      
          this%i_avg=this%i_avg+1
          this%nsample = max(this%i_avg-this%startsave+1,1)
          this%reset_step = .TRUE.
        else
          this%reset_step=.FALSE.
        endif

      end subroutine update_compteur

      subroutine initialize_compteur(this,nseg,startsave)
        implicit none
        class(compteur_type), intent(inout) :: this
        integer, intent(in) :: startsave,nseg

        this%reset_step = .TRUE.
        this%i=1
        this%i_avg = 0
        this%nsample = 1
        if(startsave<=0) then
          this%startsave=1
        else
          this%startsave = startsave
        endif
        this%nseg = nseg

      end subroutine initialize_compteur

      subroutine correct_window_spectrum(spec)
        use glassman
        implicit none
        real(r_p), intent(inout) :: spec(3*nseg)
        complex(t_p), allocatable :: s_fft(:),workfft(:)
        integer :: k
        integer :: win,winmax
        real(r_p) :: corr(nseg),hann

        win=nint(3.d0*nseg/4.d0)+1
        allocate(s_fft(3*nseg))
        allocate(workfft(3*nseg))

        s_fft=cmplx(spec,0.d0)
        call SPCFFT(s_fft,3*nseg,-1,workfft)
        corr(1)=1.d0
        do k=2,nseg
          corr(k) = real(nseg,r_p)/real(nseg-k+1,r_p)
        enddo
        do k=1,win
          hann=sin(0.5*acos(-1.d0)*(1.d0+real(k,r_p)
     &       /real(win,r_p)))**2
          corr(nseg-win+k)=corr(nseg-win+k)*hann + 1.d0 - hann
        enddo

        do k=2,nseg
          s_fft(k)=s_fft(k)*corr(k)
          s_fft(3*nseg-k+2)=s_fft(3*nseg-k+2)*corr(k)
        enddo
        call SPCFFT(s_fft,3*nseg,1,workfft)
        spec(:)=real(s_fft,r_p)/real(3*nseg,r_p)

      end subroutine correct_window_spectrum

      subroutine compute_dipole(dip,dipind,full)
      use atmlst
      use atomsmirror
      use boxes
      use charge
      use domdec
      use iounit
      use mpole
      use polar
      use potent
      use tinheader ,only:ti_p,re_p
      use units
      use mpi
      implicit none
      integer i,iipole,iglob,ierr
      integer iichg
      real(r_p), intent(inout) ::  dip(3),dipind(3)
      real(r_p), save ::  dipx,dipy,dipz
      real(r_p), save::  dipindx,dipindy,dipindz
      logical, intent(in) :: full
      real(t_p) q,xr,yr,zr
      real(t_p) mux,muy,muz,mudx,mudy,mudz,mupx,mupy,mupz
      logical, save :: first_in=.TRUE.
c
      if (first_in) then
!$acc enter data create(dipx,dipy,dipz,dipindx,dipindy,dipindz)
        first_in=.FALSE.
      endif
!$acc data present(dipx,dipy,dipz,dipindx,dipindy,dipindz
!$acc&             ,dip, dipind) async

!$acc serial async
      dipx = 0.0_re_p
      dipy = 0.0_re_p
      dipz = 0.0_re_p
      dipindx = 0.0_re_p
      dipindy = 0.0_re_p
      dipindz = 0.0_re_p
!$acc end serial
      if (use_mpole) then
!$acc parallel loop default(present) async
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = ipole(iipole)
          q = rpole(1,iipole)
          mux = q*x(iglob) !+ rpole(2,iipole)
          muy = q*y(iglob) !+ rpole(3,iipole)
          muz = q*z(iglob) !+ rpole(4,iipole)
          if (full) then
            mudx = uind(1,iipole)
            mudy = uind(2,iipole)
            mudz = uind(3,iipole)
            mupx = uinp(1,iipole)
            mupy = uinp(2,iipole)
            mupz = uinp(3,iipole)
            dipindx = dipindx + 0.5*(mudx+mupx) + rpole(2,iipole)
            dipindy = dipindy + 0.5*(mudy+mupy) + rpole(3,iipole)
            dipindz = dipindz + 0.5*(mudz+mupz) + rpole(4,iipole)
          endif
          dipx = dipx + mux
          dipy = dipy + muy
          dipz = dipz + muz
        end do
      else if (use_charge) then
!$acc parallel loop default(present) async
        do i = 1, nionloc
          iichg = chgglob(i)
          iglob = iion(iichg)
          xr = x(iglob)
          yr = y(iglob)
          zr = z(iglob)
          q = pchg(iichg)
          dipx = dipx + q*xr 
          dipy = dipy + q*yr 
          dipz = dipz + q*zr 
        end do
      end if

!$acc serial async
      dip(1) = dipx
      dip(2) = dipy
      dip(3) = dipz
      dipind(1) = dipindx
      dipind(2) = dipindy
      dipind(3) = dipindz
!$acc end serial
!$acc end data
      !dip = [dipx,dipy,dipz]
      !dipind = [dipindx,dipindy,dipindz]
      end subroutine compute_dipole

      end module spectra
