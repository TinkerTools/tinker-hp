
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Tex_sizeas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine qtbinit  --  initialize a dynamics trajectory ##
!     ##                                                           ##
!     ###############################################################
!
!     Initialisation of our random values and Ht arrays to create
!     colored noise 
!
!     Literature reference: 
!     J-L Barrat & D. Rodney JStatPhys (2011)
!     and H Dammak PRL 103, 190601 (2009

      subroutine qtbinit(dt)
      use atmtyp
       use adqtb
       use atoms
       use bath
       use cutoff
       use domdec
       use energi
       use freeze
       use katoms
       use langevin
       use math
       use mdstuf
       use moldyn
       use sizes
       use timestat
       use usage
       use units
       use mpi

      implicit none
      integer i,j,k,l,iglob
      integer compteur_lines
      real*8  C_omega   !sinus cardinal
      real*8 f_omega   !fermi dira! distribution
      real*8  theta_tilde
      real*8 dt
      !real*8 hbar                         !planck constant in good units
      real*8 k_b                          !bolzmann in good units
      real*8 omega
      real*8 t
      real*8 normal
c      real*8, allocatable :: Htilde(:,:)
      real*8, allocatable :: r_memory(:)
      real*8 :: freq
      logical :: adqtb_restart
      character*120 gamma_restart_file
      character*3  numero
      real*8 Tseg_file, gamma_file, omegacut_file
c
c     perform dynamic allocation of some pointer arrays
c
      call alloc_shared_qtb

      if (allocated(rt)) deallocate (rt)
      allocate (rt(3,nloc,nseg))


      if (allocated(repartnoise)) deallocate (repartnoise)
      allocate(repartnoise(n))
      allocate(ntype(maxtyp))
      allocate(adqtb_type(1:n))
      adqtb_type=0
      ntype=0
      typemax=0
            
      do i=1,n
        k=type(i)
          if((k.gt.0).and.(atomic(i)/=0)) then
            ntype(k)=ntype(k)+1 
          endif
      enddo

      do k=1,maxtyp
        if((ntype(k).ne.0)) then
               typemax=typemax+1
               ntype(k)=typemax
        endif
      enddo

      do i=1,n
        k = type(i)
          if((k.gt.0).and.(atomic(i)/=0)) then
            adqtb_type(i)=ntype(k)
          endif
      enddo 

!
!     constant in Kcal/mol and ps
!
      !hbar=(planck*1.d11*avogadro)/(2*pi)
      k_b=boltzmann


!     Arguments
!
      domega = (2.*pi)/(3.*nseg*dt)
!
!     allocate Ht and noise arrays
!
c      allocate(noise(3,n,-nseg:2*nseg-1))
      if (allocated(Htilde)) deallocate (Htilde)
      allocate (Htilde(0:3*nseg-1,1:typemax))
!
c      allocate (rt(3,n,1:nseg))

!     create Htilde
!     
c      open(1,file='H_omega.out')
c      write(1,500)
c 500  format(8x, 'Frequency',9x, 2x, 'Htilde',6x) 
!


      Htilde=0.0
      do i=1,typemax
        do k=0, (3*nseg)/2
          omega=k*domega
            if  (k .eq. 0) then
              Htilde(k,i)=sqrt(k_b*kelvin)
c              write(1,'(2e20.5)') omega,Htilde(k,i)
            else if (noQTB) then
              Htilde(k,i)=sqrt(k_b*kelvin)
              Htilde(3*nseg-k,i)=sqrt(k_b*kelvin)
c              write(1,'(2e20.5)') omega, Htilde(k,i)
            else
              C_omega=(1-2*exp(-gamma*dt)*cos(omega*dt)+exp(-2*gamma*dt)
     &                )/((gamma**2+omega**2)*dt**2)
              theta_tilde=hbar*abs(omega)*(1.0/2.0+1.0
     &                   /(exp(hbar*abs(omega)/(k_b*kelvin))-1))
              f_omega=1.0/(1+exp((abs(omega)-omegacut)/omegasmear))
              Htilde(k,i)=sqrt(theta_tilde*f_omega*C_omega)
              Htilde(3*nseg-k,i)=sqrt(theta_tilde*f_omega*C_omega)
c                  write(1,'(2e20.5)') omega, Htilde(k,i)
c                  write(1,'(6e20.5)') omega, C_omega, theta
c     &                                _tilde, f_omega,Htilde(k)
c     &                  ,ns*dt*C_omega*theta_tilde*f_omega
            endif
          enddo
      enddo
!
c      close(1)

!     Initialisation for adQTB-r
      if (adaptive) then
        nad=int(omegacut/domega)

          if (allocated(vad)) deallocate (vad)
          if (allocated(fad)) deallocate (fad)
          allocate(vad(3,nloc,nseg))
          allocate(fad(3,nloc,nseg))

          
          allocate(mCvv_average_type(0:nad-1,1:typemax))
          allocate(Cvf_average_type(0:nad-1,1:typemax))
          allocate(dFDR_average_type(0:nad-1,1:typemax))

          allocate(gamma_type(0:nad-1,1:typemax))



          allocate(vad_piston(nseg))
          allocate(fad_piston(nseg))

          allocate(gamma_piston(0:nad-1))

          allocate(mCvv_average_piston(0:nad-1))
          allocate(Cvf_average_piston(0:nad-1))
          allocate(dFDR_average_piston(0:nad-1))


          do i=1,typemax
            do k=0,nad-1
              gamma_type(k,i)=gamma
            enddo
          enddo


          do k=0,nad-1
            gamma_piston(k)=gammapiston
          enddo
          
          allocate(corr_pot_ratio(0:nad-1))
               
          corr_pot_ratio(:)=0d0
            
          if(corr_pot) then
            open(315,file='corr_pot.out')
            read(315,*) Tseg_file, gamma_file,omegacut_file
            if((Tseg_file.ne.Tseg).or.(omegacut_file.ne.omegacut) 
     &          .or.(gamma_file.ne.gamma)) then
                write(*,*) 'ERROR, NOT THE SAME PARAMETERS WHILE USING
     & THE POTENTIAL CORRECTION'
                call fatal
            endif  
            do i=0,nad-1
               read(315,*) omega,corr_pot_ratio(i)
            enddo
          close(315)
          endif

    

          inquire (file='gamma_restart.out', exist=adqtb_restart)
          if (adqtb_restart) then
            skipseg=1
            open(98, file='gamma_restart.out')
              do j=0,nad-1
                read(98,*) freq,gamma_type(j,:)
              enddo
            close(98)
          endif
          call adHnoise(dt)
          else

      endif
      do i=1,nloc
        iglob = glob(i)
          do j=1,3
            do k=1,3*nseg
              noise(j,iglob,k)=normal()
            enddo
        enddo
      enddo

      call convolseg

      end      
c
c     subroutine alloc_shared_qtb : allocate shared memory pointers for qtb
c     parameter arrays
c
      subroutine alloc_shared_qtb
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use atoms
      use adqtb
      use domdec
      use qtb
      use mpi
      implicit none
      integer :: win,win2
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(3)
c
      if (associated(noise)) deallocate(noise)
c      if (associated(rt)) deallocate(rt)
c      if (associated(vad)) deallocate(vad)
c      if (associated(fad)) deallocate(fad)
c
c     noise
c
      arrayshape=(/3,n,3*nseg/)
      if (hostrank == 0) then
        windowsize = int(3*n*3*nseg,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, win, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,noise,arrayshape)
cc
cc     rt
cc
c      arrayshape=(/3,n,nseg/)
c      if (hostrank == 0) then
c        windowsize = int(3*n*nseg,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
c      else
c        windowsize = 0_MPI_ADDRESS_KIND
c      end if
c      disp_unit = 1
cc
cc    allocation
cc
c      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
c     $  hostcomm, baseptr, win, ierr)
c      if (hostrank /= 0) then
c        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
c     $  baseptr, ierr)
c      end if
cc
cc    association with fortran pointer
cc
c      CALL C_F_POINTER(baseptr,rt,arrayshape)

c      if (adaptive) then
cc
cc     vad
cc
c      arrayshape=(/3,n,nseg/)
c      if (hostrank == 0) then
c        windowsize = int(3*n*nseg,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
c      else
c        windowsize = 0_MPI_ADDRESS_KIND
c      end if
c      disp_unit = 1
cc
cc    allocation
cc
c      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
c     $  hostcomm, baseptr, win, ierr)
c      if (hostrank /= 0) then
c        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
c     $  baseptr, ierr)
c      end if
cc
cc    association with fortran pointer
cc
c      CALL C_F_POINTER(baseptr,vad,arrayshape)
cc
cc     fad
cc
c      arrayshape=(/3,n,nseg/)
c      if (hostrank == 0) then
c        windowsize = int(3*n*nseg,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
c      else
c        windowsize = 0_MPI_ADDRESS_KIND
c      end if
c      disp_unit = 1
cc
cc    allocation
cc
c      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
c     $  hostcomm, baseptr, win, ierr)
c      if (hostrank /= 0) then
c        CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit,
c     $  baseptr, ierr)
c      end if
cc
cc    association with fortran pointer
cc
c      CALL C_F_POINTER(baseptr,fad,arrayshape)
c      end if
c
      return
      end
