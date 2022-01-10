      subroutine fftw(istep,dt)

c      use, intrinsic :: iso_c_binding

      use atmtyp
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
      use adqtb
      use timestat
      use usage
      use units
      use mpi
      implicit none
      integer i,j,k,l,ierr
      integer iglob,istep
      real*8 dt
      real*8,  allocatable :: mCvv_type(:,:), Cvf_type(:,:)
c      real*8, allocatable :: Cff_type(:,:)
      real*8, allocatable :: dFDR_type(:,:)
      double complex, allocatable :: s_in(:), s_out_v(:), s_out_r(:)
      integer*8 planv, planr,est,s
      character*3  numero
      character*120 Cvf_file
c      character*120 Cff_file
      character*120 mCvv_file
      character*120 deltaFDR_file
c      character*120 gamma_history
      INTERFACE
         function int_to_str(value) result(str)
           IMPLICIT NONE
           INTEGER, INTENT(in) :: value
           CHARACTER(:), ALLOCATABLE :: str
         end
      END INTERFACE

c     Allocation of the different arrays used for the fft
      allocate(s_in(3*nseg))   
      allocate(s_out_v(3*nseg))
      allocate(s_out_r(3*nseg))

c     Allocation of the different used for the adqtb
      allocate(mCvv_type(0:nad-1,1:typemax))
      allocate(Cvf_type(0:nad-1,1:typemax))
      allocate(dFDR_type(0:nad-1,1:typemax))
c      allocate(Cff_type(0:nad-1,1:typemax))
c

      s = 3*nseg
      est=1

      do i=1,typemax
            do j=0,nad-1
                  mCvv_type(j,i)=0d0
                  Cvf_type(j,i)=0d0
c                  Cff_type(j,i)=0d0
                  dFDR_type(j,i)=0d0
            enddo
      enddo
c
      if (compteur .le. startsavespec) then
            mCvv_average_type=0d0
            Cvf_average_type=0d0
            dFDR_average_type=0d0
c            Cff_average_type=0d0
      endif

      Cvf_file='Cvf.out'

      mCvv_file='mCvv.out'

      DeltaFDR_file='DeltaFDT.out'

      open(68, file=Cvf_file)
c      open(108, file=Cff_file)
      open(78, file=mCvv_file)
      open(88, file=DeltaFDR_file)
      open(70, file='gamma_history.out',access='append')
c      open(98,file='gamma_restart.out')
c
      write(68,499,advance='no') 'Frequency' 
c      write(108,499,advance='no') 'Frequency' 
      write(78,499,advance='no') 'Frequency'
      write(88,499,advance='no') 'Frequency'
      write(70,499,advance='no') 'Frequency'
c      write(98,499,advance='no') 'Frequency'
  499 format(6x,A)
            do i=1,typemax
                  write(68,500,advance='no') name(i)
c                  write(108,500,advance='no') name(i)
                  write(78,500,advance='no') name(i)
                  write(88,500,advance='no') name(i)
                  write(70,500,advance='no') name(i)
c                 write(98,500,advance='no') name(i)
  500             format (6x,'Type_',A)
            enddo
      close(70)
      close(68)
      close(78)
      close(88)
      close(98)
      open(68, file=Cvf_file,access='append')
c      open(108, file=Cff_file)
      open(78, file=mCvv_file,access='append')
      open(88, file=DeltaFDR_file,access='append')
      open(70, file='gamma_history.out',access='append')
      open(98, file='gamma_restart.out')

      call dfftw_plan_dft_1d(planv,s,s_in,s_out_v,1,est)
      call dfftw_plan_dft_1d(planr,s,s_in,s_out_r,1,est)
      do i=1,nloc
         iglob=glob(i)
         if (use(iglob) .AND. (atomic(iglob) /= 0)) then
            do j=1,3
                  s_in(:)=dcmplx(0d0,0d0)
                  do k = 1, nseg
                        s_in(k)=dcmplx(vad(j,i,k),0d0)
                  end do
                  call dfftw_execute(planv,s_in,s_out_v) 
                   s_in(:)=dcmplx(0d0,0d0)
                  do k = 1, nseg
                        s_in(k)=dcmplx(fad(j,i,k),0d0)
                  end do
                  call dfftw_execute(planr,s_in,s_out_r)
                  mCvv_type(:,adqtb_type(iglob))
     &                =mCvv_type(:, adqtb_type(iglob))
     &                   +1.0/(3.0*n) *mass(iglob)
     &                   *abs(s_out_v(1:nad))**2/nseg
                  Cvf_type(:,adqtb_type(iglob))
     &                =Cvf_type(:,adqtb_type(iglob))
     &                   +1.0/(3.0*n)*real(s_out_v(1:nad)
     &                   *conjg(s_out_r(1:nad)/nseg))
             enddo
         endif
      enddo
c
c     the master gets the sums of Cvv and Cvf contributions
c
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,mCvv_type,typemax*nad,MPI_REAL8
     &,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
      else
        call mPI_REDUCE(mCvv_type,mCvv_type,typemax*nad,MPI_REAL8
     &,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
      end if
      if (rank.eq.0) then
        call MPI_REDUCE(MPI_IN_PLACE,Cvf_type,typemax*nad,MPI_REAL8
     &,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
      else
        call MPI_REDUCE(Cvf_type,Cvf_type,typemax*nad,MPI_REAL8
     &,MPI_SUM,0,
     $     MPI_COMM_WORLD,ierr)
      end if
c      if (rank.eq.0) then
c        call MPI_REDUCE(MPI_IN_PLACE,Cff_type,typemax*nad,MPI_REAL8
c     &,MPI_SUM,0,
c     $     MPI_COMM_WORLD,ierr)
c      else
c        call MPI_REDUCE(Cff_type,Cff_type,typemax*nad,MPI_REAL8
c     &,MPI_SUM,0,
c     $     MPI_COMM_WORLD,ierr)
c      end if

      if (rank.eq.0) then 

        do i=1,typemax
         do j=0,nad-1 
           dFDR_type(:,i)=mCvv_type(:,i)
     &                              *gamma_type(:,i)
     &                              -Cvf_type(:,i)
         enddo
        enddo
c
c        if (compteur .lt. 100) then
c             dFDR_type(80:120)=0.
c        endif
         do i=1,typemax
            do j=0,nad-1
                  gamma_type(j,i)=max(
     &                  gamma_type(j,i)-A_gamma*dt
     &                  *nseg*dFDR_type(j,i)*gamma
     &                  /sqrt(sum(abs(dFDR_type(:,i))**2))
     &                  ,0.)
            enddo
        enddo  
        mCvv_average_type=mCvv_average_type+mCvv_type
        Cvf_average_type=Cvf_average_type+Cvf_type
        dFDR_average_type=dFdr_average_type+dFdr_type

        do l=0,nad-1
              write(68,'('//int_to_str(1+typemax)//'E16.8)') 
     &              domega*l,Cvf_average_type(l,:)
     &              /max(compteur-startsavespec+1,1)
              write(78,'('//int_to_str(1+typemax)//'E16.8)')
     &               domega*l,mCvv_average_type(l,:)
     &              /max(compteur-startsavespec+1,1)
              write(88,'('//int_to_str(1+typemax)//'E16.8)')
     &               domega*l,dFDR_average_type(l,:)
     &              /max(compteur-startsavespec+1,1)
         write(70,'('//int_to_str(1+typemax)//'E16.8)') 
     &              domega*l,gamma_type(l,:)
         write(98,'('//int_to_str(1+typemax)//'E16.8)') 
     &              domega*l,gamma_type(l,:)
        enddo
      end if
c
c     the master sends the adapted gamma_r to everybody
c
      call MPI_BCAST(gamma_type,typemax*nad,MPI_REAL8,0,MPI_COMM_WORLD
     &,ierr)
      close(68)
      close(78)
      close(88)
      close (70)
      close(98)
c     close(108)

      call dfftw_destroy_plan(planv)
      call dfftw_destroy_plan(planr)
      
      end
      
      function int_to_str(value) result(str)
           IMPLICIT NONE
           INTEGER, INTENT(in) :: value
           CHARACTER(:), ALLOCATABLE :: str
           INTEGER :: n_char
           CHARACTER(10) :: n_char_char

           if(value==0) then
               n_char=1
           else
               n_char=int(log10(real(value)))+1
           endif
           !write(0,*) n_char
           allocate(character(n_char) :: str)
           write(n_char_char,'(i10.10)') n_char
          ! write(0,*) n_char_char
           write(str,'(i'//trim(n_char_char)//')') value
       end function int_to_str
