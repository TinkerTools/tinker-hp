
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Tex_sizeas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine convoltseg-- Does a convolution at each Tseg  ##
!     ##                                                           ##
!     ###############################################################
!

      subroutine convolseg
      use adqtb
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
      use qtb
      use timestat
      use usage
      use units
      use mpi
      implicit none
      integer i,j,k,l,istep,ierr
      integer iglob
      real*8 dt
      real*8 normal
      double complex, allocatable :: s_in(:), s_out(:)
      double complex, allocatable :: s_out_inv(:) 
      integer*8 plan_f,plan_b,est,s


      allocate(s_in(3*nseg))   
      allocate(s_out(3*nseg))
      allocate(s_out_inv(3*nseg))   

      s=3*nseg
      est=1
      rt=0d0
      repartnoise = 0
      call dfftw_plan_dft_1d(plan_f,s,s_in,s_out,-1,est)
      call dfftw_plan_dft_1d(plan_b,s,s_out,s_out_inv,1,est)

      do i=1,nloc
            iglob=glob(i)
            repartnoise(iglob) = rank
            if(atomic(iglob)==0) CYCLE
            do j=1,3
                  noise(j,iglob,:)=cshift(noise(j,iglob,:),nseg)
                  do k=2*nseg+1,3*nseg
                        noise(j,iglob,k)=normal()
                  enddo
            enddo
      enddo
      
     

      do i=1,nloc
            iglob = glob(i)
            if(atomic(iglob)==0) CYCLE
            do j=1,3
                  s_in(:)=dcmplx(noise(j,iglob,:),0d0)
                  call dfftw_execute(plan_f,s_in,s_out) 
                  s_out(:)=s_out(:)*Htilde(:,adqtb_type(iglob))
                  call dfftw_execute(plan_b,s_out,s_out_inv)
                  rt(j,i,:)=real(s_out_inv(nseg+1:2*nseg))/(3*nseg)
            enddo
      enddo

      call MPI_BARRIER(COMM_BEAD,ierr)
c
c     communicate rank that generated qtb noise
c 
c      
      call MPI_ALLREDUCE(MPI_IN_PLACE,repartnoise,n,MPI_INT,MPI_SUM,
     $   COMM_BEAD,ierr) 


c      call dfftw_destroy(plan_f)
c      call dfftw_destroy(plan_b)
      end      

