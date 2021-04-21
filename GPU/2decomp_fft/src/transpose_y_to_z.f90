!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains the routines that transpose data from Y to Z pencil

  subroutine transpose_y_to_z_real(src, dst, opt_decomp)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#ifdef SHM
    real(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%ROW_INFO%SND_P
    call mem_split_yz_real(src, s1, s2, s3, work1, dims(2), &
         decomp%y2dist, decomp)
#else
    call mem_split_yz_real(src, s1, s2, s3, work1_r, dims(2), &
         decomp%y2dist, decomp)
#endif

    ! define receive buffer
#ifdef SHM
    work2_p = decomp%ROW_INFO%RCV_P
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y2cnts_s, decomp%y2disp_s, &
            real_type, work2, decomp%z2cnts_s, decomp%z2disp_s, &
            real_type, decomp%ROW_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    if (decomp%even) then
       call MPI_ALLTOALL(work1_r, decomp%y2count, &
            real_type, dst, decomp%z2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
    else
       call MPI_ALLTOALL(work1_r, decomp%y2count, &
            real_type, work2_r, decomp%z2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
    end if
#else
    call MPI_ALLTOALLV(work1_r, decomp%y2cnts, decomp%y2disp, &
         real_type, dst, decomp%z2cnts, decomp%z2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    call mem_merge_yz_real(work2, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
#else
#ifdef EVEN
    if (.not. decomp%even) then
       call mem_merge_yz_real(work2_r, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
    end if
#else
    ! note the receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed
#endif
#endif
    
    return
  end subroutine transpose_y_to_z_real


#ifdef OCC
  subroutine transpose_y_to_z_real_start(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if
    
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)

    ! rearrange source array as send buffer
    call mem_split_yz_real(src, s1, s2, s3, sbuf, dims(2), &
         decomp%y2dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%y2count, real_type, &
         rbuf, decomp%z2count, real_type, &
         DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y2cnts, decomp%y2disp, real_type, &
         rbuf, decomp%z2cnts, decomp%z2disp, real_type, &
         DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_z_real_start


  subroutine transpose_y_to_z_real_wait(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    call NBC_WAIT(handle, ierror)

    dst = rbuf

    return
  end subroutine transpose_y_to_z_real_wait
#endif


  subroutine transpose_y_to_z_complex(src, dst, opt_decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    complex comp
    integer i,j,k

#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%ROW_INFO%SND_P_c
    call mem_split_yz_complex(src, s1, s2, s3, work1, dims(2), &
         decomp%y2dist, decomp)
#else
    call mem_split_yz_complex(src, s1, s2, s3, work1_c, dims(2), &
         decomp%y2dist, decomp)
#endif
    
    ! define receive buffer
#ifdef SHM
    work2_p = decomp%ROW_INFO%RCV_P_c
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y2cnts_s, decomp%y2disp_s, &
            complex_type, work2, decomp%z2cnts_s, decomp%z2disp_s, &
            complex_type, decomp%ROW_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    if (decomp%even) then
       call MPI_ALLTOALL(work1_c, decomp%y2count, &
            complex_type, dst, decomp%z2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    else
       call MPI_ALLTOALL(work1_c, decomp%y2count, &
            complex_type, work2_c, decomp%z2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    end if
#else
    call MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp, &
         complex_type, dst, decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    call mem_merge_yz_complex(work2, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
#else
#ifdef EVEN
    if (.not. decomp%even) then
       call mem_merge_yz_complex(work2_c, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
    end if
#else
    ! note the receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed
#endif
#endif

    return
  end subroutine transpose_y_to_z_complex


#ifdef _OPENACC
  subroutine cutranspose_y_to_z_complex(src, dst, opt_decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN), contiguous  :: src
    complex(mytype), dimension(:,:,:), intent(OUT), contiguous :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    real(mytype) :: cre(3)

#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif

    integer :: s1,s2,s3,d1,d2,d3,i,j,k
    integer :: ierror
    !print*,'trans_yz'

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%ROW_INFO%SND_P_c
    call cumem_split_yz_complex(src, s1, s2, s3, work1, size(work1), dims(2), &
         decomp%y2dist, decomp)
#else
    call cumem_split_yz_complex(src, s1, s2, s3, work1_c, size(work1_c), dims(2), &
         decomp%y2dist, decomp)
#endif

    ! define receive buffer
#ifdef SHM
    work2_p = decomp%ROW_INFO%RCV_P_c
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif

#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
!$acc host_data use_device(work1, work2)
       call MPI_ALLTOALLV(work1, decomp%y2cnts_s, decomp%y2disp_s, &
            complex_type, work2, decomp%z2cnts_s, decomp%z2disp_s, &
            complex_type, decomp%ROW_INFO%SMP_COMM, ierror)
!$acc end host_data
    end if
#else

#ifdef EVEN
!$acc host_data use_device(work1_c, work2_c,dst)
    if (decomp%even) then
       call MPI_ALLTOALL(work1_c, decomp%y2count, &
            complex_type, dst, decomp%z2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    else
       call MPI_ALLTOALL(work1_c, decomp%y2count, &
            complex_type, work2_c, decomp%z2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    end if
!$acc end host_data
#else
    if (nproc.ne.1) then
       call Decomp2d_MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp,complex_type, &
                                       dst, decomp%z2cnts, decomp%z2disp,complex_type, &
                                       DECOMP_2D_COMM_ROW, ierror)
!      !$acc wait(rec_queue)
!      !$acc host_data use_device(work1_c, dst)
!      call MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp, &
!           complex_type, dst, decomp%z2cnts, decomp%z2disp, &
!           complex_type, DECOMP_2D_COMM_ROW, ierror)
!      !$acc end host_data
    else
       s1 = decomp%y2cnts(0)
       !$acc parallel loop present(dst,work1_c) async(rec_queue)
       do i=1,s1
          dst(i,1,1)=work1_c(i)
       end do
    end if

#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    call cumem_merge_yz_complex(work2, size(work2), d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
#else
#ifdef EVEN
    if (.not. decomp%even) then
       call cumem_merge_yz_complex(work2_c, size(work2_c), d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
    end if
#else
    ! note the receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed
#endif
#endif

    return
  end subroutine cutranspose_y_to_z_complex
#endif

  subroutine make_copy(b_in,b_out,ncount)
  implicit none
  complex(mytype), intent(in) :: b_in(*)
  complex(mytype), intent(out):: b_out(*)
  integer, intent(in) :: ncount
  integer i

  !$acc parallel loop present(b_in,b_out)
  do i=1,ncount
     b_out(i) = b_in(i)
  end do

  end subroutine


#ifdef OCC
  subroutine transpose_y_to_z_complex_start(handle, src, dst, sbuf, &
       rbuf, opt_decomp)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if
    
    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)

    ! rearrange source array as send buffer
    call mem_split_yz_complex(src, s1, s2, s3, sbuf, dims(2), &
         decomp%y2dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%y2count, &
         complex_type, rbuf, decomp%z2count, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y2cnts, decomp%y2disp, &
         complex_type, rbuf,decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_z_complex_start


  subroutine transpose_y_to_z_complex_wait(handle, src, dst, sbuf, &
       rbuf, opt_decomp)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    call NBC_WAIT(handle, ierror)

    dst = rbuf

    return
  end subroutine transpose_y_to_z_complex_wait
#endif


  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(n1,n2,n3), intent(IN) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%y2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y2count + 1
#else
       pos = decomp%y2disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_yz_real


  subroutine mem_split_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%y2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y2count + 1
#else
       pos = decomp%y2disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_yz_complex

  subroutine cumem_split_yz_complex(in,n1,n2,n3,out,no,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3,no
    complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
    complex(mytype), dimension(no), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos
    !print*, 'cumem_split_yz',no

!$acc data present(in,out)
    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%y2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y2count + 1
#else
       pos = decomp%y2disp(m) + 1
#endif
#endif

!$acc parallel loop collapse(3) async(rec_queue)
       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos+(i-1)+(j-i1)*n1+(k-1)*n1*(i2-i1+1)) = in(i,j,k)
             end do
          end do
       end do
    end do
!$acc end data
    return
  end subroutine cumem_split_yz_complex

  subroutine mem_merge_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    real(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%z2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%z2count + 1
#else
       pos = decomp%z2disp(m) + 1
#endif
#endif

       do k=i1,i2
          do j=1,n2
             do i=1,n1
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yz_real


  subroutine mem_merge_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)
    
    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%z2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%z2count + 1
#else
       pos = decomp%z2disp(m) + 1
#endif
#endif

       do k=i1,i2
          do j=1,n2
             do i=1,n1
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yz_complex

  subroutine cumem_merge_yz_complex(in,no,n1,n2,n3,out,iproc,dist,decomp)
    
    implicit none
    
    integer, intent(IN) :: n1,n2,n3,no
    complex(mytype), dimension(no), intent(IN) :: in
    complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos
    !print*, 'cumem_merge_yz',n1,n2,n3

!$acc data present(in,out)
    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%z2disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%z2count + 1
#else
       pos = decomp%z2disp(m) + 1
#endif
#endif

!$acc parallel loop collapse(3) async(rec_queue)
       do k=i1,i2
          do j=1,n2
             do i=1,n1
                out(i,j,k) = in(pos+(i-1)+(j-1)*n1+(k-i1)*n1*n2)
             end do
          end do
       end do
    end do
!$acc end data

    return
  end subroutine cumem_merge_yz_complex
