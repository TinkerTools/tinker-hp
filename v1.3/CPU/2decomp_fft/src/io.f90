!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This module provides parallel IO facilities for applications based on
! 2D decomposition.

module decomp_2d_io

  use decomp_2d
  use MPI
#ifdef T3PIO
  use t3pio
#endif

  implicit none

  private        ! Make everything private unless declared public

  public :: decomp_2d_write_one, decomp_2d_read_one, &
       decomp_2d_write_var, decomp_2d_read_var, &
       decomp_2d_write_scalar, decomp_2d_read_scalar, &
       decomp_2d_write_plane, decomp_2d_write_every

  ! Generic interface to handle multiple data types and decompositions

  interface decomp_2d_write_one
     module procedure write_one_real
     module procedure write_one_complex
  end interface decomp_2d_write_one

  interface decomp_2d_read_one
     module procedure read_one_real
     module procedure read_one_complex
  end interface decomp_2d_read_one

  interface decomp_2d_write_var
     module procedure write_var_real
     module procedure write_var_complex
  end interface decomp_2d_write_var

  interface decomp_2d_read_var
     module procedure read_var_real
     module procedure read_var_complex
  end interface decomp_2d_read_var

  interface decomp_2d_write_scalar
     module procedure write_scalar_real
     module procedure write_scalar_complex
     module procedure write_scalar_integer
  end interface decomp_2d_write_scalar

  interface decomp_2d_read_scalar
     module procedure read_scalar_real
     module procedure read_scalar_complex
     module procedure read_scalar_integer
  end interface decomp_2d_read_scalar

  interface decomp_2d_write_plane
     module procedure write_plane_3d_real
     module procedure write_plane_3d_complex
!     module procedure write_plane_2d
  end interface decomp_2d_write_plane

  interface decomp_2d_write_every
     module procedure write_every_real
     module procedure write_every_complex
  end interface decomp_2d_write_every
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to write a single 3D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_one_real(ipencil,var,filename,opt_decomp)
    
    implicit none
    
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type

    data_type = real_type

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

! This file contain common code to be included by subroutines 
! 'mpiio_write_one_...' in io.f90

    ! Using MPI-IO to write a distributed 3D array into a file

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if
    
    ! determine subarray parameters
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif

#ifdef T3PIO
    call MPI_INFO_CREATE(info, ierror)
    gs = ceiling(real(sizes(1),mytype)*real(sizes(2),mytype)* &
         real(sizes(3),mytype)/1024./1024.)
    call t3pio_set_info(MPI_COMM_WORLD, info, "./", ierror, &
         GLOBAL_SIZE=gs, factor=1)
#endif

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
#ifdef T3PIO
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, info, fh, ierror)
#else
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
#endif
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
#ifdef T3PIO
    call MPI_INFO_FREE(info,ierror)
#endif

    return
  end subroutine write_one_real


  subroutine write_one_complex(ipencil,var,filename,opt_decomp)
    
    implicit none
    
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type

    data_type = complex_type

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

! This file contain common code to be included by subroutines 
! 'mpiio_write_one_...' in io.f90

    ! Using MPI-IO to write a distributed 3D array into a file

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if
    
    ! determine subarray parameters
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif

#ifdef T3PIO
    call MPI_INFO_CREATE(info, ierror)
    gs = ceiling(real(sizes(1),mytype)*real(sizes(2),mytype)* &
         real(sizes(3),mytype)/1024./1024.)
    call t3pio_set_info(MPI_COMM_WORLD, info, "./", ierror, &
         GLOBAL_SIZE=gs, factor=1)
#endif

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
#ifdef T3PIO
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, info, fh, ierror)
#else
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
#endif
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
#ifdef T3PIO
    call MPI_INFO_FREE(info,ierror)
#endif
    
    return
  end subroutine write_one_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to read from a file a single 3D array
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_one_real(ipencil,var,filename,opt_decomp)
    
    implicit none
    
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(INOUT) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type

    data_type = real_type
    
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

! This file contain common code to be included by subroutines 
! 'mpiio_read_one_...' in io.f90

    ! Using MPI-IO to write a distributed 3D array into a file

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if
    
    ! determine subarray parameters
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine read_one_real


  subroutine read_one_complex(ipencil,var,filename,opt_decomp)
    
    implicit none
    
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(INOUT) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type
    
    data_type = complex_type

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

! This file contain common code to be included by subroutines 
! 'mpiio_read_one_...' in io.f90

    ! Using MPI-IO to write a distributed 3D array into a file

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if
    
    ! determine subarray parameters
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine read_one_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the writing
  !  operation to prepare the writing of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_var_real(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(IN) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = real_type

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

! This file contain common code to be included by subroutines 
! 'write_var_...' in io.f90

  ! Using MPI-IO to write a distributed 3D variable to a file. File 
  ! operations (open/close) need to be done in calling application. This
  ! allows multiple variables to be written to a single file. Together 
  ! with the corresponding read operation, this is the perfect solution
  ! for applications to perform restart/checkpointing.

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    ! Create file type and set file view
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for the next write operation
    disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    if (data_type == complex_type) then
       disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    end if

    return
  end subroutine write_var_real


  subroutine write_var_complex(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = complex_type

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

! This file contain common code to be included by subroutines 
! 'write_var_...' in io.f90

  ! Using MPI-IO to write a distributed 3D variable to a file. File 
  ! operations (open/close) need to be done in calling application. This
  ! allows multiple variables to be written to a single file. Together 
  ! with the corresponding read operation, this is the perfect solution
  ! for applications to perform restart/checkpointing.

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    ! Create file type and set file view
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for the next write operation
    disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    if (data_type == complex_type) then
       disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    end if

    return
  end subroutine write_var_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_var_real(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(INOUT) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = real_type

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

! This file contain common code to be included by subroutines 
! 'read_var_...' in io.f90

  ! Using MPI-IO to read a distributed 3D variable from a file. File 
  ! operations (open/close) need to be done in calling application. This 
  ! allows multiple variables to be read from a single file. Together 
  ! with the corresponding write operation, this is the perfect solution
  ! for applications to perform restart/checkpointing.

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    ! Create file type and set file view
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for the next read operation
    disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    if (data_type == complex_type) then
       disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    end if

    return
  end subroutine read_var_real


  subroutine read_var_complex(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(INOUT) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = complex_type

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

! This file contain common code to be included by subroutines 
! 'read_var_...' in io.f90

  ! Using MPI-IO to read a distributed 3D variable from a file. File 
  ! operations (open/close) need to be done in calling application. This 
  ! allows multiple variables to be read from a single file. Together 
  ! with the corresponding write operation, this is the perfect solution
  ! for applications to perform restart/checkpointing.

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    ! Create file type and set file view
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for the next read operation
    disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    if (data_type == complex_type) then
       disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    end if

    return
  end subroutine read_var_complex

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write scalar variables as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_scalar_real(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(IN) :: n              ! number of scalars
    real(mytype), dimension(n), &
         intent(IN) :: var                ! array of scalars

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         real_type,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n ! only one rank needs to write
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, real_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes

    return
  end subroutine write_scalar_real


  subroutine write_scalar_complex(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    complex(mytype), dimension(n), intent(IN) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,complex_type, &
         complex_type,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, complex_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes*2

    return
  end subroutine write_scalar_complex


  subroutine write_scalar_integer(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    integer, dimension(n), intent(IN) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, MPI_INTEGER, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_INTEGER,m,ierror)
    disp = disp + n*m

    return
  end subroutine write_scalar_integer


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read scalar variables as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_scalar_real(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(IN) :: n              ! number of scalars
    real(mytype), dimension(n), &
         intent(INOUT) :: var             ! array of scalars

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         real_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, real_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes

    return
  end subroutine read_scalar_real
  

  subroutine read_scalar_complex(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    complex(mytype), dimension(n), intent(INOUT) :: var

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,complex_type, &
         complex_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, complex_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes*2

    return
  end subroutine read_scalar_complex


  subroutine read_scalar_integer(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    integer, dimension(n), intent(INOUT) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, MPI_INTEGER, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_INTEGER,m,ierror)
    disp = disp + n*m

    return
  end subroutine read_scalar_integer


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 2D slice of the 3D data to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_plane_3d_real(ipencil,var,iplane,n,filename, &
       opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    real(mytype), allocatable, dimension(:,:,:) :: wk2d
    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, data_type

    data_type = real_type

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

! This file contain common code to be included by subroutines 
! 'mpiio_write_plane_3d_...' in io.f90

    ! It is much easier to implement if all mpi ranks participate I/O.
    ! Transpose the 3D data if necessary.

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    if (iplane==1) then
       allocate(wk(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
       if (ipencil==1) then
          wk = var
       else if (ipencil==2) then
          call transpose_y_to_x(var,wk,decomp)
       else if (ipencil==3) then
          allocate(wk2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          call transpose_z_to_y(var,wk2,decomp)
          call transpose_y_to_x(wk2,wk,decomp)
          deallocate(wk2)
       end if
       allocate(wk2d(1,decomp%xsz(2),decomp%xsz(3)))
       do k=1,decomp%xsz(3)
          do j=1,decomp%xsz(2)
             wk2d(1,j,k)=wk(n,j,k)
          end do
       end do
       sizes(1) = 1
       sizes(2) = decomp%ysz(2)
       sizes(3) = decomp%zsz(3)
       subsizes(1) = 1
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = 0
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1

    else if (iplane==2) then
       allocate(wk(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
       if (ipencil==1) then
          call transpose_x_to_y(var,wk,decomp)
       else if (ipencil==2) then
          wk = var
       else if (ipencil==3) then
          call transpose_z_to_y(var,wk,decomp)
       end if
       allocate(wk2d(decomp%ysz(1),1,decomp%ysz(3)))
       do k=1,decomp%ysz(3)
          do i=1,decomp%ysz(1)
             wk2d(i,1,k)=wk(i,n,k)
          end do
       end do
       sizes(1) = decomp%xsz(1)
       sizes(2) = 1
       sizes(3) = decomp%zsz(3)
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = 1
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = 0
       starts(3) = decomp%yst(3)-1

    else if (iplane==3) then
       allocate(wk(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))
       if (ipencil==1) then
          allocate(wk2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          call transpose_x_to_y(var,wk2,decomp)
          call transpose_y_to_z(wk2,wk,decomp)
          deallocate(wk2)
       else if (ipencil==2) then
          call transpose_y_to_z(var,wk,decomp)
       else if (ipencil==3) then
          wk = var
       end if
       allocate(wk2d(decomp%zsz(1),decomp%zsz(2),1))
       do j=1,decomp%zsz(2)
          do i=1,decomp%zsz(1) 
             wk2d(i,j,1)=wk(i,j,n)
          end do
       end do
       sizes(1) = decomp%xsz(1)
       sizes(2) = decomp%ysz(2)
       sizes(3) = 1
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = 1
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = 0
    end if

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, wk2d, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    deallocate(wk,wk2d)

    return
  end subroutine write_plane_3d_real


  subroutine write_plane_3d_complex(ipencil,var,iplane,n, &
       filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    complex(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    complex(mytype), allocatable, dimension(:,:,:) :: wk2d
    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, data_type

    data_type = complex_type

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

! This file contain common code to be included by subroutines 
! 'mpiio_write_plane_3d_...' in io.f90

    ! It is much easier to implement if all mpi ranks participate I/O.
    ! Transpose the 3D data if necessary.

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    if (iplane==1) then
       allocate(wk(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
       if (ipencil==1) then
          wk = var
       else if (ipencil==2) then
          call transpose_y_to_x(var,wk,decomp)
       else if (ipencil==3) then
          allocate(wk2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          call transpose_z_to_y(var,wk2,decomp)
          call transpose_y_to_x(wk2,wk,decomp)
          deallocate(wk2)
       end if
       allocate(wk2d(1,decomp%xsz(2),decomp%xsz(3)))
       do k=1,decomp%xsz(3)
          do j=1,decomp%xsz(2)
             wk2d(1,j,k)=wk(n,j,k)
          end do
       end do
       sizes(1) = 1
       sizes(2) = decomp%ysz(2)
       sizes(3) = decomp%zsz(3)
       subsizes(1) = 1
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = 0
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1

    else if (iplane==2) then
       allocate(wk(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
       if (ipencil==1) then
          call transpose_x_to_y(var,wk,decomp)
       else if (ipencil==2) then
          wk = var
       else if (ipencil==3) then
          call transpose_z_to_y(var,wk,decomp)
       end if
       allocate(wk2d(decomp%ysz(1),1,decomp%ysz(3)))
       do k=1,decomp%ysz(3)
          do i=1,decomp%ysz(1)
             wk2d(i,1,k)=wk(i,n,k)
          end do
       end do
       sizes(1) = decomp%xsz(1)
       sizes(2) = 1
       sizes(3) = decomp%zsz(3)
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = 1
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = 0
       starts(3) = decomp%yst(3)-1

    else if (iplane==3) then
       allocate(wk(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))
       if (ipencil==1) then
          allocate(wk2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          call transpose_x_to_y(var,wk2,decomp)
          call transpose_y_to_z(wk2,wk,decomp)
          deallocate(wk2)
       else if (ipencil==2) then
          call transpose_y_to_z(var,wk,decomp)
       else if (ipencil==3) then
          wk = var
       end if
       allocate(wk2d(decomp%zsz(1),decomp%zsz(2),1))
       do j=1,decomp%zsz(2)
          do i=1,decomp%zsz(1) 
             wk2d(i,j,1)=wk(i,j,n)
          end do
       end do
       sizes(1) = decomp%xsz(1)
       sizes(2) = decomp%ysz(2)
       sizes(3) = 1
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = 1
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = 0
    end if

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, wk2d, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    deallocate(wk,wk2d)

    return
  end subroutine write_plane_3d_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 2D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !************** TO DO ***************
!  subroutine write_plane_2d(ipencil,var,filename)
!    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
!    real(mytype), dimension(:,:), intent(IN) :: var ! 2D array
!    character(len=*), intent(IN) :: filename
!
!    if (ipencil==1) then
!       ! var should be defined as var(xsize(2)
!
!    else if (ipencil==2) then
!
!    else if (ipencil==3) then
!
!    end if
!
!    return
!  end subroutine write_plane_2d


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write 3D array data for every specified mesh point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_every_real(ipencil,var,iskip,jskip,kskip, &
       filename, from1)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iskip,jskip,kskip 
    character(len=*), intent(IN) :: filename
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
                                  ! .false. - save n,2n,3n...

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, key,color,newcomm, data_type
    integer, dimension(3) :: xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    data_type = real_type

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

! This file contain common code to be included by subroutines 
! 'write_every_...' in io.f90

  ! To write every few points of a 3D array to a file

    ! work out the distribution parameters, which may be different from 
    ! the default distribution used by the decomposition library
    !  For exmample if nx=17 and p_row=4
    !    distribution is: 4 4 4 5

    ! If writing from the 1st element
    !  If saving every 3 points, then 5 points to be saved (17/3)
    !    default distribution would be 1 1 1 2
    !    However, 1st block (1-4) contains the 3rd point
    !             2nd block (5-8) contains the 6th point
    !             3rd block (9-12) contains the 9th and 12th point
    !             4th block (13-17) contains then 15th point
    !    giving a 1 1 2 1 distribution
    !    So cannot use the base decomposition library for such IO

    ! If writing from the n-th element (n=?skip)
    !  If saving every 3 points, then 6 points to be saved
    !    However, 1st block (1-4) contains the 1st & 4th point
    !             2nd block (5-8) contains the 7th point
    !             3rd block (9-12) contains the 10th point
    !             4th block (13-17) contains then 12th & 15th point
    !    giving a 1 2 2 1 distribution

    skip(1)=iskip
    skip(2)=jskip
    skip(3)=kskip

    do i=1,3
       if (from1) then
          xst(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xst(i)=xst(i)+1
          xen(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xst(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xst(i)=xst(i)+1
          xen(i) = xend(i)/skip(i)
       end if
       xsz(i) = xen(i)-xst(i)+1
    end do
       
    do i=1,3
       if (from1) then
          yst(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) yst(i)=yst(i)+1
          yen(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          yst(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) yst(i)=yst(i)+1
          yen(i) = yend(i)/skip(i)
       end if
       ysz(i) = yen(i)-yst(i)+1
    end do

    do i=1,3
       if (from1) then
          zst(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zst(i)=zst(i)+1
          zen(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zst(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zst(i)=zst(i)+1
          zen(i) = zend(i)/skip(i)
       end if
       zsz(i) = zen(i)-zst(i)+1
    end do

    ! if 'skip' value is large it is possible that some ranks do not 
    ! contain any points to be written. Subarray constructor requires 
    ! nonzero size so it is not possible to use MPI_COMM_WORLD for IO.
    ! Create a sub communicator for this...
    color = 1
    key = 0  ! rank order doesn't matter
    if (ipencil==1) then
       if (xsz(1)==0 .or. xsz(2)==0 .or. xsz(3)==0) then
          color = 2
       end if
    else if (ipencil==2) then
       if (ysz(1)==0 .or. ysz(2)==0 .or. ysz(3)==0) then
          color = 2
       end if
    else if (ipencil==3) then
       if (zsz(1)==0 .or. zsz(2)==0 .or. zsz(3)==0) then
          color = 2
       end if
    end if
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,newcomm,ierror)

    if (color==1) then ! only ranks in this group do IO collectively
       
       ! generate subarray information
       sizes(1) = xsz(1)
       sizes(2) = ysz(2)
       sizes(3) = zsz(3)
       if (ipencil==1) then
          subsizes(1) = xsz(1)
          subsizes(2) = xsz(2)
          subsizes(3) = xsz(3)
          starts(1) = xst(1)-1
          starts(2) = xst(2)-1
          starts(3) = xst(3)-1
       else if (ipencil==2) then
          subsizes(1) = ysz(1)
          subsizes(2) = ysz(2)
          subsizes(3) = ysz(3)
          starts(1) = yst(1)-1
          starts(2) = yst(2)-1
          starts(3) = yst(3)-1
       else if (ipencil==3) then
          subsizes(1) = zsz(1)
          subsizes(2) = zsz(2)
          subsizes(3) = zsz(3)
          starts(1) = zst(1)-1
          starts(2) = zst(2)-1
          starts(3) = zst(3)-1
       end if
       
       ! copy data from original array
       ! needs a copy of original array in global coordinate 
       if (ipencil==1) then
          allocate(wk(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
          allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
          wk2=var
          if (from1) then
             do k=xst(3),xen(3)
                do j=xst(2),xen(2)
                   do i=xst(1),xen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=xst(3),xen(3)
                do j=xst(2),xen(2)
                   do i=xst(1),xen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if   
       else if (ipencil==2) then
          allocate(wk(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))
          allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
          wk2=var
          if (from1) then
             do k=yst(3),yen(3)
                do j=yst(2),yen(2)
                   do i=yst(1),yen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=yst(3),yen(3)
                do j=yst(2),yen(2)
                   do i=yst(1),yen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if
       else if (ipencil==3) then
          allocate(wk(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))
          allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
          wk2=var
          if (from1) then
             do k=zst(3),zen(3)
                do j=zst(2),zen(2)
                   do i=zst(1),zen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=zst(3),zen(3)
                do j=zst(2),zen(2)
                   do i=zst(1),zen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if
       end if
       deallocate(wk2)

       ! MPI-IO
       call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
            MPI_ORDER_FORTRAN, data_type, newtype, ierror)
       call MPI_TYPE_COMMIT(newtype,ierror)
       call MPI_FILE_OPEN(newcomm, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_VIEW(fh,disp,data_type, &
            newtype,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_WRITE_ALL(fh, wk, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            data_type, MPI_STATUS_IGNORE, ierror)
       call MPI_FILE_CLOSE(fh,ierror)
       call MPI_TYPE_FREE(newtype,ierror)

       deallocate(wk)

    end if ! color==1

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    
    return
  end subroutine write_every_real


  subroutine write_every_complex(ipencil,var,iskip,jskip,kskip, &
       filename, from1)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iskip,jskip,kskip 
    character(len=*), intent(IN) :: filename
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
                                  ! .false. - save n,2n,3n...

    complex(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, key,color,newcomm, data_type
    integer, dimension(3) :: xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    data_type = complex_type

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

! This file contain common code to be included by subroutines 
! 'write_every_...' in io.f90

  ! To write every few points of a 3D array to a file

    ! work out the distribution parameters, which may be different from 
    ! the default distribution used by the decomposition library
    !  For exmample if nx=17 and p_row=4
    !    distribution is: 4 4 4 5

    ! If writing from the 1st element
    !  If saving every 3 points, then 5 points to be saved (17/3)
    !    default distribution would be 1 1 1 2
    !    However, 1st block (1-4) contains the 3rd point
    !             2nd block (5-8) contains the 6th point
    !             3rd block (9-12) contains the 9th and 12th point
    !             4th block (13-17) contains then 15th point
    !    giving a 1 1 2 1 distribution
    !    So cannot use the base decomposition library for such IO

    ! If writing from the n-th element (n=?skip)
    !  If saving every 3 points, then 6 points to be saved
    !    However, 1st block (1-4) contains the 1st & 4th point
    !             2nd block (5-8) contains the 7th point
    !             3rd block (9-12) contains the 10th point
    !             4th block (13-17) contains then 12th & 15th point
    !    giving a 1 2 2 1 distribution

    skip(1)=iskip
    skip(2)=jskip
    skip(3)=kskip

    do i=1,3
       if (from1) then
          xst(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xst(i)=xst(i)+1
          xen(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xst(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xst(i)=xst(i)+1
          xen(i) = xend(i)/skip(i)
       end if
       xsz(i) = xen(i)-xst(i)+1
    end do
       
    do i=1,3
       if (from1) then
          yst(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) yst(i)=yst(i)+1
          yen(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          yst(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) yst(i)=yst(i)+1
          yen(i) = yend(i)/skip(i)
       end if
       ysz(i) = yen(i)-yst(i)+1
    end do

    do i=1,3
       if (from1) then
          zst(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zst(i)=zst(i)+1
          zen(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zst(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zst(i)=zst(i)+1
          zen(i) = zend(i)/skip(i)
       end if
       zsz(i) = zen(i)-zst(i)+1
    end do

    ! if 'skip' value is large it is possible that some ranks do not 
    ! contain any points to be written. Subarray constructor requires 
    ! nonzero size so it is not possible to use MPI_COMM_WORLD for IO.
    ! Create a sub communicator for this...
    color = 1
    key = 0  ! rank order doesn't matter
    if (ipencil==1) then
       if (xsz(1)==0 .or. xsz(2)==0 .or. xsz(3)==0) then
          color = 2
       end if
    else if (ipencil==2) then
       if (ysz(1)==0 .or. ysz(2)==0 .or. ysz(3)==0) then
          color = 2
       end if
    else if (ipencil==3) then
       if (zsz(1)==0 .or. zsz(2)==0 .or. zsz(3)==0) then
          color = 2
       end if
    end if
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,newcomm,ierror)

    if (color==1) then ! only ranks in this group do IO collectively
       
       ! generate subarray information
       sizes(1) = xsz(1)
       sizes(2) = ysz(2)
       sizes(3) = zsz(3)
       if (ipencil==1) then
          subsizes(1) = xsz(1)
          subsizes(2) = xsz(2)
          subsizes(3) = xsz(3)
          starts(1) = xst(1)-1
          starts(2) = xst(2)-1
          starts(3) = xst(3)-1
       else if (ipencil==2) then
          subsizes(1) = ysz(1)
          subsizes(2) = ysz(2)
          subsizes(3) = ysz(3)
          starts(1) = yst(1)-1
          starts(2) = yst(2)-1
          starts(3) = yst(3)-1
       else if (ipencil==3) then
          subsizes(1) = zsz(1)
          subsizes(2) = zsz(2)
          subsizes(3) = zsz(3)
          starts(1) = zst(1)-1
          starts(2) = zst(2)-1
          starts(3) = zst(3)-1
       end if
       
       ! copy data from original array
       ! needs a copy of original array in global coordinate 
       if (ipencil==1) then
          allocate(wk(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
          allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
          wk2=var
          if (from1) then
             do k=xst(3),xen(3)
                do j=xst(2),xen(2)
                   do i=xst(1),xen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=xst(3),xen(3)
                do j=xst(2),xen(2)
                   do i=xst(1),xen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if   
       else if (ipencil==2) then
          allocate(wk(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))
          allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
          wk2=var
          if (from1) then
             do k=yst(3),yen(3)
                do j=yst(2),yen(2)
                   do i=yst(1),yen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=yst(3),yen(3)
                do j=yst(2),yen(2)
                   do i=yst(1),yen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if
       else if (ipencil==3) then
          allocate(wk(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))
          allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
          wk2=var
          if (from1) then
             do k=zst(3),zen(3)
                do j=zst(2),zen(2)
                   do i=zst(1),zen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=zst(3),zen(3)
                do j=zst(2),zen(2)
                   do i=zst(1),zen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if
       end if
       deallocate(wk2)

       ! MPI-IO
       call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
            MPI_ORDER_FORTRAN, data_type, newtype, ierror)
       call MPI_TYPE_COMMIT(newtype,ierror)
       call MPI_FILE_OPEN(newcomm, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_VIEW(fh,disp,data_type, &
            newtype,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_WRITE_ALL(fh, wk, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            data_type, MPI_STATUS_IGNORE, ierror)
       call MPI_FILE_CLOSE(fh,ierror)
       call MPI_TYPE_FREE(newtype,ierror)

       deallocate(wk)

    end if ! color==1

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    
    return
  end subroutine write_every_complex


end module decomp_2d_io
