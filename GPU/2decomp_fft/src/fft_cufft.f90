!=======================================================================
! This (was) part of the 2DECOMP&FFT library
!
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil)
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
! This file is a modified file from 2DECOMP&FFT that uses cuda FFT for GPU
! it has been modified specifically for TinkerHP :
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!=======================================================================

#ifndef OVERWRITE
#warning(" 2DECOMP&FFT's 'OVERWRITE' mode always on in this file ")
#endif

! Define generic names for cudaFFT subroutines/variables for simple/double precision
#ifdef DOUBLE_PREC
#define CUFFT_X2X CUFFT_Z2Z
#define CUFFT_EXEC_X2X cufftExecZ2Z 
#else
#define CUFFT_X2X CUFFT_C2C
#define CUFFT_EXEC_X2X cufftExecC2C 
#endif


module decomp_2d_cufft
#ifdef _OPENACC
  use decomp_2d  ! 2D decomposition module
  use decomp_2d_fft
  use openacc
  use cufft      ! cuda_fft library module
  implicit none
  private

  integer :: cuplan_xy  ! Plan for 2D xy FFT
  integer :: cuplan_z   ! Plan for 1D z FFT
  integer :: cuplan_xyz ! Plan for 3D xyz FFT
  integer(4) :: ierr
  integer, parameter :: CUFFT_ERROR = 43, CUFFT_NOT_IMPLEMENTED = 44

  ! Declare generic interfaces to handle different inputs
  ! (~ interface from lib2decomp)
  public :: init_cufft_engine, finalize_cufft_engine, &
            get_rec_dir_queue, decomp_2d_cufft_3d, &
            decomp2d_cufftGetSize

  interface decomp_2d_cufft_3d
    module procedure cufft_3d_c2c
    ! Only c2c is used in TinkerHP
  end interface

contains

  ! Return a cudaFFT plan for multiple 1D c2c FFTs in Z direction
  subroutine cuc2c_1m_z_plan(plan1, decomp, isign)
    integer, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign
    integer(int_ptr_kind()) s_work

    ierr = cufftPlanMany(plan1, 1, decomp%zsz(3),    &
      decomp%zsz(3), decomp%zsz(1)*decomp%zsz(2), 1, &
      decomp%zsz(3), decomp%zsz(1)*decomp%zsz(2), 1, &
      CUFFT_X2X, decomp%zsz(1)*decomp%zsz(2))

    ierr = ierr + cufftGetSizeMany(plan1, 1, decomp%zsz(3),    &
      decomp%zsz(3), decomp%zsz(1)*decomp%zsz(2), 1, &
      decomp%zsz(3), decomp%zsz(1)*decomp%zsz(2), 1, &
      CUFFT_X2X, decomp%zsz(1)*decomp%zsz(2), s_work)
    cufft_worksize = cufft_worksize + s_work
    if (ierr /= 0) &
      call decomp_2d_abort(CUFFT_ERROR, "Cannot create plan for Z FFTs batch (cufftPlanMany)"  )
  end subroutine cuc2c_1m_z_plan

  ! Return a cudaFFT plan for multiple 2D c2c FFTs in xy direction
  subroutine cuc2c_1m_xy_plan(plan1, decomp, isign)
    implicit none
    integer(4), intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign
    integer :: n(2), nembed(2), stride, dist, batch, rank
    integer, pointer:: embed
    integer(int_ptr_kind()) s_work

    rank = 2
    n(1) = decomp%xsz(2); n(2)=decomp%xsz(1)
    nembed(:) = n(:)
    !nullify(embed)
    stride = 1
    dist = n(1)*n(2)
    batch = decomp%xsz(3)
    ierr = cufftPlanMany( plan1,   rank,    n, &
                         nembed, stride, dist, &
                         nembed, stride, dist, &
                      CUFFT_X2X, batch)
    ierr = ierr + cufftGetSizeMany( plan1,   rank,    n, &
                         nembed, stride, dist, &
                         nembed, stride, dist, &
                      CUFFT_X2X, batch, s_work)
    cufft_worksize = cufft_worksize + s_work
    if (ierr /= 0) &
      call decomp_2d_abort(CUFFT_ERROR, "Cannot create plan for XY FFTs batch (cufftPlanMany)"  )
  end subroutine

  ! Return a cudaFFT plan for 3D c2c FFTs in xyz direction
  subroutine cuc2c_1m_3D_plan(plan1, decomp, isign)
    implicit none
    integer(4), intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign
    integer :: n(3), nembed(3), stride, dist, batch, rank
    integer, pointer:: embed
    integer(int_ptr_kind()) s_work

    rank = 3
    n = decomp%xsz(:)
    nembed(:) = n(:)
    !nullify(embed)
    stride = 1
    dist   = n(1)*n(2)*n(3)
    batch  = 1

    ierr = cufftPlan3d(plan1,decomp%xsz(3),decomp%xsz(2),decomp%xsz(1) &
                      ,CUFFT_X2X)
    ierr = ierr + cufftGetSize3d(plan1,decomp%xsz(3),decomp%xsz(2) &
                      ,decomp%xsz(1) ,CUFFT_X2X, s_work)
!   ierr = cufftPlanMany( plan1,   rank,    n, &
!                        nembed, stride, dist, &
!                        nembed, stride, dist, &
!                     CUFFT_X2X, batch)
    cufft_worksize = cufft_worksize + s_work
    if (ierr /= 0) &
      call decomp_2d_abort(CUFFT_ERROR, "Cannot create plan for 3D FFTs batch (cufftPlan3D)"  )
  end subroutine

  ! This routine performs one-time initialisations for the FFT engine
  subroutine init_cufft_engine
    implicit none
    integer :: ierr

    if (nrank==0) then
      write(*,*) ' '
      write(*,*) '***** Using the CUFFT engine *****'
      write(*,*) ' '
    end if

    if (format /= PHYSICAL_IN_X) &
      call decomp_2d_abort(CUFFT_NOT_IMPLEMENTED, "format /= PHYSICAL_IN_X unsupported with cufft")

    ierr = 0

    if (nproc.eq.1) then
       call cuc2c_1m_3D_plan(cuplan_xyz, ph, CUFFT_FORWARD )
       ierr = ierr + cufftSetStream(cuplan_xyz,acc_get_cuda_stream(rec_queue))
    else
       call cuc2c_1m_z_plan (cuplan_z, ph, CUFFT_FORWARD )
       call cuc2c_1m_xy_plan(cuplan_xy, ph, CUFFT_FORWARD )
       ierr = ierr + cufftSetStream(cuplan_z, acc_get_cuda_stream(rec_queue))
       ierr = ierr + cufftSetStream(cuplan_xy,acc_get_cuda_stream(rec_queue))
    end if

    if (ierr.ne.0) &
      call decomp_2d_abort(CUFFT_ERROR, "Assigning cuda stream to plan failed (cufftSetStream)"  )
  end subroutine init_cufft_engine

  ! This routine performs one-time finalisations for the FFT engine
  subroutine finalize_cufft_engine
    if (nproc.eq.1) then
       ierr = cufftDestroy(cuplan_xyz)
    else
       ierr = cufftDestroy(cuplan_xy)
       ierr = cufftDestroy(cuplan_z)
    end if
    if (ierr.ne.0) &
      call decomp_2d_abort(CUFFT_ERROR, "Destroy cufft plan failed (cufftDestroy)" )
  end subroutine finalize_cufft_engine

  ! Following routines calculate multiple one/two-dimensional FFTs to form
  ! the basis of three-dimensional FFTs.
  ! ------------------------------------------------

  ! c2c transform, multiple 1D FFTs in z direction
  subroutine cuc2c_1m_z(inout, isign, plan1)
    complex(mytype), dimension(:,:,:), intent(INOUT), contiguous :: inout
    integer, intent(IN) :: isign
    integer(4), intent(IN) :: plan1

    !$acc host_data use_device(inout)
    ierr = CUFFT_EXEC_X2X(plan1, inout, inout, isign)
    !$acc end host_data
    if (ierr.ne.0) &
      call decomp_2d_abort(CUFFT_ERROR, "Z FFT execution failed (cufftExecZ2Z/C2C)" )
  end subroutine cuc2c_1m_z


  ! c2c transform, multiple 2D FFTs in xy direction
  subroutine cuc2c_1m_xy(inout, isign, plan1)
    complex(mytype), dimension(:,:,:), intent(INOUT), contiguous :: inout
    integer, intent(IN) :: isign
    integer(4), intent(IN) :: plan1
    !$acc host_data use_device(inout)
    ierr = CUFFT_EXEC_X2X(plan1, inout, inout, isign)
    !$acc end host_data
    if (ierr.ne.0) &
      call decomp_2d_abort(CUFFT_ERROR, "XY FFT execution failed (cufftExecZ2Z/C2C)" )
  end subroutine

  ! c2c transform, 3D FFT
  subroutine cuc2c_1m_3D(in, out, isign, plan1)
    complex(mytype), dimension(:,:,:), intent(INOUT), contiguous :: in
    complex(mytype), dimension(:,:,:), intent(INOUT), contiguous :: out
    integer   , intent(IN) :: isign
    integer(4), intent(IN) :: plan1

    !$acc host_data use_device(in,out)
    ierr = CUFFT_EXEC_X2X(plan1, in, out, isign)
    !$acc end host_data
    if (ierr.ne.0) &
      call decomp_2d_abort(CUFFT_ERROR, "3D FFT execution failed (cufftExecZ2Z/C2C)" )
  end subroutine


  ! End of 1D/2D routines
  ! -------------------------------------------

  ! 3D FFT - complex to complex
  subroutine cufft_3d_c2c(in, out, isign)
    complex(mytype), dimension(:,:,:), intent(INOUT) :: in
    complex(mytype), dimension(:,:,:), intent(OUT) :: out
    integer, intent(IN) :: isign

    if(format/=PHYSICAL_IN_X) &
      call decomp_2d_abort(CUFFT_NOT_IMPLEMENTED, "format /= PHYSICAL_IN_X unsupported with cufft")
    if (dims(1)>1) &
      call decomp_2d_abort(CUFFT_NOT_IMPLEMENTED, "Unsupported : xy not contiguous")

    ! Serial Mode
    if (nproc.eq.1) then
       if (isign == DECOMP_2D_FFT_FORWARD .or. isign == DECOMP_2D_FFT_BACKWARD) then
          call cuc2c_1m_3D(in,out,isign,cuplan_xyz)
          return
       else
          call decomp_2d_abort(CUFFT_ERROR, "Invalid argument : isign should be DECOMP_2D_FFT_FORWARD/BACKWARD" )
       end if
    end if

    ! MPI Mode
    if (isign == DECOMP_2D_FFT_FORWARD ) then
      ! 2D FFT xy
      call cuc2c_1m_xy(in,isign,cuplan_xy)
      ! Communicate to get all Z on local proc
      call cutranspose_y_to_z(in,out,ph)
      ! 1D FFT z
      call cuc2c_1m_z(out,isign,cuplan_z)
    else if (isign == DECOMP_2D_FFT_BACKWARD ) then
      ! 1D FFT z
      call cuc2c_1m_z(in,isign,cuplan_z)
      ! Communicate to get all XY on local proc
      call cutranspose_z_to_y(in,out,ph)
      ! 2D FFT xy
      call cuc2c_1m_xy(out,isign,cuplan_xy)
    else
      call decomp_2d_abort(CUFFT_ERROR, "Invalid argument : isign should be DECOMP_2D_FFT_FORWARD/BACKWARD" )
    end if
    decomp2d_mpi_fcall=1
  end subroutine

  ! get direct and reciproqual queue from tinker
  subroutine get_rec_dir_queue(recip,direc)
    implicit none
    integer, intent(in) :: recip,direc
    dir_queue = direc
    rec_queue = recip
  end subroutine

  function decomp2d_cufftGetSize() result(ss)
  integer(int_ptr_kind()) ss
  ss = cufft_worksize
  end function

#endif
end module decomp_2d_cufft
