c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################
c     ##                                                   ##
c     ##  module tinMemory  --  Application Memory manager ##
c     ##                                                   ##
c     #######################################################
c
c     s_shmem   shared memory size
c     sd_shmem  device shared memory size
c     sd_ddmem  device duplicated memory size
c     s_prmem   private memory size
c     sd_prmem  device private memory size
c     s_driver  CUDA Driver memory size
c     s_sfWork  Scaling factor workspace size
c     s_tinWork Tinker-HP workSpace Estimation

#include "tinker_macro.h"

      module nvshmem
      use,intrinsic:: iso_c_binding
#ifdef _CUDA
      use cudafor
#endif
      integer(c_int),bind(C):: mype, npes
      integer(4) COMM_NVSHMEM
#ifdef USE_NVSHMEM
      interface
        subroutine nvshmem_init( mype, npes, mpi_comm )
     &             bind(C,name="nvshmem_init_f")
        import c_int
        integer(c_int) :: mype,npes,mpi_comm
        end subroutine
        subroutine nvshmem_finalize ()
     &             bind(C)
        end subroutine
      end interface

      interface
        type(c_devptr) function nvshmem_malloc ( size )
     &                 bind(C)
          import c_size_t
          integer(c_size_t),value:: size
        end function
        type(c_devptr) function nvshmem_calloc ( size )
     &                 bind(C,name="nvshmem_calloc_f")
          import c_size_t
          integer(c_size_t),value:: size
        end function
        subroutine nvshmem_free ( ptr ) bind(C)
          import c_devptr
          type(c_devptr),value :: ptr
        end subroutine

        type(c_devptr) function nvshmem_ptr( ptr,pe )
     &                 bind(C)
        import c_int,c_devptr
        type(c_devptr)   ,value:: ptr
        integer(c_int),value:: pe
        end function
      end interface
#endif
      interface value_pe
        integer(4) function value_pe_i4( input )
        integer(4),intent(in)::input
        end function
        integer(8) function value_pe_i8( input )
        integer(8),intent(in)::input
        end function
      end interface
      contains

      integer(4) function value_pe_i4( input ) result(valu)
      integer(4),intent(in)::input
      valu = ceiling(1.0*input/npes)
      end function

      integer(8) function value_pe_i8( input ) result(valu)
      integer(8),intent(in)::input
      valu = ceiling(1.0*input/npes,kind=8)
      end function

      subroutine init_nvshmem( mpi_comm )
      implicit none
      integer(4),intent(in):: mpi_comm
      npes = 1
      mype = 0
      COMM_NVSHMEM = mpi_comm
#ifdef USE_NVSHMEM
      call nvshmem_init(mype,npes,COMM_NVSHMEM)
#endif
      end subroutine

      end module 

      module tinMemory
#ifdef USE_NVSHMEM
      use tinTypes,only: r2dDPC=>  Real2dDevPointerContainer
     &            ,      r3dDPC=>  Real3dDevPointerContainer
     &            ,        rDPC=>    RealDevPointerContainer
     &            ,       c8DPC=>   Char8DevPointerContainer
     &            ,      i2dDPC=>   Int2dDevPointerContainer
     &            ,        iDPC=>     IntDevPointerContainer
     &            ,        lDPC=> LogicalDevPointerContainer
     &            ,       l1DPC=>Logical1DevPointerContainer
#endif

      implicit none
      integer(4),parameter      :: mipk=int_ptr_kind()
      integer(int_ptr_kind()) s_shmem,sd_shmem,sd_ddmem
      integer(int_ptr_kind()) s_prmem,sd_prmem

      integer(int_ptr_kind()) s_cufft,s_curand,s_nvshmem
      integer(int_ptr_kind()) s_driver,s_sfwork,s_tinWork

      integer(int_ptr_kind()) szoi1,szoi,szoi8,szoip
      integer(int_ptr_kind()) szor4,szor8,szoTp,szoRp
      integer(int_ptr_kind()) szol1,szol
      integer(int_ptr_kind()) szoc8
      real(int_ptr_kind()),protected::Mio,Mo

      integer(1),private:: i1
      integer(8),private:: i8
      integer(4)   ,private:: i4
      real(t_p) ,private:: real_tp
      real(r_p) ,private:: real_rp
      real(4)   ,private:: r4
      real(8)   ,private:: r8
      logical   ,private:: lgc
      logical(1),private:: lgc1
      character*8,private:: char8

      enum, bind(C)
        enumerator :: memhost=00
        enumerator :: memacc
        enumerator :: memcuda
        enumerator :: memnvsh
        enumerator :: memfree
        enumerator :: memrealloc
        enumerator :: memsmartrealloc
      end enum

      integer:: debMem=0

      integer(4) mhostonly, mhostacc, mhostnvsh, mhostaccnvsh
     &      , mnvshonly, macconly, mfhostacc

      ! Extra allocation for parallel run
      logical   extra_alloc
      integer(4)  ,protected :: s_alloc
      real(t_p),parameter :: mem_inc=0.10

      integer(4),protected :: n_realloc_i
      integer(4),protected :: n_realloc_r

      parameter(
     &   szoi         = sizeof(i4)                 ,
     &   szoi1        = sizeof(i1)                 ,
     &   szoi8        = sizeof(i8)                 ,
     &   szoip        = sizeof(szoi)               ,
     &   szoTp        = sizeof(real_tp)            ,
     &   szoRp        = sizeof(real_rp)            ,
     &   szor4        = sizeof(r4)                 ,
     &   szor8        = sizeof(r8)                 ,
     &   szol1        = sizeof(lgc1)               ,
     &   szol         = sizeof(lgc)                ,
     &   szoc8        = sizeof(char8)
     &   )
      parameter(
     &   mhostonly    =     2**(memhost)           ,
     &   macconly     =     2**(memacc)            ,
     &   mnvshonly    =     2**(memnvsh)           ,
     &   mhostacc     = ior(2**(memacc) ,mhostonly),
     &   mhostnvsh    = ior(2**(memnvsh),mhostonly),
     &   mhostaccnvsh = ior(2**(memacc) ,mhostnvsh),
     &   mfhostacc    = ior(2**(memfree),mhostacc ),
     &   Mio          = 1024*1024                  ,
     &   Mo           = 1d6
     &   )

#ifdef USE_NVSHMEM_CUDA
      ! One pointer element of each Pointer container to serve
      ! as a temporary storage
      ! It will mainly be use to reassign device pointers from shared
      ! memory
      type(  iDPC),device,pointer::   iDPC_expl(:)
      type(  rDPC),device,pointer::   rDPC_expl(:)
      type(i2dDPC),device,pointer:: i2dDPC_expl(:)
      type(r2dDPC),device,pointer:: r2dDPC_expl(:)
#endif

      ! Memory request procedure for shared memory spaces
      interface shmem_request
        module subroutine shmem_int_req
     &          (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &           nshArray,d_nshArray,
#endif
     &           config)
          integer(4),pointer:: shArray(:)
          integer(4) request_shape(:)
          integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
          type(iDPC),   allocatable,optional:: nshArray(:)
          type(iDPC),device,pointer,optional:: d_nshArray(:)
#endif
          integer(4),optional::config
        end subroutine
        module subroutine shmem_int_req2
     &          (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &           nshArray,d_nshArray,
#endif
     &           config)
          integer(4),pointer:: shArray(:,:)
          integer(4) request_shape(:)
          integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
          type(i2dDPC),   allocatable,optional:: nshArray(:)
          type(i2dDPC),device,pointer,optional:: d_nshArray(:)
#endif
          integer(4),optional::config
        end subroutine
        module subroutine shmem_int1_req
     &    (shArray,winarray,request_shape,config)
          integer(1),pointer:: shArray(:)
          integer(4) request_shape(:)
          integer(4) winarray
          integer(4),optional::config
        end subroutine

        module subroutine shmem_logical_req
     &           (shArray,winarray,request_shape,
#ifdef USE_NVSHMEM_CUDA
     &           nshArray,d_nshArray,
#endif  
     &           config,start)
          implicit none
          logical,pointer:: shArray(:)
          integer(4) request_shape(:)
          integer(4) winarray
#ifdef USE_NVSHMEM_CUDA
          type(lDPC),   allocatable,optional::   nshArray(:)
          type(lDPC),device,pointer,optional:: d_nshArray(:)
#endif  
          integer(4),optional::config,start
        end subroutine
        module subroutine shmem_logical1_req
     &            (shArray,winarray,request_shape,
#ifdef USE_NVSHMEM_CUDA
     &             nshArray,d_nshArray,
#endif   
     &             config)
          implicit none
          logical(1),pointer:: shArray(:)
          integer(4) request_shape(:)
          integer(4) winarray
#ifdef USE_NVSHMEM_CUDA
          type(l1DPC)   ,allocatable,optional::   nshArray(:)
          type(l1DPC),device,pointer,optional:: d_nshArray(:)
#endif
          integer(4),optional::config
        end subroutine
        module subroutine shmem_char8_req
     &            (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &             nshArray,d_nshArray,
#endif
     &             config)
          implicit none
          character(8),pointer:: shArray(:)
          integer(4) request_shape(:)
          integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
          type(c8DPC),allocatable,optional:: nshArray(:)
          type(c8DPC),device,pointer,optional:: d_nshArray(:)
#endif
          integer(4)  ,optional::config
        end subroutine

        module subroutine shmem_real_req
     &          (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &           nshArray,d_nshArray,
#endif
     &           config)
          real(t_p),pointer:: shArray(:)
          integer(4) request_shape(:)
          integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
          type(rDPC),   allocatable,optional:: nshArray(:)
          type(rDPC),device,pointer,optional:: d_nshArray(:)
#endif
          integer(4),optional::config
        end subroutine
        module subroutine shmem_real_req2
     &          (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &           nshArray,d_nshArray,
#endif
     &           config)
          real(t_p),pointer:: shArray(:,:)
          integer(4) request_shape(:)
          integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
          type(r2dDPC),   allocatable,optional:: nshArray(:)
          type(r2dDPC),device,pointer,optional:: d_nshArray(:)
#endif
          integer(4),optional::config
        end subroutine
        module subroutine shmem_real_req3
     &          (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &           nshArray,d_nshArray,
#endif
     &           config)
          real(t_p),pointer:: shArray(:,:,:)
          integer(4) request_shape(:)
          integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
          type(r3dDPC),   allocatable,optional:: nshArray(:)
          type(r3dDPC),device,pointer,optional:: d_nshArray(:)
#endif
          integer(4),optional::config
        end subroutine
      end interface

      interface shmem_requestm
        module subroutine shmem_realm_req
     &          (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &           nshArray,d_nshArray,
#endif
     &           config)
          real(r_p),pointer:: shArray(:)
          integer(4) request_shape(:)
          integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
          type(rDPC),   allocatable,optional:: nshArray(:)
          type(rDPC),device,pointer,optional:: d_nshArray(:)
#endif
          integer(4),optional::config
        end subroutine
      end interface

      interface
        module integer(4) function shmem_index(iglob,asize,locpe)
     &                 result(ind)
        implicit none
        integer(4),intent(in)::iglob
        integer(4),intent(in)::asize
        integer(4),intent(out)::locpe
        end function
      end interface

      interface
      module subroutine nvshmem_get_HeapSize
      end subroutine
      end interface

      ! Display amount of memory used by the application
      interface
        module subroutine print_memory_usage()
        end subroutine
      end interface

      ! reallocating procedures on private memory
      interface prmem_request
      module subroutine prmem_int_req( array, n, async, queue, config )
        integer(4), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
      end subroutine
      module subroutine prmem_int_req1( array, sz_array, n, 
     &                  async, queue, config )
        integer(4), allocatable, intent(inout) :: array(:)
        integer(8), intent(in) :: n
        integer(8), intent(inout) :: sz_array
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
      end subroutine
      module subroutine prmem_4int1a_req( array, n,
     &                  async, queue, config )
        integer(4), allocatable, intent(inout) :: array(:)
        integer(8), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
      end subroutine
      module subroutine prmem_pint_req( array, n, async, config )
        integer(4), pointer, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: config
      end subroutine
      module subroutine prmem_int1_req( array, n, async )
        integer(1), allocatable, intent(inout) :: array(:)
        integer(4)   , intent(in) :: n
        logical   , optional, intent(in) :: async
      end subroutine
      module subroutine prmem_int_req2( array, nl, nc, async, queue,
     &                  config, nlst, ncst, cdim )
        integer(4), allocatable, intent(inout) :: array(:,:)
        integer(4), intent(in) :: nc, nl
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
        integer(4), optional, intent(in) :: nlst, ncst
        logical, optional, intent(in) :: cdim
      end subroutine
      module subroutine prmem_int1_req2( array, nl, nc, async )
        integer(1), allocatable, intent(inout) :: array(:,:)
        integer(4), intent(in) :: nc, nl
        logical, optional, intent(in) :: async
      end subroutine
      module subroutine prmem_logi_req( array, n, async, queue, config )
        logical, allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
      end subroutine
      module subroutine prmem_real_req( array, n, async, queue, config,
     &                  nst )
        real(t_p), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue
        integer(4), optional, intent(in) :: config,nst
      end subroutine
      module subroutine prmem_preal_req( array, n, async )
        real(t_p), pointer, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
      end subroutine
      module subroutine prmem_real_req2( array, nl, nc, async, queue,
     &                  config, nlst, ncst )
        real(t_p), allocatable, intent(inout) :: array(:,:)
        integer(4)  , intent(in) :: nc, nl
        logical  , optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
        integer(4), optional, intent(in) :: nlst, ncst
      end subroutine
      module subroutine prmem_real_req3(array, nx,ny,nz, async, queue_,
     &                  config)
        real(t_p), allocatable, intent(inout) :: array(:,:,:)
        integer(4)  , intent(in) :: nx, ny, nz
        logical  , optional, intent(in) :: async
        integer(4)  , optional, intent(in) :: queue_, config
      end subroutine
      module subroutine prmem_real_req4(array, nx, ny, nz, nc, async)
        real(t_p), allocatable, intent(inout) :: array(:,:,:,:)
        integer(4)  , intent(in) :: nx, ny, nz, nc
        logical  , optional, intent(in) :: async
      end subroutine
      end interface

      !Behave as same as --prmem_request-- but operates on
      ! real(r_p) type 
      interface prmem_requestm
      module subroutine prmem_realm_req( array, n, async )
        real(r_p), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
      end subroutine
      module subroutine prmem_realm_req1( array, n, async )
        real(r_p), allocatable, intent(inout) :: array(:)
        integer(mipk), intent(in) :: n
        logical, optional, intent(in) :: async
      end subroutine
      module subroutine prmem_prealm_req( array, n, async )
        real(r_p), pointer, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
      end subroutine
      module subroutine prmem_realm_req2( array, nl, nc, async )
        real(r_p), allocatable, intent(inout) :: array(:,:)
        integer(4)  , intent(in) :: nc, nl
        logical  , optional, intent(in) :: async
      end subroutine
      module subroutine prmem_realpm_req2( array, nl, nc, async )
        real(r_p), pointer, intent(inout) :: array(:,:)
        integer(4)  , intent(in) :: nc, nl
        logical  , optional, intent(in) :: async
      end subroutine
      module subroutine prmem_int8_req( array, n, async, queue, config )
        integer(8), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
      end subroutine
      module subroutine prmem_int8_req1(array, n, async, queue, config )
        integer(8), allocatable, intent(inout) :: array(:)
        integer(8), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
      end subroutine
      module subroutine prmem_int8_req2( array, nl, nc, async, queue,
     &                  config, nlst, ncst )
        integer(8), allocatable, intent(inout) :: array(:,:)
        integer(4), intent(in) :: nc, nl
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
        integer(4), optional, intent(in) :: nlst, ncst
      end subroutine
      module subroutine prmem_8int2p_req( array, nl, nc, async, queue,
     &                  config, nlst, ncst )
        implicit none
        integer(8), pointer, intent(inout) :: array(:,:)
        integer(4)   , intent(in) :: nc, nl
        logical   , optional   , intent(in) :: async
        integer(4)   , optional   , intent(in) :: queue, config
        integer(4)   , optional   , intent(in) :: nlst, ncst
      end subroutine
      module subroutine prmem_realm_req3(array, nx,ny,nz, async, queue_,
     &                  config)
        real(r_p), allocatable, intent(inout) :: array(:,:,:)
        integer(4)  , intent(in) :: nx, ny, nz
        logical  , optional, intent(in) :: async
        integer(4)  , optional, intent(in) :: queue_, config
      end subroutine
      end interface

      !! move alloc routine
      interface prmem_mvrequest
      module subroutine prmem_int_mvreq( array, n, async,queue,config )
        implicit none
        integer(4), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
      end subroutine
      end interface

      interface prmemGetAllocSize
      integer(4) module function get_prmem_alloc_size_i4(n) result(res)
      integer,intent(in):: n
      end function
      end interface

      ! Update routine for shared memory
      interface shmem_update_device
        module subroutine shmem_update_device_int(src,n,async_q,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &             dst,nd,
#endif
     &             config)
        implicit none
        integer(4) ,intent(in):: n
        integer(4) ,intent(in):: src(*)
        integer(4) ,optional,intent(in):: async_q
#if defined(_CUDA) && defined(USE_NVSHMEM)
        integer(4) ,optional,intent(in):: nd
        integer(4),device,optional:: dst(*)
#endif
        integer(4) ,optional,intent(in):: config
        end subroutine
        module subroutine shmem_update_device_real(src,n,async_q,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &             dst,nd,
#endif
     &             config)
        implicit none
        integer(4) ,intent(in):: n
        real(t_p),intent(in):: src(*)
        integer(4) ,optional,intent(in):: async_q
#if defined(_CUDA) && defined(USE_NVSHMEM)
        integer(4) ,optional,intent(in):: nd
        real(t_p),device,optional:: dst(*)
#endif
        integer(4) ,optional,intent(in):: config
        end subroutine
      end interface

#ifdef _CUDA
      interface nvshmem_check_data
      module subroutine nvshmem_int_check( sol, dat, n, npe, config )
      integer(4) sol(*)
      integer(4),device::dat(*)
      integer(4) n,npe
      integer(4),optional::config
      end subroutine
      module subroutine nvshmem_real_check( sol, dat, n, npe, config )
      implicit none
      real(t_p) sol(*)
      real(t_p),device::dat(*)
      integer(4) n,npe
      integer(4),optional::config
      end subroutine
      end interface
#endif

      interface mem_get
      module procedure mem_get_int
      module procedure mem_get_real
      end interface

      contains

#ifdef USE_NVSHMEM_CUDA
      attributes(global) subroutine DPC_test (devptr,npes)
      implicit none
      type(i2dDPC),device::devptr(:)
      integer(4),value::npes
      integer(4) i
      print*,'test',size(devptr),size(devptr(1)%pel)
      if (npes.gt.1) print*,'test1',size(devptr(2)%pel)
c     do i = 1,3
c        print*,'t',i,devptr(1)%pel(1,i)
c     end do
      end subroutine
#endif

      subroutine init_memhandle()
        implicit none
        s_shmem  = 0
        s_prmem  = 0
        sd_shmem = 0
        sd_prmem = 0
        sd_ddmem = 0
#ifdef _OPENACC
        s_driver = 443*Mio
#else
        s_driver = 0
#endif
        s_cufft  = 0
        s_curand = 0
#ifdef USE_NVSHMEM_CUDA
        call nvshmem_get_HeapSize
#else
        s_nvshmem= 0
#endif
        s_sfWork = 0
        s_tinWork= 0
        n_realloc_i = 0
        n_realloc_r = 0
        extra_alloc = .false.
        s_alloc  = 0
      end subroutine

      subroutine mem_get_int(hmem,dmem)
      implicit none
      integer(int_ptr_kind()),intent(out):: hmem,dmem
      hmem =  s_prmem +  s_shmem
      dmem = sd_prmem + sd_shmem + sd_ddmem
      end subroutine
      subroutine mem_get_real(hmem,dmem)
      implicit none
      real(8),intent(out):: hmem,dmem
      hmem = ( s_prmem +  s_shmem )/Mio
      dmem = (sd_prmem + sd_shmem + sd_ddmem)/Mio
      end subroutine

      end module
