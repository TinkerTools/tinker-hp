c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################
c     ##                                                   ##
c     ##  module subMemory  --  Application Memory manager ##
c     ##                                                   ##
c     #######################################################
c
c     tinMemory submodule contains the implementation part of
c     of tinMemory
c
#include "tinker_macro.h"

      submodule(tinMemory) subMemory
      use atoms  ,only: natoms=>n
      use domdec ,only: rank,hostrank,hostcomm,nproc
     &           ,nbloc,nrec_recep
      use fft    ,only: n1mpimax,n2mpimax,n3mpimax
      use,intrinsic:: iso_c_binding
      use mpi
      use mpole  ,only: npolelocnl,npoleloc,npolerecloc
     &           ,npolebloc
      use nvshmem
#ifdef _OPENACC
      use openacc
#endif
      use utilgpu,only: mem_set,rec_stream,rec_queue
      use vdw    ,only: nvdwlocnl
      integer(4),private:: mv_realloc_i=0

      contains

      module subroutine shmem_int_req
     &                 (shArray,winarray,request_shape,
#ifdef USE_NVSHMEM_CUDA
     &                  nshArray,d_nshArray,
#endif
     &                  config)
      implicit none
      integer(4),pointer:: shArray(:)
      integer(4) request_shape(:)
      integer(4) winarray
#ifdef USE_NVSHMEM_CUDA
      type(iDPC),   allocatable,optional::   nshArray(:)
      type(iDPC),device,pointer,optional:: d_nshArray(:)
#endif
      integer(4),optional::config

      integer(4) request_size,configure,i,istat
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4) :: disp_unit,ierr
      type(c_ptr) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
      integer(c_size_t) sh_size
#endif
 13   format('shmem_int_req::deallocate array',I4)
 14   format('shmem_int_req::allocate array',I7,I4)
      if (present(config)) then
         configure = config
      else
         configure = mhostonly
      end if

c
c     Deallocation
c
      if (associated(shArray)) then
         if (debMem.and.hostrank.eq.0) print 13, configure
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - size(shArray)*szoi
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - size(shArray)*szoi
            CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            CALL MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if
#ifdef USE_NVSHMEM_CUDA
      ! Deallocate nvshmem shared memory data structure
      if (present(nshArray).and.btest(configure,memnvsh)) then
         if (allocated(nshArray)) then
            do i = 0,npes-1
               if (i.eq.mype) then
                  sd_shmem = sd_shmem - size(nshArray(mype)%pel)*szoi
                  peptr = c_devloc(nshArray(mype)%pel)
                  call nvshmem_free(peptr)
               else
                  nullify(nshArray(i)%pel)
               end if
            end do
            deallocate(nshArray)
            deallocate(d_nshArray)
         end if
      end if
#endif

      ! Fetch for size per process
      request_size = request_shape(1)
      if (request_size.eq.0) return
      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                4_MPI_ADDRESS_KIND
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1

c
c     allocation
c
      if (debMem.and.hostrank.eq.0) print 14,request_size,rank
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape)
        s_shmem = s_shmem + size(shArray)*szoi
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + size(shArray)*szoi
      end if

#ifdef USE_NVSHMEM_CUDA
      if (present(nshArray).and.btest(configure,memnvsh)) then
         allocate(nshArray(0:npes-1),stat=istat)
         call ERROR_CHECK(istat,__FILE__,__LINE__,"Allocate")
         allocate(d_nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,"device Allocate")
         sh_size = value_pe(request_size)
         cuptr = nvshmem_malloc( int(sh_size*sizeof(sh_size),c_size_t) )
         call c_f_pointer(cuptr, nshArray(mype)%pel, [sh_size])
         sd_shmem = sd_shmem + size(nshArray(mype)%pel)*szoi
         do i = 0, npes-1
            if (i.ne.npes) then
               peptr = nvshmem_ptr( cuptr,i )
               call c_f_pointer(peptr, nshArray(i)%pel, [sh_size])
            end if
         end do
         istat = cudamemcpy(c_devloc(d_nshArray(0)),c_loc(nshArray(0)),
     &                      sizeof(nshArray))
         CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,"cudamemcpy")
      end if
#endif
      end

      module subroutine shmem_int_req2
     &                 (shArray,winarray,request_shape,
#ifdef USE_NVSHMEM_CUDA
     &                  nshArray,d_nshArray,
#endif
     &                  config)
      implicit none
      integer(4),pointer:: shArray(:,:)
      integer(4) request_shape(:)
      integer(4) winarray
#ifdef USE_NVSHMEM_CUDA
      type(i2dDPC)   ,allocatable,optional::   nshArray(:)
      type(i2dDPC),device,pointer,optional:: d_nshArray(:)
#endif
      integer(4),optional::config

      integer(mipk) request_size
      integer(4) configure,i,istat
      integer(4) stride,nelem
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4) :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
      integer(c_size_t) sh_size
#endif

 13   format('shmem_int_req2::deallocate array',I4)
 14   format('shmem_int_req2::allocate array',2I8,I10,I4)
      if (present(config)) then
         configure = config
      else
         configure = mhostonly
      end if
c
c     Deallocation
c
      if (associated(shArray)) then
         if (debMem.and.hostrank.eq.0) print 13, configure
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - (size(shArray,1)*szoi)*size(shArray,2)
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - (size(shArray,1)*szoi)*size(shArray,2)
            CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            CALL MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if
#ifdef USE_NVSHMEM_CUDA
      ! Deallocate nvshmem shared memory data structure
      if (present(nshArray).and.btest(configure,memnvsh)) then
         if (allocated(nshArray)) then
            do i = 0,npes-1
               if (i.eq.mype) then
                  sd_shmem = sd_shmem - size(nshArray(mype)%pel)*szoi
                  peptr = c_devloc(nshArray(mype)%pel)
                  call nvshmem_free(peptr)
               else
                  nullify(nshArray(i)%pel)
               end if
            end do
            deallocate(nshArray)
            deallocate(d_nshArray)
         end if
      end if
#endif

      ! Fetch for size per process
      request_size = int(request_shape(1),mipk)*request_shape(2)
      if (request_size.eq.0) return

      stride = request_shape(1)
      nelem  = request_shape(2)
      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                4_MPI_ADDRESS_KIND
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1

c
c     allocation
c
      if (debMem.and.hostrank.eq.0) print 14,request_shape
     &                              ,request_size,rank
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape)
        s_shmem = s_shmem + request_size*szoi
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + request_size*szoi
      end if

#ifdef USE_NVSHMEM_CUDA
      if (present(nshArray).and.present(d_nshArray).and.
     &    btest(configure,memnvsh)) then

         allocate(nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,"Allocate")
         allocate(d_nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,"device Allocate")
         sh_size = value_pe(nelem)
         cuptr   = nvshmem_calloc( int(stride*sh_size*szoi,c_size_t) )
         call c_f_pointer(cuptr, nshArray(mype)%pel, [stride,sh_size])
         sd_shmem = sd_shmem + stride*sh_size*szoi
         do i = 0, npes-1
            if (i.ne.npes) then
               peptr = nvshmem_ptr( cuptr,i )
               call c_f_pointer(peptr,nshArray(i)%pel,[stride,sh_size])
            end if
         end do
         !d_nshArray = nshArray
         istat = cudamemcpy(c_devloc(d_nshArray(0)),c_loc(nshArray(0)),
     &                      sizeof(nshArray))
         CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,"cudamemcpy")
      end if
#endif
      end

      module subroutine shmem_int1_req
     &                 (shArray,winarray,request_shape,config)
      implicit none
      integer(1),pointer:: shArray(:)
      integer(4) request_shape(:)
      integer(4) winarray
      integer(4),optional::config

      integer(4) request_size,configure,i,istat
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4)     :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
      integer(c_size_t) sh_size
#endif

      if (present(config)) then
         configure = config
      else
         configure = mhostonly
      end if

c
c     Deallocation
c
      if (associated(shArray)) then
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - size(shArray)*szoi1
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - size(shArray)*szoi1
            call MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            call MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if

      ! Fetch for size per process
      request_size = request_shape(1)
      if (request_size.eq.0) return
      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                1_MPI_ADDRESS_KIND
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1

c
c     allocation
c
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape)
        s_shmem = s_shmem + size(shArray)*szoi1
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + size(shArray)*szoi1
      end if
      end

      module subroutine shmem_logical_req
     &                 (shArray,winarray,request_shape,
#ifdef USE_NVSHMEM_CUDA
     &                  nshArray,d_nshArray,
#endif
     &                  config,start)
      implicit none
      logical,pointer:: shArray(:)
      integer(4) request_shape(:)
      integer(4) winarray
#ifdef USE_NVSHMEM_CUDA
      type(lDPC)   ,allocatable,optional::   nshArray(:)
      type(lDPC),device,pointer,optional:: d_nshArray(:)
#endif
      integer(4),optional::config,start

      integer(4) request_size,configure,i,istat,start_
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4) :: disp_unit,ierr
      type(c_ptr) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
      integer(c_size_t) sh_size
#endif

      ! Default config
      configure = mhostonly
      start_    = 1
      ! Configurate options
      if (present(config)) configure = config
      if (present(start )) start_    = start
c
c     Deallocation
c
      if (associated(shArray)) then
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - size(shArray)*szol
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - size(shArray)*szol
            CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            CALL MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if
#if defined(_CUDA) && defined(USE_NVSHMEM)
      ! Deallocate nvshmem shared memory data structure
      if (present(nshArray).and.btest(configure,memnvsh)) then
         if (allocated(nshArray)) then
            do i = 0,npes-1
               if (i.eq.mype) then
                  sd_shmem = sd_shmem - size(nshArray(mype)%pel)*szol
                  peptr = c_devloc(nshArray(mype)%pel)
                  call nvshmem_free(peptr)
               else
                  nullify(nshArray(i)%pel)
               end if
            end do
            deallocate(nshArray)
            deallocate(d_nshArray)
         end if
      end if
#endif

      ! Fetch for size per process
      request_size = request_shape(1)
      if (request_size.eq.0) return
      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                4_MPI_ADDRESS_KIND
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1

c
c     allocation
c
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape)
        s_shmem = s_shmem + sizeof(shArray)

        ! Reshape pointer if necessary
        if (start_.ne.1)
     &  shArray(start_:start_+request_size-1) => shArray
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + sizeof(shArray)
      end if

#if defined(_CUDA) && defined(USE_NVSHMEM)
      if (present(nshArray).and.btest(configure,memnvsh)) then
         allocate(nshArray(0:npes-1),stat=istat)
         call ERROR_CHECK(istat,__FILE__,__LINE__,"Allocate")
         allocate(d_nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,"device Allocate")
         sh_size = value_pe(request_size)
         cuptr = nvshmem_malloc( int(sh_size*sizeof(sh_size),c_size_t) )
         call c_f_pointer(cuptr, nshArray(mype)%pel, [sh_size])
         sd_shmem = sd_shmem + sizeof(nshArray(mype)%pel)
         do i = 0, npes-1
            if (i.ne.npes) then
               peptr = nvshmem_ptr( cuptr,i )
               call c_f_pointer(peptr, nshArray(i)%pel, [sh_size])
            end if
         end do
         istat = cudamemcpy(c_devloc(d_nshArray(0)),c_loc(nshArray(0)),
     &                      sizeof(nshArray))
         CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,"cudamemcpy")
      end if
#endif
      end

      module subroutine shmem_logical1_req
     &                 (shArray,winarray,request_shape,
#ifdef USE_NVSHMEM_CUDA
     &                  nshArray,d_nshArray,
#endif
     &                  config)
      implicit none
      logical(1),pointer:: shArray(:)
      integer(4) request_shape(:)
      integer(4) winarray
#ifdef USE_NVSHMEM_CUDA
      type(l1DPC)   ,allocatable,optional::   nshArray(:)
      type(l1DPC),device,pointer,optional:: d_nshArray(:)
#endif
      integer(4),optional::config

      integer(4) request_size,configure,i,istat
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4) :: disp_unit,ierr
      type(c_ptr) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
      integer(c_size_t) sh_size
#endif

      if (present(config)) then
         configure = config
      else
         configure = mhostonly
      end if

c
c     Deallocation
c
      if (associated(shArray)) then
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - sizeof(shArray)
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - sizeof(shArray)
            CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            CALL MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if
#if defined(_CUDA) && defined(USE_NVSHMEM)
      ! Deallocate nvshmem shared memory data structure
      if (present(nshArray).and.btest(configure,memnvsh)) then
         if (allocated(nshArray)) then
            do i = 0,npes-1
               if (i.eq.mype) then
                  sd_shmem = sd_shmem - sizeof(nshArray(mype)%pel)
                  peptr = c_devloc(nshArray(mype)%pel)
                  call nvshmem_free(peptr)
               else
                  nullify(nshArray(i)%pel)
               end if
            end do
            deallocate(nshArray)
            deallocate(d_nshArray)
         end if
      end if
#endif

      ! Fetch for size per process
      request_size = request_shape(1)
      if (request_size.eq.0) return
      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                1_MPI_ADDRESS_KIND
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1

c
c     allocation
c
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape)
        s_shmem = s_shmem + sizeof(shArray)
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + sizeof(shArray)
      end if

#if defined(_CUDA) && defined(USE_NVSHMEM)
      if (present(nshArray).and.btest(configure,memnvsh)) then
         allocate(nshArray(0:npes-1),stat=istat)
         call ERROR_CHECK(istat,__FILE__,__LINE__,"Allocate")
         allocate(d_nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,"device Allocate")
         sh_size = value_pe(request_size)
         cuptr = nvshmem_malloc( int(sh_size,c_size_t) )
         call c_f_pointer(cuptr, nshArray(mype)%pel, [sh_size])
         sd_shmem = sd_shmem + sizeof(nshArray(mype)%pel)
         do i = 0, npes-1
            if (i.ne.npes) then
               peptr = nvshmem_ptr( cuptr,i )
               call c_f_pointer(peptr, nshArray(i)%pel, [sh_size])
            end if
         end do
         istat = cudamemcpy(c_devloc(d_nshArray(0)),c_loc(nshArray(0)),
     &                      sizeof(nshArray))
         CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,"cudamemcpy")
      end if
#endif
      end

      module subroutine shmem_char8_req
     &                 (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &                  nshArray,d_nshArray,
#endif
     &                  config)
      implicit none
      character(8),pointer:: shArray(:)
      integer(4) request_shape(:)
      integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
      type(c8DPC),allocatable,optional:: nshArray(:)
      type(c8DPC),device,pointer,optional:: d_nshArray(:)
#endif
      integer(4)  ,optional::config

      integer(4) request_size,configure,i,istat
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4) :: disp_unit,ierr
      character(8) adata
      TYPE(C_PTR) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
      integer(c_size_t) sh_size
#endif

      if (present(config)) then
         configure = config
      else
         configure = mhostonly
      end if
c
c     Deallocation
c
      if (associated(shArray)) then
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - size(shArray)*szoc8
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - size(shArray)*szoc8
            CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            CALL MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if
#if defined(_CUDA) && defined(USE_NVSHMEM)
      ! Deallocate nvshmem shared memory data structure
      if (present(nshArray).and.btest(configure,memnvsh)) then
         if (allocated(nshArray)) then
            do i = 0,npes-1
               if (i.eq.mype) then
                  sd_shmem = sd_shmem - size(nshArray(mype)%pel)*szoc8
                  peptr = c_devloc(nshArray(mype)%pel)
                  call nvshmem_free(peptr)
               else
                  nullify(nshArray(i)%pel)
               end if
            end do
            deallocate(nshArray)
            deallocate(d_nshArray)
         end if
      end if
#endif

      ! Fetch for size per process
      request_size = request_shape(1)
      if (request_size.eq.0) return
      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                8_MPI_ADDRESS_KIND
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c     allocation
c
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape)
        s_shmem = s_shmem + size(shArray)*szoc8
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + size(shArray)*szoc8
      end if

#if defined(USE_NVSHMEM) && defined(_CUDA)
      if (present(nshArray).and.present(d_nshArray).and.
     &    btest(configure,memnvsh)) then

         allocate(nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,"Allocate")
         allocate(d_nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,"device Allocate")
         sh_size = value_pe(request_size)
         cuptr   = nvshmem_calloc(
     &              sh_size*sizeof(adata) )
         call c_f_pointer(cuptr, nshArray(mype)%pel, [sh_size])
         sd_shmem = sd_shmem + sizeof(nshArray(mype)%pel)
         do i = 0, npes-1
            if (i.ne.npes) then
               peptr = nvshmem_ptr( cuptr,i )
               call c_f_pointer(peptr,nshArray(i)%pel,[sh_size])
            end if
         end do
         !d_nshArray = nshArray
         istat = cudamemcpy(c_devloc(d_nshArray(0)),c_loc(nshArray(0)),
     &                      sizeof(nshArray))
         CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,"cudamemcpy")
      end if
#endif
      end

      module subroutine shmem_real_req
     &                 (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &                  nshArray,d_nshArray,
#endif
     &                  config)
      implicit none
      real(t_p),pointer:: shArray(:)
      integer(4) request_shape(:)
      integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
      type(rDPC),allocatable,optional:: nshArray(:)
      type(rDPC),device,pointer,optional:: d_nshArray(:)
#endif
      integer(4)  ,optional::config

      integer(8) request_size
      integer(4) configure,i,istat
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4) :: disp_unit,ierr
      real(t_p) adata
      TYPE(C_PTR) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
      integer(c_size_t) sh_size
#endif

      if (present(config)) then
         configure = config
      else
         configure = mhostonly
      end if
c
c     Deallocation
c
      if (associated(shArray)) then
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - size(shArray)*szoTp
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - size(shArray)*szoTp
            CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            CALL MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if
#if defined(_CUDA) && defined(USE_NVSHMEM)
      ! Deallocate nvshmem shared memory data structure
      if (present(nshArray).and.btest(configure,memnvsh)) then
         if (allocated(nshArray)) then
            do i = 0,npes-1
               if (i.eq.mype) then
                  sd_shmem = sd_shmem - size(nshArray(mype)%pel)*szoTp
                  peptr = c_devloc(nshArray(mype)%pel)
                  call nvshmem_free(peptr)
               else
                  nullify(nshArray(i)%pel)
               end if
            end do
            deallocate(nshArray)
            deallocate(d_nshArray)
         end if
      end if
#endif

      ! Fetch for size per process
      request_size = request_shape(1)
      if (request_size.eq.0) return

      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                MPI_SHARED_TYPE
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c     allocation
c
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape)
        s_shmem = s_shmem + size(shArray)*szoTp
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + size(shArray)*szoTp
      end if

#if defined(USE_NVSHMEM) && defined(_CUDA)
      if (present(nshArray).and.present(d_nshArray).and.
     &    btest(configure,memnvsh)) then

         allocate(nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "Allocate")
         allocate(d_nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "device Allocate")
         sh_size = value_pe(request_size)
         cuptr   = nvshmem_calloc(sh_size*szoTp )
         call c_f_pointer(cuptr, nshArray(mype)%pel, [sh_size])
         sd_shmem = sd_shmem + size(nshArray(mype)%pel)*szoTp
         do i = 0, npes-1
            if (i.ne.npes) then
               peptr = nvshmem_ptr( cuptr,i )
               call c_f_pointer(peptr,nshArray(i)%pel,[sh_size])
            end if
         end do
         !d_nshArray = nshArray
         istat = cudamemcpy(c_devloc(d_nshArray(0)),c_loc(nshArray(0)),
     &                      sizeof(nshArray))
         CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "cudamemcpy")
      end if
#endif
      end

      module subroutine shmem_realm_req
     &                 (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &                  nshArray,d_nshArray,
#endif
     &                  config)
      implicit none
      real(r_p),pointer:: shArray(:)
      integer(4) request_shape(:)
      integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
      type(rDPC),allocatable,optional:: nshArray(:)
      type(rDPC),device,pointer,optional:: d_nshArray(:)
#endif
      integer(4)  ,optional::config

      integer(8) request_size
      integer(4) configure,i,istat
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4) :: disp_unit,ierr
      real(t_p) adata
      TYPE(C_PTR) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
      integer(c_size_t) sh_size
#endif

      if (present(config)) then
         configure = config
      else
         configure = mhostonly
      end if
c
c     Deallocation
c
      if (associated(shArray)) then
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - size(shArray)*szoRp
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - size(shArray)*szoRp
            CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            CALL MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if
#if defined(_CUDA) && defined(USE_NVSHMEM)
      ! Deallocate nvshmem shared memory data structure
      if (present(nshArray).and.btest(configure,memnvsh)) then
         if (allocated(nshArray)) then
            do i = 0,npes-1
               if (i.eq.mype) then
                  sd_shmem = sd_shmem - size(nshArray(mype)%pel)*szoTp
                  peptr = c_devloc(nshArray(mype)%pel)
                  call nvshmem_free(peptr)
               else
                  nullify(nshArray(i)%pel)
               end if
            end do
            deallocate(nshArray)
            deallocate(d_nshArray)
         end if
      end if
#endif

      ! Fetch for size per process
      request_size = request_shape(1)
      if (request_size.eq.0) return

      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                MPI_SHMRED_TYPE
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c     allocation
c
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape)
        s_shmem = s_shmem + size(shArray)*szoRp
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + size(shArray)*szoRp
      end if

#if defined(USE_NVSHMEM) && defined(_CUDA)
      if (present(nshArray).and.present(d_nshArray).and.
     &    btest(configure,memnvsh)) then

         allocate(nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "Allocate")
         allocate(d_nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "device Allocate")
         sh_size = value_pe(request_size)
         cuptr   = nvshmem_calloc(sh_size*szoTp )
         call c_f_pointer(cuptr, nshArray(mype)%pel, [sh_size])
         sd_shmem = sd_shmem + size(nshArray(mype)%pel)*szoTp
         do i = 0, npes-1
            if (i.ne.npes) then
               peptr = nvshmem_ptr( cuptr,i )
               call c_f_pointer(peptr,nshArray(i)%pel,[sh_size])
            end if
         end do
         !d_nshArray = nshArray
         istat = cudamemcpy(c_devloc(d_nshArray(0)),c_loc(nshArray(0)),
     &                      sizeof(nshArray))
         CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "cudamemcpy")
      end if
#endif
      end

      module subroutine shmem_real_req2
     &                 (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &                  nshArray,d_nshArray,
#endif
     &                  config)
      implicit none
      real(t_p),pointer:: shArray(:,:)
      integer(4) request_shape(:)
      integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
      type(r2dDPC),allocatable,optional:: nshArray(:)
      type(r2dDPC),device,pointer,optional:: d_nshArray(:)
#endif
      integer(4)  ,optional::config

      integer(8) request_size
      integer(4) configure,i,istat
      integer(4) stride,nelem
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4) :: disp_unit,ierr
      real(t_p) adata
      TYPE(C_PTR) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
#endif
      integer(c_size_t) sh_size

      if (present(config)) then
         configure = config
      else
         configure = mhostonly
      end if

c
c     Deallocation
c
      if (associated(shArray)) then
         sh_size = size(shArray,1)*size(shArray,2)*szoTp
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - sh_size
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - sh_size
            CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            CALL MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if
#if defined(_CUDA) && defined(USE_NVSHMEM)
      ! Deallocate nvshmem shared memory data structure
      if (present(nshArray).and.btest(configure,memnvsh)) then
         if (allocated(nshArray)) then
            sh_size = size(nshArray(mype)%pel)*szoTp
            do i = 0,npes-1
               if (i.eq.mype) then
                  sd_shmem = sd_shmem - sh_size
                  peptr = c_devloc(nshArray(mype)%pel)
                  call nvshmem_free(peptr)
               else
                  nullify(nshArray(i)%pel)
               end if
            end do
            deallocate(nshArray)
            deallocate(d_nshArray)
         end if
      end if
#endif

      ! Fetch for size per process
      request_size = int(request_shape(1),8)*request_shape(2)
      if (request_size.eq.0) return

      stride = request_shape(1)
      nelem  = request_shape(2)
      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                MPI_SHARED_TYPE
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c     allocation
c
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape(1:2))
        s_shmem = s_shmem + request_size*szoTp
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + request_size*szoTp
      end if

#if defined(USE_NVSHMEM) && defined(_CUDA)
      if (present(nshArray).and.present(d_nshArray).and.
     &    btest(configure,memnvsh)) then

         allocate(nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "Allocate")
         allocate(d_nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "device Allocate")
         sh_size = value_pe(nelem)
         request_size = stride*sh_size*szoTp
         cuptr   = nvshmem_malloc( request_size )
         call c_f_pointer(cuptr, nshArray(mype)%pel, [stride,sh_size])
         sd_shmem = sd_shmem + request_size
         do i = 0, npes-1
            if (i.ne.npes) then
               peptr = nvshmem_ptr( cuptr,i )
               call c_f_pointer(peptr,nshArray(i)%pel,[stride,sh_size])
            end if
         end do
         !d_nshArray = nshArray
         istat = cudamemcpy(c_devloc(d_nshArray(0)),c_loc(nshArray(0)),
     &                      sizeof(nshArray))
         CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "cudamemcpy")
      end if
#endif
      end

      module subroutine shmem_real_req3
     &                 (shArray,winarray,request_shape,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &                  nshArray,d_nshArray,
#endif
     &                  config)
      implicit none
      real(t_p),pointer:: shArray(:,:,:)
      integer(4) request_shape(:)
      integer(4) winarray
#if defined(_CUDA) && defined(USE_NVSHMEM)
      type(r3dDPC),allocatable,optional:: nshArray(:)
      type(r3dDPC),device,pointer,optional:: d_nshArray(:)
#endif
      integer(4)  ,optional::config

      integer(8) request_size
      integer(4) configure,i,istat
      integer(4) stride,stride2,nelem
      integer(kind=MPI_ADDRESS_KIND) :: windowsize
      integer(4) :: disp_unit,ierr
      real(t_p) adata
      TYPE(C_PTR) :: baseptr
#ifdef _CUDA
      type(c_devptr):: peptr,cuptr
#endif
      integer(c_size_t) sh_size

      if (present(config)) then
         configure = config
      else
         configure = mhostonly
      end if

c
c     Deallocation
c
      if (associated(shArray)) then
         sh_size = size(shArray,1)*size(shArray,2)*
     &             size(shArray,3)*szoTp
         if (btest(configure,memacc)) then
            sd_ddmem = sd_ddmem - sh_size
!$acc exit data delete(shArray)
         end if
         if (btest(configure,memhost)) then
            s_shmem = s_shmem - sh_size
            CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $           disp_unit, baseptr, ierr)
            CALL MPI_Win_free(winarray,ierr)
            nullify(shArray)
         end if
      end if
#if defined(_CUDA) && defined(USE_NVSHMEM)
      ! Deallocate nvshmem shared memory data structure
      if (present(nshArray).and.btest(configure,memnvsh)) then
         if (allocated(nshArray)) then
            sh_size = size(nshArray(mype)%pel,1)*szoTp*
     &      size(nshArray(mype)%pel,2)*size(nshArray(mype)%pel,3)
            do i = 0,npes-1
               if (i.eq.mype) then
                  sd_shmem = sd_shmem - sh_size
                  peptr = c_devloc(nshArray(mype)%pel)
                  call nvshmem_free(peptr)
               else
                  nullify(nshArray(i)%pel)
               end if
            end do
            deallocate(nshArray)
            deallocate(d_nshArray)
         end if
      end if
#endif

      ! Fetch for size per process
      request_size = int(request_shape(1),8)*request_shape(2)*
     &               request_shape(3)
      if (request_size.eq.0) return

      stride = request_shape(1)
      stride2= request_shape(2)
      nelem  = request_shape(3)
      if (hostrank == 0) then
         windowsize = int(request_size,MPI_ADDRESS_KIND)*
     &                MPI_SHARED_TYPE
      else
         windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c     allocation
c
      if (btest(configure,memhost)) then
        CALL MPI_Win_allocate_shared(windowsize, disp_unit,
     $       MPI_INFO_NULL, hostcomm, baseptr, winarray, ierr)

        if (hostrank.ne.0) then
           CALL MPI_Win_shared_query(winarray, 0, windowsize,
     $          disp_unit, baseptr, ierr)
        end if

        !association with fortran pointer
        CALL C_F_POINTER(baseptr,shArray,request_shape(1:3))
        s_shmem = s_shmem + request_size*szoTp
      end if

      if (btest(configure,memacc)) then
!$acc enter data create(shArray)
         sd_ddmem = sd_ddmem + request_size*szoTp
      end if

#if defined(USE_NVSHMEM) && defined(_CUDA)
      if (present(nshArray).and.present(d_nshArray).and.
     &    btest(configure,memnvsh)) then

         allocate(nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,"Allocate")
         allocate(d_nshArray(0:npes-1),stat=istat)
         CALL ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "device Allocate")
         sh_size = value_pe(nelem)
         request_size = stride*stride2*sh_size*szoTp
         cuptr   = nvshmem_malloc( request_size )
         call c_f_pointer(cuptr, nshArray(mype)%pel, 
     &                     [stride,stride2,sh_size])
         sd_shmem = sd_shmem + request_size
         do i = 0, npes-1
            if (i.ne.npes) then
               peptr = nvshmem_ptr( cuptr,i )
               call c_f_pointer(peptr,nshArray(i)%pel,
     &                        [stride,stride2,sh_size])
            end if
         end do
         !d_nshArray = nshArray
         istat = cudamemcpy(c_devloc(d_nshArray(0)),c_loc(nshArray(0)),
     &                      sizeof(nshArray))
         CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,
     &        "cudamemcpy")
      end if
#endif
      end

      module subroutine shmem_update_device_int(src,n,async_q,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &           dst,nd,
#endif
     &           config)
      implicit none
      integer(4) ,intent(in):: n
      integer(4) ,intent(in):: src(*)
      integer(4) ,optional,intent(in):: async_q
#if defined(_CUDA) && defined(USE_NVSHMEM)
      integer(4) ,optional,intent(in):: nd
      integer(4),device,optional:: dst(*)
#endif
      integer(4) ,optional,intent(in):: config

      integer(4) configure,nd_
      integer(4) l_beg,l_end,istat
      integer(4) length

      if (present(config)) then
         configure = config
      else
         configure = ishft(1,memacc)
      end if

      if (btest(configure,memacc)) then
         if (present(async_q)) then
!$acc update device(src(1:n)) async(async_q)
         else
!$acc update device(src(1:n))
         end if
      end if
#if defined(_CUDA) && defined(USE_NVSHMEM)
      if (btest(configure,memnvsh)) then
         if (present(nd).and.present(dst)) then
            l_beg  =        mype *nd + 1
            l_end  = min((1+mype)*nd , n)
            length = l_end-l_beg+1
            !dst(1:length) = src(l_beg:l_end)
            istat  = cudamemcpyasync(dst,src(l_beg),length)
            CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,
     &           "cudamemcpy")
         else
12    format("WARNING shmem_update_device !!!",/,
     &       "Destination or size is missing to perform copy",/,
     &       " *** Wrong way to call for the routine;",
     &       " source code alteration *** ")
            print 12
         end if
      end if
#endif
      end subroutine

      module subroutine shmem_update_device_real(src,n,async_q,
#if defined(_CUDA) && defined(USE_NVSHMEM)
     &           dst,nd,
#endif
     &           config)
      implicit none
      integer(4) ,intent(in):: n
      real(t_p) ,intent(in):: src(*)
      integer(4) ,optional,intent(in):: async_q
#if defined(_CUDA) && defined(USE_NVSHMEM)
      integer(4) ,optional,intent(in):: nd
      real(t_p),device,optional:: dst(*)
#endif
      integer(4) ,optional,intent(in):: config

      integer(4) configure,nd_
      integer(4) l_beg,l_end,istat
      integer(4) length

      if (present(config)) then
         configure = config
      else
         configure = ishft(1,memacc)
      end if

      if (btest(configure,memacc)) then
         if (present(async_q)) then
!$acc update device(src(1:n)) async(async_q)
         else
!$acc update device(src(1:n))
         end if
      end if
#if defined(_CUDA) && defined(USE_NVSHMEM)
      if (btest(configure,memnvsh)) then
         if (present(nd).and.present(dst)) then
            l_beg  =        mype *nd + 1
            l_end  = min((1+mype)*nd , n)
            length = l_end-l_beg+1
            !dst(1:length) = src(l_beg:l_end)
            istat  = cudamemcpyasync(dst,src(l_beg),length)
            CALL CUDA_ERROR_CHECK(istat,__FILE__,__LINE__,
     &           "cudamemcpy")
         else
12    format("WARNING shmem_update_device_real !!!",/,
     &       "Destination or size is missing to perform copy",/,
     &       " *** Wrong way to call for the routine;",
     &       " source code alteration *** ")
            print 12
         end if
      end if
#endif
      end subroutine

      ! Check for Error
      subroutine ERROR_CHECK(istat,filename,line,error_type)
      implicit none
      integer(4) ,intent(in) :: istat
      character(*),intent(in)::filename
      integer(4),intent(in)  :: line
      character(*),intent(in):: error_type
      integer(4) ,parameter :: SUCCESS=0

65    format ("FORTRAN ASSERT: ",A,A,1x,6I,/,15x,A,1x,"line",5I)

      if (istat.ne.SUCCESS) then
         write(0,65) error_type," Error :",istat,filename,line
         call fatal
      end if
      end subroutine

      logical function IsNump(chara)
      implicit none
      character(1),intent(in):: chara
      select case (chara)
         case ('.')   ; IsNump = .true.
         case ('0')   ; IsNump = .true.
         case ('1')   ; IsNump = .true.
         case ('2')   ; IsNump = .true.
         case ('3')   ; IsNump = .true.
         case ('4')   ; IsNump = .true.
         case ('5')   ; IsNump = .true.
         case ('6')   ; IsNump = .true.
         case ('7')   ; IsNump = .true.
         case ('8')   ; IsNump = .true.
         case ('9')   ; IsNump = .true.
         case default ; IsNump = .false.
      end select
      end function

#ifdef _CUDA
      ! Check for CUDA Error
      subroutine CUDA_ERROR_CHECK(istat,filename,line,error_type)
      implicit none
      integer(4) ,intent(in) :: istat
      character(*),intent(in)::filename
      integer(4),intent(in)  :: line
      character(*),intent(in):: error_type
      integer(4) ,parameter :: SUCCESS=0

65    format ("CUDA ASSERT: ",A,A,1x,6I,A," at",2x,"line",5I,/,15x,A)

      if (istat.ne.cudasuccess) then
         write(0,65) error_type," Error :",istat,filename,line,
     &               cudageterrorstring(istat)
         call fatal
      end if
      end subroutine
#endif

      module subroutine nvshmem_get_HeapSize
      implicit none
      integer(4) ierr, i, length, istat
      character(64) hsize
      real(8) nv_size_f
      integer(4),parameter::sucess=0

      call get_environment_variable("NVSHMEM_SYMMETRIC_SIZE",
     &     hsize,length,status=ierr)

12    format("ERROR",I6," Reading NVSHMEM_SYMMETRIC_SIZE variable ",A,
     &  /,"Switching to default value")

      if (ierr.eq.sucess) then
         i = 1
         do while (i.le.length.and.(IsNump(hsize(i:i))
     &         .or.hsize(i:i).ne.'k'.or.hsize(i:i).ne.'K'
     &         .or.hsize(i:i).ne.'m'.or.hsize(i:i).ne.'M'
     &         .or.hsize(i:i).ne.'g'.or.hsize(i:i).ne.'G'
     &         .or.hsize(i:i).ne.'t'.or.hsize(i:i).ne.'T'))
            i = i + 1
         end do
         i = i-1

         if      (hsize(i:i).eq.'k'.or.hsize(i:i).eq.'K') then
            read (hsize(1:i-1),*) nv_size_f
            s_nvshmem = nv_size_f*2**10
         else if (hsize(i:i).eq.'m'.or.hsize(i:i).eq.'M') then
            read (hsize(1:i-1),*) nv_size_f
            s_nvshmem = nv_size_f*Mio
         else if (hsize(i:i).eq.'g'.or.hsize(i:i).eq.'G') then
            read (hsize(1:i-1),*) nv_size_f
            s_nvshmem = nv_size_f*2**30
         else if (hsize(i:i).eq.'t'.or.hsize(i:i).eq.'T') then
            read (hsize(1:i-1),*) nv_size_f
            s_nvshmem = nv_size_f*2**40
         else
            read(hsize(1:i),*,iostat=istat) s_nvshmem
            if (istat.ne.sucess)
     &         read(hsize(1:i-1),*,iostat=istat) s_nvshmem
            if (istat.ne.sucess) then
               print 12, istat, hsize(1:i)
               s_nvshmem = 30*Mio
            end if
         end if
      else
         s_nvshmem = 1030*Mio
      end if

      end subroutine

      ! Return nvshmem array index location in data struct
      module integer(4) function shmem_index(iglob,asize,locpe)
     &       result(ind)
      implicit none
      integer(4),intent(in)::iglob
      integer(4),intent(in)::asize
      integer(4),intent(out)::locpe
      locpe = (iglob-1)/asize
      ind   = mod(iglob-1,asize) + 1
      end function shmem_index

      module subroutine print_memory_usage()
        implicit none
        integer(4) i,ierr, hsize
        real(8) sm,pm,sdm,pdm,ddm,odm,wor,tdm,tm,adm
        real(8),save:: diff=0
        integer:: cin=1
 13     format(
     &  "iRank"," |","Host--",2x,"Shr",7x,"Pvt",5x,"Total"," |",
     &  "Device--",2x,"Shr",7x,"Pvt",6x,"Libs",5x,"Total"," |",
     &  1x,"WorkSpace",6x,"diff",2x,"Memory (MiB)")
 14     format(,I5," |",1x,3F10.3," |",3x,4F10.3," |",2F10.3)
 15     format(
     &  " Rank",I4,4x,"On host",F12.3,4x,"On device",F12.3)

        if (cin.gt.1) call WorkSpaceEstimation

        if (hostrank.eq.0) then
           sm = s_shmem/Mio
        else
           sm = 0
        end if

        pm  =   s_prmem/Mio
        sdm =  sd_shmem/Mio
        ddm =  sd_ddmem/Mio
        pdm =  sd_prmem/Mio
        wor = s_tinWork/Mio
        odm = (s_curand+s_nvshmem+s_driver+s_cufft)/Mio
        tdm = ddm + pdm + wor + odm
        tm  = sm + pm + wor
#ifdef _OPENACC
        diff= tdm - diff
#else
        diff= tm - diff
#endif

        if (rank.eq.0) write(0,13)
        call mpi_Comm_size(hostcomm,hsize,ierr)
        call mpi_barrier(hostcomm,ierr)
        do i = 0,hsize-1
           if (rank.eq.i)
     &     write(0,14) rank,sm,pm,tm,sdm+ddm,pdm,odm,tdm,wor,diff
           call mpi_barrier(hostcomm,ierr)
        end do
        cin = cin+1

#ifdef _OPENACC
        diff= tdm
#else
        diff= tm
#endif
      end subroutine

      subroutine WorkSpaceEstimation
      use potent
      implicit none
      integer(mipk) wev,wem,wec,wep,wps,wsf

      if (nproc.gt.1) then
         wsf = s_sfWork
      else
         wsf = 0
      end if
      wev = 0
      wem = 0
      wep = 0
      wec = 0
      wps = 0

c     if (use_mpole) then
c     wem = 3*npoleloc*szoTp + max(9*npolelocnl+3*nbloc,1)*szoTp
c     end if
c     if (use_polar) then
c     wps = 12*npolebloc*szoTp + 6*npolerecloc*szoTp + 6*nproc*szoi
c    &    + 10*npoleloc*szoTp
c    &    + 18*npoleloc*szoTp + 12*npolebloc*szoTp
c    &    + 6*npolerecloc*szoTp
c    &    + 3*npoleloc*szoTp
c     wep = 3*npoleloc*szoTp
c    &    + max (12*npolelocnl*szoTp, 40*npolerecloc*szoTp
c    &      + 4*n1mpimax*n2mpimax*n3mpimax*nrec_recep*szoTp)
c     end if

      s_tinWork = maxval( [wsf,wem,wep,wps,wec,wev] )
      end subroutine
c
#ifdef USE_NVSHMEM_CUDA
      module subroutine nvshmem_int_check( sol, dat, n, npe, config )
      implicit none
      integer(4) sol(*)
      integer(4),device::dat(*)
      integer(4) n,npe
      integer(4),optional::config

      integer(4) i,j,st,en,cfg,diff,diff1

      if (present(config)) then
         cfg = config
      else
         cfg = 0
      end if

      if (cfg.eq.1) then
!$acc enter data copyin(sol(1:n))
      end if

      diff = 0
      st = mype*npe +1
      en = min((mype+1)*npe,n)

!$acc parallel loop num_gangs(1) vector_length(512)
!$acc&         present(sol(1:n)) copy(diff)
      do i = st,en
         j = dat(i-st+1)
         if (sol(i).ne.j) then
!$acc atomic capture
         diff  = diff + 1
         diff1 = diff
!$acc end atomic
         if (diff1.lt.20) then
         print*,"diff",rank,i,sol(i),j
         end if
         end if
      end do

      if (cfg.eq.1) then
!$acc exit data delete(sol(1:n))
      end if

      call mpi_barrier(hostcomm,i)
 32   format('I Pe',I4," find",I8," diff(s) [int] ",
     &       " over",I8," elements")
      write(*,32) mype,diff,npe

      end subroutine
c
      module subroutine nvshmem_real_check( sol, dat, n, npe, config )
      implicit none
      real(t_p) sol(*)
      real(t_p),device::dat(*)
      integer(4) n,npe
      integer(4),optional::config

      integer(4) i,j,st,en,cfg,diff,diff1
      real(t_p) cap

      if (present(config)) then
         cfg = config
      else
         cfg = 0
      end if

      if (cfg.eq.1) then
!$acc enter data copyin(sol(1:n))
      end if

      diff = 0
      st = mype*npe +1
      en = min((mype+1)*npe,n)

!$acc parallel loop num_gangs(1) vector_length(256)
!$acc&         present(sol(1:n)) copy(diff)
      do i = st,en
         cap = dat(i-st+1)
         if (sol(i).ne.cap) then
!$acc atomic capture
         diff  = diff + 1
         diff1 = diff
!$acc end atomic
         if (diff1.lt.25) then
         print*,"diff",rank,i,sol(i),cap
         end if
         end if
      end do

      if (cfg.eq.1) then
!$acc exit data delete(sol(1:n))
      end if

      call mpi_barrier(hostcomm,i)
 32   format('I Pe',I4," find",I8," diff(s) [real] ",
     &       " over",I8," elements")
      write(*,32) mype,diff,npe

      end subroutine
#endif

      ! Reallocate `array` if it's size is smaller than `size`
      ! Array is also reallocated on GPU using `exit/enter data` directives
      ! `array` must NOT be declare create
      ! If async is true,
      module subroutine prmem_int_req( array, n, async, queue, config )
        implicit none
        integer(4), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config

        integer(4) cfg
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if ( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
 13     format(A13,L3,1I10,I4)
 14     format(A13,L3,2I10,I4)
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        if (.not. allocated(array)) then
           if (n.eq.0) return
           s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &              ,extra_alloc.and.n.ne.natoms)
           allocate(array(s_alloc))
           s_array = s_alloc*szoi
           s_prmem =  s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print 13,'alloc i',n.ne.s_alloc,n,rank
        else if (btest(cfg,memfree).or.n.eq.0) then
           s_array = size(array)*szoi
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           if(debMem) print*, 'dealloc i ',rank
           s_prmem = s_prmem - s_array
           deallocate(array)
        else
           if ( n > size(array) .or. n < 4*size(array)/5 ) then

              ! Inc realloca counter
              n_realloc_i = n_realloc_i + 1
              call enable_extra_alloc(0)
              ! Compute extra-reallocation if necessary
              s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                    ,extra_alloc.and.n.ne.natoms)
              if(debMem) print 14,'realloc i',n.ne.s_alloc
     &                           ,n,size(array),rank

              s_array = size(array)*szoi
              s_prmem = s_prmem - s_array
              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              deallocate(array)
              allocate(array(s_alloc))
              s_array = s_alloc*szoi
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if
           end if
        end if

        if (size(array).eq.0) print*,'tro prmem_int_req',n,rank

      end subroutine

      subroutine enable_extra_alloc(cfg)
      implicit none
      integer(4),intent(in)::cfg
      if (cfg.eq.0.and.n_realloc_i.gt.20) then
         if(debMem) print*, 'enable extra allocation i',rank
         extra_alloc = .true.
         n_realloc_i = 0
      else if(cfg.eq.1.and.n_realloc_r.gt.20) then
         if(debMem) print*, 'enable extra allocation r',rank
         extra_alloc = .true.
         n_realloc_r = 0
      end if
      end subroutine

      integer(4) module function get_prmem_alloc_size_i4(n) result(res)
      implicit none
      integer,intent(in):: n
      res = merge(n+int(real(n,8)*mem_inc),n
     &           ,extra_alloc.and.n.ne.natoms)
      end function

      module subroutine prmem_int8_req( array, n, async, queue, config )
        implicit none
        integer(8), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config

        integer(4) cfg
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if ( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
 13     format(A13,L3,1I14,I4)
 14     format(A13,L3,2I14,I4)
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        if (.not. allocated(array)) then
           s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                    ,extra_alloc.and.n.ne.natoms)
           allocate(array(s_alloc))
           s_array = s_alloc*szoi8
           s_prmem = s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print 13,'alloc i8 ',n.ne.s_alloc,n,rank
        else if (btest(cfg,memfree).or.n.eq.0) then
           s_array = size(array)*szoi8
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           !print*, 'deallocate array'
           s_prmem = s_prmem - s_array
           deallocate(array)
        else
           if ( n > size(array) .or. n < 4*size(array)/5 ) then
              if(debMem) print 14,'realloc i8',n.ne.s_alloc
     &                           ,n,size(array),rank
              s_array = size(array)*szoi8
              s_prmem = s_prmem - s_array
              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              deallocate(array)
              allocate(array(n))
              s_array = n*szoi8
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if
           end if
        end if

        !if (size(array).eq.0) print*,'trouble prmem_int8_req',n,rank

      end subroutine

      module subroutine prmem_int8_req1(array, n, async, queue, config)
        implicit none
        integer(8), allocatable, intent(inout) :: array(:)
        integer(8), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config

        integer(4) cfg
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if ( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
 13     format(A13,L3,1I14,I4)
 14     format(A13,L3,2I14,I4)
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        if (.not. allocated(array)) then
           s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                    ,extra_alloc.and.n.ne.natoms)
           allocate(array(s_alloc))
           s_array = s_alloc*szoi8
           s_prmem = s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print 13,'alloc i8 ',n.ne.s_alloc,n,rank
        else if (btest(cfg,memfree).or.n.eq.0) then
           s_array = size(array)*szoi8
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           !print*, 'deallocate array'
           s_prmem = s_prmem - s_array
           deallocate(array)
        else
           if ( n > size(array) .or. n < 4*size(array)/5 ) then
              if(debMem) print 14,'realloc i8',n.ne.s_alloc
     &                           ,n,size(array),rank
              s_array = size(array)*szoi8
              s_prmem = s_prmem - s_array
              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              deallocate(array)
              allocate(array(n))
              s_array = n*szoi8
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if
           end if
        end if

        !if (size(array).eq.0) print*,'trouble prmem_int8_req',n,rank

      end subroutine

      module subroutine prmem_int_req1( array, sz_array, n, 
     &                  async, queue, config )
        implicit none
        integer(4)   , allocatable, intent(inout) :: array(:)
        integer(8), intent(in)    :: n
        integer(8), intent(inout) :: sz_array
        logical   , optional, intent(in) :: async
        integer(4)   , optional, intent(in) :: queue, config

        integer(4)   cfg
        integer(8) s_alloc8
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if ( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
 13     format(A13,L3,1I14,I4)
 14     format(A13,L3,2I14,I4)
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        !TODO  Report to PGI
        !size(array,kind=8) !is not working with pgi
        if (.not. allocated(array)) then

           s_alloc8 = merge(n+int(real(n,8)*mem_inc,8),n
     &                     ,extra_alloc.and.n.ne.natoms)
           allocate(array(s_alloc8))
          sz_array = s_alloc8
           s_array = sz_array*szoi
           s_prmem = s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print*,'alloc i o1',n.ne.s_alloc8,n,rank

        else if (btest(cfg,memfree).or.n.eq.0) then

           s_array = sz_array*szoi
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           s_prmem = s_prmem - s_array
           deallocate(array)
          sz_array = 0

        else

           if ( n > sz_array .or. n < 4*sz_array/5 ) then
              s_alloc8 = merge(n+int(real(n,8)*mem_inc,8),n
     &                        ,extra_alloc.and.n.ne.natoms)
              if(debMem) print 14,'realloc i o1',n.ne.s_alloc8
     &                           ,n,sz_array,rank

              s_array = sz_array*szoi
              s_prmem = s_prmem - s_array
              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              deallocate(array)
              allocate(array(s_alloc8))
             sz_array = s_alloc8
              s_array = sz_array*szoi
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if
           end if
        end if
        !if(size(array).eq.0) print*,'trouble prmem_int_req1',n,rank

      end subroutine

      module subroutine prmem_4int1a_req( array, n,
     &                  async, queue, config )
        implicit none
        integer(4)   , allocatable, intent(inout) :: array(:)
        integer(8), intent(in)    :: n
        logical   , optional, intent(in) :: async
        integer(4)   , optional, intent(in) :: queue, config

        integer(4)   cfg
        integer(8) s_alloc8
        integer(int_ptr_kind()) s_array,sz_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if ( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
 13     format(A13,L3,1I14,I4)
 14     format(A13,L3,2I14,I4)
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        !TODO  Report to PGI
        !size(array,kind=8) !is not working with pgi
        if (.not. allocated(array)) then

           s_alloc8 = merge(n+int(real(n,8)*mem_inc,8),n
     &                     ,extra_alloc.and.n.ne.natoms)
           allocate(array(s_alloc8))
           s_array = s_alloc8*szoi
           s_prmem = s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print*,'alloc i o1',n.ne.s_alloc8,n,rank

        else if (btest(cfg,memfree).or.n.eq.0) then

           s_array = size(array)*szoi
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           s_prmem = s_prmem - s_array
           deallocate(array)

        else

           sz_array = size(array)
           if ( n > sz_array .or. n < 3*sz_array/5 ) then
              s_alloc8 = merge(n+int(real(n,8)*mem_inc,8),n
     &                        ,extra_alloc.and.n.ne.natoms)
              if(debMem) print 14,'realloc i o1',n.ne.s_alloc8
     &                           ,n,size(array),rank

              s_array = sz_array*szoi
              s_prmem = s_prmem - s_array
              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              deallocate(array)
              allocate(array(s_alloc8))
             sz_array = s_alloc8
              s_array = sz_array*szoi
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if
           end if
        end if
        !if(size(array).eq.0) print*,'trouble prmem_int_req1',n,rank

      end subroutine

      module subroutine prmem_pint_req( array, n, async, config )
        implicit none
        integer(4), pointer, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: config
        integer(4) cfg

        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if ( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
#endif
 14     format(A13,L3,2I10,2x,I4)

        cfg = mhostacc
        if (present(config)) cfg = config

        if (.not. associated(array)) then
           s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                    ,extra_alloc.and.n.ne.natoms)
           allocate(array(s_alloc))
           s_array = s_alloc*szoi
           s_prmem = s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
           sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print*,'alloc pi ',n.ne.s_alloc,n,rank
        else
           if ( n>size(array) .or. n<4*size(array)/5 ) then

              n_realloc_i = n_realloc_i + 1
              call enable_extra_alloc(0)
              s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                       ,extra_alloc.and.n.ne.natoms)
              if(debMem) print 14,'realloc pi ',n.ne.s_alloc
     &                           ,n,size(array),rank

              s_array = size(array)*szoi
              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              s_prmem = s_prmem - s_array
              deallocate(array)
              nullify(array)
              allocate(array(s_alloc))
              s_array = s_alloc*szoi
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if
           end if
        end if
        !if(size(array).eq.0) print*,'trouble prmem_pint_req',n,rank

      end subroutine
      module subroutine prmem_int1_req( array, n, async )
        implicit none
        integer(1), allocatable, intent(inout) :: array(:)
        integer(4)   , intent(in) :: n
        logical   , optional, intent(in) :: async

        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if ( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
#endif

        if (.not. allocated(array)) then
           s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &              ,extra_alloc.and.n.ne.natoms)
           allocate(array(s_alloc))
            s_array = s_alloc*szoi1
!$acc enter data create(array) async( async_queue )
            s_prmem =  s_prmem + s_array
           sd_prmem = sd_prmem + s_array
           if(debMem) print*,'alloc i1',n.ne.s_alloc,n,rank
        else
           if ( n > size(array) .or. n < 4*size(array)/5 ) then
              if(debMem) print*,'realloc i1  ',n.ne.s_alloc
     &                         ,n,size(array),rank
               s_array = size(array)*szoi1
               s_prmem =  s_prmem - s_array
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              deallocate(array)
              allocate(array(n))
              s_array  = n*szoi1
!$acc enter data create(array) async( async_queue )
               s_prmem =  s_prmem + s_array
              sd_prmem = sd_prmem + s_array
           end if
        end if
        !if(size(array).eq.0) print*,'trouble prmem_int1_req',n,rank

      end subroutine
      ! Request heap memory on 2d integer(4) array
      module subroutine prmem_int_req2( array, nl, nc, async, queue,
     &                  config, nlst, ncst, cdim )
        implicit none
        integer(4), allocatable, intent(inout) :: array(:,:)
        integer(4), intent(in) :: nc, nl
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config
        integer(4), optional, intent(in) :: nlst, ncst
        logical, optional, intent(in) :: cdim

        integer(4) ashape(2)
        integer(4) cfg, nlstr, ncstr
        integer(int_ptr_kind()) s_array
        logical f_col
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
 13     format(A13,L3,2I10,2x,I4)
 14     format(A13,L3,3I10,2x,I4)
 18     format("ERROR: prmem_int_req2 nl different from array",/,
     &         5x,"Cannot Procede to reallocation",/,
     &         5x,"old shape",2I10,1x,"request shape",2I10)

        ! Optional parameters
        cfg = mhostacc
        nlstr = 1; ncstr = 1; f_col = .true.
        if (present(config)) cfg = config
        if (present(nlst)) nlstr = nlst
        if (present(ncst)) ncstr = ncst
        if (present(cdim)) f_col = cdim

        if (.not.allocated(array)) then

           if (f_col) then
             s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &                      ,extra_alloc.and.nc.ne.natoms)
             allocate(array(nlstr:nl,ncstr:s_alloc))
             s_array = (nl-nlstr+1)*(s_alloc-ncstr+1)*szoi
           else
             s_alloc = merge(nl+int(real(nl,t_p)*mem_inc),nl
     &                      ,extra_alloc.and.nl.ne.natoms)
             allocate(array(nlstr:s_alloc,ncstr:nc))
             s_array = (s_alloc-nlstr+1)*(nc-ncstr+1)*szoi
           end if

           s_prmem = s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem =sd_prmem + s_array
           end if
           if(debMem) print 13,'alloc i 2',extra_alloc,nl,nc,rank

        else if (btest(cfg,memfree).or.nc*nl.eq.0) then
           ashape = shape(array)
           s_array = ashape(1)*ashape(2)*szoi
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           s_prmem = s_prmem - s_array
           deallocate(array)
        else
           ashape = shape(array)
           if ( nc>ashape(2).and.f_col ) then
              if (nl.ne.ashape(1)) then
                 print 18, nl,nc,ashape
                 call fatal
              end if

              n_realloc_i = n_realloc_i + 1
              s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &                       ,extra_alloc.and.nc.ne.natoms)
              if(debMem) print 14,'realloc i 2 ',nc.ne.s_alloc
     &                           ,nc,ashape,rank

              s_array = ashape(1)*ashape(2)*szoi
              s_prmem = s_prmem - s_array
              if (btest(cfg,memacc)) then
             sd_prmem =sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              deallocate(array)
              allocate(array(nlstr:nl,ncstr:s_alloc))
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              end if
              s_array = (nl-nlstr+1)*(s_alloc-ncstr+1)*szoi
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) sd_prmem =sd_prmem + s_array
           else if ( nl>ashape(1).and..not.f_col ) then

              n_realloc_i = n_realloc_i + 1
              s_alloc = merge(nl+int(real(nl,t_p)*mem_inc),nl
     &                       ,extra_alloc.and.nl.ne.natoms)
              if(debMem) print 14,'realloc i 2 ',nl.ne.s_alloc
     &                           ,nl,ashape,rank

              s_array = ashape(1)*ashape(2)*szoi
              s_prmem = s_prmem - s_array
              if (btest(cfg,memacc)) then
             sd_prmem =sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              deallocate(array)
              allocate(array(nlstr:s_alloc,ncstr:nc))
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              end if
              s_array = (s_alloc-nlstr+1)*(nc-ncstr+1)*szoi
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) sd_prmem =sd_prmem + s_array
           end if
        end if

c       if(size(array).eq.0) print*,'trouble prmem_int_req2',nl,nc,rank
      end subroutine
      module subroutine prmem_int1_req2( array, nl, nc, async )
        implicit none
        integer(1), allocatable, intent(inout) :: array(:,:)
        integer(4), intent(in) :: nc, nl
        logical, optional, intent(in) :: async

        integer(4) ashape(2)
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
#endif
  14    format (A13,L3,3I10,2x,I4)
  66    format ('error MOD_tinMemory:prmem_int1_req2 ',
     &          'improper array shape !! ', 2I4,' over ', 2I4,/,
     &          'cannot procede to reallocation ')

        if (.not.allocated(array)) then
           s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &                    ,extra_alloc.and.nc.ne.natoms)
           allocate(array(nl,s_alloc))
!$acc enter data create(array) async( async_queue )
           s_prmem = s_prmem + nl*s_alloc*szoi1
          sd_prmem =sd_prmem + nl*s_alloc*szoi1
        else
           ashape = shape(array)
           if ( nl.ne.ashape(1) ) then
              print 66, ashape, nl,nc
              call fatal
           end if
           if ( nc > ashape(2) ) then

              s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &                       ,extra_alloc.and.nc.ne.natoms)
              if(debMem) print 14,'realloc i1 2',nc.ne.s_alloc
     &                           ,nc,ashape,rank

              s_prmem = s_prmem - ashape(1)*ashape(2)*szoi1
             sd_prmem =sd_prmem - ashape(1)*ashape(2)*szoi1
!$acc exit data delete(array) async( async_queue )
              deallocate(array)
              allocate(array(nl,s_alloc))
!$acc enter data create(array) async( async_queue )
              s_prmem = s_prmem + nl*s_alloc*szoi1
             sd_prmem =sd_prmem + nl*s_alloc*szoi1
           end if
        end if
        !if(size(array).eq.0) print*,'trouble prmem_int1_req2',nl,nc,rank

      end subroutine

      module subroutine prmem_int8_req2( array, nl, nc, async, queue,
     &                  config, nlst, ncst )
        implicit none
        integer(8), allocatable, intent(inout) :: array(:,:)
        integer(4)   , intent(in) :: nc, nl
        logical   , optional   , intent(in) :: async
        integer(4)   , optional   , intent(in) :: queue, config
        integer(4)   , optional   , intent(in) :: nlst, ncst

        integer(4) ashape(2)
        integer(4) cfg, nlstr, ncstr
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        nlstr = 1; ncstr = 1;
        if (present(nlst)) nlstr = nlst
        if (present(ncst)) ncstr = ncst

        if (.not.allocated(array)) then
           s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &              ,extra_alloc.and.nc.ne.natoms)
           allocate(array(nlstr:nl,ncstr:s_alloc))
           s_array = (nl-nlstr+1)*(s_alloc-ncstr+1)*szoi8
!$acc enter data create(array) async( async_queue )
           s_prmem = s_prmem + s_array
          sd_prmem =sd_prmem + s_array
          !print*,'alloc i8 2',nl,nc,s_array
        else if (btest(cfg,memfree).or.nc*nl.eq.0) then
            ashape = shape(array)
           s_array = ashape(1)*ashape(2)*szoi8
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           s_prmem = s_prmem - s_array
           deallocate(array)
        else
           ashape = shape(array)
           if ( nc>ashape(2) .or. nl.ne.ashape(1) ) then
              if(debMem) print*,'realloc i8 2',nc,nl,ashape
              s_array = ashape(1)*ashape(2)*szoi8
              s_prmem = s_prmem - s_array
             sd_prmem =sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              deallocate(array)
              allocate(array(nlstr:nl,ncstr:nc))
!$acc enter data create(array) async( async_queue )
              s_array = (nl-nlstr+1)*(nc-ncstr+1)*szoi8
              s_prmem = s_prmem + s_array
             sd_prmem =sd_prmem + s_array
           end if
        end if
        !if(size(array).eq.0) print*,'trouble prmem_int8_req2',nl,nc,rank

      end subroutine

      module subroutine prmem_8int2p_req( array, nl, nc, async, queue,
     &                  config, nlst, ncst )
        implicit none
        integer(8), pointer, intent(inout) :: array(:,:)
        integer(4)   , intent(in) :: nc, nl
        logical   , optional   , intent(in) :: async
        integer(4)   , optional   , intent(in) :: queue, config
        integer(4)   , optional   , intent(in) :: nlst, ncst

        integer(4) ashape(2)
        integer(4) cfg, nlstr, ncstr
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        nlstr = 1; ncstr = 1;
        if (present(nlst)) nlstr = nlst
        if (present(ncst)) ncstr = ncst

        if (.not.associated(array)) then
           s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &              ,extra_alloc.and.nc.ne.natoms)
           allocate(array(nlstr:nl,ncstr:s_alloc))
           s_array = (nl-nlstr+1)*(s_alloc-ncstr+1)*szoi8
!$acc enter data create(array) async( async_queue )
           s_prmem = s_prmem + s_array
          sd_prmem =sd_prmem + s_array
          !print*,'alloc i8 2',nl,nc,s_array
        else if (btest(cfg,memfree).or.nc*nl.eq.0) then
            ashape = shape(array)
           s_array = ashape(1)*ashape(2)*szoi8
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           s_prmem = s_prmem - s_array
           deallocate(array)
        else
           ashape = shape(array)
           if ( nc>ashape(2) .or. nl.ne.ashape(1) ) then
              if(debMem) print*,'realloc i8 2',nc,nl,ashape
              s_array = ashape(1)*ashape(2)*szoi8
              s_prmem = s_prmem - s_array
             sd_prmem =sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              deallocate(array)
              allocate(array(nlstr:nl,ncstr:nc))
!$acc enter data create(array) async( async_queue )
              s_array = (nl-nlstr+1)*(nc-ncstr+1)*szoi8
              s_prmem = s_prmem + s_array
             sd_prmem =sd_prmem + s_array
           end if
        end if
        !if(size(array).eq.0) print*,'trouble prmem_8int2p_req',nl,nc,rank

      end subroutine

      module subroutine prmem_logi_req( array, n, async, queue, config )
        implicit none
        logical, allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config

        integer(4) cfg
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if ( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        if (.not. allocated(array)) then
           s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                    ,extra_alloc.and.n.ne.natoms)
           allocate(array(s_alloc))
           s_prmem =  s_prmem + sizeof(array)
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + sizeof(array)
           end if
           if (debMem) print*,'alloc l  ',n.ne.s_alloc,n,rank
        else if (btest(cfg,memfree).or.n.eq.0) then
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - sizeof(array)
!$acc exit data delete(array) async( async_queue )
           end if
           s_prmem = s_prmem - sizeof(array)
           deallocate(array)
        else
           if ( n > size(array) .or. n < 4*size(array)/5 ) then
              if(debMem) print*,'realloc l  ',n.ne.s_alloc
     &                         ,n,size(array),rank
              s_prmem = s_prmem - sizeof(array)
              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - sizeof(array)
!$acc exit data delete(array) async( async_queue )
              end if
              deallocate(array)
              allocate(array(n))
              s_prmem = s_prmem + sizeof(array)
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + sizeof(array)
              end if
           end if
        end if
        !if(size(array).eq.0) print*,'trouble prmem_logi_req',n,rank

      end subroutine

      ! Request heap memory on tinker real data
      module subroutine prmem_real_req( array, n, async, queue, config,
     &                  nst )
        implicit none
        real(t_p), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue
        integer(4), optional, intent(in) :: config, nst

        integer(4) nstr, cfg
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
 13     format(A13,L3,1I10,I4)
 14     format(A13,L3,2I10,I4)

        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if
        nstr = 1
        if (present(nst)) nstr=nst

        if (.not. allocated(array)) then
            s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &              ,extra_alloc.and.n.ne.natoms)
            if (s_alloc.eq.0) return
            allocate(array(nstr:s_alloc))
            if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
#ifdef _OPENACC
               if ( async_queue.ne.acc_async_sync ) then
               call mem_set(array,0.0_ti_p,int(s_alloc-nstr+1,mipk)
     &                     ,rec_stream)
               else
#endif
               call mem_set(array,0.0_ti_p,int(s_alloc-nstr+1,mipk)
     &                     ,int(0,mipk))
#ifdef _OPENACC
               end if
#endif
            end if
            s_array = s_alloc*szoTp
            s_prmem = s_prmem + s_array
           sd_prmem =sd_prmem + s_array
           if(debMem) print 13,'alloc r ',n.ne.s_alloc,n,rank
        else
            if( n-nstr+1 > size(array).or.nstr.ne.1 ) then

              n_realloc_r = n_realloc_r + 1
              call enable_extra_alloc(1)
              s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                       ,extra_alloc.and.n.ne.natoms)
              if(debMem) print 14,'realloc r ',n.ne.s_alloc
     &                         ,n,size(array),rank

              s_array = size(array)*szoTp
              if (btest(cfg,memacc)) then
                 sd_prmem =sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              s_prmem = s_prmem - s_array
              deallocate(array)
              if (s_alloc.eq.0) return
              allocate(array(nstr:s_alloc))
              s_array = s_alloc*szoTp
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
#ifdef _OPENACC
                if ( async_queue.ne.acc_async_sync ) then
                call mem_set(array,0.0_ti_p,int(s_alloc-nstr+1,mipk)
     &                      ,rec_stream)
                else
#endif
                call mem_set(array,0.0_ti_p,int(s_alloc-nstr+1,mipk)
     &                      ,int(0,mipk))
#ifdef _OPENACC
                end if
#endif
              sd_prmem =sd_prmem + s_array
              end if
            endif
        endif
        !if(size(array).eq.0) print*,'trouble prmem_real_req',n,rank

      end subroutine

      module subroutine prmem_realm_req( array, n, async )
        implicit none
        real(r_p), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
#endif
 14     format(A13,L3,2I10,I4)

        if(.not. allocated(array)) then
            s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                     ,extra_alloc.and.n.ne.natoms)
            allocate(array(s_alloc))
            s_array = s_alloc*szoRp
!$acc enter data create(array) async( async_queue )
#ifdef _OPENACC
            if ( async_queue.ne.acc_async_sync ) then
            call mem_set(array,0.0_re_p,int(s_alloc,mipk)
     &                  ,rec_stream)
            else
#endif
            call mem_set(array,0.0_re_p,int(s_alloc,mipk)
     &                  ,int(0,mipk))
#ifdef _OPENACC
            end if
#endif
            s_prmem =  s_prmem + s_array
           sd_prmem = sd_prmem + s_array
           if(debMem) print*,'alloc rm ',n,s_alloc
        else
            if( n > size(array) ) then

               n_realloc_r = n_realloc_r + 1
               call enable_extra_alloc(1)
               s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                        ,extra_alloc.and.n.ne.natoms)
               if(debMem) print 14,'realloc rm ',n.ne.s_alloc
     &                            ,n,size(array),rank

                s_array = size(array)*szoRp
                s_prmem =  s_prmem - s_array
               sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
               deallocate(array)
               allocate(array(s_alloc))
               s_array = s_alloc*szoRp
!$acc enter data create(array) async( async_queue )
#ifdef _OPENACC
               if ( async_queue.ne.acc_async_sync ) then
               call mem_set(array,0.0_re_p,int(s_alloc,mipk)
     &                     ,rec_stream)
               else
#endif
               call mem_set(array,0.0_re_p,int(s_alloc,mipk)
     &                     ,int(0,mipk))
#ifdef _OPENACC
               end if
#endif
                s_prmem =  s_prmem + s_array
               sd_prmem = sd_prmem + s_array
            endif
        end if
        !if(size(array).eq.0) print*,'trouble prmem_realm_req',n,rank

      end subroutine

      module subroutine prmem_realm_req1( array, n, async )
        implicit none
        real(r_p), allocatable, intent(inout) :: array(:)
        integer(mipk), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
#endif
 14     format(A13,L3,2I10,I4)

        if(.not. allocated(array)) then
            s_alloc = merge(n+int(real(n,8)*mem_inc,mipk),n
     &                     ,extra_alloc.and.n.ne.int(natoms,mipk))
            allocate(array(s_alloc))
            s_array = s_alloc*szoRp
!$acc enter data create(array) async( async_queue )
#ifdef _OPENACC
            if ( async_queue.ne.acc_async_sync ) then
            call mem_set(array,0.0_re_p,int(s_alloc,mipk)
     &                  ,rec_stream)
            else
#endif
            call mem_set(array,0.0_re_p,int(s_alloc,mipk)
     &                  ,int(0,mipk))
#ifdef _OPENACC
            end if
#endif
            s_prmem =  s_prmem + s_array
           sd_prmem = sd_prmem + s_array
           if(debMem) print*,'alloc rm ',n,s_alloc
        else
            if( n > size(array) ) then

               n_realloc_r = n_realloc_r + 1
               call enable_extra_alloc(1)
               s_alloc = merge(n+int(real(n,8)*mem_inc,mipk),n
     &                        ,extra_alloc.and.n.ne.natoms)
               if(debMem) print 14,'realloc rm ',n.ne.s_alloc
     &                            ,n,size(array),rank

                s_array = size(array)*szoRp
                s_prmem =  s_prmem - s_array
               sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
               deallocate(array)
               allocate(array(s_alloc))
               s_array = s_alloc*szoRp
!$acc enter data create(array) async( async_queue )
#ifdef _OPENACC
               if ( async_queue.ne.acc_async_sync ) then
               call mem_set(array,0.0_re_p,int(s_alloc,mipk)
     &                     ,rec_stream)
               else
#endif
               call mem_set(array,0.0_re_p,int(s_alloc,mipk)
     &                     ,int(0,mipk))
#ifdef _OPENACC
               end if
#endif
                s_prmem =  s_prmem + s_array
               sd_prmem = sd_prmem + s_array
            endif
        end if
        !if(size(array).eq.0) print*,'trouble prmem_realm_req',n,rank

      end subroutine

      module subroutine prmem_prealm_req( array, n, async )
        implicit none
        real(r_p), pointer, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
#endif
 13     format(A13,L3,1I10,I4)
 14     format(A13,L3,2I10,I4)

        if(.not. associated(array)) then
            s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                     ,extra_alloc.and.n.ne.natoms)
            allocate(array(s_alloc))
            s_array = s_alloc*szoRp
!$acc enter data create(array) async( async_queue )
            s_prmem =  s_prmem + s_array
           sd_prmem = sd_prmem + s_array
           if(debMem) print*,'alloc prm ',n,s_alloc
        else
            if( n > size(array) ) then

               n_realloc_r = n_realloc_r + 1
               s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                        ,extra_alloc.and.n.ne.natoms)
               if(debMem) print 14,'realloc pm ',n.ne.s_alloc
     &                            ,n,size(array),rank

                s_array = size(array)*szoRp
                s_prmem =  s_prmem - s_array
               sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
               deallocate(array)
               allocate(array(s_alloc))
               s_array = s_alloc*szoRp
!$acc enter data create(array) async( async_queue )
                s_prmem =  s_prmem + s_array
               sd_prmem = sd_prmem + s_array
            endif
        endif
        !if(size(array).eq.0) print*,'trouble prmem_prealm_req',n,rank

      end subroutine

      module subroutine prmem_preal_req( array, n, async )
        implicit none
        real(t_p), pointer, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
#endif
 13     format(A13,L3,1I10,I4)
 14     format(A13,L3,2I10,I4)

        if(.not. associated(array)) then
            s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                     ,extra_alloc.and.n.ne.natoms)
            allocate(array(s_alloc))
!$acc enter data create(array) async( async_queue )
            s_prmem =  s_prmem + sizeof(array)
           sd_prmem = sd_prmem + sizeof(array)
           if(debMem) print 13,'alloc pr ',n.ne.s_alloc,n,rank
        else
            if( n > size(array) ) then

               n_realloc_r = n_realloc_r + 1
               s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                        ,extra_alloc.and.n.ne.natoms)
               if(debMem) print 14,'realloc pr ',n.ne.s_alloc
     &                            ,n,size(array),rank

                s_prmem = s_prmem - sizeof(array)
               sd_prmem = sd_prmem - sizeof(array)
!$acc exit data delete(array) async( async_queue )
               deallocate(array)
               nullify(array)
               allocate(array(s_alloc))
!$acc enter data create(array) async( async_queue )
                s_prmem =  s_prmem + s_alloc*szoTp
               sd_prmem = sd_prmem + s_alloc*szoTp
            endif
        endif
        !if(size(array).eq.0) print*,'trouble prmem_preal_req',n,rank

      end subroutine
      ! Request heap memory on tinker 2D real data
      module subroutine prmem_real_req2( array, nl, nc, async, queue,
     &                  config, nlst, ncst )
        implicit none
        real(t_p), allocatable, intent(inout) :: array(:,:)
        integer(4)  , intent(in) :: nc, nl
        logical  , optional, intent(in) :: async
        integer(4)  , optional, intent(in) :: queue, config
        integer(4)  , optional, intent(in) :: nlst, ncst

        integer(4) ashape(2)
        integer(4) cfg, nlstr, ncstr
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        nlstr = 1; ncstr = 1;
        if (present(nlst)) nlstr = nlst
        if (present(ncst)) ncstr = ncst

  66    format ('error MOD_memory.f:prmem_real_req2 ',
     &          'improper array shape !! ', 2I4,' over ', 2I4,/,
     &          'cannot procede to reallocation ')

        if (.not.allocated(array)) then
           s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &              ,extra_alloc.and.nc.ne.natoms)
           allocate(array(nlstr:nl,ncstr:s_alloc))
           s_array = (nl-nlstr+1)*(s_alloc-ncstr+1)*szoTp
           s_prmem = s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
           sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print*,'alloc r 2',nc.ne.s_alloc,nl,nc
        else
           ashape  = shape(array)
           s_array = ashape(1)*ashape(2)*szoTp
           if (nl.ne.ashape(1)) then
              print 66, ashape, nl,nc
              call fatal
           end if
           if ( nc>ashape(2) .or. nl.ne.ashape(1) ) then
 12           format(A13,L3,3I10,2x,I5)
              n_realloc_r = n_realloc_r + 1
              call enable_extra_alloc(1)
              s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &                 ,extra_alloc.and.nc.ne.natoms)
              if(debMem) print 12,'realloc r 2',nc.ne.s_alloc
     &                         ,nc,ashape,rank

              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              s_prmem = s_prmem - s_array
              deallocate(array)
              allocate(array(nlstr:nl,ncstr:s_alloc))
              s_array = (nl-nlstr+1)*(s_alloc-ncstr+1)*szoTp
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if
           end if
        end if
        !if(size(array).eq.0) print*,'trouble prmem_real_req2',nl,nc,rank

      end subroutine
      ! Request heap memory on tinker 2D real mixed precision data
      module subroutine prmem_realm_req2( array, nl, nc, async )
        implicit none
        real(r_p), allocatable, intent(inout) :: array(:,:)
        integer(4)  , intent(in) :: nc, nl
        logical  , optional, intent(in) :: async
        integer(4) ashape(2)
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
#endif
  66    format ('error MOD_subMemory.f:prmem_realm_req2 ',
     &          'improper array shape !! ', 2I4,' over ', 2I4,/,
     &          'cannot procede to reallocation ')

        if (.not.allocated(array)) then
           s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &              ,extra_alloc.and.nc.ne.natoms)
           allocate(array(nl,s_alloc))
           s_array = nl*s_alloc*szoRp
!$acc enter data create(array) async( async_queue )
            s_prmem =  s_prmem + s_array
           sd_prmem = sd_prmem + s_array
           if(debMem) print*,'alloc rm 2',nc.ne.s_alloc,nl,nc
        else
           ashape = shape(array)
           if (nl.ne.ashape(1)) then
              print 66, ashape, nl,nc
              call fatal
           end if
           if ( nc > ashape(2) ) then
 12           format(A13,L3,3I10,2x,I5)
              n_realloc_r = n_realloc_r + 1
              s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &                       ,extra_alloc.and.nc.ne.natoms)
              if(debMem) print 12,'realloc rm 2',nc.ne.s_alloc
     &                         ,nc,ashape,rank

              s_array  = ashape(1)*ashape(2)*szoRp
               s_prmem =  s_prmem - s_array
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              deallocate(array)
              allocate  (array(nl,s_alloc))
              s_array  = nl*s_alloc*szoRp
               s_prmem =  s_prmem + s_array
              sd_prmem = sd_prmem + s_array
!$acc enter data create(array) async( async_queue )
           end if
        end if
      !if(size(array).eq.0) print*,'trouble prmem_realm_req2',nl,nc,rank

      end subroutine
      ! Request heap memory on tinker 2D real pointer mixed precision data
      module subroutine prmem_realpm_req2( array, nl, nc, async )
        implicit none
        real(r_p), pointer, intent(inout) :: array(:,:)
        integer(4)  , intent(in) :: nc, nl
        logical  , optional, intent(in) :: async
        integer(4) ashape(2)
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
#endif
  66    format ('error MOD_subMemory.f:prmem_realm_req2 ',
     &          'improper array shape !! ', 2I4,' over ', 2I4,/,
     &          'cannot procede to reallocation ')

        if (.not.associated(array)) then
           s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &              ,extra_alloc.and.nc.ne.natoms)
           allocate(array(nl,s_alloc))
           s_array = nl*s_alloc*szoRp
!$acc enter data create(array) async( async_queue )
            s_prmem =  s_prmem + s_array
           sd_prmem = sd_prmem + s_array
           if(debMem) print*,'alloc rm 2',nc.ne.s_alloc,nl,nc
        else
           ashape = shape(array)
           if (nl.ne.ashape(1)) then
              print 66, ashape, nl,nc
              call fatal
           end if
           if ( nc > ashape(2) ) then
 12           format(A13,L3,3I10,2x,I5)
              n_realloc_r = n_realloc_r + 1
              s_alloc = merge(nc+int(real(nc,t_p)*mem_inc),nc
     &                       ,extra_alloc.and.nc.ne.natoms)
              if(debMem) print 12,'realloc rm 2',nc.ne.s_alloc
     &                         ,nc,ashape,rank

              s_array  = ashape(1)*ashape(2)*szoRp
               s_prmem =  s_prmem - s_array
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              deallocate(array)
              allocate  (array(nl,s_alloc))
              s_array  = nl*s_alloc*szoRp
               s_prmem =  s_prmem + s_array
              sd_prmem = sd_prmem + s_array
!$acc enter data create(array) async( async_queue )
           end if
        end if
      !if(size(array).eq.0) print*,'trouble prmem_realm_req2',nl,nc,rank

      end subroutine

      ! Request heap memory on tinker 3D real data
      module subroutine prmem_real_req3(array, nx,ny,nz, async, queue_,
     &                  config)
        implicit none
        real(t_p), allocatable, intent(inout) :: array(:,:,:)
        integer(4)  , intent(in) :: nx, ny, nz
        logical  , optional, intent(in) :: async
        integer(4)  , optional, intent(in) :: queue_, config
        integer(4) ashape(3)
        integer(4) cfg
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
        if( present(queue_) ) async_queue = queue_
#endif
  66    format ('error MOD_subMemory.f:prmem_real_req3 ',
     &          'improper array shape !! ', 3I4,' over ', 3I4,/,
     &          'cannot procede to reallocation ')

        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        if (.not.allocated(array)) then
           s_alloc = merge(nz+int(real(nz,t_p)*mem_inc),nz
     &              ,extra_alloc.and.nz.ne.natoms)
           allocate(array(nx,ny,s_alloc))
           s_array = nx*ny*s_alloc*szoTp
           s_prmem = s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
           sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print*,'alloc r 3',nz.ne.s_alloc,nx,ny,nz
        else if (btest(cfg,memfree).or.nx*ny*nz.eq.0) then
           ashape = shape(array)
           s_array= ashape(1)*ashape(2)*ashape(3)*szoTp
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           s_prmem = s_prmem - s_array
           deallocate(array)
        else
           ashape = shape(array)
           if (.not.btest(cfg,memrealloc).and.
     &        (nx.ne.ashape(1) .or. ny.ne.ashape(2))) then
              print 66, ashape, nx,ny,nz
              call fatal
           end if
           if ( btest(cfg,memrealloc).or.nz > ashape(3) ) then

 13           format(A13,L3,4I10,2x,I5)
              n_realloc_r = n_realloc_r + 1
              s_alloc = merge(nz+int(real(nz,t_p)*mem_inc),nz
     &                       ,extra_alloc.and.nz.ne.natoms)
              if(debMem) print 13,'realloc r 3',nz.ne.s_alloc
     &                           ,nz,ashape,rank

              s_array = ashape(1)*ashape(2)*ashape(3)*szoTp
              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              s_prmem = s_prmem - s_array
              deallocate(array)
              allocate(array(nx,ny,s_alloc))
              s_array = nx*ny*s_alloc*szoTp
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if
           end if
        end if
      !if(size(array).eq.0) print*,'trouble prmem_real_req3',ny,nz,rank

      end subroutine

      module subroutine prmem_realm_req3(array, nx,ny,nz, async, queue_,
     &                  config)
        implicit none
        real(r_p), allocatable, intent(inout) :: array(:,:,:)
        integer(4)  , intent(in) :: nx, ny, nz
        logical  , optional, intent(in) :: async
        integer(4)  , optional, intent(in) :: queue_, config
        integer(4) ashape(3)
        integer(4) cfg
        integer(int_ptr_kind()) s_array,so_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
        if( present(queue_) ) async_queue = queue_
#endif
  66    format ('error MOD_subMemory.f:prmem_real_req3 ',
     &          'improper array shape !! ', 3I4,' over ', 3I4,/,
     &          'cannot procede to reallocation ')

        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

        if (.not.allocated(array)) then
           s_alloc = merge(nz+int(real(nz,t_p)*mem_inc),nz
     &              ,extra_alloc.and.nz.ne.natoms)
           allocate(array(nx,ny,s_alloc))
           s_array = nx*ny*s_alloc*szoTp
           s_prmem = s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
           sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print*,'alloc rm 3',nz.ne.s_alloc,nz,rank
        else if (btest(cfg,memfree).or.nx*ny*nz.eq.0) then
           ashape = shape(array)
           s_array= ashape(1)*ashape(2)*ashape(3)*szoTp
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           s_prmem = s_prmem - s_array
           deallocate(array)
        else
           ashape = shape(array)
c          if (.not.btest(cfg,memrealloc).and.
c    &        (nx.ne.ashape(1) .or. ny.ne.ashape(2))) then
c             print 66, ashape, nx,ny,nz
c             call fatal
c          end if
           so_array= ashape(1)*ashape(2)*ashape(3)*szoTp
           s_array = nx*ny*nz*szoTp
           if ( btest(cfg,memrealloc).or.s_array > so_array ) then

 13           format(A13,3I10,3x,3I10,2x,I4)
              if(debMem) print 13,'realloc rm 3',nx,ny,nz,ashape,rank

              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - so_array
!$acc exit data delete(array) async( async_queue )
              end if
              s_prmem = s_prmem - so_array
              deallocate(array)
              allocate(array(nx,ny,nz))
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if
           end if
        end if
      !if(size(array).eq.0) print*,'trouble prmem_realm_req3',ny,nz,rank

      end subroutine

      ! Request heap memory on tinker 4D real data
      module subroutine prmem_real_req4(array, nx, ny, nz, nc, async)
        implicit none
        real(t_p), allocatable, intent(inout) :: array(:,:,:,:)
        integer(4)  , intent(in) :: nx, ny, nz, nc
        logical  , optional, intent(in) :: async
        integer(4) ashape(4)
        integer(int_ptr_kind()) s_array
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if( present(async) ) then
            if( async ) async_queue = acc_async_noval
        end if
#endif
  66    format ('error MOD_subMemory.f:prmem_real_req4 ',/,
     &          'improper array shape !! ', 4I4, ' over ', 4I4,/,
     &          'cannot procede to reallocation ')

        if (.not.allocated(array)) then
           s_alloc = merge(nz+int(real(nz,t_p)*mem_inc),nz
     &                  ,extra_alloc.and.nz.ne.natoms)
           allocate(array(nx,ny,s_alloc,nc))
!$acc enter data create(array) async( async_queue )
           s_array  = nx*ny*s_alloc*nc*szoTp
            s_prmem =  s_prmem + s_array
           sd_prmem = sd_prmem + s_array
           if(debMem) print*,'alloc r 4',nx,ny,nz,nc,rank
        else
           ashape = shape(array)
           if (nx.ne.ashape(1) .or. ny.ne.ashape(2) .or.
     &         nc.ne.ashape(4) ) then
              print 66, ashape, nx,ny,nz,nc
              call fatal
           end if
           if (nz<ashape(3)) then
        !TODO  Add warning when nz<size(array,3)
        !      array should be pass as an assume shape arg to a routine
           end if
           if ( nz > ashape(3) ) then

 13           format(A13,L3,I10,2x,2I4,I10,I4,2x,I5)
              n_realloc_r = n_realloc_r + 1
              s_alloc = merge(nz+int(real(nz,t_p)*mem_inc),nz
     &                       ,extra_alloc.and.nz.ne.natoms)
              if(debMem) print 13,'realloc r 4',nz.ne.s_alloc
     &                           ,nz,ashape,rank

              s_array = ashape(1)*ashape(2)*ashape(3)*ashape(4)*szoTp
              s_prmem =  s_prmem - s_array
             sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              deallocate(array)
              allocate(array(nx,ny,s_alloc,nc))
              s_array = nx*ny*s_alloc*nc*szoTp
              s_prmem =  s_prmem + s_array
             sd_prmem = sd_prmem + s_array
!$acc enter data create(array) async( async_queue )
           end if
        end if

      end subroutine

c
c     Move allocation functions
c
      module subroutine prmem_int_mvreq( array, n, async,queue,config )
        implicit none
        integer(4), allocatable, intent(inout) :: array(:)
        integer(4), intent(in) :: n
        logical, optional, intent(in) :: async
        integer(4), optional, intent(in) :: queue, config

        integer(4) cfg,sz_array,i
        integer(int_ptr_kind()) s_array
        integer(4),allocatable:: buffer(:)
#ifdef _OPENACC
        integer(4) :: async_queue
        async_queue = acc_async_sync
        if ( present(async) ) then
           if( async ) async_queue = acc_async_noval
        end if
        if( present(queue) ) then
            async_queue = queue
        end if
#endif
        if (present(config)) then
           cfg = config
        else
           cfg = mhostacc
        end if

#ifdef _OPENACC
        ! Check mistake in configuation
        if (.not.btest(cfg,memacc)) then
11      format(" OpenACC ERROR : prmem_int_mvreq",/,
     &  "This function only operate on device when compiled with",
     &  " OpenACC options",/,
     &  "  --- Fetch for his call in source code to fix the issue")
           call fatal
        end if
#endif

        if (.not. allocated(array)) then
           s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                  ,extra_alloc.and.n.ne.natoms)
           allocate(array(s_alloc))
           s_array = s_alloc*szoi
           s_prmem =  s_prmem + s_array
           if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
           end if
           if(debMem) print*,'mvalloc i',n.ne.s_alloc,n,s_array
        else if (btest(cfg,memfree).or.n.eq.0) then
           s_array = size(array)*szoi
           if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
           end if
           !print*, 'deallocate array'
           s_prmem = s_prmem - s_array
           deallocate(array)
        else
           if ( n > size(array) ) then

 16           format(A,L3,2I10,2x,2I5)
 17           format(A,L3,3I10,2x,I5)
              n_realloc_i  = n_realloc_i + 1
              mv_realloc_i = mv_realloc_i + 1
              s_alloc = merge(n+int(real(n,t_p)*mem_inc),n
     &                  ,extra_alloc.and.n.ne.natoms)

              ! prmem_int_mvreq special extra allocation
              if (mv_realloc_i.gt.5) then
                 s_alloc = n+int(real(n,t_p)*mem_inc)
                 mv_realloc_i=0
              end if
              if(debMem) print 16,'mvrealloc i ',n.ne.s_alloc
     &                           ,n,size(array),n_realloc_i,rank

              sz_array = size(array)
              allocate(buffer(sz_array))
!$acc data create(buffer) async( async_queue )

!$acc parallel loop async( async_queue ) present(array)
              do i = 1,sz_array  ! Save array to buffer
                 buffer(i) = array(i)
              end do

              s_array = sz_array*szoi
              s_prmem = s_prmem - s_array
              if (btest(cfg,memacc)) then
              sd_prmem = sd_prmem - s_array
!$acc exit data delete(array) async( async_queue )
              end if
              deallocate(array)
              allocate(array(s_alloc))
              s_array = s_alloc*szoi
              s_prmem = s_prmem + s_array
              if (btest(cfg,memacc)) then
!$acc enter data create(array) async( async_queue )
              sd_prmem = sd_prmem + s_array
              end if

!$acc parallel loop async( async_queue ) present(array)
              do i = 1,sz_array
                 array(i) = buffer(i)
              end do

!$acc end data
              deallocate(buffer)
           end if
        end if

      end subroutine

      end submodule
