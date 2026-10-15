  module thrust
    use iso_c_binding
    use cudafor

    private
    public :: thrust_sort,thrust_inclusive_scan,thrust_exclusive_scan
    public :: thrust_stable_sort_by_key,thrust_remove_zero,thrust_remove
    public :: thrust_cache_init,thrust_clear_cache_memory

    integer(cuda_count_kind),private :: numbytes=0
    logical                 ,private :: mem_alloc=.false.
    logical                 ,private :: cache_initialized=.false.
    type(c_devptr)          ,private :: thrust_cache_ptr=c_null_ptr

    interface thrust_sort

       subroutine thrust_sort_int    (input, N) &
                  bind(C,name="sort_int_wrapper")
         import  c_int
         integer(c_int),device :: input(*)
         integer(c_int),value  :: N
!DIR$ ignore_tkr (r) input
       end subroutine

       subroutine sort_int_async (input, N, stream) &
                  bind(C,name="sort_int_async_wrapper")
         import  c_int,cuda_stream_kind
         integer(c_int),device :: input(*)
         integer(c_int),value  :: N
         integer(cuda_stream_kind),value :: stream
!DIR$ ignore_tkr (r) input
       end subroutine

       subroutine thrust_sort_float  (input, N) &
                  bind(C,name="sort_float_wrapper")
         import  c_float,c_int
         real   (c_float),device :: input(*)
         integer(c_int  ),value  :: N
!DIR$ ignore_tkr (r) input
       end subroutine

       subroutine thrust_sort_double (input, N) &
                  bind(C,name="sort_double_wrapper")
         import  c_double,c_int
         real   (c_double),device :: input(*)
         integer(c_int   ),value  :: N
!DIR$ ignore_tkr (r) input
       end subroutine

     end interface

    interface thrust_inclusive_scan

       subroutine thrust_inclusive_scan_int  (in, N, out) &
                  bind(C,name="inclusive_scan_int_wrapper")
         import  c_int
         integer(c_int),device :: in(*)
         integer(c_int),value  :: N
         integer(c_int),device :: out(*)
!DIR$ ignore_tkr (r) in, (r) out
       end subroutine

       subroutine thrust_inclusive_scan_int_async (in, N, out, s) &
                  bind(C,name="inclusive_scan_int_async_wrapper")
         import  c_int,cuda_stream_kind
         integer(c_int),device :: in(*)
         integer(c_int),value  :: N
         integer(c_int),device :: out(*)
         integer(cuda_stream_kind),value :: s
!DIR$ ignore_tkr (r) in, (r) out
       end subroutine

    end interface

    interface thrust_exclusive_scan

       subroutine thrust_exclusive_scan_int  (in, N, out) &
                  bind(C,name="exclusive_scan_int_wrapper")
         import  c_int
         integer(c_int),device :: in(*)
         integer(c_int),value  :: N
         integer(c_int),device :: out(*)
!DIR$ ignore_tkr (r) in, (r) out
       end subroutine
       subroutine thrust_exclusive_scan_int_async  (in, N, out, s) &
                  bind(C,name="exclusive_scan_int_async_wrapper")
         import  c_int,cuda_stream_kind
         integer(c_int),device :: in(*)
         integer(c_int),value  :: N
         integer(cuda_stream_kind),value  :: s
         integer(c_int),device :: out(*)
!DIR$ ignore_tkr (r) in, (r) out
       end subroutine

    end interface

    interface thrust_stable_sort_by_key

       subroutine thrust_stable_sort_by_key_int  (in, N, val) &
                  bind(C,name="stable_sort_by_key_int_wrapper")
         import  c_int
         integer(c_int),device :: in(*)
         integer(c_int),value  :: N
         integer(c_int),device :: val(*)
!DIR$ ignore_tkr (r) in
       end subroutine
       subroutine thrust_stable_sort_by_key_int_async  (in, N, val, s) &
                  bind(C,name="stable_sort_by_key_int_async_wrapper")
         import  c_int,cuda_stream_kind
         integer(c_int),device :: in(*)
         integer(c_int),value  :: N
         integer(cuda_stream_kind),value  :: s
         integer(c_int),device :: val(*)
!DIR$ ignore_tkr (r) in
       end subroutine

    end interface

    interface thrust_remove_zero

       subroutine thrust_remove_zero_int (in, N) &
                  bind(C,name="remove_zero_wrapper")
         import c_int
         integer(c_int),device :: in(*)
         integer(c_int),value  :: N
!DIR$ ignore_tkr (r) in
       end subroutine
       subroutine thrust_remove_zero_async_int (in, N, s) &
                  bind(C,name="remove_zero_async_wrapper")
         import c_int,cuda_stream_kind
         integer(c_int),device :: in(*)
         integer(c_int),value  :: N
         integer(cuda_stream_kind),value  :: s
!DIR$ ignore_tkr (r) in
       end subroutine

    end interface

    interface thrust_remove

       subroutine thrust_remove_async_wrapper (in, N, val, new_last, s) &
                  bind(C,name="remove_async_wrapper")
          import c_int,cuda_stream_kind
          integer(c_int),device:: in(*)
          integer(c_int),value :: N,val
          integer(c_int)       :: new_last
          integer(cuda_stream_kind),value :: s
!DIR$ ignore_tkr (r) in
       end subroutine

    end interface

    interface
      subroutine thrust_cache_alloc ( ptr,n ) &
                 bind(C)
        import cuda_count_kind,c_devptr
        integer(cuda_count_kind),value::n
        type(c_devptr), ptr
      end subroutine
      subroutine thrust_cache_dealloc ( ) &
                 bind(C)
      end subroutine
    end interface

    contains

    subroutine thrust_cache_init(n, device_memory)
    implicit none
    integer,intent(in) :: n
    integer(8),intent(inout) :: device_memory
    integer :: size_integers

    if (.not.cache_initialized) then
       size_integers = max(250000,n*7)
       call thrust_alloc_cache_memory(size_integers,device_memory)
       device_memory = device_memory + size_integers*sizeof(n)
       cache_initialized = .true.
    end if
    end subroutine

    subroutine thrust_alloc_cache_memory(n, device_memory)
    use iso_c_binding ,only: c_associated
    implicit none
    integer,intent(in)::n
    integer(8),intent(inout)::device_memory
    integer ierr
    type(c_ptr)::ptrc
 20 format('Error deallocating device memory on thurst_module.f90 : return error ',I10,/,A40)
 30 format('Error allocating device memory on thurst_module.f90 : return error ',I10,/,A40)

    if (mem_alloc) call thrust_clear_cache_memory(device_memory)

    ierr = 0
    numbytes = n*sizeof(n)
    ierr     = cudaMalloc(thrust_cache_ptr, numbytes)
    if (ierr.ne.cudasuccess) print 30, ierr,cudageterrorstring(ierr)
    call thrust_cache_alloc( thrust_cache_ptr,numbytes )
    mem_alloc = .true.

    end subroutine

    subroutine thrust_clear_cache_memory(device_memory)
    implicit none
    integer(8),intent(inout) :: device_memory
    integer ierr
 20 format('Error deallocating device memory on thurst_module.f90 : return error ',I10,/,A40)

    ierr = 0
    if (mem_alloc) ierr = cudaFree(thrust_cache_ptr)
    if (ierr.ne.cudasuccess) print 20, ierr,cudageterrorstring(ierr)
    call thrust_cache_dealloc
    mem_alloc=.false.
    device_memory = device_memory - numbytes
    numbytes = 0
    cache_initialized = .false.
    end subroutine

  end module
