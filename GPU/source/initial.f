c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine initial  --  initial values and program setup  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "initial" sets up original values for some parameters
c     and variables that might not otherwise get initialized
c
c
#include "tinker_macro.h"
      subroutine initial
      use atoms
      use bath
      use bound
      use cell
      use deriv
      use domdec
      use energi
      use files
      use inform
      use inter
      use iounit
      use keys
      use linmin
      use minima
      use moldyn
      use neigh
      use output
      use params
      use precis ,only:p_tiny=>tiny,small
     &           ,p_huge=>huge
      use random_mod
      use timestat
      use tinMemory
      use tinheader
      use utilgpu
      use mpi
      use virial
      implicit none
!$    integer omp_get_num_procs
      integer ierr
      integer ::cuda_success=0
      real(t_p) precise
      integer :: status,hostnm,getpid,device_type,i=0

c
c     cores, thread count and options for OpenMP
c
      nproc = 1
      nthread = 1
c     call omp_set_num_threads (nthread)
      nrec = 0
c
c     default unit numbers for input and output
c
      input = 5
      iout  = 6
c
c     command line arguments to the program
c
      call command
c
c     Number of MPI processes and rank of the current MPI process
c
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproctot,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,ranktot,ierr)
c     CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
c    $                         MPI_INFO_NULL, hostcomm,ierr)
c     CALL MPI_Comm_rank(hostcomm,hostrank,ierr)
c
c     Display copyright notice and license number
c
      if (ranktot.eq.0) call promo

c
c     values of machine precision constants
c
      p_tiny = precise (1)
      small  = precise (2)
      p_huge = precise (3)
c
c     number of lines in the keyfile
c
      nkey = 0
c
c     number of lines in the parameter file
c
      nprm = 0
c
c     number of atoms in the system
c
      n = 0
c
c     number of molecules in the system
c
c     nmol = 0
c
c     number of unit cell replicates
c
c     ncell = 0
c
c     number of mutated atoms in the system
c
      nprior = 0
c
c     flags for information levels within the program
c
      silent = .false.
      verbose = .false.
      debug = .false.
      abort = .false.
c
c     flags for periodic boundaries
c
      use_bounds = .false.
c      use_replica = .false.
      use_polymer = .false.
c
c     Stats on energy
c
      etot_ave = 0.0
      etot_std = 0.0
c
c     flags for temperature and pressure baths
c
      isothermal = .false.
      isobaric   = .false.
c
c     flags for rebuilding of neighbor lists
c
      auto_lst = .true.
      dovlst = .true.
      domlst = .true.
c
c     type of coordinates file
c
      coordtype = 'NONE'
c
c     atomic symbols for elements
c
      call initatom
c
c     names of biopolymer residue types
c
      call initres
c
c     Initiate memory manager
c
      call init_memhandle
c
c     Default moldyn data
c
      nalt     = i_init
      nalt2    = i_init
      stepfast = i_init
      stepint  = i_init
c
c     default parameters used by optimizations
c
      fctmin   = 0.0_ti_p
      maxiter  = 0
      nextiter = 0
      iprint   = -1
      iwrite   = -1
      n_fwriten= 0
      stpmax   = 0.0_ti_p
c
c     initialize timer values
c
      call initiate_timers()
      call timer_enter(timer_prog)

      dotstgrad      = .false.
      calc_e         = .true.

      end
c
c     subroutine initmpi
c
c     initialize local communicator and generate replicas if necessary
c
      subroutine initmpi
      use domdec
      use inform
      use mpi
      use tinMemory,only: extra_alloc
      implicit none
      integer ierr,iproc,rank_beadloc

      integer, allocatable :: bead_rank(:)
c
      allocate (bead_rank(nproctot))
c
c     not dealing with replicas for now
c
      nproc = nproctot
      bead_rank = 0
c
      do iproc = 0, nproctot-1
        rank_beadloc = int(iproc/nproc)
        bead_rank(iproc+1) = rank_beadloc
      end do

      CALL MPI_Comm_split(MPI_COMM_WORLD,bead_rank(ranktot+1),
     $     ranktot,COMM_TINKER,ierr)

      call MPI_COMM_SIZE(COMM_TINKER,nproc,ierr)
      call MPI_COMM_RANK(COMM_TINKER,rank,ierr)
      CALL MPI_Comm_split_type(COMM_TINKER, MPI_COMM_TYPE_SHARED, 0,
     $     MPI_INFO_NULL, hostcomm,ierr)
      CALL MPI_Comm_rank(hostcomm,hostrank,ierr)

      call initDebugEnv

      call initDevice

      if (nproc.gt.1) extra_alloc=.true.

      deallocate (bead_rank)
      end
c
c     Initialize GPU
c
      subroutine initDevice
      use deriv  ,only: delambdae,delambdav
      use domdec ,only: rank,hostcomm
      use energi
      use inter  ,only: einter
      use inform
      use mamd   ,only: aMDwattype
      use nvshmem
      use utilgpu
#ifdef _OPENACC
      use interfaces,only:C_init_env
      use utilcu    ,only: copy_data_to_cuda_env
#endif
      use sizes  ,only: tinkerdebug
      use virial
      !use tinMemory
      !use utils
      implicit none
      integer,volatile:: ierr,i
c
c     Inititialize gpu stuff
c
      ngpus         = 0
      nSMP          = 1
      nSPcores      = 1
      cores_SMP     = 1

#ifdef _OPENACC
      call selectDevice
#endif
      call init_nvshmem(hostcomm)
#ifdef _OPENACC
      call initDevice_stream
#endif
      !call mHalte
c
c     For Mpi debug purpose
c
      if (btest(tinkerdebug,tindGdb)) then
         if (rank.eq.0) then
            print*
            write(*,*) "*** Enabled Debug with [cuda-]Gdb ***"
            i = 0
            do while (i.eq.0)
               call sleep(2)
            end do
         end if
         call MPI_BARRIER(hostcomm,ierr)
      end if

      call set_warning
c
c     For blocking computation of direct space
c
      Ndir_block        = 1
      Ndir_async_block  = 0
!$acc enter data create(aMDwattype)
c
c     create virial component on device
c
      call create_vir_device
c
c     create energy data on device
c
      call create_energi_on_device()
!$acc enter data create(einter,delambdae,delambdav)
c
c     set to zero reductions buffers for (action,energy and virial)
c
!$acc enter data create(ered_buff,ered_buf1,lam_buff
!$acc&          ,vred_buff,vred_buf1,nred_buff)
      call zero_mdred_buffers
c
c     mapping parallelism for gpu execution
c
#ifdef _OPENACC
      call C_init_env(devicenum,nproc,rank,tinkerdebug)
      call copy_data_to_cuda_env(nproc)
#endif

      ngangs_rec  = 20
      gpu_gangs   = 80
      mod_gangs   = .true.
      gpu_workers = 1
      gpu_vector  = 32

!$acc update device(rank,ngpus,devicenum,gpu_gangs
!$acc&  ,gpu_workers,gpu_vector)
      end subroutine

      subroutine initmpi_reps()
      use domdec
      use iounit
      use inform
      use replicas
      use mpi
      use tinMemory,only: extra_alloc
      implicit none
      integer ierr,color,ncomm

      if (nproctot.lt.nreps) then
        nproc = 1
        if (ranktot.eq.0) then
          write(iout,*) 'each process should deal with max 1 replica'
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call fatal
      else if (mod(nproctot,nreps).ne.0) then
        if (ranktot.eq.0) then
          write(iout,*) 'Error: inconsistent number 
     &     of process for parallelism'
          write(iout,*) 'the total number of processors
     &      should be lower
     &     or a multiple of nreps'
          call fatal
        end if
      else
        nproc = nproctot/nreps
      end if
c
c     if nreps > 1, nproc is defined earlier (in dynamic_rep.f)
c     else, we specify that we use standard parallelization (nproc=nproctot)
c
      if(nreps==1) nproc=nproctot
      rank_reploc = int(ranktot/nproc)
      ncomm = int(nproctot/nproc)
      if ((ncomm-nproc*nproctot).gt.0) ncomm = ncomm+1

      CALL MPI_Comm_split(MPI_COMM_WORLD,rank_reploc,
     $     ranktot,COMM_TINKER,ierr)

      call MPI_COMM_SIZE(COMM_TINKER,nproc,ierr)
      call MPI_COMM_RANK(COMM_TINKER,rank,ierr)
      CALL MPI_Comm_split_type(COMM_TINKER, MPI_COMM_TYPE_SHARED, 0,
     $     MPI_INFO_NULL, hostcomm,ierr)
      CALL MPI_Comm_rank(hostcomm,hostrank,ierr)
      rank_reploc = int(ranktot/nproc)

      call initDebugEnv
c
      call initDevice
c
      if (nproc.gt.1) extra_alloc=.true.
c
c     create inter root communicator
c
      color = 1
      if (rank.eq.0) color = 0
      call MPI_Comm_split(MPI_COMM_WORLD,color,ranktot,COMM_ROOT2ROOT,
     $  ierr)
      
      end subroutine initmpi_reps

