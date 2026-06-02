!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine initial  --  initial values and program setup  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "initial" sets up original values for some parameters
!     and variables that might not otherwise get initialized
!
!
#include "tinker_macro.h"
subroutine initial
   use atoms
   use beads
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
   use precis ,only:p_tiny=>tiny,small&
              ,p_huge=>huge
   use random_mod
   use timestat
   use tinMemory
   use tinheader
   use utilgpu
   use uprior
   use mpi
   use virial
   implicit none
!$ integer omp_get_num_procs
   integer ierr
   integer ::cuda_success=0
   real(t_p) precise
   integer :: status,hostnm,getpid,device_type,i=0

!
!     cores, thread count and options for OpenMP
!
   nbeads = 1
   nbeads_ctr = 0
   nproc = 1
   nthread = 1
!     call omp_set_num_threads (nthread)
   nrec = 0
!
!     default unit numbers for input and output
!
   input = 5
   iout  = 6
!
!     command line arguments to the program
!
   call command
!
!     Number of MPI processes and rank of the current MPI process
!
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproctot,ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,ranktot,ierr)
!     CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
!    $                         MPI_INFO_NULL, hostcomm,ierr)
!     CALL MPI_Comm_rank(hostcomm,hostrank,ierr)
!
!     Display copyright notice and license number
!
   if (ranktot.eq.0) call promo

!
!     values of machine precision constants
!
   p_tiny = precise (1)
   small  = precise (2)
   p_huge = precise (3)
!
!     number of lines in the keyfile
!
   nkey = 0
!
!     number of lines in the parameter file
!
   nprm = 0
!
!     number of atoms in the system
!
   n = 0
!
!     number of molecules in the system
!
!     nmol = 0
!
!     number of unit cell replicates
!
!     ncell = 0
!
!     number of mutated atoms in the system
!
   nprior = 0
!
!     flags for information levels within the program
!
   silent = .false.
   verbose = .false.
   debug = .false.
   abort = .false.
!
!     flags for periodic boundaries
!
   use_bounds = .false.
!      use_replica = .false.
   use_polymer = .false.
!
!     Stats on energy
!
   etot_ave = 0.0
   etot_std = 0.0
!
!     flags for temperature and pressure baths
!
   isothermal = .false.
   isobaric   = .false.
!
!     flags for rebuilding of neighbor lists
!
   auto_lst = .true.
   dovlst = .true.
   domlst = .true.
!
!     type of coordinates file
!
   coordtype = 'NONE'
!
!     atomic symbols for elements
!
   call initatom
!
!     Flag for predictor-Corrector
!
   polpred  = 'ASPC'
   use_pred = merge(.true.,.false.,app_id.eq.dynamic_a)
!
!     names of biopolymer residue types
!
   call initres
!
!     Initiate memory manager
!
   call init_memhandle
!
!     Default moldyn data
!
   nalt     = i_init
   nalt2    = i_init
   stepfast = i_init
   stepint  = i_init
!
!     default parameters used by optimizations
!
   fctmin   = 0.0_ti_p
   maxiter  = 0
   nextiter = 0
   iprint   = -1
   iwrite   = -1
   n_fwriten= 0
   stpmax   = 0.0_ti_p
!
!     initialize timer values
!
   call initiate_timers()
   call timer_enter(timer_prog)

   dotstgrad      = .false.
   calc_e         = .true.

end

subroutine set_nproc()
   use domdec
   use iounit
   use replicas
   use mpi
   use beads
   implicit none
   integer ierr,nbeads_para

   if (path_integral_md .AND. use_reps) then
      if (ranktot.eq.0) then
         write(iout,*) 'Error: PIMD and classical replicas&
            &     cannot be used together'
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call fatal
   end if

   if (use_reps) then
      ! check consistency of replica numbers for classical replicas
      ! and define nproc
      if(nreps<=0) then
         if (ranktot.eq.0) then
            write(iout,*) 'Error: number of replicas should be > 0'
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call fatal
      endif
      if (nproctot.lt.nreps) then
         nproc = 1
         if (ranktot.eq.0) then
            write(iout,*) 'each process should deal with max 1 replica'
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call fatal
      else if (mod(nproctot,nreps).ne.0) then
         if (ranktot.eq.0) then
            write(iout,*) 'Error: inconsistent number&
               &     of process for parallelism'
            write(iout,*) 'the total number of processors&
               &      should be lower&
               &     or a multiple of nreps'
            call fatal
         end if
      else
         nproc = nproctot/nreps
      end if
      if (ranktot==0) then
         write(iout,*) 'Using multiple replicas, nreps = ',nreps
      endif

   elseif(path_integral_md) then
      ! check consistency of replica numbers for PIMD
      if(contract .and. nbeads_ctr>1) then
         nbeads_para = nbeads_ctr
      else
         nbeads_para = nbeads
      endif
      if(nbeads_para<=0) then
         if (ranktot.eq.0) then
            write(iout,*) 'Error: number of beads should be > 0'
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call fatal
      endif


      if (nproctot.lt.nbeads_para) then
         nproc = 1
      else if (mod(nproctot,nbeads_para).ne.0) then
         if (ranktot.eq.0) then
            write(iout,*) 'Error: inconsistent number&
               &     of process for parallelism'
            write(iout,*) 'the total number of processors&
               &      should be lower&
               &     or a multiple of nbeads'
            call fatal
         end if
      else
         nproc = nproctot/nbeads_para
      end if

   else ! standard parallelization (no replicas)
      nproc = nproctot
   end if

end subroutine set_nproc
!
!     subroutine initmpi
!
!     initialize local communicator and generate replicas if necessary
!
subroutine initmpi
   use domdec
   use inform
   use mpi
   use tinMemory,only: extra_alloc
   use beads
   use replicas
   implicit none
   integer ierr,iproc,ncomm,color

   call set_nproc()

   ! compute rank for spatial communicator
   rank_reploc = int(ranktot/nproc)
   ncomm       = int(nproctot/nproc)
   if ((ncomm-nproc*nproctot).gt.0) ncomm = ncomm+1

   ! split MPI_COMM_WORLD into spatial communicator COMM_TINKER
   CALL MPI_Comm_split(MPI_COMM_WORLD,rank_reploc&
           ,ranktot,COMM_TINKER,ierr)
   call MPI_COMM_SIZE(COMM_TINKER,nproc,ierr)
   call MPI_COMM_RANK(COMM_TINKER,rank,ierr)

   CALL MPI_Comm_split_type(COMM_TINKER, MPI_COMM_TYPE_SHARED,0&
           ,MPI_INFO_NULL, hostcomm,ierr)
   CALL MPI_Comm_rank(hostcomm,hostrank,ierr)

   call initDebugEnv

   call initDevice

   COMM_ROOT2ROOT = MPI_COMM_NULL

   if (path_integral_md) then
      ! create polymer communicator
      CALL MPI_Comm_split(MPI_COMM_WORLD,rank,ranktot,COMM_POLYMER,ierr)
      call MPI_COMM_SIZE(COMM_POLYMER,nproc_polymer,ierr)
      call MPI_COMM_RANK(COMM_POLYMER,rank_polymer,ierr)

   else !if(use_reps) then
      ! create inter root communicator
      color = 1
      if (rank.eq.0) color = 0
      call MPI_Comm_split(MPI_COMM_WORLD,color,ranktot,COMM_ROOT2ROOT,ierr)

   endif

   if (nproc.gt.1) extra_alloc=.true.

end
!
!     Initialize GPU
!
subroutine initDevice
   use deriv  ,only: delambdae,delambdav,delambdaesave,&
      &delambdavsave
   use domdec ,only: rank,ranktot,hostcomm
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
!
!     Inititialize gpu stuff
!
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
!
!     For Mpi debug purpose
!
   if (btest(tinkerdebug,tindGdb)) then
      if (ranktot.eq.0) then
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
!
!     For blocking computation of direct space
!
   Ndir_block        = 1
   Ndir_async_block  = 0
!$acc enter data create(aMDwattype)
!
!     create virial component on device
!
   call create_vir_device
!
!     create energy data on device
!
   call create_energi_on_device()
!$acc enter data create(einter,delambdae,delambdav,delambdaesave, &
!$acc delambdavsave)
!
!     set to zero reductions buffers for (action,energy and virial)
!
!$acc enter data create(ered_buff,ered_buf1,lam_buff &
!$acc          ,vred_buff,vred_buf1,nred_buff)
   call zero_mdred_buffers
!
!     mapping parallelism for gpu execution
!
#ifdef _OPENACC
   call C_init_env(devicenum,nproc,rank,tinkerdebug)
   call copy_data_to_cuda_env(nproc)
#endif

   ngangs_rec  = 20
   gpu_gangs   = 80
   mod_gangs   = .true.
   gpu_workers = 1
   gpu_vector  = 32

!$acc update device(rank,ngpus,devicenum,gpu_gangs &
!$acc  ,gpu_workers,gpu_vector)
end subroutine

