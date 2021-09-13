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
      subroutine initial
      use atoms
      use bath
      use bound
      use cell
      use deriv
      use domdec
      use files
      use inform
      use iounit
      use keys
      use linmin
      use minima
      use neigh
      use output
      use params
      use precis
      use timestat
      use mpi
      implicit none
!$    integer omp_get_num_procs
      integer ierr
      real*8 precise
c
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
      iout = 6
c
c     command line arguments to the program
c
      call command
c
c     Number of MPI processes and rank of the current MPI process
c
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproctot,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,ranktot,ierr)
c
c     Display copyright notice and license number
      if (ranktot.eq.0) call promo

c
c     values of machine precision constants
c
      tiny = precise (1)
      small = precise (2)
      huge = precise (3)
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
c      nmol = 0
cc
cc     number of unit cell replicates
cc
c      ncell = 0
c
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
c     flags for temperature and pressure baths
c
      isothermal = .false.
      isobaric = .false.
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
c     default parameters used by optimizations
c
      fctmin = 0.0d0
      maxiter = 0
      nextiter = 0
      iprint = -1
      iwrite = -1
      stpmax = 0.0d0
c
c     initialize timer values
c
      timestep = 0d0
      timeinte = 0.0d0
      timereneig = 0.0d0
      timecommpos = 0.0d0
      timeparam = 0.0d0
      timenl = 0d0
      timegrad = 0.0d0
      timered = 0.0d0
      timetp = 0.0d0
      timecommforces = 0.0d0
      timebonded = 0d0
      timevdw = 0d0
      timeelec = 0d0
      timepolar = 0d0
      timecleargrad = 0d0
      timereal = 0d0
      timerec = 0d0
      timecommforcesrec = 0d0
      timecommforcesreal = 0d0
      timegrid = 0d0
      timefft = 0d0
      timescalar = 0d0
      timegrid2 = 0d0
      timerecreccomm = 0d0
      dotstgrad = .false.
      return
      end
c
c     subroutine initmpi
c
c     initialize local communicator and generate replicas if necessary
c
      subroutine initmpi
      use domdec
      use mpi
      implicit none
      integer ierr,iproc,rank_beadloc

      integer, allocatable :: bead_rank(:)
c
      allocate (bead_rank(nproctot))
c
c     not dealing with replicas for now
c
      bead_rank = 0d0
c
      do iproc = 0, nproctot-1
        rank_beadloc = int(iproc/nproc)
        bead_rank(iproc+1) = rank_beadloc
      end do


      CALL MPI_Comm_split(MPI_COMM_WORLD,bead_rank(ranktot+1),
     $  ranktot,COMM_TINKER,ierr)

      call MPI_COMM_SIZE(COMM_TINKER,nproc,ierr)
      call MPI_COMM_RANK(COMM_TINKER,rank,ierr)
      CALL MPI_Comm_split_type(COMM_TINKER, MPI_COMM_TYPE_SHARED, 0,
     $  MPI_INFO_NULL, hostcomm,ierr)
      CALL MPI_Comm_rank(hostcomm,hostrank,ierr)
      deallocate (bead_rank)
c
      return
      end
