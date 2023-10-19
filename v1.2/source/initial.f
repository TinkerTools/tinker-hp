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
      use beads, only: nbeads,nbeads_ctr
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
      use replicas
      use timestat
      use mpi
      implicit none
!$    integer omp_get_num_procs
      integer ierr
      real*8 precise
c
c     replicas
c
      nreps = 0
c
c     cores, thread count and options for OpenMP
c
      nbeads = 1
      nbeads_ctr = 0
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
          write(iout,*) 'Error: PIMD and classical replicas 
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
        if (ranktot==0) then
          write(iout,*) 'Using multiple replicas, nreps = ',nreps
        endif

      elseif(path_integral_md) then
        ! check consistency of replica numbers for PIMD
        if(centroid_recip) then
          nbeads_para = 1
        elseif(contract .and. nbeads_ctr>1) then
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
            write(iout,*) 'Error: inconsistent number 
     &     of process for parallelism'
            write(iout,*) 'the total number of processors
     &      should be lower
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

      subroutine initmpi
      use domdec
      use inform
      use mpi
      use replicas
      use beads
      implicit none
      integer ierr
      integer ncomm,color

      call set_nproc()

      ! compute rank for spatial communicator
      rank_reploc = int(ranktot/nproc)
      ncomm = int(nproctot/nproc)
      if ((ncomm-nproc*nproctot).gt.0) ncomm = ncomm+1

      ! split MPI_COMM_WORLD into spatial communicator COMM_TINKER
      CALL MPI_Comm_split(MPI_COMM_WORLD,rank_reploc,
     $     ranktot,COMM_TINKER,ierr)
      call MPI_COMM_SIZE(COMM_TINKER,nproc,ierr)
      call MPI_COMM_RANK(COMM_TINKER,rank,ierr)


      CALL MPI_Comm_split_type(COMM_TINKER, MPI_COMM_TYPE_SHARED, 0,
     $  MPI_INFO_NULL, hostcomm,ierr)
      CALL MPI_Comm_rank(hostcomm,hostrank,ierr)

      !!! replica-specific communicators
      if(path_integral_md) then
        ! create polymer communicator
        CALL MPI_Comm_split(MPI_COMM_WORLD,rank,
     $    ranktot,COMM_POLYMER,ierr)
        call MPI_COMM_SIZE(COMM_POLYMER,nproc_polymer,ierr)
        call MPI_COMM_RANK(COMM_POLYMER,rank_polymer,ierr)

      elseif(use_reps) then
        ! create inter root communicator
        color = 1
        if (rank.eq.0) color = 0
        call MPI_Comm_split(MPI_COMM_WORLD,color,ranktot,COMM_ROOT2ROOT,
     $    ierr)

      endif

      end
