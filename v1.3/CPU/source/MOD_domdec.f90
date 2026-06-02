!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ############################################################################################
!     ##                                                                                        ##
!     ##  module domdec  --  system parameters for OpenMP/MPI domain decomposition computation  ##
!     ##                                                                                        ##
!     ############################################################################################
!
!     nproctot  !<total number of MPI process (within MPI_COMM_WORLD)
!     ranktot  !<total rank of the MPI process within MPI_COMM_WORLD
!     nxdd = !<number of subdivisions along the x axis
!     nydd = !<number of subdivisions along the y axis
!     nzdd = !<number of subdivisions along the z axis
!
!     COMM_TINKER !<local MPI communicator in which a dynamic, analyze, testgrad or minimize run
!      will take place
!     nproc     !<number of MPI processes during a dynamic, analyze, testgrad or minimize run
!     rank      !<rank of the current MPI process within COMM_TINKER
!     rank_bis  !<rank of the current MPI process within comm_dir or comm_rec
!     nrec      !<number of processes assigned to the computation of reciprocal space contribution
!     ndir      !<number of processes assigned to the computation of direct space contribution
!     comm_rec  !<MPI group communicator associated to the reciprocal space
!     comm_dir  !<MPI group communicator associated to the direct space
!     nthread   !<number of threads to be used with OpenMP
!     hostcomm  !<MPI group communicator associated to processes within a node
!     hostrank !<rank of the current MPI process within hostcomm
!
!     n_recep1  !<number of MPI process to receive positions from to compute electrostatic interactions
!     n_send1  !<number of MPI process to send positions to to compute electrostatic interactions
!     n_recep2  !<number of MPI process to receive positions from to compute vdw interactions
!     n_send2  !<number of MPI process to send positions to to compute vdw interactions
!
!     n_recepshort1  !<number of MPI process to receive positions from to compute short range electrostatic interactions
!     n_sendshort1  !<number of MPI process to send positions to to compute short range electrostatic interactions
!     n_recepshort2  !<number of MPI process to receive positions from to compute short range vdw interactions
!     n_sendshort2  !<number of MPI process to send positions to to compute short range vdw interactions
!
!     nrec_recep  !<number of MPI process to receive positions from to compute reciprocal interactions
!     (recip-recip communications)
!     nrec_send  !<number of MPI process to send positions to to compute reciprocal interactions
!     (recip-recip communications)
!     nrec_recep1  !<number of MPI process to receive positions from to compute reciprocal interactions
!     polarization only, no torques (recip-recip communications)
!     nrec_send1  !<number of MPI process to send positions to to compute reciprocal interactions
!     polarization only, no torques (recip-recip communications)
!     nrecdir_recep  !<number of MPI process to receive positions from to compute reciprocal interactions
!     (recip-direct communications)
!     nrecdir_send  !<number of MPI process to send positions to to compute reciprocal interactions
!     (recip-direct communications)
!     nrecdir_recep2  !<number of MPI process to receive positions from to compute reciprocal interactions, without proc already in precdir_recep1
!     (recip-direct communications)
!     nrecdir_send2  !<number of MPI process to send positions to to compute reciprocal interactions, without proc already in precdir_send1
!     (recip-direct communications)
!     (recip-direct communications)
!     (recip-direct communications)
!
!
!     nbig_recep  !<number of MPI process to receive positions from to compute largest non bonded interactions
!     nbig_send  !<number of MPI process to send positions to to compute largest non bonded interactions
!     nbigshort_recep  !<number of MPI process to receive positions from to compute largest short range non bonded interactions
!     nbigshort_send  !<number of MPI process to send positions to to compute largest short range non bonded interactions
!     nneig_recep  !<number of MPI process to receive positions from to compute bonded interactions
!     nneig_send  !<number of MPI process to send positions to to compute bonded interactions
!
!     p*_recep*  !<list of the corresponding processes
!     p*_send*  !<list of the corresponding processes
!
!     nloc  !<local number of atoms
!     nbloc !<local + neighbors number of atoms
!     nblocloop !<local + neighbors number of atoms : nbloc if nbloc is a multiple of 16, or the first one greater
!     nlocrec  !<local reciprocal number of atoms
!     nlocrec2  !<local + reciprocal neighbors number of atoms
!     nlocnl !<local nl number of atoms
!     nblocrecdir !<local + neighbors direct+reciprocal number of atoms
!
!     domlen !<number of atoms in the domains
!     domlenrec !<number of reciprocal atoms in the reciprocal domains
!     domlen !<number of multipoles in the domains
!     domlenrec !<number of reciprocal multipoles in the reciprocal domains
!
!     glob !<local-global correspondance
!     loc !<global-local correspondance
!     globrec !<local-global reciprocal correspondance
!     locrec !<global-local reciprocal correspondance
!     repart !<global index-domain correspondance
!     repartrec !<global index-domain reciprocal correspondance
!
!     bufbeg* !<index of the first atom concerned by each process
!     buflen1,buflen2,buf1,buf2,bufbeg1,bufbeg2 explicit direct-reciprocal atomic correspondance,
!     for polarization solvers :
!      - buflen* : number of atoms involved
!      - bufbeg* : index of the first atom concerned by each process
!      - buf*    : global index of the atoms involved
!
!     nx_box : !<size of each subdomain along the x-axis
!     ny_box : !<size of each subdomain along the y-axis
!     nz_box : !<size of each subdomain along the z-axis
!     scala : !<factor for triclinic domain decomposition along the a-axis
!     scalb : !<factor for triclinic domain decomposition along the b-axis
!     scalc : !<factor for triclinic domain decomposition along the c-axis
!     xbegproc : !<x coordinate of the beginning of each domain
!     ybegproc : !<y coordinate of the beginning of each domain
!     zbegproc : !<z coordinate of the beginning of each domain
!     xendproc : !<x coordinate of the ending of each domain
!     yendproc : !<y coordinate of the ending of each domain
!     zendproc : !<z coordinate of the ending of each domain
!     nxdd,nydd,nzdd : number of divisions along the axes, for domain decomposition
!     abegproc : !<a coordinate of the beginning of each domain
!     bbegproc : !<b coordinate of the beginning of each domain
!     cbegproc : !<c coordinate of the beginning of each domain
!     aendproc : !<a coordinate of the ending of each domain
!     bendproc : !<b coordinate of the ending of each domain
!     cendproc : !<c coordinate of the ending of each domain
!     nadd,nbdd,ncdd : number of divisions along the axes, for domain decomposition
!
module domdec
   implicit none
   integer :: nxdd !<number of subdivisions along the x axis
   integer :: nydd !<number of subdivisions along the y axis
   integer :: nzdd !<number of subdivisions along the z axis
   integer :: nproctot !<total number of MPI process (within MPI_COMM_WORLD)
   integer :: ranktot !<total rank of the MPI process within MPI_COMM_WORLD
   integer, target :: COMM_TINKER !<local MPI communicator in which a dynamic, analyze, testgrad or minimize run will take place
   integer, target :: nproc !<number of MPI processes during a dynamic, analyze, testgrad or minimize run
   integer, target :: rank !<rank of the current MPI process within COMM_TINKER
   integer :: rank_bis !<rank of the current MPI process within comm_dir or comm_rec
   integer :: nthread !<number of threads to be used with OpenMP
   integer :: nrec !<number of processes assigned to the computation of reciprocal space contribution
   integer :: ndir !<number of processes assigned to the computation of direct space contribution
   integer :: comm_rec !<MPI group communicator associated to the reciprocal space
   integer :: comm_dir !<MPI group communicator associated to the direct space
   integer :: COMM_POLYMER !<MPI group communicator associated to the PIMD parallelism
   integer :: rank_polymer !<MPI rank associated to COMM_POLYMER
   integer :: nproc_polymer !<number of procs within COMM_POLYMER
   integer :: hostrank !<MPI rank within each node
   integer :: hostcomm !<MPI communicator within each node
   integer :: n_recep1 !<number of MPI process to receive positions from to compute electrostatic interactions
   integer :: n_send1 !<number of MPI process to send positions to to compute electrostatic interactions
   integer :: nrec_recep !<number of MPI process to receive positions from to compute reciprocal interactions
   integer :: nrec_send !<number of MPI process to send positions to to compute reciprocal interactions
   integer :: n_recep2 !<number of MPI process to receive positions from to compute vdw interactions
   integer :: n_send2 !<number of MPI process to send positions to to compute vdw interactions
   integer :: nrecdir_recep !<number of MPI process to receive positions from to compute reciprocal interactions (recip-direct communications)
   integer :: nrecdir_send !<number of MPI process to send positions to to compute reciprocal interactions (recip-direct communications)
   integer :: nrecdir_recep2 !<number of MPI process to receive positions from to compute reciprocal interactions, without proc already in precdir_recep1
   integer :: nrecdir_send2 !<number of MPI process to send positions to to compute reciprocal interactions, without proc already in precdir_send1
   integer :: n_recepshort1 !<number of MPI process to receive positions from to compute short range electrostatic interactions
   integer :: n_sendshort1 !<number of MPI process to send positions to to compute short range electrostatic interactions
   integer :: n_recepshort2 !<number of MPI process to receive positions from to compute short range vdw interactions
   integer :: n_sendshort2 !<number of MPI process to send positions to to compute short range vdw interactions
   integer :: nneig_recep !<number of MPI process to receive positions from to compute bonded interactions
   integer :: nneig_send !<number of MPI process to send positions to to compute bonded interactions
   integer :: nrecdir_recep1 !<number of MPI process to receive positions from to compute reciprocal interactions
   integer :: nrecdir_send1 !<number of MPI process to send positions to to compute reciprocal interactions
   integer :: nbig_recep !<number of MPI process to receive positions from to compute largest non bonded interactions
   integer :: nbig_send !<number of MPI process to send positions to to compute largest non bonded interactions
   integer :: nbigshort_recep !<number of MPI process to receive positions from to compute largest short range non bonded interactions
   integer :: nbigshort_send !<number of MPI process to send positions to to compute largest short range non bonded interactions
   integer :: nbloc !<local + neighbors number of atoms
   integer :: nloc !<local number of atoms
   integer :: nlocrec !<local reciprocal number of atoms
   integer :: nlocrec2 !<local + reciprocal neighbors number of atoms
   integer :: nlocnl !<local nl number of atoms
   integer :: nblocrecdir !<local + neighbors direct+reciprocal number of atoms
   integer :: nblocloop !<variable number of atoms
   integer, allocatable,target:: domlen(:) !<number of atoms in the domains
   integer, allocatable,target:: domlenrec(:) !<number of reciprocal atoms in the reciprocal domains
   integer, allocatable,target:: domlenpole(:) !<number of multipoles in the domains
   integer, allocatable,target:: domlenpolerec(:) !<number of multipoles in the reciprocal domains
   integer, allocatable,target:: p_recep1(:) !<list of processes associated with n_recep1
   integer, allocatable,target:: p_send1(:) !<list of processes associated with n_send1
   integer, allocatable,target:: p_recep2(:) !<list of processes associated with n_recep2
   integer, allocatable,target:: p_send2(:) !<list of processes associated with n_send2
   integer, allocatable,target:: p_recepshort1(:) !<list of processes associated with n_recepshort1
   integer, allocatable,target:: p_sendshort1(:) !<list of processes associated with n_sendshort1
   integer, allocatable,target:: p_recepshort2(:) !<list of processes associated with n_recepshort2
   integer, allocatable,target:: p_sendshort2(:) !<list of processes associated with n_sendshort2
   integer, allocatable,target:: pneig_recep(:) !<list of processes associated with nneig_recep
   integer, allocatable,target:: pneig_send(:) !<list of processes associated with nneig_send
   integer, allocatable,target:: precdir_recep(:) !<list of processes associated with nrecdir_recep
   integer, allocatable,target:: precdir_send(:) !<list of processes associated with nrecdir_send
   integer, allocatable,target:: precdir_recep1(:) !<list of processes associated with nrecdir_recep1
   integer, allocatable,target:: precdir_send1(:) !<list of processes associated with nrecdir_send1
   integer, allocatable,target:: precdir_recep2(:) !<list of processes associated with nrecdir_recep2
   integer, allocatable,target:: precdir_send2(:) !<list of processes associated with nrecdir_send2
   integer, allocatable,target:: pbig_recep(:) !<list of processes associated with nbig_recep
   integer, allocatable,target:: pbig_send(:) !<list of processes associated with nbig_send
   integer, allocatable,target:: pbigshort_recep(:) !<list of processes associated with nbigshort_recep
   integer, allocatable,target:: pbigshort_send(:) !<list of processes associated with nbigshort_send
   integer, allocatable,target:: glob(:) !<local-global correspondance
   integer, allocatable,target:: loc(:) !<global-local correspondance
   integer, allocatable,target:: globrec(:) !<local-global reciprocal correspondance
   integer, allocatable,target:: locrec(:) !<global-local reciprocal correspondance
   integer, allocatable,target:: prec_send(:) !<list of processes associated with nrec_send
   integer, allocatable,target:: prec_recep(:) !<list of processes associated with nrec_recep
   integer, allocatable,target:: repartrec(:) !<global index-domain reciprocal correspondance
   integer, allocatable,target:: repart(:) !<global index-domain correspondance
   integer, allocatable,target:: bufbeg(:) !<index of the first atom concerned by each process for direct reciprocal atomic correspondance for polarization solvers
   integer, allocatable,target:: bufbegpole(:) !<index of the first multipole concerned by each process for direct reciprocal atomic correspondance for polarization solvers
   integer, allocatable,target:: bufbegrec(:) !<index of the first reciprocal multipole concerned by each process for direct reciprocal atomic correspondance for polarization solvers
   integer, allocatable,target:: buflen1(:) !<number of atoms concerned by each process for direct reciprocal atomic correspondance for polarization solvers
   integer, allocatable,target:: buflen2(:) !<number of atoms concerned by each process for direct reciprocal atomic correspondance for polarization solvers
   integer, allocatable,target:: bufbeg1(:) !<index of the first atom concerned by each process for direct reciprocal atomic correspondance for polarization solvers
   integer, allocatable,target:: bufbeg2(:) !<index of the first atom concerned by each process for direct reciprocal atomic correspondance for polarization solvers
   integer, allocatable,target:: buf1(:) !<list of atoms concerned by each process for direct reciprocal atomic correspondance for polarization solvers
   integer, allocatable,target:: buf2(:) !<list of atoms concerned by each process for direct reciprocal atomic correspondance for polarization solvers
   real*8 :: na_box !<size of each subdomain along the a-axis
   real*8 :: nb_box !<size of each subdomain along the b-axis
   real*8 :: nc_box !<size of each subdomain along the c-axis
   real*8 :: scala !<factor for triclinic domain decomposition along the a-axis
   real*8 :: scalb !<factor for triclinic domain decomposition along the b-axis
   real*8 :: scalc !<factor for triclinic domain decomposition along the c-axis
   real*8, allocatable :: zbegproc(:) !<z coordinate of the beginning of each domain
   real*8, allocatable :: zendproc(:) !<z coordinate of the end of each domain
   real*8, allocatable :: ybegproc(:) !<y coordinate of the beginning of each domain
   real*8, allocatable :: yendproc(:) !<y coordinate of the end of each domain
   real*8, allocatable :: xbegproc(:) !<x coordinate of the beginning of each domain
   real*8, allocatable :: xendproc(:) !<x coordinate of the end of each domain
   real*8, allocatable :: abegproc(:) !<a coordinate of the beginning of each domain
   real*8, allocatable :: aendproc(:) !<a coordinate of the end of each domain
   real*8, allocatable :: bbegproc(:) !<b coordinate of the beginning of each domain
   real*8, allocatable :: bendproc(:) !<b coordinate of the end of each domain
   real*8, allocatable :: cbegproc(:) !<c coordinate of the beginning of each domain
   real*8, allocatable :: cendproc(:) !<c coordinate of the end of each domain
   save
end
