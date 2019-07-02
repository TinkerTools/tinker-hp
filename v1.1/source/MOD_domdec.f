c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################################################
c     ##                                                                                        ##
c     ##  module domdec  --  system parameters for OpenMP/MPI domain decomposition computation  ##
c     ##                                                                                        ##
c     ############################################################################################
c
c
c     nproc     number of MPI processes
c     rank      rank of the current MPI process within MPI_COMM_WORLD
c     rank_bis  rank of the current MPI process within comm_dir or comm_rec
c     nrec      number of processes assigned to the computation of reciprocal space contribution
c     ndir      number of processes assigned to the computation of direct space contribution
c     comm_rec  MPI group communicator associated to the reciprocal space
c     comm_dir  MPI group communicator associated to the direct space
c     nthread   number of threads to be used with OpenMP
c     hostcomm  MPI group communicator associated to processes within a node
c     hostrank rank of the current MPI process within hostcomm
c
c     n_recep1  number of MPI process to receive positions from to compute electrostatic interactions
c     n_send1  number of MPI process to send positions to to compute electrostatic interactions
c     n_recep2  number of MPI process to receive positions from to compute vdw interactions
c     n_send2  number of MPI process to send positions to to compute vdw interactions
c
c     nrec_recep  number of MPI process to receive positions from to compute reciprocal interactions
c     (recip-recip communications)
c     nrec_send  number of MPI process to send positions to to compute reciprocal interactions
c     (recip-recip communications)
c     nrec_recep1  number of MPI process to receive positions from to compute reciprocal interactions
c     polarization only, no torques (recip-recip communications)
c     nrec_send1  number of MPI process to send positions to to compute reciprocal interactions
c     polarization only, no torques (recip-recip communications)
c     nrecdir_recep  number of MPI process to receive positions from to compute reciprocal interactions
c     (recip-direct communications)
c     nrecdir_send  number of MPI process to send positions to to compute reciprocal interactions
c     (recip-direct communications)
c
c
c     ntorque_recep  number of MPI process to receive positions from to compute electrostatic interactions + associated torques
c     ntorque_send  number of MPI process to send positions to to compute electrostatic interactions + associated torques
c     nbig_recep  number of MPI process to receive positions from to compute largest non bonded interactions
c     nbig_send  number of MPI process to send positions to to compute largest non bonded interactions
c     nneig_recep  number of MPI process to receive positions from to compute bonded interactions
c     nneig_send  number of MPI process to send positions to to compute bonded interactions
c
c     p*_recep*  list of the corresponding processes
c     p*_send*  list of the corresponding processes
c
c     nloc  local number of atoms
c     nbloc local + neighbors number of atoms
c     nlocrec  local reciprocal number of atoms
c     nblocrec local + neighbors reciprocal number of atoms
c     nlocnl local nl number of atoms
c     nblocrecdir local + neighbors direct+reciprocal number of atoms
c
c     domlen number of atoms in the domains
c     domlenrec number of reciprocal atoms in the reciprocal domains
c     domlen number of multipoles in the domains
c     domlenrec number of reciprocal multipoles in the reciprocal domains
c
c     glob local-global correspondance
c     loc global-local correspondance
c     globrec local-global reciprocal correspondance
c     locrec global-local reciprocal correspondance
c     globrec1 local-global direct-reciprocal correspondance
c     locrec1 global-local direct-reciprocal correspondance
c     repart global index-domain correspondance
c     repartrec global index-domain reciprocal correspondance
c
c     bufbeg* index of the first atom concerned by each process
c     buflen1,buflen2,buf1,buf2,bufbeg1,bufbeg2 explicit direct-reciprocal atomic correspondance, 
c     for polarization solvers :
c      - buflen* : number of atoms involved
c      - bufbeg* : index of the first atom concerned by each process
c      - buf*    : global index of the atoms involved
c
c     nx_box : size of each subdomain along the x-axis
c     ny_box : size of each subdomain along the y-axis
c     nz_box : size of each subdomain along the z-axis
c     xbegproc : x coordinate of the beginning of each domain
c     ybegproc : y coordinate of the beginning of each domain
c     zbegproc : z coordinate of the beginning of each domain
c     xendproc : x coordinate of the ending of each domain
c     yendproc : y coordinate of the ending of each domain
c     zendproc : z coordinate of the ending of each domain
c
      module domdec
      implicit none
      integer nproc,rank,rank_bis,nthread,nrec,ndir,comm_rec,comm_dir
      integer hostrank,hostcomm
      integer n_recep1, n_send1, nrec_recep,nrec_send
      integer n_recep2, n_send2, nrecdir_recep,nrecdir_send
      integer ntorque_recep,ntorque_send
      integer nneig_recep,nneig_send
      integer nrecdir_recep1,nrecdir_send1
      integer nbig_recep,nbig_send
      integer nbloc,nloc,nlocrec,nblocrec
      integer nlocnl,nblocrecdir
      integer, allocatable :: domlen(:), domlenrec(:)
      integer, allocatable :: domlenpole(:), domlenpolerec(:)
      integer, allocatable :: p_recep1(:), p_send1(:)
      integer, allocatable :: p_recep2(:), p_send2(:)
      integer, allocatable :: pneig_recep(:), pneig_send(:)
      integer, allocatable :: precdir_recep(:), precdir_send(:)
      integer, allocatable :: precdir_recep1(:), precdir_send1(:)
      integer, allocatable :: ptorque_recep(:), ptorque_send(:)
      integer, allocatable :: pbig_recep(:), pbig_send(:)
      integer, allocatable :: glob(:),loc(:)
      integer, allocatable :: globrec(:), locrec(:)
      integer, allocatable :: globrec1(:), locrec1(:)
      integer, allocatable :: prec_send(:), prec_recep(:)
      integer, allocatable :: repartrec(:),repart(:)
      integer, allocatable :: bufbeg(:),bufbegpole(:)
      integer, allocatable :: bufbegrec(:)
      integer, allocatable :: buflen1(:),buflen2(:)
      integer, allocatable :: bufbeg1(:),bufbeg2(:)
      integer, allocatable :: buf1(:),buf2(:)
      real*8 nx_box,ny_box,nz_box
      real*8, allocatable :: zbegproc(:),zendproc(:)
      real*8, allocatable :: ybegproc(:),yendproc(:)
      real*8, allocatable :: xbegproc(:),xendproc(:)
      save
      end
