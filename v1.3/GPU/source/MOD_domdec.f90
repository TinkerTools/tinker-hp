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
!     nproctot  total number of MPI process (within MPI_COMM_WORLD)
!     ranktot  total rank of the MPI process within MPI_COMM_WORLD
!     nxdd = number of subdivisions along the x axis
!     nydd = number of subdivisions along the y axis
!     nzdd = number of subdivisions along the z axis
!
!     COMM_TINKER local MPI communicator in which a dynamic, analyze, testgrad or minimize run
!      will take place
!     nproc     number of MPI processes during a dynamic, analyze, testgrad or minimize run
!     rank      rank of the current MPI process within COMM_TINKER
!     rank_bis  rank of the current MPI process within comm_dir or comm_rec
!     nrec      number of processes assigned to the computation of reciprocal space contribution
!     ndir      number of processes assigned to the computation of direct space contribution
!     comm_rec  MPI group communicator associated to the reciprocal space
!     comm_dir  MPI group communicator associated to the direct space
!     nthread   number of threads to be used with OpenMP
!     hostcomm  MPI group communicator associated to processes within a node
!     hostrank rank of the current MPI process within hostcomm
!
!     n_recep1  number of MPI process to receive positions from to compute electrostatic interactions
!     n_send1  number of MPI process to send positions to to compute electrostatic interactions
!     n_recep2  number of MPI process to receive positions from to compute vdw interactions
!     n_send2  number of MPI process to send positions to to compute vdw interactions
!
!     n_recepshort1  number of MPI process to receive positions from to compute short range electrostatic interactions
!     n_sendshort1  number of MPI process to send positions to to compute short range electrostatic interactions
!     n_recepshort2  number of MPI process to receive positions from to compute short range vdw interactions
!     n_sendshort2  number of MPI process to send positions to to compute short range vdw interactions
!
!     nrec_recep  number of MPI process to receive positions from to compute reciprocal interactions
!     (recip-recip communications)
!     nrec_send  number of MPI process to send positions to to compute reciprocal interactions
!     (recip-recip communications)
!     nrec_recep1  number of MPI process to receive positions from to compute reciprocal interactions
!     polarization only, no torques (recip-recip communications)
!     nrec_send1  number of MPI process to send positions to to compute reciprocal interactions
!     polarization only, no torques (recip-recip communications)
!     nrecdir_recep  number of MPI process to receive positions from to compute reciprocal interactions
!     (recip-direct communications)
!     nrecdir_send  number of MPI process to send positions to to compute reciprocal interactions
!     (recip-direct communications)
!     nrecdir_recep2  number of MPI process to receive positions from to compute reciprocal interactions, without proc already in precdir_recep1
!     (recip-direct communications)
!     nrecdir_send2  number of MPI process to send positions to to compute reciprocal interactions, without proc already in precdir_send1
!     (recip-direct communications)
!
!
!     ntorque_recep  number of MPI process to receive positions from to compute electrostatic interactions + associated torques
!     ntorque_send  number of MPI process to send positions to to compute electrostatic interactions + associated torques
!     nbig_recep  number of MPI process to receive positions from to compute largest non bonded interactions
!     nbig_send  number of MPI process to send positions to to compute largest non bonded interactions
!     nbigshort_recep  number of MPI process to receive positions from to compute largest short range non bonded interactions
!     nbigshort_send  number of MPI process to send positions to to compute largest short range non bonded interactions
!     ntorqueshort_recep  number of MPI process to receive positions from to compute electrostatic interactions + associated torques
!     ntorqueshort_send  number of MPI process to send positions to to compute electrostatic interactions + associated torques
!     nneig_recep  number of MPI process to receive positions from to compute bonded interactions
!     nneig_send  number of MPI process to send positions to to compute bonded interactions
!
!     p*_recep*  list of the corresponding processes
!     p*_send*  list of the corresponding processes
!
!     nloc  local number of atoms
!     nbloc local + neighbors number of atoms
!     nlocrec  local reciprocal number of atoms
!     nlocrec2  local + reciprocal neighbors number of atoms
!     nlocnl local nl number of atoms
!     nlocnlb     first multiple of BLOCK_SIZE after nlocnl
!     nblocrecdir local + neighbors direct+reciprocal number of atoms
!
!     domlen number of atoms in the domains
!     domlenrec number of reciprocal atoms in the reciprocal domains
!     domlen number of multipoles in the domains
!     domlenrec number of reciprocal multipoles in the reciprocal domains
!
!     glob local-global correspondance
!     loc global-local correspondance
!     globrec local-global reciprocal correspondance
!     locrec global-local reciprocal correspondance
!     repart global index-domain correspondance
!     repartrec global index-domain reciprocal correspondance
!
!     bufbeg* index of the first atom concerned by each process
!     buflen1,buflen2,buf1,buf2,bufbeg1,bufbeg2 explicit direct-reciprocal atomic correspondance,
!     for polarization solvers :
!      - buflen* : number of atoms involved
!      - bufbeg* : index of the first atom concerned by each process
!      - buf*    : global index of the atoms involved
!
!     nx_box : size of each subdomain along the x-axis
!     ny_box : size of each subdomain along the y-axis
!     nz_box : size of each subdomain along the z-axis
!     scala : factor for triclinic domain decomposition along the a-axis
!     scalb : factor for triclinic domain decomposition along the b-axis
!     scalc : factor for triclinic domain decomposition along the c-axis
!     xbegproc : x coordinate of the beginning of each domain
!     ybegproc : y coordinate of the beginning of each domain
!     zbegproc : z coordinate of the beginning of each domain
!     xendproc : x coordinate of the ending of each domain
!     yendproc : y coordinate of the ending of each domain
!     zendproc : z coordinate of the ending of each domain
!     abegproc : a coordinate of the beginning of each domain
!     bbegproc : b coordinate of the beginning of each domain
!     cbegproc : c coordinate of the beginning of each domain
!     aendproc : a coordinate of the ending of each domain
!     bendproc : b coordinate of the ending of each domain
!     cendproc : c coordinate of the ending of each domain
!     nxdd,nydd,nzdd : number of divisions along the axes, for domain decomposition
!
#include "tinker_macro.h"
module domdec
   implicit none
   logical Bdecomp1d,Bdecomp2d,Bdecomp3d
   integer,parameter:: masterRank=0
   integer nxdd,nydd,nzdd
   integer nproctot,ranktot
   integer,target:: COMM_TINKER,nproc,rank
   integer nproc_polymer,rank_polymer,COMM_POLYMER
   integer rank_bis,nthread,nrec,ndir,comm_rec,comm_dir
   integer hostrank,hostcomm
   integer n_recep1, n_send1, nrec_recep,nrec_send
   integer n_recep2, n_send2, nrecdir_recep,nrecdir_send
   integer nrecdir_recep2,nrecdir_send2
   integer nrecdir_recep3,nrecdir_send3
   integer n_recepshort1,n_sendshort1,n_recepshort2,n_sendshort2
   integer ntorque_recep,ntorque_send
   integer ntorqueshort_recep,ntorqueshort_send
   integer nneig_recep,nneig_send
   integer nrecdir_recep1,nrecdir_send1
   integer nbig_recep,nbig_send
   integer nbigshort_recep,nbigshort_send
   integer nbloc,nloc,nlocrec,nlocrec2
   integer nlocnl,nblocrecdir
   integer nlocnlb
   integer nblocloop

   integer,allocatable,target:: domlen(:),domlenrec(:)&
          , domlenpole(:),domlenpolerec(:)
   integer,allocatable,target:: p_recep1(:),p_send1(:)&
          , p_recep2(:),p_send2(:),p_recepshort1(:),p_sendshort1(:)&
          , p_recepshort2(:),p_sendshort2(:)
   integer,allocatable,target:: pneig_recep(:), pneig_send(:)&
          , precdir_recep(:), precdir_send(:)&
          , precdir_recep1(:), precdir_send1(:)&
          , precdir_recep2(:), precdir_send2(:)
   integer,allocatable,target:: ptorque_recep(:), ptorque_send(:)&
          , ptorqueshort_recep(:),ptorqueshort_send(:)
   integer,allocatable,target:: pbig_recep(:), pbig_send(:)&
          , pbigshort_recep(:), pbigshort_send(:)&
          , prec_send(:), prec_recep(:)
   integer,allocatable,target:: glob(:),loc(:),globrec(:),locrec(:)&
          , globrec1(:),locrec1(:),repartrec(:),repart(:)
   integer,allocatable,target:: bufbeg(:),bufbegpole(:),bufbegrec(:)&
          , bufbeg1(:),bufbeg2(:)
   integer,allocatable :: buf1(:), buf2(:), buflen1(:), buflen2(:)

   real(r_p) nx_box,ny_box,nz_box
   real(r_p) na_box,nb_box,nc_box
   real(r_p) scala,scalb,scalc
   real(t_p), allocatable :: xbegproc(:),xendproc(:)&
      &, ybegproc(:),yendproc(:),zbegproc(:),zendproc(:)
   real(t_p), allocatable :: abegproc(:),aendproc(:)&
            , bbegproc(:),bendproc(:)&
            , cbegproc(:),cendproc(:)

!$acc declare create(rank,rank_bis)
!$acc declare create(prec_send)
!$acc declare create(bufbeg1,bufbeg2)
!$acc declare create(xbegproc,xendproc,ybegproc,yendproc,zbegproc,zendproc)

end
