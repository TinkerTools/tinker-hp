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
c     nproctot  total number of MPI process (within MPI_COMM_WORLD)
c     ranktot  total rank of the MPI process within MPI_COMM_WORLD 
c     nxdd = number of subdivisions along the x axis
c     nydd = number of subdivisions along the y axis
c     nzdd = number of subdivisions along the z axis
c
c     COMM_TINKER local MPI communicator in which a dynamic, analyze, testgrad or minimize run
c      will take place
c     nproc     number of MPI processes during a dynamic, analyze, testgrad or minimize run
c     rank      rank of the current MPI process within COMM_TINKER
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
c     n_recepshort1  number of MPI process to receive positions from to compute short range electrostatic interactions
c     n_sendshort1  number of MPI process to send positions to to compute short range electrostatic interactions
c     n_recepshort2  number of MPI process to receive positions from to compute short range vdw interactions
c     n_sendshort2  number of MPI process to send positions to to compute short range vdw interactions
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
c     nrecdir_recep2  number of MPI process to receive positions from to compute reciprocal interactions, without proc already in precdir_recep1
c     (recip-direct communications)
c     nrecdir_send2  number of MPI process to send positions to to compute reciprocal interactions, without proc already in precdir_send1
c     (recip-direct communications)
c
c
c     ntorque_recep  number of MPI process to receive positions from to compute electrostatic interactions + associated torques
c     ntorque_send  number of MPI process to send positions to to compute electrostatic interactions + associated torques
c     nbig_recep  number of MPI process to receive positions from to compute largest non bonded interactions
c     nbig_send  number of MPI process to send positions to to compute largest non bonded interactions
c     nbigshort_recep  number of MPI process to receive positions from to compute largest short range non bonded interactions
c     nbigshort_send  number of MPI process to send positions to to compute largest short range non bonded interactions
c     ntorqueshort_recep  number of MPI process to receive positions from to compute electrostatic interactions + associated torques
c     ntorqueshort_send  number of MPI process to send positions to to compute electrostatic interactions + associated torques
c     nneig_recep  number of MPI process to receive positions from to compute bonded interactions
c     nneig_send  number of MPI process to send positions to to compute bonded interactions
c
c     p*_recep*  list of the corresponding processes
c     p*_send*  list of the corresponding processes
c
c     nloc  local number of atoms
c     nbloc local + neighbors number of atoms
c     nlocrec  local reciprocal number of atoms
c     nlocrec2  local + reciprocal neighbors number of atoms
c     nlocnl local nl number of atoms
c     nlocnlb     first multiple of BLOCK_SIZE after nlocnl
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
c     nxdd,nydd,nzdd : number of divisions along the axes, for domain decomposition
c
#include "tinker_macro.h"
      module domdec
      implicit none
      logical Bdecomp1d,Bdecomp2d,Bdecomp3d
      integer,parameter:: masterRank=0
      integer nxdd,nydd,nzdd
      integer nproctot,ranktot
      integer,target:: COMM_TINKER,nproc,rank
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

      integer,allocatable,target:: domlen(:),domlenrec(:)
     &       , domlenpole(:),domlenpolerec(:)
      integer,allocatable,target:: p_recep1(:),p_send1(:)
     &       , p_recep2(:),p_send2(:),p_recepshort1(:),p_sendshort1(:)
     &       , p_recepshort2(:),p_sendshort2(:)
      integer,allocatable,target:: pneig_recep(:), pneig_send(:)
     &       , precdir_recep(:), precdir_send(:)
     &       , precdir_recep1(:), precdir_send1(:)
     &       , precdir_recep2(:), precdir_send2(:)
      integer,allocatable,target:: ptorque_recep(:), ptorque_send(:)
     &       , ptorqueshort_recep(:),ptorqueshort_send(:)
      integer,allocatable,target:: pbig_recep(:), pbig_send(:)
     &       , pbigshort_recep(:), pbigshort_send(:)
     &       , prec_send(:), prec_recep(:)
      integer,allocatable,target:: glob(:),loc(:),globrec(:),locrec(:)
     &       , globrec1(:),locrec1(:),repartrec(:),repart(:)
      integer,allocatable,target:: bufbeg(:),bufbegpole(:),bufbegrec(:)
     &       , bufbeg1(:),bufbeg2(:)
      integer,allocatable :: buf1(:),buf2(:)
      integer,allocatable :: buflen1(:),buflen2(:)

      real(t_p) nx_box,ny_box,nz_box
!DIR$ ATTRIBUTES ALIGN:64:: xbegproc,xendproc,ybegproc,yendproc
!DIR$ ATTRIBUTES ALIGN:64:: zbegproc,zendproc
      real(t_p), allocatable :: xbegproc(:),xendproc(:)
     &         , ybegproc(:),yendproc(:),zbegproc(:),zendproc(:)

!$acc declare create(rank,rank_bis)
!$acc declare create(nrec_send)
!$acc declare create(p_recep1,p_recep2)
!$acc declare create(prec_send,prec_recep)
!$acc declare create(bufbeg1,bufbeg2)
!$acc declare create(xbegproc,xendproc)
!$acc declare create(ybegproc,yendproc)
!$acc declare create(zbegproc,zendproc)

      end
