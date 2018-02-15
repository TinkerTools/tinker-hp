c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  openmp.i  --  system parameters for OpenMP/MPI computation  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     nproc     number of MPI processes
c     rank      rank of the current MPI process within MPI_COMM_WORLD
c     rank_bis  rank of the current MPI process within comm_dir or comm_rec
c     nrec      number of processes assigned to the computation of reciprocal space
c      contribution
c     ndir      number of processes assigned to the computation of direct space
c      contribution
c     nthread   number of threads to be used with OpenMP
c     comm_rec   MPI group communicator associated to the reciprocal space
c     comm_dir   MPI group communicator associated to the direct space
c
c
      integer nproc,rank,rank_bis,nthread,nrec,ndir,comm_rec,comm_dir
      integer hostrank,hostcomm
      integer n_recep1, n_send1, nrec_recep,nrec_send
      integer n_recep2, n_send2, nrecdir_recep,nrecdir_send
      integer nbloc,nloc,nlocold,nblocold,nlocrec,nblocrec
      integer nlocnl,nmoleloc,nblocrecdir
      integer ntorque_recep,ntorque_send
      integer nneig_recep,nneig_send
      integer nrecdir_recep1,nrecdir_send1
      integer nbig_recep,nbig_send
      real*8, pointer :: zbegproc(:),zendproc(:)
      real*8, pointer :: ybegproc(:),yendproc(:)
      real*8, pointer :: xbegproc(:),xendproc(:)
      integer, pointer :: domlen(:), domlenrec(:)
      integer, pointer :: p_recep1(:), p_send1(:)
      integer, pointer :: p_recep2(:), p_send2(:)
      integer, pointer :: pneig_recep(:), pneig_send(:)
      integer, pointer :: precdir_recep(:), precdir_send(:)
      integer, pointer :: precdir_recep1(:), precdir_send1(:)
      integer, pointer :: ptorque_recep(:), ptorque_send(:)
      integer, pointer :: pbig_recep(:), pbig_send(:)
      integer, pointer :: glob(:)
      integer, pointer :: loc(:)
      integer, pointer :: globrec(:), locrec(:)
      integer, pointer :: prec_send(:), prec_recep(:)
      integer, pointer :: repartrec(:),repart(:)
      integer, pointer :: bufbeg(:),bufbegpole(:)
      integer, pointer :: bufbegrec(:),locrec1(:),globrec1(:)
      integer, pointer :: buflen1(:),buflen2(:),bufbeg1(:),bufbeg2(:)
      integer, pointer :: buf1(:),buf2(:)
      real*8 nx_box,ny_box,nz_box
      common /openmp/ nx_box,ny_box,nz_box,
     $  nproc,nthread,rank,nrecdir_recep,nrecdir_send,
     $  hostrank,hostcomm,
     $  rank_bis,nrec,comm_rec,comm_dir,
     $  ndir,nlocnl,
     $  p_recep1,p_send1,n_recep1,n_send1,
     $  p_recep2,p_send2,n_recep2,n_send2,
     $  nrec_send,nrec_recep,prec_recep,prec_send,
     $  repart,repartrec,bufbeg,buf1,buf2,bufbeg1,bufbeg2,
     $  buflen1,buflen2,nbloc,glob,loc,
     $  zbegproc,zendproc,domlen,domlenrec,globrec,
     $  locrec,bufbegpole,nloc,nlocrec,nlocold,nmoleloc,
     $  ntorque_recep,ntorque_send,ptorque_recep,ptorque_send,
     $  nneig_recep,nneig_send,pneig_recep,pneig_send,
     $  precdir_recep,precdir_send,xbegproc,xendproc,ybegproc,
     $  yendproc,nblocrec,precdir_recep1,precdir_send1,nrecdir_recep1,
     $  nrecdir_send1,nbig_recep,nbig_send,pbig_recep,pbig_send,
     $  bufbegrec,locrec1,globrec1,nblocrecdir
