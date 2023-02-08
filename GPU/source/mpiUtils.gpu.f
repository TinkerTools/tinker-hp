c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################
c     ##                                                  ##
c     ##  mpiUtils subModule  -- Tinker-HP                ##
c     ##                                                  ##
c     ######################################################
c
c
#include "tinker_macro.h"
      submodule(utilcomm) mpiutils
      use domdec
      use inform
      use iso_c_binding
      use mpi
      use potent
      use timestat
      use tinheader
      use tinMemory
      use utilgpu
      use sizes    ,only: tinkerdebug
      implicit none
      private

      public :: commDDd_ext_c,commDDd_add_c,commDDrd_add_c

      contains
#include "convert.f.inc"

      subroutine setv_opt(opt_,n_loc,n_bloc,nsend,nrecv
     &               ,p_send,p_recv,p_beg,d_len,name)
      integer,intent(in ):: opt_
      character(*),intent(in) :: name
      integer,intent(out):: nsend,nrecv,n_loc,n_bloc
      integer,pointer:: p_send(:),p_recv(:),p_beg(:),d_len(:)

      if      (opt_.eq.ucBig) then
         n_loc     =  nloc
         n_bloc    =  nbloc
         nsend     =  nbig_send
         nrecv     =  nbig_recep
         p_send(1:nproc) => pbig_send (:)
         p_recv(1:nproc) => pbig_recep(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else if (opt_.eq.ucNeig) then
         n_loc     =  nloc
         n_bloc    =  nbloc
         nsend     =  nneig_send
         nrecv     =  nneig_recep
         p_send(1:nproc) => pneig_send (:)
         p_recv(1:nproc) => pneig_recep(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else if (opt_.eq.ucShort) then
         n_loc     =  nloc
         n_bloc    =  nbloc
         nsend     =  nbigshort_send
         nrecv     =  nbigshort_recep
         p_send(1:nproc) => pbigshort_send (:)
         p_recv(1:nproc) => pbigshort_recep(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else if (opt_.eq.ucRec) then
         n_loc     =  nlocrec
         n_bloc    =  nlocrec2
         nsend     =  nrec_send
         nrecv     =  nrec_recep
         p_send(1:nproc) => prec_send (:)
         p_recv(1:nproc) => prec_recep(:)
         p_beg (1:nproc) => bufbegrec(:)
         d_len (1:nproc) => domlenrec(:)
      else if (opt_.eq.ucDirRec) then
         n_loc     =  nloc
         n_bloc    =  nblocrecdir
         nsend     =  nrecdir_send2
         nrecv     =  nrecdir_recep2
         p_send(1:nproc) => precdir_send2 (:)
         p_recv(1:nproc) => precdir_recep2(:)
         p_beg (1:nproc) => bufbeg(:)
         d_len (1:nproc) => domlen(:)
      else
 16      format(' ERROR! routine ',A,
     &       ,/,' --- Available options ',3I3
     &       ,/,' ---   UNKNOWN OPTION  ',I3   )
         write(0,16) name,ucNeig,ucShort,ucBig,opt_
         __TINKER_FATAL__
      end if
      end subroutine

      module subroutine commDDd_ext_c(extb,dim,opt)
      implicit  none
      real(t_p),intent(inout):: extb(*)
      integer  ,optional     :: dim,opt

      integer   i,j,k,tag,ierr,siz,opt_,dim_,task
      integer   nsend,nrecv,psend,precv,pbeg,dlen,idxr,commloc
      integer   n_loc,n_bloc,nr_loc,nr_bloc
      integer   status(MPI_STATUS_SIZE)
      integer  ,pointer:: p_send(:),p_recv(:),p_beg(:),d_len(:)
      parameter (tag=0)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
c
      task = ucComm
      opt_ = ucBig
      dim_ = 3
      if (present(opt)) then
         task = opt
         opt_ = iand(opt,15)
      end if
      if (present(dim)) dim_=dim
c
      commloc = COMM_TINKER
      call setv_opt(opt_,n_loc,n_bloc,nsend,nrecv
     &         ,p_send,p_recv,p_beg,d_len,'commDDd_ext_c')
      if (nsend.eq.0.and.nrecv.eq.0) return
c
 21   format(3x,a,' dim(',I0,') opt&task(',2I5,')')
      if (deb_Path) write(*,21) '>> commDDd_ext_c',dim_,opt_,task
!$acc host_data use_device(extb)
c
c     MPI : begin reception
c
      if (iand(task,ucRecv).ne.0) then
!$acc wait(def_queue)
      do i = 1, nrecv
         precv = p_recv(i)
         dlen  = dim_* d_len(precv+1)
         idxr  = dim_*(p_beg(precv+1)-1)+1
         call MPI_IRECV(extb(idxr),dlen,MPI_TPREC,precv,tag
     &                 ,COMM_TINKER,reqr_dirdir(i),ierr)
      end do
      end if
c
c     MPI : begin sending
c
      if (iand(task,ucSend).ne.0) then
!$acc wait(def_queue)
      do i = 1, nsend
         call MPI_ISEND(extb,dim_*nloc,MPI_TPREC,p_send(i),tag
     &           ,COMM_TINKER,reqs_dirdir(i),ierr)
      end do
      end if
!$acc end host_data
c
c     MPI : begin wait
c
      if (iand(task,ucWait).ne.0) then
      do i= 1,nsend; call MPI_WAIT(reqs_dirdir(i),status,ierr); end do
      do i= 1,nrecv; call MPI_WAIT(reqr_dirdir(i),status,ierr); end do
 31   format(3x,a)
      if (btest(tinkerdebug,tindPath)) call MPI_BARRIER(hostcomm,ierr)
      if (deb_Path) write(*,31) '<< commDDd_ext_c',opt_
      end if

      end subroutine

c     module subroutine commDDd_ext_m
c     end subroutine

c     module subroutine commDDd_ext_f
c     end subroutine

      module subroutine commDDd_add_c(addb,dim,opt)
      implicit  none
      real(t_p),intent(inout):: addb(*)
      integer  ,optional     :: dim,opt

      integer   i,j,k,tag,ierr,siz,opt_,dim_,dimn,task
      integer   nsend,nrecv,psend,precv,pbeg,dlen,idxs,idxr,commloc
      integer   n_loc,n_bloc,nr_loc,nr_bloc
      integer   status(MPI_STATUS_SIZE)
      integer  ,pointer:: p_send(:),p_recv(:),p_beg(:),d_len(:)
      real(t_p),pointer:: buffer(:)
      type(c_ptr) cptr
      parameter (tag=0)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
c
      task = ucComm
      opt_ = ucBig
      dim_ = 3
      if (present(opt)) then
         task = opt
         opt_ = iand(opt,15)
      end if
      if (present(dim)) dim_=dim
c
      commloc = COMM_TINKER
      call setv_opt(opt_,n_loc,n_bloc,nsend,nrecv
     &         ,p_send,p_recv,p_beg,d_len,'commDDd_add_c')
      if (nsend.eq.0.and.nrecv.eq.0) return
c
 21   format(3x,a,' dim(',I0,') opt&task(',2I5,')')
      if (deb_Path) write(*,21) '>> commDDd_add_c',dim_,opt_,task
      ! Create exchange buffer
      siz = dim_*max(1,nloc)*nrecv
      call prmem_request(ucP0%buffMpi_r1,siz,async=.false.)
      cptr = c_loc(ucP0%buffMpi_r1)
      call c_f_pointer(cptr,buffer,[siz])
c
      !MPI : Start reception in buffer
!$acc host_data use_device(buffer,addb)
      if (iand(task,ucRecv).ne.0) then
!$acc wait(def_queue)
      do i = 1, nrecv
         precv = p_recv(i)
         idxr  = dim_*nloc*(i-1)+1
         dlen  = dim_*nloc
         call MPI_IRECV(buffer(idxr),dlen,MPI_TPREC,precv,tag
     &           ,COMM_TINKER,reqr_dird_a(i),ierr)
      end do
      end if
c
      !MPI : Start sending
      if (iand(task,ucSend).ne.0) then
!$acc wait(def_queue)
      do i = 1, nsend
         psend = p_send(i)
         idxs  = dim_*(p_beg(psend+1)-1)+1
         dlen  = dim_* d_len(psend+1)
         call MPI_ISEND(addb(idxs),dlen,MPI_TPREC,psend,tag
     &           ,COMM_TINKER,reqs_dird_a(i),ierr)
      end do
      end if
!$acc end host_data
c
      ! Wait for Exchanges to end
      if (iand(task,ucWait).ne.0) then
      dimn = dim_*nloc
      do i = 1,nrecv; call MPI_WAIT(reqr_dird_a(i),status,ierr); end do;
c
      !MPI : move in global arrays
!$acc parallel loop gang vector async(def_queue) present(addb,buffer)
      do j = 1,dim_*nloc; do i = 0,nsend-1;
         addb(j) = addb(j) + buffer(j+dimn*i)
      end do; end do
c
      do i = 1,nsend; call MPI_WAIT(reqs_dird_a(i),status,ierr); end do;
 31   format(3x,a)
      if (btest(tinkerdebug,tindPath)) call MPI_BARRIER(hostcomm,ierr)
      if (deb_Path) write(*,31) '<< commDDd_add_c'
      end if
c
      nullify(buffer)
      end subroutine

c     module subroutine commDDd_add_m
c     end subroutine

c     module subroutine commDDd_add_f
c     end subroutine

      module subroutine commDDrd_add_c(dirb,recb,dim,opt)
      implicit  none
      real(t_p),intent(inout):: dirb(*)
      real(t_p),intent(inout):: recb(*)
      integer  ,optional     :: dim,opt

      integer   i,j,k,tag,ierr,siz,siz1,siz2,rankloc,dim_,dimn,task
      integer   psend,precv,pbeg,dlen,idxs,idxr,commloc
      integer   jloc,jglob,jlocrec,iglob,iloc
      integer   status(MPI_STATUS_SIZE)
      integer  ,pointer:: p_send(:),p_recv(:),p_beg(:),d_len(:)
      real(t_p),pointer:: buffer(:),buffers(:)
      type(c_ptr) cptr0,cptr1
      parameter (tag=0)
c
      if (nproc.eq.1) then
!$acc parallel loop async(def_queue) default(present)
         do j = 1, nlocrec*dim_
            dirb(j) = dirb(j) + recb(j)
         end do
         return
      end if
c
      task = ucComm
      dim_ = 3
      if (present(opt)) then
         task = opt
      end if
      if (present(dim)) dim_=dim
c
      ! Reserve buffers space
      rankloc = merge(rank_bis,rank,use_pmecore)
      siz1    = dim_*max(1,nloc)*nrecdir_send1
      siz2    = dim_*max(1,nblocrecdir)
      dimn    = dim_*nloc
      call prmem_request(ucP0%buffMpi_r1,siz1)
      call prmem_request(ucP0%buffMpi_r2,siz2)
      cptr0   = c_loc(ucP0%buffMpi_r1)
      cptr1   = c_loc(ucP0%buffMpi_r2)
      call c_f_pointer(cptr0,buffer ,[siz1])
      call c_f_pointer(cptr1,buffers,[siz2])
 
      if (iand(task,ucRecv).ne.0) then
 21   format(3x,a,' dim(',I0,') opt&task(',2I5,')')
      if (deb_Path) write(*,21) '>> commDDrd_add_c',dim_,task
      !MPI : begin reception in buffer
!$acc wait(def_queue)
!$acc host_data use_device(buffer)
      do i = 1, nrecdir_send1
         precv = precdir_send1(i) 
         if (precv.ne.rank) then
           call MPI_IRECV(buffer(1+(i-1)*dimn),dimn,MPI_TPREC
     &             ,precv,tag,COMM_TINKER,reqr_recdir(i),ierr)
         end if
      end do
!$acc end host_data
      end if
 
      if (iand(task,ucSend).ne.0) then
      !Move in buffers
      do i = 1, nrecdir_recep1
         psend = precdir_recep1(i) 
         if (psend.ne.rank) then
            pbeg = bufbeg(psend+1)-1 
            dlen = domlen(psend+1)
!$acc parallel loop collapse(2) async(def_queue)
!$acc&      present(repartrec,locrec,glob,recb,buffers)
            do j = 1, dlen; do k = 1,dim_
               jloc  = pbeg+j
               jglob = glob(jloc)
               if (repartrec(jglob).eq.rankloc) then
                 jlocrec = locrec(jglob)
                 buffers(k+(jloc-1)*dim_) = recb(k+(jlocrec-1)*dim_)
               else
                 buffers(k+(jloc-1)*dim_) = 0
               end if
            end do; end do
         end if
      end do
!$acc wait(def_queue)
c
      ! MPI : Start Sending
!$acc host_data use_device(buffers)
      do i = 1, nrecdir_recep1
         psend =     precdir_recep1(i)
         if (psend.ne.rank) then
            idxs  = 1 + (bufbeg(psend+1)-1)*dim_
            dlen  = dim_*domlen(psend+1)
            call MPI_ISEND(buffers(idxs),dlen,MPI_TPREC
     &              ,psend,tag,COMM_TINKER,reqs_recdir(i),ierr)
         end if
      end do
!$acc end host_data
      end if
c
      if (iand(task,ucWait).ne.0) then
!$acc wait(def_queue)
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank)
     &      call MPI_WAIT(reqr_recdir(i),status,ierr)
      end do
c
c     MPI : move in output
c
      do i = 1, nrecdir_send1
         if (precdir_send1(i).ne.rank) then
!$acc parallel loop async(def_queue) default(present)
            do j = 1,dimn
               dirb(j) = dirb(j) + (buffer(j+(i-1)*dimn))
            end do
         else
!$acc parallel loop collapse(2) async(def_queue) default(present)
            do j = 1,nlocrec; do k = 1,dim_
               iglob = globrec(j)
               iloc  = loc(iglob)
               if (repart(iglob).eq.rank) then
                  dirb(k+(iloc-1)*dim_) = dirb(k+(iloc-1)*dim_)
     &                                  + recb(k+(   j-1)*dim_)
               end if
            end do; end do
         end if
      end do

      do i = 1, nrecdir_recep1
         if (precdir_recep1(i).ne.rank)
     &      call MPI_WAIT(reqs_recdir(i),status,ierr)
      end do
 31   format(3x,a)
      if (btest(tinkerdebug,tindPath)) call MPI_BARRIER(hostcomm,ierr)
      if (deb_Path) write(*,31) '<< commDDrd_add_c'
      end if

      nullify(buffer)
      nullify(buffers)
      end subroutine

      end module
