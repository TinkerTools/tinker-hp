!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ######################################################
!     ##                                                  ##
!     ##  mpiUtils subModule  -- Tinker-HP                ##
!     ##                                                  ##
!     ######################################################
!
!
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

contains
#include "convert.inc.f90"

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
!
      if (nproc.eq.1) then
!$acc parallel loop async(def_queue) default(present)
         do j = 1, nlocrec*dim_
            dirb(j) = dirb(j) + recb(j)
         end do
         return
      end if
!
      task = ucComm
      dim_ = 3
      if (present(opt)) then
         task = opt
      end if
      if (present(dim)) dim_=dim
!
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
21       format(3x,a,' dim(',I0,') opt&task(',2I5,')')
         if (deb_Path) write(*,21) '>> commDDrd_add_c',dim_,task
         !MPI : begin reception in buffer
!$acc wait(def_queue)
!$acc host_data use_device(buffer)
         do i = 1, nrecdir_send1
            precv = precdir_send1(i)
            if (precv.ne.rank) then
               call MPI_IRECV(buffer(1+(i-1)*dimn),dimn,MPI_TPREC&
                  &,precv,tag,COMM_TINKER,reqr_recdir(i),ierr)
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
!$acc parallel loop collapse(2) async(def_queue) &
!$acc      present(repartrec,locrec,glob,recb,buffers)
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
!
         ! MPI : Start Sending
!$acc host_data use_device(buffers)
         do i = 1, nrecdir_recep1
            psend =     precdir_recep1(i)
            if (psend.ne.rank) then
               idxs  = 1 + (bufbeg(psend+1)-1)*dim_
               dlen  = dim_*domlen(psend+1)
               call MPI_ISEND(buffers(idxs),dlen,MPI_TPREC&
                  &,psend,tag,COMM_TINKER,reqs_recdir(i),ierr)
            end if
         end do
!$acc end host_data
      end if
!
      if (iand(task,ucWait).ne.0) then
!$acc wait(def_queue)
         do i = 1, nrecdir_send1
            if (precdir_send1(i).ne.rank)&
               &call MPI_WAIT(reqr_recdir(i),status,ierr)
         end do
!
!     MPI : move in output
!
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
                        dirb(k+(iloc-1)*dim_) = dirb(k+(iloc-1)*dim_)&
                           &+ recb(k+(   j-1)*dim_)
                     end if
                  end do; end do
            end if
         end do

         do i = 1, nrecdir_recep1
            if (precdir_recep1(i).ne.rank)&
               &call MPI_WAIT(reqs_recdir(i),status,ierr)
         end do
31       format(3x,a)
         if (btest(tinkerdebug,tindPath)) call MPI_BARRIER(hostcomm,ierr)
         if (deb_Path) write(*,31) '<< commDDrd_add_c'
      end if

      nullify(buffer)
      nullify(buffers)
   end subroutine

end submodule
