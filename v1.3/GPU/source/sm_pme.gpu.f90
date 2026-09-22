!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  submodule pme_lifecycle  --  PME owner lifecycle           ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
submodule (pme) pme_lifecycle
   use mpi, only: MPI_SUCCESS, MPI_Type_free, MPI_Waitall, MPI_STATUSES_IGNORE
   use pme1
   use tinMemory, only: prmem_request,mfhostacc,s_prmem,sd_prmem,szoTp
#ifdef _CUDA
   use pmestuffcu, only: igrid_t,thetai1_t,thetai2_t,thetai3_t,&
      &qgridin_t,qgrid2in_t
#endif
   implicit none

contains

   module subroutine initialize_pme_device(queue)
      implicit none
      integer,intent(in) :: queue

      if (pme_device_ready) return
!$acc enter data create(pme_ftc,pme_ctf) async(queue)
      pme_device_ready=.true.
   end subroutine initialize_pme_device

   module subroutine release_pme_device(queue, cleanup_status)
      implicit none
      integer,intent(in) :: queue
      integer,intent(out),optional :: cleanup_status
      integer :: ierr, first_error

      first_error = MPI_SUCCESS

      pme_ftc_ready=.false.
      pme_ctf_ready=.false.
      if (pme_device_ready) then
!$acc wait(queue)
!$acc exit data delete(pme_ftc,pme_ctf)
      end if

      if (allocated(front_send_req)) then
         call MPI_Waitall(size(front_send_req),front_send_req,MPI_STATUSES_IGNORE,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
      end if
      if (allocated(front_recv_req)) then
         call MPI_Waitall(size(front_recv_req),front_recv_req,MPI_STATUSES_IGNORE,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
      end if
      if (allocated(front_split_snd)) then
         call MPI_Waitall(size(front_split_snd),front_split_snd,MPI_STATUSES_IGNORE,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
      end if
      if (allocated(front_split_rcv)) then
         call MPI_Waitall(size(front_split_rcv),front_split_rcv,MPI_STATUSES_IGNORE,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
      end if
      if (allocated(back_send_req)) then
         call MPI_Waitall(size(back_send_req),back_send_req,MPI_STATUSES_IGNORE,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
      end if
      if (allocated(back_recv_req)) then
         call MPI_Waitall(size(back_recv_req),back_recv_req,MPI_STATUSES_IGNORE,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
      end if
      if (allocated(back_split_snd)) then
         call MPI_Waitall(size(back_split_snd),back_split_snd,MPI_STATUSES_IGNORE,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
      end if
      if (allocated(back_split_rcv)) then
         call MPI_Waitall(size(back_split_rcv),back_split_rcv,MPI_STATUSES_IGNORE,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
      end if

#ifdef _CUDA
      if (associated(igrid_t)) nullify(igrid_t)
      if (associated(thetai1_t)) nullify(thetai1_t)
      if (associated(thetai2_t)) nullify(thetai2_t)
      if (associated(thetai3_t)) nullify(thetai3_t)
      if (associated(qgridin_t)) nullify(qgridin_t)
      if (associated(qgrid2in_t)) nullify(qgrid2in_t)
#endif
      if (allocated(qgridin_2d).or.allocated(qgridout_2d).or.&
         &allocated(qgrid2in_2d).or.allocated(qgrid2out_2d)) then
!$acc exit data detach(qgridin_p,qgrid2in_p,qgridout_p,qgrid2out_p)
      end if
      if (associated(thetai1_p)) nullify(thetai1_p)
      if (associated(thetai2_p)) nullify(thetai2_p)
      if (associated(thetai3_p)) nullify(thetai3_p)
      if (associated(qgridin_p)) nullify(qgridin_p)
      if (associated(qgridout_p)) nullify(qgridout_p)
      if (associated(qgrid2in_p)) nullify(qgrid2in_p)
      if (associated(qgrid2out_p)) nullify(qgrid2out_p)

      if (allocated(thetai1)) call prmem_request(thetai1,0,0,0,config=mfhostacc)
      if (allocated(thetai2)) call prmem_request(thetai2,0,0,0,config=mfhostacc)
      if (allocated(thetai3)) call prmem_request(thetai3,0,0,0,config=mfhostacc)
      if (allocated(qfac_2d)) call prmem_request(qfac_2d,0,0,0,config=mfhostacc)
      if (allocated(cphirec)) call prmem_request(cphirec,0,0,config=mfhostacc)
      if (allocated(fphirec)) call prmem_request(fphirec,0,0,config=mfhostacc)
      if (allocated(igrid)) call prmem_request(igrid,0,0,config=mfhostacc)

      if (allocated(qgridin_2d)) then
!$acc exit data delete(qgridin_2d)
         s_prmem=s_prmem-size(qgridin_2d)*szoTp
#ifdef _OPENACC
         sd_prmem=sd_prmem-size(qgridin_2d)*szoTp
#endif
         deallocate(qgridin_2d)
      end if
      if (allocated(qgridout_2d)) then
!$acc exit data delete(qgridout_2d)
         s_prmem=s_prmem-size(qgridout_2d)*szoTp
#ifdef _OPENACC
         sd_prmem=sd_prmem-size(qgridout_2d)*szoTp
#endif
         deallocate(qgridout_2d)
      end if
      if (allocated(qgrid2in_2d)) then
!$acc exit data delete(qgrid2in_2d)
         s_prmem=s_prmem-size(qgrid2in_2d)*szoTp
#ifdef _OPENACC
         sd_prmem=sd_prmem-size(qgrid2in_2d)*szoTp
#endif
         deallocate(qgrid2in_2d)
      end if
      if (allocated(qgrid2out_2d)) then
!$acc exit data delete(qgrid2out_2d)
         s_prmem=s_prmem-size(qgrid2out_2d)*szoTp
#ifdef _OPENACC
         sd_prmem=sd_prmem-size(qgrid2out_2d)*szoTp
#endif
         deallocate(qgrid2out_2d)
      end if
      if (allocated(qgridmpi)) then
!$acc exit data delete(qgridmpi)
         s_prmem=s_prmem-size(qgridmpi)*szoTp
#ifdef _OPENACC
         sd_prmem=sd_prmem-size(qgridmpi)*szoTp
#endif
         deallocate(qgridmpi)
      end if
      if (Is_MPI_GRID1D) then
         call MPI_Type_free(MPI_GRID1D,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
         Is_MPI_GRID1D=.false.
      end if
      if (Is_MPI_GRID2D) then
         call MPI_Type_free(MPI_GRID2D,ierr)
         if (ierr /= MPI_SUCCESS .and. first_error == MPI_SUCCESS) first_error = ierr
         Is_MPI_GRID2D=.false.
      end if
      if (allocated(front_send_req)) deallocate(front_send_req)
      if (allocated(front_recv_req)) deallocate(front_recv_req)
      if (allocated(front_split_snd)) deallocate(front_split_snd)
      if (allocated(front_split_rcv)) deallocate(front_split_rcv)
      if (allocated(back_send_req)) deallocate(back_send_req)
      if (allocated(back_recv_req)) deallocate(back_recv_req)
      if (allocated(back_split_snd)) deallocate(back_split_snd)
      if (allocated(back_split_rcv)) deallocate(back_split_rcv)
      pme_device_ready=.false.
      if (present(cleanup_status)) cleanup_status = first_error
   end subroutine release_pme_device

end submodule pme_lifecycle
