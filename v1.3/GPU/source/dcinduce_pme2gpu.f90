!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     Implemented by Dominique Nocito and Gregory Beran
!     solve the polarization equations with a DC/JI-DIIS algorithm
!
!     literature reference:
!     "Fast Divide-and-Conquer Scheme for Evaluating Polarization
!     in Classical Force Fields" D. Nocito and G. Beran. J. Chem. Phys. 146,
!     114103 (2017)
!
#include "tinker_macro.h"
subroutine dcinduce_pme2gpu
   use atmlst
   use divcon
   use domdec
   use ewald
   use iounit
   use inform    ,only: deb_Path
   use interfaces,only: otf_dc_efld0_directgpu_p
   use math
   use mpole
   use parallel, only: polar_comm_t, polar_post_comm, polar_finish, &
      &polar_post_dip_comm, polar_finish_dip_comm, domain_workspace_tp2d, &
      &domain_workspace_tp3d
   use parallel, only: commrecdirfields, commrecdirfieldsgpu, commrecdirsolv, commrecdirsolvgpu
   use parallel, only: commfield, commfieldgpu, commfieldshort
   use pme
   use polar
   use polpot
   use potent
   use tinMemory ,only: prmem_request
   use units
   use uprior
   use utils
   use mdevice, only: after, direct_work, reciprocal_work, record_platform_error
   use status, only: failed, status_t
   use mpi
   implicit none
!
!     without separate cores for reciprocal part
!
   integer i, j, k, nrhs, proc
!
!     MPI
!
   integer iglob, ierr, iipole
   integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
   type(polar_comm_t) :: comm_ticket
   real(t_p) time0,time1
   type(status_t) :: dep_status
!
   parameter (nrhs=2)
   real(t_p)  wtime0, wtime1, wtime2, udsum, upsum
   real(t_p)  term
   real(t_p), allocatable :: ef(:,:,:), mu(:,:,:),&
      &cphi(:,:), murec(:,:,:)
   real(t_p), pointer :: buffermpi1(:,:),buffermpi2(:,:)
   real(t_p), pointer :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
!
!     Get direct space induced field for divide and conquer methods
!
   external pc_dc_tmatxb_pmegpu
   external otf_dc_tmatxb_pmegpu
!
   if (deb_Path)&
      &write(*,'(2x,a)') 'dcinduce_pme2gpu'
   if (.not.use_polar) return
#ifdef _OPENACC
   if (clsttype .eq. 1) then
      call record_platform_error(-9, 'kmean clustering is not supported with OpenACC')
      return
   end if
   if (precomp .eq. 1) then
      call record_platform_error(-9, 'dcinduce pc_dc is not supported with OpenACC')
      return
   end if
#endif
!
1000 format(' illegal polalg in newinduce.')
1010 format(' time for the ',a,F14.5)
1020 format(' total elapsed time in newinduce: ',F14.5)
!
!     allocate some memory and clear the arrays:
!
   call domain_workspace_tp3d(1,3,nrhs,max(npoleloc,1),buffermpimu1)
   call domain_workspace_tp3d(2,3,nrhs,max(npolerecloc,1),buffermpimu2)
   call domain_workspace_tp2d(1,10,max(npoleloc,1),buffermpi1)
   call domain_workspace_tp2d(2,10,max(npolerecloc,1),buffermpi2)
!
   allocate (mu   (3,nrhs,max(1,npolebloc)))
   allocate (murec(3,nrhs,max(1,npolerecloc)))
   allocate (ef   (3,nrhs,max(1,npolebloc)))
   allocate (cphi (10,max(npoleloc,1)))

   if (.not.use_mpole) then
      j = max(npolerecloc,1)
      call prmem_request(cphirec,10,j,async=.true.)
      call prmem_request(fphirec,20,j,async=.true.)
   end if
!
!$acc data create(mu,murec,ef,cphi) &
!$acc     present(cphirec,fphirec,rpole,uind,uinp,polarity, &
!$acc  upalt,udalt,bpred,poleloc,poleglob,polerecglob,ipole, &
!$acc  repart) async

   allocate (reqrecdirsend(nproc))
   allocate (reqrecdirrec(nproc))

   ! set arrays to zero
   call set_to_zero2(mu,ef       ,3*nrhs*npolebloc,reciprocal_work%queue)
   call set_to_zero1(buffermpimu1,3*nrhs*npoleloc ,reciprocal_work%queue)
   call set_to_zero2(murec,buffermpimu2,3*nrhs*npolerecloc,&
      &reciprocal_work%queue)
   call set_to_zero2(cphi   ,buffermpi1,10*npoleloc   ,reciprocal_work%queue)
   call set_to_zero2(cphirec,buffermpi2,10*npolerecloc,reciprocal_work%queue)
   call set_to_zero1(fphirec,           20*npolerecloc,reciprocal_work%queue)
!
!     compute the electric fields:
!
   wtime0 = mpi_wtime()
!
   if(clsttype .eq. 1) then
      !kmeans clustering
      call cluster1
   else if(clsttype .eq. 2) then
      !segmentation clustering
      call cluster2
   end if

#ifdef _OPENACC
   if (direct_work%queue.ne.reciprocal_work%queue)& !Start direct_work%stream cover
      &call after(direct_work,reciprocal_work,dep_status)
      if (failed(dep_status)) __TINKER_FATAL__
#endif

   if(precomp .eq. 1) then
      call pc_dc_efld0_direct(nrhs,ef)
   else
      call otf_dc_efld0_directgpu_p(nrhs,ef)
   end if
!
!     Factorization of Z has been interleaved with passing fields
!
   call commfield2gpu(nrhs,ef)
!
!    compute the reciprocal space contribution (fields)
!
   call efld0_recipgpu(cphi,ef)
!
!     Begin the reception of the reciprocal fields
!
   call commrecdirfieldsgpu(0,cphirec,cphi,buffermpi1,buffermpi2,&
      &reqrecdirrec,reqrecdirsend)
   call commrecdirfieldsgpu(1,cphirec,cphi,buffermpi1,buffermpi2,&
      &reqrecdirrec,reqrecdirsend)
   call commrecdirfieldsgpu(2,cphirec,cphi,buffermpi1,buffermpi2,&
      &reqrecdirrec,reqrecdirsend)
!
   ! Direct exchange is posted after the initial field is ready.
!
!     Add direct and reciprocal fields
!
   term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
#ifdef _OPENACC
   if (direct_work%queue.ne.reciprocal_work%queue) then !close 'direct_work%stream' recover
      call after(reciprocal_work,direct_work,dep_status)
      if (failed(dep_status)) __TINKER_FATAL__
   end if
#endif
!$acc parallel loop collapse(3) async(reciprocal_work%queue)
   do i = 1,npoleloc
      do j = 1,nrhs
         do k = 1,3
            iipole = poleglob(i)
            ef(k,j,i) = ef(k,j,i) + term*rpole(k+1,iipole)&
               &- cphi(k+1,i)
         end do
      end do
   end do
   wtime1 = mpi_wtime()

!
!     guess the dipoles:
!
!     predicted values for always stable predictor-corrector method
!
   if (use_pred .and. nualt.eq.maxualt) then
      call pred_InitField(mu)
   else if (polgsf.eq.0) then
!$acc parallel loop collapse(3) async(reciprocal_work%queue)
      do i = 1, npoleloc
         do k = 1, nrhs
            do j = 1, 3
               iipole = poleglob(i)
               mu(j,k,i) = polarity(iipole)*ef(j,k,i)
            end do
         end do
      end do
   else
!$acc parallel loop collapse(3) async(reciprocal_work%queue)
      do i = 1, npoleloc
         do j = 1,nrhs
            do k = 1, 3
               iipole = poleglob(i)
               if (btest(j,0)) then !j==1
                  mu(k,j,i) = uind(k,iipole)
               else  !j/=1
                  mu(k,j,i) = uinp(k,iipole)
               end if
            end do
         end do
      end do
   end if

   call polar_post_comm(comm_ticket,nrhs,mu)
   call polar_finish(comm_ticket,nrhs,mu)
   call polar_post_dip_comm(comm_ticket,nrhs,murec,mu,buffermpimu1,buffermpimu2)
   call polar_finish_dip_comm(comm_ticket,nrhs,murec,mu,buffermpimu1,buffermpimu2)
!
!     now, call the proper solver.
!
   if(precomp .eq. 1) then
      call inducedc_pme2gpu( pc_dc_tmatxb_pmegpu,nrhs,.true.,ef,mu,&
         &murec)
   else
      call inducedc_pme2gpu(otf_dc_tmatxb_pmegpu,nrhs,.true.,ef,mu,&
         &murec)
   end if


   wtime2 = mpi_wtime()
   if (polprt.ge.1.and.rank.eq.0) then
      if (polprt.ge.2) then
         write(iout,1010) 'fields:  ', wtime1-wtime0
         write(iout,1010) 'dipoles: ', wtime2-wtime1
      end if
      write(iout,1020) wtime2 - wtime0
   end if
!
!     move the computed dipoles in the common block.
!
!$acc parallel loop collapse(2) async(reciprocal_work%queue)
   do i = 1, npolebloc
      do j = 1, 3
         iipole = poleglob(i)
         uind(j,iipole) = mu(j,1,i)
         uinp(j,iipole) = mu(j,2,i)
      end do
   end do
!$acc parallel loop collapse(2) async(reciprocal_work%queue)
   do i = 1, npolerecloc
      do j = 1, 3
         iipole = polerecglob(i)
         iglob  = ipole(iipole)
         if (repart(iglob).ne.rank) then
            uind(j,iipole) = murec(j,1,i)
            uinp(j,iipole) = murec(j,2,i)
         else
            uind(j,iipole) = mu(j,1,poleloc(iipole))
            uinp(j,iipole) = mu(j,2,poleloc(iipole))
         end if
      end do
   end do
!$acc end data
!
!     update the lists of previous induced dipole values
!
   if (use_pred) call pred_SetAlt

   deallocate (reqrecdirrec)
   deallocate (reqrecdirsend)
   deallocate (ef)
   deallocate (mu)
   deallocate (murec)
   deallocate (cphi)
end

subroutine inducedc_pme2gpu(matvec,nrhs,dodiis,ef,mu,murec)
   use atmlst
   use divcon
   use domdec
   use ewald
   use iounit
   use inform     ,only: deb_Path
   use interfaces
   use math
   use mpole
   use parallel, only: polar_comm_t, polar_post_comm, polar_finish, &
      &polar_post_dip_comm, polar_finish_dip_comm, domain_workspace_tp3d
   use parallel, only: commrecdirfields, commrecdirfieldsgpu, commrecdirsolv, commrecdirsolvgpu
   use parallel, only: commfield, commfieldgpu, commfieldshort
   use polar
   use polpot
   use timestat
   use tinMemory  ,only: prmem_request
   use mpi
   use sizes
   use utils
   use mdevice, only: after, direct_work, reciprocal_work
   use status, only: failed, status_t
   implicit none
!
!     solves the polarizable dipoles linear equations by jacobi
!     iterations coupled to the DIIS extrapolation.
!
   integer  ,intent(in   ):: nrhs
   logical  ,intent(in   ):: dodiis
   real(t_p),intent(in   ):: ef   (3,nrhs,npolebloc)
   real(t_p),intent(inout):: mu   (3,nrhs,npolebloc)&
      &, murec(3,nrhs,npolerecloc)
   external matvec

   integer info, proc, knd, pos, mdim, ierror
   integer i, ii, j, k, k1, it, ind, ndismx, nmat, lenb
   integer iglob, iipole, ierr
   integer kst, zst, read0, write0
   real(t_p) zero, one, rnorm(2), rr
   real(t_p) rnorm1,rnorm2
   real(t_p) term
   real(t_p) time0,time1,time2,time3
   parameter (ndismx=25,zero=0,one=1)
   type(status_t) :: dep_status

   integer ipiv(ndismx+1)
   real(t_p), allocatable :: munew(:,:,:),h(:,:,:), cex(:)&
      &, xdiis(:,:), ediis(:,:), bmat(:,:), bloc(:,:)
   real(t_p), allocatable :: dcfield1(:),dcfield2(:)
   real(t_p), allocatable :: dipfield(:,:,:)&
      &, dipfieldbis(:,:,:)
   real(t_p), pointer :: buffermpi1(:,:,:),buffermpi2(:,:,:)
   real(t_p), pointer :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
!
!     MPI
!
   integer reqnorm,reqdiis(2*ndismx+1)
   integer status(MPI_STATUS_SIZE)
   integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
   type(polar_comm_t) :: comm_ticket, dip_ticket
!
!
1000 format(' itsolv converged after ',I3,' iterations.',/,&
      &' final residual norm = ',3D14.7)
1010 format(' residual norm at iteration ',I3,':',3D12.2)
1020 format(' Jacobi/DIIS solver: induced dipoles',/,&
      &' ipole       mux         muy         muz')
1021 format(' Jacobi/DIIS solver: induced p-dipoles',/,&
      &' ipole       mux         muy         muz')
1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
!
   if (deb_Path)&
      &write(*,'(3x,a)') 'inducedc_pme2gpu'

   call domain_workspace_tp3d(1,3,nrhs,max(npoleloc,1),buffermpi1)
   call domain_workspace_tp3d(2,3,nrhs,max(npoleloc,1),buffermpimu1)
   call domain_workspace_tp3d(3,3,nrhs,max(npolerecloc,1),buffermpi2)
   call domain_workspace_tp3d(4,3,nrhs,max(npolerecloc,1),buffermpimu2)
!
   allocate (munew      (3,nrhs,max(1,npolebloc  )))
   allocate (h          (3,nrhs,max(1,npolebloc  )))
   allocate (dipfield   (3,nrhs,max(1,npoleloc   )))
   allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
   allocate (dcfield1   (maxdim))
   allocate (dcfield2   (maxdim))
   allocate (reqrecdirsend(nproc))
   allocate (reqrecdirrec (nproc))
!$acc enter data create(munew,h,dipfield,dipfieldbis,dcfield1,dcfield2, &
!$acc    rnorm1,rnorm2,ipiv)  async
!
#ifdef _OPENACC
   k1 = maxdim**2
   call prmem_request(zmatDn,k1,async=.true.)
#endif
   if (dodiis) then
      nmat = 1
      lenb = ndismx + 1
      allocate (xdiis(3*nrhs*max(npoleloc,1),ndismx))
      allocate (ediis(3*nrhs*max(npoleloc,1),ndismx))
      allocate (bmat(lenb,lenb))
      allocate (bloc(lenb,lenb))
      allocate (cex(lenb))
!$acc enter data create(xdiis,ediis,bmat,bloc,cex) async
      call set_to_zero1(bmat,lenb*lenb,reciprocal_work%queue)
   end if

   !Initiate arrays
   call set_to_zero3(buffermpi1,buffermpimu1,dipfield&
      &,size(dipfield),reciprocal_work%queue)
   call set_to_zero3(buffermpi2,buffermpimu2,dipfieldbis&
      &,size(dipfieldbis),reciprocal_work%queue)
   call set_to_zero2(munew,h,size(h),reciprocal_work%queue)
   call set_to_zero2(dcfield2,dcfield2,size(dcfield2),reciprocal_work%queue)
#ifdef _OPENACC
   !start Async overlapping with matvec
   if (direct_work%queue.ne.reciprocal_work%queue)&
      &call after(direct_work,reciprocal_work,dep_status)
      if (failed(dep_status)) __TINKER_FATAL__
#endif
!$acc data present(klst,mu,poleglob,ipole,atmofst,dipfield,dipfieldbis, &
!$acc   dcfield2,dcfield1,munew,ediis,xdiis) async
!
!     main loop:
!
   do it = 1, politer
      rnorm(1) = zero
      rnorm(2) = zero
!$acc serial async present(rnorm1,rnorm2)
      rnorm1   = zero; rnorm2 = zero;
!$acc end serial
!
      time0 = mpi_wtime()
      call matvec(nrhs,.false.,mu,h)
      time1 = mpi_wtime()
      if (it.eq.1) timerealdip = timerealdip + time1-time0
      call tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
      time2 = mpi_wtime()
      if (it.eq.1) timerecdip = timerecdip + time2-time1

      call commfieldgpu(nrhs,h)

      call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,&
         &buffermpi2,reqrecdirrec,reqrecdirsend)
!
      call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,&
         &buffermpi2,reqrecdirrec,reqrecdirsend)
!
!
!     jacobi step:
!
      call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,&
         &buffermpi2,reqrecdirrec,reqrecdirsend)
      time2 = mpi_wtime()
      term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi

      !print*,'mean mdim',3*sum(npergrp)/km
      do k = 1,km
         mdim = npergrp(k)*3

         if (mdim .gt. 0) then

            !ind = 0
            !Build the field for each block
!$acc parallel loop collapse(3) async
            do i = 1, npergrp(k)
               do k1 = 1,2
                  do j = 1,3
                     !ind = ind + 1
                     ii     = klst(i,k)
                     iipole = poleglob(ii)
                     iglob  = ipole(iipole)
                     pos    = (atmofst(iglob)-1)*3
                     if (btest(k1,0)) then ! if (k==2)
                        dcfield2(pos+j) = ef(j,2,ii) -&
                           &(h(j,2,ii) + dipfield(j,2,ii)-term*mu(j,2,ii))
                     else                  ! if (k==1)
                        dcfield1(pos+j) = ef(j,1,ii) -&
                           &(h(j,1,ii) + dipfield(j,1,ii)-term*mu(j,1,ii))
                     end if
                  end do
               end do
            end do

            knd = kofst(k) + mdim*(mdim+1)/2
#ifdef _OPENACC
            kst = kofst(k)+1
!$acc parallel loop vector_length(32) async
            do i = 0,mdim-1
               zst  = i*mdim - (i*(i-1))/2
!$acc loop vector
               do j = i,mdim-1
                  read0  = zst    + j-i
                  write0 = i*mdim + j+1
                  zmatDn(write0) = zmat(kst+read0)
               end do
            end do
!$acc host_data use_device(zmatDn,dcfield1,dcfield2)
            call cuPOTRS(mdim,zmatDn,mdim,dcfield1,mdim,reciprocal_work%stream)
            call cuPOTRS(mdim,zmatDn,mdim,dcfield2,mdim,reciprocal_work%stream)
!$acc end host_data
#else
            call M_PPTRS("L",mdim,1,zmat(kofst(k)+1:knd),&
               &dcfield1(1:mdim),mdim,ierror)
            if(ierror .ne. 0) print*,"error with solving dipoles"
            call M_PPTRS("L",mdim,1,zmat(kofst(k)+1:knd),&
               &dcfield2(1:mdim),mdim,ierror)
            if(ierror .ne. 0) print*,"error with solving dipoles"
#endif
!
!  Update dipoles
!
!$acc parallel loop collapse(3) async
            do i = 1, npergrp(k); do k1=1,2; do j=1,3
               ii    = klst(i,k)
               iglob = ipole(poleglob(ii))
               pos   = (atmofst(iglob)-1)*3
               if (btest(k1,0)) then
                  munew(j,2,ii) = dcfield2(pos+j)
               else
                  munew(j,1,ii) = dcfield1(pos+j)
               end if
            end do; end do; end do

         end if
      end do

!$acc parallel loop collapse(3) present(rnorm1,rnorm2) async
      do i = 1,npoleloc; do k=1,nrhs; do j=1,3
         if (btest(k,0)) then
            rnorm2 = rnorm2 + (munew(j,k,i)-mu(j,k,i))**2
         else
            rnorm1 = rnorm1 + (munew(j,k,i)-mu(j,k,i))**2
         end if
      end do; end do; end do
!$acc update host(rnorm1,rnorm2) async

      if (nproc>1) then
!$acc wait
         rnorm(1)=rnorm1; rnorm(2)=rnorm2;
         call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_TPREC,&
            &MPI_SUM,COMM_TINKER,reqnorm,ierr)
      end if

      if (dodiis) then

!$acc parallel loop collapse(3) async
         do i = 1,npoleloc; do k = 1,nrhs; do j = 1,3
            ind = (i-1)*6 + (k-1)*3 + j
            xdiis(ind,nmat) = munew(j,k,i)
            ediis(ind,nmat) = munew(j,k,i) - mu(j,k,i)
         end do; end do; end do
!
!          Compute Pulay's Matrix and extrapolate
!
         if(nocomdiis .eq. 1)then
            call no_com_diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,&
               &bmat,nmat)
         else
            reqdiis(:) = MPI_REQUEST_NULL
            call diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,&
               &bmat,nmat,reqdiis,COMM_TINKER)
            !TODO This loop is no longer useful; MPI_IAllReduce is not CUDA-Aware
            do i = 1, 2*nmat-3
               call MPI_WAIT(reqdiis(i),status,ierr)
            end do
         end if

         call utils_amove(lenb*lenb,bmat,bloc,reciprocal_work%queue)  !copy bmat in bloc
         call set_to_zero1(cex,lenb,reciprocal_work%queue)
         call set_to_zero1(munew,size(munew),reciprocal_work%queue)
!$acc serial async
         cex(1) = one
!$acc end serial

#ifdef _OPENACC
!$acc host_data use_device(bloc,ipiv,cex)
         call cuGESV(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,reciprocal_work%stream)
!$acc end host_data
#else
         call M_gesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#endif
         call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)

      end if
      time2 = mpi_wtime()
!
      call polar_post_comm(comm_ticket,nrhs,mu,.false.,munew)
      call polar_finish(comm_ticket,nrhs,mu)
      call polar_post_dip_comm(dip_ticket,nrhs,murec,munew,buffermpimu1,buffermpimu2)
      call polar_finish_dip_comm(dip_ticket,nrhs,murec,munew,buffermpimu1,buffermpimu2)
      call utils_amove(3*nrhs*npoleloc,munew,mu,reciprocal_work%queue)   !copy munew in mu
!
!     compute the norm of the increment.
!
      if (nproc>1) then
         call MPI_WAIT(reqnorm,status,ierr)
      else
!$acc wait
         rnorm(1)=rnorm1; rnorm(2)=rnorm2
      end if

      rr = zero
      do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/real(3*npolar,t_p))
         rr = max(rnorm(k),rr)
      end do
      if (rr.lt.poleps .and. it .ge. 3) then
         if (polprt.ge.1.and.rank.eq.0)&
            &write(6,1000) it, (rnorm(k), k = 1, nrhs)
         goto 10
      end if
      if (polprt.ge.1.and.rank.eq.0)&
         &write(6,1010) it, (rnorm(k), k = 1, nrhs)
   end do

10 continue
!$acc end data

   if (polprt.ge.3) then
      write(iout,1020)
!$acc update host(poleglob)
      do i = 1, npoleloc
         iipole = poleglob(i)
         iglob = ipole(iipole)
         write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
      end do
   end if
   if (polprt.ge.4) then
      write(iout,1021)
!$acc update host(poleglob)
      do i = 1, npoleloc
         iipole = poleglob(i)
         iglob = poleglob(iipole)
         write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
      end do
   end if
!
!     free the memory.
!
!$acc exit data delete(munew,dipfield,dipfieldbis,dcfield1,dcfield2,h) &
!$acc     async
   deallocate (reqrecdirsend)
   deallocate (reqrecdirrec)
   deallocate (munew)
   deallocate (dipfield)
   deallocate (dipfieldbis)
   deallocate (dcfield1)
   deallocate (dcfield2)
   deallocate (h)
   if (dodiis) then
!$acc exit data delete(xdiis,ediis,bmat,bloc,cex) async
      deallocate (xdiis)
      deallocate (ediis)
      deallocate (bmat)
      deallocate (bloc)
      deallocate (cex)
   end if

end
