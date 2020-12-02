c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c      
c     Implemented by Dominique Nocito and Gregory Beran
c     solve the polarization equations with a DC/JI-DIIS algorithm
c
c     literature reference:
c     "Fast Divide-and-Conquer Scheme for Evaluating Polarization 
c     in Classical Force Fields" D. Nocito and G. Beran. J. Chem. Phys. 146,
c     114103 (2017)
c
#include "tinker_precision.h"
      subroutine dcinduce_pme2gpu
      use atmlst
      use divcon
      use domdec
      use ewald
      use iounit
      use interfaces,only: otf_dc_efld0_directgpu_p
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use units
      use uprior
      use utilcomm ,buffermpi1=>buffermpi2d,buffermpi2=>buffermpi2d1
      use utils
      use utilgpu
      use mpi
      implicit none
c
c     without separate cores for reciprocal part
c      
      integer i, j, k, nrhs, proc
c
c     MPI
c
      integer iglob, ierr, iipole
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      real(t_p) time0,time1
c
      parameter (nrhs=2)
      real(t_p)  wtime0, wtime1, wtime2, udsum, upsum
      real(t_p)  term, xx(1)
      real(t_p), allocatable :: ef(:,:,:), mu(:,:,:),
     &                          cphi(:,:), murec(:,:,:)
c
c     Get direct space induced field for divide and conquer methods
c
      external pc_dc_tmatxb_pmegpu
      external otf_dc_tmatxb_pmegpu
c
      if (rank.eq.0.and.tinkerdebug)
     &   write(*,'(2x,a)') 'dcinduce_pme2gpu'
      if (.not.use_polar) return
c
 1000 format(' illegal polalg in newinduce.')
 1010 format(' time for the ',a,F14.5)
 1020 format(' total elapsed time in newinduce: ',F14.5)
c
c     allocate some memory and clear the arrays:
c
      call prmem_request(buffermpimu1,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(buffermpimu2,3,nrhs,max(npolerecloc,1),
     &     async=.true.)
      call prmem_request(buffermpi1,10,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(buffermpi2,10,max(npolerecloc,1),
     &     async=.true.)
c
      allocate (mu   (3,nrhs,max(1,npolebloc)))
      allocate (murec(3,nrhs,max(1,npolerecloc)))
      allocate (ef   (3,nrhs,max(1,npolebloc)))
      allocate (cphi (10,max(npoleloc,1)))

      if (.not.use_mpole) then
         j = max(npolerecloc,1)
         call prmem_request(cphirec,10,j,async=.true.)
         call prmem_request(fphirec,20,j,async=.true.)
      end if
c
!$acc data create(mu,murec,ef,cphi)
!$acc&     present(cphirec,fphirec,rpole,uind,uinp,polarity,
!$acc&  upalt,udalt,bpred,poleloc,poleglob,polerecglob,ipole,
!$acc&  repart) async

      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))

      ! set arrays to zero
      call set_to_zero2(mu,ef       ,3*nrhs*npolebloc,rec_queue)
      call set_to_zero1(buffermpimu1,3*nrhs*npoleloc ,rec_queue)
      call set_to_zero2(murec,buffermpimu2,3*nrhs*npolerecloc,
     &                 rec_queue)
      call set_to_zero2(cphi   ,buffermpi1,10*npoleloc   ,rec_queue)
      call set_to_zero2(cphirec,buffermpi2,10*npolerecloc,rec_queue)
      call set_to_zero1(fphirec,           20*npolerecloc,rec_queue)
c
c     compute the electric fields:
c
      wtime0 = mpi_wtime()
c
      if(clsttype .eq. 1) then
        !kmeans clustering
#ifdef _OPENACC
        call openacc_abort("kmean clustering not yet supported with"
     &       //" OpenAcc !!")
#endif
        call cluster1
      else if(clsttype .eq. 2) then
        !segmentation clustering
        call cluster2
      end if

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) !Start dir_stream cover
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)
#endif

      if(precomp .eq. 1) then
#ifdef _OPENACC
         call openacc_abort("dcinduce pc_dc not supported with OpenAcc")
#endif
         call pc_dc_efld0_directgpu(nrhs,ef)
      else
         call otf_dc_efld0_directgpu_p(nrhs,ef)
      end if
c
c     Factorization of Z has been interleaved with passing fields
c
      call commfield2gpu(nrhs,ef)
c
c    compute the reciprocal space contribution (fields)
c
      call efld0_recipgpu(cphi,ef)
c
c     Begin the reception of the reciprocal fields
c
      call commrecdirfieldsgpu(0,cphirec,cphi,buffermpi1,buffermpi2,
     $     reqrecdirrec,reqrecdirsend)
      call commrecdirfieldsgpu(1,cphirec,cphi,buffermpi1,buffermpi2,
     $     reqrecdirrec,reqrecdirsend)
      call commrecdirfieldsgpu(2,cphirec,cphi,buffermpi1,buffermpi2,
     $     reqrecdirrec,reqrecdirsend)
c
      call commdirdirgpu(nrhs,0,mu,reqrec,reqsend)
c
c     Add direct and reciprocal fields
c
      term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
      def_queue = rec_queue
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then !close 'dir_stream' recover
         call stream_wait_async(dir_stream,rec_stream,dir_event)
      end if
#endif
!$acc parallel loop collapse(3) async(def_queue)
      do i = 1,npoleloc
         do j = 1,nrhs
            do k = 1,3
               iipole = poleglob(i)
               ef(k,j,i) = ef(k,j,i) + term*rpole(k+1,iipole)
     &                               - cphi(k+1,i)
            end do
         end do
      end do
      wtime1 = mpi_wtime()

c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
      if (use_pred .and. nualt.eq.maxualt) then
         call pred_InitField(mu)
      else if (polgsf.eq.0) then
!$acc parallel loop collapse(3) async(def_queue)
         do i = 1, npoleloc
            do k = 1, nrhs
               do j = 1, 3
                  iipole = poleglob(i)
                  mu(j,k,i) = polarity(iipole)*ef(j,k,i)
               end do
            end do
         end do
      else
!$acc parallel loop collapse(3) async(def_queue)
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

      call commdirdirgpu(nrhs,1,mu,reqrec,reqsend)
      call commdirdirgpu(nrhs,2,mu,reqrec,reqsend)
      call commrecdirdipgpu(nrhs,0,murec,mu,buffermpimu1,buffermpimu2,
     $     req2rec,req2send)
      call commrecdirdipgpu(nrhs,1,murec,mu,buffermpimu1,buffermpimu2,
     $     req2rec,req2send)
      call commrecdirdipgpu(nrhs,2,murec,mu,buffermpimu1,buffermpimu2,
     $     req2rec,req2send)
c
c     now, call the proper solver.
c
      if(precomp .eq. 1) then
        call inducedc_pme2gpu( pc_dc_tmatxb_pmegpu,nrhs,.true.,ef,mu,
     &                        murec)
      else
        call inducedc_pme2gpu(otf_dc_tmatxb_pmegpu,nrhs,.true.,ef,mu,
     &                        murec)
      end if


      wtime2 = mpi_wtime()
      if (polprt.ge.1.and.rank.eq.0) then
        if (polprt.ge.2) then
          write(iout,1010) 'fields:  ', wtime1-wtime0
          write(iout,1010) 'dipoles: ', wtime2-wtime1
        end if
        write(iout,1020) wtime2 - wtime0
      end if
c
c     move the computed dipoles in the common block.
c
!$acc parallel loop collapse(2) async(def_queue)
      do i = 1, npolebloc
        do j = 1, 3
          iipole = poleglob(i)
          uind(j,iipole) = mu(j,1,i)
          uinp(j,iipole) = mu(j,2,i)
        end do
      end do
!$acc parallel loop collapse(2) async(def_queue)
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
c
c     update the lists of previous induced dipole values
c
      if (use_pred) call pred_SetAlt

      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
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
      use interfaces
      use math
      use mpole
      use polar
      use polpot
      use timestat
      use mpi
      use sizes
      use utils
      use utilcomm ,buffermpi1=>buffermpi3d,buffermpi2=>buffermpi3d1
      use utilgpu
      implicit none
c
c     solves the polarizable dipoles linear equations by jacobi
c     iterations coupled to the DIIS extrapolation.
c
      integer  ,intent(in   ):: nrhs
      logical  ,intent(in   ):: dodiis
      real(t_p),intent(in   ):: ef   (3,nrhs,npolebloc)
      real(t_p),intent(inout):: mu   (3,nrhs,npolebloc)
     &                        , murec(3,nrhs,npolerecloc)
      external matvec

      integer info, proc, knd, pos, mdim, ierror
      integer i, ii, j, k, k1, it, ind, ndismx, nmat, lenb
      integer iglob, iipole, ierr
      integer kst, zst, read0, write0
      real(t_p) zero, one, rnorm(2), rr, xx
      real(t_p) rnorm1,rnorm2
      real(t_p) term
      real(t_p) time0,time1,time2,time3
      parameter (ndismx=25,zero=0,one=1,xx=0)

      integer ipiv(ndismx+1)
      real(t_p), allocatable :: munew(:,:,:),h(:,:,:), cex(:)
     &         , xdiis(:,:), ediis(:,:), bmat(:,:), bloc(:,:)
      real(t_p), allocatable :: dcfield1(:),dcfield2(:)
      real(t_p), allocatable :: dipfield(:,:,:)
     &         , dipfieldbis(:,:,:)
c
c     MPI
c
      integer reqnorm,reqdiis(2*ndismx+1)
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
c
c
 1000 format(' itsolv converged after ',I3,' iterations.',/,
     $       ' final residual norm = ',3D14.7)
 1010 format(' residual norm at iteration ',I3,':',3D12.2)
 1020 format(' Jacobi/DIIS solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Jacobi/DIIS solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
c
      if (rank.eq.0.and.tinkerdebug)
     &   write(*,'(3x,a)') 'inducedc_pme2gpu'

      call prmem_request(buffermpi1,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(buffermpimu1,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(buffermpi2,3,nrhs,max(npolerecloc,1),
     &     async=.true.)
      call prmem_request(buffermpimu2,3,nrhs,max(npolerecloc,1),
     &     async=.true.)
c
      allocate (munew      (3,nrhs,max(1,npolebloc  )))
      allocate (h          (3,nrhs,max(1,npolebloc  )))
      allocate (dipfield   (3,nrhs,max(1,npoleloc   )))
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      allocate (dcfield1   (maxdim))
      allocate (dcfield2   (maxdim))
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec (nproc))
      allocate (reqrec       (nproc))
      allocate (reqsend      (nproc))
      allocate (req2rec      (nproc))
      allocate (req2send     (nproc))
!$acc enter data create(munew,h,dipfield,dipfieldbis,dcfield1,dcfield2,
!$acc&    rnorm1,rnorm2,ipiv)  async
c
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
        call set_to_zero1(bmat,lenb*lenb,rec_queue)
      end if

      !Initiate arrays
      call set_to_zero3(buffermpi1,buffermpimu1,dipfield
     &                 ,size(dipfield),rec_queue)
      call set_to_zero3(buffermpi2,buffermpimu2,dipfieldbis
     &                 ,size(dipfieldbis),rec_queue)
      call set_to_zero2(munew,h,size(h),rec_queue)
      call set_to_zero2(dcfield2,dcfield2,size(dcfield2),rec_queue)
#ifdef _OPENACC
      !start Async overlapping with matvec
      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)
#endif
!$acc data present(klst,mu,poleglob,ipole,atmofst,dipfield,dipfieldbis,
!$acc&   dcfield2,dcfield1,munew,ediis,xdiis) async
c
c     main loop:
c
      do it = 1, politer
         rnorm(1) = zero
         rnorm(2) = zero
!$acc serial async present(rnorm1,rnorm2)
         rnorm1   = zero; rnorm2 = zero;
!$acc end serial
c
         time0 = mpi_wtime()
         call matvec(nrhs,.false.,mu,h)
         time1 = mpi_wtime()
         if (it.eq.1) timerealdip = timerealdip + time1-time0
         call tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
         time2 = mpi_wtime()
         if (it.eq.1) timerecdip = timerecdip + time2-time1

         call commfieldgpu(nrhs,h)

         call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $                       buffermpi2,reqrecdirrec,reqrecdirsend)
c        
         call commdirdirgpu(nrhs,0,mu,reqrec,reqsend)
         call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $                       buffermpi2,reqrecdirrec,reqrecdirsend)
c
         call commrecdirdipgpu(nrhs,0,murec,munew,buffermpimu1,
     $                      buffermpimu2,req2rec,req2send)
c
c     jacobi step:
c
         call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $                       buffermpi2,reqrecdirrec,reqrecdirsend)
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
                       dcfield2(pos+j) = ef(j,2,ii) - 
     $                  (h(j,2,ii) + dipfield(j,2,ii)-term*mu(j,2,ii))
                    else                  ! if (k==1)
                       dcfield1(pos+j) = ef(j,1,ii) - 
     $                  (h(j,1,ii) + dipfield(j,1,ii)-term*mu(j,1,ii))
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
            call cuPOTRS(mdim,zmatDn,mdim,dcfield1,mdim,rec_stream)
            call cuPOTRS(mdim,zmatDn,mdim,dcfield2,mdim,rec_stream)
!$acc end host_data
#else
            call M_PPTRS("L",mdim,1,zmat(kofst(k)+1:knd),
     &                   dcfield1(1:mdim),mdim,ierror)
            if(ierror .ne. 0) print*,"error with solving dipoles"
            call M_PPTRS("L",mdim,1,zmat(kofst(k)+1:knd),
     &                   dcfield2(1:mdim),mdim,ierror)
            if(ierror .ne. 0) print*,"error with solving dipoles"
#endif
c
c  Update dipoles
c
!$acc parallel loop collapse(3) async
            do i = 1, npergrp(k)
               do k1 = 1, 2
                  do j = 1, 3
                     ii    = klst(i,k)
                     iglob = ipole(poleglob(ii))
                     pos   = (atmofst(iglob)-1)*3
                     if (btest(k1,0)) then
                        munew(j,2,ii) = dcfield2(pos+j)
                     else
                        munew(j,1,ii) = dcfield1(pos+j)
                     end if
                  end do
               end do
            end do

            end if
         end do

!$acc parallel loop collapse(3) present(rnorm1,rnorm2) async
         do i = 1,npoleloc
            do k = 1, nrhs
               do j = 1,3
                  if (btest(k,0)) then
                     rnorm2 = rnorm2 + (munew(j,k,i)-mu(j,k,i))**2
                  else
                     rnorm1 = rnorm1 + (munew(j,k,i)-mu(j,k,i))**2
                  end if
               end do
            end do
         end do
!$acc update host(rnorm1,rnorm2) async

         if (nproc>1) then
!$acc wait
            rnorm(1)=rnorm1; rnorm(2)=rnorm2;
            call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_TPREC,
     $                         MPI_SUM,COMM_TINKER,reqnorm,ierr)
         end if

         if (dodiis) then

!$acc parallel loop collapse(3) async
           do i = 1, npoleloc
             do k = 1, nrhs
               do j = 1, 3
                 ind = (i-1)*6 + (k-1)*3 + j
                 xdiis(ind,nmat) = munew(j,k,i)
                 ediis(ind,nmat) = munew(j,k,i) - mu(j,k,i)
               end do
             end do
           end do
c
c          Compute Pulay's Matrix and extrapolate
c
           if(nocomdiis .eq. 1)then
             call no_com_diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $            bmat,nmat)
           else
             reqdiis(:) = MPI_REQUEST_NULL
             call diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $            bmat,nmat,reqdiis,COMM_TINKER)
             !TODO This loop is no longer useful; MPI_IAllReduce is not CUDA-Aware
             do i = 1, 2*nmat-3
                call MPI_WAIT(reqdiis(i),status,ierr)
             end do
           end if

           call utils_amove(lenb*lenb,bmat,bloc,rec_queue)  !copy bmat in bloc
           call set_to_zero1(cex,lenb,rec_queue)
           call set_to_zero1(munew,size(munew),rec_queue)
!$acc serial async
           cex(1) = one
!$acc end serial

#ifdef _OPENACC
!$acc host_data use_device(bloc,ipiv,cex)
           call cuGESV(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,rec_stream)
!$acc end host_data
#else
           call M_dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#endif
           call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)

         end if
         time2 = mpi_wtime()
c
         call commdirdirgpu(nrhs,1,munew,reqrec,reqsend)
         call commrecdirdipgpu(nrhs,1,murec,munew,buffermpimu1,
     $        buffermpimu2,req2rec,req2send)
         call commdirdirgpu(nrhs,2,mu,reqrec,reqsend)
         call commrecdirdipgpu(nrhs,2,murec,munew,buffermpimu1,
     $        buffermpimu2,req2rec,req2send)
         call utils_amove(3*nrhs*npoleloc,munew,mu,rec_queue)   !copy munew in mu
c
c     compute the norm of the increment.
c
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
            if (polprt.ge.1.and.rank.eq.0)
     $         write(6,1000) it, (rnorm(k), k = 1, nrhs)
            goto 10
         end if
         if (polprt.ge.1.and.rank.eq.0)
     $      write(6,1010) it, (rnorm(k), k = 1, nrhs)
      end do

 10   continue
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
c
c     free the memory.
c
!$acc exit data delete(munew,dipfield,dipfieldbis,dcfield1,dcfield2,h)
!$acc&     async
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
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
