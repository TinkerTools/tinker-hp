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
c
#include "tinker_precision.h"
      subroutine dcinduce_pmegpu
      use atmlst
      use divcon
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use units
      use uprior
      use mpi
      use tinheader
      implicit none
c
c     with separate cores for reciprocal part
c
      integer i, j, k, nrhs
c
c     MPI
c
      integer ierr, iipole, proc,tag
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer status(MPI_STATUS_SIZE)
c
      parameter (nrhs=2)
      real(t_p)  wtime0, wtime1, wtime2, udsum, upsum
      real(t_p)  dtime0,dtime1,ctime
      real(t_p)  term, xx(1)
      real(t_p), allocatable :: ef(:,:,:), mu(:,:,:),
     &                                  murec(:,:,:)
      real(t_p), allocatable :: cphi(:,:)
c
      real(t_p), allocatable :: buffermpi(:,:),
     &                      buffermpimu(:,:,:)
c
c     Get direct space induced field for divide and conquer methods
c
      external pc_dc_tmatxb_pme
      external otf_dc_tmatxb_pme
c
      if (.not.use_polar) return
c
 1000 format(' illegal polalg in newinduce.')
 1010 format(' time for the ',a,F14.5)
 1020 format(' total elapsed time in newinduce: ',F14.5)
c
c     allocate some memory and clear the arrays:
c
      call MPI_BARRIER(COMM_TINKER,ierr)
c
      if (rank.le.ndir-1) then
        allocate (mu(3,nrhs,max(1,npolebloc)))
        mu = 0.0_ti_p
        allocate (murec(3,nrhs,max(npolerecloc,1)))
        murec = 0.0_ti_p
        allocate (ef(3,nrhs,max(1,npolebloc)))
        ef = 0.0_ti_p
        allocate (cphi(10,max(1,npoleloc)))
        cphi = 0.0_ti_p
        if (allocated(cphirec)) deallocate (cphirec)
        allocate (cphirec(10,max(npolerecloc,1)))
        cphirec = 0.0_ti_p
        allocate (buffermpi(10,max(npoleloc,1)))
        buffermpi = 0.0_ti_p
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0_ti_p
      else
        allocate (mu(3,nrhs,max(1,npolebloc)))
        mu = 0.0_ti_p
        allocate (murec(3,nrhs,max(npolerecloc,1)))
        murec = 0.0_ti_p
        if (allocated(fphirec)) deallocate (fphirec)
        allocate (fphirec(20,max(npolerecloc,1)))
        fphirec = 0.0_ti_p
        allocate (cphi(10,max(1,npoleloc)))
        cphi = 0.0_ti_p
        if (allocated(cphirec)) deallocate (cphirec)
        allocate (cphirec(10,max(npolerecloc,1)))
        cphirec = 0.0_ti_p
        allocate (buffermpi(10,max(npolerecloc,1)))
        buffermpi = 0.0_ti_p
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
        buffermpimu = 0.0_ti_p
      end if
c
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c     compute the electric fields:
c
      wtime0 = mpi_wtime()
c
c     MPI : begin reception
c
      ctime = 0.0_ti_p
      if (rank.le.ndir-1) then
        call commdirdir(nrhs,0,mu,reqrec,reqsend)
      end if
c
      if (rank.le.ndir-1) then
        call commrecdirfields(0,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
      else
c
c    compute the reciprocal space contribution (fields)
c
        call efld0_recip(cphi)
        call commrecdirfields(1,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)
c
      end if
c
c    The real space processes compute the real fields and  add them to the recip ones
c
      if (rank.le.ndir-1) then
        term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi

        if(clsttype .eq. 1) then
C           kmeans clustering
          call cluster1
        else if(clsttype .eq. 2) then
C           segmentation clustering
          call cluster2
        end if

        if(precomp .eq. 1) then
          call pc_dc_efld0_direct(nrhs,ef)
        else
          call otf_dc_efld0_direct(nrhs,ef)
        end if


C      ctime = ctime + dtime1 - dtime0
c
c     Factorization of Z has been interleaved with passing fields
c

        call commfield2(nrhs,ef)

c
        call commrecdirfields(2,cphirec,cphi,buffermpi,buffermpi,
     $   reqrecdirrec,reqrecdirsend)

c
c       Add direct and reciprocal fields
c
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            ef(j,1,i)  = ef(j,1,i) - cphi(j+1,i) +
     $        term*rpole(j+1,iipole)
            ef(j,2,i)  = ef(j,2,i) - cphi(j+1,i) +
     $        term*rpole(j+1,iipole)
          end do
        end do
      end if
      wtime1 = mpi_wtime()
c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
      if (use_pred.and.nualt.eq.maxualt) then
         call pred_InitField(mu)
      end if

      if (rank.le.ndir-1) then

       if (.not.(use_pred .and. nualt.eq.maxualt)) then

       if (polgsf.eq.0) then
         do i = 1, npoleloc
           iipole = poleglob(i)
           do k = 1, nrhs
             do j = 1, 3
               mu(j,k,i) = polarity(iipole)*ef(j,k,i)
             end do
           end do
         end do
       else
         do i = 1, npoleloc
           iipole = poleglob(i)
           do j = 1, 3
             mu(j,1,i) = uind(j,iipole)
             mu(j,2,i) = uinp(j,iipole)
           end do
         end do
       end if

       end if
       call commdirdir(nrhs,1,mu,reqrec,reqsend)
c
       call commrecdirdip(nrhs,1,murec,mu,buffermpimu,buffermpimu,
     $      req2rec,req2send)
c
       call commdirdir(nrhs,2,mu,reqrec,reqsend)
       call commrecdirdip(nrhs,2,murec,mu,buffermpimu,buffermpimu,
     $      req2rec,req2send)

      else

        call commrecdirdip(nrhs,0,murec,mu,buffermpimu,buffermpimu,
     $       req2rec,req2send)
        call commrecdirdip(nrhs,2,murec,mu,buffermpimu,buffermpimu,
     $       req2rec,req2send)

      end if
c
c     now, call the proper solver.
c
      if(precomp .eq. 1) then
        if (rank.le.ndir-1) then
          call inducedc_pme(pc_dc_tmatxb_pme,nrhs,.true.,ef,mu,xx)
        else
          call inducedc_pme(pc_dc_tmatxb_pme,nrhs,.true.,xx,xx,murec)
        end if
      else
        if (rank.le.ndir-1) then
          call inducedc_pme(otf_dc_tmatxb_pme,nrhs,.true.,ef,mu,xx)
        else
          call inducedc_pme(otf_dc_tmatxb_pme,nrhs,.true.,xx,xx,murec)
        end if
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
      if (rank.le.ndir-1) then
        do i = 1, npolebloc
          iipole = poleglob(i)
          do j = 1, 3
            uind(j,iipole) = mu(j,1,i)
            uinp(j,iipole) = mu(j,2,i)
          end do
        end do
      else
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          do j = 1, 3
            uind(j,iipole) = murec(j,1,i)
            uinp(j,iipole) = murec(j,2,i)
          end do
        end do
      end if
      deallocate (buffermpi)
      deallocate (buffermpimu)
      if (rank.le.ndir-1) then
        deallocate (ef)
        deallocate (cphi)
        deallocate (mu)
      else
c        deallocate (fphi)
c        deallocate (cphirec)
        deallocate (murec)
      end if
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
      deallocate (reqrecdirrec)
      deallocate (reqrecdirsend)
c
c     update the lists of previous induced dipole values
c
      if (use_pred) call pred_SetAlt

      end


      subroutine inducedc_pmegpu(matvec,nrhs,dodiis,ef,mu,murec)
      use atmlst
      use divcon
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use timestat
      use mpi
      use tinheader
      implicit none
c
c     solves the polarizable dipoles linear equations by jacobi
c     iterations coupled to the DIIS extrapolation.
c
      integer fnd,mdim,inf,iglob,pos,knd,nrhs,info
      real(t_p)  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      real(t_p), allocatable :: munew(:,:,:), h(:,:,:)
      real(t_p), allocatable :: xdiis(:,:), ediis(:,:),
     $   bmat(:,:),bloc(:,:), cex(:)
      integer i,ii, j, k, iipole, it, ind, ndismx, nmat, lenb
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real(t_p)  zero, one, rnorm(2), rr, xx(1)
      real(t_p) term
      real(t_p) time0,time1,time2,time3
      real(t_p)  dtime0,dtime1,ctime
      save    zero, one, xx
      external matvec
      real(t_p), allocatable :: dcfield1(:),dcfield2(:)
      real(t_p), allocatable :: dipfield(:,:,:)
      real(t_p), allocatable :: dipfieldbis(:,:,:)
c
c     MPI
c
      real(t_p), allocatable ::buffermpi1(:,:,:),
     &                         buffermpi2(:,:,:)
      real(t_p), allocatable ::buffermpimu(:,:,:)
      integer ierr, tag, proc
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer reqnorm,reqdiis(2*ndismx+1)
      integer status(MPI_STATUS_SIZE)
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
c
      zero   = 0.0d0
      one    = 1.0d0
      xx(1)  = 0.0d0

      if (rank.le.ndir-1) then
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0_ti_p
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0_ti_p
        allocate (dipfield(3,nrhs,max(1,npoleloc)))
        dipfield = 0.0_ti_p
        allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
        dipfieldbis = 0.0_ti_p
        allocate (dcfield1(maxdim))
        allocate (dcfield2(maxdim))
        dcfield1 = 0.0_ti_p
        dcfield2 = 0.0_ti_p
        allocate (h(3,nrhs,npolebloc))
        h = 0.0_ti_p
        if (dodiis) then
          nmat = 1
          lenb = ndismx + 1
          allocate (xdiis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (ediis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (bmat(lenb,lenb))
          allocate (bloc(lenb,lenb))
          allocate (cex(lenb))
          bmat = 0.0_ti_p
        end if
      else
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0_ti_p
        allocate (buffermpimu(3,nrhs,max(1,npolerecloc)))
        buffermpimu = 0.0_ti_p
        allocate (dipfield(3,nrhs,max(1,npoleloc)))
        dipfield = 0.0_ti_p
        allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
        dipfieldbis = 0.0_ti_p
      end if
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0.0_ti_p
      allocate (buffermpi2(3,nrhs,max(1,npolerecloc)))
      buffermpi2 = 0.0_ti_p
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
c
c     main loop:
c
      do it = 1, politer
C        ctime = 0.0_ti_p
        rnorm = 0.0_ti_p
c
c     MPI : begin reception
c
        if (rank.le.ndir-1) then
          call commdirdir(nrhs,0,mu,reqrec,reqsend)
        else
          call commrecdirdip(nrhs,0,dipfieldbis,dipfield,
     $      buffermpimu,buffermpimu,req2rec,req2send)
        end if
c
c      Begin the reception of the reciprocal fields
c
        if (rank.le.ndir-1) then
           call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $     buffermpi2,reqrecdirrec,reqrecdirsend)
c
        else
          time0 = mpi_wtime()
          call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
          time1 = mpi_wtime()
          if (it.eq.1) timerecdip = timerecdip + time1-time0
          call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
          call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)

        end if
c
        if (rank.le.ndir-1) then
c     jacobi step:
c
c    The real space processes extract the recip fields, compute the real fields
c    and add them
c
         term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
         time0 = mpi_wtime()
         call matvec(nrhs,.false.,mu,h)
         time1 = mpi_wtime()
         if (it.eq.1) timerealdip = timerealdip + time1-time0
         call commfield(nrhs,h)
c
         call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $      buffermpi2,reqrecdirrec,reqrecdirsend)
c
          do k = 1,km
            mdim = npergrp(k)*3
            if(mdim .gt. 0)then
            do i = 1, npergrp(k)
              ii = klst(i,k)
              iipole = poleglob(ii)
              iglob = ipole(iipole)
C
C  Build the field for each block
C
              do j = 1,3
                pos=(atmofst(iglob)-1)*3
                dcfield1(pos+j) = ef(j,1,ii) - 
     $(h(j,1,ii) + dipfield(j,1,ii)-term*mu(j,1,ii))
                dcfield2(pos+j) = ef(j,2,ii) - 
     $(h(j,2,ii) + dipfield(j,2,ii)-term*mu(j,2,ii))
              end do
            end do

            knd = kofst(k) + mdim*(mdim+1)/2 
#ifdef _OPENACC
            call fatal_device("inducedc_pmegpu")
#else
            call DPPTRS("L",mdim,1,zmat(kofst(k)+1:knd),
     &    dcfield1(1:mdim),mdim,inf)
            if(inf .ne. 0) print*,"error with solving dipoles"
            call DPPTRS("L",mdim,1,zmat(kofst(k)+1:knd),
     &    dcfield2(1:mdim),mdim,inf)
            if(inf .ne. 0) print*,"error with solving dipoles"
#endif

C
C  Update dipoles
C
            do i = 1, npergrp(k)
              ii = klst(i,k)
              iipole = poleglob(ii)
              iglob = ipole(iipole)
              do j = 1,3
                pos=(atmofst(iglob)-1)*3
                munew(j,1,ii) = dcfield1(pos+j)
                munew(j,2,ii) = dcfield2(pos+j)
              end do
            end do

            end if
          end do

          do k = 1, nrhs
            rnorm(k)=sum((munew(:,k,1:npoleloc)-mu(:,k,1:npoleloc))**2)
          end do
c
        end if
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_TPREC,
     $    MPI_SUM,COMM_TINKER,reqnorm,ierr)

        if (rank.le.ndir-1) then
          if (dodiis) then
            ind = 0
            do i = 1, npoleloc
              do k = 1, nrhs
                do j = 1, 3
                  ind = ind + 1
                  xdiis(ind,nmat) = munew(j,k,i)
                  ediis(ind,nmat) = munew(j,k,i) - mu(j,k,i)
                end do
              end do
            end do
c
c           Compute Pulay's Matrix and extrapolate
c
            if(nocomdiis .eq. 1)then
              call no_com_diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $        bmat,nmat)
            else
              call diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $        bmat,nmat,reqdiis,comm_dir)
              do i = 1, 2*nmat-3
                call MPI_WAIT(reqdiis(i),status,ierr)
              end do
            end if
            bloc = bmat
            cex = 0.0_ti_p
            cex(1) = one
#ifdef _OPENACC
            call fatal_device("inducedc_pmegpu")
#else
            call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#endif
            munew = 0.0_ti_p
            call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
          end if
          time2 = mpi_wtime()
          call commdirdir(nrhs,1,munew,reqrec,reqsend)
          call commrecdirdip(nrhs,1,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
c
          time0 = mpi_wtime()
          call commdirdir(nrhs,2,mu,reqrec,reqsend)
          call commrecdirdip(nrhs,2,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
          mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)
          time1 = mpi_wtime()
          time3 = mpi_wtime()
          mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)

        else
          call commrecdirdip(nrhs,2,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
        end if
c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/real(3*npolar,t_p))
c          rnorm(k) = sqrt(rnorm(k))
          rr = max(rnorm(k),rr)
        end do
        if (polprt.ge.2.and.rank.eq.0)
     $         write(6,1010) it, (rnorm(k), k = 1, nrhs)
        if (rr.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1000) it, (rnorm(k), k = 1, nrhs)
          goto 10
        end if
      end do
 10   continue
c
      if (polprt.ge.3) then
        write(iout,1020)
        do i = 1, npoleloc
          iipole = poleglob(i)
          write(iout,1030) iipole, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npoleloc
          iipole = poleglob(i)
          write(iout,1030) iipole, (mu(j,2,i), j = 1, 3)
        end do
      end if
c
c     free the memory.
c
      deallocate (buffermpimu)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqrecdirsend)
      deallocate (reqrecdirrec)
      if (rank.le.ndir-1) then
        deallocate (buffermpi1)
        deallocate (munew)
        deallocate (h)
        deallocate (dipfield)
        if (dodiis) then
          deallocate (xdiis)
          deallocate (ediis)
          deallocate (bmat)
          deallocate (bloc)
          deallocate (cex)
        end if
      else
        deallocate (buffermpi2)
        deallocate (dipfieldbis)
      end if
      return
      end


      subroutine dc_factorgpu
c
c     compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use divcon
      use domdec
#ifdef _OPENACC
      use interfaces ,only: cuPOTRF
      use utilgpu    ,only: rec_stream
#endif
      use utilgpu    ,only: prmem_request
      use utils
      implicit none
      integer k,mdim,knd,kbeg,inf
c     loop over the k blocks and perform a cholesky factorization
#ifdef _OPENACC
      integer i, j, read0, zst

      call prmem_request(zmatDn,maxdim**2,async=.true.)
#endif

!$acc data present(zmatDn,zmat,npergrp,kofst) async
      do k = 1, km
        mdim = npergrp(k)*3
        if (mdim.gt.0) then
          knd  = kofst(k) + mdim*(mdim+1)/2 
          kbeg = kofst(k) + 1
#ifdef _OPENACC
!$acc parallel loop vector_length(32) async
          do i = 0,mdim-1         ! Zmat packed to zmatDn Dense
             zst  = i*mdim - (i*(i-1))/2
!$acc loop vector
             do j = i,mdim-1
                read0  = zst + j-i
                !print*, i,j,mdim,kbeg,read0
                zmatDn(i*mdim + j+1) = zmat(kbeg+read0)
             end do
          end do

!$acc host_data use_device(zmatDn)
          call cuPOTRF(mdim,zmatDn,mdim,rec_stream)
!$acc end host_data
c         call DPOTRF("L",mdim,zmatDn,mdim,inf)

!$acc parallel loop vector_length(32) async
          do i = 0,mdim-1           ! Write result back in packed format of zmat
             zst = mdim*i - (i*(i-1))/2
!$acc loop vector
             do j = i,mdim-1
                zmat(kbeg+zst + j-i) = zmatDn(i*mdim + j+1)
             end do
          end do
#else
          call M_PPTRF("L",mdim,zmat(kbeg:knd),inf)
          if(inf .ne. 0)then
          print*,rank,k,"error ",inf,"with factorization of zmat"
          end if
#endif
        end if
      end do
!$acc end data
      end

      subroutine pc_dc_tmatxb_pmegpu(nrhs,dodiag,mu,efi)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use atmlst
      use atoms
      use divcon
      use domdec
      use ewald
      use group
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use sizes
      use shunt
      use tinheader
      implicit none
      integer i,nrhs,iglob,kglob,kbis,iipole,kkpole,kkpoleloc
      integer ipoleloc,l
      real(t_p)  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      logical dodiag
      integer j, ii, kkk, irhs, tpos
      real(t_p) dx, dy, dz, d2
      real(t_p) dukx, duky, dukz, pukx, puky, pukz
      real(t_p) duix, duiy, duiz, puix, puiy, puiz
      real(t_p)  zero, one, f50
      real(t_p)  cutoff2
      save    zero, one, f50
      character*10 mode

      if (rank.eq.0.and.tinkerdebug)
     &   write(*,'(4x,a)') 'pc_dc_tmatxb_pmegpu'

      zero = 0.0d0
      one  = 1.0d0
      f50  = 50.d0
c
c     initialize the result vector
c
      do i = 1, npolebloc
        do irhs = 1, nrhs
          do j = 1, 3
            efi(j,irhs,i) = zero
          end do
        end do
      end do

c
c     gather some parameters, then set up the damping factors.
c
      mode = 'EWALD'
      call switch (mode)
      cutoff2 = cut2
      tpos = 1
c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        l = grplst(iglob)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        if (i.eq.0) cycle
        duix = mu(1,1,ipoleloc)
        duiy = mu(2,1,ipoleloc)
        duiz = mu(3,1,ipoleloc)
        puix = mu(1,2,ipoleloc)
        puiy = mu(2,2,ipoleloc)
        puiz = mu(3,2,ipoleloc)
c
        do kkk = 1,nelst(ii)
          kglob = elst(kkk,ii)
c
c     only build induced field for interactons outside of a block
c     check if the atoms are in the same block
c
          if(l .ne. grplst(kglob) .or. l .lt. 0)then
          kkpole = pollist(kglob)
          kbis = loc(kglob)
          kkpoleloc = poleloc(kkpole)
          if (kkpoleloc.eq.0) cycle
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          call image(dx,dy,dz)
          d2 = dx*dx + dy*dy + dz*dz
          if (d2.le.cutoff2) then
            dukx = mu(1,1,kkpoleloc)
            duky = mu(2,1,kkpoleloc)
            dukz = mu(3,1,kkpoleloc)
            pukx = mu(1,2,kkpoleloc)
            puky = mu(2,2,kkpoleloc)
            pukz = mu(3,2,kkpoleloc)

            efi(1,1,ipoleloc) = efi(1,1,ipoleloc)
     $  + ytab(tpos)*dukx + ytab(tpos+1)*duky + ytab(tpos+2)*dukz
            efi(2,1,ipoleloc) = efi(2,1,ipoleloc)
     $  + ytab(tpos+1)*dukx + ytab(tpos+3)*duky + ytab(tpos+4)*dukz
            efi(3,1,ipoleloc) = efi(3,1,ipoleloc)
     $  + ytab(tpos+2)*dukx + ytab(tpos+4)*duky + ytab(tpos+5)*dukz

            efi(1,2,ipoleloc) = efi(1,2,ipoleloc)
     $  + ytab(tpos)*pukx + ytab(tpos+1)*puky + ytab(tpos+2)*pukz
            efi(2,2,ipoleloc) = efi(2,2,ipoleloc)
     $  + ytab(tpos+1)*pukx + ytab(tpos+3)*puky + ytab(tpos+4)*pukz
            efi(3,2,ipoleloc) = efi(3,2,ipoleloc)
     $  + ytab(tpos+2)*pukx + ytab(tpos+4)*puky + ytab(tpos+5)*pukz
c
            efi(1,1,kkpoleloc) = efi(1,1,kkpoleloc)
     $  + ytab(tpos)*duix + ytab(tpos+1)*duiy + ytab(tpos+2)*duiz
            efi(2,1,kkpoleloc) = efi(2,1,kkpoleloc) 
     $  + ytab(tpos+1)*duix + ytab(tpos+3)*duiy + ytab(tpos+4)*duiz
            efi(3,1,kkpoleloc) = efi(3,1,kkpoleloc)
     $  + ytab(tpos+2)*duix + ytab(tpos+4)*duiy + ytab(tpos+5)*duiz

            efi(1,2,kkpoleloc) = efi(1,2,kkpoleloc)
     $  + ytab(tpos)*puix + ytab(tpos+1)*puiy + ytab(tpos+2)*puiz
            efi(2,2,kkpoleloc) = efi(2,2,kkpoleloc)
     $  + ytab(tpos+1)*puix + ytab(tpos+3)*puiy + ytab(tpos+4)*puiz
            efi(3,2,kkpoleloc) = efi(3,2,kkpoleloc)
     $  + ytab(tpos+2)*puix + ytab(tpos+4)*puiy + ytab(tpos+5)*puiz
            tpos = tpos + 6
            end if
          end if
        end do

      end do
      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
        do i = 1, npoleloc
          iipole = poleglob(i)
          if (polarity(iipole) .eq. 0_ti_p) then
             cycle
          else
            do irhs = 1, nrhs
              do j = 1, 3
                efi(j,irhs,i) = efi(j,irhs,i) +
     $             mu(j,irhs,i)/polarity(iipole)
              end do
            end do
          end if
        end do
      end if

      end

c
c     subroutine commfield : communicate some direct fields (Newton third law)
c
      subroutine commfield2gpu(nrhs,ef)
      use atoms
      use domdec
      use mpole
      use potent
      use mpi
      use utilgpu ,only: prmem_request,dir_queue
      use utilcomm,only: buff_field
      implicit none
      integer nrhs,i,j,k,l,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer req(nproc*nproc)
      real(t_p) ef(3,nrhs,*)
c
      if (use_pmecore) then
         commloc = comm_dir
      else
         commloc = COMM_TINKER
      end if
c
      if (n_send1>0)
     &   call prmem_request(buff_field,3,nrhs,
     &            max(npoleloc,1),n_send1,async=.true.)
c
!$acc host_data use_device(buff_field)
      do i = 1, n_send1
        tag = nproc*rank + p_send1(i) + 1
        call MPI_IRECV(buff_field(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,
     $       p_send1(i),tag,commloc,req(tag),ierr)
      end do
!$acc end host_data
c
!$acc host_data use_device(ef)
      do i = 1, n_recep1
!$acc wait(dir_queue)
        tag = nproc*p_recep1(i) + rank + 1
        call MPI_ISEND(ef(1,1,bufbegpole(p_recep1(i)+1)),
     $       3*nrhs*domlen(p_recep1(i)+1),MPI_TPREC,p_recep1(i),
     $       tag,commloc,req(tag),ierr)
      end do
!$acc end host_data

      call dc_factorgpu
c
      do i = 1, n_send1
        tag = nproc*rank + p_send1(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, n_recep1
        tag = nproc*p_recep1(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
!$acc parallel loop collapse(4) present(ef,buff_field) async
      do i = 1, n_send1
         do j = 1, npoleloc
            do k = 1, nrhs
               do l = 1,3
                  ef(l,k,j) = ef(l,k,j) + buff_field(l,k,j,i)
               end do
            end do
         end do
      end do
      end

c
c     subroutine commfield : communicate some direct fields (Newton third law)
c
      subroutine commfield2shortgpu(nrhs,ef)
      use atoms
      use domdec
      use mpole
      use potent
      use mpi
      use utilgpu ,only: prmem_request,dir_queue
      use utilcomm,only: buff_field
      implicit none
      integer nrhs,i,j,k,l,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer req(nproc*nproc)
      real(t_p) ef(3,nrhs,*)
      real(t_p), allocatable :: buffer(:,:,:,:)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if

      if (n_send1>0)
     &   call prmem_request(buff_field,3,nrhs,
     &          max(npoleloc,1),n_send1,async=.true.)
c
      do i = 1, n_sendshort1
        tag = nproc*rank + p_sendshort1(i) + 1
        call MPI_IRECV(buff_field(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,
     $       p_sendshort1(i),tag,commloc,req(tag),ierr)
      end do
c
      do i = 1, n_recepshort1
        tag = nproc*p_recepshort1(i) + rank + 1
        call MPI_ISEND(ef(1,1,bufbegpole(p_recepshort1(i)+1)),
     $       3*nrhs*domlen(p_recepshort1(i)+1),MPI_TPREC,
     $       p_recepshort1(i),tag,commloc,req(tag),ierr)
      end do

      call dc_factorgpu
c
      do i = 1, n_sendshort1
        tag = nproc*rank + p_sendshort1(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, n_recepshort1
        tag = nproc*p_recepshort1(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
      do i = 1, n_sendshort1
        do j = 1, npoleloc
          do k = 1, nrhs
            do l = 1,3
               ef(l,k,j) = ef(l,k,j) + buff_field(l,k,j,i)
            end do
          end do
        end do
      end do
c
      end

