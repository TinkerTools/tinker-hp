c      
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c      
      subroutine dcinduce_pme2
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
      implicit none
c
c     evaluate induced dipole moments and the polarization energy
c     using a DC-JI/DIIS algorithm
c      
      integer i, j, k, nrhs, proc
c
c     MPI
c
      integer iglob, ierr, iipole
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer tag,status(MPI_STATUS_SIZE)
      real*8 time0,time1
c
      parameter (nrhs=2)
      real*8  wtime0, wtime1, wtime2, omp_get_wtime, udsum, upsum
      real*8  term, xx(1)
      real*8, allocatable :: ef(:,:,:), mu(:,:,:), murec(:,:,:)
      real*8, allocatable :: cphi(:,:)
c
      real*8, allocatable :: buffermpi1(:,:),buffermpi2(:,:)
      real*8, allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
c
c     Get direct space induced field for divide and conquer methods
c
      external pc_dc_tmatxb_pme2
      external otf_dc_tmatxb_pme2
c
      if (.not.use_polar) return
c
 1000 format(' illegal polalg in newinduce.')
 1010 format(' time for the ',a,F14.5)
 1020 format(' total elapsed time in newinduce: ',F14.5)
c
c     allocate some memory and clear the arrays:
c
      allocate (mu(3,nrhs,max(1,npolebloc)))
      mu = 0d0
      allocate (murec(3,nrhs,max(1,npolerecloc)))
      murec = 0d0
c
      allocate (buffermpi1(10,max(npoleloc,1)))
      buffermpi1 = 0d0
      allocate (buffermpi2(10,max(npolerecloc,1)))
      buffermpi2 = 0d0
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0d0
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0d0
c
      allocate (ef(3,nrhs,max(1,npolebloc)))
      ef = 0d0
      allocate (cphi(10,max(npoleloc,1)))
      cphi = 0d0
      if (allocated(cphirec)) deallocate (cphirec)
      allocate (cphirec(10,max(npolerecloc,1)))
      cphirec = 0d0
      if (allocated(fphirec)) deallocate(fphirec)
      allocate (fphirec(20,max(npolerecloc,1)))
      fphirec = 0d0
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
c$    wtime0 = omp_get_wtime()
c
c    compute the reciprocal space contribution (fields)
c
      call efld0_recip2(cphi)
c
c     Begin the reception of the reciprocal fields
c
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_IRECV(buffermpi1(1,bufbeg1(proc+1)),
     $      10*buflen1(proc+1),
     $      MPI_REAL8,proc,tag,COMM_BEAD,reqrecdirrec(i),ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          do j = 0, buflen2(proc+1)-1
            call amove(10,cphirec(1,polerecloc(
     $        buf2(bufbeg2(proc+1)+j))),buffermpi2(1,bufbeg2(proc+1)+j))
          end do
        end if
      end do
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_ISEND(buffermpi2(1,bufbeg2(proc+1)),
     $     10*buflen2(proc+1),MPI_REAL8,proc,tag,COMM_BEAD,
     $     reqrecdirsend(i),ierr)
        end if
      end do
c
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          call MPI_WAIT(reqrecdirsend(i),status,ierr)
        end if
      end do
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          call MPI_WAIT(reqrecdirrec(i),status,ierr)
          do j = 0, buflen1(proc+1)-1
            call amove(10,buffermpi1(1,bufbeg1(proc+1)+j),
     $        cphi(1,poleloc(buf1(bufbeg1(proc+1)+j))))
          end do
        end if
      end do
c
      if(clsttype .eq. 1) then
C         kmeans clustering
        call cluster1
      else if(clsttype .eq. 2) then
C         segmentation clustering
        call cluster2
      end if

      if(precomp .eq. 1) then
        call pc_dc_efld0_direct(nrhs,ef)
      else
        call otf_dc_efld0_direct(nrhs,ef)
      end if
c
c     Factorization of Z has been interleaved with passing fields
c
      call commfield2(nrhs,ef)
c
c     MPI : begin reception
c
      do i = 1, n_recep1
        tag = nproc*rank + p_recep1(i) + 1
        call MPI_IRECV(mu(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),tag,
     $   COMM_BEAD,reqrec(i),ierr)
      end do
c
c       Add direct and reciprocal fields
c
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npoleloc
        iipole = poleglob(i)
        do j = 1, 3
          ef(j,1,i)  = ef(j,1,i) - cphi(j+1,i) +
     $       term*rpole(j+1,iipole)
          ef(j,2,i)  = ef(j,2,i) - cphi(j+1,i) +
     $       term*rpole(j+1,iipole)
        end do
      end do
c$    wtime1 = omp_get_wtime()

c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
      if (use_pred .and. nualt.eq.maxualt) then
        call ulspred
        do i = 1, npoleloc
          iipole = poleglob(i)
          do j = 1, 3
            udsum = 0.0d0
            upsum = 0.0d0
            do k = 1, nualt - 1
              udsum = udsum + bpred(k)*udalt(k,j,iipole)
              upsum = upsum + bpred(k)*upalt(k,j,iipole)
            end do
            mu(j,1,i) = udsum
            mu(j,2,i) = upsum
          end do
        end do
      else if (polgsf.eq.0) then
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
c
c     MPI : begin sending
c
      do i = 1, n_send1
        tag = nproc*p_send1(i) + rank + 1
        call MPI_ISEND(mu,3*nrhs*npoleloc,
     $   MPI_REAL8,p_send1(i),tag,COMM_BEAD,
     $   reqsend(i),ierr)
      end do
c
      do i = 1, n_recep1
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, n_send1
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     Begin reception of mu for PME
c
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),
     $      3*nrhs*buflen2(proc+1),MPI_REAL8,proc,tag,COMM_BEAD,
     $      req2rec(i),ierr)
        end if
      end do
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          do j = 0, buflen1(proc+1)-1
            call amove(3*nrhs,mu(1,1,poleloc(buf1(bufbeg1(proc+1)+j))),
     $        buffermpimu1(1,1,bufbeg1(proc+1)+j))
          end do
        end if
      end do
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_ISEND(buffermpimu1(1,1,bufbeg1(proc+1)),
     $     3*nrhs*buflen1(proc+1),MPI_REAL8,proc,tag,COMM_BEAD,
     $     req2send(i),ierr)
        end if
      end do
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          call MPI_WAIT(req2send(i),status,ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          call MPI_WAIT(req2rec(i),status,ierr)
          do j = 0, buflen2(proc+1)-1
            call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),
     $        murec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
          end do
        end if
      end do

c
c     now, call the proper solver.
c
      if(precomp .eq. 1) then
        call inducedc_pme2(pc_dc_tmatxb_pme2,nrhs,.true.,ef,mu,murec)
      else
        call inducedc_pme2(otf_dc_tmatxb_pme2,nrhs,.true.,ef,mu,murec)
      end if


c$    wtime2 = omp_get_wtime()
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
      do i = 1, npolebloc
        iipole = poleglob(i)
        do j = 1, 3
          uind(j,iipole) = mu(j,1,i)
          uinp(j,iipole) = mu(j,2,i)
        end do
      end do
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        do j = 1, 3
          if (repart(iglob).ne.rank) then
            uind(j,iipole) = murec(j,1,i)
            uinp(j,iipole) = murec(j,2,i)
          else
            uind(j,iipole) = mu(j,1,poleloc(iipole))
            uinp(j,iipole) = mu(j,2,poleloc(iipole))
          end if
        end do
      end do
c
c     update the lists of previous induced dipole values
c
      if (use_pred) then
         nualt = min(nualt+1,maxualt)
         do i = 1, npolebloc
           iipole = poleglob(i)
            do j = 1, 3
               do k = nualt, 2, -1
                  udalt(k,j,iipole) = udalt(k-1,j,iipole)
                  upalt(k,j,iipole) = upalt(k-1,j,iipole)
               end do
               udalt(1,j,iipole) = uind(j,iipole)
               upalt(1,j,iipole) = uinp(j,iipole)
             end do
         end do
      end if
      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
      deallocate (reqrecdirrec)
      deallocate (reqrecdirsend)
      deallocate (buffermpi1)
      deallocate (buffermpi2)
      deallocate (buffermpimu1)
      deallocate (buffermpimu2)
      deallocate (ef)
c      deallocate (fphi)
      deallocate (mu)
      deallocate (murec)
      deallocate (cphi)
c      deallocate (cphirec)
      return
      end

      subroutine inducedc_pme2(matvec,nrhs,dodiis,ef,mu,murec)
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
      implicit none
c
c     solves the polarizable dipoles linear equations by jacobi
c     iterations coupled to the DIIS extrapolation.
c
      integer nrhs, info, proc,knd,pos,mdim,inf
      real*8  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      real*8, allocatable :: munew(:,:,:),h(:,:,:)
      real*8, allocatable :: xdiis(:,:), ediis(:,:),
     $   bmat(:,:),bloc(:,:), cex(:)
      integer i,ii, j, k, it, ind, ndismx, nmat, lenb, tag
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real*8  zero, one, rnorm(2), rr, xx(1)
      real*8 term
      real*8 time0,time1,time2,time3
      save    zero, one, xx
      data    zero/0.0d0/, one/1.0d0/, xx/0.0d0/
      external matvec
      real*8, allocatable :: dcfield1(:),dcfield2(:)
      real*8, allocatable :: dipfield(:,:,:)
      real*8, allocatable :: dipfieldbis(:,:,:)
c
c     MPI
c
      real*8, allocatable :: buffermpi1(:,:,:),buffermpi2(:,:,:)
      real*8, allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
      integer iglob, iipole, ierr
      integer reqnorm,reqdiis(2*ndismx+1)
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
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
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0d0
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0d0
      allocate (buffermpi2(3,nrhs,max(npolerecloc,1)))
      buffermpi2 = 0d0
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0d0
c
      allocate (munew(3,nrhs,max(1,npolebloc)))
      munew = 0d0
      allocate (dipfield(3,nrhs,max(1,npoleloc)))
      dipfield = 0d0
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      dipfieldbis = 0d0
      allocate (dcfield1(maxdim))
      allocate (dcfield2(maxdim))
      dcfield1 = 0.0d0
      dcfield2 = 0.0d0
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
c
      allocate (h(3,nrhs,max(1,npolebloc)))
      h = 0d0
      if (dodiis) then
        nmat = 1
        lenb = ndismx + 1
        allocate (xdiis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (ediis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (bmat(lenb,lenb))
        allocate (bloc(lenb,lenb))
        allocate (cex(lenb))
        bmat = 0d0
      end if

c
c     main loop:
c
      do it = 1, politer
        do k = 1, nrhs
          rnorm(k) = zero
        end do
c
        time0 = mpi_wtime()
        call tmatxbrecip2(mu,murec,nrhs,dipfield,dipfieldbis)
        time1 = mpi_wtime()
        if (it.eq.1) timerecdip = timerecdip + time1-time0
        call matvec(nrhs,.false.,mu,h)
        time2 = mpi_wtime()
        if (it.eq.1) timerealdip = timerealdip + time2-time1
        call commfield(nrhs,h)
c
c     Begin the reception of the reciprocal fields
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
           tag = nproc*rank + proc + 1
           call MPI_IRECV(buffermpi1(1,1,bufbeg1(proc+1)),
     $       3*nrhs*buflen1(proc+1),
     $       MPI_REAL8,proc,tag,COMM_BEAD,reqrecdirrec(i),ierr)
          end if
        end do
c
c     MPI : begin reception
c
        do i = 1, n_recep1
          tag = nproc*rank + p_recep1(i) + 1
        call MPI_IRECV(mu(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),tag,
     $   COMM_BEAD,reqrec(i),ierr)
        end do
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            do j = 0, buflen2(proc+1)-1
              call amove(3*nrhs,dipfieldbis(1,1,polerecloc(
     $          buf2(bufbeg2(proc+1)+j))),
     $          buffermpi2(1,1,bufbeg2(proc+1)+j))
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpi2(1,1,bufbeg2(proc+1)),
     $       3*nrhs*buflen2(proc+1),
     $       MPI_REAL8,proc,tag,COMM_BEAD,reqrecdirsend(i),ierr)
          end if
        end do
c
c     Begin reception of mu for PME
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),3*nrhs*
     $       buflen2(proc+1),MPI_REAL8,proc,tag,COMM_BEAD,
     $       req2rec(i),ierr)
          end if
        end do
c
c     jacobi step:
c
c       Wait for the reciprocal fields
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirsend(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(reqrecdirrec(i),status,ierr)
            do j = 0, buflen1(proc+1)-1
              call amove(3*nrhs,buffermpi1(1,1,bufbeg1(proc+1)+j),
     $          dipfield(1,1,poleloc(buf1(bufbeg1(proc+1)+j))))
            end do
          end if
        end do
        time2 = mpi_wtime()
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi

          do k = 1,km
            mdim = npergrp(k)*3
            if(mdim .gt. 0)then
            ind = 0
            do i = 1, npergrp(k)
              ii = klst(i,k)
              iipole = poleglob(ii)
              iglob = ipole(iipole)
C
C  Build the field for each block
C
              do j = 1,3
                ind = ind + 1
                pos=(atmofst(iglob)-1)*3
                dcfield1(pos+j) = ef(j,1,ii) - 
     $(h(j,1,ii) + dipfield(j,1,ii)-term*mu(j,1,ii))
                dcfield2(pos+j) = ef(j,2,ii) - 
     $(h(j,2,ii) + dipfield(j,2,ii)-term*mu(j,2,ii))
              end do
            end do

            knd = kofst(k) + mdim*(mdim+1)/2 
            call DPPTRS("L",mdim,1,zmat(kofst(k)+1:knd),
     &    dcfield1(1:mdim),mdim,inf)
            if(inf .ne. 0) print*,"error with solving dipoles"
            call DPPTRS("L",mdim,1,zmat(kofst(k)+1:knd),
     &    dcfield2(1:mdim),mdim,inf)
            if(inf .ne. 0) print*,"error with solving dipoles"

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
          rnorm(k) = sum((munew(:,k,1:npoleloc)-mu(:,k,1:npoleloc))**2)
        end do


        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_REAL8,
     $    MPI_SUM,COMM_BEAD,reqnorm,ierr)

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
c         Compute Pulay's Matrix and extrapolate
c
          if(nocomdiis .eq. 1)then
            call no_com_diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $      bmat,nmat)
          else
            call diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $        bmat,nmat,reqdiis,comm_dir)
            do i = 1, 2*nmat-3
              call MPI_WAIT(reqdiis(i),status,ierr)
            end do
          end if
          bloc = bmat
          cex = 0d0
          cex(1) = one
          call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
          munew = 0d0
          call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
        end if
        time2 = mpi_wtime()
c
c     MPI : begin sending
c
        do i = 1, n_send1
          tag = nproc*p_send1(i) + rank + 1
          call MPI_ISEND(munew,3*nrhs*npoleloc,
     $     MPI_REAL8,p_send1(i),tag,COMM_BEAD,
     $     reqsend(i),ierr)
        end do
c
c     MPI : begin sending of mu for PME
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            do j = 0, buflen1(proc+1)-1
              call amove(3*nrhs,munew(1,1,poleloc(
     $          buf1(bufbeg1(proc+1)+j))),
     $          buffermpimu1(1,1,bufbeg1(proc+1)+j))
            end do
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpimu1(1,1,bufbeg1(proc+1)),3*nrhs*
     $       buflen1(proc+1),MPI_REAL8,proc,tag,COMM_BEAD,
     $       req2send(i),ierr)
          end if
        end do
       do i = 1, n_recep1
          call MPI_WAIT(reqrec(i),status,ierr)
        end do
        do i = 1, n_send1
          call MPI_WAIT(reqsend(i),status,ierr)
        end do
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_WAIT(req2send(i),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            call MPI_WAIT(req2rec(i),status,ierr)
            do j = 0, buflen2(proc+1)-1
              call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),
     $          murec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
            end do
          end if
        end do
        mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)
c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/dble(3*npolar))
          rr = max(rnorm(k),rr)
        end do
        if (rr.lt.poleps .and. it .ge. 3) then
          if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1000) it, (rnorm(k), k = 1, nrhs)
          goto 10
        end if
        if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1010) it, (rnorm(k), k = 1, nrhs)
      end do
 10   continue
      if (polprt.ge.3) then
        write(iout,1020)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = ipole(iipole)
          write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = poleglob(iipole)
          write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
        end do
      end if
c
c     free the memory.
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqrecdirsend)
      deallocate (reqrecdirrec)
      deallocate (buffermpi1)
      deallocate (buffermpi2)
      deallocate (buffermpimu1)
      deallocate (buffermpimu2)
      deallocate (munew)
      deallocate (dipfield)
      deallocate (dipfieldbis)
      deallocate (dcfield1)
      deallocate (dcfield2)
      deallocate (h)
      if (dodiis) then
        deallocate (xdiis)
        deallocate (ediis)
        deallocate (bmat)
        deallocate (bloc)
        deallocate (cex)
      end if

      return
      end

      subroutine pc_dc_tmatxb_pme2(nrhs,dodiag,mu,efi)
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
      use shunt
      implicit none
      integer i,nrhs,iglob,kglob,kbis,iipole,kkpole,kkpoleloc
      integer ipoleloc,l
      real*8  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      logical dodiag
      integer j, ii, kkk, irhs, inl ,tpos
      real*8 dx, dy, dz, d, d2
      real*8 dukx, duky, dukz, pukx, puky, pukz
      real*8 duix, duiy, duiz, puix, puiy, puiz
      real*8  zero, one, f50
      real*8  erfc, cutoff2
      save    zero, one, f50
      data    zero/0.d0/, one/1.d0/, f50/50.d0/
      character*6 mode

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
          do irhs = 1, nrhs
            do j = 1, 3
              efi(j,irhs,i) = efi(j,irhs,i) +
     $           mu(j,irhs,i)/polarity(iipole)
            end do
          end do
        end do
      end if

      return
      end

      subroutine otf_dc_tmatxb_pme2(nrhs,dodiag,mu,efi)
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
      use shunt
      implicit none
      integer i,nrhs,iglob,kglob,kbis,iipole,kkpole,kkpoleloc
      integer ipoleloc,l
      real*8  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      logical dodiag
      integer j, ii, kkk, irhs, inl
      real*8  dx, dy, dz, d, d2, damp, expdamp, pgamma,
     $  scale3, scale5, pdi, pti 
      real*8 ralpha, alsq2, alsq2n, bfac, exp2a
      real*8 rr3, rr5, dukx, duky, dukz, pukx, puky, pukz, puir, pukr,
     $  duir, dukr
      real*8 duix, duiy, duiz, puix, puiy, puiz
      real*8 bn(0:3), fid(3), fip(3), fimd(3), fimp(3)
      real*8 fkd(3), fkp(3), fkmd(3), fkmp(3)
      real*8  zero, one, f50
      real*8, allocatable :: dscale(:)
      real*8  erfc, cutoff2
      save    zero, one, f50
      data    zero/0.d0/, one/1.d0/, f50/50.d0/
      character*6 mode
      external erfc
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
      allocate (dscale(n))
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
c     gather some parameters, then set up the damping factors.
c
      mode = 'EWALD'
      call switch (mode)
      cutoff2 = cut2
c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        l = grplst(iglob)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        if (i.eq.0) cycle
        pdi = pdamp(iipole)
        pti = thole(iipole)
        duix = mu(1,1,ipoleloc)
        duiy = mu(2,1,ipoleloc)
        duiz = mu(3,1,ipoleloc)
        puix = mu(1,2,ipoleloc)
        puiy = mu(2,2,ipoleloc)
        puiz = mu(3,2,ipoleloc)
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = u1scale
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = u2scale
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = u3scale
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = u4scale
        end do
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
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
            d  = sqrt(d2)
            dukx = mu(1,1,kkpoleloc)
            duky = mu(2,1,kkpoleloc)
            dukz = mu(3,1,kkpoleloc)
            pukx = mu(1,2,kkpoleloc)
            puky = mu(2,2,kkpoleloc)
            pukz = mu(3,2,kkpoleloc)
c
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0d0 * aewald**2
            alsq2n = 0.0d0
            if (aewald .gt. 0.0d0)
     &        alsq2n = 1.0d0 / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 2
              bfac = dble(j+j-1)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / d2
            end do
c
            scale3 = dscale(kglob)
            scale5 = dscale(kglob)
            damp = pdi*pdamp(kkpole)
            if (damp.ne.zero) then
              pgamma = min(pti,thole(kbis))
              damp = -pgamma*(d/damp)**3
              if (damp .gt. -f50) then
                expdamp = exp(damp)
                scale3 = scale3 * (one - expdamp)
                scale5 = scale5 * (one - expdamp*(one - damp))
              end if
            end if
c
c     compute the field.
c
            rr3 = (1.0d0-scale3) / (d*d2)
            rr5 = 3.0d0 * (1.0d0-scale5) / (d*d2*d2)
            duir = dx*duix + dy*duiy + dz*duiz
            dukr = dx*dukx + dy*duky + dz*dukz
            puir = dx*puix + dy*puiy + dz*puiz
            pukr = dx*pukx + dy*puky + dz*pukz
            fimd(1) = -bn(1)*dukx + bn(2)*dukr*dx
            fimd(2) = -bn(1)*duky + bn(2)*dukr*dy
            fimd(3) = -bn(1)*dukz + bn(2)*dukr*dz
            fkmd(1) = -bn(1)*duix + bn(2)*duir*dx
            fkmd(2) = -bn(1)*duiy + bn(2)*duir*dy
            fkmd(3) = -bn(1)*duiz + bn(2)*duir*dz
            fimp(1) = -bn(1)*pukx + bn(2)*pukr*dx
            fimp(2) = -bn(1)*puky + bn(2)*pukr*dy
            fimp(3) = -bn(1)*pukz + bn(2)*pukr*dz
            fkmp(1) = -bn(1)*puix + bn(2)*puir*dx
            fkmp(2) = -bn(1)*puiy + bn(2)*puir*dy
            fkmp(3) = -bn(1)*puiz + bn(2)*puir*dz
            fid(1) = -rr3*dukx + rr5*dukr*dx
            fid(2) = -rr3*duky + rr5*dukr*dy
            fid(3) = -rr3*dukz + rr5*dukr*dz
            fkd(1) = -rr3*duix + rr5*duir*dx
            fkd(2) = -rr3*duiy + rr5*duir*dy
            fkd(3) = -rr3*duiz + rr5*duir*dz
            fip(1) = -rr3*pukx + rr5*pukr*dx
            fip(2) = -rr3*puky + rr5*pukr*dy
            fip(3) = -rr3*pukz + rr5*pukr*dz
            fkp(1) = -rr3*puix + rr5*puir*dx
            fkp(2) = -rr3*puiy + rr5*puir*dy
            fkp(3) = -rr3*puiz + rr5*puir*dz
            efi(1,1,ipoleloc) = efi(1,1,ipoleloc) - fimd(1) + fid(1)
            efi(2,1,ipoleloc) = efi(2,1,ipoleloc) - fimd(2) + fid(2)
            efi(3,1,ipoleloc) = efi(3,1,ipoleloc) - fimd(3) + fid(3)
            efi(1,2,ipoleloc) = efi(1,2,ipoleloc) - fimp(1) + fip(1)
            efi(2,2,ipoleloc) = efi(2,2,ipoleloc) - fimp(2) + fip(2)
            efi(3,2,ipoleloc) = efi(3,2,ipoleloc) - fimp(3) + fip(3)
c
            efi(1,1,kkpoleloc) = efi(1,1,kkpoleloc) - fkmd(1) + fkd(1)
            efi(2,1,kkpoleloc) = efi(2,1,kkpoleloc) - fkmd(2) + fkd(2)
            efi(3,1,kkpoleloc) = efi(3,1,kkpoleloc) - fkmd(3) + fkd(3)
            efi(1,2,kkpoleloc) = efi(1,2,kkpoleloc) - fkmp(1) + fkp(1)
            efi(2,2,kkpoleloc) = efi(2,2,kkpoleloc) - fkmp(2) + fkp(2)
            efi(3,2,kkpoleloc) = efi(3,2,kkpoleloc) - fkmp(3) + fkp(3)
            end if
          end if
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = 1.0d0
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = 1.0d0
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = 1.0d0
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = 1.0d0
        end do
      end do
      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
        do i = 1, npoleloc
          iipole = poleglob(i)
          do irhs = 1, nrhs
            do j = 1, 3
              efi(j,irhs,i) = efi(j,irhs,i) +
     $           mu(j,irhs,i)/polarity(iipole)
            end do
          end do
        end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end
