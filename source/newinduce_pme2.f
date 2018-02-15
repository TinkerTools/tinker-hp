c      
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c      
      subroutine newinduce_pme2
      implicit none
c
c     evaluate induced dipole moments and the polarization energy
c     using either a (preconditioned) conjugate gradient algorithm or
c     Jacobi iterations coupled with DIIS extrapolation.
c
      include 'sizes.i'
      include 'atmlst.i'
      include 'ewald.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'math.i'
      include 'openmp.i'
      include 'pme.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
      include 'uprior.i'
      include 'mpif.h'
      integer i, j, k, nrhs, proc
c
c     MPI
c
      integer iglob, ierr, iipole
      integer reqrecdir(nproc*nproc)
      integer req(nproc*nproc),req2(nproc*nproc), tag
      integer status(MPI_STATUS_SIZE)
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
      external tmatxb_pme2
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
      if (associated(cphirec)) deallocate (cphirec)
      allocate (cphirec(10,max(npolerecloc,1)))
      cphirec = 0d0
      if (associated(fphirec)) deallocate(fphirec)
      allocate (fphirec(20,max(npolerecloc,1)))
      fphirec = 0d0
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
     $      MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir(tag),ierr)
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
     $     10*buflen2(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $     reqrecdir(tag),ierr)
        end if
      end do
c
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_WAIT(reqrecdir(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_WAIT(reqrecdir(tag),status,ierr)
          do j = 0, buflen1(proc+1)-1
            call amove(10,buffermpi1(1,bufbeg1(proc+1)+j),
     $        cphi(1,poleloc(buf1(bufbeg1(proc+1)+j))))
          end do
        end if
      end do
c
      call efld0_direct2(nrhs,ef)
c
      call commfield(nrhs,ef)
c
c     MPI : begin reception
c
      do i = 1, n_recep1
        tag = nproc*rank + p_recep1(i) + 1
        call MPI_IRECV(mu(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c       Add direct and reciprocal fields
c
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c$omp parallel do default(shared) private(i,j,iipole)
      do i = 1, npoleloc
        iipole = poleglob(i)
        do j = 1, 3
          ef(j,1,i)  = ef(j,1,i) - cphi(j+1,i) +
     $       term*rpole(j+1,iipole)
          ef(j,2,i)  = ef(j,2,i) - cphi(j+1,i) +
     $       term*rpole(j+1,iipole)
        end do
      end do
c$omp end parallel do
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
     $   MPI_REAL8,p_send1(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, n_recep1
        tag = nproc*rank + p_recep1(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, n_send1
        tag = nproc*p_send1(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     Begin reception of mu for PME
c
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),
     $      3*nrhs*buflen2(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $      req2(tag),ierr)
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
     $     3*nrhs*buflen1(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $     req2(tag),ierr)
        end if
      end do
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_WAIT(req2(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_WAIT(req2(tag),status,ierr)
          do j = 0, buflen2(proc+1)-1
            call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),
     $        murec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
          end do
        end if
      end do
c
c     now, call the proper solver.
c
      if (polalg.eq.1) then
        call inducepcg_pme2(tmatxb_pme2,nrhs,.true.,ef,mu,murec)
      else if (polalg.eq.2) then
        call inducejac_pme2(tmatxb_pme2,nrhs,.true.,ef,mu,murec)
      else
         if (rank.eq.0) write(iout,1000) 
         call fatal
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
c
      subroutine inducepcg_pme2(matvec,nrhs,precnd,ef,mu,murec)
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'ewald.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'iounit.i'
      include 'units.i'
      include 'openmp.i'
      include 'math.i'
      include 'timestat.i'
      include 'mpif.h'
c
c     solves the polarizable dipoles linear equations by preconditioned
c     conjugate gradient. A diagonal preconditioner is used when precnd
c     is true, otherwise the standard conjugate gradient algorithm is
c     recovered by setting the preconditioner to one.
c
      integer nrhs
      real*8  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      logical precnd
      real*8, allocatable :: res(:,:,:), h(:,:,:), pp(:,:,:), zr(:,:,:),
     $  diag(:)
      integer i, it, j, k
      real*8  ggold(2), ggnew(2), gnorm(2), gg(2), alphacg(2), ene(2)
      real*8  zero, pt5, one, resnrm, term
      save    zero, pt5, one
      data    zero/0.0d0/, pt5/0.50d0/, one/1.0d0/
      external matvec
      real*8, allocatable :: dipfield1(:,:),dipfield2(:,:)
      real*8, allocatable :: dipfield1bis(:,:),dipfield2bis(:,:)
      real*8 time0,time1,time2
c
c     MPI
c
      real*8, allocatable :: buffermpi11(:,:),buffermpi12(:,:)
      real*8, allocatable :: buffermpi21(:,:),buffermpi22(:,:)
      real*8, allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
      integer iglob, iipole, ierr, tag, proc
      integer req1, req2, req3, req4
      integer status(MPI_STATUS_SIZE)
      integer reqbis(nproc*nproc),req(nproc*nproc),reqend(nproc*nproc)
      integer reqendbis(nproc*nproc)
      integer,allocatable :: reqrecdir1(:),reqrecdir2(:)
c
 1000 format(' cgiter converged after ',I3,' iterations.',/,
     $       ' final energy        = ',2D14.7,/,
     $       ' final residual norm = ',2D14.7)
 1010 format(' energy and residual norm at iteration ',I3,':',4D12.2)
 1020 format(' Conjugate gradient solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Conjugate gradient solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
 1040 format(' Using a diagonal preconditioner.')
c
      allocate (buffermpi11(3,max(npoleloc,1)))
      buffermpi11 = 0d0
      allocate (buffermpi12(3,max(npoleloc,1)))
      buffermpi12 = 0d0
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0d0
      allocate (buffermpi21(3,max(npolerecloc,1)))
      buffermpi21 = 0d0
      allocate (buffermpi22(3,max(npolerecloc,1)))
      buffermpi22 = 0d0
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0d0
c
      allocate (dipfield1(3,max(1,npoleloc)))
      allocate (dipfield2(3,max(1,npoleloc)))
      dipfield1 = 0d0
      dipfield2 = 0d0
      allocate (dipfield1bis(3,max(1,npolerecloc)))
      allocate (dipfield2bis(3,max(1,npolerecloc)))
      dipfield1bis = 0d0
      dipfield2bis = 0d0
      allocate (reqrecdir1(nproc*nproc))
      allocate (reqrecdir2(nproc*nproc))
c
c     allocate some memory and setup the preconditioner:
c
      allocate (res(3,nrhs,max(1,npoleloc)))
      allocate (h(3,nrhs,max(1,npolebloc)))
      allocate (pp(3,nrhs,max(1,npolebloc)))
      allocate (zr(3,nrhs,max(1,npoleloc)))
      allocate (diag(npoleloc))
      if (precnd) then
        do i = 1, npoleloc
          iipole = poleglob(i)
          diag(i) = polarity(iipole)
        end do
        if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
        do i = 1, npoleloc
          diag(i) = one
        end do
      end if
c
c     initialize
c
      res = 0d0
      pp = 0d0
      zr = 0d0
      h = 0d0
c
c     now, compute the initial direction
c
      ggold = 0d0
c
      time0 = mpi_wtime()
      call tmatxbrecip2(mu,murec,nrhs,dipfield1,dipfield2,dipfield1bis,
     $  dipfield2bis)
      time1 = mpi_wtime()
      timerecdip = timerecdip + time1-time0
      call matvec(nrhs,.true.,mu,h)
      time2 = mpi_wtime()
      timerealdip = timerealdip + time2-time1
      call commfield(nrhs,h)
c
c     Begin the reception of the reciprocal fields
c
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
         tag = nproc*rank + proc + 1
         call MPI_IRECV(buffermpi11(1,bufbeg1(proc+1)),
     $     3*buflen1(proc+1),
     $     MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir1(tag),ierr)
         call MPI_IRECV(buffermpi12(1,bufbeg1(proc+1)),
     $     3*buflen1(proc+1),
     $     MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir2(tag),ierr)
        end if
      end do
c
c     MPI : begin reception
c
      do i = 1, n_recep1
        tag = nproc*rank + p_recep1(i) + 1
        call MPI_IRECV(pp(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
      end do
c
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          do j = 0, buflen2(proc+1)-1
            call amove(3,dipfield1bis(1,polerecloc(
     $        buf2(bufbeg2(proc+1)+j))),
     $        buffermpi21(1,bufbeg2(proc+1)+j))
            call amove(3,dipfield2bis(1,polerecloc(
     $        buf2(bufbeg2(proc+1)+j))),
     $        buffermpi22(1,bufbeg2(proc+1)+j))
          end do
        end if
      end do
c
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_ISEND(buffermpi21(1,bufbeg2(proc+1)),
     $     3*buflen2(proc+1),
     $     MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir1(tag),ierr)
          call MPI_ISEND(buffermpi22(1,bufbeg2(proc+1)),
     $     3*buflen2(proc+1),
     $     MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir2(tag),ierr)
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
     $     buflen2(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $     reqbis(tag),ierr)
        end if
      end do
c
c       Wait for the reciprocal fields
c
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_WAIT(reqrecdir1(tag),status,ierr)
          call MPI_WAIT(reqrecdir2(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_WAIT(reqrecdir1(tag),status,ierr)
          call MPI_WAIT(reqrecdir2(tag),status,ierr)
          do j = 0, buflen1(proc+1)-1
            call amove(3,buffermpi11(1,bufbeg1(proc+1)+j),
     $        dipfield1(1,poleloc(buf1(bufbeg1(proc+1)+j))))
            call amove(3,buffermpi12(1,bufbeg1(proc+1)+j),
     $        dipfield2(1,poleloc(buf1(bufbeg1(proc+1)+j))))
          end do
        end if
      end do
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c$omp parallel do default(shared) private(i,j)
      do i = 1, npoleloc
          do j = 1,3
            h(j,1,i) = h(j,1,i)+dipfield1(j,i)-
     $        term*mu(j,1,i)
            h(j,2,i) = h(j,2,i)+dipfield2(j,i)-
     $        term*mu(j,2,i)
          end do
        end do
c$omp end parallel do
c$omp parallel do default(shared) private(i,j,k)
c$omp+ reduction(+:ggold)
      do i = 1, npoleloc
        do k = 1, nrhs
          do j = 1, 3
            res(j,k,i) = ef(j,k,i) - h(j,k,i)
            zr(j,k,i)  = diag(i)*res(j,k,i)
            pp(j,k,i)  = zr(j,k,i)
            ggold(k)   = ggold(k)  + res(j,k,i)*zr(j,k,i)
          end do
        end do
      end do
c$omp end parallel do
      call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_REAL8,MPI_SUM,
     $     MPI_COMM_WORLD,req1,ierr)
c
c     MPI : begin sending
c
      do i = 1, n_send1
        tag = nproc*p_send1(i) + rank + 1
        call MPI_ISEND(pp,3*nrhs*npoleloc,
     $   MPI_REAL8,p_send1(i),tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
c     MPI : begin sending of mu for PME
c
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          do j = 0, buflen1(proc+1)-1
            call amove(3*nrhs,pp(1,1,poleloc(buf1(bufbeg1(proc+1)+j))),
     $        buffermpimu1(1,1,bufbeg1(proc+1)+j))
          end do
        end if
      end do
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_ISEND(buffermpimu1(1,1,bufbeg1(proc+1)),3*nrhs*
     $     buflen1(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $     reqbis(tag),ierr)
        end if
      end do
      do i = 1, n_recep1
        tag = nproc*rank + p_recep1(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, n_send1
        tag = nproc*p_send1(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_WAIT(reqbis(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_WAIT(reqbis(tag),status,ierr)
          do j = 0, buflen2(proc+1)-1
            call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),
     $        murec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
          end do
        end if
      end do
      call MPI_WAIT(req1,status,ierr)
c
c     now, start the main loop:
c
      do it = 1, politer
        do k = 1, nrhs
          gg(k) = zero
        end do
c
        call tmatxbrecip2(pp,murec,nrhs,dipfield1,dipfield2,
     $   dipfield1bis,dipfield2bis)
        call matvec(nrhs,.true.,pp,h)
        call commfield(nrhs,h)
c
c     Begin the reception of the reciprocal fields
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
           tag = nproc*rank + proc + 1
           call MPI_IRECV(buffermpi11(1,bufbeg1(proc+1)),
     $       3*buflen1(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir1(tag),ierr)
           call MPI_IRECV(buffermpi12(1,bufbeg1(proc+1)),
     $       3*buflen1(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir2(tag),ierr)
          end if
        end do
c
c     MPI : begin reception
c
        do i = 1, n_recep1
          tag = nproc*rank + p_recep1(i) + 1
        call MPI_IRECV(pp(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
        end do
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            do j = 0, buflen2(proc+1)-1
              call amove(3,dipfield1bis(1,polerecloc(
     $          buf2(bufbeg2(proc+1)+j))),
     $          buffermpi21(1,bufbeg2(proc+1)+j))
              call amove(3,dipfield2bis(1,polerecloc(
     $          buf2(bufbeg2(proc+1)+j))),
     $          buffermpi22(1,bufbeg2(proc+1)+j))
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpi21(1,bufbeg2(proc+1)),
     $       3*buflen2(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir1(tag),ierr)
            call MPI_ISEND(buffermpi22(1,bufbeg2(proc+1)),
     $       3*buflen2(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir2(tag),ierr)
          end if
        end do
c
c     Begin reception of mu for PME
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),
     $      3*nrhs*buflen2(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $      reqbis(tag),ierr)
          end if
        end do
c
c       Wait for the reciprocal fields
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_WAIT(reqrecdir1(tag),status,ierr)
            call MPI_WAIT(reqrecdir2(tag),status,ierr)
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_WAIT(reqrecdir1(tag),status,ierr)
            call MPI_WAIT(reqrecdir2(tag),status,ierr)
            do j = 0, buflen1(proc+1)-1
              call amove(3,buffermpi11(1,bufbeg1(proc+1)+j),
     $          dipfield1(1,poleloc(buf1(bufbeg1(proc+1)+j))))
              call amove(3,buffermpi12(1,bufbeg1(proc+1)+j),
     $          dipfield2(1,poleloc(buf1(bufbeg1(proc+1)+j))))
            end do
          end if
        end do
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c$omp parallel do default(shared) private(i,j)
        do i = 1, npoleloc
            do j = 1,3
              h(j,1,i) = h(j,1,i)+dipfield1(j,i)-
     $          term*pp(j,1,i)
              h(j,2,i) = h(j,2,i)+dipfield2(j,i)-
     $          term*pp(j,2,i)
            end do
          end do
c$omp end parallel do
c
c$omp parallel do default(shared) private(i,j,k)
c$omp+ reduction(+:gg)
        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              gg(k)  = gg(k)  + pp(j,k,i)*h(j,k,i)
            end do
          end do
        end do
c$omp end parallel do
        call MPI_IALLREDUCE(MPI_IN_PLACE,gg(1),nrhs,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,req2,ierr)
        call MPI_WAIT(req2,status,ierr)
        do k = 1, nrhs
          if (gg(k).eq.zero) return
          alphacg(k)  = ggold(k)/gg(k)
          ggnew(k)  = zero
          ene(k)    = zero
        end do
c$omp parallel do default(shared) private(i,j,k)
c$omp+ reduction(+:ggnew) reduction(-:ene)
        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              mu(j,k,i)  = mu(j,k,i) + alphacg(k)*pp(j,k,i)
              res(j,k,i) = res(j,k,i) - alphacg(k)*h(j,k,i)
              zr(j,k,i)  = diag(i)*res(j,k,i)
              ggnew(k)   = ggnew(k) + res(j,k,i)*zr(j,k,i)
              ene(k)     = ene(k) - pt5*mu(j,k,i)*
     $           (res(j,k,i)+ef(j,k,i))
            end do
          end do
        end do
c$omp end parallel do
        call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_REAL8,
     $    MPI_SUM,MPI_COMM_WORLD,req3,ierr)
        call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_REAL8,MPI_SUM,
     $    MPI_COMM_WORLD,req4,ierr)
        call MPI_WAIT(req3,status,ierr)
        call MPI_WAIT(req4,status,ierr)
        resnrm = zero
        do k = 1, nrhs
          gnorm(k) = sqrt(ggnew(k)/dble(3*npolar))
          resnrm   = max(resnrm,gnorm(k))
        end do
        if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     $    it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)
c$omp parallel do default(shared) private(i,j,k)
        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
              pp(j,k,i) = zr(j,k,i) + ggnew(k)/ggold(k)*pp(j,k,i)
            end do
          end do
        end do
c$omp end parallel do
c
c     MPI : begin sending
c
        do i = 1, n_send1
          tag = nproc*p_send1(i) + rank + 1
          call MPI_ISEND(pp,3*nrhs*
     $     npoleloc,MPI_REAL8,p_send1(i),tag,MPI_COMM_WORLD,
     $     req(tag),ierr)
        end do
c
c
c     MPI : begin sending of mu for PME
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            do j = 0, buflen1(proc+1)-1
              call amove(3*nrhs,pp(1,1,poleloc(
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
     $       buflen1(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $       reqbis(tag),ierr)
          end if
        end do
        do i = 1, n_recep1
          tag = nproc*rank + p_recep1(i) + 1
          call MPI_WAIT(req(tag),status,ierr)
        end do
        do i = 1, n_send1
          tag = nproc*p_send1(i) + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
        end do
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_WAIT(reqbis(tag),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_WAIT(reqbis(tag),status,ierr)
            do j = 0, buflen2(proc+1)-1
              call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),
     $          murec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
            end do
          end if
        end do
c
        do k = 1, nrhs
          ggold(k) = ggnew(k)
        end do
        if (resnrm.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     $      (ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
          goto 10
        end if
      end do
 10   continue
c
c     MPI : begin reception
c
      do i = 1, n_recep1
        tag = nproc*rank + p_recep1(i) + 1
        call MPI_IRECV(mu(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),tag,
     $   MPI_COMM_WORLD,reqend(tag),ierr)
      end do
c
c     Begin reception of mu for PME
c
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_IRECV(buffermpimu2(1,1,bufbeg2(proc+1)),3*nrhs*
     $     buflen2(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $     reqendbis(tag),ierr)
        end if
      end do
c
c     MPI : begin sending
c
      do i = 1, n_send1
        tag = nproc*p_send1(i) + rank + 1
        call MPI_ISEND(mu,3*nrhs*
     $   npoleloc,MPI_REAL8,p_send1(i),tag,MPI_COMM_WORLD,
     $   reqend(tag),ierr)
      end do
c
c
c     MPI : begin sending of mu for PME
c
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
          call MPI_ISEND(buffermpimu1(1,1,bufbeg1(proc+1)),3*nrhs*
     $     buflen1(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $     reqendbis(tag),ierr)
        end if
      end do
      do i = 1, n_send1
        tag = nproc*p_send1(i) + rank + 1
        call MPI_WAIT(reqend(tag),status,ierr)
      end do
      do i = 1, n_recep1
        tag = nproc*rank + p_recep1(i) + 1
        call MPI_WAIT(reqend(tag),status,ierr)
      end do
c
      do i = 1, nrecdir_send1
        proc = precdir_send1(i)
        if (proc.ne.rank) then
          tag = nproc*proc + rank + 1
          call MPI_WAIT(reqendbis(tag),status,ierr)
        end if
      end do
      do i = 1, nrecdir_recep1
        proc = precdir_recep1(i)
        if (proc.ne.rank) then
          tag = nproc*rank + proc + 1
          call MPI_WAIT(reqendbis(tag),status,ierr)
          do j = 0, buflen2(proc+1)-1
            call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),
     $        murec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
          end do
        end if
      end do
c
c     finalize and debug printing:
c
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
          iglob = ipole(iipole)
          write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
        end do
      end if
      deallocate (reqrecdir1)
      deallocate (reqrecdir2)
      deallocate (buffermpi11)
      deallocate (buffermpi12)
      deallocate (buffermpimu1)
      deallocate (buffermpimu2)
      deallocate (buffermpi21)
      deallocate (buffermpi22)
      deallocate (dipfield1)
      deallocate (dipfield2)
      deallocate (res)
      deallocate (h)
      deallocate (pp)
      deallocate (zr)
      deallocate (diag)
      return
      end
c
      subroutine inducejac_pme2(matvec,nrhs,dodiis,ef,mu,murec)
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'ewald.i'
      include 'polar.i'
      include 'mpole.i'
      include 'math.i'
      include 'polpot.i'
      include 'iounit.i'
      include 'openmp.i'
      include 'timestat.i'
      include 'mpif.h'
c
c     solves the polarizable dipoles linear equations by jacobi
c     iterations coupled to the DIIS extrapolation.
c
      integer nrhs, info, proc
      real*8  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      real*8, allocatable :: munew(:,:,:), h(:,:,:)
      real*8, allocatable :: xdiis(:,:), ediis(:,:), bmat(:,:),
     $   bloc(:,:), cex(:)
      integer i, j, k, it, ind, ndismx, nmat, lenb, tag
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real*8  zero, one, rnorm(2), rr, xx(1)
      real*8 term
      real*8 time0,time1,time2,time3
      save    zero, one, xx
      data    zero/0.0d0/, one/1.0d0/, xx/0.0d0/
      external matvec
      real*8, allocatable :: dipfield1(:,:),dipfield2(:,:)
      real*8, allocatable :: dipfield1bis(:,:),dipfield2bis(:,:)
c
c     MPI
c
      real*8, allocatable :: buffermpi11(:,:),buffermpi12(:,:)
      real*8, allocatable :: buffermpi21(:,:),buffermpi22(:,:)
      real*8, allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
      integer iglob, iipole, ierr
      integer reqnorm,reqdiis(2*ndismx+1)
      integer status(MPI_STATUS_SIZE)
      integer req(nproc*nproc),req2(nproc*nproc)
      integer,allocatable :: reqrecdir1(:),reqrecdir2(:)
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
      allocate (buffermpi11(3,max(npoleloc,1)))
      buffermpi11 = 0d0
      allocate (buffermpi12(3,max(1,npoleloc)))
      buffermpi12 = 0d0
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = 0d0
      allocate (buffermpi21(3,max(npolerecloc,1)))
      buffermpi21 = 0d0
      allocate (buffermpi22(3,max(npolerecloc,1)))
      buffermpi22 = 0d0
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = 0d0
c
      allocate (munew(3,nrhs,npolebloc))
      munew = 0d0
      allocate (dipfield1(3,max(1,npoleloc)))
      allocate (dipfield2(3,max(1,npoleloc)))
      dipfield1 = 0d0
      dipfield2 = 0d0
      allocate (dipfield1bis(3,max(1,npolerecloc)))
      allocate (dipfield2bis(3,max(1,npolerecloc)))
      dipfield1bis = 0d0
      dipfield2bis = 0d0
      allocate (reqrecdir1(nproc*nproc))
      allocate (reqrecdir2(nproc*nproc))
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
       call tmatxbrecip2(mu,murec,nrhs,dipfield1,dipfield2,dipfield1bis,
     $   dipfield2bis)
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
           call MPI_IRECV(buffermpi11(1,bufbeg1(proc+1)),
     $       3*buflen1(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir1(tag),ierr)
           call MPI_IRECV(buffermpi12(1,bufbeg1(proc+1)),
     $       3*buflen1(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir2(tag),ierr)
          end if
        end do
c
c     MPI : begin reception
c
        do i = 1, n_recep1
          tag = nproc*rank + p_recep1(i) + 1
        call MPI_IRECV(mu(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),tag,
     $   MPI_COMM_WORLD,req(tag),ierr)
        end do
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            do j = 0, buflen2(proc+1)-1
              call amove(3,dipfield1bis(1,polerecloc(
     $          buf2(bufbeg2(proc+1)+j))),
     $          buffermpi21(1,bufbeg2(proc+1)+j))
              call amove(3,dipfield2bis(1,polerecloc(
     $          buf2(bufbeg2(proc+1)+j))),
     $          buffermpi22(1,bufbeg2(proc+1)+j))
            end do
          end if
        end do
c
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_ISEND(buffermpi21(1,bufbeg2(proc+1)),
     $       3*buflen2(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir1(tag),ierr)
            call MPI_ISEND(buffermpi22(1,bufbeg2(proc+1)),
     $       3*buflen2(proc+1),
     $       MPI_REAL8,proc,tag,MPI_COMM_WORLD,reqrecdir2(tag),ierr)
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
     $       buflen2(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $       req2(tag),ierr)
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
            tag = nproc*proc + rank + 1
            call MPI_WAIT(reqrecdir1(tag),status,ierr)
            call MPI_WAIT(reqrecdir2(tag),status,ierr)
          end if
        end do
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_WAIT(reqrecdir1(tag),status,ierr)
            call MPI_WAIT(reqrecdir2(tag),status,ierr)
            do j = 0, buflen1(proc+1)-1
              call amove(3,buffermpi11(1,bufbeg1(proc+1)+j),
     $          dipfield1(1,poleloc(buf1(bufbeg1(proc+1)+j))))
              call amove(3,buffermpi12(1,bufbeg1(proc+1)+j),
     $          dipfield2(1,poleloc(buf1(bufbeg1(proc+1)+j))))
            end do
          end if
        end do
        time2 = mpi_wtime()
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c$omp parallel do default(shared) private(i,j)
        do i =  1, npoleloc
          do j = 1,3
            h(j,1,i) = h(j,1,i)+dipfield1(j,i)-
     $        term*mu(j,1,i)
            h(j,2,i) = h(j,2,i)+dipfield2(j,i)-
     $        term*mu(j,2,i)
          end do
        end do
c$omp end parallel do
c
c$omp parallel do default(shared) private(i,iipole,j,k)
c$omp+ reduction(+:rnorm)
        do i = 1, npoleloc
          iipole = poleglob(i)
          do k = 1, nrhs
            do j = 1, 3
              munew(j,k,i) = polarity(iipole)*(ef(j,k,i) - h(j,k,i))
              rnorm(k) = rnorm(k) + (munew(j,k,i) - mu(j,k,i))**2
            end do
          end do
        end do
c$omp end parallel do
c
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_REAL8,
     $    MPI_SUM,MPI_COMM_WORLD,reqnorm,ierr)
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
          call diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $      bmat,nmat,reqdiis,MPI_COMM_WORLD)
c
          do i = 1, 2*nmat-3
            call MPI_WAIT(reqdiis(i),status,ierr)
          end do
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
     $     MPI_REAL8,p_send1(i),tag,MPI_COMM_WORLD,
     $     req(tag),ierr)
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
     $       buflen1(proc+1),MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $       req2(tag),ierr)
          end if
        end do
       do i = 1, n_recep1
          tag = nproc*rank + p_recep1(i) + 1
          call MPI_WAIT(req(tag),status,ierr)
        end do
        do i = 1, n_send1
          tag = nproc*p_send1(i) + rank + 1
          call MPI_WAIT(req(tag),status,ierr)
        end do
c
        do i = 1, nrecdir_send1
          proc = precdir_send1(i)
          if (proc.ne.rank) then
            tag = nproc*proc + rank + 1
            call MPI_WAIT(req2(tag),status,ierr)
          end if
        end do
        do i = 1, nrecdir_recep1
          proc = precdir_recep1(i)
          if (proc.ne.rank) then
            tag = nproc*rank + proc + 1
            call MPI_WAIT(req2(tag),status,ierr)
            do j = 0, buflen2(proc+1)-1
              call amove(3*nrhs,buffermpimu2(1,1,bufbeg2(proc+1)+j),
     $          murec(1,1,polerecloc(buf2(bufbeg2(proc+1)+j))))
            end do
          end if
        end do
        call amove(3*nrhs*npoleloc,munew,mu)
c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/dble(3*npolar))
          rr = max(rnorm(k),rr)
        end do
        if (rr.lt.poleps) then
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
      deallocate (reqrecdir1)
      deallocate (reqrecdir2)
      deallocate (buffermpi11)
      deallocate (buffermpi12)
      deallocate (buffermpimu1)
      deallocate (buffermpimu2)
      deallocate (buffermpi21)
      deallocate (buffermpi22)
      deallocate (munew)
      deallocate (dipfield1)
      deallocate (dipfield2)
      deallocate (dipfield1bis)
      deallocate (dipfield2bis)
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
c
      subroutine efld0_direct2(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'couple.i'
      include 'ewald.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'openmp.i'
      integer i,iglob,kglob,kbis,nrhs,inl
      real*8  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole, ipoleloc
      real*8  dx, dy, dz, d, d2, d3, d5, pdi, pti, ck,
     $  dkx, dky, dkz, dkr, qkxx, qkyy, qkzz, qkxy, qkxz, qkyz,
     $  qkx, qky, qkz, qkr
      real*8 dix, diy, diz, dir, qixx, qiyy, qizz, qixy, qixz, qiyz,
     $  qix, qiy, qiz, qir, ci
      real*8 damp, pgamma, expdamp, scale3, scale5, scale7
      real*8 drr3,drr5,drr7
      real*8 prr3,prr5,prr7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 zero, pt6, one, f50
      real*8 ralpha, alsq2, alsq2n, bfac, exp2a
      real*8 bn(0:3), fim(3), fid(3), fip(3)
      real*8 fkm(3), fkd(3), fkp(3)
      real*8 erfc, cutoff2
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      save   zero, pt6, one, f50
      data   zero/0.d0/, pt6/0.6d0/, one/1.d0/, f50/50.d0/
      character*6 mode
c
 1000 format(' Warning, system moved too much since last neighbor list
     $   update, try lowering nlupdate')
c
      mode = 'EWALD'
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      dscale = 1.0d0
      pscale = 1.0d0
c
c
c$omp parallel do default(shared) private(i,ii,iipole,iglob,pdi,pti,ci)
c$omp+ private(diy,diz,qixx,qixy,qixz,qiyy,qiyz,qizz,j,kkk,k,kglob,kbis)
c$omp+ private(dx,dy,dz,d2,d,ralpha,bn,alsq2,alsq2n,exp2a,bfac,d3,dix)
c$omp+ private(d5,ck,dkx,dky,dkz,qkxx,qkxy,qkxz,qkyy,qkyz,qkzz,damp)
c$omp+ private(pgamma,expdamp,scale3,scale5,scale7,dsc3,dsc5,dsc7)
c$omp+ private(psc3,psc5,psc7,drr3,drr5,drr7,prr3,prr5,prr7)
c$omp+ private(dkr,qkx,qky,qkz,qkr,fim,fid,fip,kkpole)
c$omp+ private(dir,qix,qiy,qiz,fkm,fkd,fkp,inl,ipoleloc)
c$omp+ reduction(+:ef)
c$omp+ firstprivate(pscale,dscale)
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        inl = polelocnl(iipole)
        if (inl.eq.0) cycle
        pdi = pdamp(iipole)
        pti = thole(iipole)
        ci   = rpole(1,iipole)
        dix  = rpole(2,iipole)
        diy  = rpole(3,iipole)
        diz  = rpole(4,iipole)
        qixx = rpole(5,iipole)
        qixy = rpole(6,iipole)
        qixz = rpole(7,iipole)
        qiyy = rpole(9,iipole)
        qiyz = rpole(10,iipole)
        qizz = rpole(13,iipole)
        do j = 1, n12(iglob)
           pscale(i12(j,iglob)) = p2scale
        end do
        do j = 1, n13(iglob)
           pscale(i13(j,iglob)) = p3scale
        end do
        do j = 1, n14(iglob)
           pscale(i14(j,iglob)) = p4scale
           do k = 1, np11(iglob)
              if (i14(j,iglob) .eq. ip11(k,iglob))
     &           pscale(i14(j,iglob)) = p4scale * p41scale
           end do
        end do
        do j = 1, n15(iglob)
           pscale(i15(j,iglob)) = p5scale
        end do
        do j = 1, np11(iglob)
           dscale(ip11(j,iglob)) = d1scale
        end do
        do j = 1, np12(iglob)
           dscale(ip12(j,iglob)) = d2scale
        end do
        do j = 1, np13(iglob)
           dscale(ip13(j,iglob)) = d3scale
        end do
        do j = 1, np14(iglob)
           dscale(ip14(j,iglob)) = d4scale
        end do
c
        do kkk = 1, nelst(inl)
          kglob = elst(kkk,inl)
          kkpole = pollist(kglob)
          kbis = poleloc(kglob)
          if (kbis.eq.0) then
            write(iout,1000)
            cycle
          end if
          dx = x(kglob) - x(iglob)
          dy = y(kglob) - y(iglob)
          dz = z(kglob) - z(iglob)
          call image(dx,dy,dz)
          d2 = dx*dx + dy*dy + dz*dz
          if (d2.le.cutoff2) then
            d  = sqrt(d2)
c
c     calculate the error function damping terms
c
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0d0 * aewald**2
            alsq2n = 0.0d0
            if (aewald .gt. 0.0d0)
     &        alsq2n = 1.0d0 / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 3
              bfac = dble(j+j-1)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / d2
            end do
c
            d3   = d*d2
            d5   = d3*d2
            ck   = rpole(1,kkpole)
            dkx  = rpole(2,kkpole)
            dky  = rpole(3,kkpole)
            dkz  = rpole(4,kkpole)
            qkxx = rpole(5,kkpole)
            qkxy = rpole(6,kkpole)
            qkxz = rpole(7,kkpole)
            qkyy = rpole(9,kkpole)
            qkyz = rpole(10,kkpole)
            qkzz = rpole(13,kkpole)
            damp = pdi*pdamp(kkpole)
            scale3 = one
            scale5 = one
            scale7 = one
            if (damp.ne.zero) then
              pgamma = min(pti,thole(kkpole))
              damp = -pgamma*(d/damp)**3
              if (damp.gt.-F50) then
                expdamp = exp(damp)
                scale3 = one - expdamp
                scale5 = one - expdamp*(one-damp)
                scale7 = one - expdamp
     &                      *(one-damp + pt6*damp**2)
              end if
            end if
            dsc3 = scale3 * dscale(kglob)
            dsc5 = scale5 * dscale(kglob)
            dsc7 = scale7 * dscale(kglob)
            psc3 = scale3 * pscale(kglob)
            psc5 = scale5 * pscale(kglob)
            psc7 = scale7 * pscale(kglob)
            drr3 = (1.0d0-dsc3) / (d*d2)
            drr5 = 3.0d0 * (1.0d0-dsc5) / (d*d2*d2)
            drr7 = 15.0d0 * (1.0d0-dsc7) / (d*d2*d2*d2)
            prr3 = (1.0d0-psc3) / (d*d2)
            prr5 = 3.0d0 * (1.0d0-psc5) / (d*d2*d2)
            prr7 = 15.0d0 * (1.0d0-psc7) / (d*d2*d2*d2)
c
c     compute some intermediate quantities
c
            dir = dix*dx + diy*dy + diz*dz
            qix = qixx*dx + qixy*dy + qixz*dz
            qiy = qixy*dx + qiyy*dy + qiyz*dz
            qiz = qixz*dx + qiyz*dy + qizz*dz
            qir = qix*dx + qiy*dy + qiz*dz
            dkr = dkx*dx + dky*dy + dkz*dz
            qkx = qkxx*dx + qkxy*dy + qkxz*dz
            qky = qkxy*dx + qkyy*dy + qkyz*dz
            qkz = qkxz*dx + qkyz*dy + qkzz*dz
            qkr = qkx*dx + qky*dy + qkz*dz
c
            fim(1) = -dx*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                   - bn(1)*dkx + 2.0d0*bn(2)*qkx
            fim(2) = -dy*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                   - bn(1)*dky + 2.0d0*bn(2)*qky
            fim(3) = -dz*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                   - bn(1)*dkz + 2.0d0*bn(2)*qkz
            fkm(1) = dx*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                  - bn(1)*dix - 2.0d0*bn(2)*qix
            fkm(2) = dy*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                  - bn(1)*diy - 2.0d0*bn(2)*qiy
            fkm(3) = dz*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                  - bn(1)*diz - 2.0d0*bn(2)*qiz
            fid(1) = -dx*(drr3*ck-drr5*dkr+drr7*qkr)
     &                   - drr3*dkx + 2.0d0*drr5*qkx
            fid(2) = -dy*(drr3*ck-drr5*dkr+drr7*qkr)
     &                  - drr3*dky + 2.0d0*drr5*qky
            fid(3) = -dz*(drr3*ck-drr5*dkr+drr7*qkr)
     &                   - drr3*dkz + 2.0d0*drr5*qkz
            fkd(1) = dx*(drr3*ci+drr5*dir+drr7*qir)
     &                  - drr3*dix - 2.0d0*drr5*qix
            fkd(2) = dy*(drr3*ci+drr5*dir+drr7*qir)
     &                  - drr3*diy - 2.0d0*drr5*qiy
            fkd(3) = dz*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diz - 2.0d0*drr5*qiz
            fip(1) = -dx*(prr3*ck-prr5*dkr+prr7*qkr)
     &                   - prr3*dkx + 2.0d0*prr5*qkx
            fip(2) = -dy*(prr3*ck-prr5*dkr+prr7*qkr)
     &                   - prr3*dky + 2.0d0*prr5*qky
            fip(3) = -dz*(prr3*ck-prr5*dkr+prr7*qkr)
     &                   - prr3*dkz + 2.0d0*prr5*qkz
            fkp(1) = dx*(prr3*ci+prr5*dir+prr7*qir)
     &                  - prr3*dix - 2.0d0*prr5*qix
            fkp(2) = dy*(prr3*ci+prr5*dir+prr7*qir)
     &                  - prr3*diy - 2.0d0*prr5*qiy
            fkp(3) = dz*(prr3*ci+prr5*dir+prr7*qir)
     &                  - prr3*diz - 2.0d0*prr5*qiz
            ef(1,1,ipoleloc) = ef(1,1,ipoleloc) + fim(1) - fid(1)
            ef(2,1,ipoleloc) = ef(2,1,ipoleloc) + fim(2) - fid(2)
            ef(3,1,ipoleloc) = ef(3,1,ipoleloc) + fim(3) - fid(3)
            ef(1,2,ipoleloc) = ef(1,2,ipoleloc) + fim(1) - fip(1)
            ef(2,2,ipoleloc) = ef(2,2,ipoleloc) + fim(2) - fip(2)
            ef(3,2,ipoleloc) = ef(3,2,ipoleloc) + fim(3) - fip(3)
c
            ef(1,1,kbis) = ef(1,1,kbis) + fkm(1) - fkd(1)
            ef(2,1,kbis) = ef(2,1,kbis) + fkm(2) - fkd(2)
            ef(3,1,kbis) = ef(3,1,kbis) + fkm(3) - fkd(3)
            ef(1,2,kbis) = ef(1,2,kbis) + fkm(1) - fkp(1)
            ef(2,2,kbis) = ef(2,2,kbis) + fkm(2) - fkp(2)
            ef(3,2,kbis) = ef(3,2,kbis) + fkm(3) - fkp(3)
          end if
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
        do j = 1, n12(iglob)
           pscale(i12(j,iglob)) = 1.0d0
        end do
        do j = 1, n13(iglob)
           pscale(i13(j,iglob)) = 1.0d0
        end do
        do j = 1, n14(iglob)
           pscale(i14(j,iglob)) = 1.0d0
        end do
        do j = 1, n15(iglob)
           pscale(i15(j,iglob)) = 1.0d0
        end do
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
c$omp end parallel do
c
      deallocate (dscale)
      deallocate (pscale)
c
      return
      end
c
      subroutine tmatxb_pme2(nrhs,dodiag,mu,efi)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'openmp.i'
      integer i,iipole,nrhs,iglob,kglob,kkpole,kkpoleloc
      integer ipoleloc
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
      allocate (dscale(n))
      do i = 1, n
         dscale(i) = 1.0d0
      end do
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
c
c$omp parallel do default(shared)
c$omp+ private(ii,i,iipole,iglob,pdi,pti,kkk,kglob,dx,dy,dz,d2,d)
c$omp+ private(dukx,duky,dukz,pukx,puky,pukz,ralpha,bn,exp2a,damp)
c$omp+ private(expdamp,scale3,scale5,rr3,rr5,dukr,pukr,fimd,fimp,fid)
c$omp+ private(puir,duix,duiy,duiz,puix,puiy,puiz,duir,fkmd,fkmp,fkd)
c$omp+ private(bfac,alsq2,alsq2n,fip,j,kkpole,kkpoleloc,inl,fkp,pgamma)
c$omp+ private(ipoleloc)
c$omp+ schedule(dynamic)
c$omp+ firstprivate(dscale)
c$omp+ reduction(+:efi)
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        inl = polelocnl(iipole)
        if (inl.eq.0) cycle
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
        do kkk = 1,nelst(inl)
          kglob = elst(kkk,inl)
          kkpole = pollist(kglob)
          kkpoleloc = poleloc(kkpole)
          if (kkpoleloc.eq.0) then
            goto 10
          end if
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
              pgamma = min(pti,thole(kkpole))
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
            do j = 1,3
              efi(j,1,ipoleloc) = efi(j,1,ipoleloc) - fimd(j) + fid(j)
              efi(j,2,ipoleloc) = efi(j,2,ipoleloc) - fimp(j) + fip(j)
              efi(j,1,kkpoleloc) = efi(j,1,kkpoleloc) - fkmd(j) + fkd(j)
              efi(j,2,kkpoleloc) = efi(j,2,kkpoleloc) - fkmp(j) + fkp(j)
            end do
          end if
 10       continue
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
c$omp end parallel do
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
c
      return
      end
c
      subroutine efld0_recip2(cphi)
c
c     Compute the reciprocal space contribution to the electric field due to the permanent 
c     multipoles
c
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'bound.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'fft.i'
      include 'math.i'
      include 'mpole.i'
      include 'pme.i'
      include 'openmp.i'
      include 'mpif.h'
      integer ierr,iglob,iloc
      integer status(MPI_STATUS_SIZE),tag
      integer i,j,k,iipole
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer proc
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 cmp(10),fmp(10), cphi(10,*)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: req(:),reqbcast(:),reqinit(:)
c
      allocate (req(nproc*nproc))
      allocate (reqbcast(nproc*nproc))
      allocate (reqinit(nproc))
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     zero out the PME charge grid
c
      qgridin_2d = 0d0
c
c     fill the pme grid, loop over the multipoles sites
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         cmp(1) = rpole(1,iipole)
         cmp(2) = rpole(2,iipole)
         cmp(3) = rpole(3,iipole)
         cmp(4) = rpole(4,iipole)
         cmp(5) = rpole(5,iipole)
         cmp(6) = rpole(9,iipole)
         cmp(7) = rpole(13,iipole)
         cmp(8) = 2.0d0 * rpole(6,iipole)
         cmp(9) = 2.0d0 * rpole(7,iipole)
         cmp(10) = 2.0d0 * rpole(10,iipole)
c
c     compute B-spline coefficients
c
        call bspline_fill_site(iglob,i)
c
c     convert Cartesian multipoles to fractional coordinates
c
        call cmp_to_fmp_site(cmp,fmp)
c
c     assign PME grid
c
        call grid_mpole_site(iglob,i,fmp)
      end do
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nproc*rank + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nproc*prec_send(i) + rank + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nrec_recep
        tag = nproc*rank + prec_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nrec_send
        tag = nproc*prec_send(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
        call aadd(2*n1mpimax*n2mpimax*n3mpimax,
     $   qgridin_2d(1,1,1,1,1),qgridmpi(1,1,1,1,i),
     $   qgridin_2d(1,1,1,1,1))
      end do
c
c     Perform 3-D FFT forward transform
c
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     make the scalar summation over reciprocal lattice
c
      if ((istart2(rank+1).eq.1).and.(jstart2(rank+1).eq.1).and.
     $   (kstart2(rank+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0d0
      end if
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c$omp parallel do default(shared) private(i,k3,j,k2,k1,m1,m2,m3)
c$omp+ private(r1,r2,r3,h1,h2,h3,hsq,term,expterm,denom)
      do k3 = kstart2(rank+1),kend2(rank+1)
        do k2 = jstart2(rank+1),jend2(rank+1)
          do k1 = istart2(rank+1),iend2(rank+1)
            m1 = k1 - 1
            m2 = k2 - 1
            m3 = k3 - 1
            if (k1 .gt. nf1)  m1 = m1 - nfft1
            if (k2 .gt. nf2)  m2 = m2 - nfft2
            if (k3 .gt. nf3)  m3 = m3 - nfft3
            if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
            r1 = dble(m1)
            r2 = dble(m2)
            r3 = dble(m3)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = 0.0d0
            if ((term .gt. -50.0d0)) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               end if
            qfac_2d(k1-istart2(rank+1)+1,k2-jstart2(rank+1)+1,k3-
     $        kstart2(rank+1)+1) = expterm
            end if
 10         continue
          end do
        end do
      end do
c$omp end parallel do
c
c     account for the zeroth grid point for a finite system
c
      if ((istart2(rank+1).eq.1).and.(jstart2(rank+1).eq.1).and.
     $   (kstart2(rank+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5d0 * pi / xbox
           qfac_2d(1,1,1) = expterm
        end if
      end if
c
c     complete the transformation of the charge grid
c
      do k = 1, ksize2(rank+1)
         do j = 1, jsize2(rank+1)
           do i = 1, isize2(rank+1)
             term = qfac_2d(i,j,k)
             qgridout_2d(1,i,j,k) = term*qgridout_2d(1,i,j,k)
             qgridout_2d(2,i,j,k) = term*qgridout_2d(2,i,j,k)
           end do
         end do
      end do
c
c     perform 3-D FFT backward transform
c
      call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     MPI : Begin reception
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nproc*rank + prec_send(i) + 1
        call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,MPI_COMM_WORLD,reqbcast(tag),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nproc*prec_recep(i) + rank + 1
        call MPI_ISEND(qgridin_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,MPI_COMM_WORLD,reqbcast(tag),ierr)
      end do
c
      do i = 1, nrec_send
        tag = nproc*rank + prec_send(i) + 1
        call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
      do i = 1, nrec_recep
        tag = nproc*prec_recep(i) + rank + 1
        call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
c
c     get field
c
c$omp parallel do default(shared) private(i,iipole,iglob,iloc)
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        iloc  = poleloc(iipole)
        call fphi_mpole_site(iglob,i)
        call fphi_to_cphi_site(fphirec(1,i),cphirec(1,i))
        if (repart(iglob).eq.rank) then
          call amove(10,cphirec(1,i),cphi(1,iloc))
        end if
      end do
c$omp end parallel do
      deallocate (qgridmpi)
      deallocate (reqinit)
      deallocate (reqbcast)
      deallocate (req)
      return
      end
c
c
      subroutine tmatxbrecip2(mu,murec,nrhs,dipfield1,dipfield2,
     $  dipfield1bis,dipfield2bis)
c
c     Compute the reciprocal space contribution to the electric field due to the current  
c     value of the induced dipoles
c
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'fft.i'
      include 'math.i'
      include 'mpole.i'
      include 'pme.i'
      include 'openmp.i'
      include 'mpif.h'
      integer ierr,iglob,iloc
      integer status(MPI_STATUS_SIZE),tag
      integer nrhs,itot,iipole
      integer temp
      integer i,j,k
      integer nf1,nf2,nf3
      integer maxk,proc
      real*8 fuind(3),fuinp(3)
      real*8 term
      real*8 a(3,3)
      real*8 fdip_phi1(10), fdip_phi2(10), fdip_sum_phi(20)
      real*8 dipfield1(3,*),dipfield2(3,*)
      real*8 dipfield1bis(3,*),dipfield2bis(3,*)
      real*8 mu(3,nrhs,*),murec(3,nrhs,*)
      integer, allocatable :: req(:),reqbcast(:)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
c
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (req(nproc*nproc))
      allocate (reqbcast(nproc*nproc))
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     zero out the PME charge grid
c
      qgrid2in_2d = 0d0
c     fill the pme grid, loop over the induced dipoles sites
      do j = 1, 3
        a(1,j) = dble(nfft1) * recip(j,1)
        a(2,j) = dble(nfft2) * recip(j,2)
        a(3,j) = dble(nfft3) * recip(j,3)
      end do
c
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        iloc  = poleloc(iipole)
c
c       Convert cartesian dipoles to fractional coordinates
c
        if (repart(iglob).ne.rank) then
          do k = 1, 3
             fuind(k) = a(k,1)*murec(1,1,i) + a(k,2)*murec(2,1,i)
     &                       + a(k,3)*murec(3,1,i)
             fuinp(k) = a(k,1)*murec(1,2,i) + a(k,2)*murec(2,2,i)
     &                       + a(k,3)*murec(3,2,i)
          end do
        else
          do k = 1, 3
             fuind(k) = a(k,1)*mu(1,1,iloc) + a(k,2)*mu(2,1,iloc)
     &                       + a(k,3)*mu(3,1,iloc)
             fuinp(k) = a(k,1)*mu(1,2,iloc) + a(k,2)*mu(2,2,iloc)
     &                       + a(k,3)*mu(3,2,iloc)
          end do
        end if
c
c     assign PME grid
c
        call grid_uind_site(iglob,i,fuind,fuinp,qgrid2in_2d)
      end do
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nproc*rank + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,MPI_COMM_WORLD,req(tag),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nproc*prec_send(i) + rank + 1
        call MPI_ISEND(qgrid2in_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,MPI_COMM_WORLD,
     $   req(tag),ierr)
      end do
c
      do i = 1, nrec_recep
        tag = nproc*rank + prec_recep(i) + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
      do i = 1, nrec_send
        tag = nproc*prec_send(i) + rank + 1
        call MPI_WAIT(req(tag),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
        call aadd(2*n1mpimax*n2mpimax*n3mpimax,
     $   qgrid2in_2d(1,1,1,1,1),qgridmpi(1,1,1,1,i),
     $   qgrid2in_2d(1,1,1,1,1))
      end do
c
c     Perform 3-D FFT forward transform
c
      call fft2d_frontmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     complete the transformation of the charge grid
c
c$omp parallel do default(shared) private(k,j,i,term)
      do k = 1, ksize2(rank+1)
         do j = 1, jsize2(rank+1)
           do i = 1, isize2(rank+1)
             term = qfac_2d(i,j,k)
             qgrid2out_2d(1,i,j,k) = term*qgrid2out_2d(1,i,j,k)
             qgrid2out_2d(2,i,j,k) = term*qgrid2out_2d(2,i,j,k)
           end do
         end do
      end do
c$omp end parallel do
c
c     perform 3-D FFT backward transform
c
      call fft2d_backmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     MPI : Begin reception
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nproc*rank + prec_send(i) + 1
        call MPI_IRECV(qgrid2in_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $   prec_send(i),tag,MPI_COMM_WORLD,reqbcast(tag),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nproc*prec_recep(i) + rank + 1
        call MPI_ISEND(qgrid2in_2d,
     $   2*n1mpimax*n2mpimax*n3mpimax,
     $   MPI_REAL8,prec_recep(i),tag,MPI_COMM_WORLD,reqbcast(tag),ierr)
      end do
c
      do i = 1, nrec_send
        tag = nproc*rank + prec_send(i) + 1
        call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
      do i = 1, nrec_recep
        tag = nproc*prec_recep(i) + rank + 1
        call MPI_WAIT(reqbcast(tag),status,ierr)
      end do
c
c     get fields
c
c$omp parallel do default(shared) private(i,iipole,iloc,iglob,k)
c$omp+ private(fdip_phi1,fdip_phi2,fdip_sum_phi)
      do i = 1, npolerecloc
        iipole = polerecglob(i)
        iglob = ipole(iipole)
        iloc  = poleloc(iipole)
        call fphi_uind_site(iglob,i,fdip_phi1,fdip_phi2,fdip_sum_phi)
        if (repart(iglob).ne.rank) then
c
c     convert the dipole fields from fractional to Cartesian
c
          do k = 1, 3
             dipfield1bis(k,i) = a(k,1)*fdip_phi1(2)
     &                           + a(k,2)*fdip_phi1(3)
     &                           + a(k,3)*fdip_phi1(4)
             dipfield2bis(k,i) = a(k,1)*fdip_phi2(2)
     &                           + a(k,2)*fdip_phi2(3)
     $                           + a(k,3)*fdip_phi2(4)
          end do
        else
          do k = 1, 3
             dipfield1(k,iloc) = a(k,1)*fdip_phi1(2)
     &                           + a(k,2)*fdip_phi1(3)
     &                           + a(k,3)*fdip_phi1(4)
             dipfield2(k,iloc) = a(k,1)*fdip_phi2(2)
     &                           + a(k,2)*fdip_phi2(3)
     $                           + a(k,3)*fdip_phi2(4)
          end do
        end if
      end do
c$omp end parallel do
      deallocate (qgridmpi)
      deallocate (req)
      deallocate (reqbcast)
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ulspred  --  induced dipole prediction coeffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ulspred" uses standard extrapolation or a least squares fit
c     to set coefficients of an induced dipole predictor polynomial
c
c     literature references:
c
c     J. Kolafa, "Time-Reversible Always Stable Predictor-Corrector
c     Method for Molecular Dynamics of Polarizable Molecules", Journal
c     of Computational Chemistry, 25, 335-342 (2004)
c
c     W. Wang and R. D. Skeel, "Fast Evaluation of Polarizable Forces",
c     Journal of Chemical Physics, 123, 164107 (2005)
c
c
      subroutine ulspred
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'uprior.i'
      integer i,j,k,m
      real*8 coeff,udk,upk
      real*8 amax,apmax
      real*8 b(maxualt)
      real*8 bp(maxualt)
      real*8 a(maxualt*(maxualt+1)/2)
      real*8 ap(maxualt*(maxualt+1)/2)
      real*8 c(maxualt,maxualt)
      real*8 cp(maxualt,maxualt)
c
c
c     set the Gear predictor binomial coefficients
c
      if (polpred .eq. 'GEAR') then
         do i = 1, nualt
            coeff = gear(i)
            bpred(i) = coeff
            bpredp(i) = coeff
            bpreds(i) = coeff
            bpredps(i) = coeff
         end do
c
c     set always stable predictor-corrector (ASPC) coefficients
c
      else if (polpred .eq. 'ASPC') then
         do i = 1, nualt
            coeff = aspc(i)
            bpred(i) = coeff
            bpredp(i) = coeff
            bpreds(i) = coeff
            bpredps(i) = coeff
         end do
c
c     derive normal equations corresponding to least squares fit
c
      else
         do k = 1, nualt
            b(k) = 0.0d0
            bp(k) = 0.0d0
            do m = k, nualt
               c(k,m) = 0.0d0
               cp(k,m) = 0.0d0
            end do
         end do
         do i = 1, npole
            do j = 1, 3
               do k = 1, nualt
                  udk = udalt(k,j,i)
                  upk = upalt(k,j,i)
                  do m = k, nualt
                     c(k,m) = c(k,m) + udk*udalt(m,j,i)
                     cp(k,m) = cp(k,m) + upk*upalt(m,j,i)
                  end do
               end do
            end do
         end do
         i = 0
         do k = 2, nualt
            b(k-1) = c(1,k)
            bp(k-1) = cp(1,k)
            do m = k, nualt
               i = i + 1
               a(i) = c(k,m)
               ap(i) = cp(k,m)
            end do
         end do
c
c     check for nonzero coefficients and solve normal equations
c
         k = nualt - 1
         amax = 0.0d0
         apmax = 0.0d0
         do i = 1, k*(k+1)/2
            amax = max(amax,a(i))
            apmax = max(apmax,ap(i))
         end do
         if (amax .ne. 0.0d0)  call cholesky (k,a,b)
         if (apmax .ne. 0.0d0)  call cholesky (k,ap,bp)
c
c     transfer the final solution to the coefficient vector
c
         do k = 1, nualt-1
            bpred(k) = b(k)
            bpredp(k) = bp(k)
            bpreds(k) = b(k)
            bpredps(k) = bp(k)
         end do
         bpred(nualt) = 0.0d0
         bpredp(nualt) = 0.0d0
         bpreds(nualt) = 0.0d0
         bpredps(nualt) = 0.0d0
      end if
      return
      end
c
c     subroutine commfield : communicate some direct fields (Newton's third law)
c
      subroutine commfield(nrhs,ef)
      implicit none
      include 'atoms.i'
      include 'openmp.i'
      include 'mpole.i'
      include 'potent.i'
      include 'mpif.h'
      integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer, allocatable :: req(:)
      real*8 ef(3,nrhs,*)
      real*8, allocatable :: buffer(:,:,:,:)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = MPI_COMM_WORLD
      end if
c
      allocate (buffer(3,nrhs,max(npoleloc,1),n_send1))
      allocate (req(nproc*nproc))
c      write(*,*) 'n_recep1=',n_recep1,'p_recep1=',p_recep1
c      write(*,*) 'n_send1=',n_send1,'p_send1=',p_send1
c
      do i = 1, n_send1
        tag = nproc*rank + p_send1(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_REAL8,
     $   p_send1(i),tag,commloc,req(tag),ierr)
      end do
c
      do i = 1, n_recep1
        tag = nproc*p_recep1(i) + rank + 1
        call MPI_ISEND(ef(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_REAL8,p_recep1(i),
     $  tag,commloc,req(tag),ierr)
      end do
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
      do i = 1, n_send1
        do j = 1, npoleloc
          do k = 1, nrhs
            ef(1,k,j) = ef(1,k,j) + buffer(1,k,j,i)
            ef(2,k,j) = ef(2,k,j) + buffer(2,k,j,i)
            ef(3,k,j) = ef(3,k,j) + buffer(3,k,j,i)
          end do
        end do
      end do
c
      deallocate (buffer)
      deallocate (req)
      return
      end
