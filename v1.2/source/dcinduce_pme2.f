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
      subroutine dcinduce_pme2
      use atmlst
      use divcon
      use domdec
      use ewald
      use inform
      use iounit
      use math
      use mpole
      use mutant
      use pme
      use polar
      use polpot
      use potent
      use units
      use uprior
      use mpi
      implicit none
c
c     without separate cores for reciprocal part
c      
      integer i, j, k, nrhs
c
c     MPI
c
      integer iglob, iipole
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)


c
      parameter (nrhs=2)
      real*8  wtime0, wtime1, wtime2, udsum, upsum
      real*8  term
      real*8 plambda,dplambdadelambdae
      real*8, allocatable :: ef(:,:,:), mu(:,:,:), murec(:,:,:)
      real*8, allocatable :: cphi(:,:)
c
      real*8, allocatable :: buffermpi1(:,:),buffermpi2(:,:)
      real*8, allocatable :: buffermpimu1(:,:,:),buffermpimu2(:,:,:)
c
c     Get direct space induced field for divide and conquer methods
c
      external pc_dc_tmatxb_pme
      external otf_dc_tmatxb_pme
c
      if (.not.use_polar) return
c
c1000 format(' illegal polalg in newinduce.')
 1010 format(' time for the ',a,F14.5)
 1020 format(' total elapsed time in dcinduce_pme2: ',F14.5)
c
      if (deb_Path) write(iout,*), 'dcinduce_pme2 '
c
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
      wtime0 = mpi_wtime()
c
c    compute the reciprocal space contribution (fields)
c
      call efld0_recip(cphi)
c
c     Begin the reception of the reciprocal fields
c
      call commrecdirfields(0,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      call commrecdirfields(1,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
      call commrecdirfields(2,cphirec,cphi,buffermpi1,buffermpi2,
     $ reqrecdirrec,reqrecdirsend)
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
c     get derivative (wrt elambda) of total permanent field for lambda dynamics
c
      if (use_lambdadyn) then
c
        call commfield(nrhs,deflambda)
      end if
c
      call commdirdir(nrhs,0,mu,reqrec,reqsend)
c
c       Add direct and reciprocal fields
c
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npoleloc
        iipole = poleglob(i)
        iglob = ipole(iipole)
        do j = 1, 3
          ef(j,1,i)  = ef(j,1,i) - cphi(j+1,i) +
     $       term*rpole(j+1,iipole)
          ef(j,2,i)  = ef(j,2,i) - cphi(j+1,i) +
     $       term*rpole(j+1,iipole)
        end do
        if (use_lambdadyn.and.mut(iglob).and.elambda.gt.0) then
          do j = 1, 3
            deflambda(j,1,i)  = deflambda(j,1,i) +
     $       term*rpole(j+1,iipole)/elambda
            deflambda(j,2,i)  = deflambda(j,2,i) +
     $       term*rpole(j+1,iipole)/elambda
          end do
        end if
      end do
      wtime1 = mpi_wtime()

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

      call commdirdir(nrhs,1,mu,reqrec,reqsend)
      call commdirdir(nrhs,2,mu,reqrec,reqsend)
      call commrecdirdip(nrhs,0,murec,mu,buffermpimu1,buffermpimu2,
     $ req2rec,req2send)
      call commrecdirdip(nrhs,1,murec,mu,buffermpimu1,buffermpimu2,
     $ req2rec,req2send)
      call commrecdirdip(nrhs,2,murec,mu,buffermpimu1,buffermpimu2,
     $ req2rec,req2send)
c
c     now, call the proper solver.
c
      if(precomp .eq. 1) then
        call inducedc_pme2(pc_dc_tmatxb_pme,nrhs,.true.,ef,mu,murec)
      else
        call inducedc_pme2(otf_dc_tmatxb_pme,nrhs,.true.,ef,mu,murec)
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
      if ((use_pred).and..not.(use_lambdadyn)) then
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
      deallocate (mu)
      deallocate (murec)
      deallocate (cphi)
      return
      end

      subroutine inducedc_pme2(matvec,nrhs,dodiis,ef,mu,murec)
      use atmlst
      use divcon
      use domdec
      use ewald
      use inform
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
      integer nrhs, info, knd,pos,mdim,inf
      real*8  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      real*8, allocatable :: munew(:,:,:),h(:,:,:)
      real*8, allocatable :: xdiis(:,:), ediis(:,:),
     $   bmat(:,:),bloc(:,:), cex(:)
      integer i,ii, j, k, it, ind, ndismx, nmat, lenb
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real*8  zero, one, rnorm(2), rr, xx(1)
      real*8 term
      save    zero, one, xx
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
      if (deb_Path) write(iout,*), 'inducedc_pme2 '
c
c
      zero = 0.0d0
      one  = 1.0d0
      xx   = 0.0d0
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
        call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
        call matvec(nrhs,.false.,mu,h)
        call commfield(nrhs,h)

        call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c
        call commdirdir(nrhs,0,mu,reqrec,reqsend)
        call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c
        call commrecdirdip(nrhs,0,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
c
c     jacobi step:
c
        call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
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
     $    MPI_SUM,COMM_TINKER,reqnorm,ierr)

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
     $        bmat,nmat,reqdiis,COMM_TINKER)
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
c
        call commdirdir(nrhs,1,munew,reqrec,reqsend)
        call commrecdirdip(nrhs,1,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
        call commdirdir(nrhs,2,mu,reqrec,reqsend)
        call commrecdirdip(nrhs,2,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
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
