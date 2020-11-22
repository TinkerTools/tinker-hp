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
      subroutine dcinduce_pme
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
      real*8  wtime0, wtime1, wtime2, udsum, upsum
      real*8  dtime0,dtime1,ctime
      real*8  term, xx(1)
      real*8, allocatable :: ef(:,:,:), mu(:,:,:), murec(:,:,:)
      real*8, allocatable :: cphi(:,:)
c
      real*8, allocatable :: buffermpi(:,:), buffermpimu(:,:,:)
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
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
      if (rank.le.ndir-1) then
        allocate (mu(3,nrhs,max(1,npolebloc)))
        mu = 0.0d0
        allocate (murec(3,nrhs,max(npolerecloc,1)))
        murec = 0.0d0
        allocate (ef(3,nrhs,max(1,npolebloc)))
        ef = 0.0d0
        allocate (cphi(10,max(1,npoleloc)))
        cphi = 0.0d0
        if (allocated(cphirec)) deallocate (cphirec)
        allocate (cphirec(10,max(npolerecloc,1)))
        cphirec = 0.0d0
        allocate (buffermpi(10,max(npoleloc,1)))
        buffermpi = 0.0d0
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0d0
      else
        allocate (mu(3,nrhs,max(1,npolebloc)))
        mu = 0.0d0
        allocate (murec(3,nrhs,max(npolerecloc,1)))
        murec = 0.0d0
        if (allocated(fphirec)) deallocate (fphirec)
        allocate (fphirec(20,max(npolerecloc,1)))
        fphirec = 0.0d0
        allocate (cphi(10,max(1,npoleloc)))
        cphi = 0.0d0
        if (allocated(cphirec)) deallocate (cphirec)
        allocate (cphirec(10,max(npolerecloc,1)))
        cphirec = 0.0d0
        allocate (buffermpi(10,max(npolerecloc,1)))
        buffermpi = 0.0d0
        allocate (buffermpimu(3,nrhs,max(npolerecloc,1)))
        buffermpimu = 0.0d0
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
      ctime = 0.0d0
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
        term = (4.0d0/3.0d0) * aewald**3 / sqrtpi

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
      if (rank.le.ndir-1) then
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
c
        call commrecdirdip(nrhs,1,murec,mu,buffermpimu,buffermpimu,
     $   req2rec,req2send)
c
       call commdirdir(nrhs,2,mu,reqrec,reqsend)
       call commrecdirdip(nrhs,2,murec,mu,buffermpimu,buffermpimu,
     $  req2rec,req2send)

      else

        call commrecdirdip(nrhs,0,murec,mu,buffermpimu,buffermpimu,
     $   req2rec,req2send)
        call commrecdirdip(nrhs,2,murec,mu,buffermpimu,buffermpimu,
     $   req2rec,req2send)
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
      if (use_pred) then
         if (rank.le.ndir-1) then
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
      end if
      return
      end


      subroutine inducedc_pme(matvec,nrhs,dodiis,ef,mu,murec)
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
      integer fnd,mdim,inf,iglob,pos,knd,nrhs,info
      real*8  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      real*8, allocatable :: munew(:,:,:), h(:,:,:)
      real*8, allocatable :: xdiis(:,:), ediis(:,:),
     $   bmat(:,:),bloc(:,:), cex(:)
      integer i,ii, j, k, iipole, it, ind, ndismx, nmat, lenb
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real*8  zero, one, rnorm(2), rr, xx(1)
      real*8 term
      real*8 time0,time1,time2,time3
      real*8  dtime0,dtime1,ctime
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
      real*8, allocatable :: buffermpimu(:,:,:)
      integer ierr, tag, proc
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer, allocatable :: reqrecdir(:)
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
      if (rank.le.ndir-1) then
        allocate (buffermpimu(3,nrhs,max(npoleloc,1)))
        buffermpimu = 0.0d0
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0d0
        allocate (dipfield(3,nrhs,max(1,npoleloc)))
        dipfield = 0.0d0
        allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
        dipfieldbis = 0.0d0
        allocate (dcfield1(maxdim))
        allocate (dcfield2(maxdim))
        dcfield1 = 0.0d0
        dcfield2 = 0.0d0
        allocate (h(3,nrhs,npolebloc))
        h = 0.0d0
        if (dodiis) then
          nmat = 1
          lenb = ndismx + 1
          allocate (xdiis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (ediis(3*nrhs*max(1,npoleloc),ndismx))
          allocate (bmat(lenb,lenb))
          allocate (bloc(lenb,lenb))
          allocate (cex(lenb))
          bmat = 0.0d0
        end if
      else
        allocate (munew(3,nrhs,max(1,npolebloc)))
        munew = 0.0d0
        allocate (buffermpimu(3,nrhs,max(1,npolerecloc)))
        buffermpimu = 0.0d0
        allocate (dipfield(3,nrhs,max(1,npoleloc)))
        dipfield = 0.0d0
        allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
        dipfieldbis = 0.0d0
      end if
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = 0.0d0
      allocate (buffermpi2(3,nrhs,max(1,npolerecloc)))
      buffermpi2 = 0.0d0
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
C        ctime = 0.0d0
        rnorm = 0.0d0
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
         term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
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
            rnorm(k)=sum((munew(:,k,1:npoleloc)-mu(:,k,1:npoleloc))**2)
          end do
c
        end if
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_REAL8,
     $    MPI_SUM,MPI_COMM_WORLD,reqnorm,ierr)

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
            cex = 0.0d0
            cex(1) = one
            call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
            munew = 0.0d0
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
         rnorm(k) = sqrt(rnorm(k)/dble(3*npolar))
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


      subroutine pc_dc_efld0_direct(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
      use couple
      use divcon
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,iglob,kglob,kbis,nrhs,inl,ipoleloc
      integer atii,atkk,l,cofst1,cofst2,cofst3,rofst
      integer knd,kbeg,inf,maxrow,mdim,tpos
      real*8  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole
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
      real*8 fkd(3), fkp(3), fkm(3)
      real*8 bcn(2),invpol
      real*8 erfc, cutoff2
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: uscale(:)
      save   zero, pt6, one, f50
      data   zero/0.d0/, pt6/0.6d0/, one/1.d0/,f50/50.d0/
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $  ' update, try lowering nlupdate')
c
      mode = 'EWALD'
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (uscale(n))
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do

      if(allocated(ytab)) deallocate(ytab)
      allocate(ytab(npolelocnl*maxelst*6))
      tpos = 1
c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        l = grplst(iglob)
c diagonal of Z mat.
        if(l .gt. 0) then
          rofst = (atmofst(iglob) - 1)*3
          if (polarity(iipole).ne.0d0)  then
            invpol = 1.0d0/polarity(iipole)
          else
            invpol = 1000.0d0
          end if
          maxrow=npergrp(l)*3
          cofst1 = rofst + 1
          cofst2 = rofst + 2
          cofst3 = rofst + 3
          cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
          cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
          cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
          zmat(rofst+1+cofst1+kofst(l)) = invpol
          zmat(rofst+2+cofst2+kofst(l)) = invpol
          zmat(rofst+3+cofst3+kofst(l)) = invpol
        end if

        if (i.eq.0) then
          write(iout,1000)
          goto 20
        end if
        pdi = pdamp(iipole)
        pti = thole(iipole)
        ci = rpole(1,iipole)
        dix = rpole(2,iipole)
        diy = rpole(3,iipole)
        diz = rpole(4,iipole)
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

        do j = 1, np11(iglob)
           uscale(ip11(j,iglob)) = u1scale
        end do
        do j = 1, np12(iglob)
           uscale(ip12(j,iglob)) = u2scale
        end do
        do j = 1, np13(iglob)
           uscale(ip13(j,iglob)) = u3scale
        end do
        do j = 1, np14(iglob)
           uscale(ip14(j,iglob)) = u4scale
        end do
c
        do kkk = 1, nelst(ii)
          kglob = elst(kkk,ii)
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

            bcn(1) = bn(1) - 
     & (1.0d0 - scale3*uscale(kglob))/ (d*d2)
            bcn(2) = bn(2) - 
     & 3.0d0*(1.0d0 - scale5*uscale(kglob))/ (d*d2*d2)

            if(l .eq. grplst(kglob) .and. l .ne. -1 ) then
              atii = (atmofst(iglob) - 1)*3
              atkk = (atmofst(kglob) - 1)*3
              maxrow=npergrp(l)*3
              if(atii .lt. atkk) then
                cofst1 = atii + 1
                cofst2 = atii + 2
                cofst3 = atii + 3
                rofst = atkk
              else
                cofst1 = atkk + 1
                cofst2 = atkk + 2
                cofst3 = atkk + 3
                rofst = atii
              end if

              cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
              cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
              cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
              zmat(rofst+1+cofst1+kofst(l)) = bcn(1) - bcn(2)*dx*dx
              zmat(rofst+2+cofst1+kofst(l)) = -bcn(2)*dx*dy
              zmat(rofst+3+cofst1+kofst(l)) = -bcn(2)*dx*dz
              zmat(rofst+1+cofst2+kofst(l)) = -bcn(2)*dx*dy
              zmat(rofst+2+cofst2+kofst(l)) = bcn(1) - bcn(2)*dy*dy
              zmat(rofst+3+cofst2+kofst(l)) = -bcn(2)*dy*dz
              zmat(rofst+1+cofst3+kofst(l)) = -bcn(2)*dx*dz
              zmat(rofst+2+cofst3+kofst(l)) = -bcn(2)*dy*dz
              zmat(rofst+3+cofst3+kofst(l)) = bcn(1) - bcn(2)*dz*dz    
              else
              ytab(tpos) = bcn(1) - bcn(2)*dx*dx
              tpos = tpos + 1
              ytab(tpos) = -bcn(2)*dx*dy
              tpos = tpos + 1
              ytab(tpos) = -bcn(2)*dx*dz
              tpos = tpos + 1
              ytab(tpos) = bcn(1) - bcn(2)*dy*dy
              tpos = tpos + 1
              ytab(tpos) = -bcn(2)*dy*dz
              tpos = tpos + 1
              ytab(tpos) = bcn(1) - bcn(2)*dz*dz 
              tpos = tpos + 1          
            end if

          end if
        end do
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
        do j = 1, np11(iglob)
           uscale(ip11(j,iglob)) = 1.0d0
        end do
        do j = 1, np12(iglob)
           uscale(ip12(j,iglob)) = 1.0d0
        end do
        do j = 1, np13(iglob)
           uscale(ip13(j,iglob)) = 1.0d0
        end do
        do j = 1, np14(iglob)
           uscale(ip14(j,iglob)) = 1.0d0
        end do
 20     continue
      end do


c
      deallocate (dscale)
      deallocate (pscale)
c
      return
      end

      subroutine pc_dc_efld0_shortreal(nrhs,ef)
c
c     Compute the short range direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
      use couple
      use divcon
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,iglob,kglob,kbis,nrhs,inl,ipoleloc
      integer atii,atkk,l,cofst1,cofst2,cofst3,rofst
      integer knd,kbeg,inf,maxrow,mdim,tpos
      real*8  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole
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
      real*8 fkd(3), fkp(3), fkm(3)
      real*8 bcn(2),invpol
      real*8 erfc, cutoff2
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: uscale(:)
      save   zero, pt6, one, f50
      data   zero/0.d0/, pt6/0.6d0/, one/1.d0/,f50/50.d0/
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $  ' update, try lowering nlupdate')
c
      mode = 'SHORTEWALD'
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (uscale(n))
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do

      if(allocated(ytab)) deallocate(ytab)
      allocate(ytab(npolelocnl*maxelst*6))
      tpos = 1
c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        l = grplst(iglob)
c diagonal of Z mat.
        if(l .gt. 0) then
          rofst = (atmofst(iglob) - 1)*3
          if (polarity(iipole).ne.0d0)  then
            invpol = 1.0d0/polarity(iipole)
          else
            invpol = 1000.0d0
          end if
          maxrow=npergrp(l)*3
          cofst1 = rofst + 1
          cofst2 = rofst + 2
          cofst3 = rofst + 3
          cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
          cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
          cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
          zmat(rofst+1+cofst1+kofst(l)) = invpol
          zmat(rofst+2+cofst2+kofst(l)) = invpol
          zmat(rofst+3+cofst3+kofst(l)) = invpol
        end if

        if (i.eq.0) then
          write(iout,1000)
          goto 20
        end if
        pdi = pdamp(iipole)
        pti = thole(iipole)
        ci = rpole(1,iipole)
        dix = rpole(2,iipole)
        diy = rpole(3,iipole)
        diz = rpole(4,iipole)
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

        do j = 1, np11(iglob)
           uscale(ip11(j,iglob)) = u1scale
        end do
        do j = 1, np12(iglob)
           uscale(ip12(j,iglob)) = u2scale
        end do
        do j = 1, np13(iglob)
           uscale(ip13(j,iglob)) = u3scale
        end do
        do j = 1, np14(iglob)
           uscale(ip14(j,iglob)) = u4scale
        end do
c
        do kkk = 1, nelst(ii)
          kglob = elst(kkk,ii)
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

            bcn(1) = bn(1) - 
     & (1.0d0 - scale3*uscale(kglob))/ (d*d2)
            bcn(2) = bn(2) - 
     & 3.0d0*(1.0d0 - scale5*uscale(kglob))/ (d*d2*d2)

            if(l .eq. grplst(kglob) .and. l .ne. -1 ) then
              atii = (atmofst(iglob) - 1)*3
              atkk = (atmofst(kglob) - 1)*3
              maxrow=npergrp(l)*3
              if(atii .lt. atkk) then
                cofst1 = atii + 1
                cofst2 = atii + 2
                cofst3 = atii + 3
                rofst = atkk
              else
                cofst1 = atkk + 1
                cofst2 = atkk + 2
                cofst3 = atkk + 3
                rofst = atii
              end if

              cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
              cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
              cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
              zmat(rofst+1+cofst1+kofst(l)) = bcn(1) - bcn(2)*dx*dx
              zmat(rofst+2+cofst1+kofst(l)) = -bcn(2)*dx*dy
              zmat(rofst+3+cofst1+kofst(l)) = -bcn(2)*dx*dz
              zmat(rofst+1+cofst2+kofst(l)) = -bcn(2)*dx*dy
              zmat(rofst+2+cofst2+kofst(l)) = bcn(1) - bcn(2)*dy*dy
              zmat(rofst+3+cofst2+kofst(l)) = -bcn(2)*dy*dz
              zmat(rofst+1+cofst3+kofst(l)) = -bcn(2)*dx*dz
              zmat(rofst+2+cofst3+kofst(l)) = -bcn(2)*dy*dz
              zmat(rofst+3+cofst3+kofst(l)) = bcn(1) - bcn(2)*dz*dz    
              else
              ytab(tpos) = bcn(1) - bcn(2)*dx*dx
              tpos = tpos + 1
              ytab(tpos) = -bcn(2)*dx*dy
              tpos = tpos + 1
              ytab(tpos) = -bcn(2)*dx*dz
              tpos = tpos + 1
              ytab(tpos) = bcn(1) - bcn(2)*dy*dy
              tpos = tpos + 1
              ytab(tpos) = -bcn(2)*dy*dz
              tpos = tpos + 1
              ytab(tpos) = bcn(1) - bcn(2)*dz*dz 
              tpos = tpos + 1          
            end if

          end if
        end do
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
        do j = 1, np11(iglob)
           uscale(ip11(j,iglob)) = 1.0d0
        end do
        do j = 1, np12(iglob)
           uscale(ip12(j,iglob)) = 1.0d0
        end do
        do j = 1, np13(iglob)
           uscale(ip13(j,iglob)) = 1.0d0
        end do
        do j = 1, np14(iglob)
           uscale(ip14(j,iglob)) = 1.0d0
        end do
 20     continue
      end do


c
      deallocate (dscale)
      deallocate (pscale)
c
      return
      end

      subroutine dc_factor
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use divcon
      use domdec
      implicit none
      integer k,mdim,knd,kbeg,inf
c     loop over the k blocks and perform a cholesky factorization
      do k = 1, km
        mdim = npergrp(k)*3
        if(mdim .gt. 0)then
          knd = kofst(k) + mdim*(mdim+1)/2 
          kbeg=kofst(k) + 1

          call DPPTRF("L",mdim,zmat(kbeg:knd),inf)

          if(inf .ne. 0)then
          print*,rank,k,"error ",inf,"with factorization of zmat"
          end if
        end if
      end do
      return
      end

      subroutine pc_dc_tmatxb_pme(nrhs,dodiag,mu,efi)
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
      character*10 mode

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
          if (polarity(iipole) .eq. 0d0) then
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

      return
      end

c
      subroutine otf_dc_tmatxb_pme(nrhs,dodiag,mu,efi)
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
      character*10 mode
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
        inl = polelocnl(iipole)
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
          if (polarity(iipole) .eq. 0d0) then
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
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end

      subroutine otf_dc_efld0_direct(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
      use couple
      use divcon
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,iglob,kglob,kbis,nrhs,inl,ipoleloc
      integer atii,atkk,l,cofst1,cofst2,cofst3,rofst
      integer knd,kbeg,inf,maxrow,mdim
      real*8  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole
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
      real*8 fkd(3), fkp(3), fkm(3)
      real*8 bcn(2),invpol
      real*8 erfc, cutoff2
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      save   zero, pt6, one, f50
      data   zero/0.d0/, pt6/0.6d0/, one/1.d0/,f50/50.d0/
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $  ' update, try lowering nlupdate')
c
      mode = 'EWALD'
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do

c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        l = grplst(iglob)
c diagonal of Z mat.
        if(l .gt. 0) then
          rofst = (atmofst(iglob) - 1)*3
          if (polarity(iipole).ne.0d0) then
            invpol = 1.0d0/polarity(iipole)
          else
            invpol = 1000.0d0
          end if
          maxrow=npergrp(l)*3
          cofst1 = rofst + 1
          cofst2 = rofst + 2
          cofst3 = rofst + 3
          cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
          cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
          cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
          zmat(rofst+1+cofst1+kofst(l)) = invpol
          zmat(rofst+2+cofst2+kofst(l)) = invpol
          zmat(rofst+3+cofst3+kofst(l)) = invpol
        end if

        if (i.eq.0) then
          write(iout,1000)
          goto 20
        end if
        pdi = pdamp(iipole)
        pti = thole(iipole)
        ci = rpole(1,iipole)
        dix = rpole(2,iipole)
        diy = rpole(3,iipole)
        diz = rpole(4,iipole)
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
        do kkk = 1, nelst(ii)
          kglob = elst(kkk,ii)
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

            if(l .eq. grplst(kglob) .and. l .ne. -1 ) then

              bcn(1) = bn(1) - (1.0d0 - scale3)/ (d*d2)
              bcn(2) = bn(2) - 3.0d0*(1.0d0 - scale5)/ (d*d2*d2)
              atii = (atmofst(iglob) - 1)*3
              atkk = (atmofst(kglob) - 1)*3
              maxrow=npergrp(l)*3
              if(atii .lt. atkk) then
                cofst1 = atii + 1
                cofst2 = atii + 2
                cofst3 = atii + 3
                rofst = atkk
              else
                cofst1 = atkk + 1
                cofst2 = atkk + 2
                cofst3 = atkk + 3
                rofst = atii
              end if

              cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
              cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
              cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
              zmat(rofst+1+cofst1+kofst(l)) = bcn(1) - bcn(2)*dx*dx
              zmat(rofst+2+cofst1+kofst(l)) = -bcn(2)*dx*dy
              zmat(rofst+3+cofst1+kofst(l)) = -bcn(2)*dx*dz
              zmat(rofst+1+cofst2+kofst(l)) = -bcn(2)*dx*dy
              zmat(rofst+2+cofst2+kofst(l)) = bcn(1) - bcn(2)*dy*dy
              zmat(rofst+3+cofst2+kofst(l)) = -bcn(2)*dy*dz
              zmat(rofst+1+cofst3+kofst(l)) = -bcn(2)*dx*dz
              zmat(rofst+2+cofst3+kofst(l)) = -bcn(2)*dy*dz
              zmat(rofst+3+cofst3+kofst(l)) = bcn(1) - bcn(2)*dz*dz              
            end if

          end if
        end do
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
 20     continue
      end do

c
      deallocate (dscale)
      deallocate (pscale)
c
      return
      end

      subroutine otf_dc_efld0_shortreal(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
      use couple
      use divcon
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,iglob,kglob,kbis,nrhs,inl,ipoleloc
      integer atii,atkk,l,cofst1,cofst2,cofst3,rofst
      integer knd,kbeg,inf,maxrow,mdim
      real*8  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole
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
      real*8 fkd(3), fkp(3), fkm(3)
      real*8 bcn(2),invpol
      real*8 erfc, cutoff2
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      save   zero, pt6, one, f50
      data   zero/0.d0/, pt6/0.6d0/, one/1.d0/,f50/50.d0/
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $  ' update, try lowering nlupdate')
c
      mode = 'SHORTEWALD'
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do

c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        l = grplst(iglob)
c diagonal of Z mat.
        if(l .gt. 0) then
          rofst = (atmofst(iglob) - 1)*3
          if (polarity(iipole).ne.0d0) then
            invpol = 1.0d0/polarity(iipole)
          else
            invpol = 1000.0d0
          end if
          maxrow=npergrp(l)*3
          cofst1 = rofst + 1
          cofst2 = rofst + 2
          cofst3 = rofst + 3
          cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
          cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
          cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
          zmat(rofst+1+cofst1+kofst(l)) = invpol
          zmat(rofst+2+cofst2+kofst(l)) = invpol
          zmat(rofst+3+cofst3+kofst(l)) = invpol
        end if

        if (i.eq.0) then
          write(iout,1000)
          goto 20
        end if
        pdi = pdamp(iipole)
        pti = thole(iipole)
        ci = rpole(1,iipole)
        dix = rpole(2,iipole)
        diy = rpole(3,iipole)
        diz = rpole(4,iipole)
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
        do kkk = 1, nelst(ii)
          kglob = elst(kkk,ii)
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

            if(l .eq. grplst(kglob) .and. l .ne. -1 ) then

              bcn(1) = bn(1) - (1.0d0 - scale3)/ (d*d2)
              bcn(2) = bn(2) - 3.0d0*(1.0d0 - scale5)/ (d*d2*d2)
              atii = (atmofst(iglob) - 1)*3
              atkk = (atmofst(kglob) - 1)*3
              maxrow=npergrp(l)*3
              if(atii .lt. atkk) then
                cofst1 = atii + 1
                cofst2 = atii + 2
                cofst3 = atii + 3
                rofst = atkk
              else
                cofst1 = atkk + 1
                cofst2 = atkk + 2
                cofst3 = atkk + 3
                rofst = atii
              end if

              cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
              cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
              cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
              zmat(rofst+1+cofst1+kofst(l)) = bcn(1) - bcn(2)*dx*dx
              zmat(rofst+2+cofst1+kofst(l)) = -bcn(2)*dx*dy
              zmat(rofst+3+cofst1+kofst(l)) = -bcn(2)*dx*dz
              zmat(rofst+1+cofst2+kofst(l)) = -bcn(2)*dx*dy
              zmat(rofst+2+cofst2+kofst(l)) = bcn(1) - bcn(2)*dy*dy
              zmat(rofst+3+cofst2+kofst(l)) = -bcn(2)*dy*dz
              zmat(rofst+1+cofst3+kofst(l)) = -bcn(2)*dx*dz
              zmat(rofst+2+cofst3+kofst(l)) = -bcn(2)*dy*dz
              zmat(rofst+3+cofst3+kofst(l)) = bcn(1) - bcn(2)*dz*dz              
            end if

          end if
        end do
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
 20     continue
      end do

c
      deallocate (dscale)
      deallocate (pscale)
c
      return
      end

c
c     subroutine commfield : communicate some direct fields (Newton's third law)
c
      subroutine commfield2(nrhs,ef)
      use atoms
      use domdec
      use mpole
      use potent
      use mpi
      implicit none
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

      call dc_factor
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

c
c     subroutine commfield : communicate some direct fields (Newton's third law)
c
cnon  subroutine commfield2short(nrhs,ef)
cnon  use atoms
cnon  use domdec
cnon  use mpole
cnon  use potent
cnon  use mpi
cnon  implicit none
cnon  integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
cnon  integer ierr,commloc
cnon  integer, allocatable :: req(:)
cnon  real*8 ef(3,nrhs,*)
cnon  real*8, allocatable :: buffer(:,:,:,:)
cnon
cnon  if (use_pmecore) then
cnon    commloc = comm_dir
cnon  else
cnon    commloc = MPI_COMM_WORLD
cnon  end if
cnon
cnon  allocate (buffer(3,nrhs,max(npoleloc,1),n_send1))
cnon  allocate (req(nproc*nproc))
cnon   write(*,*) 'n_recep1=',n_recep1,'p_recep1=',p_recep1
cnon   write(*,*) 'n_send1=',n_send1,'p_send1=',p_send1
cnon
cnon  do i = 1, n_sendshort1
cnon    tag = nproc*rank + p_sendshort1(i) + 1
cnon    call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_REAL8,
cnon $   p_sendshort1(i),tag,commloc,req(tag),ierr)
cnon  end do
cnon
cnon  do i = 1, n_recepshort1
cnon    tag = nproc*p_recepshort1(i) + rank + 1
cnon    call MPI_ISEND(ef(1,1,bufbegpole(p_recepshort1(i)+1)),
cnon $   3*nrhs*domlen(p_recepshort1(i)+1),MPI_REAL8,p_recepshort1(i),
cnon $  tag,commloc,req(tag),ierr)
cnon  end do
cnon
cnon  call dc_factor
cnon
cnon  do i = 1, n_sendshort1
cnon    tag = nproc*rank + p_sendshort1(i) + 1
cnon    call MPI_WAIT(req(tag),status,ierr)
cnon  end do
cnon  do i = 1, n_recepshort1
cnon    tag = nproc*p_recepshort1(i) + rank + 1
cnon    call MPI_WAIT(req(tag),status,ierr)
cnon  end do
cnon
cnon  do i = 1, n_sendshort1
cnon    do j = 1, npoleloc
cnon      do k = 1, nrhs
cnon        ef(1,k,j) = ef(1,k,j) + buffer(1,k,j,i)
cnon        ef(2,k,j) = ef(2,k,j) + buffer(2,k,j,i)
cnon        ef(3,k,j) = ef(3,k,j) + buffer(3,k,j,i)
cnon      end do
cnon    end do
cnon  end do
cnon
cnon  deallocate (buffer)
cnon  deallocate (req)
cnon  return
cnon

