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
      if (rank.le.ndir-1) then
       if (use_pred .and. nualt.eq.maxualt) then
         call ulspred
         do i = 1, npoleloc
           iipole = poleglob(i)
           do j = 1, 3
             udsum = 0.0_ti_p
             upsum = 0.0_ti_p
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
      if ((use_pred).and..not.(use_lambdadyn)) then
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
      use tinheader
      use timestat
      use mpi
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
          call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
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
         call matvec(nrhs,.false.,mu,h)
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
            call fatal_device("inducedc_pme")
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
            call fatal_device("inducedc_pme")
#else
            call M_gesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#endif
            munew = 0.0_ti_p
            call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
          end if
          call commdirdir(nrhs,1,munew,reqrec,reqsend)
          call commrecdirdip(nrhs,1,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
c
          call commdirdir(nrhs,2,mu,reqrec,reqsend)
          call commrecdirdip(nrhs,2,murec,munew,buffermpimu,
     $    buffermpimu,req2rec,req2send)
          mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)
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


      subroutine pc_dc_efld0_direct(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
      use bound
      use chgpen
      use couple
      use cutoff
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
      use potent
      use shunt
      use mpi
      implicit none
      integer i,iglob,kglob,nrhs,ipoleloc,nnelst
      real(t_p)  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole,kbis
      integer atii,atkk,l,cofst1,cofst2,cofst3,rofst
      integer maxrow,tpos
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1,rr2
      real(t_p) rr3,rr5,rr7
      real(t_p) rr3i,rr5i,rr7i
      real(t_p) rr3k,rr5k,rr7k
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qiyy,qizz
      real(t_p) qixy,qixz,qiyz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkyy,qkzz
      real(t_p) qkxy,qkxz,qkyz
      real(t_p) dir,dkr
      real(t_p) qix,qiy,qiz,qir
      real(t_p) qkx,qky,qkz,qkr
      real(t_p) fid(3), fip(3)
      real(t_p) fkd(3), fkp(3)
      real(t_p) fiu(9)
      real(t_p) invpol
      real(t_p) cutoff2
      real(t_p) corei,corek
      real(t_p) vali,valk
      real(t_p) alphai,alphak
      real(t_p) dmp3,dmp5,dmp7
      real(t_p) dmp3u,dmp5u
      real(t_p) dmpi(7),dmpk(7)
      real(t_p) dmpik(7),dmpe(7)
      real(t_p) scalek,scaleu
      real(t_p), allocatable :: dscale(:)
      real(t_p), allocatable :: pscale(:)
      real(t_p), allocatable :: uscale(:)
      real(t_p), allocatable :: wscale(:)
      logical shortrange
      character*11 mode
      character*80 :: RoutineName
c
 1000 format(' Warning, system moved too much since last neighbor list
     $  update, try lowering nlupdate')
c
      shortrange = use_polarshortreal
      if (shortrange) then 
         RoutineName='pc_dc_efld0_shortreal'
         mode = 'SHORTEWALD'
      else
         RoutineName='pc_dc_efld0_direct'
         mode = 'EWALD'
      endif
c
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      dscale = 1.0d0
      pscale = 1.0d0
      uscale = 1.0d0
      wscale = 1.0d0

      if(allocated(ytab)) deallocate(ytab)
      allocate(ytab(npolelocnl*maxelst*6))
      tpos = 1

c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        l        = grplst    (iglob)
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

        if ((i.eq.0).or.(i.gt.nbloc)) then
          write(iout,1000)
          cycle
        end if
        ci  = rpole(1,iipole)
        dix = rpole(2,iipole)
        diy = rpole(3,iipole)
        diz = rpole(4,iipole)
        qixx = rpole(5,iipole)
        qixy = rpole(6,iipole)
        qixz = rpole(7,iipole)
        qiyy = rpole(9,iipole)
        qiyz = rpole(10,iipole)
        qizz = rpole(13,iipole)
        if (use_chgpen) then
           corei = pcore(iipole)
           vali = pval(iipole)
           alphai = palpha(iipole)
        end if
c
c     set exclusion coefficients for connected atoms
c
        if (dpequal) then
           do j = 1, n12(iglob)
              pscale(i12(j,iglob)) = p2scale
              do k = 1, np11(iglob)
                 if (i12(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i12(j,iglob)) = p2iscale
              end do
              dscale(i12(j,iglob)) = pscale(i12(j,iglob))
           end do
           do j = 1, n13(iglob)
              pscale(i13(j,iglob)) = p3scale
              do k = 1, np11(iglob)
                 if (i13(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i13(j,iglob)) = p3iscale
              end do
              dscale(i13(j,iglob)) = pscale(i13(j,iglob))
           end do
           do j = 1, n14(iglob)
              pscale(i14(j,iglob)) = p4scale
              do k = 1, np11(iglob)
                  if (i14(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i14(j,iglob)) = p4iscale
              end do
              dscale(i14(j,iglob)) = pscale(i14(j,iglob))
           end do
           do j = 1, n15(iglob)
              pscale(i15(j,iglob)) = p5scale
              do k = 1, np11(iglob)
                 if (i15(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i15(j,iglob)) = p5iscale
              end do
              dscale(i15(j,iglob)) = pscale(i15(j,iglob))
           end do
        else
           do j = 1, n12(iglob)
              pscale(i12(j,iglob)) = p2scale
              do k = 1, np11(iglob)
                 if (i12(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i12(j,iglob)) = p2iscale
              end do
           end do
           do j = 1, n13(iglob)
              pscale(i13(j,iglob)) = p3scale
              do k = 1, np11(iglob)
                 if (i13(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i13(j,iglob)) = p3iscale
              end do
           end do
           do j = 1, n14(iglob)
              pscale(i14(j,iglob)) = p4scale
              do k = 1, np11(iglob)
                  if (i14(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i14(j,iglob)) = p4iscale
              end do
           end do
           do j = 1, n15(iglob)
              pscale(i15(j,iglob)) = p5scale
              do k = 1, np11(iglob)
                 if (i15(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i15(j,iglob)) = p5iscale
              end do
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
        end if
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
        do j = 1, n12(iglob)
           wscale(i12(j,iglob)) = w2scale
        end do
        do j = 1, n13(iglob)
           wscale(i13(j,iglob)) = w3scale
        end do
        do j = 1, n14(iglob)
           wscale(i14(j,iglob)) = w4scale
        end do
        do j = 1, n15(iglob)
           wscale(i15(j,iglob)) = w5scale
        end do
c
        if (shortrange) then
          nnelst = nshortelst(ii)
        else
          nnelst = nelst(ii)
        end if
        do kkk = 1, nnelst
          if (shortrange) then
            kkpole = shortelst(kkk,ii)
          else
            kkpole = elst(kkk,ii)
          end if
          kbis = poleloc(kkpole)
          kglob = ipole(kkpole)
          if ((kbis.eq.0).or.(kbis.gt.npolebloc)) then
            write(iout,1000)
            cycle
          end if
          xr = x(kglob) - x(iglob)
          yr = y(kglob) - y(iglob)
          zr = z(kglob) - z(iglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr* yr + zr*zr
          if (r2.le.cutoff2) then
            r = sqrt(r2)
            rr1 = 1.0d0 / r
            rr2 = rr1 * rr1
            rr3 = rr2 * rr1
            rr5 = 3.0d0 * rr2 * rr3
            rr7 = 5.0d0 * rr2 * rr5
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
c
c     intermediates involving moments and separation distance
c
            dir = dix*xr + diy*yr + diz*zr
            qix = qixx*xr + qixy*yr + qixz*zr
            qiy = qixy*xr + qiyy*yr + qiyz*zr
            qiz = qixz*xr + qiyz*yr + qizz*zr
            qir = qix*xr + qiy*yr + qiz*zr
            dkr = dkx*xr + dky*yr + dkz*zr
            qkx = qkxx*xr + qkxy*yr + qkxz*zr
            qky = qkxy*xr + qkyy*yr + qkyz*zr
            qkz = qkxz*xr + qkyz*yr + qkzz*zr
            qkr = qkx*xr + qky*yr + qkz*zr
c
c     calculate real space Ewald error function damping
c
            call dampewald (7,r,r2,1.0d0,dmpe)
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               call dampthole (iipole,kkpole,7,r,dmpik)
               scalek = dscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dmp7 = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
               fid(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkx + 2.0d0*dmp5*qkx
               fid(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dky + 2.0d0*dmp5*qky
               fid(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkz + 2.0d0*dmp5*qkz
               fkd(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*dix - 2.0d0*dmp5*qix
               fkd(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diy - 2.0d0*dmp5*qiy
               fkd(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diz - 2.0d0*dmp5*qiz
               scalek = pscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dmp7 = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
               fip(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkx + 2.0d0*dmp5*qkx
               fip(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dky + 2.0d0*dmp5*qky
               fip(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkz + 2.0d0*dmp5*qkz
               fkp(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*dix - 2.0d0*dmp5*qix
               fkp(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diy - 2.0d0*dmp5*qiy
               fkp(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diz - 2.0d0*dmp5*qiz

               call dampthole2 (iipole,kkpole,5,r,dmpik)
               scaleu = uscale(kglob)
               dmp3u = dmpe(3) - (1.0d0-scaleu*dmpik(3))*rr3
               dmp5u = dmpe(5) - (1.0d0-scaleu*dmpik(5))*rr5
               fiu(1) = dmp3u - dmp5u*xr*xr
               fiu(2) = -dmp5u*xr*yr
               fiu(3) = -dmp5u*xr*zr
               fiu(4) = -dmp5u*xr*yr
               fiu(5) = dmp3u - dmp5u*yr*yr
               fiu(6) = -dmp5u*yr*zr
               fiu(7) = -dmp5u*xr*zr
               fiu(8) = -dmp5u*yr*zr
               fiu(9) = dmp3u - dmp5u*zr*zr
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               corek = pcore(kkpole)
               valk = pval(kkpole)
               alphak = palpha(kkpole)
               call dampdir (r,alphai,alphak,dmpi,dmpk)
               scalek = dscale(kglob)
               rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
               rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
               rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
               rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
               rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
               rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
               rr3 = dmpe(3) - (1.0d0-scalek)*rr3
               fid(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx
               fid(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fid(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkd(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkd(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkd(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz
               scalek = pscale(kglob)
               rr3 = rr2 * rr1
               rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
               rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
               rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
               rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
               rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
               rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
               rr3 = dmpe(3) - (1.0d0-scalek)*rr3
               fip(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx
               fip(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fip(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkp(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkp(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkp(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz

               call dampmut (r,alphai,alphak,dmpik)
               scaleu = wscale(kglob)
               rr3 = rr2 * rr1
               dmp3u = dmpe(3) - (1.0d0-scaleu*dmpik(3))*rr3
               dmp5u = dmpe(5) - (1.0d0-scaleu*dmpik(5))*rr5
               fiu(1) = dmp3u - dmp5u*xr*xr
               fiu(2) = -dmp5u*xr*yr
               fiu(3) = -dmp5u*xr*zr
               fiu(4) = -dmp5u*xr*yr
               fiu(5) = dmp3u - dmp5u*yr*yr
               fiu(6) = -dmp5u*yr*zr
               fiu(7) = -dmp5u*xr*zr
               fiu(8) = -dmp5u*yr*zr
               fiu(9) = dmp3u - dmp5u*zr*zr
            end if
c
c     increment the field at each site due to this interaction
c
            do j = 1, 3
               ef(j,1,ipoleloc) = ef(j,1,ipoleloc) + fid(j)
               ef(j,1,kbis) = ef(j,1,kbis) + fkd(j)
               ef(j,2,ipoleloc) = ef(j,2,ipoleloc) + fip(j)
               ef(j,2,kbis) = ef(j,2,kbis) + fkp(j)
            end do

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

              zmat(rofst+1+cofst1+kofst(l)) = fiu(1)
              zmat(rofst+2+cofst1+kofst(l)) = fiu(2)
              zmat(rofst+3+cofst1+kofst(l)) = fiu(3)
              zmat(rofst+1+cofst2+kofst(l)) = fiu(4)
              zmat(rofst+2+cofst2+kofst(l)) = fiu(5)
              zmat(rofst+3+cofst2+kofst(l)) = fiu(6)
              zmat(rofst+1+cofst3+kofst(l)) = fiu(7)
              zmat(rofst+2+cofst3+kofst(l)) = fiu(8)
              zmat(rofst+3+cofst3+kofst(l)) = fiu(9)
              else
              ytab(tpos) = fiu(1)
              tpos = tpos + 1
              ytab(tpos) = fiu(2)
              tpos = tpos + 1
              ytab(tpos) = fiu(3)
              tpos = tpos + 1
              ytab(tpos) = fiu(5)
              tpos = tpos + 1
              ytab(tpos) = fiu(6)
              tpos = tpos + 1
              ytab(tpos) = fiu(9)
              tpos = tpos + 1
            end if
          end if
        end do
c
c     reset exclusion coefficients for connected atoms
c
        if (dpequal) then
           do j = 1, n12(iglob)
              pscale(i12(j,iglob)) = 1.0d0
              dscale(i12(j,iglob)) = 1.0d0
           end do
           do j = 1, n13(iglob)
              pscale(i13(j,iglob)) = 1.0d0
              dscale(i13(j,iglob)) = 1.0d0
           end do
           do j = 1, n14(iglob)
              pscale(i14(j,iglob)) = 1.0d0
              dscale(i14(j,iglob)) = 1.0d0
           end do
           do j = 1, n15(iglob)
              pscale(i15(j,iglob)) = 1.0d0
              dscale(i15(j,iglob)) = 1.0d0
           end do
        else
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
        end if
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
        do j = 1, n12(iglob)
           wscale(i12(j,iglob)) = 1.0d0
        end do
        do j = 1, n13(iglob)
           wscale(i13(j,iglob)) = 1.0d0
        end do
        do j = 1, n14(iglob)
           wscale(i14(j,iglob)) = 1.0d0
        end do
        do j = 1, n15(iglob)
           wscale(i15(j,iglob)) = 1.0d0
        end do
      end do
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (uscale)
      deallocate (wscale)

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

#ifdef _OPENACC
          call fatal_device("dc_factor")
#else
          call DPPTRF("L",mdim,zmat(kbeg:knd),inf)
#endif

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
      use sizes
      use atoms
      use atmlst
      use divcon
      use domdec
      use ewald
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use mpi
      implicit none
      integer i,nrhs,iglob,kglob,kbis,iipole,kkpole,kkpoleloc,nnelst
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

      logical shortrange
      character*11 mode
      character*80 :: RoutineName



      zero = 0.0d0
      one  = 1.0d0
      f50  = 50.d0

      shortrange = use_polarshortreal
      if (shortrange) then 
         RoutineName='pc_dc_tmatxb_shortreal'
         mode = 'SHORTEWALD'
      else
         RoutineName='pc_dc_tmatxb_pme'
         mode = 'EWALD'
      endif
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
      call switch (mode)

      cutoff2 = cut2
      tpos = 1
c
      do ii = 1, npolelocnl
        iipole   = poleglobnl(ii)
        iglob    = ipole  (iipole)
        i        = loc    (iglob)
        l        = grplst (iglob)
        ipoleloc = poleloc(iipole)
        if (i.eq.0) cycle
        duix = mu(1,1,ipoleloc)
        duiy = mu(2,1,ipoleloc)
        duiz = mu(3,1,ipoleloc)
        puix = mu(1,2,ipoleloc)
        puiy = mu(2,2,ipoleloc)
        puiz = mu(3,2,ipoleloc)
c
        if (shortrange) then
          nnelst = nshortelst(ii)
        else
          nnelst = nelst(ii)
        end if
        if (shortrange) then
          nnelst = nshortelst(ii)
        else
          nnelst = nelst(ii)
        end if
        do kkk = 1, nnelst
          if (shortrange) then
            kkpole = shortelst(kkk,ii)
          else
            kkpole = elst(kkk,ii)
          end if
c
c     only build induced field for interactons outside of a block
c     check if the atoms are in the same block
c
          kglob = ipole(kkpole)
          if(l .ne. grplst(kglob) .or. l .lt. 0)then
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
      use bound
      use chgpen
      use couple
      use divcon
      use domdec
      use ewald
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use mpi
      implicit none
      integer i,nrhs,iglob,kglob,iipole,kkpole,kkpoleloc,nnelst
      integer ipoleloc,l
      real(t_p)  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1,rr2
      real(t_p) rr3,rr5
      real(t_p) rr3ik,rr5ik
      real(t_p) scalek
      real(t_p) dmp3,dmp5
      real(t_p) fid(3),fkd(3)
      real(t_p) fip(3),fkp(3)
      real(t_p) dmpik(7),dmpe(7)
      real(t_p) dlocal(6)
      real(t_p) duix,duiy,duiz
      real(t_p) puix,puiy,puiz
      real(t_p) dukx,duky,dukz
      real(t_p) pukx,puky,pukz
      real(t_p) alphai,alphak
      real(t_p) pol
      real(t_p), allocatable :: uscale(:)
      real(t_p), allocatable :: wscale(:)
      logical dodiag
      logical shortrange
      integer   j, ii, kkk, irhs
      real(t_p)  cutoff2
      character*11 mode
      character*80 :: RoutineName

      shortrange = use_polarshortreal
      if (shortrange) then 
         RoutineName='tmatxb_shortreal'
         mode = 'SHORTEWALD'
      else
         RoutineName='tmatxb_pme'
         mode = 'EWALD'
      endif
c
c     initialize the result vector
c
      efi = 0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
      allocate (wscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      uscale = 1.0d0
      wscale = 1.0d0
c
c     gather some parameters, then set up the damping factors.
c
      call switch (mode)

      cutoff2 = cut2
c
      do ii = 1, npolelocnl
        iipole   = poleglobnl(ii)
        iglob    = ipole     (iipole)
        i        = loc       (iglob)
        ipoleloc = poleloc(iipole)
        l        = grplst (iglob)
        if (i.eq.0) cycle
        duix = mu(1,1,ipoleloc)
        duiy = mu(2,1,ipoleloc)
        duiz = mu(3,1,ipoleloc)
        puix = mu(1,2,ipoleloc)
        puiy = mu(2,2,ipoleloc)
        puiz = mu(3,2,ipoleloc)
        if (use_chgpen) then
           alphai = palpha(iipole)
        end if
c
c
c     set exclusion coefficients for connected atoms
c
        do j = 1, n12(iglob)
           wscale(i12(j,iglob)) = w2scale
        end do
        do j = 1, n13(iglob)
           wscale(i13(j,iglob)) = w3scale
        end do
        do j = 1, n14(iglob)
           wscale(i14(j,iglob)) = w4scale
        end do
        do j = 1, n15(iglob)
           wscale(i15(j,iglob)) = w5scale
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
        if (shortrange) then
          nnelst = nshortelst(ii)
        else
          nnelst = nelst(ii)
        end if
        do kkk = 1, nnelst
          if (shortrange) then
            kkpole = shortelst(kkk,ii)
          else
            kkpole = elst(kkk,ii)
          end if
          kglob = ipole(kkpole)
c
c     only build induced field for interactions outside of a block
c     check if the atoms are in the same block
c

          if(l .eq. grplst(kglob) .and. l .ge. 0)cycle
          kkpoleloc = poleloc(kkpole)
          if (kkpoleloc.eq.0) cycle
          xr = x(kglob) - x(iglob)
          yr = y(kglob) - y(iglob)
          zr = z(kglob) - z(iglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr* yr + zr*zr
          if (r2 .le. off2) then
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
            r = sqrt(r2)
            rr1 = 1.0d0 / r
            rr2 = rr1 * rr1
            rr3 = rr2 * rr1
            rr5 = 3.0d0 * rr2 * rr3
            dukx = mu(1,1,kkpoleloc)
            duky = mu(2,1,kkpoleloc)
            dukz = mu(3,1,kkpoleloc)
            pukx = mu(1,2,kkpoleloc)
            puky = mu(2,2,kkpoleloc)
            pukz = mu(3,2,kkpoleloc)
            if (use_chgpen) then
               alphak = palpha(kkpole)
            end if
c
c     calculate real space Ewald error function damping
c
            call dampewald (7,r,r2,1.0d0,dmpe)
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               call dampthole2 (iipole,kkpole,5,r,dmpik)
               scalek = uscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dlocal(1) = -dmp3 + dmp5*xr*xr
               dlocal(2) = dmp5*xr*yr
               dlocal(3) = dmp5*xr*zr
               dlocal(4) = -dmp3 + dmp5*yr*yr
               dlocal(5) = dmp5*yr*zr
               dlocal(6) = -dmp3 + dmp5*zr*zr
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               call dampmut (r,alphai,alphak,dmpik)
               scalek = wscale(kglob)
               rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dlocal(1) = -rr3ik + rr5ik*xr*xr
               dlocal(2) = rr5ik*xr*yr
               dlocal(3) = rr5ik*xr*zr
               dlocal(4) = -rr3ik + rr5ik*yr*yr
               dlocal(5) = rr5ik*yr*zr
               dlocal(6) = -rr3ik + rr5ik*zr*zr
            end if
            fid(1) = dlocal(1)*dukx+dlocal(2)*duky+dlocal(3)*dukz
            fid(2) = dlocal(2)*dukx+dlocal(4)*duky+dlocal(5)*dukz
            fid(3) = dlocal(3)*dukx+dlocal(5)*duky+dlocal(6)*dukz
            fkd(1) = dlocal(1)*duix+dlocal(2)*duiy+dlocal(3)*duiz
            fkd(2) = dlocal(2)*duix+dlocal(4)*duiy+dlocal(5)*duiz
            fkd(3) = dlocal(3)*duix+dlocal(5)*duiy+dlocal(6)*duiz

            fip(1) = dlocal(1)*pukx+dlocal(2)*puky+dlocal(3)*pukz
            fip(2) = dlocal(2)*pukx+dlocal(4)*puky+dlocal(5)*pukz
            fip(3) = dlocal(3)*pukx+dlocal(5)*puky+dlocal(6)*pukz
            fkp(1) = dlocal(1)*puix+dlocal(2)*puiy+dlocal(3)*puiz
            fkp(2) = dlocal(2)*puix+dlocal(4)*puiy+dlocal(5)*puiz
            fkp(3) = dlocal(3)*puix+dlocal(5)*puiy+dlocal(6)*puiz
            do j = 1, 3
               efi(j,1,ipoleloc) = efi(j,1,ipoleloc) - fid(j)
               efi(j,1,kkpoleloc) = efi(j,1,kkpoleloc) - fkd(j)
               efi(j,2,ipoleloc) = efi(j,2,ipoleloc) - fip(j)
               efi(j,2,kkpoleloc) = efi(j,2,kkpoleloc) - fkp(j)
            end do
          end if
        end do
c
c     reset interaction scaling coefficients for connected atoms
c
        do j = 1, n12(iglob)
           wscale(i12(j,iglob)) = 1.0d0
        end do
        do j = 1, n13(iglob)
           wscale(i13(j,iglob)) = 1.0d0
        end do
        do j = 1, n14(iglob)
           wscale(i14(j,iglob)) = 1.0d0
        end do
        do j = 1, n15(iglob)
           wscale(i15(j,iglob)) = 1.0d0
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
      end do
      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
        do i = 1, npoleloc
          iipole = poleglob(i)
c
c     if no polarisability, take a negligeable value to allow convergence
c
          if (polarity(iipole).eq.0.0d0) then
             pol = tinypol ** (-1)
          else
             pol  = polarity(iipole) ** (-1)
          endif
          do irhs = 1, nrhs
            do j = 1, 3
              efi(j,irhs,i) = efi(j,irhs,i) +
     $           mu(j,irhs,i)*pol
            end do
          end do
        end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
      deallocate (wscale)
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
      use bound
      use chgpen
      use couple
      use cutoff
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
      use potent
      use shunt
      use mpi
      implicit none
      integer i,iglob,kglob,nrhs,ipoleloc,nnelst
      integer atii,atkk,l,cofst1,cofst2,cofst3,rofst
      integer maxrow
      real(t_p)  ef(3,nrhs,npolebloc)
      integer ii, j, k, kkk, iipole, kkpole,kbis
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1,rr2
      real(t_p) rr3,rr5,rr7
      real(t_p) rr3i,rr5i,rr7i
      real(t_p) rr3k,rr5k,rr7k
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qiyy,qizz
      real(t_p) qixy,qixz,qiyz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkyy,qkzz
      real(t_p) qkxy,qkxz,qkyz
      real(t_p) dir,dkr
      real(t_p) qix,qiy,qiz,qir
      real(t_p) qkx,qky,qkz,qkr
      real(t_p) fid(3), fip(3)
      real(t_p) fkd(3), fkp(3)
      real(t_p) fiu(9)
      real(t_p) cutoff2
      real(t_p) corei,corek
      real(t_p) vali,valk
      real(t_p) alphai,alphak
      real(t_p) dmp3,dmp5,dmp7
      real(t_p) dmp3u,dmp5u
      real(t_p) dmpi(7),dmpk(7)
      real(t_p) dmpik(7),dmpe(7)
      real(t_p) invpol
      real(t_p) scalek
      real(t_p), allocatable :: dscale(:)
      real(t_p), allocatable :: pscale(:)
      real(t_p), allocatable :: uscale(:)
      real(t_p), allocatable :: wscale(:)
      logical shortrange
      character*11 mode
      character*80 :: RoutineName
c
 1000 format(' Warning, system moved too much since last neighbor list
     $  update, try lowering nlupdate')
c
      shortrange = use_polarshortreal
      if (shortrange) then 
         RoutineName='otf_dc_efld0_shortreal'
         mode = 'SHORTEWALD'
      else
         RoutineName='otf_dc_efld0_direct'
         mode = 'EWALD'
      endif
c
      call switch (mode)
      cutoff2 = cut2
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      dscale = 1.0d0
      pscale = 1.0d0
      uscale = 1.0d0
      wscale = 1.0d0

c
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob = ipole(iipole)
        i = loc(iglob)
        ipoleloc = poleloc(iipole)
        l        = grplst    (iglob)
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

        if ((i.eq.0).or.(i.gt.nbloc)) then
          write(iout,1000)
          cycle
        end if
        ci  = rpole(1,iipole)
        dix = rpole(2,iipole)
        diy = rpole(3,iipole)
        diz = rpole(4,iipole)
        qixx = rpole(5,iipole)
        qixy = rpole(6,iipole)
        qixz = rpole(7,iipole)
        qiyy = rpole(9,iipole)
        qiyz = rpole(10,iipole)
        qizz = rpole(13,iipole)
        if (use_chgpen) then
           corei = pcore(iipole)
           vali = pval(iipole)
           alphai = palpha(iipole)
        end if
c
c     set exclusion coefficients for connected atoms
c
        if (dpequal) then
           do j = 1, n12(iglob)
              pscale(i12(j,iglob)) = p2scale
              do k = 1, np11(iglob)
                 if (i12(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i12(j,iglob)) = p2iscale
              end do
              dscale(i12(j,iglob)) = pscale(i12(j,iglob))
           end do
           do j = 1, n13(iglob)
              pscale(i13(j,iglob)) = p3scale
              do k = 1, np11(iglob)
                 if (i13(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i13(j,iglob)) = p3iscale
              end do
              dscale(i13(j,iglob)) = pscale(i13(j,iglob))
           end do
           do j = 1, n14(iglob)
              pscale(i14(j,iglob)) = p4scale
              do k = 1, np11(iglob)
                  if (i14(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i14(j,iglob)) = p4iscale
              end do
              dscale(i14(j,iglob)) = pscale(i14(j,iglob))
           end do
           do j = 1, n15(iglob)
              pscale(i15(j,iglob)) = p5scale
              do k = 1, np11(iglob)
                 if (i15(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i15(j,iglob)) = p5iscale
              end do
              dscale(i15(j,iglob)) = pscale(i15(j,iglob))
           end do
        else
           do j = 1, n12(iglob)
              pscale(i12(j,iglob)) = p2scale
              do k = 1, np11(iglob)
                 if (i12(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i12(j,iglob)) = p2iscale
              end do
           end do
           do j = 1, n13(iglob)
              pscale(i13(j,iglob)) = p3scale
              do k = 1, np11(iglob)
                 if (i13(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i13(j,iglob)) = p3iscale
              end do
           end do
           do j = 1, n14(iglob)
              pscale(i14(j,iglob)) = p4scale
              do k = 1, np11(iglob)
                  if (i14(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i14(j,iglob)) = p4iscale
              end do
           end do
           do j = 1, n15(iglob)
              pscale(i15(j,iglob)) = p5scale
              do k = 1, np11(iglob)
                 if (i15(j,iglob) .eq. ip11(k,iglob))
     &              pscale(i15(j,iglob)) = p5iscale
              end do
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
        end if
        do j = 1, n12(iglob)
           wscale(i12(j,iglob)) = w2scale
        end do
        do j = 1, n13(iglob)
           wscale(i13(j,iglob)) = w3scale
        end do
        do j = 1, n14(iglob)
           wscale(i14(j,iglob)) = w4scale
        end do
        do j = 1, n15(iglob)
           wscale(i15(j,iglob)) = w5scale
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
        if (shortrange) then
          nnelst = nshortelst(ii)
        else
          nnelst = nelst(ii)
        end if
        do kkk = 1, nnelst
          if (shortrange) then
            kkpole = shortelst(kkk,ii)
          else
            kkpole = elst(kkk,ii)
          end if
          kbis = poleloc(kkpole)
          kglob = ipole(kkpole)
          if ((kbis.eq.0).or.(kbis.gt.npolebloc)) then
            write(iout,1000)
            cycle
          end if
          xr = x(kglob) - x(iglob)
          yr = y(kglob) - y(iglob)
          zr = z(kglob) - z(iglob)
          if (use_bounds)  call image (xr,yr,zr)
          r2 = xr*xr + yr* yr + zr*zr
          if (r2.le.cutoff2) then
            r = sqrt(r2)
            rr1 = 1.0d0 / r
            rr2 = rr1 * rr1
            rr3 = rr2 * rr1
            rr5 = 3.0d0 * rr2 * rr3
            rr7 = 5.0d0 * rr2 * rr5
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
c
c     intermediates involving moments and separation distance
c
            dir = dix*xr + diy*yr + diz*zr
            qix = qixx*xr + qixy*yr + qixz*zr
            qiy = qixy*xr + qiyy*yr + qiyz*zr
            qiz = qixz*xr + qiyz*yr + qizz*zr
            qir = qix*xr + qiy*yr + qiz*zr
            dkr = dkx*xr + dky*yr + dkz*zr
            qkx = qkxx*xr + qkxy*yr + qkxz*zr
            qky = qkxy*xr + qkyy*yr + qkyz*zr
            qkz = qkxz*xr + qkyz*yr + qkzz*zr
            qkr = qkx*xr + qky*yr + qkz*zr
c
c     calculate real space Ewald error function damping
c
            call dampewald (7,r,r2,1.0d0,dmpe)
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               call dampthole (iipole,kkpole,7,r,dmpik)
               scalek = dscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dmp7 = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
               fid(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkx + 2.0d0*dmp5*qkx
               fid(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dky + 2.0d0*dmp5*qky
               fid(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkz + 2.0d0*dmp5*qkz
               fkd(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*dix - 2.0d0*dmp5*qix
               fkd(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diy - 2.0d0*dmp5*qiy
               fkd(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diz - 2.0d0*dmp5*qiz
               scalek = pscale(kglob)
               dmp3 = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5 = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               dmp7 = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
               fip(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkx + 2.0d0*dmp5*qkx
               fip(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dky + 2.0d0*dmp5*qky
               fip(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)
     &                     - dmp3*dkz + 2.0d0*dmp5*qkz
               fkp(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*dix - 2.0d0*dmp5*qix
               fkp(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diy - 2.0d0*dmp5*qiy
               fkp(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)
     &                     - dmp3*diz - 2.0d0*dmp5*qiz

               call dampthole2 (iipole,kkpole,5,r,dmpik)
               scalek = uscale(kglob)
               dmp3u = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5u = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               fiu(1) = dmp3u - dmp5u*xr*xr
               fiu(2) = -dmp5u*xr*yr
               fiu(3) = -dmp5u*xr*zr
               fiu(4) = -dmp5u*xr*yr
               fiu(5) = dmp3u - dmp5u*yr*yr
               fiu(6) = -dmp5u*yr*zr
               fiu(7) = -dmp5u*xr*zr
               fiu(8) = -dmp5u*yr*zr
               fiu(9) = dmp3u - dmp5u*zr*zr
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               corek = pcore(kkpole)
               valk = pval(kkpole)
               alphak = palpha(kkpole)
               call dampdir (r,alphai,alphak,dmpi,dmpk)
               scalek = dscale(kglob)
               rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
               rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
               rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
               rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
               rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
               rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
               rr3 = dmpe(3) - (1.0d0-scalek)*rr3
               fid(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx
               fid(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fid(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkd(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkd(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkd(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz
               scalek = pscale(kglob)
               rr3 = rr2 * rr1
               rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
               rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
               rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
               rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
               rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
               rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
               rr3 = dmpe(3) - (1.0d0-scalek)*rr3
               fip(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx
               fip(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fip(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkp(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkp(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkp(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz
               call dampmut (r,alphai,alphak,dmpik)
               scalek = wscale(kglob)
               rr3 = rr2 * rr1
               dmp3u = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
               dmp5u = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
               fiu(1) = dmp3u - dmp5u*xr*xr
               fiu(2) = -dmp5u*xr*yr
               fiu(3) = -dmp5u*xr*zr
               fiu(4) = -dmp5u*xr*yr
               fiu(5) = dmp3u - dmp5u*yr*yr
               fiu(6) = -dmp5u*yr*zr
               fiu(7) = -dmp5u*xr*zr
               fiu(8) = -dmp5u*yr*zr
               fiu(9) = dmp3u - dmp5u*zr*zr
            end if
c
c     increment the field at each site due to this interaction
c
            do j = 1, 3
               ef(j,1,ipoleloc) = ef(j,1,ipoleloc) + fid(j)
               ef(j,1,kbis) = ef(j,1,kbis) + fkd(j)
               ef(j,2,ipoleloc) = ef(j,2,ipoleloc) + fip(j)
               ef(j,2,kbis) = ef(j,2,kbis) + fkp(j)
            end do

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
              zmat(rofst+1+cofst1+kofst(l)) = fiu(1)
              zmat(rofst+2+cofst1+kofst(l)) = fiu(2)
              zmat(rofst+3+cofst1+kofst(l)) = fiu(3)
              zmat(rofst+1+cofst2+kofst(l)) = fiu(4)
              zmat(rofst+2+cofst2+kofst(l)) = fiu(5)
              zmat(rofst+3+cofst2+kofst(l)) = fiu(6)
              zmat(rofst+1+cofst3+kofst(l)) = fiu(7)
              zmat(rofst+2+cofst3+kofst(l)) = fiu(8)
              zmat(rofst+3+cofst3+kofst(l)) = fiu(9)
            end if

          end if
        end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(iglob)
               pscale(i12(j,iglob)) = 1.0d0
               dscale(i12(j,iglob)) = 1.0d0
            end do
            do j = 1, n13(iglob)
               pscale(i13(j,iglob)) = 1.0d0
               dscale(i13(j,iglob)) = 1.0d0
            end do
            do j = 1, n14(iglob)
               pscale(i14(j,iglob)) = 1.0d0
               dscale(i14(j,iglob)) = 1.0d0
            end do
            do j = 1, n15(iglob)
               pscale(i15(j,iglob)) = 1.0d0
               dscale(i15(j,iglob)) = 1.0d0
            end do
         else
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
         end if
         do j = 1, n12(iglob)
            wscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            wscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            wscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            wscale(i15(j,iglob)) = 1.0d0
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
      end do
c
      deallocate (dscale)
      deallocate (pscale)
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
      real(t_p) ef(3,nrhs,*)
      real(t_p), allocatable :: buffer(:,:,:,:)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
c
      allocate (buffer(3,nrhs,max(npoleloc,1),n_send1))
      allocate (req(nproc*nproc))
c      write(*,*) 'n_recep1=',n_recep1,'p_recep1=',p_recep1
c      write(*,*) 'n_send1=',n_send1,'p_send1=',p_send1
c
      do i = 1, n_send1
        tag = nproc*rank + p_send1(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,
     $   p_send1(i),tag,commloc,req(tag),ierr)
      end do
c
      do i = 1, n_recep1
        tag = nproc*p_recep1(i) + rank + 1
        call MPI_ISEND(ef(1,1,bufbegpole(p_recep1(i)+1)),
     $   3*nrhs*domlen(p_recep1(i)+1),MPI_TPREC,p_recep1(i),
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
      subroutine commfield2short(nrhs,ef)
      use atoms
      use domdec
      use mpole
      use potent
      use mpi
      implicit none
      integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc
      integer, allocatable :: req(:)
      real(t_p) ef(3,nrhs,*)
      real(t_p), allocatable :: buffer(:,:,:,:)
c
      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if

      allocate (buffer(3,nrhs,max(npoleloc,1),n_send1))
      allocate (req(nproc*nproc))
c      write(*,*) 'n_recep1=',n_recep1,'p_recep1=',p_recep1
c      write(*,*) 'n_send1=',n_send1,'p_send1=',p_send1
c
      do i = 1, n_sendshort1
        tag = nproc*rank + p_sendshort1(i) + 1
        call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npoleloc,MPI_TPREC,
     $   p_sendshort1(i),tag,commloc,req(tag),ierr)
      end do
c
      do i = 1, n_recepshort1
        tag = nproc*p_recepshort1(i) + rank + 1
        call MPI_ISEND(ef(1,1,bufbegpole(p_recepshort1(i)+1)),
     $   3*nrhs*domlen(p_recepshort1(i)+1),MPI_TPREC,p_recepshort1(i),
     $  tag,commloc,req(tag),ierr)
      end do

      call dc_factor
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

