
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     evaluate induced dipole moments and the polarization energy
!     of a group of atom without PBC
!     using a (preconditioned) conjugate gradient algorithm
!     or a JI/DIIS solver
!     useful to define the polarization to be retrieved from the total
!     one in "group-scaled" interactions
!
!     literature reference:
!     "Scalable Evaluation of Polarization Energy and Associated Forces
!     in Polarizable Molecular Dynamics: II. Toward Massively Parallel
!     Computations Using Smooth Particle Mesh Ewald",L. Lagardere et al.,
!     J. Chem. Theory Comput., 2015, 11 (6), pp 2589–2599
!
!> @brief 
!> evaluate induced dipole moments and the polarization energy
!> of a group of atom without PBC
!> using a (preconditioned) conjugate gradient algorithm
!> or a JI/DIIS solver
!> useful to define the polarization to be retrieved from the total
!> one in "group-scaled" interactions
!> @param no params
subroutine newinduce_group
   use atmlst
   use domdec
   use group
   use inform
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
!
   integer i, j, k, nrhs
   integer iglob,iglobgroup
!
!     MPI
!
   integer iipole
   integer, allocatable :: reqsend(:),reqrec(:)
   integer, allocatable :: req2send(:),req2rec(:)
!
   parameter (nrhs=2)

   real*8, allocatable :: ef(:,:,:), mu(:,:,:)
!
   external tmatxb_group
!
   if (.not.use_polar) return
!
1000 format(' illegal polalg in newinduce.')
!
   if (deb_Path) write(iout,*), 'newinduce_group '
!
!
!
!     allocate some memory and clear the arrays:
!
   allocate (mu(3,nrhs,max(1,npolegroup)))
   mu = 0d0
!
   allocate (ef(3,nrhs,max(1,npolegroup)))
   ef = 0d0
!
   allocate (reqsend(nproc))
   allocate (reqrec(nproc))
   allocate (req2send(nproc))
   allocate (req2rec(nproc))

   call orderbuffergroup
!
!     compute the electric fields:
!
   call efld0_group(nrhs,ef)
!
   call commfieldfull(nrhs,ef)
!
!     guess the dipoles:
!
!     predicted values for always stable predictor-corrector method
!
   do i = 1, npolelocgroup
      iglobgroup = ipolegroup(poleglobgroup(i))
      iglob = globglobgroup(iglobgroup)
      iipole = pollist(iglob)
      do k = 1, nrhs
         do j = 1, 3
            mu(j,k,i) = polarity(iipole)*ef(j,k,i)
         end do
      end do
   end do
!
   call commdirdirfull(nrhs,0,mu,reqrec,reqsend)
   call commdirdirfull(nrhs,1,mu,reqrec,reqsend)
   call commdirdirfull(nrhs,2,mu,reqrec,reqsend)
!

!
!     now, call the proper solver.
!
   if (polalg.eq.1) then
      call inducepcg_group(tmatxb_group,nrhs,.true.,ef,mu)
   else if (polalg.eq.2) then
      call inducejac_group(tmatxb_group,nrhs,.true.,ef,mu)
   else
      if (rank.eq.0) write(iout,1000)
      call fatal
   end if

!
!     move the computed dipoles in the module.
!
   do i = 1, npolegroup
      iglobgroup = poleglobgroup(i)
      do j = 1, 3
         uindgroup(j,iglobgroup) = mu(j,1,i)
         uinpgroup(j,iglobgroup) = mu(j,2,i)
      end do
   end do

   deallocate (reqrec)
   deallocate (reqsend)
   deallocate (req2rec)
   deallocate (req2send)
   deallocate (ef)
   deallocate (mu)
   return
end
!
!> @brief 
!> Preconditioned Conjugate Gradient Solver for "group" polarization
!> @param[in] matvec: matrix-vector product routine
!> @param[in] nrhs: number of right hand side
!> @param[in] precnd: use of diagonal preconditioner or not
!> @param[in] ef: permanent electric field (right hand side)
!> @param[in] mu: induced dipoles (guess at entry)
subroutine inducepcg_group(matvec,nrhs,precnd,ef,mu)
   use atmlst
   use domdec
   use group
   use iounit
   use math
   use mpole
   use polar
   use polpot
   use timestat
   use units
   use mpi
   implicit none
!
!     solves the polarizable dipoles linear equations by preconditioned
!     conjugate gradient. A diagonal preconditioner is used when precnd
!     is true, otherwise the standard conjugate gradient algorithm is
!     recovered by setting the preconditioner to one.
!
   integer nrhs
   real*8  ef(3,nrhs,*), mu(3,nrhs,*)
   logical precnd
   real*8, allocatable :: res(:,:,:), h(:,:,:), pp(:,:,:), zr(:,:,:),&
   &diag(:)
   integer i, it, j, k
   real*8  ggold(2), ggnew(2), gnorm(2), gg(2), alphacg(2), ene(2)
   real*8  zero, pt5, one, resnrm
   save    zero, pt5, one
   data    zero/0.0d0/, pt5/0.50d0/, one/1.0d0/
   external matvec
!
!     MPI
!
   integer iglob, iipole, ierr
   integer iglobgroup
   integer req1, req2, req3, req4
   integer status(MPI_STATUS_SIZE)
   integer,allocatable :: reqrec(:),reqsend(:)
   integer,allocatable :: req2rec(:),req2send(:)
   integer,allocatable :: reqendrec(:),reqendsend(:)
   integer,allocatable :: req2endrec(:),req2endsend(:)
!
1000 format(' cgiter shortreal converged after ',I3,' iterations.',/,&
   &' final energy        = ',2D14.7,/,&
   &' final residual norm = ',2D14.7)
1010 format(' energy and residual norm at iteration ',I3,':',4D12.2)
1020 format(' Conjugate gradient solver: induced dipoles',/,&
   &' ipole       mux         muy         muz')
1021 format(' Conjugate gradient solver: induced p-dipoles',/,&
   &' ipole       mux         muy         muz')
1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
1040 format(' Using a diagonal preconditioner.')
!
   allocate (reqrec(nproc))
   allocate (reqsend(nproc))
   allocate (req2rec(nproc))
   allocate (req2send(nproc))
   allocate (reqendsend(nproc))
   allocate (reqendrec(nproc))
   allocate (req2endsend(nproc))
   allocate (req2endrec(nproc))
!
!     allocate some memory and setup the preconditioner:
!
   allocate (res(3,nrhs,max(1,npolelocgroup)))
   allocate (h(3,nrhs,max(1,npolegroup)))
   allocate (pp(3,nrhs,max(1,npolegroup)))
   allocate (zr(3,nrhs,max(1,npolelocgroup)))
   allocate (diag(npolelocgroup))
   if (precnd) then
      do i = 1, npolelocgroup
         iglobgroup = ipolegroup(poleglobgroup(i))
         iglob = globglobgroup(iglobgroup)
         iipole = pollist(iglob)
         if (polarity(iipole).eq.0) then
            diag(i) = tinypol
         else
            diag(i) = polarity(iipole)
         end if
      end do
      if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
   else
      do i = 1, npolelocgroup
         diag(i) = one
      end do
   end if
!
!     initialize
!
   res = 0d0
   pp = 0d0
   zr = 0d0
   h = 0d0
!
!     now, compute the initial direction
!
   ggold = 0d0
!
   call matvec(nrhs,.true.,mu,h)
   call commfieldfull(nrhs,h)
!
   res(:,:,1:npolelocgroup) = ef(:,:,1:npolelocgroup)-&
   &h(:,:,1:npolelocgroup)
   do k = 1, nrhs
      do j = 1, 3
         zr(j,k,1:npolelocgroup) = diag(1:npolelocgroup)*&
         &res(j,k,1:npolelocgroup)
      end do
   end do
   pp(:,:,1:npolelocgroup) = zr(:,:,1:npolelocgroup)
   do k = 1, nrhs
      ggold(k) = sum(res(:,k,:)*zr(:,k,:))
   end do
   call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_REAL8,MPI_SUM,&
   &COMM_TINKER,req1,ierr)
!
   call commdirdirfull(nrhs,0,pp,reqrec,reqsend)
   call commdirdirfull(nrhs,1,pp,reqrec,reqsend)
   call commdirdirfull(nrhs,2,pp,reqrec,reqsend)
   call MPI_WAIT(req1,status,ierr)
!
!     now, start the main loop:
!
   do it = 1, politer
      do k = 1, nrhs
         gg(k) = zero
      end do
!
      call matvec(nrhs,.true.,pp,h)
      call commfieldfull(nrhs,h)
!
      do k = 1, nrhs
         gg(k) = sum(pp(:,k,1:npolelocgroup)*h(:,k,1:npolelocgroup))
      end do
      call MPI_IALLREDUCE(MPI_IN_PLACE,gg(1),nrhs,MPI_REAL8,MPI_SUM,&
      &COMM_TINKER,req2,ierr)
      call MPI_WAIT(req2,status,ierr)
      do k = 1, nrhs
         if (gg(k).eq.zero) return
         alphacg(k)  = ggold(k)/gg(k)
         ggnew(k)  = zero
         ene(k)    = zero
      end do
      do k = 1, nrhs
         mu(:,k,1:npolelocgroup) = mu(:,k,1:npolelocgroup)+ alphacg(k)*&
         &pp(:,k,1:npolelocgroup)
         res(:,k,1:npolelocgroup) =res(:,k,1:npolelocgroup)-alphacg(k)*&
         &h(:,k,1:npolelocgroup)
      end do
      do k = 1, nrhs
         do j = 1, 3
            zr(j,k,1:npolelocgroup) = diag(1:npolelocgroup)&
            &*res(j,k,1:npolelocgroup)
         end do
      end do
      do k = 1, nrhs
         ggnew(k) = sum(res(:,k,:)*zr(:,k,:))
         ene(k) = -pt5*sum(mu(:,k,1:npolelocgroup)*&
         &(res(:,k,1:npolelocgroup)+ ef(:,k,1:npolelocgroup)))
      end do
      call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_REAL8,&
      &MPI_SUM,COMM_TINKER,req3,ierr)
      call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_REAL8,MPI_SUM,&
      &COMM_TINKER,req4,ierr)
      call MPI_WAIT(req3,status,ierr)
      call MPI_WAIT(req4,status,ierr)
      resnrm = zero
      do k = 1, nrhs
         gnorm(k) = sqrt(ggnew(k)/dble(3*npolegroup))
         resnrm   = max(resnrm,gnorm(k))
      end do
      if (polprt.ge.2.and.rank.eq.0) write(iout,1010)&
      &it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)
      do k = 1, nrhs
         pp(:,k,1:npolelocgroup) = zr(:,k,1:npolelocgroup)&
         &+ggnew(k)/ggold(k)*pp(:,k,1:npolelocgroup)
      end do
!
      call commdirdirfull(nrhs,0,pp,reqrec,reqsend)
      call commdirdirfull(nrhs,1,pp,reqrec,reqsend)
      call commdirdirfull(nrhs,2,pp,reqrec,reqsend)
!
      ggold = ggnew
      if (resnrm.lt.poleps) then
         if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,&
         &(ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
         goto 10
      end if
   end do
10 continue
!
   call commdirdirfull(nrhs,0,mu,reqendrec,reqendsend)
   call commdirdirfull(nrhs,1,mu,reqendrec,reqendsend)
   call commdirdirfull(nrhs,2,mu,reqendrec,reqendsend)
!
!     finalize and debug printing:
!
   if (polprt.ge.3) then
      write(iout,1020)
      do i = 1, npolelocgroup
         iglobgroup = ipolegroup(poleglobgroup(i))
         iglob = globglobgroup(iglobgroup)
         write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
      end do
   end if
   if (polprt.ge.4) then
      write(iout,1021)
      do i = 1, npolelocgroup
         iglobgroup = ipolegroup(poleglobgroup(i))
         iglob = globglobgroup(iglobgroup)
         write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
      end do
   end if
   deallocate (reqsend)
   deallocate (reqrec)
   deallocate (req2send)
   deallocate (req2rec)
   deallocate (reqendsend)
   deallocate (reqendrec)
   deallocate (req2endsend)
   deallocate (req2endrec)
   deallocate (res)
   deallocate (h)
   deallocate (pp)
   deallocate (zr)
   deallocate (diag)
   return
end
!
!> @brief 
!> Jacobi/DIIS Solver for "group" polarization
!> @param[in] matvec: matrix-vector product routine
!> @param[in] nrhs: number of right hand side
!> @param[in] dodiis: use of DIIS extrapolation or not
!> @param[in] ef: permanent electric field (right hand side)
!> @param[in] mu: induced dipoles (guess at entry)
subroutine inducejac_group(matvec,nrhs,dodiis,ef,mu)
   use atmlst
   use domdec
   use group
   use iounit
   use math
   use mpole
   use polar
   use polpot
   use potent
   use timestat
   use mpi
   implicit none
!
!     solves the polarizable dipoles linear equations by jacobi
!     iterations coupled to the DIIS extrapolation.
!
   integer nrhs, info
   real*8  ef(3,nrhs,*), mu(3,nrhs,*)
   real*8, allocatable :: munew(:,:,:), h(:,:,:)
   real*8, allocatable :: xdiis(:,:), ediis(:,:), bmat(:,:),&
   &bloc(:,:), cex(:)
   integer i, j, k, it, ind, ndismx, nmat, lenb
   logical dodiis
   parameter (ndismx=25)
   integer ipiv(ndismx+1)
   real*8  zero, one, rnorm(2), rr, xx(1)

   save    zero, one, xx
   external matvec
!
!     MPI
!
   integer iglob, iipole, ierr
   integer iglobgroup
   integer reqnorm,reqdiis(2*ndismx+1)
   integer status(MPI_STATUS_SIZE)
   integer,allocatable :: reqrec(:),reqsend(:)
   integer,allocatable :: req2rec(:),req2send(:)
!
1000 format(' itsolv shortreal converged after ',I3,' iterations.',/,&
   &' final residual norm = ',3D14.7)
1010 format(' residual norm at iteration ',I3,':',3D12.2)
1020 format(' Jacobi/DIIS solver: induced dipoles',/,&
   &' ipole       mux         muy         muz')
1021 format(' Jacobi/DIIS solver: induced p-dipoles',/,&
   &' ipole       mux         muy         muz')
1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
!
!
   zero  = 0.0d0
   one   = 1.0d0
   xx(1) = 0.0d0

   allocate (munew(3,nrhs,max(1,npolegroup)))
   munew = 0d0
   allocate (reqrec(nproc))
   allocate (reqsend(nproc))
   allocate (req2rec(nproc))
   allocate (req2send(nproc))
!
   allocate (h(3,nrhs,max(1,npolegroup)))
   h = 0d0
   if (dodiis) then
      nmat = 1
      lenb = ndismx + 1
      allocate (xdiis(3*nrhs*max(npolelocgroup,1),ndismx))
      allocate (ediis(3*nrhs*max(npolelocgroup,1),ndismx))
      allocate (bmat(lenb,lenb))
      allocate (bloc(lenb,lenb))
      allocate (cex(lenb))
      bmat = 0d0
   end if
!
!     main loop:
!
   do it = 1, politer
      h = 0d0
      rnorm = 0d0
!
      call matvec(nrhs,.false.,mu,h)
      call commfieldfull(nrhs,h)
!
!     jacobi step:
!
!
      do i = 1, npolelocgroup
         iglobgroup = ipolegroup(poleglobgroup(i))
         iglob = globglobgroup(iglobgroup)
         iipole = pollist(iglob)
         do k = 1, nrhs
            do j = 1, 3
               munew(j,k,i) = polarity(iipole)*(ef(j,k,i) - h(j,k,i))
            end do
         end do
      end do
!
      do k = 1, nrhs
         rnorm(k) = sum((munew(1:3,k,1:npolelocgroup)-&
         &mu(1:3,k,1:npolelocgroup))**2)
      end do
!
      call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_REAL8,&
      &MPI_SUM,COMM_TINKER,reqnorm,ierr)
      if (dodiis) then
         ind = 0
         do i = 1, npolelocgroup
            do k = 1, nrhs
               do j = 1, 3
                  ind = ind + 1
                  xdiis(ind,nmat) = munew(j,k,i)
                  ediis(ind,nmat) = munew(j,k,i) - mu(j,k,i)
               end do
            end do
         end do
!
!         Compute Pulay's Matrix and extrapolate
!
         call diis(ndismx,3*nrhs*npolelocgroup,xdiis,ediis,&
         &bmat,nmat,reqdiis,COMM_TINKER)
!
         do i = 1, 2*nmat-3
            call MPI_WAIT(reqdiis(i),status,ierr)
         end do
         bloc = bmat
         cex = 0d0
         cex(1) = one
         call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
         munew = 0d0
         call extrap(3*nrhs*npolelocgroup,nmat-1,xdiis,cex,munew)
      end if
!
      call commdirdirfull(nrhs,0,mu,reqrec,reqsend)
      call commdirdirfull(nrhs,1,munew,reqrec,reqsend)
      call commdirdirfull(nrhs,2,mu,reqrec,reqsend)
      mu(:,:,1:npolelocgroup) = munew(:,:,1:npolelocgroup)

!
!     compute the norm of the increment.
!
      call MPI_WAIT(reqnorm,status,ierr)
      rr = zero
      do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/dble(3*npolegroup))
         rr = max(rnorm(k),rr)
      end do
      if (rr.lt.poleps) then
         if (polprt.ge.1.and.rank.eq.0)&
         &write(6,1000) it, (rnorm(k), k = 1, nrhs)
         goto 10
      end if
      if (polprt.ge.1.and.rank.eq.0)&
      &write(6,1010) it, (rnorm(k), k = 1, nrhs)
   end do
10 continue
   if (polprt.ge.3) then
      write(iout,1020)
      do i = 1, npolelocgroup
         iglobgroup = ipolegroup(poleglobgroup(i))
         iglob = globglobgroup(iglobgroup)
         write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
      end do
   end if
   if (polprt.ge.4) then
      write(iout,1021)
      do i = 1, npolelocgroup
         iglobgroup = ipolegroup(poleglobgroup(i))
         iglob = globglobgroup(iglobgroup)
         write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
      end do
   end if
!
!     free the memory.
!
   deallocate (reqsend)
   deallocate (reqrec)
   deallocate (req2send)
   deallocate (req2rec)
   deallocate (munew)
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
!
!> @brief 
!> Compute the direct space contribution to the permanent electric field.
!>  Also compute the "p" field, which is used to
!> compute the energy according to the AMOEBA force field.
!> @param[in] nrhs: number of right hand side
!> @param[out] ef: electric field
subroutine efld0_group(nrhs,ef)
!
!     Compute the direct space contribution to the permanent electric field.
!      Also compute the "p" field, which is used to
!     compute the energy according to the AMOEBA force field.
!
   use atmlst
   use atoms
   use bound
   use chgpen
   use couple
   use cutoff
   use domdec
   use ewald
   use group
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
   integer iglob,kglob,nrhs
   integer iglobgroup,kglobgroup
   real*8  ef(3,nrhs,npolegroup)
   integer ii, j, k,kk, kkk, iipole, kkpole
   real*8 xr,yr,zr
   real*8 r,r2,rr1,rr2
   real*8 rr3,rr5,rr7
   real*8 rr3i,rr5i,rr7i
   real*8 rr3k,rr5k,rr7k
   real*8 ci,dix,diy,diz
   real*8 qixx,qiyy,qizz
   real*8 qixy,qixz,qiyz
   real*8 ck,dkx,dky,dkz
   real*8 qkxx,qkyy,qkzz
   real*8 qkxy,qkxz,qkyz
   real*8 dir,dkr
   real*8 qix,qiy,qiz,qir
   real*8 qkx,qky,qkz,qkr
   real*8 fid(3), fip(3)
   real*8 fkd(3), fkp(3)
   real*8 corei,corek
   real*8 vali,valk
   real*8 alphai,alphak
   real*8 dmp3,dmp5,dmp7
   real*8 dmpi(7),dmpk(7)
   real*8 dmpik(7)
   real*8 scalek
   real*8, allocatable :: dscale(:)
   real*8, allocatable :: pscale(:)
!
   allocate (dscale(n))
   allocate (pscale(n))
   dscale = 1.0d0
   pscale = 1.0d0

!
   do ii = 1, npolelocgroup
      iglobgroup = ipolegroup(poleglobgroup(ii))
      iglob = globglobgroup(iglobgroup)
      iipole = pollist(iglob)
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
!
!     set exclusion coefficients for connected atoms
!
      if (dpequal) then
         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = p2scale
            do k = 1, np11(iglob)
               if (i12(j,iglob) .eq. ip11(k,iglob))&
               &pscale(i12(j,iglob)) = p2iscale
            end do
            dscale(i12(j,iglob)) = pscale(i12(j,iglob))
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = p3scale
            do k = 1, np11(iglob)
               if (i13(j,iglob) .eq. ip11(k,iglob))&
               &pscale(i13(j,iglob)) = p3iscale
            end do
            dscale(i13(j,iglob)) = pscale(i13(j,iglob))
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = p4scale
            do k = 1, np11(iglob)
               if (i14(j,iglob) .eq. ip11(k,iglob))&
               &pscale(i14(j,iglob)) = p4iscale
            end do
            dscale(i14(j,iglob)) = pscale(i14(j,iglob))
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = p5scale
            do k = 1, np11(iglob)
               if (i15(j,iglob) .eq. ip11(k,iglob))&
               &pscale(i15(j,iglob)) = p5iscale
            end do
            dscale(i15(j,iglob)) = pscale(i15(j,iglob))
         end do
      else
         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = p2scale
            do k = 1, np11(iglob)
               if (i12(j,iglob) .eq. ip11(k,iglob))&
               &pscale(i12(j,iglob)) = p2iscale
            end do
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = p3scale
            do k = 1, np11(iglob)
               if (i13(j,iglob) .eq. ip11(k,iglob))&
               &pscale(i13(j,iglob)) = p3iscale
            end do
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = p4scale
            do k = 1, np11(iglob)
               if (i14(j,iglob) .eq. ip11(k,iglob))&
               &pscale(i14(j,iglob)) = p4iscale
            end do
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = p5scale
            do k = 1, np11(iglob)
               if (i15(j,iglob) .eq. ip11(k,iglob))&
               &pscale(i15(j,iglob)) = p5iscale
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
!
      do kkk = iglobgroup+1, npolegroup
         kglobgroup = ipolegroup(kkk)
         kglob = globglobgroup(kglobgroup)
         kkpole = pollist(kglob)
         kk = polelocgroup(kkk)
         xr = x(kglob) - x(iglob)
         yr = y(kglob) - y(iglob)
         zr = z(kglob) - z(iglob)
         r2 = xr*xr + yr* yr + zr*zr
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
!
!     intermediates involving moments and separation distance
!
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
!
!     find the field components for Thole polarization damping
!
         if (use_thole) then
            call damptholed (iipole,kkpole,7,r,dmpik)
            scalek = dscale(kglob)
            dmp3 = scalek*dmpik(3)*rr3
            dmp5 = scalek*dmpik(5)*rr5
            dmp7 = scalek*dmpik(7)*rr7
            fid(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)&
            &- dmp3*dkx + 2.0d0*dmp5*qkx
            fid(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)&
            &- dmp3*dky + 2.0d0*dmp5*qky
            fid(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)&
            &- dmp3*dkz + 2.0d0*dmp5*qkz
            fkd(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)&
            &- dmp3*dix - 2.0d0*dmp5*qix
            fkd(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)&
            &- dmp3*diy - 2.0d0*dmp5*qiy
            fkd(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)&
            &- dmp3*diz - 2.0d0*dmp5*qiz
            scalek = pscale(kglob)
            dmp3 = scalek*dmpik(3)*rr3
            dmp5 = scalek*dmpik(5)*rr5
            dmp7 = scalek*dmpik(7)*rr7
            fip(1) = -xr*(dmp3*ck-dmp5*dkr+dmp7*qkr)&
            &- dmp3*dkx + 2.0d0*dmp5*qkx
            fip(2) = -yr*(dmp3*ck-dmp5*dkr+dmp7*qkr)&
            &- dmp3*dky + 2.0d0*dmp5*qky
            fip(3) = -zr*(dmp3*ck-dmp5*dkr+dmp7*qkr)&
            &- dmp3*dkz + 2.0d0*dmp5*qkz
            fkp(1) = xr*(dmp3*ci+dmp5*dir+dmp7*qir)&
            &- dmp3*dix - 2.0d0*dmp5*qix
            fkp(2) = yr*(dmp3*ci+dmp5*dir+dmp7*qir)&
            &- dmp3*diy - 2.0d0*dmp5*qiy
            fkp(3) = zr*(dmp3*ci+dmp5*dir+dmp7*qir)&
            &- dmp3*diz - 2.0d0*dmp5*qiz
!
!     find the field components for charge penetration damping
!
         else if (use_chgpen) then
            corek = pcore(kkpole)
            valk = pval(kkpole)
            alphak = palpha(kkpole)
            call dampdir (r,alphai,alphak,dmpi,dmpk)
            scalek = dscale(kglob)
            rr3i = scalek*dmpi(3)*rr3
            rr5i = scalek*dmpi(5)*rr5
            rr7i = scalek*dmpi(7)*rr7
            rr3k = scalek*dmpk(3)*rr3
            rr5k = scalek*dmpk(5)*rr5
            rr7k = scalek*dmpk(7)*rr7
            rr3 = scalek*rr3
            fid(1) = -xr*(rr3*corek + rr3k*valk&
            &- rr5k*dkr + rr7k*qkr)&
            &- rr3k*dkx + 2.0d0*rr5k*qkx
            fid(2) = -yr*(rr3*corek + rr3k*valk&
            &- rr5k*dkr + rr7k*qkr)&
            &- rr3k*dky + 2.0d0*rr5k*qky
            fid(3) = -zr*(rr3*corek + rr3k*valk&
            &- rr5k*dkr + rr7k*qkr)&
            &- rr3k*dkz + 2.0d0*rr5k*qkz
            fkd(1) = xr*(rr3*corei + rr3i*vali&
            &+ rr5i*dir + rr7i*qir)&
            &- rr3i*dix - 2.0d0*rr5i*qix
            fkd(2) = yr*(rr3*corei + rr3i*vali&
            &+ rr5i*dir + rr7i*qir)&
            &- rr3i*diy - 2.0d0*rr5i*qiy
            fkd(3) = zr*(rr3*corei + rr3i*vali&
            &+ rr5i*dir + rr7i*qir)&
            &- rr3i*diz - 2.0d0*rr5i*qiz
            scalek = pscale(kglob)
            rr3 = rr2 * rr1
            rr3i = scalek*dmpi(3)*rr3
            rr5i = scalek*dmpi(5)*rr5
            rr7i = scalek*dmpi(7)*rr7
            rr3k = scalek*dmpk(3)*rr3
            rr5k = scalek*dmpk(5)*rr5
            rr7k = scalek*dmpk(7)*rr7
            rr3 = scalek*rr3
            fip(1) = -xr*(rr3*corek + rr3k*valk&
            &- rr5k*dkr + rr7k*qkr)&
            &- rr3k*dkx + 2.0d0*rr5k*qkx
            fip(2) = -yr*(rr3*corek + rr3k*valk&
            &- rr5k*dkr + rr7k*qkr)&
            &- rr3k*dky + 2.0d0*rr5k*qky
            fip(3) = -zr*(rr3*corek + rr3k*valk&
            &- rr5k*dkr + rr7k*qkr)&
            &- rr3k*dkz + 2.0d0*rr5k*qkz
            fkp(1) = xr*(rr3*corei + rr3i*vali&
            &+ rr5i*dir + rr7i*qir)&
            &- rr3i*dix - 2.0d0*rr5i*qix
            fkp(2) = yr*(rr3*corei + rr3i*vali&
            &+ rr5i*dir + rr7i*qir)&
            &- rr3i*diy - 2.0d0*rr5i*qiy
            fkp(3) = zr*(rr3*corei + rr3i*vali&
            &+ rr5i*dir + rr7i*qir)&
            &- rr3i*diz - 2.0d0*rr5i*qiz
         end if
!
!     increment the field at each site due to this interaction
!
         do j = 1, 3
            ef(j,1,ii) = ef(j,1,ii) + fid(j)
            ef(j,1,kk) = ef(j,1,kk) + fkd(j)
            ef(j,2,ii) = ef(j,2,ii) + fip(j)
            ef(j,2,kk) = ef(j,2,kk) + fkp(j)
         end do
      end do
!
!     reset exclusion coefficients for connected atoms
!
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
   end do
!
   deallocate (dscale)
   deallocate (pscale)
   return
end

!> @brief 
!> Compute the direct space contribution to the electric field due to the current value
!> of the induced dipoles
!> @param[in] nrhs: number of right hand side
!> @param[in] dodiag: compute diagonal contribution of matrix-vector product or not
!> @param[in] mu: induced dipoles
!> @param[out] efi: electrif field due to induced dipoles
subroutine tmatxb_group(nrhs,dodiag,mu,efi)
!
!     Compute the direct space contribution to the electric field due to the current value
!     of the induced dipoles
!
   use atmlst
   use atoms
   use bound
   use chgpen
   use couple
   use domdec
   use ewald
   use group
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
   integer i,nrhs,iglob,kglob,iipole,kkpole
   integer kk,iglobgroup,kglobgroup
   real*8  mu(3,nrhs,npolegroup), efi(3,nrhs,npolegroup)
   real*8 xr,yr,zr
   real*8 r,r2,rr1,rr2
   real*8 rr3,rr5
   real*8 rr3ik,rr5ik
   real*8 scalek
   real*8 dmp3,dmp5
   real*8 fid(3),fkd(3)
   real*8 fip(3),fkp(3)
   real*8 dmpik(7)
   real*8 dlocal(6)
   real*8 duix,duiy,duiz
   real*8 puix,puiy,puiz
   real*8 dukx,duky,dukz
   real*8 pukx,puky,pukz
   real*8 alphai,alphak
   real*8 pol
   real*8, allocatable :: uscale(:)
   real*8, allocatable :: wscale(:)
   logical dodiag
   integer j, ii, kkk, irhs
!
!     initialize the result vector
!
   efi = 0d0
!
!     perform dynamic allocation of some local arrays
!
   allocate (uscale(n))
   allocate (wscale(n))
!
!     set arrays needed to scale connected atom interactions
!
   uscale = 1.0d0
   wscale = 1.0d0
!
   do ii = 1, npolelocgroup
      iglobgroup = ipolegroup(poleglobgroup(ii))
      iglob = globglobgroup(iglobgroup)
      iipole = pollist(iglob)
      duix = mu(1,1,ii)
      duiy = mu(2,1,ii)
      duiz = mu(3,1,ii)
      puix = mu(1,2,ii)
      puiy = mu(2,2,ii)
      puiz = mu(3,2,ii)
      if (use_chgpen) then
         alphai = palpha(iipole)
      end if
!
!
!     set exclusion coefficients for connected atoms
!
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
      do kkk = iglobgroup+1, npolegroup
         kglobgroup = ipolegroup(kkk)
         kglob = globglobgroup(kglobgroup)
         kkpole = pollist(kglob)
         kk = polelocgroup(kkk)
         xr = x(kglob) - x(iglob)
         yr = y(kglob) - y(iglob)
         zr = z(kglob) - z(iglob)
         r2 = xr*xr + yr* yr + zr*zr
!
!     compute the distances and the scaling factors according to
!     Thole's model.
!
         r = sqrt(r2)
         rr1 = 1.0d0 / r
         rr2 = rr1 * rr1
         rr3 = rr2 * rr1
         rr5 = 3.0d0 * rr2 * rr3
         dukx = mu(1,1,kk)
         duky = mu(2,1,kk)
         dukz = mu(3,1,kk)
         pukx = mu(1,2,kk)
         puky = mu(2,2,kk)
         pukz = mu(3,2,kk)
         if (use_chgpen) then
            alphak = palpha(kkpole)
         end if
!
!     find the field components for Thole polarization damping
!
         if (use_thole) then
            call dampthole (iipole,kkpole,5,r,dmpik)
            scalek = uscale(kglob)
            dmp3 = -scalek*dmpik(3)*rr3
            dmp5 = -scalek*dmpik(5)*rr5
            dlocal(1) = -dmp3 + dmp5*xr*xr
            dlocal(2) = dmp5*xr*yr
            dlocal(3) = dmp5*xr*zr
            dlocal(4) = -dmp3 + dmp5*yr*yr
            dlocal(5) = dmp5*yr*zr
            dlocal(6) = -dmp3 + dmp5*zr*zr
!
!     find the field components for charge penetration damping
!
         else if (use_chgpen) then
            call dampmut (r,alphai,alphak,dmpik)
            scalek = wscale(kglob)
            rr3ik = -scalek*dmpik(3)*rr3
            rr5ik = -scalek*dmpik(5)*rr5
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
            efi(j,1,ii) = efi(j,1,ii) + fid(j)
            efi(j,1,kk) = efi(j,1,kk) + fkd(j)
            efi(j,2,ii) = efi(j,2,ii) + fip(j)
            efi(j,2,kk) = efi(j,2,kk) + fkp(j)
         end do
      end do
!
!     reset interaction scaling coefficients for connected atoms
!
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
!
!     if dodiag is true, also compute the "self-induced" field,
!     i.e., the diagonal portion of the matrix/vector product.
!
      do i = 1, npolelocgroup
         iglobgroup = ipolegroup(poleglobgroup(i))
         iglob = globglobgroup(iglobgroup)
         iipole = pollist(iglob)
!
!     if no polarisability, take a negligeable value to allow convergence
!
         if (polarity(iipole).eq.0.0d0) then
            pol = tinypol ** -1
         else
            pol  = polarity(iipole) ** -1
         endif
         do irhs = 1, nrhs
            do j = 1, 3
               efi(j,irhs,i) = efi(j,irhs,i) +&
               &mu(j,irhs,i)*pol
            end do
         end do
      end do
   end if
!
!     perform deallocation of some local arrays
!
   deallocate (uscale)
   deallocate (wscale)
   return
end

!
!    subroutine orderbuffergroup : get sizes and indexes of all the domains
!
!> @brief 
!> get sizes and indexes of all the domains
!> @param no params
subroutine orderbuffergroup
   use atoms
   use domdec
   use group
   use mpole
   use potent
   use mpi
   implicit none
   integer i,iproc,tag,status(MPI_STATUS_SIZE)
   integer ierr,commloc
   integer count
   integer, allocatable :: reqrec(:), reqsend(:)
!
   if (use_pmecore) then
      commloc = comm_dir
   else
      commloc = COMM_TINKER
   end if
   allocate (reqrec(nproc))
   allocate (reqsend(nproc))
!
!     deal with atoms first
!
   domlengroup(rank+1) = nlocatgroup

   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         tag = nproc*rank + iproc + 1
         call MPI_IRECV(domlengroup(iproc+1),1,MPI_INT,iproc,tag,&
         &COMM_TINKER,reqrec(iproc+1),ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         tag = nproc*iproc + rank + 1
         call MPI_ISEND(domlengroup(rank+1),1,MPI_INT,iproc,tag,&
         &COMM_TINKER,reqsend(iproc+1),ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         call MPI_WAIT(reqsend(iproc+1),status,ierr)
         call MPI_WAIT(reqrec(iproc+1),status,ierr)
      end if
   end do
!
   bufbeggroup(rank+1) = 1
   count = domlengroup(rank+1)
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         if (domlengroup(iproc+1).eq.0) then
            bufbeggroup(iproc+1) = 1
         else
            bufbeggroup(iproc+1) = count + 1
         end if
         count = count + domlengroup(iproc+1)
      end if
   end do
!
!     get the indexes of the atoms in the group frame
!
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         tag = nproc*rank + iproc + 1
         call MPI_IRECV(globgroup(bufbeggroup(iproc+1)),&
         &domlengroup(iproc+1),MPI_INT,iproc,tag,&
         &COMM_TINKER,reqrec(iproc+1),ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         tag = nproc*iproc + rank + 1
         call MPI_ISEND(globgroup,domlengroup(rank+1),MPI_INT,&
         &iproc,tag,COMM_TINKER,reqsend(iproc+1),ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         call MPI_WAIT(reqrec(iproc+1),status,ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         call MPI_WAIT(reqsend(iproc+1),status,ierr)
      end if
   end do

   do iproc = 1, nproc
      if (iproc.ne.(rank+1)) then
         do i = 1, domlengroup(iproc)
            locgroup(globgroup(bufbeggroup(iproc)+i-1)) =&
            &bufbeggroup(iproc)+i-1
         end do
      end if
   end do
!
!     then deal with multipoles
!
   domlenpolegroup(rank+1) = npolelocgroup
!
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         tag = nproc*rank + iproc + 1
         call MPI_IRECV(domlenpolegroup(iproc+1),1,MPI_INT,iproc,tag,&
         &COMM_TINKER,reqrec(iproc+1),ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         tag = nproc*iproc + rank + 1
         call MPI_ISEND(domlenpolegroup(rank+1),1,MPI_INT,iproc,tag,&
         &COMM_TINKER,reqsend(iproc+1),ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         call MPI_WAIT(reqsend(iproc+1),status,ierr)
         call MPI_WAIT(reqrec(iproc+1),status,ierr)
      end if
   end do
!
   bufbegpolegroup(rank+1) = 1
   count = domlenpolegroup(rank+1)
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         if (domlenpolegroup(iproc+1).eq.0) then
            bufbegpolegroup(iproc+1) = 1
         else
            bufbegpolegroup(iproc+1) = count + 1
         end if
         count = count + domlenpolegroup(iproc+1)
      end if
   end do
!
!     get the indexes of the mulipoles in the group frame
!
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         tag = nproc*rank + iproc + 1
         call MPI_IRECV(poleglobgroup(bufbegpolegroup(iproc+1)),&
         &domlenpolegroup(iproc+1),MPI_INT,iproc,tag,&
         &COMM_TINKER,reqrec(iproc+1),ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         tag = nproc*iproc + rank + 1
         call MPI_ISEND(poleglobgroup,domlenpolegroup(rank+1),MPI_INT,&
         &iproc,tag,COMM_TINKER,reqsend(iproc+1),ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         call MPI_WAIT(reqrec(iproc+1),status,ierr)
      end if
   end do
   do iproc = 0, nproc-1
      if (iproc.ne.rank) then
         call MPI_WAIT(reqsend(iproc+1),status,ierr)
      end if
   end do

   do iproc = 1, nproc
      if (iproc.ne.(rank+1)) then
         do i = 1, domlenpolegroup(iproc)
            polelocgroup(poleglobgroup(bufbegpolegroup(iproc)+i-1)) =&
            &bufbegpolegroup(iproc)+i-1
         end do
      end if
   end do
!
   deallocate (reqrec)
   deallocate (reqsend)
   return
end
!
!     subroutine commfieldfull : communicate some direct fields (Newton's third law)
!
!> @brief 
!> communicate some direct fields (Newton's third law)
!> @param[in] nrhs: number of right hand side
!> @param[in] ef: electric field
subroutine commfieldfull(nrhs,ef)
   use atoms
   use domdec
   use group
   use mpole
   use potent
   use mpi
   implicit none
   integer nrhs,i,j,k,tag,status(MPI_STATUS_SIZE)
   integer ierr,commloc
   integer, allocatable :: reqrec(:), reqsend(:)
   real*8 ef(3,nrhs,*)
   real*8, allocatable :: buffer(:,:,:,:)
!
   if (use_pmecore) then
      commloc = comm_dir
   else
      commloc = COMM_TINKER
   end if
!
   allocate (buffer(3,nrhs,max(npolelocgroup,1),nproc))
   allocate (reqrec(nproc))
   allocate (reqsend(nproc))
!
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         tag = nproc*rank + i
         call MPI_IRECV(buffer(1,1,1,i),3*nrhs*npolelocgroup,MPI_REAL8,&
         &i-1,tag,commloc,reqrec(i),ierr)
      end if
   end do
!
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         tag = nproc*(i-1) + rank + 1
         call MPI_ISEND(ef(1,1,bufbegpolegroup(i)),&
         &3*nrhs*domlenpolegroup(i),MPI_REAL8,i-1,&
         &tag,commloc,reqsend(i),ierr)
      end if
   end do
!
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         tag = nproc*rank + i
         call MPI_WAIT(reqrec(i),status,ierr)
      end if
   end do
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         tag = nproc*(i-1) + rank + 1
         call MPI_WAIT(reqsend(i),status,ierr)
      end if
   end do
!
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         do j = 1, npolelocgroup
            do k = 1, nrhs
               ef(1,k,j) = ef(1,k,j) + buffer(1,k,j,i)
               ef(2,k,j) = ef(2,k,j) + buffer(2,k,j,i)
               ef(3,k,j) = ef(3,k,j) + buffer(3,k,j,i)
            end do
         end do
      end if
   end do
!
   deallocate (buffer)
   deallocate (reqrec)
   deallocate (reqsend)
   return
end
!
!    subroutine commdirdirfull : deal with communications of direct dipoles
!
!    rule determines what to do:
!        - 0: start reception of direct dipoles
!        - 1: start sending of direct dipoles
!        - 2: wait for the communications to be done
!
!> @brief 
!> deal with communications of direct dipoles
!> @param[in] nrhs: number of right hand side
!> @param[in] rule: integer determining what to do:
!>        - 0: start reception of direct dipoles
!>        - 1: start sending of direct dipoles
!>        - 2: wait for the communications to be done
!> @param[in] mu: induced dipoles
!> @param[in] reqrec: array of MPI request for reception
!> @param[in] reqsend: array of MPI request for sending
subroutine commdirdirfull(nrhs,rule,mu,reqrec,reqsend)
   use domdec
   use group
   use iounit
   use mpole
   use mpi
   implicit none
   integer nrhs,rule,ierr,status(MPI_STATUS_SIZE),tag,i
   integer :: reqrec(nproc),reqsend(nproc)
   real*8 mu(3,nrhs,npolegroup)
1000 format(' illegal rule in commdirdirfull.')

   if (rule.eq.0) then
!
!     MPI : begin reception
!
      do i = 1, nproc
         if (i.ne.(rank+1)) then
            tag = nproc*rank + i
            call MPI_IRECV(mu(1,1,bufbegpolegroup(i)),&
            &3*nrhs*domlenpolegroup(i),MPI_REAL8,i-1,tag,&
            &COMM_TINKER,reqrec(i),ierr)
         end if
      end do

   else if (rule.eq.1) then
!
      do i = 1, nproc
         if (i.ne.(rank+1)) then
            tag = nproc*(i-1) + rank + 1
            call MPI_ISEND(mu,3*nrhs*npolelocgroup,&
            &MPI_REAL8,i-1,tag,COMM_TINKER,&
            &reqsend(i),ierr)
         end if
      end do
   else if (rule.eq.2) then
      do i = 1, nproc
         if (i.ne.(rank+1)) then
            call MPI_WAIT(reqrec(i),status,ierr)
         end if
      end do
      do i = 1, nproc
         if (i.ne.(rank+1)) then
            call MPI_WAIT(reqsend(i),status,ierr)
         end if
      end do
   else
      if (rank.eq.0) write(iout,1000)
      call fatal
   end if

   return
end
!
!     subroutine commforcesgroup : deal with communications of "polarization group" forces
!
!> @brief 
!> deal with communications of "polarization group" forces
!> @param no params
subroutine commforcesgroup
   use sizes
   use atoms
   use deriv
   use domdec
   use group
   use potent
   use mpi
   implicit none
   integer i,tag,ierr
   integer, allocatable :: reqsend(:),reqrec(:)
   real*8, allocatable :: buffer(:,:,:)
   integer status(MPI_STATUS_SIZE)
!
!     allocate some arrays
!
   allocate (reqsend(nproc))
   allocate (reqrec(nproc))
   allocate (buffer(3,max(1,nlocatgroup),nproc))
!      buffer = 0d0
!      buffers = 0d0
!
!     MPI : begin reception in buffer
!
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         tag = nproc*rank + i
         call MPI_IRECV(buffer(1,1,i),3*nlocatgroup,&
         &MPI_REAL8,i-1,tag,COMM_TINKER,reqrec(i),ierr)
      end if
   end do
!
!     MPI : begin sending
!
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         tag = nproc*(i-1) + rank + 1
         call MPI_ISEND(depgroup(1,bufbeggroup(i)),&
         &3*domlengroup(i),&
         &MPI_REAL8,i-1,tag,COMM_TINKER,&
         &reqsend(i),ierr)
      end if
   end do
!
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         call MPI_WAIT(reqrec(i),status,ierr)
      end if
   end do
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         call MPI_WAIT(reqsend(i),status,ierr)
      end if
   end do
!
!     MPI : move in global arrays
!
   do i = 1, nproc
      if (i.ne.(rank+1)) then
         depgroup(:,1:nlocatgroup) = depgroup(:,1:nlocatgroup) +&
         &buffer(:,1:nlocatgroup,i)
      end if
   end do
   deallocate (reqsend)
   deallocate (reqrec)
   deallocate (buffer)
   return
end
