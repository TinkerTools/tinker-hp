c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     subroutine emrecip: computes reciprocal space part of multipole and dipole
c     polarizability interaction using PME
c
      subroutine emrecip
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'atmlst.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'fft.i'
      include 'math.i'
      include 'mpole.i'
      include 'pme.i'
      include 'polar.i'
      include 'potent.i'
      include 'openmp.i'
      include 'mpif.h'
      integer ierr,iipole,proc
      integer status(MPI_STATUS_SIZE),tag
      integer iglob
      integer i,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nprocloc
      real*8 e,r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 a(3,3)
      real*8 fuind(3)
      real*8 cmp(10)
      real*8, allocatable ::  fmp(:,:)
      integer, allocatable :: req(:),reqbcast(:)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
c
c     dynamic allocation of local arrays
c
      allocate (fmp(20,npolerecloc))
      if (use_pmecore) then
        nprocloc = nrec
      else
        nprocloc = nproc
      end if
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (req(nprocloc*nprocloc))
      allocate (reqbcast(nprocloc*nprocloc))
c
      if (.not.use_polar) then
        if (associated(cphirec)) deallocate (cphirec)
        allocate (cphirec(10,max(npolerecloc,1)))
        cphirec = 0d0
        if (associated(fphirec)) deallocate(fphirec)
        allocate (fphirec(20,max(npolerecloc,1)))
        fphirec = 0d0
      end if
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
      if (.not.use_polar) then
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          iglob = ipole(iipole)
          call bspline_fill_site(iglob,i)
        end do
      end if
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
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
         call cmp_to_fmp_site(cmp,fmp(1,i))
      end do
c
      if (.not.use_polar) then
c       zero out the PME charge grid
        qgridin_2d = 0d0
c
c       MPI : Begin reception
c
        do i = 1, nrec_recep
          tag = nproc*rank + prec_recep(i) + 1
          call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $     n3mpimax,MPI_REAL8,prec_recep(i),tag,MPI_COMM_WORLD,req(tag),
     $     ierr)
        end do
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          iglob = ipole(iipole)
          call grid_mpole_site(iglob,i,fmp(1,i))
        end do
c
c     MPI : begin sending
c
        do i = 1, nrec_send
          proc = prec_send(i)
          tag = nproc*prec_send(i) + rank + 1
          call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $     2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,
     $     MPI_COMM_WORLD,req(tag),ierr)
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
c       do the reduction 'by hand'
c
        do i = 1, nrec_recep
          qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1) + 
     $      qgridmpi(:,:,:,:,i) 
        end do
c
c       Perform 3-D FFT forward transform
c

        call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $   n3mpimax)
c
c       make the scalar summation over reciprocal lattice
c
        if ((istart2(rank+1).eq.1).and.(jstart2(rank+1).eq.1).and.
     $     (kstart2(rank+1).eq.1)) then
             qfac_2d(1,1,1) = 0.0d0
        end if
        pterm = (pi/aewald)**2
        volterm = pi * volbox
        nf1 = (nfft1+1) / 2
        nf2 = (nfft2+1) / 2
        nf3 = (nfft3+1) / 2
ccc$omp parallel do default(shared) private(i,k3,j,k2,k1,m1,m2,m3)
ccc$omp+ private(r1,r2,r3,h1,h2,h3,hsq,term,expterm,denom)
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
     $          kstart2(rank+1)+1) = expterm
              end if
 10           continue
            end do
          end do
        end do
ccc$omp end parallel do
c
c       account for the zeroth grid point for a finite system
c
        if ((istart2(rank+1).eq.1).and.(jstart2(rank+1).eq.1).and.
     $     (kstart2(rank+1).eq.1)) then
          if (.not. use_bounds) then
             expterm = 0.5d0 * pi / xbox
             qfac_2d(1,1,1) = expterm
          end if
        end if
c
c       complete the transformation of the charge grid
c
c        do k = 1, ksize2(rank+1)
c           do j = 1, jsize2(rank+1)
c             do i = 1, isize2(rank+1)
c               term = qfac_2d(i,j,k)
c               qgridout_2d(1,i,j,k) = term*qgridout_2d(1,i,j,k)
c               qgridout_2d(2,i,j,k) = term*qgridout_2d(2,i,j,k)
c             end do
c           end do
c        end do
        qgridout_2d(1,:,:,:) = qfac_2d*qgridout_2d(1,:,:,:)
        qgridout_2d(2,:,:,:) = qfac_2d*qgridout_2d(2,:,:,:)
c
c       perform 3-D FFT backward transform
c
        call fft2d_backmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $   n3mpimax)
c
c       MPI : Begin reception
c
        do i = 1, nrec_send
          proc = prec_send(i)
          tag = nproc*rank + prec_send(i) + 1
          call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $     2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $     prec_send(i),tag,MPI_COMM_WORLD,reqbcast(tag),ierr)
        end do
c
c       MPI : begin sending
c
        do i = 1, nrec_recep
          tag = nproc*prec_recep(i) + rank + 1
          call MPI_ISEND(qgridin_2d,
     $     2*n1mpimax*n2mpimax*n3mpimax,
     $     MPI_REAL8,prec_recep(i),tag,MPI_COMM_WORLD,reqbcast(tag),
     $     ierr)
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
      end if
c
      e = 0.0d0
cc$omp parallel do default(shared) private(i,iipole,iglob,k)
cc$omp+ reduction(+:e)
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         if (.not.use_polar) then
           call fphi_mpole_site(iglob,i)
         end if
         do k = 1, 10
            e = e + fmp(k,i)*fphirec(k,i)
         end do
      end do
cc$omp end parallel do
c
      e = 0.5d0 * electric * e
      em = em + e
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      if (use_polar) then
         do i = 1, 3
            a(1,i) = dble(nfft1) * recip(i,1)
            a(2,i) = dble(nfft2) * recip(i,2)
            a(3,i) = dble(nfft3) * recip(i,3)
         end do
         e = 0.0d0
         do i = 1, npolerecloc
            iipole = polerecglob(i)
            iglob = ipole(iipole)
            fuind = 0d0
            do k = 1, 3
               fuind(k) = a(k,1)*uind(1,iipole) + a(k,2)*uind(2,iipole)
     &                         + a(k,3)*uind(3,iipole)
               e = e + fuind(k)*fphirec(k+1,i)
            end do
         end do
         e = 0.5d0 * electric * e
         ep = ep + e
      end if
c
c     deallocate local arrays
c
      deallocate (qgridmpi)
      deallocate (fmp)
      deallocate (req)
      deallocate (reqbcast)
c      call fftendmpi(planf,planb)
      return
      end
