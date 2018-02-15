c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c    "empole1pme" calculates the multipole and dipole polarization energy
c    and derivatives in PBC with PME with respect to cartesian coordinates
c
c
      subroutine empole1pme
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'timestat.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,ii,iglob
      integer iipole
      real*8 e,ei,eintra
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 term,f,fterm
      real*8 cii,dii,qii,uii
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 trq(3),trqi(3)
      real*8 frcx(3),frcy(3),frcz(3)
      real*8 time0,time1
c
c     zero out multipole and polarization energy and derivatives
c
      em = 0.0d0
      ep = 0.0d0
      eintra = 0.0d0
      dem = 0d0
      dep = 0d0
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induce dipole moment at each atom
c
      if (use_pmecore) then
        call newinduce_pme
      else
        call newinduce_pme2
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        time0 = mpi_wtime()
        call emrecip1
        time1 = mpi_wtime()
        timerec = timerec + time1 - time0
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        time0 = mpi_wtime()
        call ereal1d(eintra)
        time1 = mpi_wtime()
        timereal = timereal + time1 - time0
c
c    compute the Ewald self energy term over all the atoms
c
        term = 2.0d0 * aewald * aewald
        fterm = -f * aewald / sqrtpi
        do ii = 1, npoleloc
           iipole = poleglob(ii)
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
           uix = uind(1,iipole)
           uiy = uind(2,iipole)
           uiz = uind(3,iipole)
           cii = ci*ci
           dii = dix*dix + diy*diy + diz*diz
           qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &              + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
           uii = dix*uix + diy*uiy + diz*uiz
           e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
           ei = fterm * term * uii / 3.0d0
           em = em + e
           ep = ep + ei
        end do
c
c       compute the self-energy torque term due to induced dipole
c
        trq(1) = 0.0d0
        trq(2) = 0.0d0
        trq(3) = 0.0d0
        term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
        do ii = 1, npoleloc
           iipole = poleglob(ii)
           iglob = ipole(iipole)
           dix = rpole(2,iipole)
           diy = rpole(3,iipole)
           diz = rpole(4,iipole)
           uix = 0.5d0 * (uind(1,iipole)+uinp(1,iipole))
           uiy = 0.5d0 * (uind(2,iipole)+uinp(2,iipole))
           uiz = 0.5d0 * (uind(3,iipole)+uinp(3,iipole))
           trqi(1) = term * (diy*uiz-diz*uiy)
           trqi(2) = term * (diz*uix-dix*uiz)
           trqi(3) = term * (dix*uiy-diy*uix)
           call torque3 (iipole,trq,trqi,frcx,frcy,frcz,torquedir,
     $       torquedirp)
        end do
c
c       add the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
           xd = 0.0d0
           yd = 0.0d0
           zd = 0.0d0
           xu = 0.0d0
           yu = 0.0d0
           zu = 0.0d0
           do ii = 1, npoleloc
              iipole = poleglob(ii)
              i = loc(ipole(iipole))
              iglob = glob(i)
              dix = rpole(2,iipole)
              diy = rpole(3,iipole)
              diz = rpole(4,iipole)
              uix = uind(1,iipole)
              uiy = uind(2,iipole)
              uiz = uind(3,iipole)
              xd = xd + dix + rpole(1,iipole)*x(iglob)
              yd = yd + diy + rpole(1,iipole)*y(iglob)
              zd = zd + diz + rpole(1,iipole)*z(iglob)
              xu = xu + uix
              yu = yu + uiy
              zu = zu + uiz
           end do
           term = (2.0d0/3.0d0) * f * (pi/volbox)
           em = em + term*(xd*xd+yd*yd+zd*zd)
           ep = ep + term*(xd*xu+yd*yu+zd*zu)
        end if
c
c     intermolecular energy is total minus intramolecular part
c
        einter = einter + em + ep - eintra
        time1 = mpi_wtime()
      end if
      return
      end
c
c     subroutine emrecip1mpi : calculates reciprocal space part of the Ewald SUM with MPI
c     energy and derivatives
c
      subroutine emrecip1
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'bound.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'fft.i'
      include 'math.i'
      include 'mpole.i'
      include 'pme.i'
      include 'polar.i'
      include 'potent.i'
      include 'timestat.i'
      include 'virial.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,j,k,ii,tag,ierr,iglob,iipole,iloc
      integer status(MPI_STATUS_SIZE)
      integer j1,j2,j3
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 cphim(4),cphid(4)
      real*8 cphip(4)
      real*8 a(3,3),ftc(10,10)
      real*8, allocatable :: frc(:,:)
      real*8, allocatable :: trq(:,:)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: fphid(:,:)
      real*8, allocatable :: fphip(:,:)
      real*8, allocatable :: fphidp(:,:)
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
      real*8, allocatable :: qgrip(:,:,:,:)
      integer, allocatable :: req(:),reqbcast(:)
      integer nprocloc,commloc,rankloc,proc
      real*8 time0,time1
c
      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = MPI_COMM_WORLD
        rankloc  = rank
      end if
c
c     derivative indices into the fphi and fphidp arrays
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (frc(3,npolerecloc))
      allocate (trq(3,npolerecloc))
      allocate (fuind(3,npolerecloc))
      allocate (fuinp(3,npolerecloc))
      allocate (cmp(10,npolerecloc))
      allocate (fmp(10,npolerecloc))
      allocate (fphid(10,npolerecloc))
      allocate (fphip(10,npolerecloc))
      allocate (fphidp(20,npolerecloc))
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
c
      vxx = 0.0d0
      vyx = 0.0d0
      vzx = 0.0d0
      vyy = 0.0d0
      vzy = 0.0d0
      vzz = 0.0d0
c
c     increment the permanent multipole energy and gradient
c
      if (.not.use_polar) then
c
c       zero out the PME charge grid
c
        qgrid2in_2d = 0d0
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          iglob = ipole(iipole)
          call bspline_fill_site(iglob,i)
        end do
      end if
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         cmp(1,i) = rpole(1,iipole)
         cmp(2,i) = rpole(2,iipole)
         cmp(3,i) = rpole(3,iipole)
         cmp(4,i) = rpole(4,iipole)
         cmp(5,i) = rpole(5,iipole)
         cmp(6,i) = rpole(9,iipole)
         cmp(7,i) = rpole(13,iipole)
         cmp(8,i) = 2.0d0 * rpole(6,iipole)
         cmp(9,i) = 2.0d0 * rpole(7,iipole)
         cmp(10,i) = 2.0d0 * rpole(10,iipole)
         call cmp_to_fmp_site(cmp(1,i),fmp(1,i))
      end do
c
      if (.not.use_polar) then
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          iglob = ipole(iipole)
          call grid_mpole_site(iglob,i,fmp(1,i))
        end do
c
c     MPI : Begin reception
c
        do i = 1, nrec_recep
          tag = nprocloc*rankloc + prec_recep(i) + 1
          call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $     n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $     commloc,req(tag),ierr)
        end do
c
c     MPI : begin sending
c
        do i = 1, nrec_send
          proc = prec_send(i)
          tag = nprocloc*prec_send(i) + rankloc + 1
          call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $     2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,
     $     commloc,req(tag),ierr)
        end do
c
        do i = 1, nrec_recep
          tag = nprocloc*rankloc + prec_recep(i) + 1
          call MPI_WAIT(req(tag),status,ierr)
        end do
        do i = 1, nrec_send
          tag = nprocloc*prec_send(i) + rankloc + 1
          call MPI_WAIT(req(tag),status,ierr)
        end do
c
c       do the reduction 'by hand'

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
        if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $     (kstart2(rankloc+1).eq.1)) then
             qfac_2d(1,1,1) = 0.0d0
        end if
        pterm = (pi/aewald)**2
        volterm = pi * volbox
        nf1 = (nfft1+1) / 2
        nf2 = (nfft2+1) / 2
        nf3 = (nfft3+1) / 2
c$omp parallel do default(shared) private(i,k3,j,k2,k1,m1,m2,m3)
c$omp+ private(r1,r2,r3,h1,h2,h3,hsq,term,expterm,denom)
        do k3 = kstart2(rankloc+1),kend2(rankloc+1)
          do k2 = jstart2(rankloc+1),jend2(rankloc+1)
            do k1 = istart2(rankloc+1),iend2(rankloc+1)
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
              qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,
     $          k3-kstart2(rankloc+1)+1) = expterm
              end if
 10           continue
            end do
          end do
        end do
c$omp end parallel do
c
c       account for the zeroth grid point for a finite system
c
        if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $     (kstart2(rankloc+1).eq.1)) then
          if (.not. use_bounds) then
             expterm = 0.5d0 * pi / xbox
             qfac_2d(1,1,1) = expterm
          end if
        end if
c
c       complete the transformation of the charge grid
c
        do k = 1, ksize2(rankloc+1)
           do j = 1, jsize2(rankloc+1)
             do i = 1, isize2(rankloc+1)
               term = qfac_2d(i,j,k)
               qgridout_2d(1,i,j,k) = term*qgridout_2d(1,i,j,k)
               qgridout_2d(2,i,j,k) = term*qgridout_2d(2,i,j,k)
             end do
           end do
        end do
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
          tag = nprocloc*rankloc + prec_send(i) + 1
          call MPI_IRECV(qgridin_2d(1,1,1,1,i+1),
     $     2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $     prec_send(i),tag,commloc,reqbcast(tag),ierr)
        end do
c
c       MPI : begin sending
c
        do i = 1, nrec_recep
          tag = nprocloc*prec_recep(i) + rankloc + 1
          call MPI_ISEND(qgridin_2d,
     $     2*n1mpimax*n2mpimax*n3mpimax,
     $     MPI_REAL8,prec_recep(i),tag,commloc,reqbcast(tag),
     $     ierr)
        end do
c
        do i = 1, nrec_send
          tag = nprocloc*rankloc + prec_send(i) + 1
          call MPI_WAIT(reqbcast(tag),status,ierr)
        end do
        do i = 1, nrec_recep
          tag = nprocloc*prec_recep(i) + rankloc + 1
          call MPI_WAIT(reqbcast(tag),status,ierr)
        end do
c
        do i = 1, npolerecloc
          iipole = polerecglob(i)
          iglob = ipole(iipole)
          call fphi_mpole_site(iglob,i)
        end do
      end if
c
      time0 = mpi_wtime()
      e = 0.0d0
      do i = 1, npolerecloc
         do j = 1, 20
            fphirec(j,i) = electric * fphirec(j,i)
         end do
        call fphi_to_cphi_site (fphirec(1,i),cphirec(1,i))
      end do
c$omp parallel do default(shared) private(i,f1,f2,f3,k)
c$omp+ reduction(+:e)
      do i = 1, npolerecloc
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 10
            e = e + fmp(k,i)*fphirec(k,i)
            f1 = f1 + fmp(k,i)*fphirec(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphirec(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphirec(deriv3(k),i)
         end do
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
      end do
c$omp end parallel do
      time1 = mpi_wtime()
      timegrid2 = timegrid2 + time1 - time0
      e = 0.5d0 * e
      em = em + e
      do ii = 1, npolerecloc
        iipole = polerecglob(ii)
        iglob = ipole(iipole)
        i = locrec1(iglob)
        demrec(1,i) = frc(1,ii)
        demrec(2,i) = frc(2,ii)
        demrec(3,i) = frc(3,ii)
      end do
c
c     distribute torques into the permanent multipole gradient
c
c$omp parallel do default(shared) private(i)
      do i = 1, npolerecloc
         trq(1,i) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*
     $                 cphirec(4,i)
     &                 + 2.0d0*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
     &                 + cmp(9,i)*cphirec(8,i) + cmp(10,i)*
     $                 cphirec(6,i)
     &                 - cmp(8,i)*cphirec(9,i) - cmp(10,i)*
     $                 cphirec(7,i)
         trq(2,i) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*
     $                 cphirec(2,i)
     &                 + 2.0d0*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
     &                 + cmp(8,i)*cphirec(10,i) + cmp(9,i)*
     $                 cphirec(7,i)
     &                 - cmp(9,i)*cphirec(5,i) - cmp(10,i)*
     $                 cphirec(8,i)
         trq(3,i) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*
     $                 cphirec(3,i)
     &                 + 2.0d0*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
     &                 + cmp(8,i)*cphirec(5,i) + cmp(10,i)*
     $                 cphirec(9,i)
     &                 - cmp(8,i)*cphirec(6,i) - cmp(9,i)*
     $                 cphirec(10,i)
      end do
c$omp end parallel do
c
      call torque2(trq,torquerec)
c
c     permanent multipole contribution to the internal virial
c
c$omp parallel do default(shared) private(i)
c$omp+ reduction(+:vxx,vyx,vzx,vyy,vzy,vzz)
      do i = 1, npolerecloc
         vxx = vxx - cmp(2,i)*cphirec(2,i) - 2.0d0*cmp(5,i)*
     $             cphirec(5,i)
     &             - cmp(8,i)*cphirec(8,i) - cmp(9,i)*cphirec(9,i)
         vyx = vyx - 0.5d0*(cmp(3,i)*cphirec(2,i)+cmp(2,i)*
     $              cphirec(3,i))
     &             - (cmp(5,i)+cmp(6,i))*cphirec(8,i)
     &             - 0.5d0*cmp(8,i)*(cphirec(5,i)+cphirec(6,i))
     &             - 0.5d0*(cmp(9,i)*cphirec(10,i)+cmp(10,i)*
     $             cphirec(9,i))
         vzx = vzx - 0.5d0*(cmp(4,i)*cphirec(2,i)+cmp(2,i)*
     $             cphirec(4,i))
     &             - (cmp(5,i)+cmp(7,i))*cphirec(9,i)
     &             - 0.5d0*cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
     &             - 0.5d0*(cmp(8,i)*cphirec(10,i)+cmp(10,i)*
     $             cphirec(8,i))
         vyy = vyy - cmp(3,i)*cphirec(3,i) - 2.0d0*cmp(6,i)*
     $             cphirec(6,i)
     &             - cmp(8,i)*cphirec(8,i) - cmp(10,i)*
     $             cphirec(10,i)
         vzy = vzy - 0.5d0*(cmp(4,i)*cphirec(3,i)+cmp(3,i)*
     $             cphirec(4,i))
     &             - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
     &             - 0.5d0*cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
     &             - 0.5d0*(cmp(8,i)*cphirec(9,i)+cmp(9,i)*
     $             cphirec(8,i))
         vzz = vzz - cmp(4,i)*cphirec(4,i) - 2.0d0*cmp(7,i)*
     $              cphirec(7,i)
     &             - cmp(9,i)*cphirec(9,i) - cmp(10,i)*
     $             cphirec(10,i)
      end do
c$omp end parallel do
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      if (use_polar) then
c
c     zero out the PME grid
c
         qgrid2in_2d = 0d0
c
         do i = 1, 3
            a(1,i) = dble(nfft1) * recip(i,1)
            a(2,i) = dble(nfft2) * recip(i,2)
            a(3,i) = dble(nfft3) * recip(i,3)
         end do
         time0 = mpi_wtime()
         do ii = 1, npolerecloc
            iipole = polerecglob(ii)
            iglob = ipole(iipole)
            do j = 1, 3
               fuind(j,ii) = a(j,1)*uind(1,iipole)+a(j,2)*uind(2,iipole)
     &                          + a(j,3)*uind(3,iipole)
               fuinp(j,ii) = a(j,1)*uinp(1,iipole)+a(j,2)*uinp(2,iipole)
     &                          + a(j,3)*uinp(3,iipole)
            end do
            call grid_uind_site(iglob,ii,fuind(1,ii),fuinp(1,ii),
     $       qgrid2in_2d)
         end do
         time1 = mpi_wtime()
         timegrid1 = timegrid1 + time1 - time0
c
         do i = 1, nrec_recep
           tag = nprocloc*rankloc + prec_recep(i) + 1
           call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $      n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $      commloc,req(tag),ierr)
         end do
c
c     MPI : begin sending
c
         time0 = mpi_wtime()
         do i = 1, nrec_send
           tag = nprocloc*prec_send(i) + rankloc + 1
           call MPI_ISEND(qgrid2in_2d(1,1,1,1,i+1),
     $      2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $      prec_send(i),tag,commloc,req(tag),ierr)
         end do
c
         do i = 1, nrec_recep
           tag = nprocloc*rankloc + prec_recep(i) + 1
           call MPI_WAIT(req(tag),status,ierr)
         end do
         do i = 1, nrec_send
           tag = nprocloc*prec_send(i) + rankloc + 1
           call MPI_WAIT(req(tag),status,ierr)
         end do
c
c     do the reduction 'by hand'
c
         do i = 1, nrec_recep
           qgrid2in_2d(:,:,:,:,1) = qgrid2in_2d(:,:,:,:,1)+
     $      qgridmpi(:,:,:,:,i)
         end do
         time1 = mpi_wtime()
         timerecreccomm = timerecreccomm + time1 - time0
c
         time0 = mpi_wtime()
         call fft2d_frontmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     $    n3mpimax)
         time1 = mpi_wtime()
         timeffts = timeffts + time1 - time0
c
c     account for the zeroth grid point for a finite system
c
         if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     $      .and.(kstart2(rankloc+1).eq.1)) then
           if (.not. use_bounds) then
              expterm = 0.5d0 * pi / xbox
              struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
     $          qgrid2in_2d(2,1,1,1,1)**2
              e = 0.5d0 * expterm * struc2
           end if
         end if
c
c     complete the transformation of the PME grid
c
         time0 = mpi_wtime()
         do k = 1, ksize2(rankloc+1)
            do j = 1, jsize2(rankloc+1)
              do i = 1, isize2(rankloc+1)
                term = qfac_2d(i,j,k)
                qgrid2out_2d(1,i,j,k) = term*qgrid2out_2d(1,i,j,k)
                qgrid2out_2d(2,i,j,k) = term*qgrid2out_2d(2,i,j,k)
              end do
            end do
         end do
         time1 = mpi_wtime()
         timescalar = timescalar + time1 - time0
c
c     perform 3-D FFT backward transform and get potential
c
         time0 = mpi_wtime()
         call fft2d_backmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     $    n3mpimax)
         time1 = mpi_wtime()
         timeffts = timeffts + time1 - time0
c
c     MPI : Begin reception
c
          do i = 1, nrec_send
            tag = nprocloc*rankloc + prec_send(i) + 1
            call MPI_IRECV(qgrid2in_2d(1,1,1,1,i+1),
     $       2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $       prec_recep(i),tag,commloc,reqbcast(tag),ierr)
          end do
c
c     MPI : begin sending
c
          do i = 1, nrec_recep
            tag = nprocloc*prec_recep(i) + rankloc + 1
            call MPI_ISEND(qgrid2in_2d(1,1,1,1,1),
     $       2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $       prec_send(i),tag,commloc,reqbcast(tag),ierr)
          end do
c
          time0 = mpi_wtime()
          do i = 1, nrec_send
            tag = nprocloc*rankloc + prec_send(i) + 1
            call MPI_WAIT(reqbcast(tag),status,ierr)
          end do
          do i = 1, nrec_recep
            tag = nprocloc*prec_recep(i) + rankloc + 1
            call MPI_WAIT(reqbcast(tag),status,ierr)
          end do
          time1 = mpi_wtime()
          timerecreccomm = timerecreccomm + time1 - time0
c
         time0 = mpi_wtime()
c$omp parallel do default(shared) private(ii,iipole,iglob,j)
         do ii = 1, npolerecloc
           iipole = polerecglob(ii)
           iglob = ipole(iipole)
           call fphi_uind_site(iglob,ii,fphid(1,ii),fphip(1,ii),
     $          fphidp(1,ii))
           do j = 1, 10
              fphid(j,ii) = electric * fphid(j,ii)
              fphip(j,ii) = electric * fphip(j,ii)
           end do
           do j = 1, 20
              fphidp(j,ii) = electric * fphidp(j,ii)
           end do
         end do
c$omp end parallel do
c
c     increment the induced dipole energy and gradient
c
         e = 0.0d0
c$omp parallel do default(shared) private(i,k,j1,j2,j3,f1,f2,f3)
c$omp+ reduction(+:e)
         do i = 1, npolerecloc
            f1 = 0.0d0
            f2 = 0.0d0
            f3 = 0.0d0
            do k = 1, 3
               j1 = deriv1(k+1)
               j2 = deriv2(k+1)
               j3 = deriv3(k+1)
               e = e + fuind(k,i)*fphirec(k+1,i)
               f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphirec(j1,i)
     &                 + fuind(k,i)*fphip(j1,i)
     &                 + fuinp(k,i)*fphid(j1,i)
               f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphirec(j2,i)
     &                 + fuind(k,i)*fphip(j2,i)
     &                 + fuinp(k,i)*fphid(j2,i)
               f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphirec(j3,i)
     &                 + fuind(k,i)*fphip(j3,i)
     &                 + fuinp(k,i)*fphid(j3,i)
            end do
            do k = 1, 10
               f1 = f1 + fmp(k,i)*fphidp(deriv1(k),i)
               f2 = f2 + fmp(k,i)*fphidp(deriv2(k),i)
               f3 = f3 + fmp(k,i)*fphidp(deriv3(k),i)
            end do
            f1 = 0.5d0 * dble(nfft1) * f1
            f2 = 0.5d0 * dble(nfft2) * f2
            f3 = 0.5d0 * dble(nfft3) * f3
            frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
            frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
            frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         end do
c$omp end parallel do
         e = 0.5d0 * e
         ep = ep + e
         do ii = 1, npolerecloc
           iipole = polerecglob(ii)
           iglob = ipole(iipole)
           iloc = loc(iglob)
           i = locrec1(iglob)
           deprec(1,i) = frc(1,ii)
           deprec(2,i) = frc(2,ii)
           deprec(3,i) = frc(3,ii)
         end do
c
c     set the potential to be the induced dipole average
c
c$omp parallel do default(shared) private(i,k)
         do i = 1, npolerecloc
            do k = 1, 10
               fphidp(k,i) = 0.5d0 * fphidp(k,i)
            end do
            call fphi_to_cphi_site(fphidp(1,i),cphirec(1,i))
         end do
c$omp end parallel do
         time1 = mpi_wtime()
         timegrid2 = timegrid2 + time1 - time0
c
c     distribute torques into the induced dipole gradient
c
c$omp parallel do default(shared) private(i)
         do i = 1, npolerecloc
            trq(1,i) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*
     $          cphirec(4,i)
     &          + 2.0d0*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
     &          + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
     &          - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
            trq(2,i) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*
     $          cphirec(2,i)
     &          + 2.0d0*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
     &          + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
     &          - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
            trq(3,i) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*
     $          cphirec(3,i)
     &          + 2.0d0*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
     &          + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
     &          - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)
         end do
c$omp end parallel do
        do i = 1, npolerecloc
           frc(1,i) = 0.0d0
           frc(2,i) = 0.0d0
           frc(3,i) = 0.0d0
        end do
        call torque2(trq,torquerecp)
        call frac_to_cart(ftc)
c
c     induced dipole contribution to the internal virial
c
c$omp parallel do default(shared) private(i,iipole,j,k)
c$omp+ reduction(+:vxx,vyx,vzx,vzy,vyy,vzz)
         do i = 1, npolerecloc
            iipole = polerecglob(i)
            do j = 2, 4
               cphim(j) = 0.0d0
               cphid(j) = 0.0d0
               cphip(j) = 0.0d0
               do k = 2, 4
                  cphim(j) = cphim(j) + ftc(j,k)*fphirec(k,i)
                  cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
                  cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
               end do
            end do
            vxx = vxx - cphirec(2,i)*cmp(2,i)
     &            - 0.5d0*(cphim(2)*(uind(1,iipole)+uinp(1,iipole))
     &            +cphid(2)*uinp(1,iipole)+cphip(2)*uind(1,iipole))
            vyx = vyx - 0.5d0*(cphirec(2,i)*cmp(3,i)+cphirec(3,i)*
     $            cmp(2,i))
     &            - 0.25d0*(cphim(2)*(uind(2,iipole)+uinp(2,iipole))
     &            +cphim(3)*(uind(1,iipole)+uinp(1,iipole))
     &            +cphid(2)*uinp(2,iipole)+cphip(2)*uind(2,iipole)
     &            +cphid(3)*uinp(1,iipole)+cphip(3)*uind(1,iipole))
            vzx = vzx - 0.5d0*(cphirec(2,i)*cmp(4,i)+cphirec(4,i)*
     $            cmp(2,i))
     &            - 0.25d0*(cphim(2)*(uind(3,iipole)+uinp(3,iipole))
     &            +cphim(4)*(uind(1,iipole)+uinp(1,iipole))
     &            +cphid(2)*uinp(3,iipole)+cphip(2)*uind(3,iipole)
     &            +cphid(4)*uinp(1,iipole)+cphip(4)*uind(1,iipole))
            vyy = vyy - cphirec(3,i)*cmp(3,i)
     &            - 0.5d0*(cphim(3)*(uind(2,iipole)+uinp(2,iipole))
     &            +cphid(3)*uinp(2,iipole)+cphip(3)*uind(2,iipole))
            vzy = vzy - 0.5d0*(cphirec(3,i)*cmp(4,i)+cphirec(4,i)*
     $            cmp(3,i))
     &            - 0.25d0*(cphim(3)*(uind(3,iipole)+uinp(3,iipole))
     &            +cphim(4)*(uind(2,iipole)+uinp(2,iipole))
     &            +cphid(3)*uinp(3,iipole)+cphip(3)*uind(3,iipole)
     &            +cphid(4)*uinp(2,iipole)+cphip(4)*uind(2,iipole))
            vzz = vzz - cphirec(4,i)*cmp(4,i)
     &            - 0.5d0*(cphim(4)*(uind(3,iipole)+uinp(3,iipole))
     &            +cphid(4)*uinp(3,iipole)+cphip(4)*uind(3,iipole))
            vxx = vxx - 2.0d0*cmp(5,i)*cphirec(5,i) - cmp(8,i)*
     $            cphirec(8,i)
     &            - cmp(9,i)*cphirec(9,i)
            vyx = vyx - (cmp(5,i)+cmp(6,i))*cphirec(8,i)
     &            - 0.5d0*(cmp(8,i)*(cphirec(6,i)+cphirec(5,i))
     &            +cmp(9,i)*cphirec(10,i)+cmp(10,i)*cphirec(9,i))
            vzx = vzx - (cmp(5,i)+cmp(7,i))*cphirec(9,i)
     &            - 0.5d0*(cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
     &             +cmp(8,i)*cphirec(10,i)+cmp(10,i)*cphirec(8,i))
            vyy = vyy - 2.0d0*cmp(6,i)*cphirec(6,i) - cmp(8,i)*
     $            cphirec(8,i)
     &            - cmp(10,i)*cphirec(10,i)
            vzy = vzy - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
     &            - 0.5d0*(cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
     &            +cmp(8,i)*cphirec(9,i)+cmp(9,i)*cphirec(8,i))
            vzz = vzz - 2.0d0*cmp(7,i)*cphirec(7,i) -
     $                cmp(9,i)*cphirec(9,i)
     &                - cmp(10,i)*cphirec(10,i)
         end do
c$omp end parallel do
      end if
c
c     increment the internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vyx
      vir(3,1) = vir(3,1) + vzx
      vir(1,2) = vir(1,2) + vyx
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vzy
      vir(1,3) = vir(1,3) + vzx
      vir(2,3) = vir(2,3) + vzy
      vir(3,3) = vir(3,3) + vzz
c
c
      vxx = 0.0d0
      vyx = 0.0d0
      vzx = 0.0d0
      vyy = 0.0d0
      vzy = 0.0d0
      vzz = 0.0d0
c
c
c     copy multipole moments and coordinates to local storage
c
c     perform dynamic allocation of some local arrays
c
      allocate (qgrip(2,isize2(rankloc+1),jsize2(rankloc+1),
     $ ksize2(rankloc+1)))
c
c     assign permanent and induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c
      if (use_polar) then
c
c    zero out the grid
c
         qgridin_2d = 0d0
c
         time0 = mpi_wtime()
         do ii = 1, npolerecloc
            iipole = polerecglob(ii)
            iglob = ipole(iipole)
            do j = 2, 4
               cmp(j,ii) = cmp(j,ii) + uinp(j-1,iipole)
            end do
            call cmp_to_fmp_site (cmp(1,ii),fmp(1,ii))
            call grid_mpole_site(iglob,ii,fmp(1,ii))
         end do
         time1 = mpi_wtime()
         timegrid1 = timegrid1 + time1 - time0
c
c     MPI : Begin reception
c
         do i = 1, nrec_recep
           tag = nprocloc*rankloc + prec_recep(i) + 1
           call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $      n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $      commloc,req(tag),ierr)
         end do
c
c     MPI : begin sending
c
         time0 = mpi_wtime()
         do i = 1, nrec_send
           tag = nprocloc*prec_send(i) + rankloc + 1
           call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $      2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $      prec_send(i),tag,commloc,req(tag),ierr)
         end do
c
         do i = 1, nrec_recep
           tag = nprocloc*rankloc + prec_recep(i) + 1
           call MPI_WAIT(req(tag),status,ierr)
         end do
         do i = 1, nrec_send
           tag = nprocloc*prec_send(i) + rankloc + 1
           call MPI_WAIT(req(tag),status,ierr)
         end do
c
c     do the reduction 'by hand'
c
         do i = 1, nrec_recep
           qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1)+
     $      qgridmpi(:,:,:,:,i)
         end do
         time1 = mpi_wtime()
         timerecreccomm = timerecreccomm + time1 - time0
c
         time0 = mpi_wtime()
         call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $    n3mpimax)
         time1 = mpi_wtime()
         timeffts = timeffts + time1 - time0
         do k = 1, ksize2(rankloc+1)
            do j = 1, jsize2(rankloc+1)
               do i = 1, isize2(rankloc+1)
                  qgrip(1,i,j,k) = qgridout_2d(1,i,j,k)
                  qgrip(2,i,j,k) = qgridout_2d(2,i,j,k)
               end do
            end do
         end do
c
c     zero out the PME grid
c
         qgridin_2d = 0d0
         time0 = mpi_wtime()
         do ii = 1, npolerecloc
            iipole = polerecglob(ii)
            iglob = ipole(iipole)
            do j = 2, 4
               cmp(j,ii) = cmp(j,ii) + uind(j-1,iipole)-uinp(j-1,iipole)
            end do
            call cmp_to_fmp_site(cmp(1,ii),fmp(1,ii))
            call grid_mpole_site(iglob,ii,fmp(1,ii))
         end do
         time1 = mpi_wtime()
         timegrid1 = timegrid1 + time1 - time0
c
c     MPI : Begin reception
c
         do i = 1, nrec_recep
           tag = nprocloc*rankloc + prec_recep(i) + 1
           call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $      n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $      commloc,req(tag),ierr)
         end do
c
c     MPI : begin sending
c
         time0 = mpi_wtime()
         do i = 1, nrec_send
           tag = nprocloc*prec_send(i) + rankloc + 1
           call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $      2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $      prec_send(i),tag,commloc,req(tag),ierr)
         end do
c
         do i = 1, nrec_recep
           tag = nprocloc*rankloc + prec_recep(i) + 1
           call MPI_WAIT(req(tag),status,ierr)
         end do
         do i = 1, nrec_send
           tag = nprocloc*prec_send(i) + rankloc + 1
           call MPI_WAIT(req(tag),status,ierr)
         end do
c
c     do the reduction 'by hand'
c
         do i = 1, nrec_recep
           qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1)+
     $      qgridmpi(:,:,:,:,i)
         end do
         time1 = mpi_wtime()
         timerecreccomm = timerecreccomm + time1 - time0
c
         time0 = mpi_wtime()
         call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $    n3mpimax)
         time1 = mpi_wtime()
         timeffts = timeffts + time1 - time0
         do ii = 1, npolerecloc
            iipole = polerecglob(ii)
            do j = 2, 4
               cmp(j,ii) = cmp(j,ii) - uind(j-1,iipole)
            end do
         end do
      else
c
         qgridin_2d = 0d0
         time0 = mpi_wtime()
         do ii = 1, npolerecloc
            iipole = polerecglob(ii)
            iglob = ipole(iipole)
            call cmp_to_fmp_site(cmp(1,ii),fmp(1,ii))
            call grid_mpole_site(iglob,ii,fmp(1,ii))
         end do
         time1 = mpi_wtime()
         timegrid1 = timegrid1 + time1 - time0
c
c     MPI : Begin reception
c
         do i = 1, nrec_recep
           tag = nprocloc*rankloc + prec_recep(i) + 1
           call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $      n3mpimax,MPI_REAL8,prec_recep(i),tag,
     $      commloc,req(tag),ierr)
         end do
c
c     MPI : begin sending
c
         time0 = mpi_wtime()
         do i = 1, nrec_send
           tag = nprocloc*prec_send(i) + rankloc + 1
           call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $      2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,
     $      prec_send(i),tag,commloc,req(tag),ierr)
         end do
c
         do i = 1, nrec_recep
           tag = nprocloc*rankloc + prec_recep(i) + 1
           call MPI_WAIT(req(tag),status,ierr)
         end do
         do i = 1, nrec_send
           tag = nprocloc*prec_send(i) + rankloc + 1
           call MPI_WAIT(req(tag),status,ierr)
         end do
c
c     do the reduction 'by hand'

         do i = 1, nrec_recep
           qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1)+
     $      qgridmpi(:,:,:,:,i)
         end do
         time1 = mpi_wtime()
         timerecreccomm = timerecreccomm + time1 - time0
c
         time0 = mpi_wtime()
         call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $    n3mpimax)
         time1 = mpi_wtime()
         timeffts = timeffts + time1 - time0
         qgrip = qgridout_2d
         qgridin_2d = 0d0
      end if
c
c     make the scalar summation over reciprocal lattice

      time0 = mpi_wtime()
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0d0
      end if
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
cc$omp parallel do default(shared) private(i,k3,j,k2,k1,m1,m2,m3)
cc$omp+ private(r1,r2,r3,h1,h2,h3,hsq,term,expterm,denom,struc2,vterm)
cc$omp+ private(eterm)
cc$omp+ reduction(+:vxx) reduction(+:vyx) reduction(+:vzx)
cc$omp+ reduction(+:vyy) reduction(+:vzy) reduction(+:vzz)
      do k3 = kstart2(rankloc+1),kend2(rankloc+1)
        do k2 = jstart2(rankloc+1),jend2(rankloc+1)
          do k1 = istart2(rankloc+1),iend2(rankloc+1)
            m1 = k1 - 1
            m2 = k2 - 1
            m3 = k3 - 1
            if (k1 .gt. nf1)  m1 = m1 - nfft1
            if (k2 .gt. nf2)  m2 = m2 - nfft2
            if (k3 .gt. nf3)  m3 = m3 - nfft3
            if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 20
            r1 = dble(m1)
            r2 = dble(m2)
            r3 = dble(m3)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = 0.0d0
            if (term .gt. -50.0d0) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               end if
               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $          k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)*
     $          qgrip(1,k1-istart2(rankloc+1)+1,
     $          k2-jstart2(rankloc+1)+1,
     $          k3-kstart2(rankloc+1)+1)
     &          + qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-
     $          jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)*
     $          qgrip(2,k1-istart2(rankloc+1)+1,
     $          k2-jstart2(rankloc+1)+1,
     $          k3-kstart2(rankloc+1)+1)
               eterm = 0.5d0 * electric * expterm * struc2
               vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
               vxx = vxx + h1*h1*vterm - eterm
               vyx = vyx + h2*h1*vterm
               vzx = vzx + h3*h1*vterm
               vyy = vyy + h2*h2*vterm - eterm
               vzy = vzy + h3*h2*vterm
               vzz = vzz + h3*h3*vterm - eterm
            end if
            qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
     $        kstart2(rankloc+1)+1) = expterm
 20         continue
          end do
        end do
      end do
cc$omp end parallel do
      time1 = mpi_wtime()
      timescalar = timescalar + time1 - time0
c
c     increment the internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vyx
      vir(3,1) = vir(3,1) + vzx
      vir(1,2) = vir(1,2) + vyx
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vzy
      vir(1,3) = vir(1,3) + vzx
      vir(2,3) = vir(2,3) + vzy
      vir(3,3) = vir(3,3) + vzz
c
c     perform deallocation of some local arrays
c
      deallocate (qgrip)
c
      deallocate (frc)
      deallocate (trq)
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (cmp)
      deallocate (fmp)
      deallocate (fphid)
      deallocate (fphip)
      deallocate (fphidp)
      deallocate (qgridmpi)
      deallocate (req)
      deallocate (reqbcast)
      return
      end
c
c
c
c     "ereal1dmpi" evaluates the real space portion of the regular Ewald
c     summation energy and gradient due to atomic multipole interactions
c     and dipole polarizability
c
c
      subroutine ereal1d (eintra)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'virial.i'
      include 'openmp.i'
      include 'mpif.h'
      integer i,j,k,iglob,kglob,kbis
      integer ii,inl,kkk,iipole,kkpole
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5
      real*8 scale7
      real*8 temp3,temp5,temp7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 usc3,usc5
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 erl,erli
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ftm2r(3),ftm2ri(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 ttm2ri(3),ttm3ri(3)
      real*8 dixdk(3)
      real*8 dkxui(3),dixuk(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 sci(8),scip(8)
      real*8 gli(7),glip(7)
      real*8 gf(7),gfi(6)
      real*8 gfr(7),gfri(6)
      real*8 gti(6),gtri(6)
      real*8 virt(3,3),emtt,eptt
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: demtt(:,:),deptt(:,:)
      real*8, allocatable :: torquedirt(:,:),torquedirpt(:,:)
      logical dorl,dorli
      character*6 mode
      external erfc
c
 1000 format(' Warning, system moved too much since last neighbor list
     $   update, try lowering nlupdate')
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (demtt(3,nbloc))
      allocate (deptt(3,nbloc))
      allocate (torquedirt(3,nbloc))
      allocate (torquedirpt(3,nbloc))
c
c
c     set arrays needed to scale connected atom interactions
c
      mscale = 1.0d0
      pscale = 1.0d0
      dscale = 1.0d0
      uscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     initialize local variables for OpenMP calculation
c
      emtt = em
      eptt = ep
      demtt = dem
      deptt = dep
      torquedirt = torquedir
      torquedirpt = torquedirp
      virt = vir
c
c     set OpenMP directives for the major loop structure
c
c     compute the real space portion of the Ewald summation
c
c$omp parallel do default(shared) firstprivate(f)
c$omp+ private(ii,iipole,i,iglob,inl,pdi,pti,ci,di,qi,j,k,kkk,kglob)
c$omp+ private(kkpole,kbis,xr,yr,zr)
c$omp+ private(r2,r,ck,dk,qk,ralpha,alsq2,alsq2n,exp2a,bfac,bn)
c$omp+ private(rr1,rr3,rr5,rr7,rr9,rr11,scale3,scale5,scale7)
c$omp+ private(ddsc3,ddsc5,ddsc7,damp,pgamma,expdamp,temp3,temp5)
c$omp+ private(dsc3,dsc5,dsc7,psc3,psc5,psc7,usc3,usc5)
c$omp+ private(temp7,dixdk,dixuk,dkxui,dixukp,dkxuip,dixr,dkxr,qir,qkr)
c$omp+ private(qiqkr,qkqir,qixqk,rxqir,rxqkr,rxqikr,rxqkir,qkrxqir)
c$omp+ private(qidk,qkdi,qiuk,qkui,qiukp,qkuip,dixqkr,dkxqir,uixqkr)
c$omp+ private(ukxqir,uixqkrp,ukxqirp,rxqidk,rxqkdi,rxqiuk,rxqkui)
c$omp+ private(rxqiukp,rxqkuip,sc,sci,scip,gl,gli,glip,e,ei,erl,erli)
c$omp+ private(dorl,dorli,ftm2r,ftm2ri,ttm2r,ttm2ri,ttm3r,ttm3ri,gf)
c$omp+ private(ftm2,gfr,gfi,ftm2i,gfri,fridmp,findmp)
c$omp+ private(ttm2,ttm3)
c$omp+ private(gti,ttm2i,gtri,ttm3i,frcxi,frcyi,frczi,iaz,iax,iay,xiz)
c$omp+ private(frcxk,frcyk,frczk,kaz,kax,kay,xkz,ykz,zkz,xkx,ykx,zkx)
c$omp+ private(xky,yky,zky)
c$omp+ private(yiz,ziz,xix,yix,zix,xiy,yiy,ziy,vxx,vyx,vzx,vyy,vzy,vzz)
c$omp+ firstprivate(mscale,pscale,dscale,uscale)
c$omp+ reduction(+:emtt,eptt,virt,demtt,deptt)
c$omp+ reduction(+:torquedirt,torquedirpt,eintra)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         i = loc(ipole(iipole))
         iglob = glob(i)
         inl = polelocnl(iipole)
         if (inl.eq.0) then
           write(iout,1000)
           cycle
         end if
         pdi = pdamp(iipole)
         pti = thole(iipole)
         ci = rpole(1,iipole)
         di(1) = rpole(2,iipole)
         di(2) = rpole(3,iipole)
         di(3) = rpole(4,iipole)
         qi(1) = rpole(5,iipole)
         qi(2) = rpole(6,iipole)
         qi(3) = rpole(7,iipole)
         qi(4) = rpole(8,iipole)
         qi(5) = rpole(9,iipole)
         qi(6) = rpole(10,iipole)
         qi(7) = rpole(11,iipole)
         qi(8) = rpole(12,iipole)
         qi(9) = rpole(13,iipole)
c
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = m2scale
            pscale(i12(j,iglob)) = p2scale
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = m3scale
            pscale(i13(j,iglob)) = p3scale
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = m4scale
            pscale(i14(j,iglob)) = p4scale
            do k = 1, np11(iglob)
                if (i14(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i14(j,iglob)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = m5scale
            pscale(i15(j,iglob)) = p5scale
         end do
         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = d1scale
            uscale(ip11(j,iglob)) = u1scale
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = d2scale
            uscale(ip12(j,iglob)) = u2scale
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = d3scale
            uscale(ip13(j,iglob)) = u3scale
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = d4scale
            uscale(ip14(j,iglob)) = u4scale
         end do
         do kkk = 1, nelst(inl)
          kglob = elst(kkk,inl)
          kkpole = pollist(kglob)
          kbis = loc(kglob)
            if (kbis.eq.0) then
              write(iout,1000)
              cycle
            end if
            xr = x(kglob) - x(iglob)
            yr = y(kglob) - y(iglob)
            zr = z(kglob) - z(iglob)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,kkpole)
               dk(1) = rpole(2,kkpole)
               dk(2) = rpole(3,kkpole)
               dk(3) = rpole(4,kkpole)
               qk(1) = rpole(5,kkpole)
               qk(2) = rpole(6,kkpole)
               qk(3) = rpole(7,kkpole)
               qk(4) = rpole(8,kkpole)
               qk(5) = rpole(9,kkpole)
               qk(6) = rpole(10,kkpole)
               qk(7) = rpole(11,kkpole)
               qk(8) = rpole(12,kkpole)
               qk(9) = rpole(13,kkpole)
c
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(kkpole)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(kkpole))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(kglob)
               dsc5 = 1.0d0 - scale5*dscale(kglob)
               dsc7 = 1.0d0 - scale7*dscale(kglob)
               psc3 = 1.0d0 - scale3*pscale(kglob)
               psc5 = 1.0d0 - scale5*pscale(kglob)
               psc7 = 1.0d0 - scale7*pscale(kglob)
               usc3 = 1.0d0 - scale3*uscale(kglob)
               usc5 = 1.0d0 - scale5*uscale(kglob)
c
c     construct necessary auxiliary vectors
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,kkpole) - di(3)*uind(2,kkpole)
               dixuk(2) = di(3)*uind(1,kkpole) - di(1)*uind(3,kkpole)
               dixuk(3) = di(1)*uind(2,kkpole) - di(2)*uind(1,kkpole)
               dkxui(1) = dk(2)*uind(3,iipole) - dk(3)*uind(2,iipole)
               dkxui(2) = dk(3)*uind(1,iipole) - dk(1)*uind(3,iipole)
               dkxui(3) = dk(1)*uind(2,iipole) - dk(2)*uind(1,iipole)
               dixukp(1) = di(2)*uinp(3,kkpole) - di(3)*uinp(2,kkpole)
               dixukp(2) = di(3)*uinp(1,kkpole) - di(1)*uinp(3,kkpole)
               dixukp(3) = di(1)*uinp(2,kkpole) - di(2)*uinp(1,kkpole)
               dkxuip(1) = dk(2)*uinp(3,iipole) - dk(3)*uinp(2,iipole)
               dkxuip(2) = dk(3)*uinp(1,iipole) - dk(1)*uinp(3,iipole)
               dkxuip(3) = dk(1)*uinp(2,iipole) - dk(2)*uinp(1,iipole)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,kkpole) + qi(4)*uind(2,kkpole)
     &                      + qi(7)*uind(3,kkpole)
               qiuk(2) = qi(2)*uind(1,kkpole) + qi(5)*uind(2,kkpole)
     &                      + qi(8)*uind(3,kkpole)
               qiuk(3) = qi(3)*uind(1,kkpole) + qi(6)*uind(2,kkpole)
     &                      + qi(9)*uind(3,kkpole)
               qkui(1) = qk(1)*uind(1,iipole) + qk(4)*uind(2,iipole)
     &                      + qk(7)*uind(3,iipole)
               qkui(2) = qk(2)*uind(1,iipole) + qk(5)*uind(2,iipole)
     &                      + qk(8)*uind(3,iipole)
               qkui(3) = qk(3)*uind(1,iipole) + qk(6)*uind(2,iipole)
     &                      + qk(9)*uind(3,iipole)
               qiukp(1) = qi(1)*uinp(1,kkpole) + qi(4)*uinp(2,kkpole)
     &                       + qi(7)*uinp(3,kkpole)
               qiukp(2) = qi(2)*uinp(1,kkpole) + qi(5)*uinp(2,kkpole)
     &                       + qi(8)*uinp(3,kkpole)
               qiukp(3) = qi(3)*uinp(1,kkpole) + qi(6)*uinp(2,kkpole)
     &                       + qi(9)*uinp(3,kkpole)
               qkuip(1) = qk(1)*uinp(1,iipole) + qk(4)*uinp(2,iipole)
     &                       + qk(7)*uinp(3,iipole)
               qkuip(2) = qk(2)*uinp(1,iipole) + qk(5)*uinp(2,iipole)
     &                       + qk(8)*uinp(3,iipole)
               qkuip(3) = qk(3)*uinp(1,iipole) + qk(6)*uinp(2,iipole)
     &                       + qk(9)*uinp(3,iipole)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,iipole)*qkr(3) - uind(3,iipole)*qkr(2)
               uixqkr(2) = uind(3,iipole)*qkr(1) - uind(1,iipole)*qkr(3)
               uixqkr(3) = uind(1,iipole)*qkr(2) - uind(2,iipole)*qkr(1)
               ukxqir(1) = uind(2,kkpole)*qir(3) - uind(3,kkpole)*qir(2)
               ukxqir(2) = uind(3,kkpole)*qir(1) - uind(1,kkpole)*qir(3)
               ukxqir(3) = uind(1,kkpole)*qir(2) - uind(2,kkpole)*qir(1)
               uixqkrp(1) = uinp(2,iipole)*qkr(3)- uinp(3,iipole)*qkr(2)
               uixqkrp(2) = uinp(3,iipole)*qkr(1)- uinp(1,iipole)*qkr(3)
               uixqkrp(3) = uinp(1,iipole)*qkr(2)- uinp(2,iipole)*qkr(1)
               ukxqirp(1) = uinp(2,kkpole)*qir(3)- uinp(3,kkpole)*qir(2)
               ukxqirp(2) = uinp(3,kkpole)*qir(1)- uinp(1,kkpole)*qir(3)
               ukxqirp(3) = uinp(1,kkpole)*qir(2)- uinp(2,kkpole)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
c
               sci(1) = uind(1,iipole)*dk(1) + uind(2,iipole)*dk(2)
     &                     + uind(3,iipole)*dk(3) + di(1)*uind(1,kkpole)
     &                     + di(2)*uind(2,kkpole) + di(3)*uind(3,kkpole)
               sci(2) = uind(1,iipole)*uind(1,kkpole) + uind(2,iipole)*
     $                      uind(2,kkpole)+uind(3,iipole)*uind(3,kkpole)
               sci(3) = uind(1,iipole)*xr+uind(2,iipole)*yr+
     $           uind(3,iipole)*zr
               sci(4) = uind(1,kkpole)*xr+uind(2,kkpole)*yr+
     $           uind(3,kkpole)*zr
               sci(7) = qir(1)*uind(1,kkpole) + qir(2)*uind(2,kkpole)
     &                     + qir(3)*uind(3,kkpole)
               sci(8) = qkr(1)*uind(1,iipole) + qkr(2)*uind(2,iipole)
     &                     + qkr(3)*uind(3,iipole)
               scip(1) = uinp(1,iipole)*dk(1) + uinp(2,iipole)*dk(2)
     &                     + uinp(3,iipole)*dk(3) + di(1)*uinp(1,kkpole)
     &                     + di(2)*uinp(2,kkpole) + di(3)*uinp(3,kkpole)
               scip(2) =uind(1,iipole)*uinp(1,kkpole)+uind(2,iipole)*
     $                   uinp(2,kkpole)+ uind(3,iipole)*uinp(3,kkpole)+
     $                   uinp(1,iipole)*uind(1,kkpole)+ uinp(2,iipole)*
     $                   uind(2,kkpole)+uinp(3,iipole)*uind(3,kkpole)
               scip(3) = uinp(1,iipole)*xr+uinp(2,iipole)*yr+
     $                   uinp(3,iipole)*zr
               scip(4) = uinp(1,kkpole)*xr+uinp(2,kkpole)*yr+
     $                   uinp(3,kkpole)*zr
               scip(7) = qir(1)*uinp(1,kkpole) + qir(2)*uinp(2,kkpole)
     &                      + qir(3)*uinp(3,kkpole)
               scip(8) = qkr(1)*uinp(1,iipole) + qkr(2)*uinp(2,iipole)
     &                      + qkr(3)*uinp(3,iipole)
c
c     calculate the gl functions for permanent components
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
               gl(4) = sc(5)*sc(6)
               gl(5) = -4.0d0 * sc(9)
               gl(6) = sc(2)
               gl(7) = 2.0d0 * (sc(7)-sc(8))
               gl(8) = 2.0d0 * sc(10)
c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
     &                + bn(2)*(gl(2)+gl(7)+gl(8))
     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
               erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                  + rr5*(gl(2)+gl(7)+gl(8))
     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               erl = erl * (1.0d0-mscale(kglob))
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               e = e - erl
               ei = ei - erli
c
               e = f * e
               ei = f * ei
               emtt = emtt + e
               eptt = eptt + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
               if (molcule(iglob) .eq. molcule(kglob)) then
                  eintra = eintra + mscale(kglob)*erl*f
                  eintra = eintra + 0.5d0*pscale(kglob)
     &                        * (rr3*(gli(1)+gli(6))*scale3
     &                              + rr5*(gli(2)+gli(7))*scale5
     &                              + rr7*gli(3)*scale7)
               end if

c     set flags to compute components without screening
c
               dorl = .false.
               dorli = .false.
               if (mscale(kglob) .ne. 1.0d0)  dorl = .true.
               if (psc3 .ne. 0.0d0)  dorli = .true.
               if (dsc3 .ne. 0.0d0)  dorli = .true.
               if (usc3 .ne. 0.0d0)  dorli = .true.
c
c     zero out force and torque components without screening
c
               ftm2r = 0d0
               ftm2ri = 0d0
               ttm2r = 0d0
               ttm2ri = 0d0
               ttm3r = 0d0
               ttm3ri = 0d0
c
c     get the permanent force with screening
c
               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
               gf(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gf(4) = 2.0d0 * bn(2)
               gf(5) = 2.0d0 * (-ck*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
               gf(7) = 4.0d0 * bn(3)
               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))
c
c     get the permanent force without screening
c
               if (dorl) then
                  gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
     &                        + rr7*(gl(2)+gl(7)+gl(8))
     &                        + rr9*(gl(3)+gl(5)) + rr11*gl(4)
                  gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
                  gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
                  gfr(4) = 2.0d0 * rr5
                  gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
                  gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
                  gfr(7) = 4.0d0 * rr7
                  ftm2r(1) = gfr(1)*xr + gfr(2)*di(1) + gfr(3)*dk(1)
     &                          + gfr(4)*(qkdi(1)-qidk(1))
     &                          + gfr(5)*qir(1) + gfr(6)*qkr(1)
     &                          + gfr(7)*(qiqkr(1)+qkqir(1))
                  ftm2r(2) = gfr(1)*yr + gfr(2)*di(2) + gfr(3)*dk(2)
     &                          + gfr(4)*(qkdi(2)-qidk(2))
     &                          + gfr(5)*qir(2) + gfr(6)*qkr(2)
     &                          + gfr(7)*(qiqkr(2)+qkqir(2))
                  ftm2r(3) = gfr(1)*zr + gfr(2)*di(3) + gfr(3)*dk(3)
     &                          + gfr(4)*(qkdi(3)-qidk(3))
     &                          + gfr(5)*qir(3) + gfr(6)*qkr(3)
     &                          + gfr(7)*(qiqkr(3)+qkqir(3))
               end if
c
c     get the induced force with screening
c
               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     + 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,iipole)+uinp(1,iipole))
     &            + bn(2)*(sci(4)*uinp(1,iipole)+scip(4)*uind(1,iipole))
     &            + gfi(3)*(uind(1,kkpole)+uinp(1,kkpole))
     &            + bn(2)*(sci(3)*uinp(1,kkpole)+scip(3)*uind(1,kkpole))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,iipole)+uinp(2,iipole))
     &            + bn(2)*(sci(4)*uinp(2,iipole)+scip(4)*uind(2,iipole))
     &            + gfi(3)*(uind(2,kkpole)+uinp(2,kkpole))
     &            + bn(2)*(sci(3)*uinp(2,kkpole)+scip(3)*uind(2,kkpole))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,iipole)+uinp(3,iipole))
     &            + bn(2)*(sci(4)*uinp(3,iipole)+scip(4)*uind(3,iipole))
     &            + gfi(3)*(uind(3,kkpole)+uinp(3,kkpole))
     &            + bn(2)*(sci(3)*uinp(3,kkpole)+scip(3)*uind(3,kkpole))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               if (dorli) then
                  gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                               + (glip(1)+glip(6))*dsc3
     &                               + scip(2)*usc3)
     &                    + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                               + (glip(7)+glip(2))*dsc5
     &                        - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                    + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
                  gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
                  gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
                  gfri(4) = 2.0d0 * rr5
                  gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
                  gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
                  ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &          (- rr3*ck*(uind(1,iipole)*psc3+uinp(1,iipole)*dsc3)
     &          + rr5*sc(4)*(uind(1,iipole)*psc5+uinp(1,iipole)*dsc5)
     &          - rr7*sc(6)*(uind(1,iipole)*psc7+uinp(1,iipole)*dsc7))
     &          + (rr3*ci*(uind(1,kkpole)*psc3+uinp(1,kkpole)*dsc3)
     &          + rr5*sc(3)*(uind(1,kkpole)*psc5+uinp(1,kkpole)*dsc5)
     &          + rr7*sc(5)*(uind(1,kkpole)*psc7
     &          + uinp(1,kkpole)*dsc7))*0.5d0
     &          + rr5*usc5*(sci(4)*uinp(1,iipole)+scip(4)*uind(1,iipole)
     &               + sci(3)*uinp(1,kkpole)
     &               +scip(3)*uind(1,kkpole))*0.5d0
     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &               + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &               + (qkuip(1)-qiukp(1))*dsc5)
     &               + gfri(5)*qir(1) + gfri(6)*qkr(1)
                  ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &          (- rr3*ck*(uind(2,iipole)*psc3+uinp(2,iipole)*dsc3)
     &          + rr5*sc(4)*(uind(2,iipole)*psc5+uinp(2,iipole)*dsc5)
     &          - rr7*sc(6)*(uind(2,iipole)*psc7+uinp(2,iipole)*dsc7))
     &          + (rr3*ci*(uind(2,kkpole)*psc3+uinp(2,kkpole)*dsc3)
     &          + rr5*sc(3)*(uind(2,kkpole)*psc5+uinp(2,kkpole)*dsc5)
     &          + rr7*sc(5)*(uind(2,kkpole)*psc7
     &          + uinp(2,kkpole)*dsc7))*0.5d0
     &          + rr5*usc5*(sci(4)*uinp(2,iipole)+scip(4)*uind(2,iipole)
     &               + sci(3)*uinp(2,kkpole)
     &          + scip(3)*uind(2,kkpole))*0.5d0
     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &               + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &               + (qkuip(2)-qiukp(2))*dsc5)
     &               + gfri(5)*qir(2) + gfri(6)*qkr(2)
                  ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &          (- rr3*ck*(uind(3,iipole)*psc3+uinp(3,iipole)*dsc3)
     &           + rr5*sc(4)*(uind(3,iipole)*psc5+uinp(3,iipole)*dsc5)
     &           - rr7*sc(6)*(uind(3,iipole)*psc7+uinp(3,iipole)*dsc7))
     &           + (rr3*ci*(uind(3,kkpole)*psc3+uinp(3,kkpole)*dsc3)
     &           + rr5*sc(3)*(uind(3,kkpole)*psc5+uinp(3,kkpole)*dsc5)
     &           + rr7*sc(5)*(uind(3,kkpole)*psc7
     &          + uinp(3,kkpole)*dsc7))*0.5d0
     &          + rr5*usc5*(sci(4)*uinp(3,iipole)+scip(4)*uind(3,iipole)
     &               + sci(3)*uinp(3,kkpole)
     &               +scip(3)*uind(3,kkpole))*0.5d0
     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &               + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &               + (qkuip(3)-qiukp(3))*dsc5)
     &               + gfri(5)*qir(3) + gfri(6)*qkr(3)
               end if
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kglob)
     &                                +(glip(1)+glip(6))*dscale(kglob))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kglob)
     &                                +(glip(2)+glip(7))*dscale(kglob))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kglob)
     &                                +glip(3)*dscale(kglob))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(kglob) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(kglob)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     get the permanent torque with screening
c
               ttm2(1) = -bn(1)*dixdk(1) + gf(2)*dixr(1)
     &                      + gf(4)*(dixqkr(1)+dkxqir(1)
     &                              +rxqidk(1)-2.0d0*qixqk(1))
     &                      - gf(5)*rxqir(1)
     &                      - gf(7)*(rxqikr(1)+qkrxqir(1))
               ttm2(2) = -bn(1)*dixdk(2) + gf(2)*dixr(2)
     &                      + gf(4)*(dixqkr(2)+dkxqir(2)
     &                              +rxqidk(2)-2.0d0*qixqk(2))
     &                      - gf(5)*rxqir(2)
     &                      - gf(7)*(rxqikr(2)+qkrxqir(2))
               ttm2(3) = -bn(1)*dixdk(3) + gf(2)*dixr(3)
     &                      + gf(4)*(dixqkr(3)+dkxqir(3)
     &                              +rxqidk(3)-2.0d0*qixqk(3))
     &                      - gf(5)*rxqir(3)
     &                      - gf(7)*(rxqikr(3)+qkrxqir(3))
               ttm3(1) = bn(1)*dixdk(1) + gf(3)*dkxr(1)
     &                      - gf(4)*(dixqkr(1)+dkxqir(1)
     &                              +rxqkdi(1)-2.0d0*qixqk(1))
     &                      - gf(6)*rxqkr(1)
     &                      - gf(7)*(rxqkir(1)-qkrxqir(1))
               ttm3(2) = bn(1)*dixdk(2) + gf(3)*dkxr(2)
     &                      - gf(4)*(dixqkr(2)+dkxqir(2)
     &                              +rxqkdi(2)-2.0d0*qixqk(2))
     &                      - gf(6)*rxqkr(2)
     &                      - gf(7)*(rxqkir(2)-qkrxqir(2))
               ttm3(3) = bn(1)*dixdk(3) + gf(3)*dkxr(3)
     &                      - gf(4)*(dixqkr(3)+dkxqir(3)
     &                              +rxqkdi(3)-2.0d0*qixqk(3))
     &                      - gf(6)*rxqkr(3)
     &                      - gf(7)*(rxqkir(3)-qkrxqir(3))
c
c     get the permanent torque without screening
c
               if (dorl) then
                  ttm2r(1) = -rr3*dixdk(1) + gfr(2)*dixr(1)
     &                          + gfr(4)*(dixqkr(1)+dkxqir(1)
     &                                   +rxqidk(1)-2.0d0*qixqk(1))
     &                          - gfr(5)*rxqir(1)
     &                          - gfr(7)*(rxqikr(1)+qkrxqir(1))
                  ttm2r(2) = -rr3*dixdk(2) + gfr(2)*dixr(2)
     &                          + gfr(4)*(dixqkr(2)+dkxqir(2)
     &                                   +rxqidk(2)-2.0d0*qixqk(2))
     &                          - gfr(5)*rxqir(2)
     &                          - gfr(7)*(rxqikr(2)+qkrxqir(2))
                  ttm2r(3) = -rr3*dixdk(3) + gfr(2)*dixr(3)
     &                          + gfr(4)*(dixqkr(3)+dkxqir(3)
     &                                   +rxqidk(3)-2.0d0*qixqk(3))
     &                          - gfr(5)*rxqir(3)
     &                          - gfr(7)*(rxqikr(3)+qkrxqir(3))
                  ttm3r(1) = rr3*dixdk(1) + gfr(3)*dkxr(1)
     &                          - gfr(4)*(dixqkr(1)+dkxqir(1)
     &                                   +rxqkdi(1)-2.0d0*qixqk(1))
     &                          - gfr(6)*rxqkr(1)
     &                          - gfr(7)*(rxqkir(1)-qkrxqir(1))
                  ttm3r(2) = rr3*dixdk(2) + gfr(3)*dkxr(2)
     &                          - gfr(4)*(dixqkr(2)+dkxqir(2)
     &                                   +rxqkdi(2)-2.0d0*qixqk(2))
     &                          - gfr(6)*rxqkr(2)
     &                          - gfr(7)*(rxqkir(2)-qkrxqir(2))
                  ttm3r(3) = rr3*dixdk(3) + gfr(3)*dkxr(3)
     &                          - gfr(4)*(dixqkr(3)+dkxqir(3)
     &                                   +rxqkdi(3)-2.0d0*qixqk(3))
     &                          - gfr(6)*rxqkr(3)
     &                          - gfr(7)*(rxqkir(3)-qkrxqir(3))
               end if
c
c     get the induced torque with screening
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               ttm2i(1) = -0.5d0*bn(1)*(dixuk(1)+dixukp(1))
     &                       + gti(2)*dixr(1) - gti(5)*rxqir(1)
     &                       + 0.5d0*gti(4)*(ukxqir(1)+rxqiuk(1)
     &                                      +ukxqirp(1)+rxqiukp(1))
               ttm2i(2) = -0.5d0*bn(1)*(dixuk(2)+dixukp(2))
     &                       + gti(2)*dixr(2) - gti(5)*rxqir(2)
     &                       + 0.5d0*gti(4)*(ukxqir(2)+rxqiuk(2)
     &                                      +ukxqirp(2)+rxqiukp(2))
               ttm2i(3) = -0.5d0*bn(1)*(dixuk(3)+dixukp(3))
     &                       + gti(2)*dixr(3) - gti(5)*rxqir(3)
     &                       + 0.5d0*gti(4)*(ukxqir(3)+rxqiuk(3)
     &                                      +ukxqirp(3)+rxqiukp(3))
               ttm3i(1) = -0.5d0*bn(1)*(dkxui(1)+dkxuip(1))
     &                       + gti(3)*dkxr(1) - gti(6)*rxqkr(1)
     &                       - 0.5d0*gti(4)*(uixqkr(1)+rxqkui(1)
     &                                      +uixqkrp(1)+rxqkuip(1))
               ttm3i(2) = -0.5d0*bn(1)*(dkxui(2)+dkxuip(2))
     &                       + gti(3)*dkxr(2) - gti(6)*rxqkr(2)
     &                       - 0.5d0*gti(4)*(uixqkr(2)+rxqkui(2)
     &                                       +uixqkrp(2)+rxqkuip(2))
               ttm3i(3) = -0.5d0*bn(1)*(dkxui(3)+dkxuip(3))
     &                       + gti(3)*dkxr(3) - gti(6)*rxqkr(3)
     &                       - 0.5d0*gti(4)*(uixqkr(3)+rxqkui(3)
     &                                      +uixqkrp(3)+rxqkuip(3))
c
c     get the induced torque without screening
c
               if (dorli) then
                  gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
                  gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
                  gtri(4) = gfri(4)
                  gtri(5) = gfri(5)
                  gtri(6) = gfri(6)
                  ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &                           + gtri(2)*dixr(1) - gtri(5)*rxqir(1)
     &                           + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &                             +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0
                  ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &                           + gtri(2)*dixr(2) - gtri(5)*rxqir(2)
     &                           + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &                             +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0
                  ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &                           + gtri(2)*dixr(3) - gtri(5)*rxqir(3)
     &                           + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &                             +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0
                  ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &                           + gtri(3)*dkxr(1) - gtri(6)*rxqkr(1)
     &                           - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &                             +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0
                  ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &                           + gtri(3)*dkxr(2) - gtri(6)*rxqkr(2)
     &                           - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &                             +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0
                  ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &                           + gtri(3)*dkxr(3) - gtri(6)*rxqkr(3)
     &                           - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &                             +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0
               end if
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2(j) = f * (ftm2(j)-(1.0d0-mscale(kglob))*ftm2r(j))
                  ftm2i(j) = f *(ftm2i(j)-ftm2ri(j))
                  ttm2(j) = f * (ttm2(j)-(1.0d0-mscale(kglob))*ttm2r(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3(j) = f * (ttm3(j)-(1.0d0-mscale(kglob))*ttm3r(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
c
c     increment gradient due to force and torque on first site
c
               demtt(1,i) = demtt(1,i) + ftm2(1)
               demtt(2,i) = demtt(2,i) + ftm2(2)
               demtt(3,i) = demtt(3,i) + ftm2(3)
               deptt(1,i) = deptt(1,i) + ftm2i(1)
               deptt(2,i) = deptt(2,i) + ftm2i(2)
               deptt(3,i) = deptt(3,i) + ftm2i(3)
               call torque3 (iipole,ttm2,ttm2i,frcxi,frcyi,frczi,
     $          torquedirt,torquedirpt)
c
c     increment gradient due to force and torque on second site
c
               demtt(1,kbis) = demtt(1,kbis) - ftm2(1)
               demtt(2,kbis) = demtt(2,kbis) - ftm2(2)
               demtt(3,kbis) = demtt(3,kbis) - ftm2(3)
               deptt(1,kbis) = deptt(1,kbis) - ftm2i(1)
               deptt(2,kbis) = deptt(2,kbis) - ftm2i(2)
               deptt(3,kbis) = deptt(3,kbis) - ftm2i(3)
               call torque3 (kkpole,ttm3,ttm3i,frcxk,frcyk,frczk,
     $          torquedirt,torquedirpt)
c
c     increment the internal virial tensor components
c
c
               iaz = zaxis(iipole)
               if (iaz .le. 0)  iaz = iglob
               iax = xaxis(iipole)
               if (iax .le. 0)  iax = iglob
               iay = yaxis(iipole)
               if (iay .le. 0)  iay = iglob
               kaz = zaxis(kkpole)
               if (kaz .le. 0)  kaz = kglob
               kax = xaxis(kkpole)
               if (kax .le. 0)  kax = kglob
               kay = yaxis(kkpole)
               if (kay .le. 0)  kay = kglob
               xiz = x(iaz) - x(iglob)
               yiz = y(iaz) - y(iglob)
               ziz = z(iaz) - z(iglob)
               xix = x(iax) - x(iglob)
               yix = y(iax) - y(iglob)
               zix = z(iax) - z(iglob)
               xiy = x(iay) - x(iglob)
               yiy = y(iay) - y(iglob)
               ziy = z(iay) - z(iglob)
               xkz = x(kaz) - x(kglob)
               ykz = y(kaz) - y(kglob)
               zkz = z(kaz) - z(kglob)
               xkx = x(kax) - x(kglob)
               ykx = y(kax) - y(kglob)
               zkx = z(kax) - z(kglob)
               xky = x(kay) - x(kglob)
               yky = y(kay) - y(kglob)
               zky = z(kay) - z(kglob)
               vxx = -xr*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)
               virt(1,1) = virt(1,1) + vxx
               virt(2,1) = virt(2,1) + vyx
               virt(3,1) = virt(3,1) + vzx
               virt(1,2) = virt(1,2) + vyx
               virt(2,2) = virt(2,2) + vyy
               virt(3,2) = virt(3,2) + vzy
               virt(1,3) = virt(1,3) + vzx
               virt(2,3) = virt(2,3) + vzy
               virt(3,3) = virt(3,3) + vzz
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0d0
            pscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0d0
            pscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0d0
            pscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0d0
            pscale(i15(j,iglob)) = 1.0d0
         end do
         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = 1.0d0
            uscale(ip11(j,iglob)) = 1.0d0
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = 1.0d0
            uscale(ip12(j,iglob)) = 1.0d0
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = 1.0d0
            uscale(ip13(j,iglob)) = 1.0d0
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = 1.0d0
            uscale(ip14(j,iglob)) = 1.0d0
         end do
      end do
c
c     end OpenMP directives for the major loop structure
c
c$omp end parallel do
c
c     add local copies to global variables for OpenMP calculation
c
      em = emtt
      ep = eptt
      dem = demtt
      dep = deptt 
      torquedir = torquedirt 
      torquedirp = torquedirpt
      vir = virt
c
c
      deallocate (torquedirt)
      deallocate (torquedirpt)
      deallocate (demtt)
      deallocate (deptt)
      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      return
      end
