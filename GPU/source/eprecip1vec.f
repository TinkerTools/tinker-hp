c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine eprecip1  --  PME recip polarize energy & derivs  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "eprecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to dipole polarization
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
#include "tinker_macro.h"
      subroutine eprecip1vec
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use fft
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use timestat
      use tinheader ,only:ti_p,re_p
      use virial
      use mpi
      implicit none
      integer status(MPI_STATUS_SIZE),tag,ierr
      integer i,j,k,m,ii,iipole,iglob,iloc
      integer j1,j2,j3
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real(t_p) e,eterm,f
      real(t_p) r1,r2,r3
      real(t_p) h1,h2,h3
      real(t_p) f1,f2,f3
      real(t_p) vxx,vyy,vzz
      real(t_p) vxy,vxz,vyz
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) vterm,struc2
      real(t_p) trq(3),fix(3)
      real(t_p) fiy(3),fiz(3)
      real(t_p) cphim(4),cphid(4)
      real(t_p) cphip(4)
      real(t_p) a(3,3),ftc(10,10)
      real(t_p), allocatable :: cmp(:,:),fmp(:,:)
      real(t_p), allocatable :: fuind(:,:)
      real(t_p), allocatable :: fuinp(:,:)
      real(t_p), allocatable :: fphid(:,:)
      real(t_p), allocatable :: fphip(:,:)
      real(t_p), allocatable :: fphidp(:,:)
      real(t_p), allocatable :: qgrip(:,:,:,:)
      real(t_p), allocatable :: qgridmpi(:,:,:,:,:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer nprocloc,commloc,rankloc,proc
!DIR$ ATTRIBUTES ALIGN:64::fuind1vec,fuind2vec,fuind3vec
      real(t_p) fuind1vec(npolereclocloop)
      real(t_p) fuind2vec(npolereclocloop)
      real(t_p) fuind3vec(npolereclocloop)
      real(t_p) time0,time1
c
      if(rank.eq.0.and.tinkerdebug)write(*,*) 'eprecip1vec'
      call timer_enter( timer_eprecip )

      if (use_pmecore) then
        nprocloc = nrec
        commloc  = comm_rec
        rankloc  = rank_bis
      else
        nprocloc = nproc
        commloc  = COMM_TINKER
        rankloc  = rank
      end if
c
c     indices into the electrostatic field array
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      allocate (fuind(3,npolerecloc))
      allocate (fuinp(3,npolerecloc))
      allocate (cmp(10,npolerecloc))
      allocate (fmp(10,npolerecloc))
      allocate (fphid(10,npolerecloc))
      allocate (fphip(10,npolerecloc))
      allocate (fphidp(20,npolerecloc))
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0_ti_p
      vxy = 0.0_ti_p
      vxz = 0.0_ti_p
      vyy = 0.0_ti_p
      vyz = 0.0_ti_p
      vzz = 0.0_ti_p
c
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
c
c     initialize variables required for the scalar summation
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c
c     remove scalar sum virial from prior multipole 3-D FFT
c
c      if (use_mpole) then
         vxx = -vmxx
         vxy = -vmxy
         vxz = -vmxz
         vyy = -vmyy
         vyz = -vmyz
         vzz = -vmzz
cc
cc     compute the arrays of B-spline coefficients
cc
c      else
c         call bspline_fill
c         call table_fill
cc
cc     assign only the permanent multipoles to the PME grid
cc     and perform the 3-D FFT forward transformation
cc
         do i = 1, npolerecloc
            iipole = polerecglob(i)
            cmp(1,i) = rpole(1,iipole)
            cmp(2,i) = rpole(2,iipole)
            cmp(3,i) = rpole(3,iipole)
            cmp(4,i) = rpole(4,iipole)
            cmp(5,i) = rpole(5,iipole)
            cmp(6,i) = rpole(9,iipole)
            cmp(7,i) = rpole(13,iipole)
            cmp(8,i) = 2.0_ti_p * rpole(6,iipole)
            cmp(9,i) = 2.0_ti_p * rpole(7,iipole)
            cmp(10,i) = 2.0_ti_p * rpole(10,iipole)
           call cmp_to_fmp_site(cmp(1,i),fmp(1,i))
         end do
c         call cmp_to_fmp (cmp,fmp)
c         call grid_mpole (fmp)
c         call fftfront
cc
cc     make the scalar summation over reciprocal lattice
cc
c         do i = 1, ntot-1
c            k3 = i/nff + 1
c            j = i - (k3-1)*nff
c            k2 = j/nfft1 + 1
c            k1 = j - (k2-1)*nfft1 + 1
c            m1 = k1 - 1
c            m2 = k2 - 1
c            m3 = k3 - 1
c            if (k1 .gt. nf1)  m1 = m1 - nfft1
c            if (k2 .gt. nf2)  m2 = m2 - nfft2
c            if (k3 .gt. nf3)  m3 = m3 - nfft3
c            r1 = real(m1,t_p)
c            r2 = real(m2,t_p)
c            r3 = real(m3,t_p)
c            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
c            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
c            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
c            hsq = h1*h1 + h2*h2 + h3*h3
c            term = -pterm * hsq
c            expterm = 0.0_ti_p
c            if (term .gt. -50.0_ti_p) then
c               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
c               expterm = exp(term) / denom
c               if (.not. use_bounds) then
c                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
c               else if (octahedron) then
c                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0_ti_p
c               end if
c               struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
c               eterm = 0.5_ti_p * f * expterm * struc2
c               vterm = (2.0_ti_p/hsq) * (1.0_ti_p-term) * eterm
c               vxx = vxx - h1*h1*vterm + eterm
c               vxy = vxy - h1*h2*vterm
c               vxz = vxz - h1*h3*vterm
c               vyy = vyy - h2*h2*vterm + eterm
c               vyz = vyz - h2*h3*vterm
c               vzz = vzz - h3*h3*vterm + eterm
c            end if
c         end do
cc
cc     account for zeroth grid point for nonperiodic system
cc
c         qfac(1,1,1) = 0.0_ti_p
c         if (.not. use_bounds) then
c            expterm = 0.5_ti_p * pi / xbox
c            struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
c            e = f * expterm * struc2
c            qfac(1,1,1) = expterm
c         end if
cc
cc     complete the transformation of the PME grid
cc
c         do k = 1, nfft3
c            do j = 1, nfft2
c               do i = 1, nfft1
c                  term = qfac(i,j,k)
c                  qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
c                  qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
c               end do
c            end do
c         end do
cc
cc     perform 3-D FFT backward transform and get potential
cc
c         call fftback
c         call fphi_mpole (fphi)
c         do i = 1, npole
c            do j = 1, 20
c               fphi(j,i) = f * fphi(j,i)
c            end do
c         end do
c         call fphi_to_cphi (fphi,cphi)
c      end if
c
c     zero out the PME grid
c
      qgrid2in_2d = 0_ti_p
c
c     convert Cartesian induced dipoles to fractional coordinates
c
!DIR$ VECTOR ALIGNED
      do i = 1, 3
         a(1,i) = real(nfft1,t_p) * recip(i,1)
         a(2,i) = real(nfft2,t_p) * recip(i,2)
         a(3,i) = real(nfft3,t_p) * recip(i,3)
      end do
      time0 = mpi_wtime()
      do ii = 1, npolerecloc
         iipole = polerecglob(ii)
         iglob = ipole(iipole)
         do j = 1, 3
            fuind(j,ii) = a(j,1)*uind(1,iipole) + a(j,2)*uind(2,iipole)
     &                      + a(j,3)*uind(3,iipole)
            fuinp(j,ii) = a(j,1)*uinp(1,iipole) + a(j,2)*uinp(2,iipole)
     &                      + a(j,3)*uinp(3,iipole)
         end do
         call grid_uind_site(iglob,ii,fuind(1,ii),fuinp(1,ii),
     $    qgrid2in_2d)
      end do
      time1 = mpi_wtime()
      timegrid1 = timegrid1 + time1 - time0
c
c     MPI : begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_TPREC,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgrid2in_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_send(i),tag,commloc,reqsend(i),ierr)
      end do
c
      do i = 1, nrec_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
        qgrid2in_2d(:,:,:,:,1) = qgrid2in_2d(:,:,:,:,1)+
     $   qgridmpi(:,:,:,:,i)
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1 - time0
c
      time0 = mpi_wtime()
      call fft2d_frontmpi(qgrid2in_2d,qgrid2out_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1 - time0
c
c     account for zeroth grid point for nonperiodic system
c
      if (.not. use_bounds) then
         expterm = 0.5_ti_p * pi / xbox
         struc2 = qgrid2in_2d(1,1,1,1,1)**2 + qgrid2in_2d(2,1,1,1,1)**2
         e = f * expterm * struc2
         ep = ep + e
      end if
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     $   .and.(kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5_ti_p * pi / xbox
           struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
     $       qgrid2in_2d(2,1,1,1,1)**2
           e = f * expterm * struc2
           ep = ep + e
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
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1 - time0
c
c     MPI : Begin reception
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        tag = nprocloc*rankloc + prec_send(i) + 1
        call MPI_IRECV(qgrid2in_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_recep(i),tag,commloc,req2rec(i),ierr)
      end do
c
c     MPI : begin sending
c
      do i = 1, nrec_recep
        tag = nprocloc*prec_recep(i) + rankloc + 1
        call MPI_ISEND(qgrid2in_2d(1,1,1,1,1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_send(i),tag,commloc,req2send(i),ierr)
      end do
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        call MPI_WAIT(req2rec(i),status,ierr)
      end do
      do i = 1, nrec_recep
        call MPI_WAIT(req2send(i),status,ierr)
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1 - time0
       time0 = mpi_wtime()
       do ii = 1, npolerecloc
         iipole = polerecglob(ii)
         iglob = ipole(iipole)
         call fphi_uind_site(iglob,ii,fphid(1,ii),fphip(1,ii),
     $        fphidp(1,ii))
         do j = 1, 10
            fphid(j,ii) = electric * fphid(j,ii)
            fphip(j,ii) = electric * fphip(j,ii)
         end do
         do j = 1, 20
            fphidp(j,ii) = electric * fphidp(j,ii)
         end do
         do j = 1, 20
            fphirec(j,ii) = electric * fphirec(j,ii)
         end do
       end do
c
c     increment the induced dipole energy and gradient
c
      e = 0.0_ti_p
      do i = 1, npolerecloc
         f1 = 0.0_ti_p
         f2 = 0.0_ti_p
         f3 = 0.0_ti_p
         do k = 1, 3
            j1 = deriv1(k+1)
            j2 = deriv2(k+1)
            j3 = deriv3(k+1)
            e = e + fuind(k,i)*fphirec(k+1,i)
            f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphirec(j1,i)
     &              + fuind(k,i)*fphip(j1,i)
     &              + fuinp(k,i)*fphid(j1,i)
            f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphirec(j2,i)
     &              + fuind(k,i)*fphip(j2,i)
     &              + fuinp(k,i)*fphid(j2,i)
            f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphirec(j3,i)
     &              + fuind(k,i)*fphip(j3,i)
     &              + fuinp(k,i)*fphid(j3,i)
         end do
         do k = 1, 10
            f1 = f1 + fmp(k,i)*fphidp(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphidp(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphidp(deriv3(k),i)
         end do
         f1 = 0.5_ti_p * real(nfft1,t_p) * f1
         f2 = 0.5_ti_p * real(nfft2,t_p) * f2
         f3 = 0.5_ti_p * real(nfft3,t_p) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         ii = locrec1(iglob)
         deprec(1,ii) = deprec(1,ii) + h1
         deprec(2,ii) = deprec(2,ii) + h2
         deprec(3,ii) = deprec(3,ii) + h3
      end do
      e = 0.5_ti_p * e
      ep = ep + e
c
c     set the potential to be the induced dipole average
c
      do i = 1, npolerecloc
         do k = 1, 10
            fphidp(k,i) = 0.5_ti_p  * fphidp(k,i)
         end do
         call fphi_to_cphi_site(fphidp(1,i),cphirec(1,i))
      end do
c
c     distribute torques into the induced dipole gradient
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         trq(1) = cmp(4,i)*cphirec(3,i) - cmp(3,i)*cphirec(4,i)
     &               + 2.0_ti_p*(cmp(7,i)-cmp(6,i))*cphirec(10,i)
     &               + cmp(9,i)*cphirec(8,i) + cmp(10,i)*cphirec(6,i)
     &               - cmp(8,i)*cphirec(9,i) - cmp(10,i)*cphirec(7,i)
         trq(2) = cmp(2,i)*cphirec(4,i) - cmp(4,i)*cphirec(2,i)
     &               + 2.0_ti_p*(cmp(5,i)-cmp(7,i))*cphirec(9,i)
     &               + cmp(8,i)*cphirec(10,i) + cmp(9,i)*cphirec(7,i)
     &               - cmp(9,i)*cphirec(5,i) - cmp(10,i)*cphirec(8,i)
         trq(3) = cmp(3,i)*cphirec(2,i) - cmp(2,i)*cphirec(3,i)
     &               + 2.0_ti_p*(cmp(6,i)-cmp(5,i))*cphirec(8,i)
     &               + cmp(8,i)*cphirec(5,i) + cmp(10,i)*cphirec(9,i)
     &               - cmp(8,i)*cphirec(6,i) - cmp(9,i)*cphirec(10,i)

         call torque_rec (iipole,trq,fix,fiy,fiz,deprec)
      end do
c
c     induced dipole contribution to the internal virial
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         do j = 2, 4
            cphim(j) = 0.0_ti_p
            cphid(j) = 0.0_ti_p
            cphip(j) = 0.0_ti_p
            do k = 2, 4
               cphim(j) = cphim(j) + ftc(j,k)*fphirec(k,i)
               cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
               cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
            end do
         end do
         vxx = vxx - cphirec(2,i)*cmp(2,i)
     &         - 0.5_ti_p*(cphim(2)*(uind(1,iipole)+uinp(1,iipole))
     &         +cphid(2)*uinp(1,iipole)+cphip(2)*uind(1,iipole))
         vxy = vxy - 0.5_ti_p*(cphirec(2,i)*cmp(3,i)+cphirec(3,i)*
     $         cmp(2,i))
     &         - 0.25_ti_p*(cphim(2)*(uind(2,iipole)+uinp(2,iipole))
     &         +cphim(3)*(uind(1,iipole)+uinp(1,iipole))
     &         +cphid(2)*uinp(2,iipole)+cphip(2)*uind(2,iipole)
     &         +cphid(3)*uinp(1,iipole)+cphip(3)*uind(1,iipole))
         vxz = vxz - 0.5_ti_p*(cphirec(2,i)*cmp(4,i)+cphirec(4,i)*
     $         cmp(2,i))
     &         - 0.25_ti_p*(cphim(2)*(uind(3,iipole)+uinp(3,iipole))
     &         +cphim(4)*(uind(1,iipole)+uinp(1,iipole))
     &         +cphid(2)*uinp(3,iipole)+cphip(2)*uind(3,iipole)
     &         +cphid(4)*uinp(1,iipole)+cphip(4)*uind(1,iipole))
         vyy = vyy - cphirec(3,i)*cmp(3,i)
     &         - 0.5_ti_p*(cphim(3)*(uind(2,iipole)+uinp(2,iipole))
     &         +cphid(3)*uinp(2,iipole)+cphip(3)*uind(2,iipole))
         vyz = vyz - 0.5_ti_p*(cphirec(3,i)*cmp(4,i)+cphirec(4,i)*
     $         cmp(3,i))
     &         - 0.25_ti_p*(cphim(3)*(uind(3,iipole)+uinp(3,iipole))
     &         +cphim(4)*(uind(2,iipole)+uinp(2,iipole))
     &         +cphid(3)*uinp(3,iipole)+cphip(3)*uind(3,iipole)
     &         +cphid(4)*uinp(2,iipole)+cphip(4)*uind(2,iipole))
         vzz = vzz - cphirec(4,i)*cmp(4,i)
     &         - 0.5_ti_p*(cphim(4)*(uind(3,iipole)+uinp(3,iipole))
     &         +cphid(4)*uinp(3,iipole)+cphip(4)*uind(3,iipole))
         vxx = vxx - 2.0_ti_p*cmp(5,i)*cphirec(5,i) - cmp(8,i)*
     $         cphirec(8,i)
     &         - cmp(9,i)*cphirec(9,i)
         vxy = vxy - (cmp(5,i)+cmp(6,i))*cphirec(8,i)
     &         - 0.5_ti_p*(cmp(8,i)*(cphirec(6,i)+cphirec(5,i))
     &         +cmp(9,i)*cphirec(10,i)+cmp(10,i)*cphirec(9,i))
         vxz = vxz - (cmp(5,i)+cmp(7,i))*cphirec(9,i)
     &         - 0.5_ti_p*(cmp(9,i)*(cphirec(5,i)+cphirec(7,i))
     &          +cmp(8,i)*cphirec(10,i)+cmp(10,i)*cphirec(8,i))
         vyy = vyy - 2.0_ti_p*cmp(6,i)*cphirec(6,i) - cmp(8,i)*
     $         cphirec(8,i)
     &         - cmp(10,i)*cphirec(10,i)
         vyz = vyz - (cmp(6,i)+cmp(7,i))*cphirec(10,i)
     &         - 0.5_ti_p*(cmp(10,i)*(cphirec(6,i)+cphirec(7,i))
     &         +cmp(8,i)*cphirec(9,i)+cmp(9,i)*cphirec(8,i))
         vzz = vzz - 2.0_ti_p*cmp(7,i)*cphirec(7,i) -
     $             cmp(9,i)*cphirec(9,i)
     &             - cmp(10,i)*cphirec(10,i) 
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fphid)
      deallocate (fphip)
      deallocate (fphidp)
c
c     perform dynamic allocation of some local arrays
c
      allocate (qgrip(2,isize2(rankloc+1),jsize2(rankloc+1),
     $ ksize2(rankloc+1)))
c
c     assign permanent and induced multipoles to the PME grid
c     and perform the 3-D FFT forward transformation
c
c
c    zero out the grid
c
      qgridin_2d = 0_ti_p
c
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         do j = 2, 4
            cmp(j,i) = cmp(j,i) + uinp(j-1,iipole)
         end do
        call cmp_to_fmp_site (cmp(1,i),fmp(1,i))
        call grid_mpole_site(iglob,i,fmp(1,i))
      end do
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_TPREC,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_send(i),tag,commloc,reqsend(i),ierr)
      end do
c
      do i = 1, nrec_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1)+
     $   qgridmpi(:,:,:,:,i)
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1 - time0
c
      time0 = mpi_wtime()
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
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
      qgridin_2d = 0_ti_p
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         do j = 2, 4
            cmp(j,i) = cmp(j,i) + uind(j-1,iipole) - uinp(j-1,iipole)
         end do
         call cmp_to_fmp_site(cmp(1,i),fmp(1,i))
         call grid_mpole_site(iglob,i,fmp(1,i))
      end do
      time1 = mpi_wtime()
      timegrid1 = timegrid1 + time1 - time0
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_TPREC,prec_recep(i),tag,
     $   commloc,reqrec(i),ierr)
      end do
c
c     MPI : begin sending
c
      time0 = mpi_wtime()
      do i = 1, nrec_send
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_TPREC,
     $   prec_send(i),tag,commloc,reqsend(i),ierr)
      end do
c
      do i = 1, nrec_recep
        call MPI_WAIT(reqrec(i),status,ierr)
      end do
      do i = 1, nrec_send
        call MPI_WAIT(reqsend(i),status,ierr)
      end do
c
c     do the reduction 'by hand'
c
      do i = 1, nrec_recep
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1)+
     $   qgridmpi(:,:,:,:,i)
      end do
      time1 = mpi_wtime()
      timerecreccomm = timerecreccomm + time1 - time0
c
      time0 = mpi_wtime()
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
      time1 = mpi_wtime()
      timeffts = timeffts + time1 - time0
c
c     make the scalar summation over reciprocal lattice
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.
     $   (kstart2(rankloc+1).eq.1)) then
           qfac_2d(1,1,1) = 0.0_ti_p
      end if
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
           r1 = real(m1,t_p)
           r2 = real(m2,t_p)
           r3 = real(m3,t_p)
           h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
           h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
           h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
           hsq = h1*h1 + h2*h2 + h3*h3
           term = -pterm * hsq
           expterm = 0.0_ti_p
           if (term .gt. -50.0_ti_p) then
              denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
              expterm = exp(term) / denom
              if (.not. use_bounds) then
                 expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
              else if (octahedron) then
                 if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0_ti_p
              end if
              struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $         k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)*
     $         qgrip(1,k1-istart2(rankloc+1)+1,
     $         k2-jstart2(rankloc+1)+1,
     $         k3-kstart2(rankloc+1)+1)
     &         + qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-
     $         jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)*
     $         qgrip(2,k1-istart2(rankloc+1)+1,
     $         k2-jstart2(rankloc+1)+1,
     $         k3-kstart2(rankloc+1)+1)
              eterm = 0.5_ti_p * f * expterm * struc2
              vterm = (2.0_ti_p/hsq) * (1.0_ti_p-term) * eterm
              vxx = vxx + h1*h1*vterm - eterm
              vxy = vxy + h1*h2*vterm
              vxz = vxz + h1*h3*vterm
              vyy = vyy + h2*h2*vterm - eterm
              vyz = vyz + h2*h3*vterm
              vzz = vzz + h3*h3*vterm - eterm
           end if
           qfac_2d(k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,k3-
     $       kstart2(rankloc+1)+1) = expterm
 20         continue
         end do
       end do
      end do
c
c     increment the internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vyz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz
c
c     perform deallocation of some local arrays
c
      deallocate (qgrip)
c      deallocate (fuind)
c      deallocate (fuinp)
      deallocate (cmp)
      deallocate (fmp)
c      deallocate (fphid)
c      deallocate (fphip)
c      deallocate (fphidp)
      deallocate (qgridmpi)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      call timer_exit( timer_eprecip )
      end
