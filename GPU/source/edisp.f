c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine edisp  --  damped dispersion energy              ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "edisp" calculates the dispersion energy; also partitions
c     the energy among the atoms
c
c     literature reference:
c
c     J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, "A Physically
c     Grounded Damped Dispersion Model with Particle Mesh Ewald
c     Summation", Journal of Chemical Physics, 149, 084115 (2018)
c
c
      subroutine edisp
      use atoms
      use domdec
      use dsppot
      use energi
      use ewald
      use inform
      use iounit
      use potent
      implicit none
      real*8 elrc
      character*11 mode
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_dewald) then
        call edisp0d
      else
        call edisp0b
      end if
c
c     apply long range dispersion correction if desired
c
      if (use_dcorr .and. .not.use_dewald) then
         mode = 'DISP'
         call evcorr (mode,elrc)
         edsp = edsp + elrc
         if (verbose .and. elrc.ne.0.0d0) then
            if (digits .ge. 8) then
               write (iout,10)  elrc
   10          format (/,' Long-Range Dispersion :',9x,f16.8)
            else if (digits .ge. 6) then
               write (iout,20)  elrc
   20          format (/,' Long-Range Dispersion :',9x,f16.6)
            else
               write (iout,30)  elrc
   30          format (/,' Long-Range Dispersion :',9x,f16.4)
            end if
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine edisp0b  --  damp dispersion energy via list    ##
c     ##                                                             ##
c     #################################################################
c
c
c     "edisp0b" calculates the damped dispersion potential energy
c     using a pairwise neighbor list
c
c
      subroutine edisp0b
      use atmlst
      use atmtyp
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use cutoff
      use disp
      use domdec
      use dsppot
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mutant
      use neigh
      use potent
      use shunt
      use usage
      implicit none
      integer i,j,k,iglob,kglob,kbis,nnvlst
      integer ii,iidisp,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ai2
      real*8 ak,ak2
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,taper
      real*8, allocatable :: dspscale(:)
      real*8 s,ds,dispshortcut2
      logical proceed,usei
      logical muti,mutk
      logical header,huge
      logical testcut,shortrange,longrange,fullrange
      character*11 mode
      character*80 :: RoutineName
c
c     compute the short, long, or full real space part of the Ewald summation
      shortrange = use_dispshort
      longrange  = use_displong
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'edispshort0b'
         mode        = 'SHORTDISP'
      else if (longrange) then
         RoutineName = 'edisplong0b'
         mode        = 'DISP'
      else
         RoutineName = 'edisp0b'
         mode        = 'DISP'
      endif
c
c     zero out the dispersion interaction
c
      edsp = 0.0d0
      if (ndisp .eq. 0) return
c
c     perform dynamic allocation of some local arrays
c
      allocate (dspscale(n))
c
c     initialize connected atom exclusion coefficients
c
      dspscale = 1.0d0
c
c     set cutoff and switching coefficients
c
      call switch (mode)
      dispshortcut2 = (dispshortcut - shortheal) ** 2
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. ndisp.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dispersion Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     find the damped dispersion energy via neighbor list search
c
      do ii = 1, ndisplocnl
         iidisp = dispglobnl(ii)
         iglob = idisp(iidisp)
         i = loc(iglob)
         ci = csix(iidisp)
         ai = adisp(iidisp)
         usei = use(iglob)
         muti = mut(iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            dspscale(i12(j,iglob)) = dsp2scale
         end do
         do j = 1, n13(iglob)
            dspscale(i13(j,iglob)) = dsp3scale
         end do
         do j = 1, n14(iglob)
            dspscale(i14(j,iglob)) = dsp4scale
         end do
         do j = 1, n15(iglob)
            dspscale(i15(j,iglob)) = dsp5scale
         end do
c
c     decide whether to compute the current interaction
c
         nnvlst = merge(nshortvlst(ii),
     &                  nvlst     (ii),
     &                  shortrange
     &                 )
         do kkk = 1, nnvlst
            kk = merge(shortvlst(kkk,ii),
     &                    vlst     (kkk,ii),
     &                    shortrange
     &                   )
            kglob = idisp(kk)
            kbis = loc(kglob)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            ck = csix(kk)
            ak = adisp(kk)
            mutk = mut(kglob)
            proceed = (usei .or. use(kglob))
            if (proceed) then
               xr = xi - x(kglob)
               yr = yi - y(kglob)
               zr = zi - z(kglob)
               if (use_bounds) call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               testcut = merge(r2 .le. off2.and.r2.ge.dispshortcut2,
     &                      r2 .le. off2,
     &                      longrange
     &                     )
               if (testcut) then
                  r = sqrt(r2)
                  r6 = r2**3
                  e = -ci * ck / r6
c
c     find the damping factor for the dispersion interaction
c
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (ai .ne. ak) then
                     ai2 = ai * ai
                     ak2 = ak * ak
                     dk2 = dk * dk
                     dk3 = dk * dk2
                     ti = ak2 / (ak2-ai2)
                     ti2 = ti * ti
                     tk = ai2 / (ai2-ak2)
                     tk2 = tk * tk
                     damp3 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2)*expk
     &                          - 2.0d0*ti2*tk*(1.0d0+di)*expi
     &                          - 2.0d0*tk2*ti*(1.0d0+dk)*expk
                     damp5 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2
     &                                       +di3/6.0d0)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2
     &                                    +dk3/6.0d0)*expk
     &                          - 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi
     &                          - 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  e = e * damp**2 * dspscale(kglob)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        e = e * vlambda
                     else if (.not.muti .or. .not.mutk) then
                        e = e * vlambda
                     end if
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if(longrange.or.fullrange) then
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        e = e * taper
                     end if
                  end if

                  if(shortrange .or. longrange)
     &            call switch_respa(r,dispshortcut,shortheal,s,ds)

                  if(shortrange) then
                     e  =   e * s
                  else if(longrange) then
                     e  =   e * (1.0d0 - s)
                  endif
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
                  edsp = edsp + e
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 4.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Dispersion',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,e
   30                format (' Disper',4x,2(i7,'-',a3),
     &                          9x,f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            dspscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            dspscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            dspscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            dspscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c                                                                        
c     perform deallocation of some local arrays
c
      deallocate (dspscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine edisp0d  --  Ewald dispersion energy via loop    ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "edisp3c" calculates the dispersion interaction energy using
c     particle mesh Ewald summation and a double loop
c
c
      subroutine edisp0d
      use analyz
      use atmlst
      use atoms
      use disp
      use domdec
      use energi
      use ewald
      use pme
      use potent
      implicit none
      integer i,ii,iglob,iidisp
      real*8 e
c
c
c     zero out the damped dispersion energy
c
      edsp = 0.0d0
      if (ndisp .eq. 0)  return
c
c     set Ewald coefficient
c
      aewald = adewald
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        if (use_disprec) then
          call edrecip
        end if
      end if
c
c     compute the real space portion of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_dispreal) then
          call edreal0d
        end if
c
c     compute the self-energy portion of the Ewald summation
c
        if (use_dispself) then
          do ii = 1, ndisploc
             iidisp = dispglob(ii)
             iglob = idisp(iidisp)
             i = loc(iglob)
             e = csix(iidisp)**2 * aewald**6 / 12.0d0
             edsp = edsp + e
          end do
        end if
      end if
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine edrecip  --  PME recip space damped dispersion  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "edrecip" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to damped dispersion
c
c
      subroutine edrecip
      use atmlst
      use boxes
      use bound
      use disp
      use domdec
      use energi
      use ewald
      use fft
      use math
      use mpi
      use pme
      use potent
      implicit none
      integer i,iidisp,iglob
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff
      integer nf1,nf2,nf3
      real*8 e,denom
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 term,expterm
      real*8 hsq,struc2
      real*8 h,hhh,b,bfac
      real*8 term1,denom0
      real*8 fac1,fac2,fac3
      real*8 erfcterm
      integer, allocatable :: req(:),reqbcast(:)
      integer status(MPI_STATUS_SIZE),tag,ierr,proc
      integer nprocloc,rankloc,commloc
      real*8, allocatable :: qgridmpi(:,:,:,:,:)
c
      if (use_pmecore) then
        nprocloc = nrec
        rankloc = rank_bis
        commloc = comm_rec
      else
        nprocloc = nproc
        rankloc = rank
        commloc = COMM_TINKER
      end if
c
c     set Ewald coefficient
c
      aewald = adewald
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     dynamic allocation of local arrays
c
      allocate (qgridmpi(2,n1mpimax,n2mpimax,n3mpimax,nrec_recep))
      allocate (req(nproc*nproc))
      allocate (reqbcast(nproc*nproc))
c
c
      do i = 1, ndisprecloc
        iidisp = disprecglob(i)
        iglob = idisp(iidisp)
        call bspline_fill_site(iglob,i)
      end do
c
      qgridin_2d = 0d0
c
c     MPI : Begin reception
c
      do i = 1, nrec_recep
        tag = nprocloc*rankloc + prec_recep(i) + 1
        call MPI_IRECV(qgridmpi(1,1,1,1,i),2*n1mpimax*n2mpimax*
     $   n3mpimax,MPI_REAL8,prec_recep(i),tag,commloc,req(tag),
     $   ierr)
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
      do i = 1, ndisprecloc
        iidisp = disprecglob(i)
        iglob = idisp(iidisp)
        call grid_disp_site(iglob,i)
      enddo
c
c     MPI : begin sending
c
      do i = 1, nrec_send
        proc = prec_send(i)
        tag = nprocloc*prec_send(i) + rankloc + 1
        call MPI_ISEND(qgridin_2d(1,1,1,1,i+1),
     $   2*n1mpimax*n2mpimax*n3mpimax,MPI_REAL8,proc,tag,
     $   commloc,req(tag),ierr)
      end do
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
        qgridin_2d(:,:,:,:,1) = qgridin_2d(:,:,:,:,1) + 
     $    qgridmpi(:,:,:,:,i) 
      end do
c
c     perform the 3-D FFT forward transformation
c
      call fft2d_frontmpi(qgridin_2d,qgridout_2d,n1mpimax,n2mpimax,
     $ n3mpimax)
c
c     use scalar sum to get the reciprocal space energy
c
      bfac = pi / aewald
      fac1 = 2.0d0 * pi**(3.5d0)
      fac2 = aewald**3
      fac3 = -2.0d0 * aewald * pi**2
      denom0 = (6.0d0*volbox) / (pi**1.5d0)
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      nff = nfft1 * nfft2
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
            h = sqrt(hsq)
            b = h * bfac
            hhh = h * hsq
            term = -b * b
            expterm = 0.0d0
            if (term .gt. -50.0d0) then
               denom = denom0 * bsmod1(k1) * bsmod2(k2) * bsmod3(k3)
               expterm = exp(term)
               erfcterm = erfc(b)
               if (.not. use_bounds) then
                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
                  erfcterm = erfcterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
                  if (mod(m1+m2+m3,2) .ne. 0)  erfcterm = 0.0d0
               end if
               term1 = fac1*erfcterm*hhh + expterm*(fac2+fac3*hsq)
               struc2 = qgridout_2d(1,k1-istart2(rankloc+1)+1,
     $  k2-jstart2(rankloc+1)+1,k3-kstart2(rankloc+1)+1)**2 + 
     $  qgridout_2d(2,k1-istart2(rankloc+1)+1,k2-jstart2(rankloc+1)+1,
     $  k3-kstart2(rankloc+1)+1)**2
               e = -(term1/denom) * struc2
               edsp = edsp + e
            end if
 10         continue
          end do
        end do
      end do
c
c     account for the total energy correction term
c
      if (rankloc.eq.0) then
        e = -csixpr * aewald**3 / denom0
        edsp = edsp + e
      end if

      deallocate (qgridmpi)
      deallocate (req)
      deallocate (reqbcast)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine edreal0d  --  real space disp energy via list    ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "edreal0d" evaluated the real space portion of the damped
c     dispersion energy using Ewald and a neighbor list
c
c
      subroutine edreal0d
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use cutoff
      use disp
      use domdec
      use dsppot
      use energi
      use ewald
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mutant
      use neigh
      use potent
      use shunt
      use usage
      implicit none
      integer i,j,iglob,kglob,kbis,nnvlst
      integer ii,iidisp,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r6
      real*8 ai,ai2
      real*8 ak,ak2
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 ralpha2
      real*8 expa,term
      real*8 damp3,damp5
      real*8 damp,scale
      real*8, allocatable :: dspscale(:)
      real*8 s,ds,dispshortcut2
      logical proceed,usei
      logical muti,mutk
      logical header
      logical testcut,shortrange,longrange,fullrange
      character*11 mode
      character*80 :: RoutineName

c     compute the short, long, or full real space part of the Ewald summation
      shortrange = use_dispshortreal
      longrange  = use_displong
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'edrealshort0d'
         mode        = 'SHORTDEWALD'
      else if (longrange) then
         RoutineName = 'edreallong0d'
         mode        = 'DEWALD'
      else
         RoutineName = 'edreal0d'
         mode        = 'DEWALD'
      endif
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (dspscale(n))
c
c     initialize connected atom exclusion coefficients
c
      dspscale = 1.0d0

      call switch (mode)
      dispshortcut2 = (dispshortcut - shortheal) ** 2
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. ndisp.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dispersion Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisplocnl
         iidisp = dispglobnl(ii)
         iglob = idisp(iidisp)
         i = loc(iglob)
         ci = csix(iidisp)
         ai = adisp(iidisp)
         usei = use(iglob)
         muti = mut(iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            dspscale(i12(j,iglob)) = dsp2scale
         end do
         do j = 1, n13(iglob)
            dspscale(i13(j,iglob)) = dsp3scale
         end do
         do j = 1, n14(iglob)
            dspscale(i14(j,iglob)) = dsp4scale
         end do
         do j = 1, n15(iglob)
            dspscale(i15(j,iglob)) = dsp5scale
         end do
         nnvlst = merge(nshortvlst(ii),
     &                  nvlst     (ii),
     &                  shortrange
     &                 )
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nnvlst
            kk = merge(shortvlst(kkk,ii),
     &                    vlst     (kkk,ii),
     &                    shortrange
     &                   )
            kglob = idisp(kk)
            kbis = loc(kglob)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            ck = csix(kk)
            ak = adisp(kk)
            mutk = mut(kglob)
            proceed = (usei .or. use(kglob))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(kglob)
               yr = yi - y(kglob)
               zr = zi - z(kglob)
               if (use_bounds) call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               testcut = merge(r2 .le. off2.and.r2.ge.dispshortcut2,
     &                      r2 .le. off2,
     &                      longrange
     &                     )
               if (testcut) then
                  r = sqrt(r2)
                  r6 = r2**3
                  ralpha2 = r2 * aewald**2
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  expa = exp(-ralpha2) * term
c
c     find the damping factor for the dispersion interaction
c
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (ai .ne. ak) then
                     ai2 = ai * ai
                     ak2 = ak * ak
                     dk2 = dk * dk
                     dk3 = dk * dk2
                     ti = ak2 / (ak2-ai2)
                     ti2 = ti * ti
                     tk = ai2 / (ai2-ak2)
                     tk2 = tk * tk
                     damp3 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2)*expk
     &                          - 2.0d0*ti2*tk*(1.0d0+di)*expi
     &                          - 2.0d0*tk2*ti*(1.0d0+dk)*expk
                     damp5 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2
     &                                       +di3/6.0d0)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2
     &                                    +dk3/6.0d0)*expk
     &                          - 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi
     &                          - 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  scale = dspscale(kglob) * damp**2
                  if (use_group)  scale = scale * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        scale = scale * vlambda
                     else if (.not.muti .or. .not.mutk) then
                        scale = scale * vlambda
                     end if
                  end if

                  scale = scale - 1.0d0
                  e = -ci * ck * (expa+scale) / r6

                  if(shortrange .or. longrange)
     &            call switch_respa(r,dispshortcut,shortheal,s,ds)

                  if(shortrange) then
                     e  =   e * s
                  else if(longrange) then
                     e  =   e * (1.0d0 - s)
                  endif
c     
c     increment the overall dispersoin energy components
c
                   edsp = edsp + e
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            dspscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            dspscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            dspscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            dspscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c                                                                        
c     perform deallocation of some local arrays
c
      deallocate (dspscale)
      return
      end
