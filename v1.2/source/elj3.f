c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine elj3  --  Lennard-Jones vdw energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "elj3" calculates the Lennard-Jones 6-12 van der Waals energy
c     and also partitions the energy among the atoms
c
c
      subroutine elj3
      use analyz
      use atoms
      use domdec
      use energi
      use inform
      use iounit
      use potent
      use vdwpot
      implicit none
      integer i
      real*8 elrc,aelrc
c
c     choose the method for summing over pairwise interactions
c
      if (use_vdwshort) then
        call eljshort3c
      else if (use_vdwlong) then
        call eljlong3c
      else
        call elj3c
      end if
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr (elrc)
         ev = ev + elrc
         aelrc = elrc / dble(n)
         do i = 1, nbloc
            aev(i) = aev(i) + aelrc
         end do
         if (verbose .and. elrc.ne.0.0d0) then
            write (iout,10)  elrc
   10       format (/,' Long Range vdw Correction :',9x,f12.4)
         end if
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine elj3c  --  Lennard-Jones analysis via list  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "elj3c" calculates the Lennard-Jones van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine elj3c
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use couple
      use domdec
      use energi
      use inter
      use inform
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,iglob,iivdw,kglob,kbis
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e
      real*8 p6,p12,eps
      real*8 rv,rdn
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical header,huge
      character*10 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      aev = 0d0
      header = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      vscale = 1.0d0
      iv14 = 0
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
      do ii = 1, nvdwbloc
         iivdw = vdwglob(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         rdn = kred(iglob)
         xr = x(iglob) - x(iv)
         yr = y(iglob) - y(iv)
         zr = z(iglob) - z(iv)
         if (use_polymer) call image(xr,yr,zr)
         xred(i) = rdn*xr + x(iv)
         yred(i) = rdn*yr + y(iv)
         zred(i) = rdn*zr + z(iv)
      end do
c
c     find the van der Waals energy via neighbor list search
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         it = jvdw(iglob)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(iglob) .or. use(iv))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            vscale(i12(j,iglob)) = v2scale
         end do
         do j = 1, n13(iglob)
            vscale(i13(j,iglob)) = v3scale
         end do
         do j = 1, n14(iglob)
            vscale(i14(j,iglob)) = v4scale
            iv14(i14(j,iglob)) = iglob
         end do
         do j = 1, n15(iglob)
            vscale(i15(j,iglob)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvlst(ii)
            kglob = vlst(kk,ii)
            kbis = loc(kglob)
            kv = ired(kglob)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(kglob) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(kglob)
               xr = xi - xred(kbis)
               yr = yi - yred(kbis)
               zr = zi - zred(kbis)
               if (use_bounds) call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(kglob) .eq. iglob) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(kglob)
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     increment the overall van der Waals energy components
c
                  if (e .ne. 0.0d0) then
                     nev = nev + 1
                     ev = ev + e
                     aev(i) = aev(i) + 0.5d0*e
                     aev(kbis) = aev(kbis) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             20x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,20) iglob,name(iglob),kglob,name(kglob)
     &                                ,rv,sqrt(rik2),e
   20                format (' VDW-LJ',4x,2(i7,'-',a3),
     &                          13x,2f10.4,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            vscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            vscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            vscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            vscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine eljshort3c  --  short range Lennard-Jones analysis via list  ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "eljshort3c" calculates the short range Lennard-Jones van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine eljshort3c
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use couple
      use cutoff
      use domdec
      use energi
      use inter
      use inform
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,iglob,iivdw,kglob,kbis
      integer ii,iv,it
      integer kk,kv,kt
      integer nevt
      integer, allocatable :: iv14(:)
      real*8 e
      real*8 p6,p12,eps
      real*8 rv,rdn
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2

      real*8 s,ds
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical header,huge
      character*10 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      aev = 0d0
      header = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      vscale = 1.0d0
      iv14 = 0
c
c     set the coefficients for the switching function
c
      mode = 'SHORTVDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
      do ii = 1, nvdwbloc
         iivdw = vdwglob(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         rdn = kred(iglob)
         xr = x(iglob) - x(iv)
         yr = y(iglob) - y(iv)
         zr = z(iglob) - z(iv)
         if (use_polymer) call image(xr,yr,zr)
         xred(i) = rdn*xr + x(iv)
         yred(i) = rdn*yr + y(iv)
         zred(i) = rdn*zr + z(iv)
      end do
c
c     find the van der Waals energy via neighbor list search
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         it = jvdw(iglob)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(iglob) .or. use(iv))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            vscale(i12(j,iglob)) = v2scale
         end do
         do j = 1, n13(iglob)
            vscale(i13(j,iglob)) = v3scale
         end do
         do j = 1, n14(iglob)
            vscale(i14(j,iglob)) = v4scale
            iv14(i14(j,iglob)) = iglob
         end do
         do j = 1, n15(iglob)
            vscale(i15(j,iglob)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nshortvlst(ii)
            kglob = shortvlst(kk,ii)
            kbis = loc(kglob)
            kv = ired(kglob)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(kglob) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(kglob)
               xr = xi - xred(kbis)
               yr = yi - yred(kbis)
               zr = zi - zred(kbis)
               if (use_bounds) call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(kglob) .eq. iglob) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(kglob)
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
c
c     use energy switching if near the cutoff distance
c
                  rik = sqrt(rik2)
                  call switch_respa(rik,off,shortheal,s,ds)
c
c     increment the overall van der Waals energy components
c
                  if (e .ne. 0.0d0) then
                     nev = nev + 1
                     ev = ev + e*s
                     aev(i) = aev(i) + 0.5d0*e*s
                     aev(kbis) = aev(kbis) + 0.5d0*e*s
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             20x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,20) iglob,name(iglob),kglob,name(kglob)
     &                                ,rv,sqrt(rik2),e
   20                format (' VDW-LJ',4x,2(i7,'-',a3),
     &                          13x,2f10.4,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            vscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            vscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            vscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            vscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
c
c     ##########################################################################
c     ##                                                                      ##
c     ##  subroutine eljlon3c  -- long range Lennard-Jones analysis via list  ##
c     ##                                                                      ##
c     ##########################################################################
c
c
c     "eljlong3c" calculates the long range Lennard-Jones van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine eljlong3c
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use couple
      use cutoff
      use domdec
      use energi
      use inter
      use inform
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,iglob,iivdw,kglob,kbis
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e
      real*8 p6,p12,eps
      real*8 rv,rdn
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8 s,ds,vdwshortcut2
      logical proceed,usei
      logical header,huge
      character*10 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      aev = 0d0
      header = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(nbloc))
      allocate (yred(nbloc))
      allocate (zred(nbloc))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      vscale = 1.0d0
      iv14 = 0
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
      vdwshortcut2 = (vdwshortcut-shortheal)**2
c
c     apply any reduction factor to the atomic coordinates
c
      do ii = 1, nvdwbloc
         iivdw = vdwglob(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         rdn = kred(iglob)
         xr = x(iglob) - x(iv)
         yr = y(iglob) - y(iv)
         zr = z(iglob) - z(iv)
         if (use_polymer) call image(xr,yr,zr)
         xred(i) = rdn*xr + x(iv)
         yred(i) = rdn*yr + y(iv)
         zred(i) = rdn*zr + z(iv)
      end do
c
c     find the van der Waals energy via neighbor list search
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         iv = ired(iglob)
         it = jvdw(iglob)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(iglob) .or. use(iv))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            vscale(i12(j,iglob)) = v2scale
         end do
         do j = 1, n13(iglob)
            vscale(i13(j,iglob)) = v3scale
         end do
         do j = 1, n14(iglob)
            vscale(i14(j,iglob)) = v4scale
            iv14(i14(j,iglob)) = iglob
         end do
         do j = 1, n15(iglob)
            vscale(i15(j,iglob)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvlst(ii)
            kglob = vlst(kk,ii)
            kbis = loc(kglob)
            kv = ired(kglob)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(kglob) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(kglob)
               xr = xi - xred(kbis)
               yr = yi - yred(kbis)
               zr = zi - zred(kbis)
               if (use_bounds) call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if ((rik2 .le. off2).and.(rik2.ge.vdwshortcut2)) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(kglob) .eq. iglob) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(kglob)
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
c
c     use energy switching if close the cutoff distance (at short range)
c
                 rik = sqrt(rik2)
                 call switch_respa(rik,vdwshortcut,shortheal,s,ds)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     increment the overall van der Waals energy components
c
                  if (e .ne. 0.0d0) then
                     nev = nev + 1
                     ev = ev + e*(1-s)
                     aev(i) = aev(i) + 0.5d0*e*(1-s)
                     aev(kbis) = aev(kbis) + 0.5d0*e*(1-s)
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             20x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,20) iglob,name(iglob),kglob,name(kglob)
     &                                ,rv,sqrt(rik2),e
   20                format (' VDW-LJ',4x,2(i7,'-',a3),
     &                          13x,2f10.4,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(iglob)
            vscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            vscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            vscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            vscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
c
