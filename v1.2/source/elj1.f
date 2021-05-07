c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine elj1  --  Lennard-Jones energy & derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj1" calculates the Lennard-Jones 6-12 van der Waals energy
c     and its first derivatives with respect to Cartesian coordinates
c
c
      subroutine elj1
      use energi
      use potent
      use virial
      use vdwpot
      implicit none
      real*8 elrc,vlrc
c
c     choose the method for summing over pairwise interactions
c
      if (use_vdwshort) then
        call eljshort1c
      else if(use_vdwlong) then
        call eljlong1c
      else
        call elj1c
      end if
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr1 (elrc,vlrc)
         ev = ev + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine elj1c  --  Lennard-Jones vdw derivs via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj1c" calculates the Lennard-Jones 12-6 van der Waals energy
c     and its first derivatives using a pairwise neighbor list
c
c
      subroutine elj1c
      use atmlst
      use atoms
      use bound
      use couple
      use deriv
      use domdec
      use energi
      use inter
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      implicit none
      integer i,j,iglob,kglob,kbis,iivdw
      integer ii,iv,it,ivloc,kvloc
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,de,p6,p12,eps
      real*8 rv,rdn
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 taper,dtaper
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy


      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical usei
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate VDW')
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      dev = 0.0d0
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
c     find van der Waals energy and derivatives via neighbor list
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         if (i.eq.0) then
           write(iout,1000)
           cycle
         end if
         iv = ired(iglob)
         ivloc = loc(iv)
         redi = kred(iglob)
         rediv = 1.0d0 - redi
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
            kvloc = loc(kv)
            if ((kvloc.eq.0).or.(kbis.gt.nbloc).or.(kvloc.gt.nbloc)) 
     $          then
              write(iout,1000)
              cycle
            end if
c
c     compute the energy contribution for this interaction
c
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
               rik = sqrt(rik2)
               p6 = rv**6 / rik2**3
               p12 = p6 * p6
               e = eps * (p12-2.0d0*p6)
               de = eps * (p12-p6) * (-12.0d0/rik)
c
c     use energy switching if near the cutoff distance
c
                if (rik2 .gt. cut2) then
                   rik3 = rik2 * rik
                   rik4 = rik2 * rik2
                   rik5 = rik2 * rik3
                   taper = c5*rik5 + c4*rik4 + c3*rik3
     &                        + c2*rik2 + c1*rik + c0
                   dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                         + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                   de = e*dtaper + de*taper
                   e = e * taper
                end if
c
c     find the chain rule terms for derivative components
c
                de = de / rik
                dedx = de * xr
                dedy = de * yr
                dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                ev  = ev  + e
                if (iglob .eq. iv) then
                   dev(1,i) = dev(1,i) + dedx
                   dev(2,i) = dev(2,i) + dedy
                   dev(3,i) = dev(3,i) + dedz
                else
                   dev(1,i) = dev(1,i) + dedx*redi
                   dev(2,i) = dev(2,i) + dedy*redi
                   dev(3,i) = dev(3,i) + dedz*redi
                   dev(1,ivloc) = dev(1,ivloc) + dedx*rediv
                   dev(2,ivloc) = dev(2,ivloc) + dedy*rediv
                   dev(3,ivloc) = dev(3,ivloc) + dedz*rediv
                end if
                if (kglob .eq. kv) then
                   dev(1,kbis) = dev(1,kbis) - dedx
                   dev(2,kbis) = dev(2,kbis) - dedy
                   dev(3,kbis) = dev(3,kbis) - dedz
                else
                   redk = kred(kglob)
                   redkv = 1.0d0 - redk
                   dev(1,kbis) = dev(1,kbis) - dedx*redk
                   dev(2,kbis) = dev(2,kbis) - dedy*redk
                   dev(3,kbis) = dev(3,kbis) - dedz*redk
                   dev(1,kvloc) = dev(1,kvloc) - dedx*redkv
                   dev(2,kvloc) = dev(2,kvloc) - dedy*redkv
                   dev(3,kvloc) = dev(3,kvloc) - dedz*redkv
                end if
c
c     increment the internal virial tensor components
c
                vxx = xr * dedx
                vyx = yr * dedx
                vzx = zr * dedx
                vyy = yr * dedy
                vzy = zr * dedy
                vzz = zr * dedz
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
c     increment the total intermolecular energy
c
                if (molcule(iglob) .ne. molcule(kglob)) then
                   einter = einter + e
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
c     ################################################################################
c     ##                                                                            ##
c     ##  subroutine eljshort1c  --  short range Lennard-Jones vdw derivs via list  ##
c     ##                                                                            ##
c     ################################################################################
c
c
c     "eljshort1c" calculates the short range Lennard-Jones 12-6 van der Waals energy
c     and its first derivatives using a pairwise neighbor list
c
c
      subroutine eljshort1c
      use atmlst
      use atoms
      use bound
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use inter
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      implicit none
      integer i,j,iglob,kglob,kbis,iivdw
      integer ii,iv,it,ivloc,kvloc
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,de,p6,p12,eps
      real*8 rv,rdn
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2


      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy


      real*8 s,ds
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical usei
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate VDW')
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      dev = 0.0d0
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
c     find van der Waals energy and derivatives via neighbor list
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         if (i.eq.0) then
           write(iout,1000)
           cycle
         end if
         iv = ired(iglob)
         ivloc = loc(iv)
         redi = kred(iglob)
         rediv = 1.0d0 - redi
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
            kvloc = loc(kv)
            if ((kvloc.eq.0).or.(kbis.gt.nbloc).or.(kvloc.gt.nbloc)) 
     $          then
              write(iout,1000)
              cycle
            end if
c
c     compute the energy contribution for this interaction
c
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
               rik = sqrt(rik2)
               p6 = rv**6 / rik2**3
               p12 = p6 * p6
               e = eps * (p12-2.0d0*p6)
               de = eps * (p12-p6) * (-12.0d0/rik)
c
c     use energy switching if near the cutoff distance
c
               call switch_respa(rik,off,shortheal,s,ds)
               e = e * s
               de = e*ds + de*s
c
c     find the chain rule terms for derivative components
c
                de = de / rik
                dedx = de * xr
                dedy = de * yr
                dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                ev  = ev  + e
                if (iglob .eq. iv) then
                   dev(1,i) = dev(1,i) + dedx
                   dev(2,i) = dev(2,i) + dedy
                   dev(3,i) = dev(3,i) + dedz
                else
                   dev(1,i) = dev(1,i) + dedx*redi
                   dev(2,i) = dev(2,i) + dedy*redi
                   dev(3,i) = dev(3,i) + dedz*redi
                   dev(1,ivloc) = dev(1,ivloc) + dedx*rediv
                   dev(2,ivloc) = dev(2,ivloc) + dedy*rediv
                   dev(3,ivloc) = dev(3,ivloc) + dedz*rediv
                end if
                if (kglob .eq. kv) then
                   dev(1,kbis) = dev(1,kbis) - dedx
                   dev(2,kbis) = dev(2,kbis) - dedy
                   dev(3,kbis) = dev(3,kbis) - dedz
                else
                   redk = kred(kglob)
                   redkv = 1.0d0 - redk
                   dev(1,kbis) = dev(1,kbis) - dedx*redk
                   dev(2,kbis) = dev(2,kbis) - dedy*redk
                   dev(3,kbis) = dev(3,kbis) - dedz*redk
                   dev(1,kvloc) = dev(1,kvloc) - dedx*redkv
                   dev(2,kvloc) = dev(2,kvloc) - dedy*redkv
                   dev(3,kvloc) = dev(3,kvloc) - dedz*redkv
                end if
c
c     increment the internal virial tensor components
c
                vxx = xr * dedx
                vyx = yr * dedx
                vzx = zr * dedx
                vyy = yr * dedy
                vzy = zr * dedy
                vzz = zr * dedz
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
c     increment the total intermolecular energy
c
                if (molcule(iglob) .ne. molcule(kglob)) then
                   einter = einter + e
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
c     #############################################################################
c     ##                                                                         ##
c     ##  subroutine eljlong1c  -- long range Lennard-Jones vdw derivs via list  ##
c     ##                                                                         ##
c     #############################################################################
c
c
c     "eljlong1c" calculates the long range Lennard-Jones 12-6 van der Waals energy
c     and its first derivatives using a pairwise neighbor list
c
c
      subroutine eljlong1c
      use atmlst
      use atoms
      use bound
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use inter
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      implicit none
      integer i,j,iglob,kglob,kbis,iivdw
      integer ii,iv,it,ivloc,kvloc
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,de,p6,p12,eps
      real*8 rv,rdn
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 taper,dtaper
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy


      real*8 s,ds,vdwshortcut2
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical usei
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate VDW')
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      dev = 0.0d0
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
c     find van der Waals energy and derivatives via neighbor list
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         if (i.eq.0) then
           write(iout,1000)
           cycle
         end if
         iv = ired(iglob)
         ivloc = loc(iv)
         redi = kred(iglob)
         rediv = 1.0d0 - redi
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
            kvloc = loc(kv)
            if ((kvloc.eq.0).or.(kbis.gt.nbloc).or.(kvloc.gt.nbloc)) 
     $          then
              write(iout,1000)
              cycle
            end if
c
c     compute the energy contribution for this interaction
c
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
               rik = sqrt(rik2)
               p6 = rv**6 / rik2**3
               p12 = p6 * p6
               e = eps * (p12-2.0d0*p6)
               de = eps * (p12-p6) * (-12.0d0/rik)
c
c     use energy switching if close the cutoff distance (at short range)
c
               call switch_respa(rik,vdwshortcut,shortheal,s,ds)
               e = (1-s)*e
               de = -e*ds-s*de
c
c     use energy switching if near the cutoff distance
c
                if (rik2 .gt. cut2) then
                   rik3 = rik2 * rik
                   rik4 = rik2 * rik2
                   rik5 = rik2 * rik3
                   taper = c5*rik5 + c4*rik4 + c3*rik3
     &                        + c2*rik2 + c1*rik + c0
                   dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                         + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                   de = e*dtaper + de*taper
                   e = e * taper
                end if
c
c     find the chain rule terms for derivative components
c
                de = de / rik
                dedx = de * xr
                dedy = de * yr
                dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                ev  = ev  + e
                if (iglob .eq. iv) then
                   dev(1,i) = dev(1,i) + dedx
                   dev(2,i) = dev(2,i) + dedy
                   dev(3,i) = dev(3,i) + dedz
                else
                   dev(1,i) = dev(1,i) + dedx*redi
                   dev(2,i) = dev(2,i) + dedy*redi
                   dev(3,i) = dev(3,i) + dedz*redi
                   dev(1,ivloc) = dev(1,ivloc) + dedx*rediv
                   dev(2,ivloc) = dev(2,ivloc) + dedy*rediv
                   dev(3,ivloc) = dev(3,ivloc) + dedz*rediv
                end if
                if (kglob .eq. kv) then
                   dev(1,kbis) = dev(1,kbis) - dedx
                   dev(2,kbis) = dev(2,kbis) - dedy
                   dev(3,kbis) = dev(3,kbis) - dedz
                else
                   redk = kred(kglob)
                   redkv = 1.0d0 - redk
                   dev(1,kbis) = dev(1,kbis) - dedx*redk
                   dev(2,kbis) = dev(2,kbis) - dedy*redk
                   dev(3,kbis) = dev(3,kbis) - dedz*redk
                   dev(1,kvloc) = dev(1,kvloc) - dedx*redkv
                   dev(2,kvloc) = dev(2,kvloc) - dedy*redkv
                   dev(3,kvloc) = dev(3,kvloc) - dedz*redkv
                end if
c
c     increment the internal virial tensor components
c
                vxx = xr * dedx
                vyx = yr * dedx
                vzx = zr * dedx
                vyy = yr * dedy
                vzy = zr * dedy
                vzz = zr * dedz
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
c     increment the total intermolecular energy
c
                if (molcule(iglob) .ne. molcule(kglob)) then
                   einter = einter + e
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
