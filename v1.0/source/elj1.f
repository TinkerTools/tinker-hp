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
      implicit none
      real*8 elrc,vlrc
      include 'sizes.i'
      include 'energi.i'
      include 'vdwpot.i'
      include 'virial.i'
c
c
c     choose the method for summing over pairwise interactions
c
      call elj1c
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
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmlst.i'
      include 'bound.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'neigh.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'virial.i'
      include 'openmp.i'
      integer i,j,k,iglob,kglob,kbis,iivdw,inl
      integer ii,iv,it,ivloc,kvloc
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,de,p6,p12,eps
      real*8 rv,rdn,fgrp
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
      real*8 evt,eintert
      real*8 virt(3,3)
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: devt(:,:)
      logical proceed,usei
      character*6 mode
c
 1000 format(' Warning, system moved too much since last neighbor list
     $   update, try lowering nlupdate VDW')
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
      allocate (devt(3,nbloc))
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
         xred(i) = rdn*(x(iglob)-x(iv)) + x(iv)
         yred(i) = rdn*(y(iglob)-y(iv)) + y(iv)
         zred(i) = rdn*(z(iglob)-z(iv)) + z(iv)
      end do
c
c     transfer global to local copies for OpenMP calculation
c
      evt = ev
      eintert = einter
      devt = dev
      virt = vir
c
c     set OpenMP directives for the major loop structure
c
c!$OMP PARALLEL default(private) shared(nvdw,ivdw,ired,kred,
c!$OMP& jvdw,xred,yred,zred,use,nvlst,vlst,n12,n13,n14,n15,
c!$OMP& i12,i13,i14,i15,v2scale,v3scale,v4scale,v5scale,
c!$OMP& use_group,off2,radmin,epsilon,radmin4,epsilon4,
c!$OMP& cut2,c0,c1,c2,c3,c4,c5,molcule) firstprivate(vscale,iv14)
c!$OMP& shared(evt,devt,virt,eintert)
c!$OMP DO reduction(+:evt,devt,virt,eintert) schedule(dynamic)
c
c     find van der Waals energy and derivatives via neighbor list
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         inl = vdwlocnl(iivdw)
         if (inl.eq.0) then
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
            kv = ired(k)
            kvloc = loc(kv)
            if ((kvloc.eq.0).or.(kbis.gt.nbloc).or.(kvloc.gt.nbloc)) 
     $          then
              write(iout,1000)
              cycle
            end if
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
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
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
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
                  evt = evt + e
                  if (iglob .eq. iv) then
                     devt(1,i) = devt(1,i) + dedx
                     devt(2,i) = devt(2,i) + dedy
                     devt(3,i) = devt(3,i) + dedz
                  else
                     devt(1,i) = devt(1,i) + dedx*redi
                     devt(2,i) = devt(2,i) + dedy*redi
                     devt(3,i) = devt(3,i) + dedz*redi
                     devt(1,ivloc) = devt(1,ivloc) + dedx*rediv
                     devt(2,ivloc) = devt(2,ivloc) + dedy*rediv
                     devt(3,ivloc) = devt(3,ivloc) + dedz*rediv
                  end if
                  if (kglob .eq. kv) then
                     devt(1,kbis) = devt(1,kbis) - dedx
                     devt(2,kbis) = devt(2,kbis) - dedy
                     devt(3,kbis) = devt(3,kbis) - dedz
                  else
                     redk = kred(kglob)
                     redkv = 1.0d0 - redk
                     devt(1,kbis) = devt(1,kbis) - dedx*redk
                     devt(2,kbis) = devt(2,kbis) - dedy*redk
                     devt(3,kbis) = devt(3,kbis) - dedz*redk
                     devt(1,kvloc) = devt(1,kvloc) - dedx*redkv
                     devt(2,kvloc) = devt(2,kvloc) - dedy*redkv
                     devt(3,kvloc) = devt(3,kvloc) - dedz*redkv
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
                  virt(1,1) = virt(1,1) + vxx
                  virt(2,1) = virt(2,1) + vyx
                  virt(3,1) = virt(3,1) + vzx
                  virt(1,2) = virt(1,2) + vyx
                  virt(2,2) = virt(2,2) + vyy
                  virt(3,2) = virt(3,2) + vzy
                  virt(1,3) = virt(1,3) + vzx
                  virt(2,3) = virt(2,3) + vzy
                  virt(3,3) = virt(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     eintert = eintert + e
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
c     end OpenMP directives for the major loop structure
c
c!$OMP END DO
c!$OMP END PARALLEL
c
c     transfer local to global copies for OpenMP calculation
c
      ev = evt
      einter = eintert
      dev = devt
      vir = virt
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      deallocate (devt)
      return
      end
