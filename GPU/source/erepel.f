c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erepel  --  Pauli repulsion energy & analysis   ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erepel" calculates the Pauli repulsion energy 
c
c     literature reference:
c
c     J. A. Rackers and J. W. Ponder, "Classical Pauli Repulsion:
c     An Anisotropic, Atomic Multipole Model", Journal of Chemical
c     Physics, 150, 084104 (2019)
c
c
      subroutine erepel
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
      call erepel0c
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine erepel0c  --  Pauli repulsion analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "erepel3b" calculates the Pauli repulsion energy 
c
c
      subroutine erepel0c
      use atmtyp
      use atmlst
      use atoms
      use bound
      use couple
      use cutoff
      use domdec
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mpole
      use mutant
      use neigh
      use potent
      use repel
      use reppot
      use shunt
      use usage
      implicit none
      integer i,j,iipole,iglob,kglob,nnelst
      integer ii,kkk,kbis,kkpole
      real*8 e,eterm
      real*8 fgrp,taper
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr,dik,qik
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 diqk,dkqi,qiqk
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 sizi,sizk,sizik
      real*8 vali,valk
      real*8 dmpi,dmpk
      real*8 dmpik(9)
      real*8 s,ds,repshortcut2,facts
      real*8, allocatable :: rscale(:)
      logical proceed,usei
      logical header
      logical muti,mutk,mutik
      logical testcut,shortrange,longrange,fullrange
      character*11 mode
      character*80 :: RoutineName

      shortrange = use_repulsshort
      longrange  = use_repulslong
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'erepel3short3c'
         mode        = 'SHORTREPULS'
      else if (longrange) then
         RoutineName = 'erepel3long3c'
         mode        = 'REPULS'
      else
         RoutineName = 'erepel3c'
         mode        = 'REPULS'
      endif
c
c
c     zero out Pauli repulsion energy and partitioning terms
c
      er = 0.0d0
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (rscale(n))
c
c     initialize connected atom exclusion coefficients
c
      rscale = 1.0d0
c
c     set the coefficients for the switching function
c
      call switch (mode)
      repshortcut2 = (repshortcut-shortheal)**2
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Pauli Repulsion Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         sizi = sizpr(iglob)
         dmpi = dmppr(iglob)
         vali = elepr(iglob)
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
         usei = use(iglob)
         muti = mut(iglob)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            rscale(i12(j,iglob)) = r2scale
         end do
         do j = 1, n13(iglob)
            rscale(i13(j,iglob)) = r3scale
         end do
         do j = 1, n14(iglob)
            rscale(i14(j,iglob)) = r4scale
         end do
         do j = 1, n15(iglob)
            rscale(i15(j,iglob)) = r5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         nnelst = merge(nshortelst(ii),
     &                  nelst     (ii),
     &                  shortrange
     &                 )
         do kkk = 1, nnelst
            kkpole = merge(shortelst(kkk,ii),
     &                     elst     (kkk,ii),
     &                     shortrange
     &                   )
            kglob = ipole(kkpole)
            kbis = loc(kglob)
            mutk = mut(kglob)
            proceed = (use(iglob) .or. use(kglob))
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            if (proceed) then
               xr = x(kglob) - xi
               yr = y(kglob) - yi
               zr = z(kglob) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
c               if (r2 .le. off2) then
               testcut = merge(r2 .le. off2.and.r2.ge.repshortcut2,
     &                         r2 .le. off2,
     &                         longrange
     &                        )
               if (testcut) then
                  r = sqrt(r2)
                  sizk = sizpr(kglob)
                  dmpk = dmppr(kglob)
                  valk = elepr(kglob)
                  ck = rpole(1,kkpole)
                  dkx = rpole(2,kkpole)
                  dky = rpole(3,kkpole)
                  dkz = rpole(4,kkpole)
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
                  dik = dix*dkx + diy*dky + diz*dkz
                  qik = qix*qkx + qiy*qky + qiz*qkz
                  diqk = dix*qkx + diy*qky + diz*qkz
                  dkqi = dkx*qix + dky*qiy + dkz*qiz
                  qiqk = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                      + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     get reciprocal distance terms for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     get damping coefficients for the Pauli repulsion energy
c     
                  call damprep (r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
     &                             9,dmpi,dmpk,dmpik)                  
c
c     compute the Pauli repulsion energy for this interaction
c
                  term1 = vali*valk
                  term2 = valk*dir - vali*dkr + dik
                  term3 = vali*qkr + valk*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  eterm = term1*dmpik(1) + term2*dmpik(3)
     &                       + term3*dmpik(5) + term4*dmpik(7)
     &                       + term5*dmpik(9)
                  sizik = sizi * sizk
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = .false.
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        mutik = .true.
                     else if (.not.muti .or. .not.mutk) then
                        mutik = .true.
                     end if
                  end if
c
c     get interaction energy, via soft core lambda scaling as needed
c
                  if (mutik) then
                     e = vlambda * sizik * rscale(kglob) * eterm
     &                      / sqrt(1.0d0-vlambda+r2)
                  else
                     e = sizik * rscale(kglob) * eterm * rr1
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
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                    end if
                  end if
c
c     scale the interaction based on its group membership
c
                 if (use_group)  e = e * fgrp
                 if(shortrange .or. longrange)
     &            call switch_respa(r,repshortcut,shortheal,s,ds)

                 if(shortrange) then
                    facts =         s
                 else if(longrange) then
                    facts = 1.0d0 - s
                 else
                    facts = 1.0d0
                 endif

                 e  = e * facts
c
c
c     increment the overall Pauli repulsion energy component
c
                  er = er + e
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
            rscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            rscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            rscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            rscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (rscale)
      return
      end
