c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal3  --  buffered 14-7 vdw energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehal3" calculates the buffered 14-7 van der Waals energy
c     and partitions the energy among the atoms
c
c
      subroutine ehal3
      use analyz
      use atoms
      use domdec
      use energi
      use inform
      use iounit
      use potent
      use vdwpot
      use mpi
      implicit none
      integer i
      real*8 elrc,aelrc
      character*11 mode
c
c     evaluate pairwise interactions
c
      call ehal3c
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         mode = 'VDW'
         call evcorr (mode,elrc)
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal3c  --  buffered 14-7 analysis via list    ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal3c" calculates the buffered 14-7 van der Waals energy
c      with respect to Cartesian coordinates using a paiwise neighbor list

c     if longrange , calculates just the long range part
c     if shortrange, calculates just the short range part
c
      subroutine ehal3c
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
      use vdw
      use vdwpot
      use mpi
      implicit none
      integer i,j,iglob,kglob,kbis,iivdw,nnvlst
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,eps,rdn,fgrp
      real*8 rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 scal,t1,t2
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8 s,ds,vdwshortcut2,facts
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical muti,mutk,mutik
      logical header,huge
      logical testcut,shortrange,longrange,fullrange
      character*11 mode
      character*80 :: RoutineName
c

c     choose the method for summing over pairwise interactions
      shortrange = use_vdwshort
      longrange  = use_vdwlong
      fullrange  = .not.(shortrange.or.longrange)

      if (shortrange) then 
         RoutineName = 'ehalshort3c'
         mode        = 'SHORTVDW'
      else if (longrange) then
         RoutineName = 'ehallong3c'
         mode        = 'VDW'
      else
         RoutineName = 'ehal3c'
         mode        = 'VDW'
      endif

c
c     zero out the van der Waals energy and partitioning terms
c
      nev          = 0
      ev = 0.0d0
      aev = 0.0d0
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
c     find van der Waals energy via neighbor list search
c
      do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         iv    = ired(iglob)
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
         muti  = mut(iglob)
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

         if (shortrange) then
           nnvlst = nshortvlst(ii)
         else
           nnvlst = nvlst(ii)
         end if
         do kk = 1, nnvlst
            if (shortrange) then
              kglob = shortvlst(kk,ii)
            else
              kglob = vlst(kk,ii)
            end if
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            kbis = loc(kglob)
            kv = ired(kglob)
            mutk = mut(kglob)
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
            testcut = merge(rik2 .le. off2.and.rik2.ge.vdwshortcut2,
     &                      rik2 .le. off2,
     &                      longrange
     &                     )
            if (testcut) then
               rik = sqrt(rik2)
               rv  = radmin(kt,it)
               eps = epsilon(kt,it)
               if (iv14(kglob) .eq. iglob) then
                  rv = radmin4(kt,it)
                  eps = epsilon4(kt,it)
               end if
               eps = eps * vscale(kglob)

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
c     get the interaction energy, via soft core if necessary
c
               if (mutik) then
                  rho = rik / rv
                  eps = eps * vlambda**scexp
                  scal = scalpha * (1.0d0-vlambda)**2
                  t1   = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                  t2   = (1.0d0+ghal) / (scal+rho**7+ghal)
                  e    = eps * t1 * (t2-2.0d0)
               else
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + ghal*rv7
                  tau = (dhal+1.0d0) / (rik + dhal*rv)
                  e = eps * rv7 * tau**7
     &                    * ((ghal+1.0d0)*rv7/rho-2.0d0)
               end if
c
c     use energy switching if near the cutoff distance at short range
c
               if(longrange.or.fullrange) then
                  if (rik2 .gt. cut2) then
                     rik3  = rik2 * rik
                     rik4  = rik2 * rik2
                     rik5  = rik2 * rik3
                     taper =  c5*rik5 + c4*rik4 + c3*rik3
     &                      + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
               endif
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp

               if(shortrange .or. longrange)
     &            call switch_respa(rik,vdwshortcut,shortheal,s,ds)

c
c     increment the overall van der Waals energy components
c
               if(shortrange) then
                  facts =         s
               else if(longrange) then
                  facts = 1.0d0 - s
               else
                  facts = 1.0d0
               endif

               e  = e * facts

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
   10                 format (/,' Individual van der Waals',
     &                           ' Interactions :',
     &                        //,' Type',14x,'Atom Names',
     &                           20x,'Minimum',4x,'Actual',
     &                           6x,'Energy',/)
                   end if
                   write (iout,20)  iglob,name(iglob),kglob,
     $                   name(kglob),
     &                              rv,sqrt(rik2),e
   20              format (' VDW-Hal',3x,2(i7,'-',a3),
     &                        13x,2f10.4,f12.4)
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
