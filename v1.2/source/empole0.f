c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole3  --  mpole/polar energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole3" calculates the electrostatic energy due to
c     atomic multipole and dipole polarizability interactions
c
c
      subroutine empole0
      use potent
      implicit none
c
c     choose the method for summing over multipole interactions
c
      if (use_lambdadyn) then
        call elambdampole0c
      else
        call empole0c
      end if
c
      return
      end
c
c     ################################################################
c     ##                                                                  ##
c     ##  subroutine elambdampole0c  --  Ewald multipole energy via list  ##
c     ##                                                                  ##
c     ################################################################
c
c
c     "elambdampole0d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list, if lambdadyn
c     is activated
c
c
      subroutine elambdampole0c
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use mutant
      use potent
      use mpi
      implicit none
      integer i,ii,iglob,iipole,ierr
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 :: elambdatemp
      real*8 :: elambdarec0,elambdarec1
c
      elambdatemp = elambda  
c
c     zero out the total atomic multipole energy
c
      em = 0.0d0
      if (npole .eq. 0)  return
      aewald = aeewald
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        if (use_mrec) then
c
c         the reciprocal part is interpolated between 0 and 1
c
          elambda = 0d0
          call MPI_BARRIER(hostcomm,ierr)
          if (hostrank.eq.0) call altelec
          call MPI_BARRIER(hostcomm,ierr)
          call rotpole
          em = 0d0
          if (elambda.lt.1d0) then
            call emrecip
          end if
          elambdarec0  = em

          elambda = 1d0
          call MPI_BARRIER(hostcomm,ierr)
          if (hostrank.eq.0) call altelec
          call MPI_BARRIER(hostcomm,ierr)
          call rotpole
          em = 0d0
          if (elambda.gt.0d0) then
            call emrecip
          end if
          elambdarec1  = em

          elambda = elambdatemp 
          em = (1-elambda)*elambdarec0 + elambda*elambdarec1
c
c         reset lambda to initial value
c
          call MPI_BARRIER(hostcomm,ierr)
          if (hostrank.eq.0) call altelec
          call MPI_BARRIER(hostcomm,ierr)
          call rotpole
        end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_mreal) then
          call emreal0d
        end if

        if (use_mself) then
c
c     compute the self-energy portion of the Ewald summation
c
          term = 2.0d0 * aewald * aewald
          fterm = -f * aewald / sqrtpi
          do ii = 1, npoleloc
             iipole = poleglob(ii)
             iglob = ipole(iipole)
             i = loc(iglob)
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
             cii = ci*ci
             dii = dix*dix + diy*diy + diz*diz
             qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &                + qixx*qixx + qiyy*qiyy + qizz*qizz
             e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
             em = em + e
          end do
c
c       compute the cell dipole boundary correction term
c
          if (boundary .eq. 'VACUUM') then
             xd = 0.0d0
             yd = 0.0d0
             zd = 0.0d0
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                dix = rpole(2,iipole)
                diy = rpole(3,iipole)
                diz = rpole(4,iipole)
                xd = xd + dix + rpole(1,iipole)*x(iglob)
                yd = yd + diy + rpole(1,iipole)*y(iglob)
                zd = zd + diz + rpole(1,iipole)*z(iglob)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               term = (2.0d0/3.0d0) * f * (pi/volbox)
               e = term * (xd*xd+yd*yd+zd*zd)
               em = em + e
             end if
          end if
        end if
      end if
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole0c  --  Ewald multipole energy via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole0d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list
c
c
      subroutine empole0c
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use potent
      use mpi
      implicit none
      integer i,ii,iglob,iipole,ierr
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the total atomic multipole energy
c
      em = 0.0d0
      if (npole .eq. 0)  return
      aewald = aeewald
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        if (use_mrec) then
          call emrecip
        end if
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        if (use_mreal) then
          call emreal0d
        end if

        if (use_mself) then
c
c     compute the self-energy portion of the Ewald summation
c
          term = 2.0d0 * aewald * aewald
          fterm = -f * aewald / sqrtpi
          do ii = 1, npoleloc
             iipole = poleglob(ii)
             iglob = ipole(iipole)
             i = loc(iglob)
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
             cii = ci*ci
             dii = dix*dix + diy*diy + diz*diz
             qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &                + qixx*qixx + qiyy*qiyy + qizz*qizz
             e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
             em = em + e
          end do
c
c       compute the cell dipole boundary correction term
c
          if (boundary .eq. 'VACUUM') then
             xd = 0.0d0
             yd = 0.0d0
             zd = 0.0d0
             do ii = 1, npoleloc
                iipole = poleglob(ii)
                iglob = ipole(iipole)
                dix = rpole(2,iipole)
                diy = rpole(3,iipole)
                diz = rpole(4,iipole)
                xd = xd + dix + rpole(1,iipole)*x(iglob)
                yd = yd + diy + rpole(1,iipole)*y(iglob)
                zd = zd + diz + rpole(1,iipole)*z(iglob)
             end do
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_REAL8,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
               term = (2.0d0/3.0d0) * f * (pi/volbox)
               e = term * (xd*xd+yd*yd+zd*zd)
               em = em + e
             end if
          end if
        end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal0d  --  real space mpole energy via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal0d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipoles using a neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
c     if longrange, calculates just the long range part
c     if shortrange, calculates just the short range part

      subroutine emreal0d
      use sizes
      use atmlst
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use cutoff
      use domdec
      use energi
      use ewald
      use group
      use math
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      use mpi
      implicit none
      integer i,j,iglob,kglob,nnelst
      integer ii,kkk,iipole,kkpole
      real*8 e,f
      real*8 scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 rr1i,rr3i,rr5i
      real*8 rr1k,rr3k,rr5k
      real*8 rr1ik,rr3ik,rr5ik
      real*8 rr7ik,rr9ik
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
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term1i,term2i,term3i
      real*8 term1k,term2k,term3k
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9),dmpe(9)
      real*8 fgrp
      real*8 s,ds,mpoleshortcut2
      real*8 facts
      logical testcut,shortrange,longrange,fullrange
      real*8, allocatable :: mscale(:)
      character*11 mode
      character*80 :: RoutineName


c     compute the short, long, or full real space part of the summation
      shortrange = use_mpoleshortreal
      longrange  = use_mpolelong
      fullrange  = .not.(shortrange.or.longrange)
      if (shortrange) then 
         RoutineName = 'emrealshort3d'
         mode        = 'SHORTEWALD'
      else if (longrange) then
         RoutineName = 'emreallong3d'
         mode        = 'EWALD'
      else
         RoutineName = 'emreal3d'
         mode        = 'EWALD'
      endif

c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      mscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc  (iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
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
         if (use_chgpen) then
            corei = pcore(iipole)
            vali = pval(iipole)
            alphai = palpha(iipole)
         end if
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = m2scale
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = m3scale
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = m4scale
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         if (shortrange) then
           nnelst = nshortelst(ii)
         else
           nnelst = nelst(ii)
         end if
         do kkk = 1, nnelst
            if (shortrange) then
              kkpole = shortelst(kkk,ii)
            else
              kkpole = elst(kkk,ii)
            end if
            kglob = ipole(kkpole)
            if (use_group)  call groups (fgrp,iglob,kglob,0,0,0,0)
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            testcut = merge(r2 .le. off2.and.r2.ge.mpoleshortcut2,
     &                      r2 .le. off2,
     &                      longrange
     &                     )
            if (testcut) then
               r = sqrt(r2)
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
     &                   + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
               call dampewald (9,r,r2,f,dmpe)
c
c     find damped multipole intermediates and energy value
c
               if (use_chgpen) then
                  corek = pcore(kkpole)
                  valk = pval(kkpole)
                  alphak = palpha(kkpole)
                  term1 = corei*corek
                  term1i = corek*vali
                  term2i = corek*dir
                  term3i = corek*qir
                  term1k = corei*valk
                  term2k = -corei*dkr
                  term3k = corei*qkr
                  term1ik = vali*valk
                  term2ik = valk*dir - vali*dkr + dik
                  term3ik = vali*qkr + valk*qir - dir*dkr
     &                         + 2.0d0*(dkqi-diqk+qiqk)
                  term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                  term5ik = qir*qkr
                  call damppole (r,9,alphai,alphak,
     &                            dmpi,dmpk,dmpik)
                  scalek = mscale(kglob)
                  if (use_group)  scalek = scalek * fgrp
                  rr1i = dmpe(1) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
                  rr1k = dmpe(1) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
                  rr1ik = dmpe(1) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = dmpe(9) - (1.0d0-scalek*dmpik(9))*rr9
                  rr1 = dmpe(1) - (1.0d0-scalek)*rr1
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  scalek = 1.0d0 - mscale(kglob)
                  if (use_group)  then
                    scalek = 1.0d0 - mscale(kglob)*fgrp
                  end if
                  rr1 = dmpe(1) - scalek*rr1
                  rr3 = dmpe(3) - scalek*rr3
                  rr5 = dmpe(5) - scalek*rr5
                  rr7 = dmpe(7) - scalek*rr7
                  rr9 = dmpe(9) - scalek*rr9
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
               end if

               if(shortrange .or. longrange)
     &            call switch_respa(r,ewaldshortcut,shortheal,s,ds)
c
c     fix the s factor, depending on the range
c
               if(shortrange) then
                  facts  =         s
               else if(longrange) then
                  facts  = 1.0d0 - s
               else
                  facts  = 1.0d0
               endif

               em = em + facts * e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
