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
#include "tinker_precision.h"
      subroutine empole0
      implicit none
c
c     choose the method for summing over multipole interactions
c
      call empole0c
c
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
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer i,ii,iglob,iipole,ierr
      real(t_p) e,f
      real(t_p) term,fterm
      real(t_p) cii,dii,qii
      real(t_p) xd,yd,zd
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
c
c
c     zero out the total atomic multipole energy
c
      em = 0.0_re_p
      if (npole .eq. 0)  return
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
          if (use_mpoleshortreal) then
            call emrealshort0d
          else if (use_mpolelong) then
            call emreallong0d
          else
            call emreal0d
          end if
        end if

        if (use_mself) then
c
c     compute the self-energy portion of the Ewald summation
c
        term = 2.0_ti_p * aewald * aewald
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
           qii = 2.0_ti_p*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &              + qixx*qixx + qiyy*qiyy + qizz*qizz
           e = fterm * (cii + term*(dii/3.0_ti_p+
     &                              2.0_ti_p*term*qii/5.0_ti_p))
           em = em + e
        end do
c
c       compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
           xd = 0.0_ti_p
           yd = 0.0_ti_p
           zd = 0.0_ti_p
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
             call MPI_ALLREDUCE(MPI_IN_PLACE,xd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,yd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,zd,1,MPI_TPREC,MPI_SUM,
     $          COMM_TINKER,ierr)
             if (rank.eq.0) then
           term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
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
      subroutine emreal0d
      use sizes
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use mplpot
      use neigh
      use shunt
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer i,j,k,inl,iglob,kglob,kbis
      integer ii,kk,kkk,iipole,kkpole
      real(t_p) e,f,bfac
      real(t_p) alsq2,alsq2n
      real(t_p) exp2a,ralpha
      real(t_p) scalekk
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1,rr3
      real(t_p) rr5,rr7,rr9
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) dri,drk,dik
      real(t_p) qrri,qrrk
      real(t_p) qrrik,qik
      real(t_p) diqrk,dkqri
      real(t_p) term1,term2,term3
      real(t_p) term4,term5
      real(t_p) bn(0:4)
      real(t_p), allocatable :: mscale(:)
      character*10 mode
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      mscale = 1.0_ti_p
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
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
         do kkk = 1, nelst(ii)
            kkpole = elst(kkk,ii)
            kglob = ipole(kkpole)
            kbis = loc(kglob)
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
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
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0_ti_p * rr3 / r2
               rr7 = 5.0_ti_p * rr5 / r2
               rr9 = 7.0_ti_p * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0_ti_p * aewald**2
               alsq2n = 0.0_ti_p
               if (aewald .gt. 0.0_ti_p)  
     &            alsq2n = 1.0_ti_p / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = real(j+j-1,t_p)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0_ti_p*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0_ti_p*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0_ti_p*qrrik
               term5 = qrri*qrrk
c
c     modify distances to account for Ewald and exclusions
c
               scalekk = 1.0_ti_p - mscale(kglob)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0_ti_p
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c     ##################################################################################
c     ##                                                                              ##
c     ##  subroutine emrealshort0d  --  short range real space mpole energy via list  ##
c     ##                                                                              ##
c     ##################################################################################
c
c
c     "emrealshort0d" evaluates the short range real space portion of the Ewald sum
c     energy due to atomic multipoles using a neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emrealshort0d
      use sizes
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use cutoff
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use mplpot
      use neigh
      use shunt
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer i,j,iglob,kglob,kbis
      integer ii,kkk,iipole,kkpole
      real(t_p) e,f,bfac
      real(t_p) alsq2,alsq2n
      real(t_p) exp2a,ralpha
      real(t_p) scalekk
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1,rr3
      real(t_p) rr5,rr7,rr9
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) dri,drk,dik
      real(t_p) qrri,qrrk
      real(t_p) qrrik,qik
      real(t_p) diqrk,dkqri
      real(t_p) term1,term2,term3
      real(t_p) term4,term5
      real(t_p) bn(0:4)
      real(t_p) s,ds
      real(t_p), allocatable :: mscale(:)
      character*10 mode
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      mscale = 1.0_ti_p
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'SHORTEWALD'
      call switch (mode)
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
         do kkk = 1, nelst(ii)
            kkpole = elst(kkk,ii)
            kglob = ipole(kkpole)
            kbis = loc(kglob)
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
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
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0_ti_p * rr3 / r2
               rr7 = 5.0_ti_p * rr5 / r2
               rr9 = 7.0_ti_p * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0_ti_p * aewald**2
               alsq2n = 0.0_ti_p
               if (aewald .gt. 0.0_ti_p)
     &            alsq2n = 1.0_ti_p / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0_ti_p*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0_ti_p*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0_ti_p*qrrik
               term5 = qrri*qrrk
c
c     modify distances to account for Ewald and exclusions
c
               scalekk = 1.0_ti_p - mscale(kglob)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               call switch_respa(r,off,shortheal,s,ds)
               em = em + e*s
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0_ti_p
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c     ##################################################################################
c     ##                                                                              ##
c     ##  subroutine emlongshort0d  --  long range real space mpole energy via list   ##
c     ##                                                                              ##
c     ##################################################################################
c
c
c     "emrealshort0d" evaluates the short range real space portion of the Ewald sum
c     energy due to atomic multipoles using a neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emreallong0d
      use sizes
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use cutoff
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use mplpot
      use neigh
      use shunt
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
      integer i,j,iglob,kglob,kbis
      integer ii,kkk,iipole,kkpole
      real(t_p) e,f,bfac
      real(t_p) alsq2,alsq2n
      real(t_p) exp2a,ralpha
      real(t_p) scalekk
      real(t_p) xi,yi,zi
      real(t_p) xr,yr,zr
      real(t_p) r,r2,rr1,rr3
      real(t_p) rr5,rr7,rr9
      real(t_p) ci,dix,diy,diz
      real(t_p) qixx,qixy,qixz
      real(t_p) qiyy,qiyz,qizz
      real(t_p) ck,dkx,dky,dkz
      real(t_p) qkxx,qkxy,qkxz
      real(t_p) qkyy,qkyz,qkzz
      real(t_p) qrix,qriy,qriz
      real(t_p) qrkx,qrky,qrkz
      real(t_p) dri,drk,dik
      real(t_p) qrri,qrrk
      real(t_p) qrrik,qik
      real(t_p) diqrk,dkqri
      real(t_p) term1,term2,term3
      real(t_p) term4,term5
      real(t_p) bn(0:4)
      real(t_p) s,ds,mpoleshortcut2
      real(t_p), allocatable :: mscale(:)
      character*10 mode
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      mscale = 1.0_ti_p
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
      mpoleshortcut2 = (mpoleshortcut-shortheal)**2
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
         do kkk = 1, nelst(ii)
            kkpole = elst(kkk,ii)
            kglob = ipole(kkpole)
            kbis = loc(kglob)
            xr = x(kglob) - xi
            yr = y(kglob) - yi
            zr = z(kglob) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if ((r2 .le. off2).and.(r2.ge.mpoleshortcut2)) then
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
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0_ti_p * rr3 / r2
               rr7 = 5.0_ti_p * rr5 / r2
               rr9 = 7.0_ti_p * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0_ti_p * aewald**2
               alsq2n = 0.0_ti_p
               if (aewald .gt. 0.0_ti_p)
     &            alsq2n = 1.0_ti_p / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0_ti_p*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0_ti_p*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0_ti_p*qrrik
               term5 = qrri*qrrk
c
c     modify distances to account for Ewald and exclusions
c
               scalekk = 1.0_ti_p - mscale(kglob)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               call switch_respa(r,off,shortheal,s,ds)
               em = em + (1-s)*e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0_ti_p
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0_ti_p
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
