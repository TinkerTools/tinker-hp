c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3  --  charge-charge energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms
c
c
#include "tinker_macro.h"
      subroutine echarge3vec
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      call echarge3cvec
c
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine echarge3cvec  --  Ewald charge analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c
c     "echarge3cvec" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a particle
c     mesh Ewald summation
c
c
      subroutine echarge3cvec
      use action
      use analyz
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use domdec
      use energi
      use ewald
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use usage
      use mpi
      use vec
      use vec_elec
      use vec_charge
      use utilvec

      implicit none
      integer i,k,ii,kk
      integer iglob,iichg,inl
      integer nnchg,nnchg1,nnchg2,ncloop8,ncloop16
      integer iglobdefault,kkchgdefault
      real(t_p) f,fi
      real(t_p) fs
      real(t_p) xi,yi,zi,e
      real(t_p) xd,yd,zd
      real(t_p) term
      real(t_p) one,half
!DIR$ ATTRIBUTES ALIGN:64::iglobvec
      integer iglobvec(nionlocloop)
!DIR$ ATTRIBUTES ALIGN:64::ivec
      integer ivec(nionlocloop)
!DIR$ ATTRIBUTES ALIGN:64::pchgvec
      real(t_p) pchgvec(nionlocloop)
      character*10 mode
    
      if(rank.eq.0.and.tinkerdebug) write (*,*) 'echarge3cvec'
c
c     zero out the Ewald summation energy and partitioning
c
      
      nec = 0
      ec  = 0.0_ti_p
      e   = 0.0_ti_p

      if (nion .eq. 0)  return
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  then
        call ecrecip
        if (use_pmecore) return
      end if

      kkchgdefault = 0 !exclusion with mask by default
      iglobdefault = 1 !exclusion with mask by default
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / sqrtpi

!DIR$ ASSUME (mod(nionlocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1,nionlocloop
         if (k.le.nionloc) then
            pchgvec (k) = pchg(chgglob(k))
            iglobvec(k) = iion(chgglob(k))
         else
            pchgvec(k) = 0.0_ti_p
            iglobvec(k) = iglobdefault
         endif
         ivec(k) = loc(iglobvec(k))
         ec = ec + fs * pchgvec(k) ** 2
      enddo
      nec = nec + nionloc

c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then

         term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
         xd = 0.0_ti_p
         yd = 0.0_ti_p
         zd = 0.0_ti_p
!DIR$ ASSUME (mod(nionlocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1,nionlocloop
            xd = xd + pchgvec(k) * x (iglobvec(k))
            yd = yd + pchgvec(k) * y (iglobvec(k))
            zd = zd + pchgvec(k) * z (iglobvec(k))
         enddo

         ec = ec + term * (xd**2 + yd**2 + zd**2 )
         nec = nec + 1
      end if
c
      one  = 1.0_ti_p
      half = 0.5_ti_p
c
c     compute the real space portion of the Ewald summation
c
      MAINLOOP:
     &do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i = loc(iglob)
         if(i.eq.0) cycle MAINLOOP
         inl = chglocnl(iichg)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         fi = f * pchg(iichg)
c
c     fill vectors to compute interactions involving this atom
c
         nnchg = nelst(inl)

         if(nnchg.eq.0) cycle MAINLOOP
         ncloop8  = (int(nnchg / 8 ) + 1) * 8 ! First multiple of 8
         ncloop16 = (int(nnchg / 16) + 1) * 16! First multiple of 16
         ncloop8  = merge(nnchg,ncloop8 , mod(nnchg,8 ).eq.0)
         ncloop16 = merge(nnchg,ncloop16, mod(nnchg,16).eq.0)

!DIR$ ASSUME (mod(ncloop16,16).eq.0)
!DIR$ NOFUSION
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX=64
         do k = 1, ncloop16
            if (k.le.nnchg) then
               kkchgvec(k) = elst(k,inl)
            else
               kkchgvec(k) = kkchgdefault !good value to point to by default
            endif
         enddo
         kk = 0
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, ncloop16
            if (kkchgvec(k) /= 0) then
               kk = kk + 1
               kkchgvec1(kk) = kkchgvec(k)
            endif
         enddo
         nnchg1 = kk

         ncloop8  = (int(nnchg1 / 8 ) + 1) * 8 ! First multiple of 8
         ncloop16 = (int(nnchg1 / 16) + 1) * 16! First multiple of 16
         ncloop8  = merge(nnchg1,ncloop8 , mod(nnchg1,8 ).eq.0)
         ncloop16 = merge(nnchg1,ncloop16, mod(nnchg1,16).eq.0)

!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, ncloop16
            if(k.le.nnchg1) then 
               kglobvec1(k)  = iion(kkchgvec1(k))
            else
               kkchgvec1(k) = kkchgdefault !good value to point to by default
            endif
         enddo
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!!DIR$ NOFUSION
         do k = 1, ncloop16
            xposvec1(k)  = xi - x(kglobvec1(k))
            yposvec1(k)  = yi - y(kglobvec1(k))
            zposvec1(k)  = zi - z(kglobvec1(k))
         enddo
         call  image3dvec(xposvec1,yposvec1,zposvec1,ncloop16)
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, ncloop16
!!           rvec1(k) = (  xposvec1(k)**2 + yposvec1(k)**2
!!    &                  + zposvec1(k)**2) ** half
         enddo

         kk = 0
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
         do k = 1, ncloop16
           if (    (  xposvec1(k)**2 + yposvec1(k)**2
     &              + zposvec1(k)**2) ** half .le.off
     &       .and. k.le.nnchg1) then
              kk = kk + 1
              kkchgvec2(kk) = kkchgvec1(k)
              kglobvec (kk) = kglobvec1(k)
!             rvec(kk)      = rvec1(k)
              xposvec  (kk) = xposvec1 (k)
              yposvec  (kk) = yposvec1 (k)
              zposvec  (kk) = zposvec1 (k)
           endif
         enddo
         nnchg2 = kk

         ncloop8  = (int(nnchg2 / 8 ) + 1) * 8 ! First multiple of 8
         ncloop16 = (int(nnchg2 / 16) + 1) * 16! First multiple of 16
         ncloop8  = merge(nnchg2,ncloop8 , mod(nnchg2,8 ).eq.0)
         ncloop16 = merge(nnchg2,ncloop16, mod(nnchg2,16).eq.0)

c
c     set exclusion coefficients for connected atoms
c
        call setscale(iglob,kglobvec,ncloop16,'cscale',scalevec)


!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, ncloop16
            rvec(k) = (  xposvec(k)**2 + yposvec(k)**2
     &                  + zposvec(k)**2) ** half
         enddo
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, ncloop16
            if (k.le.nnchg2) then
               fikvec(k) = fi * pchg(kkchgvec2(k))
            else
               fikvec(k) = 0.0_ti_p
            endif 
            rewvec(k)    = aewald * rvec(k)
         enddo

         call vderfc(ncloop8,rewvec,erfvec)
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, ncloop16
            if (k.le.nnchg2) then
               term1vec(k) =  erfvec(k) + scalevec(k)   - one 
               term2vec(k) = (rvec  (k) + ebuffer) ** ( - one )
               evec(k)     = fikvec(k) * term1vec(k) * term2vec(k)
               e = e + evec(k)
               if(scalevec(k).ne.0.0_ti_p) nec = nec + 1
            endif
         enddo


      end do MAINLOOP
      ec =  ec + e ! use ec outside the previous loop, because it prevents
                   ! vectorization !
      return
      end
