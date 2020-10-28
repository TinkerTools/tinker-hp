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
c     atomic multipole and dipole polarizability interactions,
c     and partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      subroutine empole3vec
      implicit none
c
c     choose the method for summing over multipole interactions
c
      call empole3cvec
c
      return
      end
c
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine empole3cvec     --  Ewald multipole analysis via list  ##
c     ##                                                                    ##
c     ########################################################################
c
c
c     "empole3d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list, and
c     partitions the energy among the atoms
c
c
      subroutine empole3cvec
      use action
      use analyz
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
      use tinheader ,only:ti_p,re_p
      use vec
      use vec_polar
      use vec_mpole

      implicit none
      integer k
      integer iipoledefault
      integer iipolevec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: iglobvec
      integer iglobvec(npolelocloop)

      real(t_p) e,f
      real(t_p) term,fterm
      real(t_p) xd,yd,zd
!DIR$ ATTRIBUTES ALIGN:64:: civec
      real(t_p) civec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: ciivec
      real(t_p) ciivec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: divecx
      real(t_p) divecx(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: divecy
      real(t_p) divecy(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: divecz
      real(t_p) divecz(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: diivec
      real(t_p) diivec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qivec1
      real(t_p) qivec1(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qivec2
      real(t_p) qivec2(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qivec3
      real(t_p) qivec3(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qivec4
      real(t_p) qivec4(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qivec5
      real(t_p) qivec5(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qivec6
      real(t_p) qivec6(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qivec7
      real(t_p) qivec7(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qivec8
      real(t_p) qivec8(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qivec9
      real(t_p) qivec9(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: qiivec
      real(t_p) qiivec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xvec
      real(t_p) xvec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: yvec
      real(t_p) yvec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zvec
      real(t_p) zvec(npolelocloop)
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'empole3cvec'
c
c     zero out the multipole and polarization energies
c
      nem = 0
      em  = 0.0_ti_p
      e = 0.0_ti_p
c     aem = 0_ti_p
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
     $  call emrecip
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
         call emreal3dvec
         iipoledefault = poleglob(npoleloc)

c
c     compute the self-energy part of the Ewald summation
c
         term  = 2.0_ti_p * aewald * aewald
         fterm = -f    * aewald / sqrtpi

!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!!!        do k = 1, npolelocloop
!!!           iipolevec(k) = iipoledefault
!!!        enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, npolelocloop
            if(k.le.npoleloc)  then
               iipolevec(k) = poleglob(k)
            else
               iipolevec(k) = iipoledefault
            endif
            iglobvec (k) = ipole(iipolevec(k))
         enddo

cdeb     if(rank.eq.0.and.tinkerdebug)write(*,*)'npoleloc,npolelocloop',
cdeb &             npoleloc,   npolelocloop

!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ NOFUSION
         do k  = 1, npolelocloop
            xvec(k) = x(iglobvec(k))
            yvec(k) = y(iglobvec(k))
            zvec(k) = z(iglobvec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            civec (k) = rpole( 1, iipolevec(k))
            divecx(k) = rpole( 2, iipolevec(k))
            divecy(k) = rpole( 3, iipolevec(k))
            divecz(k) = rpole( 4, iipolevec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            qivec1(k) = rpole( 5, iipolevec(k))
            qivec2(k) = rpole( 6, iipolevec(k))
            qivec3(k) = rpole( 7, iipolevec(k))
            qivec4(k) = rpole( 6, iipolevec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            qivec5(k) = rpole( 9, iipolevec(k))
            qivec6(k) = rpole(10, iipolevec(k))
            qivec7(k) = rpole( 7, iipolevec(k))
            qivec8(k) = rpole(10, iipolevec(k))
            qivec9(k) = rpole(13, iipolevec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k  = 1, npolelocloop
            ciivec(k) = civec (k) ** 2
            diivec(k) = divecx(k) ** 2 + divecy(k) ** 2 + divecz(k) ** 2
            qiivec(k) =  2.0_ti_p * (  qivec2(k) ** 2 + qivec3(k) ** 2
     &                            + qivec6(k) ** 2 )
     &                 + qivec1(k) ** 2
     &                 + qivec5(k) ** 2
     &                 + qivec9(k) ** 2
         enddo

!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ NOFUSION
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            if (k.le.npoleloc) then
                e =  e
     &              + fterm * (  ciivec(k)
     &                         + term * (  diivec(k) / 3.0_ti_p
     &                         + 2.0_ti_p/5.0_ti_p * term * qiivec(k)
     &                                  )
     &                        )
               nem = nem + 1
            endif
         enddo
         em = em + e !em out of previous loop because it prevents
                     !vectorization !
c
c       compute the cell dipole boundary correction term
c
         if (boundary .eq. 'VACUUM') then
            xd = 0.0_ti_p
            yd = 0.0_ti_p
            zd = 0.0_ti_p
            term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)

!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
            do k  = 1, npolelocloop
               xd = xd + divecx(k) +  civec(k) * xvec(k)
               yd = yd + divecy(k) +  civec(k) * yvec(k)
               zd = zd + divecz(k) +  civec(k) * zvec(k)
            enddo

            em = em + term * (xd * xd + yd * yd + zd * zd)
            nem = nem + 1
         end if
      end if
      return
      end
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emreal3dvec  --  real space mpole analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emreal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions, and partitions
c     the energy among the atoms using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emreal3dvec
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use chgpot
      use couple
      use domdec
      use energi
      use ewald
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      use tinheader ,only:ti_p,re_p
      use utilvec
      use mpi
      use vec
      use vec_polar
      use vec_mpole

      implicit none
      integer i
      integer k
      integer ii,kk,iipole
      integer iglob,nnelst,nnelst1
      integer kpoledefault
      integer neloop8,neloop16

      real(t_p) alsq2,alsq2n
      real(t_p) ci
      real(t_p) f
      real(t_p) half,one
      real(t_p) xi,yi,zi

      character*10 mode
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'emreal3dvec'
      if (npole .eq. 0)  return
c
c     set conversion factor, cutoff and switching coefficients
c
      one  = 1.0_ti_p
      half = 0.5_ti_p
      f    = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      alsq2 = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p)  alsq2n = 1.0_ti_p / (sqrtpi*aewald)

      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc(iglob)
         xi     = x(iglob)
         yi     = y(iglob)
         zi     = z(iglob)

         ci       = rpole( 1, iipole)
         divec(1) = rpole( 2, iipole)
         divec(2) = rpole( 3, iipole)
         divec(3) = rpole( 4, iipole)
c        qivec  is qixx, qixy, qixz, qixy, qiyy, qiyz,qixz, qiyz, qizz
         qivec(1) = rpole( 5, iipole)
         qivec(2) = rpole( 6, iipole)
         qivec(3) = rpole( 7, iipole)
         qivec(4) = rpole( 6, iipole)
         qivec(5) = rpole( 9, iipole)
         qivec(6) = rpole(10, iipole)
         qivec(7) = rpole( 7, iipole)
         qivec(8) = rpole(10, iipole)
         qivec(9) = rpole(13, iipole)

c
c     evaluate all sites within the cutoff distance
c
         nnelst   = nelst(ii)

         neloop8  = (int(nnelst / 8 ) + 1) * 8 ! First multiple of 8
         neloop16 = (int(nnelst / 16) + 1) * 16! First multiple of 16
         neloop8  = merge(nnelst,neloop8 , mod(nnelst,8 ).eq.0)
         neloop16 = merge(nnelst,neloop16, mod(nnelst,16).eq.0)

         kpoledefault = elst(nnelst,ii)

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED 
!DIR LOOP COUNT MAX=80
         do k = 1, neloop16
            if(k.le.nnelst) then 
               kpolevec(k) = elst(k,ii)
            else
               kpolevec(k) = kpoledefault
            endif
            kglobvec(k) = ipole(kpolevec(k))
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR NOFUSION
         do k = 1, neloop16
            xposvec(k)  = x(kglobvec(k)) - xi
            yposvec(k)  = y(kglobvec(k)) - yi
            zposvec(k)  = z(kglobvec(k)) - zi
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         call  image3dvec(xposvec,yposvec,zposvec,neloop16)
c
c     evaluate all sites within the cutoff distance
c
         kk = 0
!DIR$ ASSUME (MOD(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop16
            if (     xposvec(k)**2 + yposvec(k)**2+ zposvec(k)**2<=off2
     &          .and.k.le.nnelst) then    ! mask for r2
               kk = kk + 1
               kpolevec1(kk) = kpolevec(k)
               kglobvec1(kk) = kglobvec(k)
               xposvec2 (kk) = xposvec (k)
               yposvec2 (kk) = yposvec (k)
               zposvec2 (kk) = zposvec (k)
            endif
         enddo
         nnelst1  = kk
         neloop8  = (int(nnelst1 / 8 ) + 1) * 8 ! First multiple of 8
         neloop16 = (int(nnelst1 / 16) + 1) * 16! First multiple of 16
         neloop8  = merge(nnelst1,neloop8 , mod(nnelst1,8 ).eq.0)
         neloop16 = merge(nnelst1,neloop16, mod(nnelst1,16).eq.0)

c
c      set exclusion coefficients for connected atoms
c
c        mscalevec = setscale2(iglob, kglobvec1,neloop16,'mscale')
         call setscale(iglob, kglobvec1,neloop16,'mscale',mscalevec)
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            ckvec (k) = rpole( 1, kpolevec1(k)) 
            dkvecx(k) = rpole( 2, kpolevec1(k))
            dkvecy(k) = rpole( 3, kpolevec1(k))
            dkvecz(k) = rpole( 4, kpolevec1(k))
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            qkvec1(k) = rpole( 5, kpolevec1(k))
            qkvec2(k) = rpole( 6, kpolevec1(k))
            qkvec3(k) = rpole( 7, kpolevec1(k))
            qkvec4(k) = rpole( 6, kpolevec1(k))
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            qkvec5(k) = rpole( 9, kpolevec1(k))
            qkvec6(k) = rpole(10, kpolevec1(k))
            qkvec7(k) = rpole( 7, kpolevec1(k))
            qkvec8(k) = rpole(10, kpolevec1(k))
            qkvec9(k) = rpole(13, kpolevec1(k))
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ LOOP COUNT MAX=64
         do k = 1, neloop8
c           r2vec(k) = xposvec2(k)**2 + yposvec2(k)**2 + zposvec2(k)**2
c           invr2vec (k) = r2vec(k) ** ( - one )
            invr2vec (k) = (  xposvec2(k)**2 + yposvec2(k)**2
     &                      + zposvec2(k)**2 ) ** ( - one )
c           rvec     (k) = r2vec(k) **     half
c           rvec     (k) = invr2vec(k) **   ( - half)
c           invrvec  (k) = r2vec(k) ** ( - half )
            invrvec  (k) = invr2vec(k) **  half 
c           ralphavec(k) = aewald   * rvec(k)
            ralphavec(k) = aewald   * invrvec(k) ** ( - one )
         enddo

c
c     get reciprocal distance terms for this interaction
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            rr1vec(k) =     f             * invrvec (k)
            rr3vec(k) =         rr1vec(k) * invr2vec(k) !       1
            rr5vec(k) = 3.0_ti_p * rr3vec(k) * invr2vec(k) !     3*1
            rr7vec(k) = 5.0_ti_p * rr5vec(k) * invr2vec(k) !   5*3*1
            rr9vec(k) = 7.0_ti_p * rr7vec(k) * invr2vec(k) ! 7*5*3*1
         enddo

c
c     calculate the real space Ewald error function terms
c
c
         call vderfc(neloop8,ralphavec,bn0vec)

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k =1, neloop8
            bn0vec  (k) = bn0vec(k)     * invrvec(k)
            exp2avec(k) = exp( - ralphavec(k) ** 2)
            bn1vec  (k) = (  1.0_ti_p      * bn0vec(k)
     &                     + alsq2      * alsq2n * exp2avec(k)
     &                    )             * invr2vec(k)
            bn2vec  (k) = (  3.0_ti_p      * bn1vec(k)
     &                     + alsq2 ** 2 * alsq2n * exp2avec(k)
     &                    )             * invr2vec(k)
            bn3vec  (k) = (  5.0_ti_p      * bn2vec(k)
     &                     + alsq2 ** 3 * alsq2n * exp2avec(k)
     &                    )             * invr2vec(k)
            bn4vec  (k) = (  7.0_ti_p      * bn3vec(k)
     &                     + alsq2 ** 4 * alsq2n * exp2avec(k)
     &                    )             * invr2vec(k)
         enddo
    
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k =1, neloop8
            bn0vec(k) = f * bn0vec(k)
            bn1vec(k) = f * bn1vec(k)
            bn2vec(k) = f * bn2vec(k)
            bn3vec(k) = f * bn3vec(k)
            bn4vec(k) = f * bn4vec(k)
         enddo

c
c     intermediates involving moments and distance separation
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            qrivecx(k) =  qivec(1) * xposvec2(k)
     &                  + qivec(2) * yposvec2(k)
     &                  + qivec(3) * zposvec2(k)
            qrivecy(k) =  qivec(4) * xposvec2(k)
     &                  + qivec(5) * yposvec2(k)
     &                  + qivec(6) * zposvec2(k)
            qrivecz(k) =  qivec(7) * xposvec2(k)
     &                  + qivec(8) * yposvec2(k)
     &                  + qivec(9) * zposvec2(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            qrkvecx(k) =  qkvec1(k) * xposvec2(k)
     &                  + qkvec2(k) * yposvec2(k)
     &                  + qkvec3(k) * zposvec2(k)
            qrkvecy(k) =  qkvec4(k) * xposvec2(k)
     &                  + qkvec5(k) * yposvec2(k)
     &                  + qkvec6(k) * zposvec2(k)
            qrkvecz(k) =  qkvec7(k) * xposvec2(k)
     &                  + qkvec8(k) * yposvec2(k)
     &                  + qkvec9(k) * zposvec2(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            qikvec(k)  =  2.0_ti_p * (  qivec(2) * qkvec2(k)
     &                             + qivec(3) * qkvec3(k)
     &                             + qivec(6) * qkvec6(k)
     &                            )
     &                  +            qivec(1) * qkvec1(k)
     &                  +            qivec(5) * qkvec5(k)
     &                  +            qivec(9) * qkvec9(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k = 1, neloop8
            drivec(k)  =  divec(1)  * xposvec2(k)
     &                  + divec(2)  * yposvec2(k)
     &                  + divec(3)  * zposvec2(k)
            drkvec(k)  =  dkvecx(k) * xposvec2(k)
     &                  + dkvecy(k) * yposvec2(k)
     &                  + dkvecz(k) * zposvec2(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            qrrivec(k) =  qrivecx(k) * xposvec2(k) !OK
     &                  + qrivecy(k) * yposvec2(k) !OK
     &                  + qrivecz(k) * zposvec2(k) !OK
            qrrkvec(k) =  qrkvecx(k) * xposvec2(k) !OK
     &                  + qrkvecy(k) * yposvec2(k) !OK
     &                  + qrkvecz(k) * zposvec2(k) !OK
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            qrrikvec(k) =  qrivecx(k) * qrkvecx(k)
     &                   + qrivecy(k) * qrkvecy(k)
     &                   + qrivecz(k) * qrkvecz(k)
            diqrkvec(k) =  divec  (1) * qrkvecx(k)
     &                   + divec  (2) * qrkvecy(k)
     &                   + divec  (3) * qrkvecz(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dkqrivec(k) =  dkvecx(k) * qrivecx(k)
     &                   + dkvecy(k) * qrivecy(k)
     &                   + dkvecz(k) * qrivecz(k)
            dikvec(k)   =  divec (1) * dkvecx(k)   !OK
     &                   + divec (2) * dkvecy(k)   !OK
     &                   + divec (3) * dkvecz(k)   !OK
         enddo
c
c     calculate intermediate terms for multipole interaction
c

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop16
            term1vec(k) =  ci       * ckvec (k)
            term2vec(k) =  ckvec(k) * drivec(k)
     &                   - ci       * drkvec(k) + dikvec(k)  
         enddo

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ NOFUSION
         do k = 1, neloop16
            term3vec(k) =  ci        * qrrkvec(k)
     &                   + ckvec (k) * qrrivec(k)
     &                   - drivec(k) * drkvec (k)
     &                   + 2.0_ti_p     * (  dkqrivec(k)
     &                                  - diqrkvec(k)
     &                                  + qikvec  (k)
     &                                 )
         enddo

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k = 1, neloop16
            term4vec(k) =  drivec(k) * qrrkvec(k)
     &                   - drkvec(k) * qrrivec(k)
     &                       - 4.0_ti_p * qrrikvec(k)
            term5vec(k) = qrrivec(k) * qrrkvec(k)
         enddo
c
c     modify error function terms to account for scaling
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            rr1vec(k) = bn0vec(k) - (1.0_ti_p - mscalevec(k))* rr1vec(k)
            rr3vec(k) = bn1vec(k) - (1.0_ti_p - mscalevec(k))* rr3vec(k)
            rr5vec(k) = bn2vec(k) - (1.0_ti_p - mscalevec(k))* rr5vec(k)
            rr7vec(k) = bn3vec(k) - (1.0_ti_p - mscalevec(k))* rr7vec(k)
            rr9vec(k) = bn4vec(k) - (1.0_ti_p - mscalevec(k))* rr9vec(k)
         enddo
c
c     compute the energy contribution for this interaction
c

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ SIMD
         do k = 1, neloop16
            if (k.le.nnelst1) then
               em = em + term1vec(k) * rr1vec(k)
     &                 + term2vec(k) * rr3vec(k)
     &                 + term3vec(k) * rr5vec(k)
     &                 + term4vec(k) * rr7vec(k)
     &                 + term5vec(k) * rr9vec(k)
c
c             count interactions number with the undamped value
c
              if(mscalevec(k) /= 0.0_ti_p) nem = nem + 1
            endif
         enddo
      end do MAINLOOP
      return
      end
