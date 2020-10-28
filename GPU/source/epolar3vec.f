c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar3  --  induced dipole energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar3" calculates the induced dipole polarization energy,
c     and partitions the energy among atoms
c
c
#include "tinker_precision.h"
      subroutine epolar3vec
      implicit none
     
c
c
c     choose the method for summing over polarization interactions
c
      call epolar3cvec
      return
      end

c     #####################################################################
c     ##                                                                 ##
c     ##  subroutine epolar3cvec  --  Ewald polarization analysis; list  ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     "epolar3c" calculates the polarization energy and analysis with
c     respect to Cartesian coordinates using particle mesh Ewald and
c     a neighbor list
c
c
      subroutine epolar3cvec
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
      use polar
      use polpot
      use potent
      use mpi
      use tinheader ,only:ti_p,re_p
      use sizes

      implicit none
      integer k,iipoledefault
      real(t_p) f
      real(t_p) term,fterm
      real(t_p) xd,yd,zd,xu,yu,zu
!DIR$ ATTRIBUTES ALIGN:64:: iipolevec
      integer iipolevec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: iglobvec
      integer iglobvec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: ivec
      real(t_p) ivec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: civec
      real(t_p) civec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: divecx
      real(t_p) divecx(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: divecy
      real(t_p) divecy(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: divecz
      real(t_p) divecz(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: uivecx
      real(t_p) uivecx(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: uivecy
      real(t_p) uivecy(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: uivecz
      real(t_p) uivecz(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: uiivec
      real(t_p) uiivec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xvec
      real(t_p) xvec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: yvec
      real(t_p) yvec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zvec
      real(t_p) zvec(npolelocloop)
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'epolar3cvec'
c
c     zero out the dipole polarization energy and components
c
      nep = 0
      ep  = 0.0_ti_p
c     aep = 0.0_ti_p
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      if (use_pmecore) then
         if (polalg.eq.5) then
            call dcinduce_pme
         else
            call newinduce_pmevec
         end if
      else
         if (polalg.eq.5) then
            call dcinduce_pme2
         else
            call newinduce_pme2vec
         end if
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     &   call eprecip
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         call epreal3dvec
         iipoledefault = poleglob(npoleloc)

c
c     set the energy unit conversion factor
c
         f = electric / dielec
         term = 2.0_ti_p * aewald * aewald
         fterm = -f * aewald / sqrtpi
c
c     compute the Ewald self-energy term over all the atoms
c

!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, npolelocloop
            if(k.le.npoleloc) then
               iipolevec(k) = poleglob(k)
            else
               iipolevec(k) = iipoledefault
            endif
            iglobvec (k) = ipole(iipolevec(k))
            xvec(k) = x(iglobvec(k))
            yvec(k) = y(iglobvec(k))
            zvec(k) = z(iglobvec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k  = 1, npolelocloop
            civec(k)  = rpole( 1, iipolevec(k))
            divecx(k) = rpole( 2, iipolevec(k))
            divecy(k) = rpole( 3, iipolevec(k))
            divecz(k) = rpole( 4, iipolevec(k))
        enddo
!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            uivecx(k) = uind ( 1, iipolevec(k))
            uivecy(k) = uind ( 2, iipolevec(k))
            uivecz(k) = uind ( 3, iipolevec(k))
         enddo

!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k  = 1, npolelocloop
            if(k.le.npoleloc) then
               uiivec(k) =  divecx(k) * uivecx(k)
     &                    + divecy(k) * uivecy(k)
     &                    + divecz(k) * uivecz(k)
            else
               uiivec(k) = 0.0_ti_p
            endif
            ep = ep + uiivec(k) * term * fterm / 3.0_ti_p
         enddo

         nep = nep + npoleloc
         
c        compute the cell dipole boundary correction term
c
         if (boundary .eq. 'VACUUM') then
            term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
            xd = 0.0_ti_p
            yd = 0.0_ti_p
            zd = 0.0_ti_p
            xu = 0.0_ti_p
            yu = 0.0_ti_p
            zu = 0.0_ti_p

!DIR$ ASSUME (mod(npolelocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
            do k  = 1, npolelocloop
               if(k.le.npoleloc) then
                  xd = xd + divecx(k) + civec(k) * xvec(k)
                  yd = yd + divecy(k) + civec(k) * yvec(k)
                  zd = zd + divecz(k) + civec(k) * zvec(k)
                  xu = xu + uivecx(k)
                  yu = yu + uivecy(k)
                  zu = zu + uivecz(k)
               endif
            enddo
            ep = ep + term *( xd * xu + yd * yu + zd * zu)
            nep = nep + 1
         end if
      end if
      return
      end
c
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine epreal3dvec  --  real space polar analysis via list  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "epreal3d" calculates the induced dipole polarization energy
c     and analysis using particle mesh Ewald and a neighbor list
c
c
      subroutine epreal3dvec
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
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use mpi
      use utilvec
      use tinheader ,only:ti_p,re_p
      use vec
      use vec_polar
      use vec_mpole
      use sizes
      use timestat

      implicit none
      integer ii,k, kk,iglob,iipole
      integer nnelst,nnelst1,neloop8,neloop16
      integer kpoledefault,ne
      real(t_p) f
      real(t_p) alsq2,alsq2n
      real(t_p) pdi,pti
      real(t_p) xi,yi,zi,e
      real(t_p) ci
      real(t_p) half,one
      character*10 mode

      if(rank.eq.0.and.tinkerdebug) write(*,*) 'epreal3dvec'
      call timer_enter( timer_epreal )

c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5_ti_p * electric / dielec
      half = 0.5_ti_p
      one = 1.0_ti_p
      e = 0.0_ti_p
      ne = 0
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the dipole polarization energy component
c
      alsq2 = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p)  alsq2n = 1.0_ti_p / (sqrtpi*aewald)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole     (iipole)
         pdi    = pdamp     (iipole)
         pti    = thole     (iipole)
         xi     = x         (iglob)
         yi     = y         (iglob)
         zi     = z         (iglob)

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

         uivec(1) = uind(1,iipole)
         uivec(2) = uind(2,iipole)
         uivec(3) = uind(3,iipole)

         nnelst = nelst(ii)

         neloop8  = (int(nnelst / 8 ) + 1) * 8 ! First multiple of 8
         neloop16 = (int(nnelst / 16) + 1) * 16! First multiple of 16
         neloop8  = merge(nnelst,neloop8 , mod(nnelst,8 ).eq.0)
         neloop16 = merge(nnelst,neloop16, mod(nnelst,16).eq.0)

         kpoledefault = elst(nnelst,ii)
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX = 160
         do k = 1, neloop8
            if(k.le.nnelst)  then 
               kpolevec(k) = elst(k,ii)
            else
               kpolevec(k) = kpoledefault
            endif
            kglobvec(k) = ipole(kpolevec(k))
            xposvec(k)  = x(kglobvec(k)) - xi
            yposvec(k)  = y(kglobvec(k)) - yi
            zposvec(k)  = z(kglobvec(k)) - zi
         enddo

         call  image3dvec(xposvec,yposvec,zposvec,neloop16)
c
c     evaluate all sites within the cutoff distance
c
         kk = 0
!DIR$ ASSUME (MOD(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop16
            if (     xposvec(k)**2 + yposvec(k)**2 + zposvec(k)**2<=off2
     &          .and.k.le.nnelst) then
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

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            ckvec (k) = rpole( 1, kpolevec1(k))
            dkvecx(k) = rpole( 2, kpolevec1(k))
            dkvecy(k) = rpole( 3, kpolevec1(k))
            dkvecz(k) = rpole( 4, kpolevec1(k))
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            qkvec1(k) = rpole( 5, kpolevec1(k))
            qkvec2(k) = rpole( 6, kpolevec1(k))
            qkvec3(k) = rpole( 7, kpolevec1(k))
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            qkvec4(k) = rpole( 6, kpolevec1(k))
            qkvec5(k) = rpole( 9, kpolevec1(k))
            qkvec6(k) = rpole(10, kpolevec1(k))
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            qkvec7(k) = rpole( 7, kpolevec1(k))
            qkvec8(k) = rpole(10, kpolevec1(k))
            qkvec9(k) = rpole(13, kpolevec1(k))
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            ukvecx(k) = uind ( 1, kpolevec1(k))
            ukvecy(k) = uind ( 2, kpolevec1(k))
            ukvecz(k) = uind ( 3, kpolevec1(k))
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tholevec (k) = thole          (kpolevec1(k))
            dampvec  (k) = pdamp          (kpolevec1(k)) * pdi
            pgammavec(k) = min(pti , thole(kpolevec1(k)))!OK
         enddo

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ NOFUSION
!DIR$ VECTOR ALIGNED
         do k = 1, neloop16
            r2vec(k) = xposvec2(k)**2 + yposvec2(k)**2 + zposvec2(k)**2
            invr2vec(k)  = r2vec(k) ** ( - one )
            rvec(k)      = r2vec(k) **     half
            ralphavec(k) = aewald   * rvec(k)
         enddo
c
c     get reciprocal distance terms for this interaction
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            invrvec(k)   = r2vec(k) ** ( - half )
            rr3vec(k) =     f * invrvec(k) * invr2vec(k)
            rr5vec(k) = 3.0_ti_p * rr3vec (k) * invr2vec(k) !     3*1
            rr7vec(k) = 5.0_ti_p * rr5vec (k) * invr2vec(k) !   5*3*1
         enddo
c
c     calculate the real space Ewald error function terms
c
c
         call vderfc(neloop8,ralphavec,bn0vec)

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            bn0vec(k)  = bn0vec(k) * invrvec(k)
            exp2avec(k) = exp(-ralphavec(k)**2)
            bn1vec(k)  = (  1.0_ti_p       * bn0vec           (k)
     &                     + alsq2      * alsq2n * exp2avec(k)
     &                   )                       * invr2vec(k)
            bn2vec(k)  = (  3.0_ti_p       * bn1vec           (k)
     &                     + alsq2 ** 2 * alsq2n * exp2avec(k)
     &                   )                       * invr2vec(k)
            bn3vec(k)  = (  5.0_ti_p       * bn2vec           (k)
     &                     + alsq2 ** 3 * alsq2n * exp2avec(k)
     &                   )                       * invr2vec(k)
             
            bn0vec(k) = f * bn0vec(k)
            bn1vec(k) = f * bn1vec(k)
            bn2vec(k) = f * bn2vec(k)
            bn3vec(k) = f * bn3vec(k)
         enddo
c
c     apply Thole polarization damping to scale factors
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            if (dampvec(k) .ne. 0.0_ti_p)
     &         invdampvec(k) = dampvec(k) ** ( - one )
               dampvec1(k) = -  pgammavec(k)
     &                        * (rvec(k) * invdampvec(k)) ** 3
         end do
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
               expdampvec1(k) = exp(dampvec1(k))
               sc3vec(k) = 1.0_ti_p -  expdampvec1(k)
               sc5vec(k) = 1.0_ti_p -  (1.0_ti_p - dampvec1(k))
     &                            * expdampvec1(k)
               sc7vec(k) = 1.0_ti_p -  (  1.0_ti_p - dampvec1(k)
     &                               + 0.6_ti_p * dampvec1(k)**2)
     &                            * expdampvec1(k)
         enddo
c
c     intermediates involving moments and distance separation
c
c     calculate intermediate terms for polarization interaction
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
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
         do k = 1, neloop8
            drivec(k) =  divec (1) * xposvec2(k)
     &                 + divec (2) * yposvec2(k)
     &                 + divec (3) * zposvec2(k)
            drkvec(k) =  dkvecx(k) * xposvec2(k)
     &                 + dkvecy(k) * yposvec2(k)
     &                 + dkvecz(k) * zposvec2(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qrrivec(k)  =  qrivecx(k) * xposvec2(k)
     &                   + qrivecy(k) * yposvec2(k)
     &                   + qrivecz(k) * zposvec2(k)
            qrrkvec(k)  =  qrkvecx(k) * xposvec2(k)
     &                   + qrkvecy(k) * yposvec2(k)
     &                   + qrkvecz(k) * zposvec2(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            urivec(k) =  uivec(1)  * xposvec2(k)
     &                 + uivec(2)  * yposvec2(k)
     &                 + uivec(3)  * zposvec2(k)
            urkvec(k) =  ukvecx(k) * xposvec2(k)
     &                 + ukvecy(k) * yposvec2(k)
     &                 + ukvecz(k) * zposvec2(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            duikvec(k)  =  divec  (1) * ukvecx(k)
     &                   + divec  (2) * ukvecy(k)
     &                   + divec  (3) * ukvecz(k)
     &                   + dkvecx (k) * uivec (1)
     &                   + dkvecy (k) * uivec (2)
     &                   + dkvecz (k) * uivec (3)
            quikvec(k)  =  qrivecx(k) * ukvecx(k)
     &                   + qrivecy(k) * ukvecy(k)
     &                   + qrivecz(k) * ukvecz(k)
     &                   - qrkvecx(k) * uivec (1)
     &                   - qrkvecy(k) * uivec (2)
     &                   - qrkvecz(k) * uivec (3)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            term1vec(k) =  ckvec(k)  * urivec (k)
     &                   - ci        * urkvec (k)
     &                   + duikvec(k)
            term2vec(k) =  2.0_ti_p     * quikvec(k)
     &                   - urivec(k) * drkvec (k)
     &                   - drivec(k) * urkvec (k)
            term3vec(k) =  urivec(k) * qrrkvec(k)
     &                   - urkvec(k) * qrrivec(k)
         enddo
c   
c     set exclusion coefficients for connected atoms
c
         pscalevec = setscale2(iglob,kglobvec1,neloop16,'pscale')
c
c     modify error function terms to account for scaling
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX=64
!DIR$ SIMD

         do k = 1, neloop8
            psr3vec(k) =  bn1vec(k)
     &                  - rr3vec(k)* (1.0_ti_p - sc3vec(k)*pscalevec(k))
            psr5vec(k) =  bn2vec(k)
     &                  - rr5vec(k)* (1.0_ti_p - sc5vec(k)*pscalevec(k))
            psr7vec(k) =  bn3vec(k)
     &                  - rr7vec(k)* (1.0_ti_p - sc7vec(k)*pscalevec(k))
         enddo
c
c     compute the energy contribution for this interaction
c

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop16
            if (k.le.nnelst1) then
               e = e + term1vec(k) * psr3vec(k)
     &               + term2vec(k) * psr5vec(k)
     &               + term3vec(k) * psr7vec(k)
               if (pscalevec(k) /= 0.0_ti_p) ne = ne + 1
            endif
         enddo
c
      end do
      nep = nep + ne ! moved outside of the oreceeding loop
      ep  = ep + e ! moved outside of the oreceeding loop
                    ! as it prevents loop vectorization !
      call timer_exit( timer_epreal )
      end
