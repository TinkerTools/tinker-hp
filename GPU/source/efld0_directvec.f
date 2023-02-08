c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_macro.h"
      subroutine efld0_directvec(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst
      use atoms
      use couple
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use vec
      use vec_elec
      use utilvec
      use mpi
      use timestat
      implicit none
      integer i,iglob,nrhs,inl,iploc,iii,kk,countsel
      integer ii,j,k,iipole
      integer nnelst,nnelst1,nnelst2,neloop8,neloop16
      integer kpoledefault,kglobdefault,kbisdefault
      real(t_p) pdi,pti,ci
      real(t_p) xi,yi,zi
      real(t_p) alsq2, alsq2n
      real(t_p) cutoff2
      real(t_p) half,one
!!DIR$ ATTRIBUTES ALIGN:64::ef
      real(t_p)  ef(3,nrhs,npolebloc)
cnew  real(t_p)  ef(max(1,npolebloc),3,nrhs)
!DIR$ ATTRIBUTES ALIGN:64::ef11
      real(t_p)  ef11(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::ef21
      real(t_p)  ef21(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::ef31
      real(t_p)  ef31(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::ef12
      real(t_p)  ef12(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::ef22
      real(t_p)  ef22(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::ef32
      real(t_p)  ef32(nblocloop)
   
      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list
     $  update, try lowering nlupdate')
c
      mode = 'EWALD'
      call switch (mode)
c
      if(rank.eq.0.and.tinkerdebug)write(*,*) 'efld0_directvec'
      call timer_enter( timer_efld0_direct )

!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ SIMD
      do k = 1, nblocloop
         if(k.le.nbloc) then
            ef11(k) = ef(1,1,k)
            ef21(k) = ef(2,1,k)
            ef31(k) = ef(3,1,k)
            ef12(k) = ef(1,2,k)
            ef22(k) = ef(2,2,k)
            ef32(k) = ef(3,2,k)
cnew        ef11(k) = ef(k,1,1)
cnew        ef21(k) = ef(k,2,1)
cnew        ef31(k) = ef(k,3,1)
cnew        ef12(k) = ef(k,1,2)
cnew        ef22(k) = ef(k,2,2)
cnew        ef32(k) = ef(k,3,2)
         endif
      enddo


      cutoff2 = cut2
      half   = 0.5_ti_p
      one    = 1.0_ti_p
      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi*aewald)

      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole   = poleglobnl(ii)
         iglob    = ipole     (iipole)
         i        = loc       (iglob)
         iploc = poleloc   (iipole)
         if ((i.eq.0).or.(i.gt.nbloc)) then
           write(iout,1000)
           cycle MAINLOOP
         end if
         countsel = 0
         pdi      = pdamp(iipole)
         pti      = thole(iipole)
         xi       = x    (iglob) 
         yi       = y    (iglob) 
         zi       = z    (iglob) 
c
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

         nnelst = nelst(ii)
         if(nnelst.eq.0) cycle MAINLOOP
c
c        default values
c
         kpoledefault = elst(nnelst,ii) ! safe value to point to
         kbisdefault  = 0          ! exclude atoms by default
         kglobdefault = ipole(1)   ! safe value to point to

         neloop8  = (int(nnelst / 8 ) + 1) * 8 ! First multiple of 8
         neloop16 = (int(nnelst / 16) + 1) * 16! First multiple of 16

         neloop8  = merge(nnelst,neloop8 , mod(nnelst,8 ).eq.0)
         neloop16 = merge(nnelst,neloop16, mod(nnelst,16).eq.0)

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k =  1, neloop16
            if(k.le.nnelst)  then
               kpolevec(k) = elst(k,ii)
               kbisvec (k) = poleloc(kpolevec(k))
               kglobvec(k) = ipole  (kpolevec(k))
            else
               kpolevec(k) = kpoledefault ! good value to point to
               kbisvec (k) = kbisdefault  ! exclude atoms by default
               kglobvec(k) = kglobdefault ! good value to point to
            endif
         enddo
c           if  some values in mask are false,
c           i.e
c               kbis = 0  somewhere
c            or kbis > npolebloc   somewhere
c           write a message

!!!!!!!!!!!!!!!        if(countsel.ne.nnelst) write(iout,1000)

c           and leave them out

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k =  1, neloop16
            xposvec1(k)  = x(kglobvec(k)) - xi
            yposvec1(k)  = y(kglobvec(k)) - yi
            zposvec1(k)  = z(kglobvec(k)) - zi
         enddo
         call image3dvec(xposvec1,yposvec1,zposvec1,neloop16)

         kk = 0
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop16
            if (      kbisvec(k) /= 0.and.kbisvec(k) <= npolebloc
     &          .and. xposvec1(k)**2 + yposvec1(k)**2 + zposvec1(k)**2
     &                <=cutoff2
     &          .and. k.le.nnelst
     &         ) then
             kk = kk + 1
             kpolevec1(kk) = kpolevec(k)
             kbisvec1 (kk) = kbisvec (k)
             kglobvec1(kk) = kglobvec(k)

             xposvec(kk) = xposvec1(k)
             yposvec(kk) = yposvec1(k)
             zposvec(kk) = zposvec1(k)
           end if
         end do
         nnelst1 = kk
         if(nnelst1.eq.0) cycle MAINLOOP ! No atoms selected

         neloop8  = (int(nnelst1 / 8 ) + 1) * 8 ! First multiple of 8
         neloop16 = (int(nnelst1 / 16) + 1) * 16! First multiple of 16

         neloop8  = merge(nnelst1,neloop8 , mod(nnelst1,8 ).eq.0)
         neloop16 = merge(nnelst1,neloop16, mod(nnelst1,16).eq.0)

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            if(k.le.nnelst1) then
               ckvec(k)  = rpole( 1, kpolevec1(k))
               dkvecx(k) = rpole( 2, kpolevec1(k))
               dkvecy(k) = rpole( 3, kpolevec1(k))
               dkvecz(k) = rpole( 4, kpolevec1(k))
            endif
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            if(k.le.nnelst1) then
               qkvec1(k) = rpole( 5, kpolevec1(k))
               qkvec2(k) = rpole( 6, kpolevec1(k))
               qkvec3(k) = rpole( 7, kpolevec1(k))
            endif
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            if(k.le.nnelst1) then
               qkvec4(k) = rpole( 6, kpolevec1(k))
               qkvec5(k) = rpole( 9, kpolevec1(k))
               qkvec6(k) = rpole(10, kpolevec1(k))
            endif
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            if(k.le.nnelst1) then
               qkvec7(k) = rpole( 7, kpolevec1(k))
               qkvec8(k) = rpole(10, kpolevec1(k))
               qkvec9(k) = rpole(13, kpolevec1(k))
            endif
         enddo

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop16
            if(k.le.nnelst1) then
               tholevec (k) = thole  (kpolevec1(k))
               dampvec  (k) = pdamp  (kpolevec1(k)) * pdi
            endif
            pgammavec(k) = min(pti,tholevec (k))
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop16
            invd2vec(k) =(xposvec(k)**2 + yposvec(k)**2 + zposvec(k)**2)
     &                   ** ( - one )
            d1vec(k)    = invd2vec(k) ** ( - half )
            invd1vec(k) = invd2vec(k) **     half 
         enddo

c
c     set exclusion coefficients for connected atom
c
         call setscale(iglob,kglobvec1,neloop16,'pscale',pscalevec)
         call setscalep(iglob,kglobvec1,neloop16,'dscale',dscalevec)
c
c     calculate the error function damping terms
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            ralphavec(k) = aewald * d1vec   (k)
            exp2avec (k) = exp( - ralphavec(k)**2)
         enddo
c
         call vderfc(neloop16,ralphavec,bn0vec)

!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1, neloop16
            bn0vec(k)  = bn0vec(k) * invd1vec(k)
            bn1vec(k)  = (  1.0_ti_p   * bn0vec (k)
     &                    + alsq2   * alsq2n * exp2avec(k)
     &                   )    * invd2vec(k)
            bn2vec(k)  = (  3.0_ti_p   * bn1vec (k)
     &                    + alsq2**2 * alsq2n * exp2avec(k)
     &                   )    * invd2vec(k)
            bn3vec(k)  = (  5.0_ti_p   * bn2vec (k)
     &                    + alsq2**3 * alsq2n * exp2avec(k)
     &                   )    * invd2vec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ NOFUSION
         do k = 1,neloop16
            if (dampvec(k) .ne. 0.0_ti_p) then
               invdampvec(k) = dampvec(k) ** ( - one )
               dampvec1(k) = -  pgammavec(k)
     &                        * (d1vec(k) * invdampvec(k)) ** 3
            end if
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
              expdampvec1(k) = exp(dampvec1(k))
              sc3vec(k) = 1.0_ti_p -  expdampvec1(k)
              sc5vec(k) = 1.0_ti_p -  expdampvec1(k)
     &                           * (1.0_ti_p - dampvec1(k))
              sc7vec(k) = 1.0_ti_p -  expdampvec1(k)
     &                           * (  1.0_ti_p - dampvec1(k)
     &                              + 0.6_ti_p * dampvec1(k)**2)
         enddo

!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            dsc3vec(k) = sc3vec(k) * dscalevec(k)
            psc3vec(k) = sc3vec(k) * pscalevec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            dsc5vec(k) = sc5vec(k) * dscalevec(k)
            psc5vec(k) = sc5vec(k) * pscalevec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            dsc7vec(k) = sc7vec(k) * dscalevec(k)
            psc7vec(k) = sc7vec(k) * pscalevec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            drr3vec(k)=          (1.0_ti_p - dsc3vec(k)) * invd1vec(k)
     &                         * invd2vec(k)
            prr3vec(k)=          (1.0_ti_p - psc3vec(k)) * invd1vec(k)
     &                         * invd2vec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            drr5vec(k)=  3.0_ti_p * (1.0_ti_p - dsc5vec(k))* invd1vec(k)
     &                         * invd2vec(k)* invd2vec(k)
            prr5vec(k)=  3.0_ti_p * (1.0_ti_p - psc5vec(k))* invd1vec(k)
     &                         * invd2vec(k)* invd2vec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            drr7vec(k)= 15.0_ti_p * (1.0_ti_p - dsc7vec(k))* invd1vec(k)
     &                         * invd2vec(k)* invd2vec(k)* invd2vec(k)
            prr7vec(k)= 15.0_ti_p * (1.0_ti_p - psc7vec(k))* invd1vec(k)
     &                         * invd2vec(k)* invd2vec(k)* invd2vec(k)
         enddo

c
c     compute some intermediate quantities
c
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            dirvec(k)  =   divec(1) * xposvec(k)
     &                   + divec(2) * yposvec(k)
     &                   + divec(3) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            qirvecx(k) =   qivec(1) * xposvec(k)
     &                   + qivec(2) * yposvec(k)
     &                   + qivec(3) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            qirvecy(k) =   qivec(4) * xposvec(k)
     &                   + qivec(5) * yposvec(k)
     &                   + qivec(6) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            qirvecz(k) =   qivec(7) * xposvec(k)
     &                   + qivec(8) * yposvec(k)
     &                   + qivec(9) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            qirrvec(k) =   qirvecx(k) * xposvec(k)
     &                   + qirvecy(k) * yposvec(k)
     &                   + qirvecz(k) * zposvec(k)

            dkrvec(k)  =   dkvecx(k)  * xposvec(k)
     &                   + dkvecy(k)  * yposvec(k)
     &                   + dkvecz(k)  * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            qkrvecx(k) =  qkvec1(k)  * xposvec(k)
     &                  + qkvec2(k)  * yposvec(k)
     &                  + qkvec3(k)  * zposvec(k)
            qkrvecy(k) =  qkvec4(k)  * xposvec(k)
     &                  + qkvec5(k)  * yposvec(k)
     &                  + qkvec6(k)  * zposvec(k)
            qkrvecz(k) =  qkvec7(k)  * xposvec(k)
     &                  + qkvec8(k)  * yposvec(k)
     &                  + qkvec9(k)  * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            qkrrvec(k)  =  qkrvecx(k) * xposvec(k)
     &                   + qkrvecy(k) * yposvec(k)
     &                   + qkrvecz(k) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fimvecx(k) = -      (  bn1vec(k) * ckvec  (k)
     &                           - bn2vec(k) * dkrvec (k)
     &                           + bn3vec(k) * qkrrvec(k)
     &                          )            * xposvec(k)
     &                   -         bn1vec(k) * dkvecx (k)
     &                   + 2.0_ti_p * bn2vec(k) * qkrvecx(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fimvecy(k) = -      (  bn1vec(k) * ckvec  (k)
     &                           - bn2vec(k) * dkrvec (k)
     &                           + bn3vec(k) * qkrrvec(k)
     &                          )            * yposvec(k)
     &                   -         bn1vec(k) * dkvecy (k)
     &                   + 2.0_ti_p * bn2vec(k) * qkrvecy(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fimvecz(k) = -      (  bn1vec(k) * ckvec  (k)
     &                           - bn2vec(k) * dkrvec (k)
     &                           + bn3vec(k) * qkrrvec(k)
     &                          )            * zposvec(k)
     &                   -         bn1vec(k) * dkvecz (k)
     &                   + 2.0_ti_p * bn2vec(k) * qkrvecz(k)
         enddo
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fkmvecx(k) =        (  bn1vec(k) * ci                
     &                           + bn2vec(k) * dirvec (k)
     &                           + bn3vec(k) * qirrvec(k)
     &                          )            * xposvec(k)
     &                   -         bn1vec(k) * divec  (1)
     &                   - 2.0_ti_p * bn2vec(k) * qirvecx(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fkmvecy(k) =        (  bn1vec(k) * ci            
     &                           + bn2vec(k) * dirvec (k)
     &                           + bn3vec(k) * qirrvec(k)
     &                          )            * yposvec(k)
     &                   -         bn1vec(k) * divec  (2)
     &                   - 2.0_ti_p * bn2vec(k) * qirvecy(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fkmvecz(k) =        (  bn1vec(k) * ci            
     &                           + bn2vec(k) * dirvec (k)
     &                           + bn3vec(k) * qirrvec(k)
     &                          )            * zposvec(k)
     &                   -         bn1vec(k) * divec  (3)
     &                   - 2.0_ti_p * bn2vec(k) * qirvecz(k)
         enddo
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fidvecx(k) = -      (  drr3vec(k) * ckvec  (k)
     &                           - drr5vec(k) * dkrvec (k)
     &                           + drr7vec(k) * qkrrvec(k)
     &                          )             * xposvec(k)
     &                   -         drr3vec(k) * dkvecx (k)
     &                   + 2.0_ti_p * drr5vec(k) * qkrvecx(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fidvecy(k) = -      (  drr3vec(k) * ckvec  (k)
     &                           - drr5vec(k) * dkrvec (k)
     &                           + drr7vec(k) * qkrrvec(k)
     &                          )             * yposvec(k)
     &                   -         drr3vec(k) * dkvecy (k)
     &                   + 2.0_ti_p * drr5vec(k) * qkrvecy(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fidvecz(k) = -      (  drr3vec(k) * ckvec  (k)
     &                           - drr5vec(k) * dkrvec (k)
     &                           + drr7vec(k) * qkrrvec(k)
     &                          )             * zposvec(k)
     &                   -         drr3vec(k) * dkvecz (k)
     &                   + 2.0_ti_p * drr5vec(k) * qkrvecz(k)
         enddo
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fkdvecx(k) =        (  drr3vec(k) * ci            
     &                           + drr5vec(k) * dirvec (k)
     &                           + drr7vec(k) * qirrvec(k)
     &                          )             * xposvec(k)
     &                   -         drr3vec(k) * divec  (1)
     &                   - 2.0_ti_p * drr5vec(k) * qirvecx(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fkdvecy(k) =        (  drr3vec(k) * ci            
     &                           + drr5vec(k) * dirvec (k)
     &                           + drr7vec(k) * qirrvec(k)
     &                          )             * yposvec(k)
     &                   -         drr3vec(k) * divec  (2)
     &                   - 2.0_ti_p * drr5vec(k) * qirvecy(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fkdvecz(k) =        (  drr3vec(k) * ci            
     &                           + drr5vec(k) * dirvec (k)
     &                           + drr7vec(k) * qirrvec(k)
     &                          )             * zposvec(k)
     &                   -         drr3vec(k) * divec  (3)
     &                   - 2.0_ti_p * drr5vec(k) * qirvecz(k)
         enddo
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fipvecx(k) = -      (  prr3vec(k) * ckvec  (k)
     &                           - prr5vec(k) * dkrvec (k)
     &                           + prr7vec(k) * qkrrvec(k)
     &                          )             * xposvec(k)
     &                   -         prr3vec(k) * dkvecx (k)
     &                   + 2.0_ti_p * prr5vec(k) * qkrvecx(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fipvecy(k) = -      (  prr3vec(k) * ckvec  (k)
     &                           - prr5vec(k) * dkrvec (k)
     &                           + prr7vec(k) * qkrrvec(k)
     &                          )             * yposvec(k)
     &                   -         prr3vec(k) * dkvecy (k)
     &                   + 2.0_ti_p * prr5vec(k) * qkrvecy(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fipvecz(k) = -      (  prr3vec(k) * ckvec  (k)
     &                           - prr5vec(k) * dkrvec (k)
     &                           + prr7vec(k) * qkrrvec(k)
     &                          )             * zposvec(k)
     &                   -         prr3vec(k) * dkvecz (k)
     &                   + 2.0_ti_p * prr5vec(k) * qkrvecz(k)
         enddo
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fkpvecx(k) =        (  prr3vec(k) * ci            
     &                           + prr5vec(k) * dirvec (k)
     &                           + prr7vec(k) * qirrvec(k)
     &                          )             * xposvec(k)
     &                   -         prr3vec(k) * divec  (1)
     &                   - 2.0_ti_p * prr5vec(k) * qirvecx(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fkpvecy(k) =        (  prr3vec(k) * ci            
     &                           + prr5vec(k) * dirvec (k)
     &                           + prr7vec(k) * qirrvec(k)
     &                          )             * yposvec(k)
     &                   -         prr3vec(k) * divec  (2)
     &                   - 2.0_ti_p * prr5vec(k) * qirvecy(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            fkpvecz(k) =        (  prr3vec(k) * ci 
     &                           + prr5vec(k) * dirvec (k)
     &                           + prr7vec(k) * qirrvec(k)
     &                          )             * zposvec(k)
     &                   -         prr3vec(k) * divec  (3)
     &                   - 2.0_ti_p * prr5vec(k) * qirvecz(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop16
            if (k.le.nnelst1) then
               ef11(iploc) = ef11(iploc) + fimvecx(k) - fidvecx(k)
               ef21(iploc) = ef21(iploc) + fimvecy(k) - fidvecy(k)
               ef31(iploc) = ef31(iploc) + fimvecz(k) - fidvecz(k)
               ef12(iploc) = ef12(iploc) + fimvecx(k) - fipvecx(k)
               ef22(iploc) = ef22(iploc) + fimvecy(k) - fipvecy(k)
               ef32(iploc) = ef32(iploc) + fimvecz(k) - fipvecz(k)
            endif
         enddo
         
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            if (k.le.nnelst1) then
               ef11(kbisvec1(k)) =  ef11(kbisvec1(k))
     &                            + fkmvecx(k) - fkdvecx(k)
               ef21(kbisvec1(k)) =  ef21(kbisvec1(k))
     &                            + fkmvecy(k) - fkdvecy(k)
               ef31(kbisvec1(k)) =  ef31(kbisvec1(k))
     &                            + fkmvecz(k) - fkdvecz(k)
               ef12(kbisvec1(k)) =  ef12(kbisvec1(k))
     &                            + fkmvecx(k) - fkpvecx(k)
               ef22(kbisvec1(k)) =  ef22(kbisvec1(k))
     &                            + fkmvecy(k) - fkpvecy(k)
               ef32(kbisvec1(k)) =  ef32(kbisvec1(k))
     &                            + fkmvecz(k) - fkpvecz(k)
            endif
         enddo

      end do MAINLOOP
c
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, nblocloop
         if (k.le.npolebloc) then
            ef(1,1,k) = ef11(k)
            ef(2,1,k) = ef21(k)
            ef(3,1,k) = ef31(k)
            ef(1,2,k) = ef12(k)
            ef(2,2,k) = ef22(k)
            ef(3,2,k) = ef32(k)
cnew        ef(k,1,1) = ef11(k)
cnew        ef(k,2,1) = ef21(k)
cnew        ef(k,3,1) = ef31(k)
cnew        ef(k,1,2) = ef12(k)
cnew        ef(k,2,2) = ef22(k)
cnew        ef(k,3,2) = ef32(k)
         endif
      enddo

      call timer_exit( timer_efld0_direct )

      return
      end
c
