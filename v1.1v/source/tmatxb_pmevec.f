c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
      subroutine tmatxb_pmevec(nrhs,dodiag,mu,efi)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use atmlst
      use atoms
      use domdec
      use ewald
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use mpi
      use vec
      use vec_polar
      use utilvec
      implicit none
      integer ii,j,k,kk,irhs
      integer i,iglob,iipole,nrhs
      integer kpoledefault,kpolelocdefault,kglobdefault
      integer nnelst,nnelst1,nnelst2,neloop8,neloop16
      integer iploc
      real*8 cutoff2
      real*8 pdi, pti ,pol
      real*8 alsq2, alsq2n,half,one
      real*8 xi,yi,zi
      real*8 time0, time1
!!DIR$ ATTRIBUTES ALIGN:64::mu,efi
      real*8  mu(3,nrhs,*), efi(3,nrhs,npolebloc)
!DIR$ ATTRIBUTES ALIGN:64::mu11,efi11
      real*8  mu11(nblocloop), efi11(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::mu21,efi21
      real*8  mu21(nblocloop), efi21(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::mu31,efi31
      real*8  mu31(nblocloop), efi31(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::mu12,efi12
      real*8  mu12(nblocloop), efi12(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::mu22,efi22
      real*8  mu22(nblocloop), efi22(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::mu32,efi32
      real*8  mu32(nblocloop), efi32(nblocloop)

      logical dodiag
      character*6 mode
    

      time0 = MPI_WTIME()

c
c     initialize the mu and the result vectors
c
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
      do k = 1, nblocloop
         if(k.le.npolebloc) then
            mu11(k) = mu(1,1,k)
            mu21(k) = mu(2,1,k)
            mu31(k) = mu(3,1,k)
            mu12(k) = mu(1,2,k)
            mu22(k) = mu(2,2,k)
            mu32(k) = mu(3,2,k)
         endif
         efi11(k) = 0.0d0
         efi12(k) = 0.0d0
         efi21(k) = 0.0d0
         efi22(k) = 0.0d0
         efi31(k) = 0.0d0
         efi32(k) = 0.0d0
      enddo

      half   = 0.5d0
      one    = 1.0d0
      alsq2  = 2.0d0 * aewald**2
      alsq2n = 0.0d0
      if (aewald .gt. 0.0d0) alsq2n = 1.0d0 / (sqrtpi*aewald)
c
c     gather some parameters, then set up the damping factors.
c
      mode    = 'EWALD'
      call switch (mode)
      cutoff2 = cut2
c
      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
c
c       skip atom if it is not polarizable
c
         if (polarity(iipole).eq.0.0d0) cycle MAINLOOP
         iglob = ipole  (iipole)
         i     = loc    (iglob)
         iploc = poleloc(iipole)
         if (i.eq.0) cycle MAINLOOP

         pdi = pdamp(iipole)
         pti = thole(iipole)
         xi  = x(iglob)
         yi  = y(iglob)
         zi  = z(iglob)
c
         nnelst   = nelst(ii)
         if(nnelst.eq.0) cycle MAINLOOP

         neloop8  = (int(nnelst / 8 ) + 1) * 8 ! First multiple of 8
         neloop16 = (int(nnelst / 16) + 1) * 16! First multiple of 16

         neloop8  = merge(nnelst,neloop8 , mod(nnelst,8 ).eq.0)
         neloop16 = merge(nnelst,neloop16, mod(nnelst,16).eq.0)

c
c        default values
c
         kpoledefault    = poleloc(nnelst) ! safe value to point to
         kpolelocdefault = 0               ! exclude value by default
         kglobdefault    = ipole(1)        ! safe value to point to

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k =  1, neloop16
            if(k.le.nnelst)  then 
               kpolevec   (k) = elst(k,ii)
               kpolelocvec(k) = poleloc(kpolevec(k))
               kglobvec   (k) = ipole  (kpolevec(k))
            else
               kpolevec   (k) = kpoledefault
               kpolelocvec(k) = kpolelocdefault
               kglobvec   (k) = kglobdefault
            endif
            xposvec1(k)  = x    (kglobvec(k)) - xi
            yposvec1(k)  = y    (kglobvec(k)) - yi
            zposvec1(k)  = z    (kglobvec(k)) - zi
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         call image3dvec(xposvec1,yposvec1,zposvec1,neloop16)

         kk = 0
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop16
            if (     kpolelocvec(k).ne.0
     &         .and. xposvec1(k)**2 + yposvec1(k)**2 + zposvec1(k)**2
     &               <=cutoff2
     &         .and.k.le.nnelst
     &         ) then
               kk = kk + 1
               kpolevec1   (kk) = kpolevec   (k)
               kglobvec1   (kk) = kglobvec   (k)
               kpolelocvec1(kk) = kpolelocvec(k)

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


!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1,neloop16
            if(k.le.nnelst1) then
               duivecx(k) = mu11(iploc)
               dukvecx(k) = mu11(kpolelocvec1(k))
               duivecy(k) = mu21(iploc)
               dukvecy(k) = mu21(kpolelocvec1(k))
               duivecz(k) = mu31(iploc)
               dukvecz(k) = mu31(kpolelocvec1(k))
            endif
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED 
!DIR$ NOFUSION
         do k = 1,neloop16
            if(k.le.nnelst1) then
               puivecx(k) = mu12(iploc)
               pukvecx(k) = mu12(kpolelocvec1(k))
               puivecy(k) = mu22(iploc)
               pukvecy(k) = mu22(kpolelocvec1(k))
               puivecz(k) = mu32(iploc)
               pukvecz(k) = mu32(kpolelocvec1(k))
            endif
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED 
!DIR$ SIMD
         do k = 1,neloop16
            if(k.le.nnelst1) then
               tholevec (k) = thole   (kpolevec1(k))
               dampvec  (k) = pdamp   (kpolevec1(k)) * pdi
            endif
            pgammavec(k) = min (pti,tholevec (k))
         enddo
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
          invd2vec(k) = (xposvec (k)**2 + yposvec(k)**2 + zposvec(k)**2)
     &                  ** ( - one )
           d1vec (k)   = invd2vec   (k) ** ( -  half )
           invd1vec(k) = invd2vec   (k) **  half 
         enddo
c
c        set exclusion coefficients for connected atoms
c
          call setscalep(iglob,kglobvec1,neloop16,'uscale',uscalevec)

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
           ralphavec(k) = aewald * d1vec   (k)
           exp2avec (k) = exp( - ralphavec(k)**2)
        enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
        call vderfc(neloop8,ralphavec,bn0vec)
!DIR$ ASSUME (mod(neloop16,16).eq.0)
        do k = 1,neloop16
         bn0vec(k)  = bn0vec(k)   * invd1vec(k)
         bn1vec(k)  = (  1.0d0    * bn0vec  (k)
     &                 + alsq2    * alsq2n     * exp2avec(k)
     &                )           * invd2vec(k)
         bn2vec(k)  = (  3.0d0    * bn1vec  (k)
     &                 + alsq2**2 * alsq2n     * exp2avec(k)
     &                )           * invd2vec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            if (dampvec(k) .ne. 0.0d0)
     &         invdampvec(k) = dampvec(k) ** ( - one )

               dampvec1(k) = -  pgammavec(k) 
     &                        * (d1vec(k) * invdampvec(k)) ** 3
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
               expdampvec1(k) = exp(dampvec1(k))
c
c     compute the field.
c
            rr3vec(k) =         (1.0d0 -
     &                           (1.0d0 - expdampvec1(k))
     &                         * uscalevec(k))
     &                         * invd1vec(k) * invd2vec(k)
            rr5vec(k) = 3.0d0 * (1.0d0 -
     &                 (1.0d0 - expdampvec1(k) * (1.0d0 - dampvec1(k)))
     &                         * uscalevec(k))
     &                         * invd1vec(k) * invd2vec(k) * invd2vec(k)
         enddo
c      
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1,neloop16
            duirvec(k) =  duivecx(k) * xposvec(k)
     &                  + duivecy(k) * yposvec(k)
     &                  + duivecz(k) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            dukrvec(k) =  dukvecx(k) * xposvec(k)
     &                  + dukvecy(k) * yposvec(k)
     &                  + dukvecz(k) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            puirvec(k) =  puivecx(k) * xposvec(k)
     &                  + puivecy(k) * yposvec(k)
     &                  + puivecz(k) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1,neloop16
            pukrvec(k) =  pukvecx(k) * xposvec(k)
     &                  + pukvecy(k) * yposvec(k)
     &                  + pukvecz(k) * zposvec(k)
         enddo
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fimdvecx(k) = - bn1vec(k) * dukvecx(k)
     &                    + bn2vec(k) * dukrvec(k) * xposvec(k)
            fimdvecy(k) = - bn1vec(k) * dukvecy(k)
     &                    + bn2vec(k) * dukrvec(k) * yposvec(k)
         enddo
!DIR$  ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fimdvecz(k) = - bn1vec(k) * dukvecz(k)
     &                    + bn2vec(k) * dukrvec(k) * zposvec(k)
                              
            fkmdvecx(k) = - bn1vec(k) * duivecx(k)
     &                    + bn2vec(k) * duirvec(k) * xposvec(k)
         enddo
!DIR$  ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            fkmdvecy(k) = - bn1vec(k) * duivecy(k)
     &                    + bn2vec(k) * duirvec(k) * yposvec(k)
            fkmdvecz(k) = - bn1vec(k) * duivecz(k)
     &                    + bn2vec(k) * duirvec(k) * zposvec(k)
        enddo
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
        do k = 1,neloop8
           fimpvecx(k) = - bn1vec(k) * pukvecx(k)
     &                   + bn2vec(k) * pukrvec(k) * xposvec(k)
           fimpvecy(k) = - bn1vec(k) * pukvecy(k)
     &                   + bn2vec(k) * pukrvec(k) * yposvec(k)
         enddo
!DIR$  ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
           fimpvecz(k) = - bn1vec(k) * pukvecz(k)
     &                   + bn2vec(k) * pukrvec(k) * zposvec(k)
           fkmpvecx(k) = - bn1vec(k) * puivecx(k)
     &                   + bn2vec(k) * puirvec(k) * xposvec(k)
         enddo
!DIR$  ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
           fkmpvecy(k) = - bn1vec(k) * puivecy(k)
     &                   + bn2vec(k) * puirvec(k) * yposvec(k)
           fkmpvecz(k) = - bn1vec(k) * puivecz(k)
     &                   + bn2vec(k) * puirvec(k) * zposvec(k)
        enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
        do k = 1,neloop8
           fidvecx(k) = - rr3vec(k) * dukvecx(k)
     &                  + rr5vec(k) * dukrvec(k) * xposvec(k)
           fidvecy(k) = - rr3vec(k) * dukvecy(k)
     &                  + rr5vec(k) * dukrvec(k) * yposvec(k)
        enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
        do k = 1,neloop8
           fidvecz(k) = - rr3vec(k) * dukvecz(k)
     &                  + rr5vec(k) * dukrvec(k) * zposvec(k)
           fkdvecx(k) = - rr3vec(k) * duivecx(k)
     &                  + rr5vec(k) * duirvec(k) * xposvec(k)
         enddo
!DIR$  ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
           fkdvecy(k) = - rr3vec(k) * duivecy(k)
     &                  + rr5vec(k) * duirvec(k) * yposvec(k)
           fkdvecz(k) = - rr3vec(k) * duivecz(k)
     &                  + rr5vec(k) * duirvec(k) * zposvec(k)
        enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
        do k = 1,neloop8
           fipvecx(k) = - rr3vec(k) * pukvecx(k)
     &                  + rr5vec(k) * pukrvec(k) * xposvec(k)
           fipvecy(k) = - rr3vec(k) * pukvecy(k)
     &                  + rr5vec(k) * pukrvec(k) * yposvec(k)
         enddo
!DIR$  ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
           fipvecz(k) = - rr3vec(k) * pukvecz(k)
     &                  + rr5vec(k) * pukrvec(k) * zposvec(k)
           fkpvecx(k) = - rr3vec(k) * puivecx(k)
     &                  + rr5vec(k) * puirvec(k) * xposvec(k)
         enddo
!DIR$  ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
           fkpvecy(k) = - rr3vec(k) * puivecy(k)
     &                  + rr5vec(k) * puirvec(k) * yposvec(k)
           fkpvecz(k) = - rr3vec(k) * puivecz(k)
     &                  + rr5vec(k) * puirvec(k) * zposvec(k)
        enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
        do k = 1, neloop8
           if(k.le.nnelst1) then
               efi11(iploc) = efi11(iploc) - fimdvecx(k) + fidvecx(k)
               efi21(iploc) = efi21(iploc) - fimdvecy(k) + fidvecy(k)
               efi31(iploc) = efi31(iploc) - fimdvecz(k) + fidvecz(k)
               efi12(iploc) = efi12(iploc) - fimpvecx(k) + fipvecx(k)
               efi22(iploc) = efi22(iploc) - fimpvecy(k) + fipvecy(k)
               efi32(iploc) = efi32(iploc) - fimpvecz(k) + fipvecz(k)
           endif
        enddo
c
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
        do k = 1, neloop16
           if(k.le.nnelst1) then
               efi11(kpolelocvec1(k)) =  efi11(kpolelocvec1(k))
     &                                 - fkmdvecx(k) + fkdvecx(k)
               efi21(kpolelocvec1(k)) =  efi21(kpolelocvec1(k))
     &                                 - fkmdvecy(k) + fkdvecy(k)
               efi31(kpolelocvec1(k)) =  efi31(kpolelocvec1(k))
     &                                 - fkmdvecz(k) + fkdvecz(k)

           endif
        enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
        do k = 1, neloop16
           if(k.le.nnelst1) then
               efi12(kpolelocvec1(k)) =  efi12(kpolelocvec1(k))
     &                                 - fkmpvecx(k) + fkpvecx(k)
               efi22(kpolelocvec1(k)) =  efi22(kpolelocvec1(k))
     &                                 - fkmpvecy(k) + fkpvecy(k)
               efi32(kpolelocvec1(k)) =  efi32(kpolelocvec1(k))
     &                                 - fkmpvecz(k) + fkpvecz(k)
           endif
        enddo

      end do MAINLOOP
      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.

!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do i = 1, npolelocloop
            if(i.le.npoleloc) then
               iipole = poleglob(i)
c
c     if no polarisability, take a negligeable value to allow convergence
c
               if (polarity(iipole).eq.0.d00) then
                  pol = tinypol ** -1
               else
                  pol  = polarity(iipole) ** -1
               endif
               efi11(i) =  efi11(i) + mu11(i) * pol
               efi21(i) =  efi21(i) + mu21(i) * pol
               efi31(i) =  efi31(i) + mu31(i) * pol
               efi12(i) =  efi12(i) + mu12(i) * pol
               efi22(i) =  efi22(i) + mu22(i) * pol
               efi32(i) =  efi32(i) + mu32(i) * pol
            endif
         end do
      endif

!DIR$ ASSUME (mod(nblocloop,16).eq.0)
      do k = 1, nblocloop
         if (k.le.npolebloc) then
            efi(1,1,k) = efi11(k)
            efi(2,1,k) = efi21(k)
            efi(3,1,k) = efi31(k)
            efi(1,2,k) = efi12(k)
            efi(2,2,k) = efi22(k)
            efi(3,2,k) = efi32(k)
         endif
      enddo


      time1 = MPI_WTIME()
      if(rank.eq.0.and.tinkertime) write(*,*) 'time tmatxbvec ',
     &      time1 - time0
      end
