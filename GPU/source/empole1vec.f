c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1" : driver for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      subroutine empole1vec
      use domdec
      use energi
      use potent
      use tinheader ,only:ti_p,re_p

      implicit none
c
c
c     choose the method for summing over multipole interactions
c
      call empole1cvec
c
c     zero out energy and derivative terms which are not in use
c
      if (.not. use_mpole) then
         em = 0.0_ti_p
      end if
      if (.not. use_polar) then
         ep = 0.0_ti_p
      end if
      return
      end
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1c  --  Ewald multipole derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list
c
c
      subroutine empole1cvec
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use potent
      use timestat
      use virial
      use mpi
      use tinheader ,only:ti_p,re_p
      use vec
      use vec_mpole

      implicit none
      integer iipole,iglob
      integer iipoledefault
      integer k
      real(t_p) f
      real(t_p) term,fterm
      real(t_p) xd,yd,zd
      real(t_p) xq,yq,zq
      real(t_p) xv,yv,zv,vterm
      real(t_p) xdfield,ydfield
      real(t_p) zdfield
      real(t_p) time0,time1
!DIR$ ATTRIBUTES ALIGN:64::iipolevec 
      integer iipolevec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: iglobvec
      integer iglobvec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: ivec
      integer ivec(npolelocloop)
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
!DIR$ ATTRIBUTES ALIGN:64:: fixvecx
      real(t_p) fixvecx(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: fixvecy
      real(t_p) fixvecy(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: fixvecz
      real(t_p) fixvecz(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: fiyvecx
      real(t_p) fiyvecx(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: fiyvecy
      real(t_p) fiyvecy(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: fiyvecz
      real(t_p) fiyvecz(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: fizvecx
      real(t_p) fizvecx(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: fizvecy
      real(t_p) fizvecy(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: fizvecz
      real(t_p) fizvecz(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: trqvecx
      real(t_p) trqvecx(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: trqvecy
      real(t_p) trqvecy(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: trqvecz
      real(t_p) trqvecz(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xvec
      real(t_p) xvec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: yvec
      real(t_p) yvec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: zvec
      real(t_p) zvec(npolelocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dem1vec
      real(t_p) dem1vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dem2vec
      real(t_p) dem2vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dem3vec
      real(t_p) dem3vec(nblocloop)
c
      if(rank.eq.0.and.tinkerdebug)write(*,*) 'empole1cvec'
      if (npole .eq. 0)  return
c
c     zero out the atomic multipole energy and derivatives
c
      em  = 0.0_ti_p
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
     &   then
        time0 = mpi_wtime()
        call emrecip1
        time1 = mpi_wtime()
        timerec = timerec + time1 - time0
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         time0 = mpi_wtime()
         call emreal1cvec4
         time1 = mpi_wtime()
         timereal = timereal + time1 - time0
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
         do k = 1, nblocloop
            dem1vec(k) = 0.0_ti_p
            dem2vec(k) = 0.0_ti_p
            dem3vec(k) = 0.0_ti_p
         enddo

c
c     compute the Ewald self-energy term over all the atoms
c
         iipoledefault = poleglob(npoleloc)
         term = 2.0_ti_p * aewald * aewald
         fterm = -f * aewald / sqrtpi
c
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
         do k = 1, npolelocloop
            trqvecx(k) = 0.0_ti_p
            trqvecy(k) = 0.0_ti_p
            trqvecz(k) = 0.0_ti_p
         enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, npolelocloop
            if(k.le.npoleloc) then
               iipolevec(k) = poleglob(k)
            else
               iipolevec(k) = iipoledefault
            endif
            iglobvec (k) = ipole(iipolevec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            xvec(k) = x(iglobvec(k))
            yvec(k) = y(iglobvec(k))
            zvec(k) = z(iglobvec(k))
            ivec      (k)  = loc(iglobvec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k  = 1, npolelocloop
            civec (k) = rpole( 1, iipolevec(k))
            divecx(k) = rpole( 2, iipolevec(k))
            divecy(k) = rpole( 3, iipolevec(k))
            divecz(k) = rpole( 4, iipolevec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            qivec1(k) = rpole( 5, iipolevec(k))
            qivec2(k) = rpole( 6, iipolevec(k))
            qivec3(k) = rpole( 7, iipolevec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k  = 1, npolelocloop
            qivec4(k) = rpole( 6, iipolevec(k))
            qivec5(k) = rpole( 9, iipolevec(k))
            qivec6(k) = rpole(10, iipolevec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            qivec7(k) = rpole( 7, iipolevec(k))
            qivec8(k) = rpole(10, iipolevec(k))
            qivec9(k) = rpole(13, iipolevec(k))
         enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            ciivec(k) = civec(k) ** 2
            diivec(k) = divecx(k)** 2 + divecy(k) ** 2 + divecz(k) ** 2
         enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k  = 1, npolelocloop
            qiivec(k) =  2.0_ti_p * (  qivec2(k) ** 2 + qivec3(k) ** 2
     &                            + qivec6(k) ** 2 )
     &                 + qivec1(k) ** 2
     &                 + qivec5(k) ** 2
     &                 + qivec9(k) ** 2
         enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocloop
            if (k.le.npoleloc)
     &         em =  em
     &             + fterm * (  ciivec(k)
     &                        + term * (  diivec(k) / 3.0_ti_p
     &                        + 2.0_ti_p/5.0_ti_p * term * qiivec(k)
     &                                 )
     &                       )
         enddo
c
c       compute the cell dipole boundary correction term
c

         if (boundary .eq. 'VACUUM') then
           xd = 0.0_ti_p
           yd = 0.0_ti_p
           zd = 0.0_ti_p
           term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
            do k  = 1, npolelocloop
               xd = xd + divecx(k) +  civec(k) * xvec(k)
               yd = yd + divecy(k) +  civec(k) * yvec(k)
               zd = zd + divecz(k) +  civec(k) * zvec(k)
            enddo
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
            do k  = 1, npolelocloop
               dem1vec(k) =  2.0_ti_p * term * civec(k) * xd
               dem2vec(k) =  2.0_ti_p * term * civec(k) * yd
               dem3vec(k) =  2.0_ti_p * term * civec(k) * zd
            enddo
            em = em + term * (xd * xd + yd * yd + zd * zd)
            xdfield = -2.0_ti_p * term * xd
            ydfield = -2.0_ti_p * term * yd
            zdfield = -2.0_ti_p * term * zd
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
            do k  = 1, npolelocloop
               trqvecx(k) = divecy(k) * zdfield - divecz(k) * ydfield
               trqvecy(k) = divecz(k) * xdfield - divecx(k) * zdfield
               trqvecz(k) = divecx(k) * ydfield - divecy(k) * xdfield
            enddo

            call torquevec2 ( iipolevec ,
     &                        trqvecx    ,
     &                        trqvecy    ,
     &                        trqvecz    ,
     &                        fixvecx    ,
     &                        fixvecy    ,
     &                        fixvecz    ,
     &                        fiyvecx    ,
     &                        fiyvecy    ,
     &                        fiyvecz    ,
     &                        fizvecx    ,
     &                        fizvecy    ,
     &                        fizvecz    ,
     &                        dem1vec,
     &                        dem2vec,
     &                        dem3vec
     &                      )
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
            do k  = 1, npolelocloop
               dem(1,ivec(k)) = dem(1,ivec(k)) + dem1vec(k)
               dem(2,ivec(k)) = dem(2,ivec(k)) + dem2vec(k)
               dem(3,ivec(k)) = dem(3,ivec(k)) + dem3vec(k)
            enddo
            xdfield = -2.0_ti_p * term * xd
            ydfield = -2.0_ti_p * term * yd
c
c       boundary correction to virial due to overall cell dipole
c
!DIR$ ASSUME (mod(npolelocloop,16).eq.0)
            do k  = 1, npolelocloop
               if(k.le.npoleloc) then
                  xd = xd + divecx(k)
                  yd = yd + divecy(k)
                  zd = zd + divecz(k)
                  xq = xq + civec(k) * xvec(k)
                  yq = yq + civec(k) * yvec(k)
                  zq = zq + civec(k) * zvec(k)
               endif
            enddo
            xv = xd * xq
            yv = yd * yq
            zv = zd * zq

            vterm = term * (  xd * xd + yd * yd + zd * zd
     &                      + 2.0_ti_p * (xv + yv + zv)
     &                      + xq * xq + yq * yq + zq * zq
     &                     )
            vir(1,1) = vir(1,1) + 2.0_ti_p * term * (xq * xq + xv)+vterm
            vir(2,1) = vir(2,1) + 2.0_ti_p * term * (xq * yq + xv)
            vir(3,1) = vir(3,1) + 2.0_ti_p * term * (xq * zq + xv)
            vir(1,2) = vir(1,2) + 2.0_ti_p * term * (yq * xq + yv)
            vir(2,2) = vir(2,2) + 2.0_ti_p * term * (yq * yq + yv)+vterm
            vir(3,2) = vir(3,2) + 2.0_ti_p * term * (yq * zq + yv)
            vir(1,3) = vir(1,3) + 2.0_ti_p * term * (zq * xq + zv)
            vir(2,3) = vir(2,3) + 2.0_ti_p * term * (zq * yq + zv)
            vir(3,3) = vir(3,3) + 2.0_ti_p * term * (zq * zq + zv)+vterm
         endif
      endif
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1c  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emreal1cvec4
      use atmlst
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use domdec
      use energi
      use ewald
      use inter
      use iounit
      use utilvec
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use shunt
      use virial
      use mpi
      use tinheader ,only:ti_p,re_p
      use vec
      use vec_polar
      use vec_mpole
      use timestat

      implicit none
      integer i,k,iglob
      integer ii,kk,iipole
      integer iipoledefault,nelstdefault,iglobdefault
      integer kpoledefault
      integer countsel
      integer nnelst,nnelst1,neloop8,neloop16
      integer iii,inl
      real(t_p) half
      real(t_p) one
      real(t_p) f
      real(t_p) alsq2,alsq2n
      real(t_p) xi,yi,zi
      real(t_p) ci
      real(t_p) vxx,vyy,vzz
      real(t_p) vxy,vxz,vyz
      real(t_p) dummy

!DIR$ ATTRIBUTES ALIGN:64:: xvec
      real(t_p)  xvec(maxelst,npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: yvec
      real(t_p)  yvec(maxelst,npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: zvec
      real(t_p)  zvec(maxelst,npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: nelstvec
      integer nelstvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: nelstvec1
      integer nelstvec1(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: kpolevecg
      integer kpolevecg(maxelst,npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iipolevec1
      integer iipolevec1(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iivec1
      integer iivec1(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: ivec
      integer ivec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: ivec1
      integer ivec1(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iipolevec
      integer iipolevec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iglobvec
      integer iglobvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iglobvec1
      integer iglobvec1(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iaxvec
      integer iaxvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iayvec
      integer iayvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iazvec
      integer iazvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: trqvecx
      real(t_p)  trqvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: trqvecy
      real(t_p)  trqvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: trqvecz
      real(t_p)  trqvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fixvecx
      real(t_p)  fixvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fixvecy
      real(t_p)  fixvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fixvecz
      real(t_p)  fixvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fiyvecx
      real(t_p)  fiyvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fiyvecy
      real(t_p)  fiyvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fiyvecz
      real(t_p)  fiyvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fizvecx
      real(t_p)  fizvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fizvecy
      real(t_p)  fizvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fizvecz
      real(t_p)  fizvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: xivecx
      real(t_p)  xivecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: xivecy
      real(t_p)  xivecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: xivecz
      real(t_p)  xivecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: yivecx
      real(t_p)  yivecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: yivecy
      real(t_p)  yivecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: yivecz
      real(t_p)  yivecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: zivecx
      real(t_p)  zivecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: zivecy
      real(t_p)  zivecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: zivecz
      real(t_p)  zivecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: temx
      real(t_p) temx(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: temy
      real(t_p) temy(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: temz
      real(t_p) temz(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dem1vec
      real(t_p) dem1vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dem2vec
      real(t_p) dem2vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dem3vec
      real(t_p) dem3vec(nblocloop)

      character*10 mode


      if(rank.eq.0.and.tinkerdebug)write(*,*) 'emreal1cvec4'

      call timer_enter( timer_emreal )
      
 1000 format(' Warning, system moved too much since last neighbor list',
     $   'update, try lowering nlupdate')

c
c     initialize torque arrays
c
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
      do k = 1,nblocloop
        dem1vec(k) = 0.0_ti_p
        dem2vec(k) = 0.0_ti_p
        dem3vec(k) = 0.0_ti_p
        temx(k)    = 0.0_ti_p
        temy(k)    = 0.0_ti_p
        temz(k)    = 0.0_ti_p
      enddo
c
c     initialize global mask for selection loop
c

c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      half     = 0.5_ti_p
      one      = 1.0_ti_p
      countsel = 0
      mode     = 'EWALD'
      call switch (mode)

      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p)  alsq2n = 1.0_ti_p / (sqrtpi*aewald)

      iipoledefault = poleglobnl(npolelocnl)! good value to point to
      nelstdefault  = nelst     (npolelocnl)! good value to point to
      iglobdefault  = ipole(iipoledefault)  ! good value to point to

c
c     compute the real space portion of the Ewald summation
c
      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         i      = loc  (iglob)
         if(i.eq.0.or.i.gt.nbloc ) then
           write(iout,1000)
           cycle MAINLOOP
         endif
         nnelst = nelst     (ii)
         if (nnelst.eq.0) cycle MAINLOOP

         neloop8  = (int(nnelst / 8 ) + 1) * 8 ! First multiple of 8
         neloop16 = (int(nnelst / 16) + 1) * 16! First multiple of 16
         neloop8  = merge(nnelst,neloop8 , mod(nnelst,8 ).eq.0)
         neloop16 = merge(nnelst,neloop16, mod(nnelst,16).eq.0)


         countsel     = 0
c        kpoledefault = elst(nnelst,ii)
         kpoledefault = 1

         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop16
            if(k.le.nnelst)  then
               kpolevec(k) = elst(k,ii)
            else
               kpolevec(k) = kpoledefault  ! good value to point to
               kpolevec1(k) = kpoledefault  ! good value to point to
            endif
            kglobvec(k)    = ipole(kpolevec(k))
            kbisvec (k)    = loc  (kglobvec(k))
         enddo

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop16
            xposvec2(k) = x(kglobvec(k)) - xi
            yposvec2(k) = y(kglobvec(k)) - yi
            zposvec2(k) = z(kglobvec(k)) - zi
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         call image3dvec(xposvec2,yposvec2,zposvec2,neloop16)

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

         kk = 0
c
c       select all sites under the cutoff
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1, neloop16
           if(      kbisvec(k).ne.0
     &        .and. xposvec2(k)**2 + yposvec2(k)**2 + zposvec2(k)**2
     &              .le.off2
     &        .and. k.le.nnelst
     &       ) then
             kk            = kk + 1
             kpolevec1(kk) = kpolevec(k)
             kglobvec1(kk) = kglobvec(k)
             kbisvec1(kk) = kbisvec(k)

             xposvec  (kk) = xposvec2   (k)
             yposvec  (kk) = yposvec2   (k)
             zposvec  (kk) = zposvec2   (k)
           end if
         end do
         nnelst1  = kk

         neloop8  = (int(nnelst1 / 8 ) + 1) * 8
         neloop16 = (int(nnelst1 / 16) + 1) * 16
         neloop8  = merge(nnelst1,neloop8 , mod(nnelst1,8 ).eq.0)
         neloop16 = merge(nnelst1,neloop16, mod(nnelst1,16).eq.0)
         
         if(nnelst1.eq.0) cycle MAINLOOP
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop16
            r2vec   (k) = xposvec(k)**2 + yposvec(k)**2 + zposvec(k)**2
            invr2vec(k) = r2vec  (k)** (- one )
            rvec    (k) = r2vec  (k)**    half
            invrvec (k) = rvec   (k)** (- one )
c
c     calculate the real space Ewald error function terms
c
            ralphavec(k) = aewald    * rvec(k)
         enddo

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop16
            if(k.le.nnelst1) then
               ckvec (k) = rpole( 1, kpolevec1(k)) !
               dkvecx(k) = rpole( 2, kpolevec1(k))
               dkvecy(k) = rpole( 3, kpolevec1(k))
               dkvecz(k) = rpole( 4, kpolevec1(k))
            endif
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop16
            if(k.le.nnelst1) then
               qkvec1(k) = rpole( 5, kpolevec1(k))
               qkvec2(k) = rpole( 6, kpolevec1(k))
               qkvec3(k) = rpole( 7, kpolevec1(k))
               qkvec4(k) = rpole( 6, kpolevec1(k))
               qkvec5(k) = rpole( 9, kpolevec1(k))
               qkvec6(k) = rpole(10, kpolevec1(k))
               qkvec7(k) = rpole( 7, kpolevec1(k))
               qkvec8(k) = rpole(10, kpolevec1(k))
               qkvec9(k) = rpole(13, kpolevec1(k))
            endif
         enddo

c
c     get reciprocal distance terms for this interaction
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            rr1vec(k)  =     f             * invrvec(k) 
            rr3vec(k)  =         rr1vec(k) * invr2vec(k)
            rr5vec(k)  = 3.0_ti_p * rr3vec(k) * invr2vec(k)
            rr7vec(k)  = 5.0_ti_p * rr5vec(k) * invr2vec(k)
            rr9vec(k)  = 7.0_ti_p * rr7vec(k) * invr2vec(k)
            rr11vec(k) = 9.0_ti_p * rr9vec(k) * invr2vec(k)
         enddo

c
         call vderfc(neloop16,ralphavec,bn0vec)

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            exp2avec(k) = exp( - ralphavec(k)**2)
            bn0vec(k)   = bn0vec(k)   * invrvec          (k)
            bn1vec(k)   = (  1.0_ti_p    * bn0vec           (k)
     &                     + alsq2    * alsq2n * exp2avec(k)
     &                    )           * invr2vec         (k)
            bn2vec(k)   = (  3.0_ti_p    * bn1vec           (k)
     &                     + alsq2**2 * alsq2n * exp2avec(k)
     &                    )           * invr2vec         (k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            bn3vec(k)   = (  5.0_ti_p    * bn2vec           (k)
     &                     + alsq2**3 * alsq2n * exp2avec(k)
     &                    )           * invr2vec         (k)
            bn4vec(k)   = (  7.0_ti_p    * bn3vec           (k)
     &                     + alsq2**4 * alsq2n * exp2avec(k)
     &                    )           * invr2vec         (k)
            bn5vec(k)   = (  9.0_ti_p    * bn4vec           (k)
     &                     + alsq2**5 * alsq2n * exp2avec(k)
     &                    )           * invr2vec         (k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            bn0vec(k) = f * bn0vec(k)
            bn1vec(k) = f * bn1vec(k)
            bn2vec(k) = f * bn2vec(k)
            bn3vec(k) = f * bn3vec(k)
            bn4vec(k) = f * bn4vec(k)
            bn5vec(k) = f * bn5vec(k)
         enddo
c
c     intermediates involving moments and distance separation
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            dikvecx(k) = divec(2) * dkvecz(k) - divec(3) * dkvecy(k)
            dikvecy(k) = divec(3) * dkvecx(k) - divec(1) * dkvecz(k)
            dikvecz(k) = divec(1) * dkvecy(k) - divec(2) * dkvecx(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            dirvecx(k) = divec(2) * zposvec(k) - divec(3) * yposvec(k)
            dirvecy(k) = divec(3) * xposvec(k) - divec(1) * zposvec(k)
            dirvecz(k) = divec(1) * yposvec(k) - divec(2) * xposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            dkrvecx(k) = dkvecy(k) * zposvec(k) - dkvecz(k) * yposvec(k)
            dkrvecy(k) = dkvecz(k) * xposvec(k) - dkvecx(k) * zposvec(k)
            dkrvecz(k) = dkvecx(k) * yposvec(k) - dkvecy(k) * xposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            drivec(k) =  divec(1)  * xposvec(k) + divec(2)  * yposvec(k)
     &                 + divec(3)  * zposvec(k)
            drkvec(k) =  dkvecx(k) * xposvec(k) + dkvecy(k) * yposvec(k)
     &                 + dkvecz(k) * zposvec(k)
            dikrvec(k) =  divec(1) * dkvecx(k)  + divec(2)  * dkvecy(k)
     &                  + divec(3) * dkvecz(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qrivecx(k) =  qivec(1) * xposvec(k) + qivec(2) * yposvec(k)
     &                  + qivec(3) * zposvec(k)
            qrivecy(k) =  qivec(4) * xposvec(k) + qivec(5) * yposvec(k)
     &                  + qivec(6) * zposvec(k)
            qrivecz(k) =  qivec(7) * xposvec(k) + qivec(8) * yposvec(k)
     &                  + qivec(9) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            qrkvecx(k) = qkvec1(k) * xposvec(k)
     &                 + qkvec2(k) * yposvec(k)
     &                 + qkvec3(k) * zposvec(k)
            qrkvecy(k) = qkvec4(k) * xposvec(k)
     &                 + qkvec5(k) * yposvec(k)
     &                 + qkvec6(k) * zposvec(k)
            qrkvecz(k) = qkvec7(k) * xposvec(k)
     &                 + qkvec8(k) * yposvec(k)
     &                 + qkvec9(k) * zposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qrrivec(k)  =  qrivecx(k) * xposvec(k)
     &                   + qrivecy(k) * yposvec(k)
     &                   + qrivecz(k) * zposvec(k)
            qrrkvec(k)  =  qrkvecx(k) * xposvec(k)
     &                   + qrkvecy(k) * yposvec(k)
     &                   + qrkvecz(k) * zposvec(k)
            qrrikvec(k) =  qrivecx(k) * qrkvecx(k)
     &                   + qrivecy(k) * qrkvecy(k)
     &                   + qrivecz(k) * qrkvecz(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8

            qikvec(k)   =  2.0_ti_p * (  qivec(2) * qkvec2(k)
     &                              + qivec(3) * qkvec3(k)
     &                              + qivec(6) * qkvec6(k)
     &                             )
     &                   +            qivec(1) * qkvec1(k)
     &                   +            qivec(5) * qkvec5(k)
     &                   +            qivec(9) * qkvec9(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qrirvecx(k) =  qrivecz(k) * yposvec(k)
     &                   - qrivecy(k) * zposvec(k)
            qrirvecy(k) =  qrivecx(k) * zposvec(k)
     &                   - qrivecz(k) * xposvec(k)
            qrirvecz(k) =  qrivecy(k) * xposvec(k)
     &                   - qrivecx(k) * yposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qrkrvecx(k) =  qrkvecz(k) * yposvec(k)
     &                   - qrkvecy(k) * zposvec(k)
            qrkrvecy(k) =  qrkvecx(k) * zposvec(k)
     &                   - qrkvecz(k) * xposvec(k)
            qrkrvecz(k) =  qrkvecy(k) * xposvec(k)
     &                   - qrkvecx(k) * yposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qrrvecx(k) =  qrkvecy(k) * qrivecz(k)
     &                  - qrkvecz(k) * qrivecy(k)
            qrrvecy(k) =  qrkvecz(k) * qrivecx(k)
     &                  - qrkvecx(k) * qrivecz(k)
            qrrvecz(k) =  qrkvecx(k) * qrivecy(k)
     &                  - qrkvecy(k) * qrivecx(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qikrvecx(k) =   qivec(1) * qrkvecx(k)
     &                    + qivec(2) * qrkvecy(k)
     &                    + qivec(3) * qrkvecz(k)
            qikrvecy(k) =   qivec(4) * qrkvecx(k)
     &                    + qivec(5) * qrkvecy(k)
     &                    + qivec(6) * qrkvecz(k)
            qikrvecz(k) =   qivec(7) * qrkvecx(k)
     &                    + qivec(8) * qrkvecy(k)
     &                    + qivec(9) * qrkvecz(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qkirvecx(k) =  qkvec1(k) * qrivecx(k)
     &                   + qkvec2(k) * qrivecy(k)
     &                   + qkvec3(k) * qrivecz(k)
            qkirvecy(k) =  qkvec4(k) * qrivecx(k)
     &                   + qkvec5(k) * qrivecy(k)
     &                   + qkvec6(k) * qrivecz(k)
            qkirvecz(k) =  qkvec7(k) * qrivecx(k)
     &                   + qkvec8(k) * qrivecy(k)
     &                   + qkvec9(k) * qrivecz(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qikrrvecx(k) =  qikrvecz(k) * yposvec(k)
     &                    - qikrvecy(k) * zposvec(k)
            qikrrvecy(k) =  qikrvecx(k) * zposvec(k)
     &                    - qikrvecz(k) * xposvec(k)
            qikrrvecz(k) =  qikrvecy(k) * xposvec(k)
     &                    - qikrvecx(k) * yposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            qkirrvecx(k) =  qkirvecz(k) * yposvec(k)
     &                    - qkirvecy(k) * zposvec(k)
            qkirrvecy(k) =  qkirvecx(k) * zposvec(k)
     &                    - qkirvecz(k) * xposvec(k)
            qkirrvecz(k) =  qkirvecy(k) * xposvec(k)
     &                    - qkirvecx(k) * yposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            diqkvecx(k) =  divec(1) * qkvec1(k)
     &                   + divec(2) * qkvec2(k)
     &                   + divec(3) * qkvec3(k)
            diqkvecy(k) =  divec(1) * qkvec4(k)
     &                   + divec(2) * qkvec5(k)
     &                   + divec(3) * qkvec6(k)
            diqkvecz(k) =  divec(1) * qkvec7(k)
     &                   + divec(2) * qkvec8(k)
     &                   + divec(3) * qkvec9(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            dkqivecx(k) =  dkvecx(k) * qivec(1)
     &                   + dkvecy(k) * qivec(2)
     &                   + dkvecz(k) * qivec(3)
            dkqivecy(k) =  dkvecx(k) * qivec(4)
     &                   + dkvecy(k) * qivec(5)
     &                   + dkvecz(k) * qivec(6)
            dkqivecz(k) =  dkvecx(k) * qivec(7)
     &                   + dkvecy(k) * qivec(8)
     &                   + dkvecz(k) * qivec(9)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            diqrkvec(k) =  divec(1) * qrkvecx(k)
     &                   + divec(2) * qrkvecy(k)
     &                   + divec(3) * qrkvecz(k)

            dkqrivec(k) =  dkvecx(k) * qrivecx(k)
     &                   + dkvecy(k) * qrivecy(k)
     &                   + dkvecz(k) * qrivecz(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            diqkrvecx(k) =  diqkvecz(k) * yposvec(k)
     &                    - diqkvecy(k) * zposvec(k)
            diqkrvecy(k) =  diqkvecx(k) * zposvec(k)
     &                    - diqkvecz(k) * xposvec(k)
            diqkrvecz(k) =  diqkvecy(k) * xposvec(k)
     &                    - diqkvecx(k) * yposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            dkqirvecx(k) =  dkqivecz(k) * yposvec(k)
     &                    - dkqivecy(k) * zposvec(k)
            dkqirvecy(k) =  dkqivecx(k) * zposvec(k)
     &                    - dkqivecz(k) * xposvec(k)
            dkqirvecz(k) =  dkqivecy(k) * xposvec(k)
     &                    - dkqivecx(k) * yposvec(k)
         enddo
!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1, neloop16
            dqiqkvecx(k) =  divec (2) * qrkvecz(k)
     &                    - divec (3) * qrkvecy(k)
     &                    + dkvecy(k) * qrivecz(k)
     &                    - dkvecz(k) * qrivecy(k)
     &                    - 2.0_ti_p     * (  qivec(2) * qkvec3(k)
     &                                   + qivec(5) * qkvec6(k)
     &                                   + qivec(8) * qkvec9(k)
     &                                   - qivec(3) * qkvec2(k)
     &                                   - qivec(6) * qkvec5(k)
     &                                   - qivec(9) * qkvec8(k)
     &                                  )
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            dqiqkvecy(k) =  divec (3) * qrkvecx(k)
     &                    - divec (1) * qrkvecz(k)
     &                    + dkvecz(k) * qrivecx(k)
     &                    - dkvecx(k) * qrivecz(k)
     &                    - 2.0_ti_p     * (  qivec(3) * qkvec1(k)
     &                                   + qivec(6) * qkvec4(k)
     &                                   + qivec(9) * qkvec7(k)
     &                                   - qivec(1) * qkvec3(k)
     &                                   - qivec(4) * qkvec6(k)
     &                                   - qivec(7) * qkvec9(k)
     &                                  )
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            dqiqkvecz(k) =  divec (1) * qrkvecy(k)
     &                    - divec (2) * qrkvecx(k)
     &                    + dkvecx(k) * qrivecy(k)
     &                    - dkvecy(k) * qrivecx(k)
     &                    - 2.0_ti_p     * (  qivec(1) * qkvec2(k)
     &                                   + qivec(4) * qkvec5(k)
     &                                   + qivec(7) * qkvec8(k)
     &                                   - qivec(2) * qkvec1(k)
     &                                   - qivec(5) * qkvec4(k)
     &                                   - qivec(8) * qkvec7(k)
     &                                  )
         enddo
c
c     calculate intermediate terms for multipole energy
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            term1vec(k) =  ci        * ckvec (k)
            term2vec(k) =  ckvec (k) * drivec(k)
     &                   - ci        * drkvec(k)
     &                   + dikrvec(k)
            term3vec(k) =  ci        * qrrkvec    (k)
     &                   + ckvec(k)  * qrrivec    (k)
     &                   - drivec(k) * drkvec     (k)
     &                   + 2.0_ti_p     * (  dkqrivec(k) - diqrkvec(k)
     &                                  + qikvec  (k) )
            term4vec(k) =  drivec(k) * qrrkvec    (k)
     &                   - drkvec(k) * qrrivec    (k)
     &                   - 4.0_ti_p     * qrrikvec   (k)
            term5vec(k) = qrrivec(k) * qrrkvec    (k)
         enddo
c
c     modify distances to account for Ewald and exclusions
c

c
c      set exclusion coefficients for connected atoms
c

!DIR$ ASSUME (mod(neloop16,16).eq.0)
        call setscale(iglob,kglobvec1,neloop16,'mscale',mscalevec )
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            rr1vec(k)  = bn0vec(k) -(1.0_ti_p - mscalevec(k))*rr1vec (k)
            rr3vec(k)  = bn1vec(k) -(1.0_ti_p - mscalevec(k))*rr3vec (k)
            rr5vec(k)  = bn2vec(k) -(1.0_ti_p - mscalevec(k))*rr5vec (k)
            rr7vec(k)  = bn3vec(k) -(1.0_ti_p - mscalevec(k))*rr7vec (k)
            rr9vec(k)  = bn4vec(k) -(1.0_ti_p - mscalevec(k))*rr9vec (k)
            rr11vec(k) = bn5vec(k) -(1.0_ti_p - mscalevec(k))*rr11vec(k)
         enddo
c
c     compute the energy contributions for this interaction
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            if(k.le.nnelst1) then
               evec(k) =  term1vec(k) * rr1vec(k)
     &                  + term2vec(k) * rr3vec(k)
     &                  + term3vec(k) * rr5vec(k)
     &                  + term4vec(k) * rr7vec(k)
     &                  + term5vec(k) * rr9vec(k)
            else
               evec(k) = 0.0_ti_p
            endif
c
c     calculate intermediate terms for force and torque
c
            devec(k) =  term1vec(k) * rr3vec (k)
     &                + term2vec(k) * rr5vec (k)
     &                + term3vec(k) * rr7vec (k)
     &                + term4vec(k) * rr9vec (k)
     &                + term5vec(k) * rr11vec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            term1vec(k) = - ckvec  (k) * rr3vec (k)
     &                    + drkvec (k) * rr5vec (k)
     &                    - qrrkvec(k) * rr7vec (k)
            term2vec(k) =  ci          * rr3vec (k)
     &                   + drivec  (k) * rr5vec (k)
     &                   + qrrivec (k) * rr7vec (k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            term3vec(k) = 2.0_ti_p                * rr5vec(k)
            term4vec(k) = 2.0_ti_p *
     &                          ( - ckvec  (k) * rr5vec(k)
     &                            + drkvec (k) * rr7vec(k)
     &                            - qrrkvec(k) * rr9vec(k)
     &                          )
            term5vec(k) = 2.0_ti_p *
     &                          ( - ci         * rr5vec(k)
     &                            - drivec(k)  * rr7vec(k)
     &                            - qrrivec(k) * rr9vec(k)
     &                          )
            term6vec(k) = 4.0_ti_p                * rr7vec(k)
         enddo
c
c     compute the force components for this interaction
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            if(k.le.nnelst1) then
            frcvecx(k) =  devec   (k) * xposvec  (k)
     &                  + term1vec(k) * divec    (1)
     &                  + term2vec(k) * dkvecx   (k)
     &                  + term3vec(k) * (diqkvecx(k) - dkqivecx(k))
     &                  + term4vec(k) * qrivecx  (k)
     &                  + term5vec(k) * qrkvecx  (k)
     &                  + term6vec(k) * (qikrvecx(k) + qkirvecx(k))

            frcvecy(k) =  devec   (k) * yposvec  (k)
     &                  + term1vec(k) * divec    (2)
     &                  + term2vec(k) * dkvecy   (k)
     &                  + term3vec(k) * (diqkvecy(k) - dkqivecy(k))
     &                  + term4vec(k) * qrivecy  (k)
     &                  + term5vec(k) * qrkvecy  (k)
     &                  + term6vec(k) * (qikrvecy(k) + qkirvecy(k))

            frcvecz(k) =  devec   (k) * zposvec  (k)
     &                  + term1vec(k) * divec    (3)
     &                  + term2vec(k) * dkvecz   (k)
     &                  + term3vec(k) * (diqkvecz(k) - dkqivecz(k))
     &                  + term4vec(k) * qrivecz  (k)
     &                  + term5vec(k) * qrkvecz  (k)
     &                  + term6vec(k) * (qikrvecz(k) + qkirvecz(k))
c
c     compute the torque components for this interaction
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
            ttmivecx(k) = - rr3vec  (k) * dikvecx   (k)
     &                    + term1vec(k) * dirvecx   (k)
     &                    + term3vec(k) * (dqiqkvecx(k) + dkqirvecx(k))
     &                    - term4vec(k) * qrirvecx  (k)
     &                    - term6vec(k) * (qikrrvecx(k) + qrrvecx (k))

            ttmivecy(k) = - rr3vec  (k) * dikvecy   (k)
     &                    + term1vec(k) * dirvecy   (k)
     &                    + term3vec(k) * (dqiqkvecy(k) + dkqirvecy(k))
     &                    - term4vec(k) * qrirvecy  (k)
     &                    - term6vec(k) * (qikrrvecy(k) + qrrvecy (k))

            ttmivecz(k) = - rr3vec  (k) * dikvecz   (k)
     &                    + term1vec(k) * dirvecz   (k)
     &                    + term3vec(k) * (dqiqkvecz(k) + dkqirvecz(k))
     &                    - term4vec(k) * qrirvecz  (k)
     &                    - term6vec(k) * (qikrrvecz(k) + qrrvecz (k))

            ttmkvecx(k) =   rr3vec  (k) * dikvecx   (k)
     &                    + term2vec(k) * dkrvecx   (k)
     &                    - term3vec(k) * (dqiqkvecx(k) + diqkrvecx(k))
     &                    - term5vec(k) * qrkrvecx  (k)
     &                    - term6vec(k) * (qkirrvecx(k) - qrrvecx (k))

            ttmkvecy(k) =   rr3vec  (k) * dikvecy   (k)
     &                    + term2vec(k) * dkrvecy   (k)
     &                    - term3vec(k) * (dqiqkvecy(k) + diqkrvecy(k))
     &                    - term5vec(k) * qrkrvecy  (k)
     &                    - term6vec(k) * (qkirrvecy(k) - qrrvecy (k))

            ttmkvecz(k) =   rr3vec  (k) * dikvecz   (k)
     &                    + term2vec(k) * dkrvecz   (k)
     &                    - term3vec(k) * (dqiqkvecz(k) + diqkrvecz(k))
     &                    - term5vec(k) * qrkrvecz  (k)
     &                    - term6vec(k) * (qkirrvecz(k) - qrrvecz (k))
            else
              frcvecx(k)  = 0.0_ti_p
              frcvecy(k)  = 0.0_ti_p
              frcvecz(k)  = 0.0_ti_p
              ttmivecx(k) = 0.0_ti_p
              ttmivecy(k) = 0.0_ti_p
              ttmivecz(k) = 0.0_ti_p
              ttmkvecx(k) = 0.0_ti_p
              ttmkvecy(k) = 0.0_ti_p
              ttmkvecz(k) = 0.0_ti_p
            endif
         enddo
c
c     increment force-based gradient and torque on first site
c
c     and compute the energy for this interaction
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
               dem1vec(i) = dem1vec(i) + frcvecx(k)
               dem2vec(i) = dem2vec(i) + frcvecy(k)
               dem3vec(i) = dem3vec(i) + frcvecz(k)

               temx   (i) = temx   (i) + ttmivecx(k)
               temy   (i) = temy   (i) + ttmivecy(k)
               temz   (i) = temz   (i) + ttmivecz(k)
               em = em + evec(k)
         enddo
c
c     increment force-based gradient and torque on second site
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            if(k.le.nnelst1) then
               dem1vec(kbisvec1(k)) = dem1vec(kbisvec1(k)) - frcvecx(k)
               dem2vec(kbisvec1(k)) = dem2vec(kbisvec1(k)) - frcvecy(k)
               dem3vec(kbisvec1(k)) = dem3vec(kbisvec1(k)) - frcvecz(k)
            endif
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            if(k.le.nnelst1) then
               temx(kbisvec1(k))    = temx(kbisvec1(k)) + ttmkvecx(k)
               temy(kbisvec1(k))    = temy(kbisvec1(k)) + ttmkvecy(k)
               temz(kbisvec1(k))    = temz(kbisvec1(k)) + ttmkvecz(k)
            endif
         enddo
c
c     zero out temporary accumulation virial components
c
         vxx = 0.0_ti_p
         vxy = 0.0_ti_p
         vxz = 0.0_ti_p
         vyy = 0.0_ti_p
         vyz = 0.0_ti_p
         vzz = 0.0_ti_p
c
c     increment the virial due to pairwise Cartesian forces
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
c              em = em + evec(k)
               vxx = vxx - xposvec(k) * frcvecx(k)
               vxy = vxy - yposvec(k) * frcvecx(k)
               vxz = vxz - zposvec(k) * frcvecx(k)
               vyy = vyy - yposvec(k) * frcvecy(k)
               vyz = vyz - zposvec(k) * frcvecy(k)
               vzz = vzz - zposvec(k) * frcvecz(k)
         enddo

         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
       end do MAINLOOP
c
c     resolve site torques then increment forces and virial
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         if (k.le.npolelocnl) then
            iipolevec(k) = poleglobnl(k)
         else
            iipolevec(k) = iipoledefault
         endif
         iglobvec (k) = ipole(iipolevec(k))
         ivec     (k) = loc  (iglobvec (k))
      enddo

!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         trqvecx(k) = temx(ivec(k))
         trqvecy(k) = temy(ivec(k))
         trqvecz(k) = temz(ivec(k))
      enddo
      call torquevec2( iipolevec,
     &                 trqvecx  ,
     &                 trqvecy  ,
     &                 trqvecz  ,
     &                 fixvecx  ,
     &                 fixvecy  ,
     &                 fixvecz  ,
     &                 fiyvecx  ,
     &                 fiyvecy  ,
     &                 fiyvecz  ,
     &                 fizvecx  ,
     &                 fizvecy  ,
     &                 fizvecz  ,
     &                 dem1vec ,
     &                 dem2vec ,
     &                 dem3vec
     &               )

!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         iaxvec(k) = xaxis(iipolevec(k))
         iayvec(k) = yaxis(iipolevec(k))
         iazvec(k) = zaxis(iipolevec(k))

         if (iaxvec(k) /= 0) then
            xivecx(k) = x(iaxvec(k)) - x(iglobvec(k))     ! xix 
            yivecx(k) = y(iaxvec(k)) - y(iglobvec(k))     ! yix 
            zivecx(k) = z(iaxvec(k)) - z(iglobvec(k))     ! zix 
         else
            xivecx(k) = 0.0_ti_p
            yivecx(k) = 0.0_ti_p
            zivecx(k) = 0.0_ti_p
         endif
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         if (iayvec(k) /= 0) then
            xivecy(k) = x(iayvec(k)) - x(iglobvec(k))     ! xiy 
            yivecy(k) = y(iayvec(k)) - y(iglobvec(k))     ! yiy 
            zivecy(k) = z(iayvec(k)) - z(iglobvec(k))     ! ziy 
         else
            xivecy(k) = 0.0_ti_p
            yivecy(k) = 0.0_ti_p
            zivecy(k) = 0.0_ti_p
         endif
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         if (iazvec(k) /= 0) then
            xivecz(k) = x(iazvec(k)) - x(iglobvec(k))     ! xiz 
            yivecz(k) = y(iazvec(k)) - y(iglobvec(k))     ! yiz 
            zivecz(k) = z(iazvec(k)) - z(iglobvec(k))     ! ziz 
         else
            xivecz(k) = 0.0_ti_p
            yivecz(k) = 0.0_ti_p
            zivecz(k) = 0.0_ti_p
         endif
      enddo
      vxx = 0.0_ti_p
      vxy = 0.0_ti_p
      vxz = 0.0_ti_p
      vyy = 0.0_ti_p
      vyz = 0.0_ti_p
      vzz = 0.0_ti_p

!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vxx = vxx + xivecx(k) * fixvecx(k)
     &             + xivecy(k) * fiyvecx(k)
     &             + xivecz(k) * fizvecx(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vxy = vxy + yivecx(k) * fixvecx(k)
     &             + yivecy(k) * fiyvecx(k)
     &             + yivecz(k) * fizvecx(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         vxz = vxz + zivecx(k) * fixvecx(k)
     &             + zivecy(k) * fiyvecx(k)
     &             + zivecz(k) * fizvecx(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         vyy = vyy + yivecx(k) * fixvecy(k)
     &             + yivecy(k) * fiyvecy(k)
     &             + yivecz(k) * fizvecy(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         vyz = vyz + zivecx(k) * fixvecy(k)
     &             + zivecy(k) * fiyvecy(k)
     &             + zivecz(k) * fizvecy(k)
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         vzz = vzz + zivecx(k) * fixvecz(k)
     &             + zivecy(k) * fiyvecz(k)
     &             + zivecz(k) * fizvecz(k)
      enddo
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vyz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz

!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ SIMD
!!     do k = 1, nblocloop
!!        if(k.le.nbloc) then
      do k = 1, nbloc
c        if(k.le.nbloc) then
            dem(1,k) = dem(1,k) + dem1vec(k)
            dem(2,k) = dem(2,k) + dem2vec(k)
            dem(3,k) = dem(3,k) + dem3vec(k)
c        endif
      enddo
      
      call timer_exit( timer_emreal )
      end
