c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar1  --  polarization energy & derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar1" calculates the induced dipole polarization energy
c     and derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      subroutine epolar1vec
      implicit none
c
c     choose the method for summing over polarization interactions
c
      call epolar1cvec
      return
      end

c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar1c  --  Ewald polarization derivs via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar1c" calculates the dipole polarization energy and
c     derivatives with respect to Cartesian coordinates using
c     particle mesh Ewald summation and a neighbor list
c
c
      subroutine epolar1cvec
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent
      use virial
      use mpi
      use sizes
      use timestat
      use tinheader ,only:ti_p,re_p
      use vec
      use vec_polar

      implicit none
      integer k
      integer iipoledefault
!DIR$ ATTRIBUTES ALIGN:64:: iipolevec
      integer iipolevec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iglobvec,ivec
      integer iglobvec(npolelocnlloop),ivec(npolelocnlloop)
      real(t_p) f,term,fterm,vterm,half
      real(t_p) xu,yu,zu,xup,yup,zup
      real(t_p) xq,yq,zq,xv,yv,zv
      real(t_p) xd,yd,zd
      real(t_p) time0,time1
!DIR$ ATTRIBUTES ALIGN:64::civec
      real(t_p) civec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: divecx,divecy
      real(t_p) divecx(npolelocnlloop),divecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: divecz
      real(t_p) divecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: uivecx,uivecy
      real(t_p) uivecx(npolelocnlloop),uivecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: uivecz
      real(t_p) uivecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: uipvecx,uipvecy
      real(t_p) uipvecx(npolelocnlloop),uipvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: uipvecz
      real(t_p) uipvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: upivecx,upivecy
      real(t_p) upivecx(npolelocnlloop),upivecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: upivecz
      real(t_p) upivecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64::depvecx,depvecy
      real(t_p) depvecx(nblocloop),depvecy(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::depvecz
      real(t_p) depvecz(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: xvec,yvec
      real(t_p) xvec(npolelocnlloop),yvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: zvec
      real(t_p) zvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fixvecx
      real(t_p) fixvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fixvecy
      real(t_p) fixvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fixvecz
      real(t_p) fixvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fiyvecx
      real(t_p) fiyvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fiyvecy
      real(t_p) fiyvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fiyvecz
      real(t_p) fiyvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fizvecx
      real(t_p) fizvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fizvecy
      real(t_p) fizvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: fizvecz
      real(t_p) fizvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: trqvecx
      real(t_p) trqvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: trqvecy
      real(t_p) trqvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: trqvecz
      real(t_p) trqvecz(npolelocnlloop)

      if(rank.eq.0.and.tinkerdebug)write(*,*) 'epolar1cvec'
!!     if(rank.eq.0.and.tinkerdebug)write(*,*) 'npolelocnlloop',
!!    &                                         npolelocnlloop

c
c
c     zero out the polarization energy and derivatives
c
      ep     = 0.0_ti_p
!     dep = 0.0_ti_p

!!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!!DIR$ VECTOR ALIGNED
!!    do k = 1, nbloc
!!        dep1(k) = 0.0_ti_p
!!        dep2(k) = 0.0_ti_p
!!        dep3(k) = 0.0_ti_p
!!     enddo

c
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
      half = 0.5_ti_p
c
c     compute the induced dipoles at each polarizable atom
c
      time0 = mpi_wtime()
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
      time1 = mpi_wtime()
      timepolarsolve = timepolarsolve + time1-time0
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $   then
         time0 = mpi_wtime()
         call eprecip1vec
         time1 = mpi_wtime()
         timerec = timerec + time1-time0
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
         time0 = mpi_wtime()
         call epreal1cvec
         time1 = mpi_wtime()
         timereal = timereal + time1-time0
         iipoledefault = poleglob(npoleloc)
c
c     compute the Ewald self-energy term over all the atoms
c
         term  =  2.0_ti_p * aewald * aewald
         fterm = -f     * aewald / sqrtpi

!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1,npolelocnlloop
            if(k.le.npoleloc) then
               iipolevec(k) = poleglob(k)
            else
               iipolevec(k) = iipoledefault
            endif
            iglobvec(k) = ipole(iipolevec(k))
            ivec(k)     = loc  (iglobvec (k))
         enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k  = 1, npolelocnlloop
            divecx (k) = rpole( 2, iipolevec(k))
            divecy (k) = rpole( 3, iipolevec(k))
            divecz (k) = rpole( 4, iipolevec(k))
            uivecx (k) = uind ( 1, iipolevec(k))
            uivecy (k) = uind ( 2, iipolevec(k))
            uivecz (k) = uind ( 3, iipolevec(k))
            uipvecx(k) = uinp ( 1, iipolevec(k))
            uipvecy(k) = uinp ( 2, iipolevec(k))
            uipvecz(k) = uinp ( 3, iipolevec(k))
         enddo

!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocnlloop
            if(k.le.npoleloc) then
               ep = ep + fterm * term / 3.0_ti_p 
     &                         * (  divecx(k) * uivecx(k) 
     &                            + divecy(k) * uivecy(k) 
     &                            + divecz(k) * uivecz(k)
     &                           ) 
            endif
         enddo
c
c        compute the self-energy torque term due to induced dipole
c
         term = (4.0_ti_p/3.0_ti_p) * f * aewald**3 / sqrtpi

!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocnlloop
            if(k.le.npoleloc) then
               upivecx(k) = half * (uivecx(k) + uipvecx(k))
               upivecy(k) = half * (uivecy(k) + uipvecy(k))
               upivecz(k) = half * (uivecz(k) + uipvecz(k))
            else
               upivecx(k) = 0.0_ti_p
               upivecy(k) = 0.0_ti_p
               upivecz(k) = 0.0_ti_p
            endif
         enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k  = 1, npolelocnlloop

               trqvecx(k) = term  * (  divecy(k) * upivecz(k)
     &                               - divecz(k) * upivecy(k))
               trqvecy(k) = term  * (  divecz(k) * upivecx(k)
     &                               - divecx(k) * upivecz(k))
               trqvecz(k) = term  * (  divecx(k) * upivecy(k)
     &                               - divecy(k) * upivecx(k))
         enddo
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         call torquevec2( iipolevec,
     &                    trqvecx  ,
     &                    trqvecy  ,
     &                    trqvecz  ,
     &                    fixvecx  ,
     &                    fixvecy  ,
     &                    fixvecz  ,
     &                    fiyvecx  ,
     &                    fiyvecy  ,
     &                    fiyvecz  ,
     &                    fizvecx  ,
     &                    fizvecy  ,
     &                    fizvecz  ,
     &                    dep(1,1:nbloc),
     &                    dep(2,1:nbloc),
     &                    dep(3,1:nbloc)
!!    &                    dep1(1:nbloc),
!!    &                    dep2(1:nbloc),
!!    &                    dep3(1:nbloc)
     &                   )

c
c       compute the cell dipole boundary correction term
c
         if (boundary .eq. 'VACUUM') then
            term     = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            do k  = 1, npolelocnlloop
               xvec(k) = x(iglobvec(k))
               yvec(k) = y(iglobvec(k))
               zvec(k) = z(iglobvec(k))
               civec(k) = rpole(1,iipolevec(k)) 
            enddo
!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            do k  = 1, npolelocnlloop
               xu  = xu  + uivecx (iipolevec(k))
               yu  = yu  + uivecy (iipolevec(k))
               zu  = zu  + uivecz (iipolevec(k))
               xup = xup + uipvecx(iipolevec(k))
               yup = yup + uipvecy(iipolevec(k))
               zup = zup + uipvecz(iipolevec(k))
               xd  = xd  + divecx(k) + civec(k) * xvec(k)
               yd  = yd  + divecy(k) + civec(k) * yvec(k)
               zd  = zd  + divecz(k) + civec(k) * zvec(k)
            enddo
            ep = ep + term * (  xd * xu + yd * yu + zd * zu)

!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
            do k  = 1, npolelocnlloop
               depvecx(k) = term * civec(k) * (xu + xup)
               depvecy(k) = term * civec(k) * (yu + yup)
               depvecz(k) = term * civec(k) * (zu + zup)

               trqvecx(k) = term * (  divecy(k) * (zu + zup)
     &                              - divecz(k) * (yu + yup))
               trqvecy(k) = term * (  divecz(k) * (xu + xup)
     &                              - divecx(k) * (zu + zup))
               trqvecz(k) = term * (  divecx(k) * (yu + yup)
     &                              - divecy(k) * (xu + xup))
            enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            do k  = 1, npolelocnlloop
               dep(1,ivec(k)) = dep(1,ivec(k)) + depvecx(k)
               dep(2,ivec(k)) = dep(2,ivec(k)) + depvecy(k)
               dep(3,ivec(k)) = dep(3,ivec(k)) + depvecz(k)
!!              dep1(ivec(k)) = dep1(ivec(k)) + depvecx(k)
!!              dep2(ivec(k)) = dep2(ivec(k)) + depvecy(k)
!!              dep3(ivec(k)) = dep3(ivec(k)) + depvecz(k)
            enddo
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
            call torquevec2 ( iipolevec,
     &                        trqvecx  ,
     &                        trqvecy  ,
     &                        trqvecz  ,
     &                        fixvecx  ,
     &                        fixvecy  ,
     &                        fixvecz  ,
     &                        fiyvecx  ,
     &                        fiyvecy  ,
     &                        fiyvecz  ,
     &                        fizvecx  ,
     &                        fizvecy  ,
     &                        fizvecz  ,
     &                        dep(1,1:nbloc),
     &                        dep(2,1:nbloc),
     &                        dep(3,1:nbloc)
!!    &                    dep1(1:nbloc),
!!    &                    dep2(1:nbloc),
!!    &                    dep3(1:nbloc)
     &                      )
c
c       boundary correction to virial due to overall cell dipole
c
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
            do k  = 1, npolelocnlloop
               xd = xd + divecx(k)
               yd = yd + divecy(k)
               zd = zd + divecz(k)
               xq = xq + civec(k) * xvec(k)
               yq = yq + civec(k) * yvec(k)
               zq = zq + civec(k) * zvec(k)
            enddo
            xv = xq * (xu + xup)
            yv = yq * (yu + yup)                        
            zv = zq * (zu + zup)                        

            vterm  =  term * (  xv + xu * xup + xd * (xu + xup)
     &                        + yv + yu * yup + yd * (xu + xup)
     &                        + zv + zu * zup + zd * (xu + xup))

            vir(1,1) = vir(1,1) + term * xv + vterm
            vir(2,1) = vir(2,1) + term * xv
            vir(3,1) = vir(3,1) + term * xv
            vir(1,2) = vir(1,2) + term * yv
            vir(2,2) = vir(2,2) + term * yv + vterm
            vir(3,2) = vir(3,2) + term * yv
            vir(1,3) = vir(1,3) + term * zv
            vir(2,3) = vir(2,3) + term * zv
            vir(3,3) = vir(3,3) + term * zv + vterm
         endif
      end if
      return
      end
c     #################################################################
c     ##                                                             ##
c     ## subroutine epreal1cvec -- Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal1d" evaluates the real space portion of the Ewald
c     summation energy and gradient due to dipole polarization
c     via a neighbor list
c
c
      subroutine epreal1cvec
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
      use math
      use molcul
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use virial
      use mpi
      use tinheader ,only:ti_p,re_p
      use vec
      use vec_polar
      use utilvec
      use timestat

      implicit none
      integer i,iglob,k,kk
      integer nnelst,nnelst1,neloop8,neloop16
      integer ii,iii,iipole
      integer nretain,countsel
      integer iipoledefault,nelstdefault
      integer kpoledefault
      real(t_p) alsq2,alsq2n
      real(t_p) exp2a,ralpha
      real(t_p) f
      real(t_p) half
      real(t_p) one
      real(t_p) pdi,pti
      real(t_p) xi,yi,zi
      real(t_p) ci
      real(t_p) t1,t2,t3
      real(t_p) t4,t5,t6
      real(t_p) vxx,vxy,vxz,vyy,vyz,vzz

!DIR$ ATTRIBUTES ALIGN:64:: iipolevec
      integer iipolevec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iglobvec
      integer iglobvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iaxvec
      integer iaxvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iayvec
      integer iayvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: iazvec
      integer iazvec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: ivec
       integer ivec(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dipolvecx
      real(t_p)  dipolvecx(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dipolvecy
      real(t_p)  dipolvecy(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: dipolvecz
      real(t_p)  dipolvecz(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: qipolvec1
      real(t_p)  qipolvec1(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: qipolvec2
      real(t_p)  qipolvec2(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: qipolvec3
      real(t_p)  qipolvec3(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: qipolvec4
      real(t_p)  qipolvec4(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: qipolvec5
      real(t_p)  qipolvec5(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: qipolvec6
      real(t_p)  qipolvec6(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: qipolvec7
      real(t_p)  qipolvec7(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: qipolvec8
      real(t_p)  qipolvec8(npolelocnlloop)
!DIR$ ATTRIBUTES ALIGN:64:: qipolvec9
      real(t_p)  qipolvec9(npolelocnlloop)
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
!DIR$ ATTRIBUTES ALIGN:64::ufldvec1
      real(t_p) ufldvec1(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::ufldvec2
      real(t_p) ufldvec2(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::ufldvec3
      real(t_p) ufldvec3(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::dufldvec1
      real(t_p) dufldvec1(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::dufldvec2
      real(t_p) dufldvec2(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::dufldvec3
      real(t_p) dufldvec3(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::dufldvec4
      real(t_p) dufldvec4(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::dufldvec5
      real(t_p) dufldvec5(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64::dufldvec6
      real(t_p) dufldvec6(nblocloop)

      character*10 mode
      if(rank.eq.0.and.tinkerdebug)write(*,*) 'epreal1cvec'

      call timer_enter( timer_epreal )

1000  format(' Warning, system moved too much since last neighbor list',
     $   ' update, try lowering nlupdate')

c
c     set arrays to store fields
c
!DIR$ ASSUME (mod(nblocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1,nblocloop
           ufldvec1(k) = 0.0_ti_p
           ufldvec2(k) = 0.0_ti_p
           ufldvec3(k) = 0.0_ti_p
          dufldvec1(k) = 0.0_ti_p
          dufldvec2(k) = 0.0_ti_p
          dufldvec3(k) = 0.0_ti_p
          dufldvec4(k) = 0.0_ti_p
          dufldvec5(k) = 0.0_ti_p
          dufldvec6(k) = 0.0_ti_p

      enddo

c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5_ti_p * electric / dielec
      mode = 'EWALD'
      call switch (mode)

      half     = 0.5_ti_p
      one      = 1.0_ti_p
      countsel = 0
      iipoledefault = poleglobnl(npolelocnl)
      nelstdefault  = nelst     (npolelocnl)

      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi * aewald)

      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole   = poleglobnl(ii)
         iglob    = ipole (iipole)
         i        = loc(iglob)
         if(i.eq.0.or.i.gt.nbloc ) then
           write(iout,1000)
           cycle MAINLOOP
         endif
         nnelst   = nelst     (ii)

         countsel = 0

c
c        No neighbours
c
         if (nnelst.eq.0) cycle MAINLOOP
         neloop8  = (int(nnelst / 8 ) + 1) * 8 ! First multiple of 8
         neloop16 = (int(nnelst / 16) + 1) * 16! First multiple of 16
         neloop8  = merge(nnelst,neloop8 , mod(nnelst,8 ).eq.0)
         neloop16 = merge(nnelst,neloop16, mod(nnelst,16).eq.0)

         xi     = x(iglob)
         yi     = y(iglob)
         zi     = z(iglob)

         kpoledefault = elst(nnelst,ii)

!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop16
            if (k.le.nnelst)  then
               kpolevec(k) = elst(k,ii)   
            else
               kpolevec(k) = kpoledefault
               kpolevec1(k) = kpoledefault
            endif
            kglobvec(k)    = ipole(kpolevec(k))
            kbisvec (k)    = loc  (kglobvec(k))
         enddo

!DIR$ ASSUME (mod(neloop16,16).eq.0)
         do k = 1, neloop16
            xposvec1(k) = x(kglobvec(k)) - xi
            yposvec1(k) = y(kglobvec(k)) - yi
            zposvec1(k) = z(kglobvec(k)) - zi
         enddo

         call image3dvec(xposvec1,yposvec1,zposvec1,neloop16)
    
         pdi = pdamp(iipole)
         pti = thole(iipole)
         ci  = rpole(1,iipole)
    
!DIR$ VECTOR ALIGNED
         divec(1) = rpole(2,iipole)
         divec(2) = rpole(3,iipole)
         divec(3) = rpole(4,iipole)
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
    
         uivec (1) = uind(1,iipole)
         uivec (2) = uind(2,iipole)
         uivec (3) = uind(3,iipole)
         uipvec(1) = uinp(1,iipole)
         uipvec(2) = uinp(2,iipole)
         uipvec(3) = uinp(3,iipole)
    
         kk = 0
!DIR$ ASSUME (mod(neloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop16
           if (      kbisvec(k) /= 0
     &         .and. xposvec1(k)**2 + yposvec1(k)**2 + zposvec1(k)**2
     &             <=off2 
     &         .and. k<=nnelst)
     &     then

             kk            = kk + 1
             kpolevec1(kk) = kpolevec(k)
             kglobvec1(kk) = kglobvec(k)
             kbisvec1 (kk) = kbisvec(k)    

             xpsvec   (kk) = xposvec1(k)
             ypsvec   (kk) = yposvec1(k)
             zpsvec   (kk) = zposvec1(k)
           end if
         end do
         nnelst1 = kk

         if(nnelst1.eq.0) cycle MAINLOOP

         neloop8  = (int(nnelst1 / 8 ) + 1) * 8
         neloop16 = (int(nnelst1 / 16) + 1) * 16
         neloop8  = merge(nnelst1,neloop8 , mod(nnelst1,8 ).eq.0)
         neloop16 = merge(nnelst1,neloop16, mod(nnelst1,16).eq.0)
c   
c     set exclusion coefficients for connected atoms
c   
         call setscale (iglob, kglobvec1,neloop16,'pscale',pscalevec)
         call setscalep(iglob, kglobvec1,neloop16,'dscale',dscalevec)
         call setscalep(iglob, kglobvec1,neloop16,'uscale',uscalevec)
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            if(k.le.nnelst1) then
               ckvec (k) = rpole( 1, kpolevec1(k))           
               dkvecx(k) = rpole( 2, kpolevec1(k))
               dkvecy(k) = rpole( 3, kpolevec1(k))
               dkvecz(k) = rpole( 4, kpolevec1(k))
            endif
         enddo

c
c        qkvec  is qkxx, qkxy, qkxz, etc  for all the sites
c    qkvec  1     2     3     4     5     6     7     8     9
c         qkxx, qkxy, qkxz, qkxy, qkyy, qkyz, qkxz, qkyz, qkzz

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            if(k.le.nnelst1) then
               qkvec1 (k) = rpole( 5, kpolevec1(k))
               qkvec2 (k) = rpole( 6, kpolevec1(k))
               qkvec3 (k) = rpole( 7, kpolevec1(k))
               qkvec4 (k) = rpole( 6, kpolevec1(k))
               qkvec5 (k) = rpole( 9, kpolevec1(k))
               qkvec6 (k) = rpole(10, kpolevec1(k))
               qkvec7 (k) = rpole( 7, kpolevec1(k))
               qkvec8 (k) = rpole(10, kpolevec1(k))
               qkvec9 (k) = rpole(13, kpolevec1(k))
            endif
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            if(k.le.nnelst1) then
               ukvecx (k) = uind ( 1, kpolevec1(k))
               ukvecy (k) = uind ( 2, kpolevec1(k))
               ukvecz (k) = uind ( 3, kpolevec1(k))

               ukpvecx(k) = uinp ( 1, kpolevec1(k))
               ukpvecy(k) = uinp ( 2, kpolevec1(k))
               ukpvecz(k) = uinp ( 3, kpolevec1(k))
            endif
         enddo

c
c     get reciprocal distance terms for this interaction
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            invr2vec(k) = (xpsvec(k)**2 + ypsvec(k)**2 + zpsvec(k)**2)
     &                   ** ( - one )
            rvec    (k) = invr2vec(k) ** ( - half )
            invrvec (k) = invr2vec(k) ** half
            ralphavec(k)  = aewald * rvec(k)               
         enddo
c
c     calculate the real space Ewald error function terms
c
c
         call vderfc(neloop16, ralphavec, bn0vec) 

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            exp2avec(k) = exp( - ralphavec(k)**2)
            bn0vec  (k) = bn0vec(k)   * invrvec          (k)
            bn1vec  (k) = (  1.0_ti_p    * bn0vec           (k)
     &                     +  alsq2   * alsq2n * exp2avec(k)
     &                    )           * invr2vec         (k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            bn2vec  (k) = (  3.0_ti_p    * bn1vec           (k)
     &                     + alsq2**2 * alsq2n * exp2avec(k)
     &                    )           * invr2vec         (k)
            bn3vec  (k) = (  5.0_ti_p    * bn2vec           (k)
     &                     + alsq2**3 * alsq2n * exp2avec(k)
     &                    )           * invr2vec         (k)
            bn4vec  (k) = (  7.0_ti_p    * bn3vec           (k)
     &                     + alsq2**4 * alsq2n * exp2avec(k)
     &                    )           * invr2vec         (k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            bn0vec(k) = f  * bn0vec(k)
            bn1vec(k) = f  * bn1vec(k)
            bn2vec(k) = f  * bn2vec(k)
            bn3vec(k) = f  * bn3vec(k)
            bn4vec(k) = f  * bn4vec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1,neloop8
            rr3vec(k) =     f * invrvec(k) * invr2vec(k)
            rr5vec(k) = 3.0_ti_p * rr3vec (k) * invr2vec(k)
            rr7vec(k) = 5.0_ti_p * rr5vec (k) * invr2vec(k)
            rr9vec(k) = 7.0_ti_p * rr7vec (k) * invr2vec(k)
         enddo
c
c     apply Thole polarization damping to scale factors
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            if(k.le.nnelst1) then
               tholevec (k) = thole        (kpolevec1(k))
               pgammavec(k) = min( pti,     tholevec (k))
               dampvec  (k) = pdi * pdamp  (kpolevec1(k))
            endif
            if (dampvec(k) .ne. 0.0_ti_p)
     &                          invdampvec(k) = dampvec(k) ** ( - one ) 
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1,neloop8
               dampvec1(k) = -  pgammavec(k)
     &                        * (rvec(k) * invdampvec(k)) ** 3
               expdampvec1(k) = exp(dampvec1(k))
               davec      (k) = dampvec1(k) * expdampvec1(k)
         enddo
c
c     intermediates involving Thole damping and scale factors
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            sc3vec (k) = 1.0_ti_p - expdampvec1(k)

            psc3vec(k) = 1.0_ti_p - sc3vec(k) * pscalevec(k)
            dsc3vec(k) = 1.0_ti_p - sc3vec(k) * dscalevec(k)
            usc3vec(k) = 1.0_ti_p - sc3vec(k) * uscalevec(k)

            psr3vec(k) = bn1vec(k) - psc3vec(k) * rr3vec(k)
            dsr3vec(k) = bn1vec(k) - dsc3vec(k) * rr3vec(k)
            usr3vec(k) = bn1vec(k) - usc3vec(k) * rr3vec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            sc5vec (k) = 1.0_ti_p - (1.0_ti_p - dampvec1(k))
     &                              *expdampvec1(k)

            psc5vec(k) = 1.0_ti_p - sc5vec(k) * pscalevec(k)
            dsc5vec(k) = 1.0_ti_p - sc5vec(k) * dscalevec(k)
            usc5vec(k) = 1.0_ti_p - sc5vec(k) * uscalevec(k)

            psr5vec(k) = bn2vec(k) - psc5vec(k) * rr5vec(k)
            dsr5vec(k) = bn2vec(k) - dsc5vec(k) * rr5vec(k)
            usr5vec(k) = bn2vec(k) - usc5vec(k) * rr5vec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            sc7vec (k) = 1.0_ti_p - (  1.0_ti_p - dampvec1(k)
     &                            + 0.6_ti_p * dampvec1(k)**2)
     &                          * expdampvec1(k)

            psc7vec(k) = 1.0_ti_p - sc7vec(k) * pscalevec(k)
            dsc7vec(k) = 1.0_ti_p - sc7vec(k) * dscalevec(k)

            psr7vec(k) = bn3vec(k) - psc7vec(k) * rr7vec(k)
            dsr7vec(k) = bn3vec(k) - dsc7vec(k) * rr7vec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1,neloop8
            rc3vecx (k) = - 3.0_ti_p * davec(k) *xpsvec(k) * invr2vec(k)
            prc3vecx(k) = rc3vecx(k) * pscalevec(k)
            drc3vecx(k) = rc3vecx(k) * dscalevec(k)
            urc3vecx(k) = rc3vecx(k) * uscalevec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1,neloop8
            rc3vecy (k) = - 3.0_ti_p * davec(k) *ypsvec(k) * invr2vec(k)
            prc3vecy(k) = rc3vecy(k) * pscalevec(k)
            drc3vecy(k) = rc3vecy(k) * dscalevec(k)
            urc3vecy(k) = rc3vecy(k) * uscalevec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1,neloop8
            rc3vecz (k) = - 3.0_ti_p  * davec(k)*zpsvec(k) * invr2vec(k)
            prc3vecz(k) =   rc3vecz(k) * pscalevec(k)
            drc3vecz(k) =   rc3vecz(k) * dscalevec(k)
            urc3vecz(k) =   rc3vecz(k) * uscalevec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            rc5vecx (k) = - dampvec1(k) * rc3vecx(k)
            prc5vecx(k) =    rc5vecx(k) * pscalevec(k)
            drc5vecx(k) =    rc5vecx(k) * dscalevec(k)
            urc5vecx(k) =    rc5vecx(k) * uscalevec(k)
!        enddo
!
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!        do k = 1,neloop8
            rc5vecy (k) = - dampvec1(k) * rc3vecy(k)
            prc5vecy(k) =    rc5vecy(k) * pscalevec(k)
            drc5vecy(k) =    rc5vecy(k) * dscalevec(k)
            urc5vecy(k) =    rc5vecy(k) * uscalevec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1,neloop8
            rc5vecz (k) = - dampvec1(k) * rc3vecz(k)
            prc5vecz(k) =    rc5vecz(k) * pscalevec(k)
            drc5vecz(k) =    rc5vecz(k) * dscalevec(k)
            urc5vecz(k) =    rc5vecz(k) * uscalevec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            rc7vecx (k) = - (0.2_ti_p + 0.6_ti_p*dampvec1(k))*rc5vecx(k)
            prc7vecx(k) = rc7vecx(k) * pscalevec(k)
            drc7vecx(k) = rc7vecx(k) * dscalevec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            rc7vecy (k) = - (0.2_ti_p + 0.6_ti_p*dampvec1(k))*rc5vecy(k)
            prc7vecy(k) = rc7vecy(k) * pscalevec(k)
            drc7vecy(k) = rc7vecy(k) * dscalevec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            rc7vecz (k) = - (0.2_ti_p + 0.6_ti_p*dampvec1(k))*rc5vecz(k)
            prc7vecz(k) = rc7vecz(k) * pscalevec(k)
            drc7vecz(k) = rc7vecz(k) * dscalevec(k)
         enddo

c
c     intermediates involving moments and distance separation
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8

            drivec(k)   =  divec (1) * xpsvec(k)
     &                   + divec (2) * ypsvec(k)
     &                   + divec (3) * zpsvec(k)
            drkvec(k)   =  dkvecx(k) * xpsvec(k)
     &                   + dkvecy(k) * ypsvec(k)
     &                   + dkvecz(k) * zpsvec(k)
            qrivecx(k)  =  qivec (1) * xpsvec(k)
     &                   + qivec (2) * ypsvec(k)
     &                   + qivec (3) * zpsvec(k)
            qrivecy(k)  =  qivec (4) * xpsvec(k)
     &                   + qivec (5) * ypsvec(k)
     &                   + qivec (6) * zpsvec(k)
            qrivecz(k)  =  qivec (7) * xpsvec(k)
     &                   + qivec (8) * ypsvec(k)
     &                   + qivec (9) * zpsvec(k)

         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            qrkvecx(k) =   qkvec1 (k) * xpsvec(k)
     &                   + qkvec2 (k) * ypsvec(k)
     &                   + qkvec3 (k) * zpsvec(k)
            qrkvecy(k) =   qkvec4 (k) * xpsvec(k)
     &                   + qkvec5 (k) * ypsvec(k)
     &                   + qkvec6 (k) * zpsvec(k)
            qrkvecz(k) =   qkvec7 (k) * xpsvec(k)
     &                   + qkvec8 (k) * ypsvec(k)
     &                   + qkvec9 (k) * zpsvec(k)
            qrrivec(k)  =  qrivecx(k) * xpsvec(k)
     &                   + qrivecy(k) * ypsvec(k)
     &                   + qrivecz(k) * zpsvec(k)
            qrrkvec(k)  =  qrkvecx(k) * xpsvec(k)
     &                   + qrkvecy(k) * ypsvec(k)
     &                   + qrkvecz(k) * zpsvec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            urivec(k)  =   uivec  (1) * xpsvec(k)
     &                   + uivec  (2) * ypsvec(k)
     &                   + uivec  (3) * zpsvec(k)
            urkvec(k)  =   ukvecx (k) * xpsvec(k)
     &                   + ukvecy (k) * ypsvec(k)
     &                   + ukvecz (k) * zpsvec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            uripvec(k) =   uipvec (1) * xpsvec(k)
     &                   + uipvec (2) * ypsvec(k)
     &                   + uipvec (3) * zpsvec(k)
            urkpvec(k) =   ukpvecx(k) * xpsvec(k)
     &                   + ukpvecy(k) * ypsvec(k)
     &                   + ukpvecz(k) * zpsvec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
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

c
c     calculate intermediate terms for polarization interaction
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            term1vec(k) =  ckvec (k) * urivec (k)
     &                   - ci        * urkvec (k) + duikvec(k)
            term2vec(k) =  2.0_ti_p     * quikvec(k)
     &                   - urivec(k) * drkvec (k)
     &                   - drivec(k) * urkvec (k)
            term3vec(k) =  urivec(k) * qrrkvec(k)
     &                   - urkvec(k) * qrrivec(k)
         enddo
c
c     compute the energy contribution for this interaction
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k = 1, neloop8
            if(k.le.nnelst1)  ep = ep + term1vec(k) * psr3vec(k)
     &                                + term2vec(k) * psr5vec(k)
     &                                + term3vec(k) * psr7vec(k)
         enddo
cDIR$ ASSUME (mod(neloop8,8).eq.0)
cDIR$ VECTOR ALIGNED
c        do k = 1, neloop8
c   
c     compute the full undamped energy for this interaction
c   
c           psr3vec(k) = rr3vec(k) * sc3vec(k) * pscalevec(k)
c           psr5vec(k) = rr5vec(k) * sc5vec(k) * pscalevec(k)
c           psr7vec(k) = rr7vec(k) * sc7vec(k) * pscalevec(k)
c           efullvec(k) =  term1vec(k) * psr3vec(k)
c    &                   + term2vec(k) * psr5vec(k)
c    &                   + term3vec(k) * psr7vec(k)
c        enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, neloop8
            qrimodvecx(k) =  qrivecy(k) * ypsvec(k)
     &                     + qrivecz(k) * zpsvec(k)
                      
            qrimodvecy(k) =  qrivecx(k) * xpsvec(k)
     &                     + qrivecz(k) * zpsvec(k)
                      
            qrimodvecz(k) =  qrivecx(k) * xpsvec(k)
     &                     + qrivecy(k) * ypsvec(k)
                      
            qrkmodvecx(k) =  qrkvecy(k) * ypsvec(k)
     &                     + qrkvecz(k) * zpsvec(k)
                      
            qrkmodvecy(k) =  qrkvecx(k) * xpsvec(k)
     &                     + qrkvecz(k) * zpsvec(k)
                      
            qrkmodvecz(k) =  qrkvecx(k) * xpsvec(k)
     &                     + qrkvecy(k) * ypsvec(k)
         enddo
c
c     get the dEd/dR terms used for direct polarization force
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm1(k) =  bn2vec(k) - dsc3vec(k) * rr5vec(k)
            dterm2(k) =  bn3vec(k) - dsc5vec(k) * rr7vec(k)
         enddo
c
c       Straight terms ( xx, yy ,zz )
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm3x(k) = - dsr3vec(k)
     &                   + dterm1 (k) * xpsvec(k)**2
     &                   - rr3vec (k) * xpsvec(k) * drc3vecx(k)
            dterm3y(k) = - dsr3vec(k)
     &                   + dterm1 (k) * ypsvec(k)**2
     &                   - rr3vec (k) * ypsvec(k) * drc3vecy(k)
            dterm3z(k) = - dsr3vec(k)
     &                   + dterm1 (k) * zpsvec(k)**2
     &                   - rr3vec (k) * zpsvec(k) * drc3vecz(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm4x(k) =  rr3vec(k)  * drc3vecx(k)
     &                  - (dterm1(k) + dsr5vec(k)) * xpsvec(k)
            dterm4y(k) =  rr3vec(k)  * drc3vecy(k)
     &                  - (dterm1(k) + dsr5vec(k)) * ypsvec(k)
            dterm4z(k) =  rr3vec(k)  * drc3vecz(k)
     &                  - (dterm1(k) + dsr5vec(k)) * zpsvec(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm5x(k) = - dsr5vec(k) + dterm2(k) * xpsvec  (k) ** 2
     &                   - rr5vec (k) * xpsvec(k) * drc5vecx(k)    
            dterm5y(k) = - dsr5vec(k) + dterm2(k) * ypsvec  (k) ** 2
     &                   - rr5vec (k) * ypsvec(k) * drc5vecy(k)
            dterm5z(k) = - dsr5vec(k) + dterm2(k) * zpsvec  (k) ** 2
     &                   - rr5vec (k) * zpsvec(k) * drc5vecz(k)    
         enddo                                                     

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm6x(k) =  (bn4vec(k) - dsc7vec(k) * rr9vec  (k))
     &                                            * xpsvec  (k) ** 2
     &                  -  bn3vec(k)                               
     &                  -  rr7vec(k) * xpsvec (k) * drc7vecx(k)    
                   
            dterm6y(k) =  (bn4vec(k) - dsc7vec(k) * rr9vec  (k))
     &                                            * ypsvec  (k) ** 2
     &                  -  bn3vec(k)                               
     &                  -  rr7vec(k) * ypsvec (k) * drc7vecy(k)    
                   
            dterm6z(k) =  (bn4vec(k) - dsc7vec(k) * rr9vec  (k))
     &                                            * zpsvec  (k) ** 2
     &                  -  bn3vec(k)                               
     &                  -  rr7vec(k) * zpsvec (k) * drc7vecz(k)    
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k = 1, neloop8
            dterm7x(k) = rr5vec   (k)         * drc5vecx(k)
     &                  - 2.0_ti_p * bn3vec(k)   * xpsvec  (k)
     &                  + (dsc5vec(k) + 1.5_ti_p * dsc7vec (k))
     &                  * rr7vec  (k)         * xpsvec  (k)
            dterm7y(k) = rr5vec   (k)         * drc5vecy(k)
     &                  - 2.0_ti_p * bn3vec(k)   * ypsvec  (k)
     &                  + (dsc5vec(k) + 1.5_ti_p * dsc7vec(k))
     &                  * rr7vec  (k)         * ypsvec  (k)
            dterm7z(k) = rr5vec   (k)         * drc5vecz(k)
     &                  - 2.0_ti_p * bn3vec(k)   * zpsvec  (k)
     &                  + (dsc5vec(k) + 1.5_ti_p * dsc7vec (k))
     &                  * rr7vec  (k)         * zpsvec  (k)
         enddo
c
c       Straight terms ( xx, yy ,zz )
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tisvecx(k) =          ci            * dterm3x(k)
     &                  +         divec    (1)  * dterm4x(k)
     &                  +         drivec   (k)  * dterm5x(k)
     &                  + 2.0_ti_p * dsr5vec  (k)  * qivec  (1)
     &                  +         qrimodvecx(k) * dsc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrivecx  (k)  * dterm7x(k)
     &                  +         qrrivec  (k)  * dterm6x(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tisvecy(k) =          ci            * dterm3y(k)
     &                  +         divec    (2)  * dterm4y(k)
     &                  +         drivec   (k)  * dterm5y(k)
     &                  + 2.0_ti_p * dsr5vec  (k)  * qivec  (5)
     &                  +         qrimodvecy(k) * dsc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrivecy  (k)  * dterm7y(k)
     &                  +         qrrivec  (k)  * dterm6y(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tisvecz(k) =          ci            * dterm3z(k)
     &                  +         divec    (3)  * dterm4z(k)
     &                  +         drivec   (k)  * dterm5z(k)
     &                  + 2.0_ti_p * dsr5vec  (k)  * qivec  (9)
     &                  +         qrimodvecz(k) * dsc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrivecz  (k)  * dterm7z(k)
     &                  +         qrrivec  (k)  * dterm6z(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tksvecx(k) =          ckvec     (k) * dterm3x(k)
     &                  -         dkvecx    (k) * dterm4x(k)
     &                  -         drkvec    (k) * dterm5x(k)
     &                  + 2.0_ti_p * dsr5vec   (k) * qkvec1 (k)
     &                  +         qrkmodvecx(k) * dsc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrkvecx   (k) * dterm7x(k)
     &                  +         qrrkvec   (k) * dterm6x(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tksvecy(k) =          ckvec     (k) * dterm3y(k)
     &                  -         dkvecy    (k) * dterm4y(k)
     &                  -         drkvec    (k) * dterm5y(k)
     &                  + 2.0_ti_p * dsr5vec   (k) * qkvec5 (k)
     &                  +         qrkmodvecy(k) * dsc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrkvecy   (k) * dterm7y(k)
     &                  +         qrrkvec   (k) * dterm6y(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tksvecz(k) =          ckvec     (k) * dterm3z(k)
     &                  -         dkvecz    (k) * dterm4z(k)
     &                  -         drkvec    (k) * dterm5z(k)
     &                  + 2.0_ti_p * dsr5vec   (k) * qkvec9 (k)
     &                  +         qrkmodvecz(k) * dsc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrkvecz   (k) * dterm7z(k)
     &                  +         qrrkvec   (k) * dterm6z(k)
         enddo                             
c
c       Crossed terms ( xy = yx , xz = zx ,yz = zy )
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tmp1vecx(k) = xpsvec(k)    * ypsvec(k)
            tmp1vecy(k) = xpsvec(k)    * zpsvec(k)
            tmp1vecz(k) = ypsvec(k)    * zpsvec(k)

            tmp2vecx(k) = drc3vecx(k)
            tmp2vecy(k) = drc3vecx(k)
            tmp2vecz(k) = drc3vecy(k)

            tmp3vecx(k) = ypsvec(k) * drc3vecx(k)
            tmp3vecy(k) = zpsvec(k) * drc3vecx(k)
            tmp3vecz(k) = zpsvec(k) * drc3vecy(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tmp4vecx(k) = xpsvec(k)
            tmp4vecy(k) = xpsvec(k)
            tmp4vecz(k) = ypsvec(k)

            tmp5vecx(k) = ypsvec(k) * drc5vecx(k)
            tmp5vecy(k) = zpsvec(k) * drc5vecx(k)
            tmp5vecz(k) = zpsvec(k) * drc5vecy(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tmp6vecx(k) = ypsvec(k) * drc7vecx(k)
            tmp6vecy(k) = zpsvec(k) * drc7vecx(k)
            tmp6vecz(k) = zpsvec(k) * drc7vecy(k)

            tmp7vecx(k) =             drc5vecx(k)
            tmp7vecy(k) =             drc5vecx(k)
            tmp7vecz(k) =             drc5vecy(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm3x(k) =  dterm1(k) * tmp1vecx(k)
     &                  - rr3vec(k) * tmp3vecx(k)
            dterm3y(k) =  dterm1(k) * tmp1vecy(k) 
     &                  - rr3vec(k) * tmp3vecy(k)
            dterm3z(k) =  dterm1(k) * tmp1vecz(k) 
     &                  - rr3vec(k) * tmp3vecz(k)
         enddo                                    
                                                  
!DIR$ ASSUME (mod(neloop8,8).eq.0)     
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8                     
            dterm4x(k) =  rr3vec(k) * tmp2vecx(k)
     &                  - dterm1(k) * tmp4vecx(k)
            dterm4y(k) =  rr3vec(k) * tmp2vecy(k)
     &                  - dterm1(k) * tmp4vecy(k)
            dterm4z(k) =  rr3vec(k) * tmp2vecz(k)
     &                  - dterm1(k) * tmp4vecz(k)

            dterm5x(k) =  dterm2(k) * tmp1vecx(k) 
     &                 -  rr5vec(k) * tmp5vecx(k)
            dterm5y(k) =  dterm2(k) * tmp1vecy(k) 
     &                 -  rr5vec(k) * tmp5vecy(k)
            dterm5z(k) =  dterm2(k) * tmp1vecz(k) 
     &                 -  rr5vec(k) * tmp5vecz(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            dterm6x(k) =  (bn4vec  (k) - dsc7vec(k) * rr9vec(k))
     &                   * tmp1vecx(k)
     &                  -  rr7vec  (k) * tmp6vecx(k)
            dterm6y(k) =  (bn4vec  (k) - dsc7vec(k) * rr9vec(k))
     &                   * tmp1vecy(k)
     &                  -  rr7vec  (k) * tmp6vecy(k)
            dterm6z(k) =  (bn4vec  (k) - dsc7vec(k) * rr9vec(k))
     &                   * tmp1vecz(k)
     &                  -  rr7vec  (k) * tmp6vecz(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8

            dterm7x(k) =  rr5vec(k)  * tmp7vecx(k)
     &                  - dterm2(k)  * tmp4vecx(k)
            dterm7y(k) =  rr5vec(k)  * tmp7vecy(k)
     &                  - dterm2(k)  * tmp4vecy(k)
            dterm7z(k) =  rr5vec(k)  * tmp7vecz(k)
     &                  - dterm2(k)  * tmp4vecz(k)
         enddo

c
c       Crossed terms ( xy = yx , xz = zx ,yz = zy )
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            ticvecx(k) =          ci         * dterm3x(k)
     &                  -         dsr5vec(k) * divec  (1)
     &                                       * ypsvec (k)
     &                  +         divec  (2) * dterm4x(k)
     &                  +         drivec (k) * dterm5x(k)
     &                  + 2.0_ti_p * dsr5vec(k) * qivec  (4)
     &                  - 2.0_ti_p * dsr7vec(k) * ypsvec (k)
     &                                       * qrivecx(k)
     &                  + 2.0_ti_p * qrivecy(k) * dterm7x(k)
     &                  +         qrrivec(k) * dterm6x(k)

            ticvecy(k) =          ci         * dterm3y(k)
     &                  -         dsr5vec(k) * divec  (1)
     &                                       * zpsvec (k)
     &                  +         divec  (3) * dterm4y(k)
     &                  +         drivec (k) * dterm5y(k)
     &                  + 2.0_ti_p * dsr5vec(k) * qivec  (7)
     &                  - 2.0_ti_p * dsr7vec(k) * zpsvec (k)
     &                                       * qrivecx(k)
     &                  + 2.0_ti_p * qrivecz(k) * dterm7y(k)
     &                  +         qrrivec(k) * dterm6y(k)

            ticvecz(k) =          ci         * dterm3z(k)
     &                  -         dsr5vec(k) * divec  (2)
     &                                       * zpsvec (k)
     &                  +         divec  (3) * dterm4z(k)
     &                  +         drivec (k) * dterm5z(k)
     &                  + 2.0_ti_p * dsr5vec(k) * qivec  (8)
     &                  - 2.0_ti_p * dsr7vec(k) * zpsvec (k)
     &                                       * qrivecy(k)
     &                  + 2.0_ti_p * qrivecz(k) * dterm7z(k)
     &                  +         qrrivec(k) * dterm6z(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tkcvecx(k) =          ckvec  (k) * dterm3x(k)
     &                  +         dsr5vec(k) * dkvecx (k)
     &                                       * ypsvec (k)
     &                  -         dkvecy (k) * dterm4x(k)
     &                  -         drkvec (k) * dterm5x(k)
     &                  + 2.0_ti_p * dsr5vec(k) * qkvec4 (k)
     &                  - 2.0_ti_p * dsr7vec(k) * ypsvec (k)
     &                                       * qrkvecx(k)
     &                  + 2.0_ti_p * qrkvecy(k) * dterm7x(k)
     &                  +         qrrkvec(k) * dterm6x(k)

            tkcvecy(k) =          ckvec  (k) * dterm3y(k)
     &                  +         dsr5vec(k) * dkvecx (k)
     &                                       * zpsvec (k)
     &                  -         dkvecz (k) * dterm4y(k)
     &                  -         drkvec (k) * dterm5y(k)
     &                  + 2.0_ti_p * dsr5vec(k) * qkvec7 (k)
     &                  - 2.0_ti_p * dsr7vec(k) * zpsvec (k)
     &                                       * qrkvecx(k)
     &                  + 2.0_ti_p * qrkvecz(k) * dterm7y(k)
     &                  +         qrrkvec(k) * dterm6y(k)

            tkcvecz(k) =          ckvec  (k) * dterm3z(k)
     &                  +         dsr5vec(k) * dkvecy (k)
     &                                       * zpsvec (k)
     &                  -         dkvecz (k) * dterm4z(k)
     &                  -         drkvec (k) * dterm5z(k)
     &                  + 2.0_ti_p * dsr5vec(k) * qkvec8 (k)
     &                  - 2.0_ti_p * dsr7vec(k) * zpsvec (k)
     &                                       * qrkvecy(k)
     &                  + 2.0_ti_p * qrkvecz(k) * dterm7z(k)
     &                  +         qrrkvec(k) * dterm6z(k)
         enddo
c
c         Construct matrixes for dot_product
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tivec1(k) = tisvecx(k) !xx   
            tivec5(k) = tisvecy(k) !yy   
            tivec9(k) = tisvecz(k) !zz   
            tivec2(k) = ticvecx(k) !xy   
            tivec4(k) = ticvecx(k) !yx   
            tivec3(k) = ticvecy(k) !xz   
            tivec7(k) = ticvecy(k) !zx   
            tivec6(k) = ticvecz(k) !yz   
            tivec8(k) = ticvecz(k) !zy   
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tkvec1(k) = tksvecx(k) !xx   
            tkvec5(k) = tksvecy(k) !yy   
            tkvec9(k) = tksvecz(k) !zz   
            tkvec2(k) = tkcvecx(k) !xy   
            tkvec4(k) = tkcvecx(k) !yx   
            tkvec3(k) = tkcvecy(k) !xz   
            tkvec7(k) = tkcvecy(k) !zx   
            tkvec6(k) = tkcvecz(k) !yz   
            tkvec8(k) = tkcvecz(k) !zy   
         enddo
c
c        do dot product
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            frcxvec(k) =  tivec1(k) * ukpvecx(k)
     &                  + tivec2(k) * ukpvecy(k)
     &                  + tivec3(k) * ukpvecz(k)
     &                  - tkvec1(k) * uipvec (1)
     &                  - tkvec2(k) * uipvec (2)
     &                  - tkvec3(k) * uipvec (3)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            frcyvec(k) =  tivec4(k) * ukpvecx(k)
     &                  + tivec5(k) * ukpvecy(k)
     &                  + tivec6(k) * ukpvecz(k)
     &                  - tkvec4(k) * uipvec (1)
     &                  - tkvec5(k) * uipvec (2)
     &                  - tkvec6(k) * uipvec (3)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            frczvec(k) =  tivec7(k) * ukpvecx(k)
     &                  + tivec8(k) * ukpvecy(k)
     &                  + tivec9(k) * ukpvecz(k)
     &                  - tkvec7(k) * uipvec (1)
     &                  - tkvec8(k) * uipvec (2)
     &                  - tkvec9(k) * uipvec (3)
         enddo
c
c     get the dEp/dR terms used for direct polarization force
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm1(k) =  bn2vec(k) - psc3vec(k) * rr5vec (k)
            dterm2(k) =  bn3vec(k) - psc5vec(k) * rr7vec (k)
c
c       Straight terms ( xx, yy ,zz )
c
            dterm3x(k) = - psr3vec(k)
     &                   + dterm1 (k) * xpsvec(k) ** 2
     &                   - rr3vec (k) * xpsvec(k) * prc3vecx(k)
                   
            dterm3y(k) = - psr3vec(k)
     &                   + dterm1 (k) * ypsvec(k) ** 2 
     &                   - rr3vec (k) * ypsvec(k) * prc3vecy(k)
                   
            dterm3z(k) = - psr3vec(k)
     &                   + dterm1 (k) * zpsvec(k) ** 2
     &                   - rr3vec (k) * zpsvec(k) * prc3vecz(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm4x(k) =   rr3vec (k) * prc3vecx(k)
     &                   - dterm1 (k) * xpsvec  (k)
     &                   - psr5vec(k) * xpsvec  (k)
            dterm4y(k) =   rr3vec (k) * prc3vecy(k)
     &                   - dterm1 (k) * ypsvec  (k)
     &                   - psr5vec(k) * ypsvec  (k)
            dterm4z(k) =   rr3vec (k) * prc3vecz(k)
     &                   - dterm1 (k) * zpsvec  (k)
     &                   - psr5vec(k) * zpsvec  (k)

            dterm5x(k) = - psr5vec(k) + dterm2(k) * xpsvec  (k) ** 2
     &                   - rr5vec( k) * xpsvec(k) * prc5vecx(k)
            dterm5y(k) = - psr5vec(k) + dterm2(k) * ypsvec  (k) ** 2
     &                   - rr5vec (k) * ypsvec(k) * prc5vecy(k)
            dterm5z(k) = - psr5vec(k) + dterm2(k) * zpsvec  (k) ** 2
     &                   - rr5vec (k) * zpsvec(k) * prc5vecz(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm6x(k) =  (bn4vec(k) - psc7vec(k) * rr9vec  (k))
     &                                            * xpsvec  (k) ** 2
     &                  - bn3vec (k)
     &                  - rr7vec (k) * xpsvec (k) * prc7vecx(k)

            dterm6y(k) =  (bn4vec(k) - psc7vec(k) * rr9vec  (k))
     &                                            * ypsvec  (k) ** 2
     &                  - bn3vec (k)
     &                  - rr7vec (k) * ypsvec (k) * prc7vecy(k)
                   
            dterm6z(k) =  (bn4vec(k) - psc7vec(k) * rr9vec  (k))
     &                                            * zpsvec  (k) ** 2
     &                  - bn3vec (k)
     &                  - rr7vec (k) * zpsvec (k) * prc7vecz(k)

            dterm7x(k) =          rr5vec(k)    * prc5vecx(k)
     &                  - 2.0_ti_p * bn3vec(k)    * xpsvec  (k)
     &                  +  (psc5vec(k) + 1.5_ti_p * psc7vec (k) )
     &                          * rr7vec(k)    * xpsvec  (k)
            dterm7y(k) =          rr5vec(k)    * prc5vecy(k)
     &                  - 2.0_ti_p * bn3vec(k)    * ypsvec  (k)
     &                  +  (psc5vec(k) + 1.5_ti_p * psc7vec (k) )
     &                          * rr7vec(k)    * ypsvec  (k)
            dterm7z(k) =          rr5vec(k)    * prc5vecz(k)
     &                  - 2.0_ti_p * bn3vec(k)    * zpsvec  (k)
     &                  +  (psc5vec(k) + 1.5_ti_p * psc7vec (k) )
     &                          * rr7vec(k)    * zpsvec  (k)
         enddo
c
c       Straight terms ( xx, yy ,zz )
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tisvecx(k) =          ci            * dterm3x(k)
     &                  +         divec     (1) * dterm4x(k)
     &                  +         drivec    (k) * dterm5x(k)
     &                  + 2.0_ti_p * psr5vec   (k) * qivec  (1)
     &                  +         qrimodvecx(k) * psc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrivecx   (k) * dterm7x(k)
     &                  +         qrrivec   (k) * dterm6x(k)

            tisvecy(k) =          ci            * dterm3y(k)
     &                  +         divec     (2) * dterm4y(k)
     &                  +         drivec    (k) * dterm5y(k)
     &                  + 2.0_ti_p * psr5vec   (k) * qivec  (5)
     &                  +         qrimodvecy(k) * psc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrivecy   (k) * dterm7y(k)
     &                  +         qrrivec   (k) * dterm6y(k)

            tisvecz(k) =          ci            * dterm3z(k)
     &                  +         divec     (3) * dterm4z(k)
     &                  +         drivec    (k) * dterm5z(k)
     &                  + 2.0_ti_p * psr5vec   (k) * qivec  (9)
     &                  +         qrimodvecz(k) * psc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrivecz   (k) * dterm7z(k)
     &                  +         qrrivec   (k) * dterm6z(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tksvecx(k) =          ckvec     (k) * dterm3x(k)
     &                  -         dkvecx    (k) * dterm4x(k)
     &                  -         drkvec    (k) * dterm5x(k)
     &                  + 2.0_ti_p * psr5vec   (k) * qkvec1 (k)
     &                  +         qrkmodvecx(k) * psc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrkvecx   (k) * dterm7x(k)
     &                  +         qrrkvec   (k) * dterm6x(k)

            tksvecy(k) =          ckvec     (k) * dterm3y(k)
     &                  -         dkvecy    (k) * dterm4y(k)
     &                  -         drkvec    (k) * dterm5y(k)
     &                  + 2.0_ti_p * psr5vec   (k) * qkvec5 (k)
     &                  +         qrkmodvecy(k) * psc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrkvecy   (k) * dterm7y(k)
     &                  +         qrrkvec   (k) * dterm6y(k)

            tksvecz(k) =          ckvec     (k) * dterm3z(k)
     &                  -         dkvecz    (k) * dterm4z(k)
     &                  -         drkvec    (k) * dterm5z(k)
     &                  + 2.0_ti_p * psr5vec   (k) * qkvec9 (k)
     &                  +         qrkmodvecz(k) * psc7vec(k)
     &                                          * rr7vec (k)
     &                  + 2.0_ti_p * qrkvecz   (k) * dterm7z(k)
     &                  +         qrrkvec   (k) * dterm6z(k)
         enddo

c
c       Crossed terms ( xy = yx , xz = zx ,yz = zy )
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tmp1vecx(k) = xpsvec(k) * ypsvec(k)
            tmp1vecy(k) = xpsvec(k) * zpsvec(k)
            tmp1vecz(k) = ypsvec(k) * zpsvec(k)

            tmp2vecx(k) = prc3vecx(k)
            tmp2vecy(k) = prc3vecx(k)
            tmp2vecz(k) = prc3vecy(k)

            tmp3vecx(k) = ypsvec(k) * prc3vecx(k)
            tmp3vecy(k) = zpsvec(k) * prc3vecx(k)
            tmp3vecz(k) = zpsvec(k) * prc3vecy(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tmp4vecx(k) = xpsvec(k)
            tmp4vecy(k) = xpsvec(k)
            tmp4vecz(k) = ypsvec(k)

            tmp5vecx(k) = ypsvec(k) * prc5vecx(k)
            tmp5vecy(k) = zpsvec(k) * prc5vecx(k)
            tmp5vecz(k) = zpsvec(k) * prc5vecy(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            tmp6vecx(k) = ypsvec(k) * prc7vecx(k)
            tmp6vecy(k) = zpsvec(k) * prc7vecx(k)
            tmp6vecz(k) = zpsvec(k) * prc7vecy(k)

            tmp7vecx(k) = prc5vecx(k)
            tmp7vecy(k) = prc5vecx(k)
            tmp7vecz(k) = prc5vecy(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
         do k = 1, neloop8
            dterm3x(k) =  dterm1(k)  * tmp1vecx(k)
     &                  - rr3vec(k)  * tmp3vecx(k)
            dterm3y(k) =  dterm1(k)  * tmp1vecy(k) 
     &                  - rr3vec(k)  * tmp3vecy(k)
            dterm3z(k) =  dterm1(k)  * tmp1vecz(k) 
     &                  - rr3vec(k)  * tmp3vecz(k)
                                                   
            dterm4x(k) =  rr3vec(k)  * tmp2vecx(k)
     &                  - dterm1(k)  * tmp4vecx(k)
            dterm4y(k) =  rr3vec(k)  * tmp2vecy(k)
     &                  - dterm1(k)  * tmp4vecy(k)
            dterm4z(k) =  rr3vec(k)  * tmp2vecz(k)
     &                  - dterm1(k)  * tmp4vecz(k)
         enddo                                     
                                                   
!DIR$ ASSUME (mod(neloop8,8).eq.0)      
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8                      
            dterm5x(k) =  dterm2(k)  * tmp1vecx(k) 
     &                  - rr5vec(k)  * tmp5vecx(k)
            dterm5y(k) =  dterm2(k)  * tmp1vecy(k) 
     &                  - rr5vec(k)  * tmp5vecy(k)
            dterm5z(k) =  dterm2(k)  * tmp1vecz(k) 
     &                  - rr5vec(k)  * tmp5vecz(k)

            dterm6x(k) =  (bn4vec(k) - psc7vec (k) * rr9vec(k))
     &                               * tmp1vecx(k)
     &                  -  rr7vec(k) * tmp6vecx(k)
            dterm6y(k) =  (bn4vec(k) - psc7vec (k) * rr9vec(k))
     &                               * tmp1vecy(k)
     &                  -  rr7vec(k) * tmp6vecy(k)
            dterm6z(k) =  (bn4vec(k) - psc7vec (k) * rr9vec(k))
     &                               * tmp1vecz(k)
     &                  -  rr7vec(k) * tmp6vecz(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm7x(k) =  rr5vec(k)  * tmp7vecx(k)
     &                  - dterm2(k)  * tmp4vecx(k)
            dterm7y(k) =  rr5vec(k)  * tmp7vecy(k)
     &                  - dterm2(k)  * tmp4vecy(k)
            dterm7z(k) =  rr5vec(k)  * tmp7vecz(k)
     &                  - dterm2(k)  * tmp4vecz(k)
         enddo

c
c       Crossed terms ( xy = yx , xz = zx ,yz = zy )
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            ticvecx(k) =          ci         * dterm3x(k)
     &                  -         psr5vec(k) * divec  (1)
     &                                       * ypsvec (k)
     &                  +         divec  (2) * dterm4x(k)
     &                  +         drivec (k) * dterm5x(k)
     &                  + 2.0_ti_p * psr5vec(k) * qivec  (4)
     &                  - 2.0_ti_p * psr7vec(k) * ypsvec (k)
     &                                       * qrivecx(k)
     &                  + 2.0_ti_p * qrivecy(k) * dterm7x(k)
     &                  +         qrrivec(k) * dterm6x(k)
                   
            ticvecy(k) =          ci         * dterm3y(k)
     &                  -         psr5vec(k) * divec  (1)
     &                                       * zpsvec (k)
     &                  +         divec  (3) * dterm4y(k)
     &                  +         drivec (k) * dterm5y(k)
     &                  + 2.0_ti_p * psr5vec(k) * qivec  (7)
     &                  - 2.0_ti_p * psr7vec(k) * zpsvec (k)
     &                                       * qrivecx(k)
     &                  + 2.0_ti_p * qrivecz(k) * dterm7y(k)
     &                  +         qrrivec(k) * dterm6y(k)

            ticvecz(k) =          ci         * dterm3z(k)
     &                  -         psr5vec(k) * divec  (2)
     &                                       * zpsvec (k)
     &                  +         divec  (3) * dterm4z(k)
     &                  +         drivec (k) * dterm5z(k)
     &                  + 2.0_ti_p * psr5vec(k) * qivec  (8)
     &                  - 2.0_ti_p * psr7vec(k) * zpsvec (k)
     &                                       * qrivecy(k)
     &                  + 2.0_ti_p * qrivecz(k) * dterm7z(k)
     &                  +         qrrivec(k) * dterm6z(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tkcvecx(k) =          ckvec  (k) * dterm3x(k)
     &                  +         psr5vec(k) * dkvecx (k)
     &                                       * ypsvec (k)
     &                  -         dkvecy (k) * dterm4x(k)
     &                  -         drkvec (k) * dterm5x(k)
     &                  + 2.0_ti_p * psr5vec(k) * qkvec4 (k)
     &                  - 2.0_ti_p * psr7vec(k) * ypsvec (k)
     &                                       * qrkvecx(k)
     &                  + 2.0_ti_p * qrkvecy(k) * dterm7x(k)
     &                  +         qrrkvec(k) * dterm6x(k)

            tkcvecy(k) =          ckvec  (k) * dterm3y(k)
     &                  +         psr5vec(k) * dkvecx (k)
     &                                       * zpsvec (k)
     &                  -         dkvecz (k) * dterm4y(k)
     &                  -         drkvec (k) * dterm5y(k)
     &                  + 2.0_ti_p * psr5vec(k) * qkvec7 (k)
     &                  - 2.0_ti_p * psr7vec(k) * zpsvec (k)
     &                                       * qrkvecx(k)
     &                  + 2.0_ti_p * qrkvecz(k) * dterm7y(k)
     &                  +         qrrkvec(k) * dterm6y(k)

            tkcvecz(k) =          ckvec  (k) * dterm3z(k)
     &                  +         psr5vec(k) * dkvecy (k)
     &                                       * zpsvec (k)
     &                  -         dkvecz (k) * dterm4z(k)
     &                  -         drkvec (k) * dterm5z(k)
     &                  + 2.0_ti_p * psr5vec(k) * qkvec8 (k)
     &                  - 2.0_ti_p * psr7vec(k) * zpsvec (k)
     &                                       * qrkvecy(k)
     &                  + 2.0_ti_p * qrkvecz(k) * dterm7z(k)
     &                  +         qrrkvec(k) * dterm6z(k)
         enddo
c
c         Construct matrixes for dot_product
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tivec1(k) = tisvecx(k) !xx   
            tivec5(k) = tisvecy(k) !yy   
            tivec9(k) = tisvecz(k) !zz   
            tivec2(k) = ticvecx(k) !xy   
            tivec4(k) = ticvecx(k) !yx   
            tivec3(k) = ticvecy(k) !xz   
            tivec7(k) = ticvecy(k) !zx   
            tivec6(k) = ticvecz(k) !yz   
            tivec8(k) = ticvecz(k) !zy   
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tkvec1(k) = tksvecx(k) !xx   
            tkvec5(k) = tksvecy(k) !yy   
            tkvec9(k) = tksvecz(k) !zz   
            tkvec2(k) = tkcvecx(k) !xy   
            tkvec4(k) = tkcvecx(k) !yx   
            tkvec3(k) = tkcvecy(k) !xz   
            tkvec7(k) = tkcvecy(k) !zx   
            tkvec6(k) = tkcvecz(k) !yz   
            tkvec8(k) = tkcvecz(k) !zy   
         enddo
c
c        Do dot product
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8

            depxvec(k) =  tivec1(k) * ukvecx(k)
     &                  + tivec2(k) * ukvecy(k)
     &                  + tivec3(k) * ukvecz(k)
     &                  - tkvec1(k) * uivec (1)
     &                  - tkvec2(k) * uivec (2)
     &                  - tkvec3(k) * uivec (3)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            depyvec(k) =  tivec4(k) * ukvecx(k)
     &                  + tivec5(k) * ukvecy(k)
     &                  + tivec6(k) * ukvecz(k)
     &                  - tkvec4(k) * uivec (1)
     &                  - tkvec5(k) * uivec (2)
     &                  - tkvec6(k) * uivec (3)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            depzvec(k) =  tivec7(k) * ukvecx(k)
     &                  + tivec8(k) * ukvecy(k)
     &                  + tivec9(k) * ukvecz(k)
     &                  - tkvec7(k) * uivec (1)
     &                  - tkvec8(k) * uivec (2)
     &                  - tkvec9(k) * uivec (3)
!        enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!        do k = 1, neloop8

            frcxvec(k) = frcxvec(k) + depxvec(k)
            frcyvec(k) = frcyvec(k) + depyvec(k)
            frczvec(k) = frczvec(k) + depzvec(k)

         enddo

c
c     get the dtau/dr terms used for mutual polarization force
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm1(k) = bn2vec(k) - usc3vec(k) * rr5vec(k)
            dterm2(k) = bn3vec(k) - usc5vec(k) * rr7vec(k)
!        enddo

c
c       Straight terms ( xx, yy ,zz )
c

            dterm3x(k) = usr5vec(k) + dterm1   (k)
            dterm3y(k) = usr5vec(k) + dterm1   (k)
            dterm3z(k) = usr5vec(k) + dterm1   (k)
            dterm4x(k) = rr3vec (k) * uscalevec(k)
            dterm4y(k) = rr3vec (k) * uscalevec(k)
            dterm4z(k) = rr3vec (k) * uscalevec(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm5x(k) = - xpsvec(k)   * dterm3x(k)
     &                   + rc3vecx(k)  * dterm4x(k)
            dterm5y(k) = - ypsvec(k)   * dterm3y(k)
     &                   + rc3vecy(k)  * dterm4y(k)
            dterm5z(k) = - zpsvec(k)   * dterm3z(k)
     &                   + rc3vecz(k)  * dterm4z(k)

            dterm6x(k) = - usr5vec(k) + xpsvec(k) ** 2 * dterm2  (k)
     &                   - rr5vec (k) * xpsvec(k)      * urc5vecx(k)
            dterm6y(k) = - usr5vec(k) + ypsvec(k) ** 2 * dterm2  (k)
     &                   - rr5vec (k) * ypsvec(k)      * urc5vecy(k)
            dterm6z(k) = - usr5vec(k) + zpsvec(k) ** 2 * dterm2  (k)
     &                   - rr5vec (k) * zpsvec(k)      * urc5vecz(k)
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tisvecx(k) =  uivec (1) * dterm5x(k)
     &                  + urivec(k) * dterm6x(k)
            tisvecy(k) =  uivec (2) * dterm5y(k)
     &                  + urivec(k) * dterm6y(k)
            tisvecz(k) =  uivec (3) * dterm5z(k)
     &                  + urivec(k) * dterm6z(k)
            tksvecx(k) =  ukvecx(k) * dterm5x(k)
     &                  + urkvec(k) * dterm6x(k)
            tksvecy(k) =  ukvecy(k) * dterm5y(k)
     &                  + urkvec(k) * dterm6y(k)
            tksvecz(k) =  ukvecz(k) * dterm5z(k)
     &                  + urkvec(k) * dterm6z(k)
         enddo

c
c       Crossed terms ( xy = yx , xz = zx ,yz = zy )
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tmp1vecx(k)  = xpsvec  (k)             
            tmp1vecy(k)  = xpsvec  (k)             
            tmp1vecz(k)  = ypsvec  (k)             

            tmp2vecx(k)  = urc3vecx(k)             
            tmp2vecy(k)  = urc3vecx(k)             
            tmp2vecz(k)  = urc3vecy(k)             

            tmp4vecx(k)  = ypsvec  (k)             
            tmp4vecy(k)  = zpsvec  (k)             
            tmp4vecz(k)  = zpsvec  (k)             

            tmp5vecx(k)  = urc5vecx(k)             
            tmp5vecy(k)  = urc5vecx(k)             
            tmp5vecz(k)  = urc5vecy(k)             

            tmp6vecx(k)  = xpsvec  (k) * dterm2(k)
            tmp6vecy(k)  = xpsvec  (k) * dterm2(k)
            tmp6vecz(k)  = ypsvec  (k) * dterm2(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            dterm4x(k) = - usr5vec (k) * tmp4vecx(k)
            dterm4y(k) = - usr5vec (k) * tmp4vecy(k)
            dterm4z(k) = - usr5vec (k) * tmp4vecz(k)
            dterm5x(k) = - tmp1vecx(k) * dterm1  (k)   
     &                   +  rr3vec (k) * tmp2vecx(k)
            dterm5y(k) = - tmp1vecy(k) * dterm1  (k)   
     &                   +  rr3vec (k) * tmp2vecy(k)
            dterm5z(k) = - tmp1vecz(k) * dterm1  (k)   
     &                   +  rr3vec (k) * tmp2vecz(k)
         enddo                                            
                                                          
!DIR$ ASSUME (mod(neloop8,8).eq.0)             
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8                             
            dterm6x(k) =   tmp6vecx(k) * tmp4vecx(k)
     &                   - rr5vec  (k) * tmp4vecx(k)
     &                                 * tmp5vecx(k)
            dterm6y(k) =   tmp6vecy(k) * tmp4vecy(k)
     &                   - rr5vec  (k) * tmp4vecy(k)
     &                                 * tmp5vecy(k)
            dterm6z(k) =   tmp6vecz(k) * tmp4vecz(k)
     &                   - rr5vec  (k) * tmp4vecz(k)
     &                                 * tmp5vecz(k)
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            ticvecx(k) =  uivec (1) * dterm4x(k)
     &                  + uivec (2) * dterm5x(k)
     &                  + urivec(k) * dterm6x(k)
            ticvecy(k) =  uivec (1) * dterm4y(k)
     &                  + uivec (3) * dterm5y(k)
     &                  + urivec(k) * dterm6y(k)
            ticvecz(k) =  uivec (2) * dterm4z(k)
     &                  + uivec (3) * dterm5z(k)
     &                  + urivec(k) * dterm6z(k)
                                             
            tkcvecx(k) =  ukvecx(k) * dterm4x(k)
     &                  + ukvecy(k) * dterm5x(k)
     &                  + urkvec(k) * dterm6x(k)
                                             
            tkcvecy(k) =  ukvecx(k) * dterm4y(k)
     &                  + ukvecz(k) * dterm5y(k)
     &                  + urkvec(k) * dterm6y(k)
                                             
            tkcvecz(k) =  ukvecy(k) * dterm4z(k)
     &                  + ukvecz(k) * dterm5z(k)
     &                  + urkvec(k) * dterm6z(k)
         enddo
c
c         Construct matrixes for dot_product
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tivec1(k) = tisvecx(k) !xx   
            tivec5(k) = tisvecy(k) !yy   
            tivec9(k) = tisvecz(k) !zz   
            tivec2(k) = ticvecx(k) !xy   
            tivec4(k) = ticvecx(k) !yx   
            tivec3(k) = ticvecy(k) !xz   
            tivec7(k) = ticvecy(k) !zx   
            tivec6(k) = ticvecz(k) !yz   
            tivec8(k) = ticvecz(k) !zy   
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            tkvec1(k) = tksvecx(k) !xx   
            tkvec5(k) = tksvecy(k) !yy   
            tkvec9(k) = tksvecz(k) !zz   
            tkvec2(k) = tkcvecx(k) !xy   
            tkvec4(k) = tkcvecx(k) !yx   
            tkvec3(k) = tkcvecy(k) !xz   
            tkvec7(k) = tkcvecy(k) !zx   
            tkvec6(k) = tkcvecz(k) !yz   
            tkvec8(k) = tkcvecz(k) !zy   
         enddo
c
c     Do dot product
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            depxvec(k) =  tivec1(k) * ukpvecx(k)
     &                  + tivec2(k) * ukpvecy(k)
     &                  + tivec3(k) * ukpvecz(k)
     &                  + tkvec1(k) * uipvec (1)
     &                  + tkvec2(k) * uipvec (2)
     &                  + tkvec3(k) * uipvec (3)
            depyvec(k) =  tivec4(k) * ukpvecx(k)
     &                  + tivec5(k) * ukpvecy(k)
     &                  + tivec6(k) * ukpvecz(k)
     &                  + tkvec4(k) * uipvec (1)
     &                  + tkvec5(k) * uipvec (2)
     &                  + tkvec6(k) * uipvec (3)
            depzvec(k) =  tivec7(k) * ukpvecx(k)
     &                  + tivec8(k) * ukpvecy(k)
     &                  + tivec9(k) * ukpvecz(k)
     &                  + tkvec7(k) * uipvec (1)
     &                  + tkvec8(k) * uipvec (2)
     &                  + tkvec9(k) * uipvec (3)


            frcxvec(k) = frcxvec(k) + depxvec(k)
            frcyvec(k) = frcyvec(k) + depyvec(k)
            frczvec(k) = frczvec(k) + depzvec(k)
         enddo

c
c     increment gradient and virial due to Cartesian forces
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            if (k.le.nnelst1) then
              dep(1,i) = dep(1,i) - frcxvec(k)
              dep(2,i) = dep(2,i) - frcyvec(k)
              dep(3,i) = dep(3,i) - frczvec(k)
!!             dep1(i) = dep1(i) - frcxvec(k)
!!             dep2(i) = dep2(i) - frcyvec(k)
!!             dep3(i) = dep3(i) - frczvec(k)
            endif
         enddo
    
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            if (k.le.nnelst1) then
               dep(1,kbisvec1(k)) = dep(1,kbisvec1(k)) + frcxvec(k)
               dep(2,kbisvec1(k)) = dep(2,kbisvec1(k)) + frcyvec(k)
               dep(3,kbisvec1(k)) = dep(3,kbisvec1(k)) + frczvec(k)
!!              dep1(kbisvec1(k)) = dep1(kbisvec1(k)) + frcxvec(k)
!!              dep2(kbisvec1(k)) = dep2(kbisvec1(k)) + frcyvec(k)
!!              dep3(kbisvec1(k)) = dep3(kbisvec1(k)) + frczvec(k)
               vir(1,1) = vir(1,1) + xpsvec(k) * frcxvec(k)
               vir(2,1) = vir(2,1) + ypsvec(k) * frcxvec(k)
               vir(3,1) = vir(3,1) + zpsvec(k) * frcxvec(k)
               vir(1,2) = vir(1,2) + ypsvec(k) * frcxvec(k)
               vir(2,2) = vir(2,2) + ypsvec(k) * frcyvec(k)
               vir(3,2) = vir(3,2) + zpsvec(k) * frcyvec(k)
               vir(1,3) = vir(1,3) + zpsvec(k) * frcxvec(k)
               vir(2,3) = vir(2,3) + zpsvec(k) * frcyvec(k)
               vir(3,3) = vir(3,3) + zpsvec(k) * frczvec(k)
            endif
         enddo
c
c     get the induced dipole field used for dipole torques
c

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k = 1, neloop8
            if (k.le.nnelst1) then
               ti3vecx(k) =    psr3vec(k) * ukvecx (k)
     &                       + dsr3vec(k) * ukpvecx(k)
               ti3vecy(k) =    psr3vec(k) * ukvecy (k)
     &                       + dsr3vec(k) * ukpvecy(k)
               ti3vecz(k) =    psr3vec(k) * ukvecz (k)
     &                       + dsr3vec(k) * ukpvecz(k)
               tk3vecx(k) =    psr3vec(k) * uivec  (1)
     &                       + dsr3vec(k) * uipvec (1)
               tk3vecy(k) =    psr3vec(k) * uivec  (2)
     &                       + dsr3vec(k) * uipvec (2)
               tk3vecz(k) =    psr3vec(k) * uivec  (3)
     &                       + dsr3vec(k) * uipvec (3)
               turi5vec(k) = - psr5vec(k) * urkvec (k)
     &                       - dsr5vec(k) * urkpvec(k)
               turk5vec(k) = - psr5vec(k) * urivec (k)
     &                       - dsr5vec(k) * uripvec(k)
            else
               ti3vecx(k) =  0.0_ti_p
               ti3vecy(k) =  0.0_ti_p
               ti3vecz(k) =  0.0_ti_p
               tk3vecx(k) =  0.0_ti_p
               tk3vecy(k) =  0.0_ti_p
               tk3vecz(k) =  0.0_ti_p
               turi5vec(k) = 0.0_ti_p
               turk5vec(k) = 0.0_ti_p
            endif
         enddo
c
c     get induced dipole field gradient used for quadrupole torques
c
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k = 1, neloop8
            if (k.le.nnelst1) then
               ti5vecx(k) =  2.0_ti_p * (  psr5vec(k) * ukvecx (k)
     &                                + dsr5vec(k) * ukpvecx(k))
               ti5vecy(k) =  2.0_ti_p * (  psr5vec(k) * ukvecy (k)
     &                                + dsr5vec(k) * ukpvecy(k))
               ti5vecz(k) =  2.0_ti_p * (  psr5vec(k) * ukvecz (k)
     &                                + dsr5vec(k) * ukpvecz(k))
               tk5vecx(k) =  2.0_ti_p * (  psr5vec(k) * uivec  (1)
     &                                + dsr5vec(k) * uipvec (1))
               tk5vecy(k) =  2.0_ti_p * (  psr5vec(k) * uivec  (2)
     &                                + dsr5vec(k) * uipvec (2))
               tk5vecz(k) =  2.0_ti_p * (  psr5vec(k) * uivec  (3)
     &                                + dsr5vec(k) * uipvec (3))
               turi7vec(k) = - psr7vec(k) * urkvec (k)
     &                       - dsr7vec(k) * urkpvec(k)
               turk7vec(k) = - psr7vec(k) * urivec (k)
     &                       - dsr7vec(k) * uripvec(k)
            else
               ti5vecx(k) =  0.0_ti_p
               ti5vecy(k) =  0.0_ti_p
               ti5vecz(k) =  0.0_ti_p
               tk5vecx(k) =  0.0_ti_p
               tk5vecy(k) =  0.0_ti_p
               tk5vecz(k) =  0.0_ti_p
               turi7vec(k) = 0.0_ti_p
               turk7vec(k) = 0.0_ti_p
            endif
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
c           if (k.le.nnelst1) then
               t1 = ti3vecx(k) + xpsvec(k) * turi5vec(k)
               t2 = ti3vecy(k) + ypsvec(k) * turi5vec(k)
               t3 = ti3vecz(k) + zpsvec(k) * turi5vec(k)

               ufldvec1(i) =  ufldvec1(i) + t1
               ufldvec2(i) =  ufldvec2(i) + t2
               ufldvec3(i) =  ufldvec3(i) + t3
c           endif
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
            if(k.le.nnelst1) then
               t1 = tk3vecx(k) + xpsvec(k) * turk5vec(k)
               t2 = tk3vecy(k) + ypsvec(k) * turk5vec(k)
               t3 = tk3vecz(k) + zpsvec(k) * turk5vec(k)

               ufldvec1(kbisvec1(k)) = ufldvec1(kbisvec1(k)) + t1
               ufldvec2(kbisvec1(k)) = ufldvec2(kbisvec1(k)) + t2
               ufldvec3(kbisvec1(k)) = ufldvec3(kbisvec1(k)) + t3
            endif
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
               t1 = xpsvec(k) * ti5vecx (k) + xpsvec(k)**2 * turi7vec(k)
               t3 = ypsvec(k) * ti5vecy (k) + ypsvec(k)**2 * turi7vec(k)
               t6 = zpsvec(k) * ti5vecz (k) + zpsvec(k)**2 * turi7vec(k)

               dufldvec1(i) =  dufldvec1(i) + t1
               dufldvec3(i) =  dufldvec3(i) + t3
               dufldvec6(i) =  dufldvec6(i) + t6
         enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, neloop8
               t2 =  xpsvec(k) * ti5vecy(k) + ypsvec(k) * ti5vecx (k)
     &             + 2.0_ti_p     * xpsvec (k) * ypsvec(k) * turi7vec(k)
               t4 =  xpsvec(k) * ti5vecz(k) + zpsvec(k) * ti5vecx (k)
     &             + 2.0_ti_p     * xpsvec (k) * zpsvec(k) * turi7vec(k)
               t5 =  ypsvec(k) * ti5vecz(k) + zpsvec(k) * ti5vecy (k)
     &             + 2.0_ti_p     * ypsvec (k) * zpsvec(k) * turi7vec(k)

               dufldvec2(i) =  dufldvec2(i) + t2
               dufldvec4(i) =  dufldvec4(i) + t4
               dufldvec5(i) =  dufldvec5(i) + t5
         enddo

!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, neloop8
            if (k.le.nnelst1) then
               t1 = xpsvec(k) * tk5vecx(k) + xpsvec(k)**2 * turk7vec(k)
               t3 = ypsvec(k) * tk5vecy(k) + ypsvec(k)**2 * turk7vec(k)
               t6 = zpsvec(k) * tk5vecz(k) + zpsvec(k)**2 * turk7vec(k)

               dufldvec1(kbisvec1(k)) =  dufldvec1(kbisvec1(k)) - t1
               dufldvec3(kbisvec1(k)) =  dufldvec3(kbisvec1(k)) - t3
               dufldvec6(kbisvec1(k)) =  dufldvec6(kbisvec1(k)) - t6
           endif
        enddo
!DIR$ ASSUME (mod(neloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
        do k = 1, neloop8
           if (k.le.nnelst1) then
               t2 =  xpsvec(k) * tk5vecy(k) + ypsvec(k) * tk5vecx (k)
     &             + 2.0_ti_p     * xpsvec (k) * ypsvec(k) * turk7vec(k)

               t4 =  xpsvec(k) * tk5vecz(k) + zpsvec(k) * tk5vecx (k)
     &             + 2.0_ti_p     * xpsvec (k) * zpsvec(k) * turk7vec(k)

               t5 =  ypsvec(k) * tk5vecz(k) + zpsvec(k) * tk5vecy (k)
     &             + 2.0_ti_p     * ypsvec(k)  * zpsvec(k) * turk7vec(k)

               dufldvec2(kbisvec1(k)) =  dufldvec2(kbisvec1(k)) - t2
               dufldvec4(kbisvec1(k)) =  dufldvec4(kbisvec1(k)) - t4
               dufldvec5(kbisvec1(k)) =  dufldvec5(kbisvec1(k)) - t5
            endif
         enddo

      end do  MAINLOOP
c
c     torque is induced field and gradient cross permanent moments
c

!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
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

!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         dipolvecx(k)  = rpole(2,iipolevec(k))
         dipolvecy(k)  = rpole(3,iipolevec(k))
         dipolvecz(k)  = rpole(4,iipolevec(k))
c    
c       qipolvec  is qixx, qixy, qixz for all the npolelocnl
c        1     2     3     4     5     6     7     8     9
c       qixx, qixy, qixz, qixy, qiyy, qiyz, qixz, qiyz, qizz

         qipolvec1(k) = rpole(  5, iipolevec(k))
         qipolvec2(k) = rpole(  6, iipolevec(k))
         qipolvec3(k) = rpole(  7, iipolevec(k))
         qipolvec4(k) = rpole(  6, iipolevec(k))
         qipolvec5(k) = rpole(  9, iipolevec(k))
         qipolvec6(k) = rpole( 10, iipolevec(k))
         qipolvec7(k) = rpole(  7, iipolevec(k))
         qipolvec8(k) = rpole( 10, iipolevec(k))
         qipolvec9(k) = rpole( 13, iipolevec(k))
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         trqvecx(k) =           dipolvecz(k)  *     ufldvec2(ivec(k)) 
     &                -         dipolvecy(k)  *     ufldvec3(ivec(k)) 
     &                +         qipolvec3(k)  *    dufldvec2(ivec(k)) 
     &                -         qipolvec2(k)  *    dufldvec4(ivec(k)) 
     &                + 2.0_ti_p * qipolvec6(k)  * (  dufldvec3(ivec(k))
     &                                           - dufldvec6(ivec(k)) )
     &                +         (qipolvec9(k) -    qipolvec5     (k))
     &                                        *    dufldvec5(ivec(k)) 
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         trqvecy(k) =           dipolvecx(k)  *     ufldvec3(ivec(k)) 
     &                -         dipolvecz(k)  *     ufldvec1(ivec(k)) 
     &                -         qipolvec6(k)  *    dufldvec2(ivec(k)) 
     &                +         qipolvec2(k)  *    dufldvec5(ivec(k)) 
     &                + 2.0_ti_p * qipolvec3(k)  * (  dufldvec6(ivec(k))
     &                                           - dufldvec1(ivec(k)) )
     &                +         (qipolvec1(k) -    qipolvec9     (k))
     &                                        *    dufldvec4(ivec(k)) 
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
         trqvecz(k) =           dipolvecy(k)  *     ufldvec1(ivec(k)) 
     &                -         dipolvecx(k)  *     ufldvec2(ivec(k)) 
     &                +         qipolvec6(k)  *    dufldvec4(ivec(k)) 
     &                -         qipolvec3(k)  *    dufldvec5(ivec(k)) 
     &                + 2.0_ti_p * qipolvec2(k)  * (  dufldvec1(ivec(k))
     &                                           - dufldvec3(ivec(k)) )
     &                +         (qipolvec5(k) -    qipolvec1     (k))
     &                                        *    dufldvec2(ivec(k)) 
      enddo

!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
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
     &                 dep(1,1:nbloc),
     &                 dep(2,1:nbloc),
     &                 dep(3,1:nbloc)
!!    &                 dep1(1:nbloc),
!!    &                 dep2(1:nbloc),
!!    &                 dep3(1:nbloc)
     &                 )

!DIR$ ASSUME (mod(npolelocnlloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, npolelocnlloop
         iaxvec(k) = xaxis(iipolevec(k))
         iayvec(k) = yaxis(iipolevec(k))
         iazvec(k) = zaxis(iipolevec(k))
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

!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         if (iaxvec(k).gt.0)  then
            xivecx(k) = x(iaxvec(k)) - x(iglobvec(k)) ! xix
            yivecx(k) = y(iaxvec(k)) - y(iglobvec(k)) ! yix
            zivecx(k) = z(iaxvec(k)) - z(iglobvec(k)) ! zix
         else
            xivecx(k) = 0.0_ti_p ! xix
            yivecx(k) = 0.0_ti_p ! yix
            zivecx(k) = 0.0_ti_p ! zix
         endif
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         if (iayvec(k).gt.0)  then
            xivecy(k) = x(iayvec(k)) - x(iglobvec(k)) ! xiy
            yivecy(k) = y(iayvec(k)) - y(iglobvec(k)) ! yiy
            zivecy(k) = z(iayvec(k)) - z(iglobvec(k)) ! ziy
         else
            xivecy(k) = 0.0_ti_p ! xiy
            yivecy(k) = 0.0_ti_p ! yiy
            zivecy(k) = 0.0_ti_p ! ziy
         endif
      enddo
!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ NOFUSION
      do k = 1, npolelocnlloop
         if (iazvec(k).gt.0) then
            xivecz(k) = x(iazvec(k)) - x(iglobvec(k)) ! xiz
            yivecz(k) = y(iazvec(k)) - y(iglobvec(k)) ! yiz
            zivecz(k) = z(iazvec(k)) - z(iglobvec(k)) ! ziz
         else
            xivecz(k) = 0.0_ti_p ! xiz
            yivecz(k) = 0.0_ti_p ! yiz
            zivecz(k) = 0.0_ti_p ! ziz
         endif
      enddo

!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
            vxx = vxx + xivecx(k) * fixvecx(k) + xivecy(k) * fiyvecx(k)
     &                + xivecz(k) * fizvecx(k)
            vxy = vxy + yivecx(k) * fixvecx(k) + yivecy(k) * fiyvecx(k)
     &                + yivecz(k) * fizvecx(k)
      enddo

!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
            vxz = vxz + zivecx(k) * fixvecx(k) + zivecy(k) * fiyvecx(k)
     &                + zivecz(k) * fizvecx(k)
            vyy = vyy + yivecx(k) * fixvecy(k) + yivecy(k) * fiyvecy(k)
     &                + yivecz(k) * fizvecy(k)
      enddo

!DIR$ ASSUME (mod(npolelocnlloop,8).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, npolelocnlloop
            vyz = vyz + zivecx(k) * fixvecy(k) + zivecy(k) * fiyvecy(k)
     &                + zivecz(k) * fizvecy(k)
            vzz = vzz + zivecx(k) * fixvecz(k) + zivecy(k) * fiyvecz(k)
     &                + zivecz(k) * fizvecz(k)
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

!        dep(1,:) = dep1
!        dep(2,:) = dep2
!        dep(3,:) = dep3
      call timer_exit( timer_epreal )
      return
      end
