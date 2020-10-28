c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_precision.h"
      module precompute_pole
        implicit none

        integer :: nprecomp
        integer, allocatable,target :: kkpole_precompute(:)
        integer, pointer :: kpoleloc_precompute(:)
        integer, allocatable,target :: iipole_precompute(:)
        integer, pointer :: iploc_precompute(:)
        real(t_p), allocatable :: rx_precompute(:) ! Precomputed distance vectors (x)
        real(t_p), allocatable :: ry_precompute(:) ! Precomputed distance vectors (y)
        real(t_p), allocatable :: rz_precompute(:) ! Precomputed distance vectors (z)
        real(t_p), pointer :: ms_precompute(:) ! m scaling factor
        real(t_p), allocatable,target :: rr3_bn1_precompute(:) ! Precomputed distance 
        real(t_p), allocatable,target :: rr5_bn2_precompute(:) ! Precomputed distance  
        logical :: precompute_solvpole=.false.
        logical :: precompute_mpole=.false.
        logical :: precompute_tmat=.true.

        include "erfcore_data.f.inc"

      contains

#include "image.f.inc"
#if defined(SINGLE) | defined(MIXED)
      include "erfcscore.f.inc"
#else
      include "erfcdcore.f.inc"
#endif

      subroutine allocate_precompute_data
      use inform ,only:deb_Path
      use mpole
      use neigh
      use utilgpu,only:dir_queue
      implicit none
      logical :: do_deallocate, do_allocate
      logical,save :: first_in=.true.
      integer i,ninteract

      ninteract = 0
!$acc parallel loop async(dir_queue)
      do i = 1,npolelocnl
         ninteract = ninteract + nelstc(i)
      end do
!$acc wait

      ! Allocate if not allocated and reallocate if smaller than needed
      if( allocated(rx_precompute) ) then
        do_deallocate = (size(rx_precompute) < ninteract)
        do_allocate   = do_deallocate
      else
        do_deallocate = .false.
        do_allocate   = .true.
      endif
      if( do_deallocate ) then
!$acc exit data
!$acc& delete(rx_precompute, ry_precompute, rz_precompute)
!$acc& delete(rr3_bn1_precompute, rr5_bn2_precompute)
!$acc& delete(iipole_precompute,kkpole_precompute)
!$acc& detach(ms_precompute,kpoleloc_precompute,iploc_precompute)

        nullify(ms_precompute)
        nullify(iploc_precompute)
        nullify(kpoleloc_precompute)
        deallocate(iipole_precompute)
        deallocate(kkpole_precompute)
        deallocate(rx_precompute)
        deallocate(ry_precompute)
        deallocate(rz_precompute)
        deallocate(rr3_bn1_precompute)
        deallocate(rr5_bn2_precompute)
      endif
      if( do_allocate ) then
        allocate( iipole_precompute(ninteract) )
        allocate( kkpole_precompute(ninteract) )
        allocate( rx_precompute(ninteract) )
        allocate( ry_precompute(ninteract) )
        allocate( rz_precompute(ninteract) ) 
        allocate( rr3_bn1_precompute(ninteract) )
        allocate( rr5_bn2_precompute(ninteract) )
        ms_precompute       => rr3_bn1_precompute
        kpoleloc_precompute => kkpole_precompute
        iploc_precompute    => iipole_precompute

        if (deb_Path) then
#if defined(SINGLE) | defined(MIXED)
        print *,"tmatxb Temporaries : 7*float*",shape(rx_precompute),
     &  ":", (7.*4.*size(rx_precompute))/(1024.*1024.*1024.), "Go"
#else
        print *,"tmatxb Temporaries : 6*double*",shape(rx_precompute),
     &  ":", (6.*8.*size(rx_precompute))/(1024.*1024.*1024.), "Go"
#endif
        end if

!$acc enter data
!$acc& create(rx_precompute, ry_precompute, rz_precompute)
!$acc& create(rr3_bn1_precompute, rr5_bn2_precompute) 
!$acc& create(iipole_precompute,kkpole_precompute)
!$acc& attach(ms_precompute,kpoleloc_precompute)
!$acc& attach(iploc_precompute)

      end if
      end subroutine

      ! Precompute multipole interactions via neigbhor list
      subroutine mpole_precompute()
      use atoms   , only : x, y, z
      use atmlst  , only : poleglobnl
      use couple  , only : allscal_n,typscal_n,numscal_n,scalbeg_n
      use chgpot  , only : dielec,electric
      use ewald   , only : aewald
      use math    , only : sqrtpi
      use utilgpu , only : def_queue
      use mpole   , only : ipole,poleloc,npolelocnl
      use mplpot  , only : m2scale,m3scale,m4scale,m5scale
      use neigh   , only : nelst, elst
      use polar   , only : polarity, thole, pdamp
      use shunt   , only : cut2
      use utilgpu , only : maxscaling1,maxscaling,dir_queue,warning
      use domdec  , only : loc, nlocnl, nbloc
      use polpot  , only : u1scale, u2scale, u3scale, u4scale
      use sizes   , only : maxelst
      use tinheader,only : ti_p
      implicit none

      integer :: j, k                
      integer :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      real(t_p) :: rx, ry, rz             ! Difference between atoms positions
      real(t_p) :: d1, d2                 ! Distance (d2 squared) between atoms      
      real(t_p) :: alsq2, alsq2n          ! ewald precomputed constants
      real(t_p) :: bn0,sc3,sc5,ralpha,exp2a
      real(t_p) :: sdamp,sdamp1,expdamp1,pgamma
      real(t_p) :: bn1,bn2,rr3,rr5
      integer   :: imscal(maxscaling)  ! Temporary array to store interaction type for neighbors
      real(t_p) :: fmscal(maxscaling)  ! Temporary array to store scaling factors for neighbors
      integer nn12,nn13,nn14,ntot,km
      real(t_p) mscale ! Scaling factor for interaction
      real(t_p) :: f
      real(t_p),save:: mscalevec(5)
      logical,save:: f_in=.true.
      integer :: nprecomp_capture

!!$acc update host(poleglobnl,polarity,loc,ipole,nelst)
      call allocate_precompute_data
      if (f_in) then
         mscalevec = [0.0_ti_p,m2scale,m3scale,m4scale,m5scale]
!$acc enter data copyin(mscalevec)
         f_in = .false.
      end if

!$acc data present(poleglobnl,ipole,poleloc)
!$acc&     present(polarity,thole,pdamp)
!$acc&     present(rx_precompute, ry_precompute, rz_precompute,
!$acc& rr3_bn1_precompute, rr5_bn2_precompute,
!$acc& ms_precompute,iipole_precompute,kkpole_precompute)
!$acc&     present(numscal_n,scalbeg_n,allscal_n,typscal_n)
!$acc&     present(mscalevec)

      ! gather some parameters, then set up the damping factors.
      nprecomp = 0

!$acc parallel loop gang vector_length(32) copy(nprecomp)
!$acc&         private(km,imscal,fmscal)
!$acc&         async(dir_queue)
      PRECOMPUTE_LOOP:
     &do ii = 1, npolelocnl
           iipole = poleglobnl(ii)
           iglob  = ipole  (iipole)
           if (loc(iglob) == 0 .or. loc(iglob)>nbloc) then
              print*,warning
              cycle PRECOMPUTE_LOOP           
           end if

           km    = 0
           if (ntot.gt.maxscaling) 
     &        print*,'scaling array too short mpole precompute'

           ntot = numscal_n(iglob)
           nn12 = scalbeg_n(iglob)
!$acc loop vector
           do j = 1,ntot
              imscal(j) = allscal_n(nn12+j)
              fmscal(j) = mscalevec(typscal_n(nn12+j))
           end do

           ! Loop on neighbors
!$acc loop vector
           do k =  1, nelst(ii)
              kpole    = elst(k,ii)
              kglob    = ipole(kpole)

              mscale   = 1
              if (km<ntot) then
!$acc loop seq
                 do j=1,ntot
                    if (imscal(j)==kglob) then
                       mscale = fmscal(j)
!$acc atomic update
                       km  = km+1
                       exit
                    end if
                 end do
              end if

              rx = x(kglob) - x(iglob)
              ry = y(kglob) - y(iglob)
              rz = z(kglob) - z(iglob)
              call image_inl(rx,ry,rz)
              d2 = rx*rx + ry*ry + rz*rz
              if (d2 > cut2) cycle

!$acc atomic capture
              nprecomp = nprecomp + 1
              nprecomp_capture = nprecomp
!$acc end atomic

              iipole_precompute(nprecomp_capture) = iipole
              kkpole_precompute(nprecomp_capture) = kpole
              rx_precompute    (nprecomp_capture) = rx
              ry_precompute    (nprecomp_capture) = ry
              rz_precompute    (nprecomp_capture) = rz
              ms_precompute    (nprecomp_capture) = 1-mscale
           enddo
      end do PRECOMPUTE_LOOP
!$acc end data
!$acc wait(def_queue)

      end subroutine

      ! Compute the direct space contribution to the electric field due to the current value of the induced dipoles
      subroutine tmatxb_precompute()
      use atoms   , only : x, y, z
      use atmlst  , only : poleglobnl
      use chgpot  , only : dielec,electric
      use ewald   , only : aewald
      use math    , only : sqrtpi
      use utilgpu , only : def_queue
      use mpole   , only : ipole,poleloc,npolelocnl
      use mplpot  , only : m2scale,m3scale,m4scale,m5scale
      use neigh   , only : nelst, elst, nelstc
      use polar   , only : polarity, thole, pdamp
      use shunt   , only : cut2
      use utilgpu , only : maxscaling1,maxscaling,dir_queue,warning
      use domdec  , only : loc, nlocnl, nbloc
      use polgrp  , only : allscal_p,typscal_p,numscal_p,scalbeg_p
      use polpot  , only : u1scale, u2scale, u3scale, u4scale
      use sizes   , only : maxelst
      use tinheader,only : ti_p
      implicit none

      integer   :: j, k
      integer   :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer   :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      real(t_p) :: rx, ry, rz               ! Difference between atoms positions
      real(t_p) :: d1, d2                   ! Distance (d2 squared) between atoms      
      real(t_p) :: alsq2, alsq2n            ! ewald precomputed constants
      real(t_p) :: bn0,sc3,sc5,ralpha,exp2a
      real(t_p) :: sdamp,sdamp1,expdamp1,pgamma
      real(t_p) :: bn1,bn2,rr3,rr5
      integer   :: iscal(maxscaling1)   ! Temporary array to store interaction type for neighbors
      real(t_p) :: fscal(maxscaling1)   ! Temporary array to store scaling factors for neighbors
      integer nnp11,nnp12,nnp13,nnp14,ki
      real(t_p),parameter:: uscale=1.0   ! Scaling factor for interaction
      real(t_p) uscalevec(5)
      real(t_p) :: f
      integer   :: nprecomp_capture
      logical,save:: f_in=.true.

      ! No need to reallocate if it has already been done precomputing
      ! mpole
      if (.not.precompute_mpole)
     &   call allocate_precompute_data
      !FIXME Scaling factor containers are not on device anymore

!$acc data present(poleglobnl,ipole,poleloc)
!$acc&     present(polarity,thole,pdamp,elst,nelst,nelstc)
!$acc&     present(iploc_precompute, kpoleloc_precompute)
!$acc&     present(rx_precompute, ry_precompute, rz_precompute,
!$acc&  rr3_bn1_precompute, rr5_bn2_precompute)
!$acc&     present(numscal_p,scalbeg_p,allscal_p,typscal_p)
c!$acc data create(uscalevec) async(dir_queue)

      ! gather some parameters, then set up the damping factors.
      call switch ('EWALD')
      alsq2    = 2 * aewald**2
      alsq2n   = 0
      f        = electric / dielec
      nprecomp = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)
c!$acc serial async(def_queue)
c      uscalevec = [u1scale,u2scale,u3scale,u4scale,0.0]
c!$acc end serial

!$acc parallel loop gang vector_length(32) copy(nprecomp)
!$acc&         async(dir_queue)
c!$acc&         private(ki,iscal,fscal)
      PRECOMPUTE_LOOP:
     &do ii = 1, npolelocnl
           iipole = poleglobnl(ii)
           !skip atom if it is not polarizable
           if (polarity(iipole) == 0) cycle PRECOMPUTE_LOOP
           iglob  = ipole  (iipole)
           iploc  = poleloc(iipole)
           if (loc(iglob) == 0) cycle PRECOMPUTE_LOOP           

c           ki     = 0
c           nnp14 = numscal_p(iglob)
c           nnp12 = scalbeg_p(iglob)
c!$acc loop vector
c           do j = 1,nnp14
c              iscal(j) = allscal_p(nnp12+j)
c              fscal(j) = uscalevec(typscal_p(nnp12+j))
c           end do

           ! Loop on neighbors
!$acc loop vector
           do k =  1, nelst(ii)
              kpole    = elst(k,ii)
              kpoleloc = poleloc(kpole)
              if (kpoleloc.eq.0) cycle
              kglob    = ipole(kpole)

c              uscale   = 1
c              if (ki<nnp14) then
c!$acc   loop seq
c                 do j=1,nnp14
c                    if (iscal(j)==kglob) then
c                       uscale  = fscal(j)
c!$acc   atomic update
c                       ki  = ki+1
c                       exit
c                    end if
c                 end do
c              end if

              rx = x(kglob) - x(iglob)
              ry = y(kglob) - y(iglob)
              rz = z(kglob) - z(iglob)
              call image_inl(rx,ry,rz)
              d2 = rx*rx + ry*ry + rz*rz
              if (d2 > cut2) cycle

              ! compute the distances and the scaling factors according to Thole's model.
              d1      = sqrt(d2)

              ralpha  = aewald * d1
              exp2a   = exp(-ralpha*ralpha)
              bn0     = erfc(ralpha)
              bn0     = bn0 / d1
              bn1     = (     bn0 +         alsq2  *alsq2n * exp2a) / d2
              bn2     = ( 3 * bn1 + alsq2 * alsq2 * alsq2n * exp2a) / d2
              sdamp   = pdamp(iipole) * pdamp(kpole)
              pgamma  = min(thole(iipole),thole(kpole))

              if( sdamp == 0 ) then
                sdamp1   = -100.0
                sc3      =   1 - exp(sdamp1) * uscale
                sc5      =   1 - exp(sdamp1) * uscale * (1 - sdamp1)
              else
                sdamp1 = - pgamma * (d1 / sdamp) ** 3
                if (sdamp1 > -50) then
                  expdamp1 = exp(sdamp1)
                  sc3      =   1 - expdamp1 * uscale
                  sc5      =   1 - expdamp1 * uscale
     &                      * (1 - sdamp1)
                else
                  sc3     = 1
                  sc5     = 1
                end if
              endif

              ! compute the field.
              rr3     =     (1 - sc3) / (d1 * d2)      
              rr5     = 3 * (1 - sc5) / (d1 * d2 * d2)

!$acc atomic capture
              nprecomp = nprecomp + 1
              nprecomp_capture = nprecomp
!$acc end atomic

              iploc_precompute (nprecomp_capture) = iploc
              kpoleloc_precompute (nprecomp_capture) = kpoleloc
              rx_precompute     (nprecomp_capture) = rx
              ry_precompute     (nprecomp_capture) = ry
              rz_precompute     (nprecomp_capture) = rz
              rr3_bn1_precompute(nprecomp_capture) = rr3 - bn1    
              rr5_bn2_precompute(nprecomp_capture) = rr5 - bn2
           enddo

      end do PRECOMPUTE_LOOP
c!$acc end data
!$acc end data
!$acc wait
      end subroutine

      subroutine tmatxb_pme_compute(mu,efi)
      use atmlst  , only :  poleglobnl
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolebloc,npolelocnl
      use neigh   , only : nelst, elst
      use polar   , only : polarity,  tinypol
      use shunt   , only : cut2
      use utilgpu , only : def_queue
      use domdec  , only : loc
      implicit none
      real(t_p) ,intent(in)   ::  mu(:,:,:)
      real(t_p) ,intent(inout):: efi(:,:,:)

      integer ::  k                
      integer :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer :: kpoleloc   ! Indexes associated to neighbor atom
      real(t_p) :: rx, ry, rz             ! Difference between atoms positions
      real(t_p) :: d2                 ! Distance (d2 squared) between atoms      
      real(t_p) :: duix, puix ! coefficients of mu(iploc)
      real(t_p) :: duiy, puiy ! coefficients of mu(iploc)
      real(t_p) :: duiz, puiz ! coefficients of mu(iploc)
      real(t_p) :: dukx, pukx ! coefficients of mu(kpoleloc)
      real(t_p) :: duky, puky ! coefficients of mu(kpoleloc)
      real(t_p) :: dukz, pukz ! coefficients of mu(kpoleloc)      
      real(t_p) :: duir, dukr, puir, pukr ! Scalar products duir = (du(i).r)
      real(t_p) :: rr3_bn1, rr5_bn2
      real(t_p) :: fidx, fipx ! Contribution of local pair to efi(iploc)
      real(t_p) :: fidy, fipy ! Contribution of local pair to efi(iploc)
      real(t_p) :: fidz, fipz ! Contribution of local pair to efi(iploc)
      real(t_p) :: fkdx, fkpx ! Contribution of local pair to efi(kpoleloc)
      real(t_p) :: fkdy, fkpy ! Contribution of local pair to efi(kpoleloc)
      real(t_p) :: fkdz, fkpz ! Contribution of local pair to efi(kpoleloc)   
      
c
c     Initiate precompute
c
      if (precompute_tmat) then
         call tmatxb_precompute()
         precompute_tmat=.false.
      end if

!$acc parallel loop async(def_queue)
!$acc&         present(poleglobnl,ipole,poleloc,polarity,mu,efi)
!$acc&         present(rx_precompute,ry_precompute,rz_precompute
!$acc&  ,rr3_bn1_precompute, rr5_bn2_precompute
!$acc&  ,iploc_precompute, kpoleloc_precompute)
      MAINLOOP:
     &do ii = 1, nprecomp

         rx = rx_precompute(ii)
         ry = ry_precompute(ii)
         rz = rz_precompute(ii)
         kpoleloc = kpoleloc_precompute(ii)
         iploc    = iploc_precompute(ii)
         rr3_bn1  = rr3_bn1_precompute(ii)   
         rr5_bn2  = rr5_bn2_precompute(ii)              

         dukx    = mu(1,1,kpoleloc)
         duky    = mu(2,1,kpoleloc)
         dukz    = mu(3,1,kpoleloc)
         pukx    = mu(1,2,kpoleloc)
         puky    = mu(2,2,kpoleloc)
         pukz    = mu(3,2,kpoleloc)
         duix    = mu(1,1,iploc)
         duiy    = mu(2,1,iploc)
         duiz    = mu(3,1,iploc)
         puix    = mu(1,2,iploc)
         puiy    = mu(2,2,iploc)
         puiz    = mu(3,2,iploc)                            

         duir    =   duix * rx + duiy * ry + duiz * rz
         dukr    =   dukx * rx + duky * ry + dukz * rz

         puir    =   puix * rx + puiy * ry + puiz * rz
         pukr    =   pukx * rx + puky * ry + pukz * rz

         fidx    = - rr3_bn1  * dukx + rr5_bn2  * dukr * rx
         fidy    = - rr3_bn1  * duky + rr5_bn2  * dukr * ry
         fidz    = - rr3_bn1  * dukz + rr5_bn2  * dukr * rz

         fkdx    = - rr3_bn1  * duix + rr5_bn2  * duir * rx
         fkdy    = - rr3_bn1  * duiy + rr5_bn2  * duir * ry
         fkdz    = - rr3_bn1  * duiz + rr5_bn2  * duir * rz

         fipx    = - rr3_bn1  * pukx + rr5_bn2  * pukr * rx
         fipy    = - rr3_bn1  * puky + rr5_bn2  * pukr * ry
         fipz    = - rr3_bn1  * pukz + rr5_bn2  * pukr * rz

         fkpx    = - rr3_bn1  * puix + rr5_bn2  * puir * rx
         fkpy    = - rr3_bn1  * puiy + rr5_bn2  * puir * ry
         fkpz    = - rr3_bn1  * puiz + rr5_bn2  * puir * rz
         
         ! increment electric field for each atoms 
!$acc atomic
         efi(1,1,iploc)    = efi(1,1,iploc) + fidx
!$acc atomic
         efi(2,1,iploc)    = efi(2,1,iploc) + fidy
!$acc atomic
         efi(3,1,iploc)    = efi(3,1,iploc) + fidz  
!$acc atomic
         efi(1,2,iploc)    = efi(1,2,iploc) + fipx
!$acc atomic
         efi(2,2,iploc)    = efi(2,2,iploc) + fipy
!$acc atomic
         efi(3,2,iploc)    = efi(3,2,iploc) + fipz
!$acc atomic
         efi(1,1,kpoleloc) = efi(1,1,kpoleloc) + fkdx
!$acc atomic
         efi(2,1,kpoleloc) = efi(2,1,kpoleloc) + fkdy
!$acc atomic
         efi(3,1,kpoleloc) = efi(3,1,kpoleloc) + fkdz 
!$acc atomic
         efi(1,2,kpoleloc) = efi(1,2,kpoleloc) + fkpx
!$acc atomic
         efi(2,2,kpoleloc) = efi(2,2,kpoleloc) + fkpy   
!$acc atomic
         efi(3,2,kpoleloc) = efi(3,2,kpoleloc) + fkpz

      end do MAINLOOP
      end subroutine
      end module
