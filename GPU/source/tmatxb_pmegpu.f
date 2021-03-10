c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_precision.h"
      module tmatxb_inl_subroutines
        use utilgpu ,only: real3,real6
        integer:: icall=1,ndec=1
        include "erfcore_data.f.inc"
        real(t_p),allocatable::mu3(:,:),mu4(:,:)
        real(t_p),allocatable::efi3(:,:),efi4(:,:)
        contains
#include "image.f.inc"
#if defined(SINGLE) | defined(MIXED)
        include "erfcscore.f.inc"
#else
        include "erfcdcore.f.inc"
#endif
        ! Compute one tmat scalar product
#include "pair_tmatxb.f.inc"
      end module

c=========================================================================================
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c=========================================================================================
      subroutine tmatxb_pmegpu(nrhs,dodiag,mu,efi)
        use atmlst  , only : poleglob
        use mpole   , only : npoleloc,npolebloc
        use polar   , only : polarity, tinypol
        use utilgpu , only : def_queue, dir_queue, rec_queue
#ifdef _OPENACC
     &              , def_stream,dir_stream
#endif
        use domdec  , only : rank
        use inform  , only : deb_Path
        use interfaces ,only : tmatxb_pme_core2,tmatxb_pme_core3,
     &                  tmatxb_pme_core_p
        use tinheader,only: ti_p
        use timestat, only : timer_tmatxb_pmegpu, timer_enter,timer_exit
        use tmatxb_inl_subroutines , only : ndec, icall
        use utils   , only : set_to_zero1
        use utilcomm, only : skpPcomm
        use potent  , only : use_polarshortreal
        use precompute_pole
        implicit none
        integer   , intent(in) :: nrhs
        logical   , intent(in) :: dodiag
        real(t_p) , intent(in) :: mu(:,:,:)
        real(t_p) , intent(out):: efi(:,:,:)
        integer   i,iipole,irhs,j
        integer npoleloc_e
        real(t_p) ipolar
        real(t_p) tmp

        if(deb_Path) write(*,'(4x,a)') 'tmatxb_pmegpu'
        call timer_enter( timer_tmatxb_pmegpu )
#ifdef _OPENACC
        def_queue  = dir_queue
        def_stream = dir_stream
#endif

        !gather some parameters, then set up the damping factors.
        if (use_polarshortreal) then
           call switch ('SHORTEWALD')
           ndec = 1
        else
           call switch ('EWALD     ')
           ! Enable Matvec kernel split
           if (dir_queue.ne.rec_queue) ndec=2
        end if

        if (icall.eq.1)
     &  call set_to_zero1(efi,3*nrhs*npolebloc,def_queue)

        ! core computation of matrice-vector product
        call tmatxb_pme_core_p(mu,efi)

        if (icall.eq.ndec) then
           if (dodiag) then
           npoleloc_e = merge(npolebloc,npoleloc,
     &             (use_polarshortreal.and.skpPcomm))
        ! if dodiag is true, also compute the "self-induced" field,
        ! i.e., the diagonal portion of the matrix/vector product.
!$acc  parallel loop collapse(3) async(def_queue)
!$acc&          present(polarity,poleglob,mu,efi)
           do i = 1, npoleloc_e
              do irhs = 1, nrhs
                 do j = 1, 3
                    ipolar = polarity(poleglob(i))
                    !if no polarisability, take a negligeable value to allow convergence
                    if (ipolar  == 0.0_ti_p) then
                       tmp = mu(j,irhs,i)*(1.0_ti_p/tinypol)
                    else
                       tmp = mu(j,irhs,i)/ipolar
                    end if
!$acc atomic
                    efi(j,irhs,i) = efi(j,irhs,i) + tmp
                 end do
              end do
           end do
           end if
        end if
        call timer_exit( timer_tmatxb_pmegpu )
      end subroutine

      subroutine incMatvecStep
      use tmatxb_inl_subroutines,only:icall
      icall = icall + 1
      end subroutine

      subroutine resetMatvecStep
      use tmatxb_inl_subroutines,only:icall,ndec
      if (icall.ne.ndec) then
13       format( ' Error Made spliting Matrix-Vector product '
     &        ,/,' Either incomplete of over calculated '
     &        ,/,' ndec  ', I4
     &        ,/,' icall ', I4 )
         write(*, 13) ndec, ncall
         call fatal
      end if
      icall = 1
      end subroutine

      subroutine tmatxb_pmegpu1(nrhs,dodiag,mu,efi)
      use atoms   , only : x, y, z
      use atmlst  , only : poleglob, poleglobnl
      use ewald   , only : aewald
      use inform  , only : deb_Path
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolebloc,npoleloc,npolelocnl
      use neigh   , only : nelst, elst,nelstc,nshortelst, shortelst
      use polar   , only : polarity, thole, pdamp, tinypol
      use shunt   , only : cut2
      use utilgpu , only : maxscaling1, def_queue, dir_queue
      use domdec  , only : rank, loc
      use polgrp  , only : np11, np12, np13, np14, ip11, ip12, ip13,ip14
      use polpot  , only : u1scale, u2scale, u3scale, u4scale
      use potent  , only : use_polarshortreal
      use timestat, only : timer_tmatxb_pmegpu, timer_enter, timer_exit
      use tmatxb_inl_subroutines, only : image_inl
      implicit none

      integer   , intent(in) :: nrhs
      logical   , intent(in) :: dodiag
      real(t_p) , intent(in) :: mu(:,:,:)
      real(t_p) , intent(out):: efi(:,:,:)

      integer :: i, j, k, irhs
      integer :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      integer,pointer :: lst(:,:),nlst(:) ! pointers to neighbor list
      real(t_p) :: rx, ry, rz             ! Difference between atoms positions
      real(t_p) :: d1, d2                 ! Distance (d2 squared) between atoms
      real(t_p) :: duix, puix ! coefficients of mu(iploc)
      real(t_p) :: duiy, puiy ! coefficients of mu(iploc)
      real(t_p) :: duiz, puiz ! coefficients of mu(iploc)
      real(t_p) :: dukx, pukx ! coefficients of mu(kpoleloc)
      real(t_p) :: duky, puky ! coefficients of mu(kpoleloc)
      real(t_p) :: dukz, pukz ! coefficients of mu(kpoleloc)
      real(t_p) :: duir, dukr, puir, pukr ! Scalar products duir = (du(i).r)
      real(t_p) :: alsq2, alsq2n ! ewald precomputed constants
      real(t_p) :: bn0,sc3,sc5,ralpha,exp2a,tmp
      real(t_p) :: sdamp,sdamp1,expdamp1,pgamma
      real(t_p) :: bn1,bn2,rr3,rr5
      real(t_p) :: rr3_bn1, rr5_bn2
      real(t_p) :: fidx, fipx ! Contribution of local pair to efi(iploc)
      real(t_p) :: fidy, fipy ! Contribution of local pair to efi(iploc)
      real(t_p) :: fidz, fipz ! Contribution of local pair to efi(iploc)
      real(t_p) :: fkdx, fkpx ! Contribution of local pair to efi(kpoleloc)
      real(t_p) :: fkdy, fkpy ! Contribution of local pair to efi(kpoleloc)
      real(t_p) :: fkdz, fkpz ! Contribution of local pair to efi(kpoleloc)

      integer :: iscal(maxscaling1)   ! Temporary array to store interaction type for neighbors
      real(t_p) :: fscal(maxscaling1) ! Temporary array to store scaling factors for neighbors
      integer nnp11,nnp12,nnp13,nnp14,ki
      real(t_p) uscale ! Scaling factor for interaction

      if(deb_Path) write(*,'(4x,a)') 'tmatxb_pmegpu1'
      call timer_enter( timer_tmatxb_pmegpu )
      def_queue = dir_queue

!$acc data present(poleglobnl,ipole,poleloc,polarity,thole,pdamp,
!$acc&  poleglob)
!$acc&     present(mu(1:3,:nrhs,:npolebloc),efi(1:3,1:nrhs,:npolebloc))

      ! initialize the result vector
!$acc kernels async(def_queue)
      efi(1:3,1:nrhs,:npolebloc) = 0
!$acc end kernels

      ! gather some parameters, then set up the damping factors.
      if (use_polarshortreal) then
         call switch ('SHORTEWALD')
          lst =>  shortelst
         nlst => nshortelst
      else
         call switch ('EWALD     ')
          lst =>  elst
         nlst => nelstc
      end if

      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)
c
!$acc parallel loop gang
!$acc&         vector_length(32)
!$acc&         private(ki,iscal,fscal)
!$acc&         present(lst,nlst)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, npolelocnl
           iipole = poleglobnl(ii)
           !skip atom if it is not polarizable
           if (polarity(iipole) == 0) cycle MAINLOOP
           iglob  = ipole  (iipole)
           iploc  = poleloc(iipole)
           if (loc(iglob) == 0) cycle MAINLOOP

           nnp11  = np11(iglob)
           nnp12  = np12(iglob) + nnp11
           nnp13  = np13(iglob) + nnp12
           nnp14  = np14(iglob) + nnp13
           !if (nnp14.gt.maxscaling1) print*,'scaling array to short in tmatxb',ii

           !fill scaling factor along kglob and interaction type
           ki = 0
!$acc loop vector
           do k = 1,nnp12
              if      (k.le.nnp11) then
                 iscal(k) = ip11(k,iglob)
                 fscal(k) = u1scale
              else
                 iscal(k) = ip12(k-nnp11,iglob)
                 fscal(k) = u2scale
              end if
           end do
!$acc loop vector
           do k = nnp12+1,nnp14
              if      (k.le.nnp13) then
                 iscal(k) = ip13(k-nnp12,iglob)
                 fscal(k) = u3scale
              else
                 iscal(k) = ip14(k-nnp13,iglob)
                 fscal(k) = u4scale
              end if
           end do

           ! Loop on neighbors
!$acc loop vector
           do k =  1, nlst(ii)
              kpole = lst(k,ii)
              kpoleloc = poleloc(kpole)
              if (kpoleloc.eq.0) cycle
              kglob = ipole(kpole)
              rx = x(kglob) - x(iglob)
              ry = y(kglob) - y(iglob)
              rz = z(kglob) - z(iglob)
              call image_inl(rx,ry,rz)
              d2 = rx*rx + ry*ry + rz*rz
              if (d2 > cut2) cycle

              ! find exclusion coefficients for connected atoms
              ! TODO : remove this ugly sequential loop
              uscale  = 1
              if (ki<nnp14) then
!$acc   loop seq
                 do j=1,nnp14
                    if (iscal(j)==kglob) then
                       uscale  = fscal(j)
!$acc   atomic update
                       ki  = ki+1
                       exit
                    end if
                 end do
              end if

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
                sdamp1  = -100.0
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
              rr3_bn1 = rr3 - bn1
              rr5_bn2 = rr5 - bn2

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
!$acc   atomic update
              efi(1,1,iploc)    = efi(1,1,iploc) + fidx
!$acc   atomic update
              efi(2,1,iploc)    = efi(2,1,iploc) + fidy
!$acc   atomic update
              efi(3,1,iploc)    = efi(3,1,iploc) + fidz
!$acc   atomic update
              efi(1,2,iploc)    = efi(1,2,iploc) + fipx
!$acc   atomic update
              efi(2,2,iploc)    = efi(2,2,iploc) + fipy
!$acc   atomic update
              efi(3,2,iploc)    = efi(3,2,iploc) + fipz
!$acc   atomic update
              efi(1,1,kpoleloc) = efi(1,1,kpoleloc) + fkdx
!$acc   atomic update
              efi(2,1,kpoleloc) = efi(2,1,kpoleloc) + fkdy
!$acc   atomic update
              efi(3,1,kpoleloc) = efi(3,1,kpoleloc) + fkdz
!$acc   atomic update
              efi(1,2,kpoleloc) = efi(1,2,kpoleloc) + fkpx
!$acc   atomic update
              efi(2,2,kpoleloc) = efi(2,2,kpoleloc) + fkpy
!$acc   atomic update
              efi(3,2,kpoleloc) = efi(3,2,kpoleloc) + fkpz
           enddo

      end do MAINLOOP

      if(dodiag) then
      ! if dodiag is true, also compute the "self-induced" field,
      ! i.e., the diagonal portion of the matrix/vector product.
!$acc parallel loop gang vector async(def_queue)
         do i = 1, npoleloc
          iipole = poleglob(i)
          !if no polarisability, take a negligeable value to allow convergence
          if (polarity(iipole)  == 0) then
             do irhs = 1, nrhs
               do j = 1, 3
                  tmp = mu(j,irhs,i)/tinypol
!$acc atomic update
                  efi(j,irhs,i) = efi(j,irhs,i) + tmp
               end do
             end do
          else
             do irhs = 1, nrhs
               do j = 1, 3
                  tmp = mu(j,irhs,i)/polarity(iipole)
!$acc atomic update
                  efi(j,irhs,i) = efi(j,irhs,i) + tmp
               end do
             end do
          end if
        end do
      end if
!$acc end data
      call timer_exit( timer_tmatxb_pmegpu )

      end subroutine

      subroutine otf_dc_tmatxb_pmegpu(nrhs,dodiag,mu,efi)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use atmlst
      use atoms
      use divcon
      use domdec
      use ewald
      use group
      use inform    ,only: deb_Path
      use interfaces,only: otf_dc_tmatxb_pme_core_p
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use sizes
      use shunt
      use tinheader ,only: ti_p
      use utils
      use utilgpu
      implicit none
      integer  ,intent(in ):: nrhs
      logical  ,intent(in ):: dodiag
      real(t_p),intent(in ):: mu(3,nrhs,npolebloc)
      real(t_p),intent(out):: efi(3,nrhs,npolebloc)

      integer i, j, irhs, iipole
      integer ipoleloc,ipolar
      real(t_p) tmp

      if (deb_Path)
     &   write(*,'(4x,a)') 'otf_dc_tmatxb_pmegpu'

#ifdef _OPENACC
      def_queue  = dir_queue
      def_stream = dir_stream
#endif
      !initialize the result vector
      call set_to_zero1(efi,size(efi),def_queue)

      call otf_dc_tmatxb_pme_core_p(mu,efi)

      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
!$acc  parallel loop collapse(3) async(def_queue)
!$acc&          present(polarity,poleglob,mu,efi)
         do i = 1, npoleloc
            do irhs = 1, nrhs
               do j = 1, 3
                  ipolar = polarity(poleglob(i))
                  !if no polarisability, take a negligeable value to allow convergence
                  if (ipolar  == 0.0_ti_p) then
                     tmp = mu(j,irhs,i)*(1.0_ti_p/tinypol)
                  else
                     tmp = mu(j,irhs,i)/ipolar
                  end if
!$acc atomic
                  efi(j,irhs,i) = efi(j,irhs,i) + tmp
               end do
            end do
         end do
      end if
      end

c======================================================================
c       -----------------------------------
c       Core routine for tmatxb
c       -----------------------------------
c======================================================================
      subroutine tmatxb_pme_core1(mu,efi)
      use atoms   , only : x, y, z
      use atmlst  , only : poleglobnl
      !use erf_mod
      use ewald   , only : aewald
      use inform  , only : deb_Path
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl
      use neigh   , only : nelst,nelstc,elst,shortelst,nshortelstc
      use polar   , only : polarity, thole, pdamp
      use shunt   , only : cut2
      use utils   , only : set_to_zero1
      use utilgpu , only : maxscaling1, def_queue, dir_queue,real3,real6
      use domdec  , only : rank, loc
      use polgrp  , only : np11, np12, np13, np14, ip11, ip12, ip13,ip14
      use polpot  , only : u1scale, u2scale, u3scale, u4scale
      use potent  , only : use_polarshortreal
      use tmatxb_inl_subroutines
      implicit none
      real(t_p) ,intent(in)   :: mu(:,:,:)
      real(t_p) ,intent(inout):: efi(:,:,:)

      integer    :: i, j, k, irhs
      integer    :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer    :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      integer,pointer :: lst(:,:),nlst(:) ! pointers to neighbor list
      real(t_p)  :: d2                     ! Distance (d2 squared) between atoms
      real(t_p)  :: alsq2, alsq2n          ! ewald precomputed constants
      real(t_p)  :: sdamp,pgamma
      real(t_p)  :: tmp
      type(real6):: dpui,dpuk  ! store coefficients of mu
      type(real3):: posi,dist  ! position of i pole and distance between i and k
      type(real3):: fid,fip    ! Contribution of local pair to efi(iploc)
      type(real3):: fkd,fkp    ! Contribution of local pair to efi(iploc)

      integer   nnp11,nnp12,nnp13,nnp14,ki
      real(t_p) uscale ! Scaling factor for interaction
      integer   :: iscal(maxscaling1) ! Temporary array to store interaction type for neighbors
      real(t_p) :: fscal(maxscaling1) ! Temporary array to store scaling factors for neighbors

      if(deb_Path)write(*,'(4x,a)') 'tmatxb_pme_core1'

      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)

      if (use_polarshortreal) then
          lst =>  shortelst
         nlst => nshortelstc
      else
          lst =>  elst
         nlst => nelstc
      end if
c
!$acc parallel loop gang
!$acc&         vector_length(32)
!$acc&         present(mu,efi,thole,pdamp,polarity,
!$acc&   poleglobnl,ipole,poleloc,x,y,z)
!$acc&         present(nlst,lst)
!$acc&         private(ki,iscal,fscal,dpui,posi)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         !skip atom if it is not polarizable
         !if (polarity(iipole) == 0) cycle MAINLOOP
         iglob  = ipole  (iipole)
         iploc  = poleloc(iipole)
         !if (loc(iglob) == 0) cycle MAINLOOP

         nnp11   = np11(iglob)
         nnp12   = np12(iglob) + nnp11
         nnp13   = np13(iglob) + nnp12
         nnp14   = np14(iglob) + nnp13
         posi%x  = x(iglob)
         posi%y  = y(iglob)
         posi%z  = z(iglob)
         dpui%x  = mu(1,1,iploc)
         dpui%y  = mu(2,1,iploc)
         dpui%z  = mu(3,1,iploc)
         dpui%xx = mu(1,2,iploc)
         dpui%yy = mu(2,2,iploc)
         dpui%zz = mu(3,2,iploc)
         if (nnp14.gt.maxscaling1)
     &      print*,'scaling array to short in tmatxb',ii

         !fill scaling factor along kglob and interaction type
         ki = 0
!$acc loop vector
         do k = 1,nnp12
            if      (k.le.nnp11) then
               iscal(k) = ip11(k,iglob)
               fscal(k) = u1scale
            else
               iscal(k) = ip12(k-nnp11,iglob)
               fscal(k) = u2scale
            end if
         end do
!$acc loop vector
         do k = nnp12+1,nnp14
            if      (k.le.nnp13) then
               iscal(k) = ip13(k-nnp12,iglob)
               fscal(k) = u3scale
            else
               iscal(k) = ip14(k-nnp13,iglob)
               fscal(k) = u4scale
            end if
         end do

         ! Loop on neighbors
!$acc loop vector private(dist,dpuk,fid,fip,fkp,fkd)
         do k =  1, nlst(ii)
            kpole    = lst(k,ii)
            kpoleloc = poleloc(kpole)
            if (kpoleloc.eq.0) cycle
            kglob  = ipole(kpole)

            dist%x = x(kglob) - posi%x
            dist%y = y(kglob) - posi%y
            dist%z = z(kglob) - posi%z
            call image_inl(dist%x,dist%y,dist%z)
            d2 = dist%x*dist%x + dist%y*dist%y + dist%z*dist%z
            if (d2 > cut2) cycle

            sdamp   = pdamp(iipole) * pdamp(kpole)
            pgamma  = min(thole(iipole),thole(kpole))

            ! find exclusion coefficients for connected atoms
            uscale  = 1
            if (ki<nnp14) then
!$acc   loop seq
               do j=1,nnp14
                  if (iscal(j)==kglob) then
                     uscale  = fscal(j)
!$acc   atomic update
                     ki  = ki+1
                     exit
                  end if
               end do
            end if

            dpuk%x   = mu(1,1,kpoleloc)
            dpuk%y   = mu(2,1,kpoleloc)
            dpuk%z   = mu(3,1,kpoleloc)
            dpuk%xx  = mu(1,2,kpoleloc)
            dpuk%yy  = mu(2,2,kpoleloc)
            dpuk%zz  = mu(3,2,kpoleloc)

            call tmatxb_couple(d2,dist,dpui,dpuk,
     &           sdamp,pgamma,aewald,alsq2,alsq2n,uscale,
     &                         fid,fip,fkd,fkp,.false.)


            ! increment electric field for each atoms
!$acc atomic update
            efi(1,1,iploc)    = efi(1,1,iploc)    + fid%x
!$acc atomic update
            efi(2,1,iploc)    = efi(2,1,iploc)    + fid%y
!$acc atomic update
            efi(3,1,iploc)    = efi(3,1,iploc)    + fid%z
!$acc atomic update
            efi(1,2,iploc)    = efi(1,2,iploc)    + fip%x
!$acc atomic update
            efi(2,2,iploc)    = efi(2,2,iploc)    + fip%y
!$acc atomic update
            efi(3,2,iploc)    = efi(3,2,iploc)    + fip%z
!$acc atomic update
            efi(1,1,kpoleloc) = efi(1,1,kpoleloc) + fkd%x
!$acc atomic update
            efi(2,1,kpoleloc) = efi(2,1,kpoleloc) + fkd%y
!$acc atomic update
            efi(3,1,kpoleloc) = efi(3,1,kpoleloc) + fkd%z
!$acc atomic update
            efi(1,2,kpoleloc) = efi(1,2,kpoleloc) + fkp%x
!$acc atomic update
            efi(2,2,kpoleloc) = efi(2,2,kpoleloc) + fkp%y
!$acc atomic update
            efi(3,2,kpoleloc) = efi(3,2,kpoleloc) + fkp%z
         enddo

      end do MAINLOOP

      end

      subroutine tmatxb_pme_core2(mu,efi)
      use atoms   , only : x, y, z
      use atmlst  , only : poleglobnl
      use domdec  , only : rank, loc
      !use erf_mod
      use ewald   , only : aewald
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl
      use neigh   , only : nelstc,elst,shortelst,nshortelstc
      use inform  , only : deb_Path
      use interfaces,only: tmatxb_correct_interactions
      use polar   , only : polarity, thole, pdamp
      use shunt   , only : cut2
      use utilgpu , only : def_queue, dir_queue,real3,real6
      !use polpot  , only : n_uscale,ucorrect_ik,ucorrect_scale
      use potent  , only : use_polarshortreal
      use tinheader, only : ti_p
      use tmatxb_inl_subroutines
      implicit none
      real(t_p) ,intent(in)   :: mu(:,:,:)
      real(t_p) ,intent(inout):: efi(:,:,:)

      integer    :: i, j, k, irhs
      integer    :: start,finla,sized        ! Indexes associated to Matvec split
      integer    :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer    :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      integer,pointer :: lst(:,:),nlst(:)    ! pointers to neighbor list
      real(t_p)  :: d2                     ! Distance (d2 squared) between atoms
      real(t_p)  :: alsq2, alsq2n          ! ewald precomputed constants
      real(t_p)  :: sdamp,pgamma
      real(t_p)  :: tmp
      type(real6):: dpui,dpuk  ! store coefficients of mu
      type(real3):: posi,dist  ! position of i pole and distance between i and k
      type(real3):: fid,fip    ! Contribution of local pair to efi(iploc)
      type(real3):: fkd,fkp    ! Contribution of local pair to efi(iploc)
      !real(t_p) uscale ! Scaling factor for interaction

      if(deb_Path)write(*,'(4x,a)') 'tmatxb_pme_core2'
      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)

      ! Split matvec kernel to ease recovering process in MPI
      if (ndec.eq.1) then
         start  = 1
         finla  = npolelocnl
      else
         sized  = npolelocnl/(ndec+1)
         start  = (icall-1)*sized + 1
         if (icall.eq.ndec) then
            finla = npolelocnl
         else
            finla = icall*sized
         end if
      end if

      if (use_polarshortreal) then
          lst =>  shortelst
         nlst => nshortelstc
      else
          lst =>  elst
         nlst => nelstc
      end if
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(mu,efi,polarity,ipole,poleloc,
!$acc&  loc,x,y,z,pdamp,thole,lst,nlst)
!$acc&         private(posi,dpui)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = start, finla
         iipole = poleglobnl(ii)
         !skip atom if it is not polarizable
         !if (polarity(iipole) == 0) cycle MAINLOOP
         iglob  = ipole  (iipole)
         iploc  = poleloc(iipole)
         !if (loc(iglob) == 0) cycle MAINLOOP

         posi%x  = x(iglob)
         posi%y  = y(iglob)
         posi%z  = z(iglob)
         dpui%x  = mu(1,1,iploc)
         dpui%y  = mu(2,1,iploc)
         dpui%z  = mu(3,1,iploc)
         dpui%xx = mu(1,2,iploc)
         dpui%yy = mu(2,2,iploc)
         dpui%zz = mu(3,2,iploc)

         ! Loop on neighbors
!$acc loop vector private(dist,dpuk,fid,fip,fkd,fkp)
         do k =  1, nlst(ii)
            kpole    = lst(k,ii)
            kpoleloc = poleloc(kpole)
            if (kpoleloc.eq.0) cycle
            kglob  = ipole(kpole)

            dist%x = x(kglob) - posi%x
            dist%y = y(kglob) - posi%y
            dist%z = z(kglob) - posi%z
            call image_inl(dist%x,dist%y,dist%z)
            d2 = dist%x*dist%x + dist%y*dist%y + dist%z*dist%z
            if (d2 > cut2) cycle

            sdamp   = pdamp(iipole) * pdamp(kpole)
            pgamma  = min(thole(iipole),thole(kpole))

            dpuk%x   = mu(1,1,kpoleloc)
            dpuk%y   = mu(2,1,kpoleloc)
            dpuk%z   = mu(3,1,kpoleloc)
            dpuk%xx  = mu(1,2,kpoleloc)
            dpuk%yy  = mu(2,2,kpoleloc)
            dpuk%zz  = mu(3,2,kpoleloc)

            call tmatxb_couple(d2,dist,dpui,dpuk,
     &           sdamp,pgamma,aewald,alsq2,alsq2n,1.0_ti_p,
     &                         fid,fip,fkd,fkp,.false.)


            ! increment electric field for each atoms
!$acc atomic update
            efi(1,1,iploc)    = efi(1,1,iploc)    + fid%x
!$acc atomic update
            efi(2,1,iploc)    = efi(2,1,iploc)    + fid%y
!$acc atomic update
            efi(3,1,iploc)    = efi(3,1,iploc)    + fid%z
!$acc atomic update
            efi(1,2,iploc)    = efi(1,2,iploc)    + fip%x
!$acc atomic update
            efi(2,2,iploc)    = efi(2,2,iploc)    + fip%y
!$acc atomic update
            efi(3,2,iploc)    = efi(3,2,iploc)    + fip%z
!$acc atomic update
            efi(1,1,kpoleloc) = efi(1,1,kpoleloc) + fkd%x
!$acc atomic update
            efi(2,1,kpoleloc) = efi(2,1,kpoleloc) + fkd%y
!$acc atomic update
            efi(3,1,kpoleloc) = efi(3,1,kpoleloc) + fkd%z
!$acc atomic update
            efi(1,2,kpoleloc) = efi(1,2,kpoleloc) + fkp%x
!$acc atomic update
            efi(2,2,kpoleloc) = efi(2,2,kpoleloc) + fkp%y
!$acc atomic update
            efi(3,2,kpoleloc) = efi(3,2,kpoleloc) + fkp%z
         enddo

      end do MAINLOOP

      call tmatxb_correct_interactions(mu,efi)

      end

      subroutine otf_dc_tmatxb_pme_core2(mu,efi)
      use atoms   , only : x, y, z
      use atmlst  , only : poleglobnl
      use domdec  , only : rank, loc
      use divcon  , only : grplst
      !use erf_mod
      use ewald   , only : aewald
      use inform  , only : deb_Path
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl
      use neigh   , only : nelstc,elst,shortelst,nshortelst
      use interfaces,only: otf_dc_tmatxb_correct_interactions
      use polar   , only : polarity, thole, pdamp
      use shunt   , only : cut2
      use utilgpu , only : def_queue, dir_queue,real3,real6
      !use polpot  , only : n_uscale,ucorrect_ik,ucorrect_scale
      use potent  , only : use_polarshortreal
      use tinheader, only : ti_p
      use tmatxb_inl_subroutines
      implicit none
      real(t_p) ,intent(in)   :: mu(:,:,:)
      real(t_p) ,intent(inout):: efi(:,:,:)

      integer    :: i, j, k, irhs, l
      integer    :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer    :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      integer,pointer :: lst(:,:),nlst(:)    ! pointers to neighbor list
      real(t_p)  :: d2                     ! Distance (d2 squared) between atoms
      real(t_p)  :: alsq2, alsq2n          ! ewald precomputed constants
      real(t_p)  :: sdamp,pgamma
      real(t_p)  :: tmp
      type(real6):: dpui,dpuk  ! store coefficients of mu
      type(real3):: posi,dist  ! position of i pole and distance between i and k
      type(real3):: fid,fip    ! Contribution of local pair to efi(iploc)
      type(real3):: fkd,fkp    ! Contribution of local pair to efi(iploc)
      !real(t_p) uscale ! Scaling factor for interaction


      if(deb_Path)
     &   write(*,'(5x,a)') 'otf_dc_tmatxb_pme_core2'
      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)

      if (use_polarshortreal) then
          lst =>  shortelst
         nlst => nshortelst
      else
          lst =>  elst
         nlst => nelstc
      end if
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(mu,efi,polarity,ipole,poleloc,
!$acc&  loc,x,y,z,pdamp,thole,lst,nlst,grplst)
!$acc&         private(posi,dpui)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole  = poleglobnl(ii)
         !skip atom if it is not polarizable
         !if (polarity(iipole) == 0) cycle MAINLOOP
         iglob   = ipole  (iipole)
         iploc   = poleloc(iipole)
         l       = grplst(iglob)

         posi%x  = x(iglob)
         posi%y  = y(iglob)
         posi%z  = z(iglob)
         dpui%x  = mu(1,1,iploc)
         dpui%y  = mu(2,1,iploc)
         dpui%z  = mu(3,1,iploc)
         dpui%xx = mu(1,2,iploc)
         dpui%yy = mu(2,2,iploc)
         dpui%zz = mu(3,2,iploc)

         ! Loop on neighbors
!$acc loop vector private(dist,dpuk,fid,fip,fkd,fkp)
         do k =  1, nlst(ii)
            kpole    = lst(k,ii)
            kpoleloc = poleloc(kpole)
            kglob    = ipole(kpole)
c
c     only build induced field for interactons outside of a block
c     check if the atoms are in the same block
c
            if (.not.(l.ne.grplst(kglob) .or. l.lt.0)) cycle

            dist%x = x(kglob) - posi%x
            dist%y = y(kglob) - posi%y
            dist%z = z(kglob) - posi%z
            call image_inl(dist%x,dist%y,dist%z)
            d2 = dist%x*dist%x + dist%y*dist%y + dist%z*dist%z
            if (d2 > cut2) cycle

            sdamp   = pdamp(iipole) * pdamp(kpole)
            pgamma  = min(thole(iipole),thole(kpole))

            dpuk%x   = mu(1,1,kpoleloc)
            dpuk%y   = mu(2,1,kpoleloc)
            dpuk%z   = mu(3,1,kpoleloc)
            dpuk%xx  = mu(1,2,kpoleloc)
            dpuk%yy  = mu(2,2,kpoleloc)
            dpuk%zz  = mu(3,2,kpoleloc)

            call tmatxb_couple(d2,dist,dpui,dpuk,
     &           sdamp,pgamma,aewald,alsq2,alsq2n,1.0_ti_p,
     &                         fid,fip,fkd,fkp,.false.)


            ! increment electric field for each atoms
!$acc atomic update
            efi(1,1,iploc)    = efi(1,1,iploc)    + fid%x
!$acc atomic update
            efi(2,1,iploc)    = efi(2,1,iploc)    + fid%y
!$acc atomic update
            efi(3,1,iploc)    = efi(3,1,iploc)    + fid%z
!$acc atomic update
            efi(1,2,iploc)    = efi(1,2,iploc)    + fip%x
!$acc atomic update
            efi(2,2,iploc)    = efi(2,2,iploc)    + fip%y
!$acc atomic update
            efi(3,2,iploc)    = efi(3,2,iploc)    + fip%z
!$acc atomic update
            efi(1,1,kpoleloc) = efi(1,1,kpoleloc) + fkd%x
!$acc atomic update
            efi(2,1,kpoleloc) = efi(2,1,kpoleloc) + fkd%y
!$acc atomic update
            efi(3,1,kpoleloc) = efi(3,1,kpoleloc) + fkd%z
!$acc atomic update
            efi(1,2,kpoleloc) = efi(1,2,kpoleloc) + fkp%x
!$acc atomic update
            efi(2,2,kpoleloc) = efi(2,2,kpoleloc) + fkp%y
!$acc atomic update
            efi(3,2,kpoleloc) = efi(3,2,kpoleloc) + fkp%z
         enddo

      end do MAINLOOP

      call otf_dc_tmatxb_correct_interactions(mu,efi)

      end

      subroutine tmatxb_pme_core3(mu,efi)
      use atoms   , only : n
      use atmlst  , only : poleglobnl
      use cell
      use domdec  , only : xbegproc,ybegproc,zbegproc, rank, loc
     &            ,nproc,rank,xendproc,yendproc,zendproc
      !use erf_mod
      use ewald   , only : aewald
      use inform  , only : deb_Path
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl,npolebloc
     &            , npolelocnlb,npolelocnlb_pair,npolelocnlb2_pair
     &            , nspnlb2=>nshortpolelocnlb2_pair
      use neigh   , only : ipole_s=>celle_glob,pglob_s=>celle_pole
     &            , ploc_s=>celle_ploc, ieblst_s=>ieblst
     &            , iseblst_s=>ishorteblst, seblst_s=>shorteblst
     &            , eblst_s=>eblst, x_s=>celle_x, y_s=>celle_y
     &            , z_s=>celle_z
      use interfaces,only: tmatxb_correct_interactions
      use polar   , only : polarity, thole, pdamp
      use potent  , only : use_polarshortreal
      use shunt   , only : cut2
      use tinheader ,only: ti_p
      use tmatxb_inl_subroutines ,only : ndec,icall
      use utilcomm, only : no_commdir
      use utilgpu , only : def_queue, dir_queue,real3,real6
     &            , maxBlock, BLOCK_SIZE
#ifdef _OPENACC
      use cudafor
      use interfaces, only : cu_tmatxb_pme
      use utilcu    , only : BLOCK_DIM,check_launch_kernel
      use utilgpu   , only : def_stream,rec_stream, nSMP
      use tmatxb_pmecu,only: tmatxb_pme_core_cu
#endif
      implicit none
      integer i
      integer start,start1,sized
      real(t_p) ,intent(in)   ::  mu(:,:,:)
      real(t_p) ,intent(inout):: efi(:,:,:)
      integer ierrSync, begin
      integer,save :: gS
      logical,save :: first_in=.true.,dyn_gS=.true.
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

      if(deb_Path) write(*,'(4x,a)') 'tmatxb_pme_core3'
      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)
      begin  = 2*npolelocnlb_pair+1

      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)
#ifdef _CUDA
      if (first_in) then
         ! Compute though occupancy the right gridSize to launch the kernel with
         first_in = .false.
         call cudaMaxGridSize("tmatxb_pme_core_cu",gS)
         if ( gS.eq.0 ) dyn_gS = .true.
         !gS = gS-nSMP
         if (deb_Path) print*, 'tmatxb blockSize ',gS
      end if

      if (use_polarshortreal) then

      if (dyn_gS) gS = min(nspnlb2/4,maxBlock)
!$acc host_data use_device(ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s
!$acc&    ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi)

c     call cu_tmatxb_pme    !Find his interface inside MOD_inteface
c    &     (ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s(begin)
c    &     ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi
c    &     ,npolelocnlb,nspnlb2,npolebloc,n,nproc
c    &     ,.not.no_commdir
c    &     ,cut2,alsq2,alsq2n,aewald
c    &     , xcell, ycell, zcell,xcell2,ycell2,zcell2
c    &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,def_stream)
      call tmatxb_pme_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s(begin)
     &     ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi
     &     ,npolelocnlb,nspnlb2,npolebloc,n
     &     ,cut2,alsq2,alsq2n,aewald
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
      call check_launch_kernel(" tmatxb_pme_core_cu")

!$acc end host_data

      else

      ! Split matvec kernel to ease recovering process in MPI
      sized  = npolelocnlb2_pair/ndec
      start  = (icall-1)*sized + 1
      start1 = begin+(start-1)*BLOCK_SIZE
      if (icall.eq.ndec) sized = npolelocnlb2_pair-start+1
      if (dyn_gS) gS = min(sized/4,maxBlock)

!$acc host_data use_device(ipole_s,pglob_s,ploc_s,ieblst_s,eblst_s
!$acc&    ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi)

c     call cu_tmatxb_pme    !Find his interface inside MOD_inteface
c    &     (ipole_s,pglob_s,ploc_s
c    &     ,ieblst_s(start),eblst_s(start1)
c    &     ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi
c    &     ,npolelocnlb,sized,npolebloc,n,nproc
c    &     ,.not.no_commdir
c    &     ,cut2,alsq2,alsq2n,aewald
c    &     , xcell, ycell, zcell,xcell2,ycell2,zcell2
c    &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,def_stream)
      call tmatxb_pme_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,ploc_s
     &     ,ieblst_s(start),eblst_s(start1)
     &     ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi
     &     ,npolelocnlb,sized,npolebloc,n
     &     ,cut2,alsq2,alsq2n,aewald
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
      call check_launch_kernel(" tmatxb_pme_core_cu")

!$acc end host_data

      end if
#else
      print 100
 100  format('tmatxb_pme_core3 is a specific device routine',/,
     &       'you are not supposed to get inside with your compile',
     &        'type.')
      call fatal
#endif
c
      call tmatxb_correct_interactions(mu,efi)

      end

      subroutine tmatxb_pme_core_v4(mu,efi)
      use atoms   , only : n
      use atmlst  , only : poleglobnl
      use cell
      use domdec  , only : xbegproc,ybegproc,zbegproc,rank,loc
     &            , nproc,rank,xendproc,yendproc,zendproc
      !use erf_mod
      use ewald   , only : aewald
      use inform  , only : deb_Path
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl,npolebloc
     &            , npolelocnlb,npolelocnlb_pair,npolelocnlb2_pair
     &            , nspnlb2=>nshortpolelocnlb2_pair
      use neigh   , only : ipole_s=>celle_glob,pglob_s=>celle_pole
     &            , ploc_s=>celle_ploc, ieblst_s=>ieblst
     &            , iseblst_s=>ishorteblst, seblst_s=>shorteblst
     &            , eblst_s=>eblst, x_s=>celle_x, y_s=>celle_y
     &            , z_s=>celle_z
      use interfaces,only: tmatxb_correct_interactions
      use polar   , only : polarity, thole, pdamp
      use potent  , only : use_polarshortreal
      use shunt   , only : cut2
      use tinheader ,only: ti_p
      use tinMemory ,only: prmem_request
      use tmatxb_inl_subroutines ,only : ndec,icall,mu3,mu4,efi3,efi4
      use utilcomm, only : no_commdir
      use utilgpu , only : def_queue, dir_queue,real3,real6
     &            , BLOCK_SIZE
#ifdef _OPENACC
      use cudafor
      use interfaces, only : cu_tmatxb_pme
      use utilcu    , only : BLOCK_DIM,check_launch_kernel
      use utilgpu   , only : def_stream,rec_stream, nSMP
      use tmatxb_pmecu,only: tmatxb_pme_core_v4_cu
#endif
      implicit none
      integer i
      integer start,start1,sized
      real(t_p) ,intent(in)   ::  mu(:,:,:)
      real(t_p) ,intent(inout):: efi(:,:,:)
      integer ierrSync, begin
      integer,save :: gS
      logical,save :: first_in=.true.,dyn_gS=.true.
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

      if(deb_Path) write(*,'(4x,a)') 'tmatxb_pme_core_v4'
      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)
      begin  = 2*npolelocnlb_pair+1

      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)
#ifdef _CUDA
      if (first_in) then
         ! Compute though occupancy the right gridSize to launch the kernel with
         first_in = .false.
         call cudaMaxGridSize("tmatxb_pme_core_cu",gS)
         if ( gS.eq.0 ) dyn_gS = .true.
         !gS = gS-nSMP
         if (deb_Path) print*, 'tmatxb blockSize ',gS
      end if

      call prmem_request(mu3 ,3,npolebloc)
      call prmem_request(mu4 ,3,npolebloc)
      call prmem_request(efi3,3,npolebloc)
      call prmem_request(efi4,3,npolebloc)

      if (use_polarshortreal) then

      if (dyn_gS) gS = nspnlb2/4
!$acc host_data use_device(ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s
!$acc&    ,x_s,y_s,z_s,pdamp,thole,polarity,mu,mu3,mu4,efi,efi3,efi4)

c     call cu_tmatxb_pme    !Find his interface inside MOD_inteface
c    &     (ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s(begin)
c    &     ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi
c    &     ,npolelocnlb,nspnlb2,npolebloc,n,nproc
c    &     ,.not.no_commdir
c    &     ,cut2,alsq2,alsq2n,aewald
c    &     , xcell, ycell, zcell,xcell2,ycell2,zcell2
c    &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,def_stream)
      call tmatxb_pme_core_v4_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s(begin)
     &     ,x_s,y_s,z_s,pdamp,thole,polarity
     &     ,mu,mu3,mu4,efi,efi3,efi4
     &     ,npolelocnlb,nspnlb2,npolebloc,n
     &     ,cut2,alsq2,alsq2n,aewald
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
      call check_launch_kernel(" tmatxb_pme_core_cu")

!$acc end host_data

      else

      if (dyn_gS) gS = npolelocnlb2_pair/4
      ! Split matvec kernel to ease recovering process in MPI
      sized  = npolelocnlb2_pair/ndec
      start  = (icall-1)*sized + 1
      start1 = begin+(start-1)*BLOCK_SIZE
      if (icall.eq.ndec) sized = npolelocnlb2_pair-start+1

!$acc host_data use_device(ipole_s,pglob_s,ploc_s,ieblst_s,eblst_s
!$acc&    ,x_s,y_s,z_s,pdamp,thole,polarity,mu,mu3,mu4,efi,efi3,efi4)

c     call cu_tmatxb_pme    !Find his interface inside MOD_inteface
c    &     (ipole_s,pglob_s,ploc_s
c    &     ,ieblst_s(start),eblst_s(start1)
c    &     ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi
c    &     ,npolelocnlb,sized,npolebloc,n,nproc
c    &     ,.not.no_commdir
c    &     ,cut2,alsq2,alsq2n,aewald
c    &     , xcell, ycell, zcell,xcell2,ycell2,zcell2
c    &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,def_stream)
      call tmatxb_pme_core_v4_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,ploc_s
     &     ,ieblst_s(start),eblst_s(start1)
     &     ,x_s,y_s,z_s,pdamp,thole,polarity
     &     ,mu,mu3,mu4,efi,efi3,efi4
     &     ,npolelocnlb,sized,npolebloc,n
     &     ,cut2,alsq2,alsq2n,aewald
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
      call check_launch_kernel(" tmatxb_pme_core_cu")

!$acc end host_data

      end if
#else
      print 100
 100  format('tmatxb_pme_core3 is a specific device routine',/,
     &       'you are not supposed to get inside with your compile',
     &        'type.')
      call fatal
#endif
c
      call tmatxb_correct_interactions(mu,efi)

      end

      subroutine otf_dc_tmatxb_pme_core3(mu,efi)
      use atoms   , only : n
      use atmlst  , only : poleglobnl
      use cell
      use divcon  , only : grplst
      use domdec  , only : xbegproc,ybegproc,zbegproc, rank, loc
     &            ,nproc,rank,xendproc,yendproc,zendproc
      !use erf_mod
      use ewald   , only : aewald
      use inform  , only : deb_Path
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl,npolebloc
     &            , npolelocnlb,npolelocnlb_pair,npolelocnlb2_pair
     &            , nspnlb2=>nshortpolelocnlb2_pair
      use neigh   , only : ipole_s=>celle_glob,pglob_s=>celle_pole
     &            , ploc_s=>celle_ploc, ieblst_s=>ieblst
     &            , iseblst_s=>ishorteblst, seblst_s=>shorteblst
     &            , eblst_s=>eblst, x_s=>celle_x, y_s=>celle_y
     &            , z_s=>celle_z
      use interfaces,only: otf_dc_tmatxb_correct_interactions
      use polar   , only : polarity, thole, pdamp
      use potent  , only : use_polarshortreal
      use shunt   , only : cut2
      use tinheader ,only: ti_p
      use utilgpu , only : def_queue, dir_queue,real3,real6
#ifdef _OPENACC
      use cudafor
      use interfaces, only : cu_tmatxb_pme
      use utilcu    , only : BLOCK_DIM,check_launch_kernel
      use utilgpu   , only : def_stream, nSMP
      use tmatxb_pmecu,only: otfdc_tmatxb_pme_core_cu
#endif
      implicit none
      real(t_p) ,intent(in)   :: mu(:,:,:)
      real(t_p) ,intent(inout):: efi(:,:,:)
      integer ierrSync, begin
      integer,save :: gS
      logical,save :: first_in=.true.,dyn_gS=.false.
      real(t_p) alsq2,alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend

      if(deb_Path)
     &   write(*,'(4x,a)') 'oft_dc_tmatxb_pme_core3'
      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)
      begin  = 2*npolelocnlb_pair+1

      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)
#ifdef _CUDA
      if (first_in) then
         ! Compute though occupancy the right gridSize to launch the kernel with
         first_in = .false.
         call cudaMaxGridSize("otf_dc_tmatxb_pme_core_cu",gS)
         if ( gS.eq.0 ) dyn_gS = .true.
         !gS = gS-nSMP
         if (deb_Path) 
     &      print*, 'otf tmatxb blockSize ',gS
      end if

      if (use_polarshortreal) then

!$acc host_data use_device(ipole_s,pglob_s,ploc_s,grplst
!$acc&    ,iseblst_s,seblst_s,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi)

      if (dyn_gS) gS = nspnlb2/4
      call otfdc_tmatxb_pme_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,ploc_s,grplst,iseblst_s,seblst_s(begin)
     &     ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi
     &     ,npolelocnlb,nspnlb2,npolebloc,n
     &     ,cut2,alsq2,alsq2n,aewald
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
      call check_launch_kernel(" otf_tmatxb_pme_core_cu")

!$acc end host_data

      else

!$acc host_data use_device(ipole_s,pglob_s,ploc_s,grplst
!$acc&    ,ieblst_s,eblst_s,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi)

      if (dyn_gS) gS = npolelocnlb2_pair/4
      call otfdc_tmatxb_pme_core_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &     (ipole_s,pglob_s,ploc_s,grplst,ieblst_s,eblst_s(begin)
     &     ,x_s,y_s,z_s,pdamp,thole,polarity,mu,efi
     &     ,npolelocnlb,npolelocnlb2_pair,npolebloc,n
     &     ,cut2,alsq2,alsq2n,aewald
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend)
      call check_launch_kernel(" otf_dc_tmatxb_pme_core_cu")

!$acc end host_data

      end if
#else
      print 100
 100  format('otf_dc_tmatxb_pme_core3 is a specific device routine',/,
     &       'you are not supposed to get inside with your compile',
     &        'type.')
      call fatal
#endif
c
      call otf_dc_tmatxb_correct_interactions(mu,efi)

      end



c======================================================================
c       --------------------------------------------------------
c       Correct scaling pairwise interactions routine for tmatxb
c       --------------------------------------------------------
c======================================================================

      subroutine tmatxb_correct_interactions(mu,efi)
      use atoms   , only : x, y, z
      use atmlst  , only : poleglobnl
      !use erf_mod
      use ewald   , only : aewald
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl,npolebloc
      use polar   , only : polarity, thole, pdamp
      use shunt   , only : cut2
      use utilgpu , only : def_queue, dir_queue,real3,real6
      use utilcomm, only : no_commdir
      use domdec  , only : loc
      use polpot  , only : n_uscale,ucorrect_ik,ucorrect_scale
      use tmatxb_inl_subroutines
      implicit none
      real(t_p) ,intent(in) :: mu (:,:,:)
      real(t_p) ,intent(out):: efi(:,:,:)

      integer    :: i, j, k, irhs
      integer    :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer    :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      real(t_p)  :: d2                     ! Distance (d2 squared) between atoms
      real(t_p)  :: alsq2, alsq2n          ! ewald precomputed constants
      real(t_p)  :: sdamp,pgamma
      type(real6):: dpui,dpuk  ! store coefficients of mu
      type(real3):: posi,dist  ! position of i pole and distance between i and k
      type(real3):: fid,fip    ! Contribution of local pair to efi(iploc)
      type(real3):: fkd,fkp    ! Contribution of local pair to efi(iploc)

      real(t_p) uscale ! Scaling factor for interaction

      if (n_uscale.eq.0) return

      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)

      ! Scaling corrections loop
!$acc parallel loop default(present)
!$acc&         private(dpui,dpuk,dist,fid,fip,fkd,fkp)
!$acc&         async(def_queue)
      do ii = 1,n_uscale
         iipole   = ucorrect_ik(2*(ii-1)+1)
         kpole    = ucorrect_ik(2*(ii-1)+2)
         uscale   = ucorrect_scale(ii)

         iglob    = ipole  (iipole)
         kglob    = ipole  (kpole)
         iploc    = poleloc(iipole)
         kpoleloc = poleloc(kpole)
         !skip atom if it is not polarizable
         !FIXME do we need to test on polarity
         if (polarity(iipole) == 0) cycle
         if (iploc.eq.0.or.iploc.gt.npolebloc) cycle
         if (kpoleloc.eq.0.or.kpoleloc.gt.npolebloc) cycle

         dpui%x   = mu(1,1,iploc)
         dpui%y   = mu(2,1,iploc)
         dpui%z   = mu(3,1,iploc)
         dpui%xx  = mu(1,2,iploc)
         dpui%yy  = mu(2,2,iploc)
         dpui%zz  = mu(3,2,iploc)

         dpuk%x   = mu(1,1,kpoleloc)
         dpuk%y   = mu(2,1,kpoleloc)
         dpuk%z   = mu(3,1,kpoleloc)
         dpuk%xx  = mu(1,2,kpoleloc)
         dpuk%yy  = mu(2,2,kpoleloc)
         dpuk%zz  = mu(3,2,kpoleloc)

         dist%x   = x(kglob) - x(iglob)
         dist%y   = y(kglob) - y(iglob)
         dist%z   = z(kglob) - z(iglob)
         call image_inl(dist%x,dist%y,dist%z)
         d2  = dist%x*dist%x + dist%y*dist%y + dist%z*dist%z
         if (d2 > cut2) cycle

         sdamp    = pdamp(iipole) * pdamp(kpole)
         pgamma   = min(thole(iipole),thole(kpole))

         call tmatxb_couple(d2,dist,dpui,dpuk,
     &        sdamp,pgamma,aewald,alsq2,alsq2n,uscale,
     &                      fid,fip,fkd,fkp,.true.)

         ! increment electric field for each atoms
!$acc atomic update
         efi(1,1,iploc)    = efi(1,1,iploc)    + fid%x
!$acc atomic update
         efi(2,1,iploc)    = efi(2,1,iploc)    + fid%y
!$acc atomic update
         efi(3,1,iploc)    = efi(3,1,iploc)    + fid%z
!$acc atomic update
         efi(1,2,iploc)    = efi(1,2,iploc)    + fip%x
!$acc atomic update
         efi(2,2,iploc)    = efi(2,2,iploc)    + fip%y
!$acc atomic update
         efi(3,2,iploc)    = efi(3,2,iploc)    + fip%z
!$acc atomic update
         efi(1,1,kpoleloc) = efi(1,1,kpoleloc) + fkd%x
!$acc atomic update
         efi(2,1,kpoleloc) = efi(2,1,kpoleloc) + fkd%y
!$acc atomic update
         efi(3,1,kpoleloc) = efi(3,1,kpoleloc) + fkd%z
!$acc atomic update
         efi(1,2,kpoleloc) = efi(1,2,kpoleloc) + fkp%x
!$acc atomic update
         efi(2,2,kpoleloc) = efi(2,2,kpoleloc) + fkp%y
!$acc atomic update
         efi(3,2,kpoleloc) = efi(3,2,kpoleloc) + fkp%z
      end do
      end

      subroutine otf_dc_tmatxb_correct_interactions(mu,efi)
      use atoms   , only : x, y, z
      use atmlst  , only : poleglobnl
      !use erf_mod
      use divcon  , only : grplst
      use ewald   , only : aewald
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl
      use polar   , only : polarity, thole, pdamp
      use shunt   , only : cut2
      use utilgpu , only : def_queue, dir_queue,real3,real6
      use domdec  , only : loc
      use polpot  , only : n_uscale,ucorrect_ik,ucorrect_scale
      use tmatxb_inl_subroutines
      implicit none
      real(t_p) ,intent(in) :: mu (:,:,:)
      real(t_p) ,intent(out):: efi(:,:,:)

      integer    :: i, j, k, irhs
      integer    :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer    :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      real(t_p)  :: d2                     ! Distance (d2 squared) between atoms
      real(t_p)  :: alsq2, alsq2n          ! ewald precomputed constants
      real(t_p)  :: sdamp,pgamma
      type(real6):: dpui,dpuk  ! store coefficients of mu
      type(real3):: posi,dist  ! position of i pole and distance between i and k
      type(real3):: fid,fip    ! Contribution of local pair to efi(iploc)
      type(real3):: fkd,fkp    ! Contribution of local pair to efi(iploc)

      real(t_p) uscale ! Scaling factor for interaction

      if (n_uscale.eq.0) return

      alsq2  = 2 * aewald**2
      alsq2n = 0
      if (aewald > 0) alsq2n = 1 / (sqrtpi*aewald)

      ! Scaling corrections loop
!$acc parallel loop default(present)
!$acc&         private(dpui,dpuk,dist,fid,fip,fkd,fkp)
!$acc&         async(def_queue)
      do ii = 1,n_uscale
         iipole   = ucorrect_ik(2*(ii-1)+1)
         kpole    = ucorrect_ik(2*(ii-1)+2)
         uscale   = ucorrect_scale(ii)

         iglob    = ipole  (iipole)
         kglob    = ipole  (kpole)
         iploc    = poleloc(iipole)
         kpoleloc = poleloc(kpole)
         if (.not.(grplst(iglob).eq.grplst(kglob).and.
     &             grplst(iglob).ne.-1)) cycle

         dpui%x   = mu(1,1,iploc)
         dpui%y   = mu(2,1,iploc)
         dpui%z   = mu(3,1,iploc)
         dpui%xx  = mu(1,2,iploc)
         dpui%yy  = mu(2,2,iploc)
         dpui%zz  = mu(3,2,iploc)

         dpuk%x   = mu(1,1,kpoleloc)
         dpuk%y   = mu(2,1,kpoleloc)
         dpuk%z   = mu(3,1,kpoleloc)
         dpuk%xx  = mu(1,2,kpoleloc)
         dpuk%yy  = mu(2,2,kpoleloc)
         dpuk%zz  = mu(3,2,kpoleloc)

         dist%x   = x(kglob) - x(iglob)
         dist%y   = y(kglob) - y(iglob)
         dist%z   = z(kglob) - z(iglob)
         call image_inl(dist%x,dist%y,dist%z)
         d2  = dist%x*dist%x + dist%y*dist%y + dist%z*dist%z
         if (d2 > cut2) cycle

         sdamp    = pdamp(iipole) * pdamp(kpole)
         pgamma   = min(thole(iipole),thole(kpole))

         call tmatxb_couple(d2,dist,dpui,dpuk,
     &        sdamp,pgamma,aewald,alsq2,alsq2n,uscale,
     &                      fid,fip,fkd,fkp,.true.)

         ! increment electric field for each atoms
!$acc atomic update
         efi(1,1,iploc)    = efi(1,1,iploc)    + fid%x
!$acc atomic update
         efi(2,1,iploc)    = efi(2,1,iploc)    + fid%y
!$acc atomic update
         efi(3,1,iploc)    = efi(3,1,iploc)    + fid%z
!$acc atomic update
         efi(1,2,iploc)    = efi(1,2,iploc)    + fip%x
!$acc atomic update
         efi(2,2,iploc)    = efi(2,2,iploc)    + fip%y
!$acc atomic update
         efi(3,2,iploc)    = efi(3,2,iploc)    + fip%z
!$acc atomic update
         efi(1,1,kpoleloc) = efi(1,1,kpoleloc) + fkd%x
!$acc atomic update
         efi(2,1,kpoleloc) = efi(2,1,kpoleloc) + fkd%y
!$acc atomic update
         efi(3,1,kpoleloc) = efi(3,1,kpoleloc) + fkd%z
!$acc atomic update
         efi(1,2,kpoleloc) = efi(1,2,kpoleloc) + fkp%x
!$acc atomic update
         efi(2,2,kpoleloc) = efi(2,2,kpoleloc) + fkp%y
!$acc atomic update
         efi(3,2,kpoleloc) = efi(3,2,kpoleloc) + fkp%z
      end do
      end
