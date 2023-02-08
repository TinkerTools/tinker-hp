c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine epolar0cvec  --  Ewald polarization derivs via list  ##
c     ##     o                                                            ##
c     ######################################################################
c
c
c     "epolar0c" calculates the dipole polarization energy with respect
c     to Cartesian coordinates using particle mesh Ewald summation and
c     a neighbor list
c
c
#include "tinker_macro.h"
      subroutine epolar0cvec
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
      use polar
      use polpot
      use potent
      use tinheader ,only:ti_p,re_p
      use mpi

      implicit none
      integer i,ii
      real(t_p) e,f
      real(t_p) term,fterm
      !DIR$ ATTRIBUTES ALIGN:64::iipolevec,iglobvec
      integer,allocatable :: iipolevec(:),iglobvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::divec,uivec,uiivec
      real(t_p),allocatable :: divec(:,:), uivec(:,:),uiivec(:)
      !DIR$ ATTRIBUTES ALIGN:64::dvec,uvec
      real(t_p),allocatable :: dvec(:), uvec(:)

      allocate(iipolevec(npole))
      allocate(iglobvec(npole))
      allocate(divec(npole,3))
      allocate(uivec(npole,3))
      allocate(uiivec(npole))
      allocate(dvec(3))
      allocate(uvec(3))

c
c
c     zero out the polarization energy and derivatives
c
      write(*,*) 'epolar0cvec'
      ep = 0.0_ti_p
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole

c     compute the induced dipoles at each polarizable atom
c
      if (use_pmecore) then
        if (polalg.eq.5) then
          call dcinduce_pme
        else
          call newinduce_pme
        end if
      else
        if (polalg.eq.5) then
          call dcinduce_pme2
        else
          call newinduce_pme2
        end if
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     $  call eprecipvec
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     $   then
        call epreal0cvec
c
c     compute the Ewald self-energy term over all the atoms
c
        term = 2.0_ti_p * aewald * aewald
        fterm = -f * aewald / sqrtpi

        iipolevec(1:npoleloc) = poleglob(1:npoleloc)
        iglobvec (1:npoleloc) = ipole(iipolevec(1:npoleloc))
        divec(1:npoleloc,1:3) =
     &            transpose (
     &                       rpole(2:4,iipolevec(1:npoleloc))
     &                      )
        uivec(1:npoleloc,1:3) =
     &            transpose (
     &                        uind(1:3,iipolevec(1:npoleloc))
     &                      )

        do ii = 1, npoleloc
           uiivec(ii) = dot_product(divec(ii,1:3),uivec(ii,1:3))
        enddo
        e = sum(uiivec(1:npoleloc)) * term * fterm / 3.0_ti_p
        ep = ep + e
c
c       compute the cell dipole boundary correction term
c
        if (boundary .eq. 'VACUUM') then
           dvec(1) = sum(   rpole(2,iipolevec(1:npoleloc))
     &                   +  rpole(1,iipolevec(1:npoleloc))
     &                    * x(iglobvec(1:npoleloc))
     &                  )
           dvec(2) = sum(   rpole(3,iipolevec(1:npoleloc))
     &                   +  rpole(1,iipolevec(1:npoleloc))
     &                    * y(iglobvec(1:npoleloc))
     &                  )
           dvec(3) = sum(   rpole(4,iipolevec(1:npoleloc))
     &                   +  rpole(1,iipolevec(1:npoleloc))
     &                    * z(iglobvec(1:npoleloc))
     &                  )
           uvec(1) = sum( uind(1,iipolevec(1:npoleloc)))
           uvec(2) = sum( uind(2,iipolevec(1:npoleloc)))
           uvec(3) = sum( uind(3,iipolevec(1:npoleloc)))

           term = (2.0_ti_p/3.0_ti_p) * f * (pi/volbox)
           ep = ep + term*dot_product(dvec,uvec)
        end if
      end if
      return
      end
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine epreal0cvec  --  real space polar energy via list  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "epreal0c" calculates the induced dipole polarization energy
c     using particle mesh Ewald summation and a neighbor list
c
c
      subroutine epreal0cvec
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
      use polar
      use polgrp
      use polpot
      use neigh
      use potent
      use shunt
      use tinheader ,only:ti_p,re_p
      use mpi

      implicit none
      integer i,j,k,inl
      integer ii,kk,iii,kkk,iipole,kkpole
      integer iglob,kglob,nnelst,nnelst1
      integer nnp11,nn12,nn13,nn1213,nn14,nn121314,nn15,ntot
      real(t_p) e,f
      real(t_p) alsq2,alsq2n
      real(t_p) pdi,pti
      real(t_p) xi,yi,zi
      real(t_p) ci
      !DIR$ ATTRIBUTES ALIGN:64::pscale
      real(t_p), allocatable :: pscale(:)
      !DIR$ ATTRIBUTES ALIGN:64::pscalevec
      integer, allocatable :: pscalevec(:)
      !DIR$ ATTRIBUTES ALIGN:64::qivec,qritmp,qrktmp
      real(t_p), allocatable :: qivec(:,:),qritmp(:),qrktmp(:)
      !DIR$ ATTRIBUTES ALIGN:64::itmp12,itmp13,itmp14,itmp15
      integer,allocatable :: itmp12(:),itmp13(:),itmp14(:),itmp15(:)
      !DIR$ ATTRIBUTES ALIGN:64:: kglobvec,kbisvec,kkpolevec
      integer, allocatable :: kglobvec(:),kbisvec(:),kkpolevec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: kglobvec1,kbisvec1,kkpolevec1
      integer, allocatable :: kglobvec1(:),kbisvec1(:),kkpolevec1(:)
      !DIR$ ATTRIBUTES ALIGN:64::posvec,posvec1
      real(t_p),allocatable :: posvec(:,:),posvec1(:),
     &                                   posvec2(:,:)
      !DIR$ ATTRIBUTES ALIGN:64::r2vec,r2vec1,rvec
      real(t_p),allocatable :: r2vec(:),r2vec1(:),rvec(:)
      !DIR$ ATTRIBUTES ALIGN:64::divec,uivec
      real(t_p),allocatable :: divec(:),uivec(:)
      !DIR$ ATTRIBUTES ALIGN:64::ckvec,dkvec,qkvec,ukvec
      real(t_p),allocatable :: ckvec(:),dkvec(:,:),qkvec(:,:,:),
     &                                  ukvec(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: mask
      logical,allocatable :: mask(:),mask1(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: drivec,drkvec,qrivec,qrkvec
      real(t_p),allocatable :: drivec(:),drkvec(:),qrivec(:,:),
     &                         qrkvec(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: qrrivec,qrrkvec
      real(t_p),allocatable :: qrrivec(:),qrrkvec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: urivec,urkvec
      real(t_p),allocatable :: urivec(:),urkvec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: duikvec,quikvec
      real(t_p),allocatable :: duikvec(:),quikvec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: rr3vec,rr5vec,rr7vec
      real(t_p),allocatable :: rr3vec(:),rr5vec(:),rr7vec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: ralphavec
      real(t_p),allocatable :: ralphavec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: bnvec
      real(t_p),allocatable :: bnvec(:,:)
      !DIR$ ATTRIBUTES ALIGN:64:: exp2avec
      real(t_p),allocatable :: exp2avec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: dampvec,dampvec1
      real(t_p),allocatable :: dampvec(:),dampvec1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: pgammavec,expdampvec1
      real(t_p),allocatable :: pgammavec(:),expdampvec1(:)
      !DIR$ ATTRIBUTES ALIGN:64:: sc3vec,sc5vec,sc7vec
      real(t_p),allocatable :: sc3vec(:),sc5vec(:),sc7vec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: psc3vec,psc5vec ,psc7vec
      real(t_p),allocatable :: psc3vec(:),psc5vec(:),psc7vec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: psr3vec,psr5vec ,psr7vec
      real(t_p),allocatable :: psr3vec(:),psr5vec(:),psr7vec(:)
      !DIR$ ATTRIBUTES ALIGN:64:: term1vec ,term2vec ,term3vec
      real(t_p),allocatable :: term1vec(:),term2vec(:),
     &                                     term3vec(:)

      character*10 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (qivec(3,3),qritmp(3),qrktmp(3))
      allocate (pscalevec(npole))
      allocate (itmp12(npole))
      allocate (itmp13(npole))
      allocate (itmp14(npole))
      allocate (itmp15(npole))
      allocate (kglobvec(maxelst),kbisvec(maxelst),kkpolevec(maxelst))
      allocate (kglobvec1(maxelst),kbisvec1(maxelst))
      allocate (kkpolevec1(maxelst))
      allocate (posvec(maxelst,3),posvec1(3*maxelst),posvec2(maxelst,3))
      allocate (r2vec(maxelst),r2vec1(maxelst))
      allocate (rvec(maxelst))
      allocate (ckvec(maxelst),dkvec(maxelst,3),qkvec(maxelst,3,3))
      allocate (divec(3),uivec(3))
      allocate (ukvec(maxelst,3))
      allocate (mask(maxelst))
      allocate (mask1(maxelst,3))
      allocate (drivec(maxelst))
      allocate (drkvec(maxelst))
      allocate (qrivec(maxelst,3))
      allocate (qrkvec(maxelst,3))
      allocate (qrrivec(maxelst),qrrkvec(maxelst))
      allocate (urivec(maxelst),urkvec(maxelst))
      allocate (duikvec(maxelst),quikvec(maxelst))
      allocate (rr3vec(maxelst))
      allocate (rr5vec(maxelst))
      allocate (rr7vec(maxelst))
      allocate (ralphavec(maxelst))
      allocate (bnvec(maxelst,0:3))
      allocate (exp2avec(maxelst))
      allocate (dampvec(maxelst),dampvec1(maxelst))
      allocate (pgammavec(maxelst),expdampvec1(maxelst))
      allocate (sc3vec(maxelst))
      allocate (sc5vec(maxelst))
      allocate (sc7vec(maxelst))
      allocate (psc3vec(maxelst))
      allocate (psc5vec(maxelst))
      allocate (psc7vec(maxelst))
      allocate (psr3vec(maxelst))
      allocate (psr5vec(maxelst))
      allocate (psr7vec(maxelst))
      allocate (term1vec(maxelst))
      allocate (term2vec(maxelst))
      allocate (term3vec(maxelst))

c
c     initialize connected atom exclusion coefficients
c
      pscale = 1.0_ti_p
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5_ti_p * electric / dielec
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
         iglob = ipole(iipole)
         i = loc(iglob)
         pdi = pdamp(iipole)
         pti = thole(iipole)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         sc3vec = 1.0_ti_p
         sc5vec = 1.0_ti_p
         sc7vec = 1.0_ti_p

cold     ci = rpole(1,iipole)
         ci = rpole(1,iipole) ! 
cold     dix = rpole(2,iipole)
cold     diy = rpole(3,iipole)
cold     diz = rpole(4,iipole)
cold     qixx = rpole(5,iipole)
cold     qixy = rpole(6,iipole)
cold     qixz = rpole(7,iipole)
cold     qiyy = rpole(9,iipole)
cold     qiyz = rpole(10,iipole)
cold     qizz = rpole(13,iipole)
c        divec is dix, diy, diz
         divec = rpole(2:4,iipole)
c    
c        qivec  is qixx, qixy, qixz
c               qixy, qiyy, qiyz
c               qixz, qiyz, qizz

            !DIR$ ASSUME_ALIGNED qivec:64
         qivec= reshape ( rpole( (/ 5,6,7,6,9,10,7,10,13 /), iipole),
     &                 (/3,3/)
     &               )
cold     uix = uind(1,iipole)
cold     uiy = uind(2,iipole)
cold     uiz = uind(3,iipole)
c        ui is  uix, uiy, uiz
         uivec = uind(1:3,iipole)


cold     do j = 1, n12(iglob)
cold        pscale(i12(j,iglob)) = p2scale
cold     end do
cold     do j = 1, n13(iglob)
cold        pscale(i13(j,iglob)) = p3scale
cold     end do
cold     do j = 1, n14(iglob)
cold        pscale(i14(j,iglob)) = p4scale
cold        do k = 1, np11(iglob)
cold            if (i14(j,iglob) .eq. ip11(k,iglob))
cold &            pscale(i14(j,iglob)) = p4scale * p41scale
cold        end do
cold     end do
cold     do j = 1, n15(iglob)
cold        pscale(i15(j,iglob)) = p5scale
cold     end do
         nnp11 = np11(iglob)
         nn12  = n12(iglob)
         nn13  = n13(iglob)
         nn14  = n14(iglob)
         nn15  = n15(iglob)

         ntot=nn12+nn13+nn14+nn15

         !DIR$ ASSUME_ALIGNED itmp12:64
         !DIR$ ASSUME_ALIGNED itmp13:64
         !DIR$ ASSUME_ALIGNED itmp14:64
         !DIR$ ASSUME_ALIGNED itmp15:64
         itmp12  = i12(1:nn12,iglob)
         itmp13  = i13(1:nn13,iglob)
         itmp14  = i14(1:nn14,iglob)
         itmp15  = i15(1:nn15,iglob)
         !DIR$ ASSUME_ALIGNED itmp12:64
         !DIR$ ASSUME_ALIGNED itmp13:64
         !DIR$ ASSUME_ALIGNED itmp14:64
         !DIR$ ASSUME_ALIGNED itmp15:64
         pscalevec(1:ntot) = (/ itmp12(1:nn12),
     &                          itmp13(1:nn13),
     &                          itmp14(1:nn14),
     &                          itmp15(1:nn15)
     &                       /)
         pscale(pscalevec(1:ntot))= (/ (p2scale,i=1,nn12),
     &                                 (p3scale,i=1,nn13),
     &                                 (p4scale,i=1,nn14),
     &                                 (p5scale,i=1,nn15)
     &                              /)
        do j = 1, nn14
            do k = 1, nnp11
                if (i14(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i14(j,iglob)) = p4scale * p41scale
            end do
         end do

c
c     evaluate all sites within the cutoff distance
c
cold    do kkk = 1, nelst(ii)
cold       kkpole = elst(kkk,ii)
cold       kglob = ipole(kkpole)
cold       xr = x(kglob) - xi
cold       yr = y(kglob) - yi
cold       zr = z(kglob) - zi
cold       if (use_bounds)  call image (xr,yr,zr)
cold       r2 = xr*xr + yr*yr + zr*zr
         nnelst =nelst(ii)

         kkpolevec(1:nnelst)    = elst           (1:nnelst,ii)
         kglobvec (1:nnelst)    = ipole(kkpolevec(1:nnelst))
         posvec   (1:nnelst,1)  = x(kglobvec     (1:nnelst)) -  xi
         posvec   (1:nnelst,2)  = y(kglobvec     (1:nnelst)) -  yi
         posvec   (1:nnelst,3)  = z(kglobvec     (1:nnelst)) -  zi

         if (use_bounds)  call imagevec2 (posvec(1:nnelst,:),nnelst)
         r2vec(1:nnelst) =  posvec(1:nnelst,1)**2
     &                    + posvec(1:nnelst,2)**2
     &                    + posvec(1:nnelst,3)**2

cold        if (r2 .le. off2) then
c
c     evaluate all sites within the cutoff distance
c
         mask(1:nnelst)       = (r2vec(1:nnelst)<=off2)    ! mask for r2
         nnelst1              = count(mask(1:nnelst))
         kkpolevec1(1:nnelst1)= pack(kkpolevec(1:nnelst),mask(1:nnelst))
         kglobvec1(1:nnelst1) = pack(kglobvec(1:nnelst),mask(1:nnelst))
         do iii = 1, 3
            posvec2(1:nnelst1,iii) = pack( posvec(1:nnelst,iii), 
     &                                     mask(1:nnelst)
     &                                   )
         enddo
         r2vec1(1:nnelst1) = posvec2(1:nnelst1,1)**2
     &                     + posvec2(1:nnelst1,2)**2
     &                     + posvec2(1:nnelst1,3)**2

cold           r = sqrt(r2)
cold           ck = rpole(1,kkpole)
cold           dkx = rpole(2,kkpole)
cold           dky = rpole(3,kkpole)
cold           dkz = rpole(4,kkpole)
cold           qkxx = rpole(5,kkpole)
cold           qkxy = rpole(6,kkpole)
cold           qkxz = rpole(7,kkpole)
cold           qkyy = rpole(9,kkpole)
cold           qkyz = rpole(10,kkpole)
cold           qkzz = rpole(13,kkpole)
cold           ukx = uind(1,kkpole)
cold           uky = uind(2,kkpole)
cold           ukz = uind(3,kkpole)
         rvec(1:nnelst1)      = sqrt(r2vec1(1:nnelst1))
         ckvec(1:nnelst1)     = rpole(1,kkpolevec1(1:nnelst1)) !
c
c        dkvec is  dkx, dky, dkz
c
         !DIR$ ASSUME_ALIGNED dkvec:64
         !DIR$ ASSUME_ALIGNED kkpolevec1:64
         dkvec(1:nnelst1,1:3) = reshape (
     &                rpole(2:4,kkpolevec1(1:nnelst1)) ,
     &                (/nnelst1,3/),ORDER=(/2,1/))

c
c            qkvec is  qkxx, qkxy, qkxz
c                      qkxy, qkyy, qkyz
c                      qkxz, qkyz, qkzz
cnew     !DIR$ ASSUME_ALIGNED qkvec:64
cnew     !DIR$ ASSUME_ALIGNED kkpolevec1:64
         qkvec(1:nnelst1,1:3,1:3)=
     &           reshape ( (/
     &                        rpole(5,kkpolevec1(1:nnelst1)),
     &                        rpole(6,kkpolevec1(1:nnelst1)),
     &                        rpole(7,kkpolevec1(1:nnelst1)),
     &                        rpole(6,kkpolevec1(1:nnelst1)),
     &                        rpole(9,kkpolevec1(1:nnelst1)),
     &                        rpole(10,kkpolevec1(1:nnelst1)),
     &                        rpole(7,kkpolevec1(1:nnelst1)),
     &                        rpole(10,kkpolevec1(1:nnelst1)),
     &                        rpole(13,kkpolevec1(1:nnelst1))
     &                     /),(/nnelst1,3,3/),ORDER=(/1,2,3/)
     &                   )
c

         ukvec(1:nnelst1,1:3) =
     &               transpose (uind(1:3,kkpolevec1(1:nnelst1))) 
c
c     get reciprocal distance terms for this interaction
c
cold           rr1 = f / r
cold           rr3 = rr1 / r2
cold           rr5 = 3.0_ti_p * rr3 / r2
cold           rr7 = 5.0_ti_p * rr5 / r2
         rr3vec(1:nnelst1) =          f / rvec(1:nnelst1)**3 !         1
         rr5vec(1:nnelst1) =  3.0_ti_p * f / rvec(1:nnelst1)**5 !     3 * 1
         rr7vec(1:nnelst1) = 15.0_ti_p * f / rvec(1:nnelst1)**7 ! 5 * 3 * 1

c
c     calculate the real space Ewald error function terms
c
cold           ralpha = aewald * r
cold           bn(0) = erfc(ralpha) / r
cold           alsq2 = 2.0_ti_p * aewald**2
cold           alsq2n = 0.0_ti_p
cold           if (aewald .gt. 0.0_ti_p)  alsq2n = 1.0_ti_p / (sqrtpi*aewald)
cold           exp2a = exp(-ralpha**2)
cold           do j = 1, 3
cold              bfac = real(j+j-1,t_p)
cold              alsq2n = alsq2 * alsq2n
cold              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
cold           end do
cold           do j = 0, 3
cold              bn(j) = f * bn(j)
cold           end do
         ralphavec(1:nnelst1)= aewald * rvec(1:nnelst1)
c
         call vderfc(nnelst1,ralphavec,bnvec(:,0))

         bnvec(1:nnelst1,0)  =  bnvec(1:nnelst1,0)  /  rvec(1:nnelst1)
         exp2avec(1:nnelst1) = exp(-ralphavec(1:nnelst1)**2)
         bnvec(1:nnelst1,1)  = (  1.0_ti_p * bnvec(1:nnelst1,0)
     &                         + alsq2 * alsq2n * exp2avec(1:nnelst1)
     &                         )    / r2vec1(1:nnelst1)
         bnvec(1:nnelst1,2)  = (  3.0_ti_p * bnvec(1:nnelst1,1)
     &                         + alsq2**2 * alsq2n * exp2avec(1:nnelst1)
     &                         )    / r2vec1(1:nnelst1)
         bnvec(1:nnelst1,3)  = (  5.0_ti_p * bnvec(1:nnelst1,2)
     &                         + alsq2**3 * alsq2n * exp2avec(1:nnelst1)
     &                         )    / r2vec1(1:nnelst1)

         bnvec(1:nnelst1,:)= f * bnvec(1:nnelst1,:)

c
c     apply Thole polarization damping to scale factors
c
         !DIR$ ASSUME_ALIGNED kkpolevec1:64
         !DIR$ ASSUME_ALIGNED dampvec:64
         !DIR$ ASSUME_ALIGNED dampvec1:64
         dampvec(1:nnelst1)  = pdi * pdamp(kkpolevec1(1:nnelst1))
         dampvec1(1:nnelst1) = dampvec(1:nnelst1)
         do iii = 1, nnelst1
             pgammavec(iii) = min(pti,thole(kkpolevec1(iii)))
         enddo
         where (dampvec(1:nnelst1) .ne. 0.0_ti_p)
                dampvec1 =   - pgammavec * (rvec/dampvec)**3
           where (dampvec1(:). gt.-50.0_ti_p)
             expdampvec1 = exp(dampvec1)
             sc3vec =   1.0_ti_p - expdampvec1
             sc5vec =   1.0_ti_p - (1.0_ti_p - dampvec1) * expdampvec1
             sc7vec =   1.0_ti_p - (1.0_ti_p - dampvec1 +
     &                                           0.6_ti_p *dampvec1**2)
     &                         * expdampvec1
            end where
         end where

c
c     intermediates involving Thole damping and scale factors
c
         psc3vec(1:nnelst1) =   1.0_ti_p 
     &                       -  sc3vec          (1:nnelst1)
     &                        * pscale(kglobvec1(1:nnelst1))
         psc5vec(1:nnelst1) =   1.0_ti_p
     &                       -  sc5vec          (1:nnelst1)
     &                        * pscale(kglobvec1(1:nnelst1))
         psc7vec(1:nnelst1) =   1.0_ti_p
     &                       -  sc7vec          (1:nnelst1)
     &                        * pscale(kglobvec1(1:nnelst1))

         psr3vec(1:nnelst1) =   bnvec           (1:nnelst1,1)
     &                       -  psc3vec         (1:nnelst1)
     &                        * rr3vec          (1:nnelst1)
         psr5vec(1:nnelst1) =   bnvec           (1:nnelst1,2)
     &                       -  psc5vec         (1:nnelst1)
     &                        * rr5vec          (1:nnelst1)
         psr7vec(1:nnelst1) =   bnvec           (1:nnelst1,3)
     &                       -  psc7vec         (1:nnelst1)
     &                        * rr7vec          (1:nnelst1)

c
c     intermediates involving moments and distance separation
c
         !DIR$ ASSUME_ALIGNED posvec2:64
         !DIR$ ASSUME_ALIGNED qivec:64
         !DIR$ ASSUME_ALIGNED qrivec:64
         !DIR$ ASSUME_ALIGNED qkvec:64
         !DIR$ ASSUME_ALIGNED qrkvec:64
         do iii = 1,nnelst1
            qritmp = matmul(qivec(:,:), posvec2(iii,:))
            qrivec(iii,:) = qritmp
            qrktmp = matmul(qkvec(iii,:,:), posvec2(iii,:))
            qrkvec(iii,:) = qrktmp
         enddo
         !DIR$ ASSUME_ALIGNED posvec2:64
         do iii = 1,nnelst1
            drivec(iii)  =  dot_product(divec(:),posvec2(iii,:))
            drkvec(iii)  =  dot_product(dkvec(iii,:),posvec2(iii,:)) ! 
            qrrivec(iii) =  dot_product(qrivec(iii,:),posvec2(iii,:)) ! 
            qrrkvec(iii) =  dot_product(qrkvec(iii,:),posvec2(iii,:))
            urivec(iii)  =  dot_product(uivec(:),posvec2(iii,:)) ! 
            urkvec(iii)  =  dot_product(ukvec(iii,:),posvec2(iii,:)) ! 
            duikvec(iii) =  dot_product(divec(:),ukvec(iii,:)) ! 
     &                    + dot_product(dkvec(iii,:),uivec(:)) ! 
            quikvec(iii) =  dot_product(qrivec(iii,:),ukvec(iii,:))
     &                    - dot_product(qrkvec(iii,:),uivec(:))
         enddo

c
c     calculate intermediate terms for polarization interaction
c
         term1vec(1:nnelst1) =  ckvec(1:nnelst1) * urivec(1:nnelst1)
     &                        - ci               * urkvec(1:nnelst1)
     &                        + duikvec(1:nnelst1)
         term2vec(1:nnelst1) =  2.0_ti_p          * quikvec(1:nnelst1)
     &                        - urivec(1:nnelst1) * drkvec(1:nnelst1)
     &                        - drivec(1:nnelst1) * urkvec(1:nnelst1)
         term3vec(1:nnelst1) =  urivec(1:nnelst1) * qrrkvec(1:nnelst1)
     &                        - urkvec(1:nnelst1) * qrrivec(1:nnelst1)
c
c     compute the energy contribution for this interaction
c
         e = sum (
     &              term1vec(1:nnelst1) * psr3vec(1:nnelst1)
     &            + term2vec(1:nnelst1) * psr5vec(1:nnelst1)
     &            + term3vec(1:nnelst1) * psr7vec(1:nnelst1)
     &           )


c
c     increment the overall polarization energy components
c
               ep = ep + e
c
c     reset exclusion coefficients for connected atoms
c
         pscale(pscalevec(1:ntot)) = 1.0_ti_p
      end do

c
c     perform deallocation of some local arrays
c
      deallocate (term3vec)
      deallocate (term2vec)
      deallocate (term1vec)
      deallocate (psr7vec)
      deallocate (psr5vec)
      deallocate (psr3vec)
      deallocate (psc7vec)
      deallocate (psc5vec)
      deallocate (psc3vec)
      deallocate (sc7vec)
      deallocate (sc5vec)
      deallocate (sc3vec)
      deallocate (pgammavec,expdampvec1)
      deallocate (dampvec,dampvec1)
      deallocate (exp2avec)
      deallocate (bnvec)
      deallocate (ralphavec)
      deallocate (rr7vec)
      deallocate (rr5vec)
      deallocate (rr3vec)
      deallocate (duikvec,quikvec)
      deallocate (urivec,urkvec)
      deallocate (qrrivec,qrrkvec)
      deallocate (qrkvec)
      deallocate (qrivec)
      deallocate (drkvec)
      deallocate (drivec)
      deallocate (mask1)
      deallocate (mask)
      deallocate (ukvec)
      deallocate (divec,uivec)
      deallocate (ckvec,dkvec,qkvec)
      deallocate (rvec)
      deallocate (r2vec,r2vec1)
      deallocate (posvec,posvec1,posvec2)
      deallocate (kkpolevec1)
      deallocate (kglobvec1,kbisvec1)
      deallocate (kglobvec,kbisvec,kkpolevec)
      deallocate (itmp15)
      deallocate (itmp14)
      deallocate (itmp13)
      deallocate (itmp12)
      deallocate (qivec,qritmp,qrktmp)
      deallocate (pscalevec)
      deallocate (pscale)
      return
      end
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine eprecipvec  --  PME recip space polarization energy  ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "eprecipvec" evaluates the reciprocal space portion of particle
c     mesh Ewald summation energy due to dipole polarization
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine eprecipvec
      use atmlst
      use atoms
      use bound
      use boxes
      use chgpot
      use domdec
      use energi
      use ewald
      use fft
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use tinheader ,only:ti_p,re_p
      use mpi

      implicit none
      integer ierr,iipole,proc
      integer status(MPI_STATUS_SIZE),tag,commloc
      integer nprocloc,rankloc
      integer i,j,k,iglob
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      real(t_p) e,r1,r2,r3
      real(t_p) f,h1,h2,h3
      real(t_p) volterm,denom
      real(t_p) hsq,expterm
      real(t_p) term,pterm
      real(t_p) struc2
      real(t_p) ftc(10,10)
      real(t_p) fuind(3)
      !DIR$ ATTRIBUTES ALIGN:64::fuindvec,iipolevec,a,nfft
      real(t_p),allocatable :: fuindvec(:,:),iipolevec(:),
     &                         a(:,:),nfft(:,:)
      !DIR$ ATTRIBUTES ALIGN:64::evec,fphirecvec,fuindvec2
      real(t_p),allocatable :: evec(:,:),fphirecvec(:,:),
     &                                    fuindvec2(:,:)

c
      if (use_pmecore) then
        nprocloc = nrec
        rankloc = rank_bis
        commloc =  comm_rec
      else
        nprocloc = nproc
        rankloc = rank
        commloc = MPI_COMM_WORLD
      end if
      allocate(iipolevec(npolerecloc))
      allocate(fuindvec(3,npolerecloc))
      allocate(fuindvec2(npolerecloc,3))
      allocate(evec(npolerecloc,3))
      allocate(fphirecvec(npolerecloc,3))
      allocate(a(3,3),nfft(3,3))
      write(*,*) 'eprecipvec'
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
cc
cc     get the fractional to Cartesian transformation matrix
cc
      call frac_to_cart (ftc)
cc
cc     initialize variables required for the scalar summation
cc
c
c     convert Cartesian induced dipoles to fractional coordinates
c
       do i = 1, 3
         a(1,i) = real(nfft1,t_p) * recip(i,1)
         a(2,i) = real(nfft2,t_p) * recip(i,2)
         a(3,i) = real(nfft3,t_p) * recip(i,3)
      end do
      e = 0.0_ti_p
      iipolevec(1:npolerecloc)=polerecglob(1:npolerecloc)
      fphirecvec(1:npolerecloc,:) = transpose(
     &                   EOSHIFT(fphirec(:,1:npolerecloc),1))
       !DIR ASSUME_ALIGNED a:64
       !DIR ASSUME_ALIGNED fuindvec:64
       do i = 1, npolerecloc
          fuindvec(:,i) = matmul(a,uind(:,polerecglob(i)))
       end do
       !DIR ASSUME_ALIGNED evec:64
       !DIR ASSUME_ALIGNED fuindvec2:64
       !DIR ASSUME_ALIGNED fuindvec:64
       fuindvec2(1:npolerecloc,:) = transpose(fuindvec(:,1:npolerecloc))

       evec(1:npolerecloc,:) =  fuindvec2(1:npolerecloc,:)
     &                        * fphirecvec(1:npolerecloc,:)
       e = sum(evec(1:npolerecloc,:))

      e  = 0.5_ti_p * electric * e
      ep = ep + e
c
c     account for zeroth grid point for nonperiodic system
c
      if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1)
     $   .and.(kstart2(rankloc+1).eq.1)) then
        if (.not. use_bounds) then
           expterm = 0.5_ti_p * pi / xbox
           struc2 = qgrid2in_2d(1,1,1,1,1)**2 +
     $       qgrid2in_2d(2,1,1,1,1)**2
           e = f * expterm * struc2
           ep = ep + e
        end if
      end if
c
      return
      end
