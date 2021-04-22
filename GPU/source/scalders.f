
      subroutine scalderfieldzmat6( a1, a2, d1, e1, d2, e2,d3,
     $                              e3, ade, adme, adte, adtb)

      !TODO 1.2 merge this file with the code precision and mpi comm

      use atmlst
      use atoms
      use bound
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
      use potent
      use shunt
      implicit none


      real*8, dimension(3,npolebloc), intent(in) :: a1, a2
      real*8, dimension(3,npolebloc), intent(in) :: e1, e2, e3, d1, d2,
     $                                               d3
      real*8, dimension(3,npolebloc), intent(out) :: adtb
      real*8, dimension(3,2,npolebloc), intent(out) :: adme, ade
      real*8, dimension(3,3,2,npolebloc), intent(out) :: adte
      
      real*8, allocatable, dimension(:,:,:) ::adE_cor,admE_cor,admE_self
      real*8, allocatable, dimension(:,:,:,:) ::adtE_cor




      integer :: ii, iipole, iglob, i, jjpole, jglob,
     $           jjj, jbis, k, j

      real*8 :: pti, pdi, rrij3, rrij5, rrij7, rrij9,
     $         urr5, urr7, 
     $scale3, scale5, 
     $          scale7, scale9, dsc3, dsc5, dsc7, dsc9, psc3, psc5,
     $          psc7, psc9, usc5, usc7, scalmuir, scalmujr, 
     $quadterm, quadterm_bis, 
     $          ci, cj, bfac, invrij2, alsq2n, cutoff2, exp2a, ralpha,
     $          d, rij2, damp, expdamp, pgamma, alsq2, term

      real*8                :: bn(0:4) 
      real*8, dimension(2)  :: srr3, srr5, srr7, srr9
      real*8, dimension(3)  :: rij, thetajr, thetair, di, dj
      real*8, dimension(3,3):: qi, qj

      real*8, allocatable, dimension(:) :: dscale, pscale, uscale
      real*8, allocatable, dimension(:,:) :: adtbcor
      

      real*8:: a1xi,a1yi,a1zi,a1xj,a1yj,a1zj
      real*8:: a2xi,a2yi,a2zi,a2xj,a2yj,a2zj
      real*8,dimension(2):: term1,term2,term3,term4,term5,term6,term7
      real*8,dimension(2):: dexdxi,dexdxj,deydyi,deydyj,dezdzi,dezdzj
      real*8,dimension(2):: dexdyi,dexdyj,dexdzi,dexdzj,deydzi,deydzj
      real*8,dimension(2):: demxdx,demxdy,demxdz,demydy,demydz,demzdz
      real*8,dimension(2):: detxdx,detxdy,detxdz,detydy,detydz,detzdz
      real*8 :: depx,depy,depz,dempx,dempy,dempz,detx,dety,detz
      real*8 :: scale1i,scale1j,scale2i,scale2j,scale3i,scale3j
      real*8::temp1jxx5,temp1jxy5,temp1jxz5,temp1jyy5,temp1jyz5,
     $     temp1jzz5
      real*8::temp1jxx7,temp1jxy7,temp1jxz7,temp1jyy7,temp1jyz7,
     $     temp1jzz7
      real*8::temp1jxx,temp1jxy,temp1jxz,temp1jyy,temp1jyz,temp1jzz
      real*8::temp1ixx5,temp1ixy5,temp1ixz5,temp1iyy5,temp1iyz5,
     $     temp1izz5
      real*8::temp1ixx7,temp1ixy7,temp1ixz7,temp1iyy7,temp1iyz7,
     $     temp1izz7
      real*8::temp1ixx,temp1ixy,temp1ixz,temp1iyy,temp1iyz,temp1izz

      real*8::temp2jxx5,temp2jxy5,temp2jxz5,temp2jyy5,temp2jyz5,
     $     temp2jzz5
      real*8::temp2jxx7,temp2jxy7,temp2jxz7,temp2jyy7,temp2jyz7,
     $     temp2jzz7
      real*8::temp2jxx,temp2jxy,temp2jxz,temp2jyy,temp2jyz,temp2jzz
      real*8::temp2ixx5,temp2ixy5,temp2ixz5,temp2iyy5,temp2iyz5,
     $     temp2izz5
      real*8::temp2ixx7,temp2ixy7,temp2ixz7,temp2iyy7,temp2iyz7,
     $     temp2izz7
      real*8::temp2ixx,temp2ixy,temp2ixz,temp2iyy,temp2iyz,temp2izz

      real*8::temp3jxx5,temp3jxy5,temp3jxz5,temp3jyy5,temp3jyz5,
     $     temp3jzz5
      real*8::temp3jxx7,temp3jxy7,temp3jxz7,temp3jyy7,temp3jyz7,
     $     temp3jzz7
      real*8::temp3jxx,temp3jxy,temp3jxz,temp3jyy,temp3jyz,temp3jzz
      real*8::temp3ixx5,temp3ixy5,temp3ixz5,temp3iyy5,temp3iyz5,
     $     temp3izz5
      real*8::temp3ixx7,temp3ixy7,temp3ixz7,temp3iyy7,temp3iyz7,
     $     temp3izz7
      real*8::temp3ixx,temp3ixy,temp3ixz,temp3iyy,temp3iyz,temp3izz

      character(10) :: mode

 1000 format(' Warning, system moved too much since last neighbor list
     $   update, try lowering nlupdate')

      

      ade = 0d0
      adme = 0d0
      adte = 0d0
      adtb = 0d0

      mode = 'EWALD'
      if (use_polarshortreal) mode = 'SHORTEWALD'
      call switch (mode)
      cutoff2 = cut2

      allocate(dscale(n), pscale(n), uscale(n))
      allocate(adtbcor(3,npolebloc))
      allocate(adE_cor(3,2,npolebloc),admE_cor(3,2,npolebloc))
      allocate(admE_self(3,2,npolebloc),adtE_cor(3,3,2,npolebloc))
      adme_cor = 0d0
      adme_self = 0d0
      adte_cor = 0d0
      ade_cor = 0d0
      dscale = 1d0
      pscale = 1d0
      uscale = 1.0d0
      adtbcor = 0d0

      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = poleloc(iglob)
         if ((i.eq.0).or.(i.gt.nbloc)) then
            write(iout,1000)
            cycle
         end if
         pdi = pdamp(iipole)
         pti = thole(iipole)

         ci      = rpole(1,iipole)
         di(1)   = rpole(2,iipole)
         di(2)   = rpole(3,iipole)
         di(3)   = rpole(4,iipole)
         qi(1,1) = rpole(5, iipole)
         qi(2,1) = rpole(6, iipole)
         qi(3,1) = rpole(7, iipole)
         qi(1,2) = rpole(6, iipole)
         qi(2,2) = rpole(9, iipole)
         qi(3,2) = rpole(10,iipole)
         qi(1,3) = rpole(7, iipole)
         qi(2,3) = rpole(10,iipole)
         qi(3,3) = rpole(13,iipole)

         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = p2scale
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = p3scale
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = p4scale
            do k = 1, np11(iglob)
               if (i14(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i14(j,iglob)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = p5scale
         end do
         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = d1scale
            uscale(ip11(j,iglob)) = u1scale
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = d2scale
            uscale(ip12(j,iglob)) = u2scale
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = d3scale
            uscale(ip13(j,iglob)) = u3scale
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = d4scale
            uscale(ip14(j,iglob)) = u4scale
         end do

         do jjj = 1, nelst(ii)
            jjpole = elst(jjj,ii)
            jglob = ipole(jjpole)
            jbis = poleloc(jjpole)
            if (jbis.eq.0) then
               write(iout,1000)
               cycle
            end if

            !Distances
            rij(1) = x(jglob) - x(iglob)
            rij(2) = y(jglob) - y(iglob)
            rij(3) = z(jglob) - z(iglob)
            call image(rij(1), rij(2), rij(3))
            rij2 = dot_product(rij,rij)
            if (rij2 .gt. cutoff2) cycle
            d  = sqrt(rij2)
            invrij2 = 1d0/rij2
            rrij3 = invrij2/d
            rrij5 = 3d0*rrij3*invrij2
            rrij7 = 5d0*rrij5*invrij2
            rrij9 = 7d0*rrij7*invrij2

            !Multipoles
            cj      = rpole(1,jjpole)
            dj(1)   = rpole(2,jjpole)
            dj(2)   = rpole(3,jjpole)
            dj(3)   = rpole(4,jjpole)
            qj(1,1) = rpole(5, jjpole)
            qj(2,1) = rpole(6, jjpole)
            qj(3,1) = rpole(7, jjpole)
            qj(1,2) = rpole(6, jjpole)
            qj(2,2) = rpole(9, jjpole)
            qj(3,2) = rpole(10,jjpole)
            qj(1,3) = rpole(7, jjpole)
            qj(2,3) = rpole(10,jjpole)
            qj(3,3) = rpole(13,jjpole)

            scalmujr = dot_product(dj,rij)!sprod(3,dj,rij)
            scalmuir = dot_product(di,rij)!sprod(3,di,rij)
            thetajr(1) = dot_product(qj(1,:),rij)
            thetajr(2) = dot_product(qj(2,:),rij)
            thetajr(3) = dot_product(qj(3,:),rij)
            thetair(1) = dot_product(qi(1,:),rij)
            thetair(2) = dot_product(qi(2,:),rij)
            thetair(3) = dot_product(qi(3,:),rij)
            quadterm     = dot_product(thetajr,rij)!sprod(3,thetajr,rij)
            quadterm_bis = dot_product(thetair,rij)!sprod(3,thetair,rij)

            scale1j = dot_product(e1(:,jbis),rij)
            scale1i = dot_product(e1(:,i),rij)
            scale2j = dot_product(e2(:,jbis),rij)
            scale2i = dot_product(e2(:,i),rij)
            scale3j = dot_product(e3(:,jbis),rij)
            scale3i = dot_product(e3(:,i),rij)

            a1xi = a1(1,i)
            a1yi = a1(2,i)
            a1zi = a1(3,i)
            a1xj = a1(1,jbis)
            a1yj = a1(2,jbis)
            a1zj = a1(3,jbis)
            a2xi = a2(1,i)
            a2yi = a2(2,i)
            a2zi = a2(3,i)
            a2xj = a2(1,jbis)
            a2yj = a2(2,jbis)
            a2zj = a2(3,jbis)

            !Errfunc damping
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0d0 * aewald**2
            alsq2n = 0.0d0
            if (aewald .gt. 0.0d0)
     &        alsq2n = 1.0d0 / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 4 
              bfac = dble(j+j-1)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a)*invrij2
            end do

            !Scalings
            damp = pdi*pdamp(jjpole)
            scale3 = 1d0
            scale5 = 1d0
            scale7 = 1d0
            scale9 = 1d0
            if (damp.ne.0d0) then
               pgamma = min(pti,thole(jjpole))
               damp = -pgamma*(d/damp)**3
               if (damp .gt. -50d0) then
                  expdamp = exp(damp)
                  scale3 = scale3 * (1.0d0 - expdamp)
                  scale5 = scale5 * (1.0d0 - expdamp*(1.0d0 - damp))
                  scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                  scale9 = 1d0 - expdamp*(1d0 - damp + (18d0*damp**2 -
     $                     9d0*damp**3)/35d0)
               end if
            end if
            dsc3 = scale3*dscale(jglob)
            dsc5 = scale5*dscale(jglob)
            dsc7 = scale7*dscale(jglob)
            dsc9 = scale9*dscale(jglob)
            psc3 = scale3*pscale(jglob)
            psc5 = scale5*pscale(jglob)
            psc7 = scale7*pscale(jglob)
            psc9 = scale9*pscale(jglob)
            !mat
            usc5 = scale5*uscale(jglob)
            usc7 = scale7*uscale(jglob)
            srr3(1) = (1d0 - dsc3)*rrij3 - bn(1)
            srr3(2) = (1d0 - psc3)*rrij3 - bn(1)
            srr5(1) = (1d0 - dsc5)*rrij5 - bn(2)
            srr5(2) = (1d0 - psc5)*rrij5 - bn(2)
            srr7(1) = (1d0 - dsc7)*rrij7 - bn(3)
            srr7(2) = (1d0 - psc7)*rrij7 - bn(3)
            srr9(1) = (1d0 - dsc9)*rrij9 - bn(4)
            srr9(2) = (1d0 - psc9)*rrij9 - bn(4)
            !mat
            urr5 = (1d0 - usc5)*rrij5 - bn(2)
            urr7 = (1d0 - usc7)*rrij7 - bn(3)

            !Actual field calculations ! 

c dexdxj (alpha=beta=1)
            term1 = rij(1)*rij(1)*srr5 - srr3 
            term2 = 2d0*srr5*rij(1) 
            term3 = -srr7*rij(1)*rij(1) + srr5
            term4 = -4d0*rij(1)*srr7
            term5 = srr9*rij(1)*rij(1) - srr7

            dexdxj = cj*term1 + dj(1)*term2
     $      + scalmujr*term3 + 2d0*qj(1,1)*srr5
     $      + thetajr(1)*term4 + quadterm*term5

            dexdxi = ci*term1 - di(1)*term2
     $      - scalmuir*term3 + 2d0*qi(1,1)*srr5
     $      + thetair(1)*term4 + quadterm_bis*term5
            
            demxdx = term1
            detxdx = 2d0*srr5-srr7*rij(1)*rij(1)

c deydyj (alpha=beta=2)
            term1 = rij(2)*rij(2)*srr5 - srr3 
            term2 = 2d0*srr5*rij(2) 
            term3 = -srr7*rij(2)*rij(2) + srr5
            term4 = -4d0*rij(2)*srr7
            term5 = srr9*rij(2)*rij(2) - srr7

            deydyj = cj*term1 + dj(2)*term2
     $      + scalmujr*term3 + 2d0*qj(2,2)*srr5
     $      + thetajr(2)*term4 + quadterm*term5

            deydyi = ci*term1 - di(2)*term2
     $      - scalmuir*term3 + 2d0*qi(2,2)*srr5
     $      + thetair(2)*term4 + quadterm_bis*term5

            demydy = term1
            detydy = 2d0*srr5-srr7*rij(2)*rij(2)

c dexdyj (alpha=1 beta=2)
            term1 = rij(1)*rij(2)*srr5 
            term2 = srr5*rij(2) 
            term3 = srr5*rij(1) 
            term4 = -srr7*rij(1)*rij(2)
            term5 = -2d0*rij(2)*srr7
            term6 = -2d0*rij(1)*srr7
            term7 = srr9*rij(1)*rij(2)

            dexdyj = cj*term1 + dj(1)*term2 + dj(2)*term3
     $      + scalmujr*term4 + 2d0*qj(1,2)*srr5
     $      + thetajr(1)*term5 +  thetajr(2)*term6
     $      + quadterm*term7

            dexdyi = ci*term1 - di(1)*term2 - di(2)*term3
     $      - scalmuir*term4 + 2d0*qi(1,2)*srr5
     $      + thetair(1)*term5 +  thetair(2)*term6
     $      + quadterm_bis*term7

            demxdy = term1
            detxdy = -srr7*rij(1)*rij(2)


c dezdzj (alpha3=beta=3)
            term1 = rij(3)*rij(3)*srr5 - srr3 
            term2 = 2d0*srr5*rij(3) 
            term3 = -srr7*rij(3)*rij(3) + srr5
            term4 = -4d0*rij(3)*srr7
            term5 = srr9*rij(3)*rij(3) - srr7

            dezdzj = cj*term1 + dj(3)*term2
     $      + scalmujr*term3 + 2d0*qj(3,3)*srr5
     $      + thetajr(3)*term4 + quadterm*term5

            dezdzi = ci*term1 - di(3)*term2
     $      - scalmuir*term3 + 2d0*qi(3,3)*srr5
     $      + thetair(3)*term4 + quadterm_bis*term5

            demzdz = term1
            detzdz = 2d0*srr5-srr7*rij(3)*rij(3)

c dexdzj (alpha=1 beta=3)
            term1 = rij(1)*rij(3)*srr5 
            term2 = srr5*rij(3) 
            term3 = srr5*rij(1) 
            term4 = -srr7*rij(1)*rij(3)
            term5 = -2d0*rij(3)*srr7
            term6 = -2d0*rij(1)*srr7
            term7 = srr9*rij(1)*rij(3)

            dexdzj = cj*term1 + dj(1)*term2 + dj(3)*term3
     $      + scalmujr*term4 + 2d0*qj(1,3)*srr5
     $      + thetajr(1)*term5 +  thetajr(3)*term6
     $      + quadterm*term7

            dexdzi = ci*term1 - di(1)*term2 - di(3)*term3
     $      - scalmuir*term4 + 2d0*qi(1,3)*srr5
     $      + thetair(1)*term5 +  thetair(3)*term6
     $      + quadterm_bis*term7

            demxdz = term1
            detxdz = -srr7*rij(1)*rij(3)


c deydzj (alpha=2 beta=3)
            term1 = rij(2)*rij(3)*srr5 
            term2 = srr5*rij(3) 
            term3 = srr5*rij(2) 
            term4 = -srr7*rij(2)*rij(3)
            term5 = -2d0*rij(3)*srr7
            term6 = -2d0*rij(2)*srr7
            term7 = srr9*rij(2)*rij(3)

            deydzj = cj*term1 + dj(2)*term2 + dj(3)*term3
     $      + scalmujr*term4 + 2d0*qj(2,3)*srr5
     $      + thetajr(2)*term5 +  thetajr(3)*term6
     $      + quadterm*term7

            deydzi = ci*term1 - di(2)*term2 - di(3)*term3
     $      - scalmuir*term4 + 2d0*qi(2,3)*srr5
     $      + thetair(2)*term5 +  thetair(3)*term6
     $      + quadterm_bis*term7

            demydz = term1
            detydz = -srr7*rij(2)*rij(3)

            depx = -a1xi*dexdxj(1) + a1xj*dexdxi(1)
     $     -a1yi*dexdyj(1) + a1yj*dexdyi(1)
     $     -a1zi*dexdzj(1) + a1zj*dexdzi(1)
            depy = -a1xi*dexdyj(1) + a1xj*dexdyi(1)
     $     -a1yi*deydyj(1) + a1yj*deydyi(1)
     $     -a1zi*deydzj(1) + a1zj*deydzi(1)
            depz = -a1xi*dexdzj(1) + a1xj*dexdzi(1)
     $     -a1yi*deydzj(1) + a1yj*deydzi(1)
     $     -a1zi*dezdzj(1) + a1zj*dezdzi(1)

            ade_cor(1,1,i) = ade_cor(1,1,i) + depx
            ade_cor(2,1,i) = ade_cor(2,1,i) + depy
            ade_cor(3,1,i) = ade_cor(3,1,i) + depz

            ade_cor(1,1,jbis) = ade_cor(1,1,jbis) - depx
            ade_cor(2,1,jbis) = ade_cor(2,1,jbis) - depy
            ade_cor(3,1,jbis) = ade_cor(3,1,jbis) - depz

            depx = -a2xi*dexdxj(2) + a2xj*dexdxi(2)
     $     -a2yi*dexdyj(2) + a2yj*dexdyi(2)
     $     -a2zi*dexdzj(2) + a2zj*dexdzi(2)
            depy = -a2xi*dexdyj(2) + a2xj*dexdyi(2)
     $     -a2yi*deydyj(2) + a2yj*deydyi(2)
     $     -a2zi*deydzj(2) + a2zj*deydzi(2)
            depz = -a2xi*dexdzj(2) + a2xj*dexdzi(2)
     $     -a2yi*deydzj(2) + a2yj*deydzi(2)
     $     -a2zi*dezdzj(2) + a2zj*dezdzi(2)

            ade_cor(1,2,i) = ade_cor(1,2,i) + depx
            ade_cor(2,2,i) = ade_cor(2,2,i) + depy
            ade_cor(3,2,i) = ade_cor(3,2,i) + depz

            ade_cor(1,2,jbis) = ade_cor(1,2,jbis) - depx
            ade_cor(2,2,jbis) = ade_cor(2,2,jbis) - depy
            ade_cor(3,2,jbis) = ade_cor(3,2,jbis) - depz

            dempx = a1xj*demxdx(1) + a1yj*demxdy(1)
     $     + a1zj*demxdz(1)
            dempy = a1xj*demxdy(1) + a1yj*demydy(1)
     $     + a1zj*demydz(1)
            dempz = a1xj*demxdz(1) + a1yj*demydz(1)
     $     + a1zj*demzdz(1)

            adme_cor(1,1,i) = adme_cor(1,1,i) + dempx
            adme_cor(2,1,i) = adme_cor(2,1,i) + dempy
            adme_cor(3,1,i) = adme_cor(3,1,i) + dempz

            dempx = a1xi*demxdx(1) + a1yi*demxdy(1)
     $     + a1zi*demxdz(1)
            dempy = a1xi*demxdy(1) + a1yi*demydy(1)
     $     + a1zi*demydz(1)
            dempz = a1xi*demxdz(1) + a1yi*demydz(1)
     $     + a1zi*demzdz(1)

            adme_cor(1,1,jbis) = adme_cor(1,1,jbis) + dempx
            adme_cor(2,1,jbis) = adme_cor(2,1,jbis) + dempy
            adme_cor(3,1,jbis) = adme_cor(3,1,jbis) + dempz

            dempx = a2xj*demxdx(2) + a2yj*demxdy(2)
     $     + a2zj*demxdz(2)
            dempy = a2xj*demxdy(2) + a2yj*demydy(2)
     $     + a2zj*demydz(2)
            dempz = a2xj*demxdz(2) + a2yj*demydz(2)
     $     + a2zj*demzdz(2)

            adme_cor(1,2,i) = adme_cor(1,2,i) + dempx
            adme_cor(2,2,i) = adme_cor(2,2,i) + dempy
            adme_cor(3,2,i) = adme_cor(3,2,i) + dempz

            dempx = a2xi*demxdx(2) + a2yi*demxdy(2)
     $     + a2zi*demxdz(2)
            dempy = a2xi*demxdy(2) + a2yi*demydy(2)
     $     + a2zi*demydz(2)
            dempz = a2xi*demxdz(2) + a2yi*demydz(2)
     $     + a2zi*demzdz(2)

            adme_cor(1,2,jbis) = adme_cor(1,2,jbis) + dempx
            adme_cor(2,2,jbis) = adme_cor(2,2,jbis) + dempy
            adme_cor(3,2,jbis) = adme_cor(3,2,jbis) + dempz

            detx = -a1xj*detxdx(1)-a1yj*detxdy(1)
     $     -a1zj*detxdz(1) 
            dety = -a1xj*detxdy(1)-a1yj*detydy(1)
     $     -a1zj*detydz(1) 
            detz = -a1xj*detxdz(1)-a1yj*detydz(1)
     $     -a1zj*detzdz(1) 

            adtE_cor(1,:,1,i) = adtE_cor(1,:,1,i) + detx*rij(:)
            adtE_cor(2,:,1,i) = adtE_cor(2,:,1,i) + dety*rij(:)
            adtE_cor(3,:,1,i) = adtE_cor(3,:,1,i) + detz*rij(:)


            detx = a1xi*detxdx(1)+a1yi*detxdy(1)
     $      + a1zi*detxdz(1) 
            dety = a1xi*detxdy(1)+a1yi*detydy(1)
     $      + a1zi*detydz(1) 
            detz = a1xi*detxdz(1)+a1yi*detydz(1)
     $      + a1zi*detzdz(1) 

            adtE_cor(1,:,1,jbis) = adtE_cor(1,:,1,jbis) + detx*rij(:)
            adtE_cor(2,:,1,jbis) = adtE_cor(2,:,1,jbis) + dety*rij(:)
            adtE_cor(3,:,1,jbis) = adtE_cor(3,:,1,jbis) + detz*rij(:)

            detx = -a2xj*detxdx(2)-a2yj*detxdy(2)
     $     -a2zj*detxdz(2) 
            dety = -a2xj*detxdy(2)-a2yj*detydy(2)
     $     -a2zj*detydz(2) 
            detz = -a2xj*detxdz(2)-a2yj*detydz(2)
     $     -a2zj*detzdz(2) 

            adtE_cor(1,:,2,i) = adtE_cor(1,:,2,i) + detx*rij(:)
            adtE_cor(2,:,2,i) = adtE_cor(2,:,2,i) + dety*rij(:)
            adtE_cor(3,:,2,i) = adtE_cor(3,:,2,i) + detz*rij(:)

            detx = a2xi*detxdx(2)+a2yi*detxdy(2)
     $      + a2zi*detxdz(2) 
            dety = a2xi*detxdy(2)+a2yi*detydy(2)
     $      + a2zi*detydz(2) 
            detz = a2xi*detxdz(2)+a2yi*detydz(2)
     $      + a2zi*detzdz(2) 

            adtE_cor(1,:,2,jbis) = adtE_cor(1,:,2,jbis) + detx*rij(:)
            adtE_cor(2,:,2,jbis) = adtE_cor(2,:,2,jbis) + dety*rij(:)
            adtE_cor(3,:,2,jbis) = adtE_cor(3,:,2,jbis) + detz*rij(:)
c
c      also compute gradient of field of dipoles "e"
c
c      alpha = beta = 1
c
            temp1jxx5 =  urr5*(scale1i + 2*rij(1)*e1(1,i))
            temp1jxx7 =  -urr7*scale1i*rij(1)*rij(1)
            temp1jxx = temp1jxx5 + temp1jxx7
            temp1ixx5 =  urr5*(scale1j + 2*rij(1)*e1(1,jbis))
            temp1ixx7 =  -urr7*scale1j*rij(1)*rij(1)
            temp1ixx = temp1ixx5 + temp1ixx7

            temp2jxx5 =  urr5*(scale2i + 2*rij(1)*e2(1,i))
            temp2jxx7 =  -urr7*scale2i*rij(1)*rij(1)
            temp2jxx = temp2jxx5 + temp2jxx7
            temp2ixx5 =  urr5*(scale2j + 2*rij(1)*e2(1,jbis))
            temp2ixx7 =  -urr7*scale2j*rij(1)*rij(1)
            temp2ixx = temp2ixx5 + temp2ixx7

            temp3jxx5 =  urr5*(scale3i + 2*rij(1)*e3(1,i))
            temp3jxx7 =  -urr7*scale3i*rij(1)*rij(1)
            temp3jxx = temp3jxx5 + temp3jxx7
            temp3ixx5 =  urr5*(scale3j + 2*rij(1)*e3(1,jbis))
            temp3ixx7 =  -urr7*scale3j*rij(1)*rij(1)
            temp3ixx = temp3ixx5 + temp3ixx7

c
c      alpha = 2, beta = 1
c
            temp1jxy5 =  urr5*(rij(2)*e1(1,i)+rij(1)*e1(2,i))
            temp1jxy7 =  -urr7*scale1i*rij(2)*rij(1)
            temp1jxy = temp1jxy5 + temp1jxy7
            temp1ixy5 =  urr5*(rij(2)*e1(1,jbis)+rij(1)*e1(2,jbis))
            temp1ixy7 =  -urr7*scale1j*rij(2)*rij(1)
            temp1ixy = temp1ixy5 + temp1ixy7

            temp2jxy5 =  urr5*(rij(2)*e2(1,i)+rij(1)*e2(2,i))
            temp2jxy7 =  -urr7*scale2i*rij(2)*rij(1)
            temp2jxy = temp2jxy5 + temp2jxy7
            temp2ixy5 =  urr5*(rij(2)*e2(1,jbis)+rij(1)*e2(2,jbis))
            temp2ixy7 =  -urr7*scale2j*rij(2)*rij(1)
            temp2ixy = temp2ixy5 + temp2ixy7

            temp3jxy5 =  urr5*(rij(2)*e3(1,i)+rij(1)*e3(2,i))
            temp3jxy7 =  -urr7*scale3i*rij(2)*rij(1)
            temp3jxy = temp3jxy5 + temp3jxy7
            temp3ixy5 =  urr5*(rij(2)*e3(1,jbis)+rij(1)*e3(2,jbis))
            temp3ixy7 =  -urr7*scale3j*rij(2)*rij(1)
            temp3ixy = temp3ixy5 + temp3ixy7
c
c      alpha = 3, beta = 1
c
            temp1jxz5 =  urr5*(rij(3)*e1(1,i)+rij(1)*e1(3,i))
            temp1jxz7 =  -urr7*scale1i*rij(3)*rij(1)
            temp1jxz = temp1jxz5 + temp1jxz7
            temp1ixz5 =  urr5*(rij(3)*e1(1,jbis)+rij(1)*e1(3,jbis))
            temp1ixz7 =  -urr7*scale1j*rij(3)*rij(1)
            temp1ixz = temp1ixz5 + temp1ixz7

            temp2jxz5 =  urr5*(rij(3)*e2(1,i)+rij(1)*e2(3,i))
            temp2jxz7 =  -urr7*scale2i*rij(3)*rij(1)
            temp2jxz = temp2jxz5 + temp2jxz7
            temp2ixz5 =  urr5*(rij(3)*e2(1,jbis)+rij(1)*e2(3,jbis))
            temp2ixz7 =  -urr7*scale2j*rij(3)*rij(1)
            temp2ixz = temp2ixz5 + temp2ixz7

            temp3jxz5 =  urr5*(rij(3)*e3(1,i)+rij(1)*e3(3,i))
            temp3jxz7 =  -urr7*scale3i*rij(3)*rij(1)
            temp3jxz = temp3jxz5 + temp3jxz7
            temp3ixz5 =  urr5*(rij(3)*e3(1,jbis)+rij(1)*e3(3,jbis))
            temp3ixz7 =  -urr7*scale3j*rij(3)*rij(1)
            temp3ixz = temp3ixz5 + temp3ixz7

c
c      alpha = beta = 2
c
            temp1jyy5 =  urr5*(scale1i + 2*rij(2)*e1(2,i))
            temp1jyy7 =  -urr7*scale1i*rij(2)*rij(2)
            temp1jyy = temp1jyy5 + temp1jyy7
            temp1iyy5 =  urr5*(scale1j + 2*rij(2)*e1(2,jbis))
            temp1iyy7 =  -urr7*scale1j*rij(2)*rij(2)
            temp1iyy = temp1iyy5 + temp1iyy7

            temp2jyy5 =  urr5*(scale2i + 2*rij(2)*e2(2,i))
            temp2jyy7 =  -urr7*scale2i*rij(2)*rij(2)
            temp2jyy = temp2jyy5 + temp2jyy7
            temp2iyy5 =  urr5*(scale2j + 2*rij(2)*e2(2,jbis))
            temp2iyy7 =  -urr7*scale2j*rij(2)*rij(2)
            temp2iyy = temp2iyy5 + temp2iyy7

            temp3jyy5 =  urr5*(scale3i + 2*rij(2)*e3(2,i))
            temp3jyy7 =  -urr7*scale3i*rij(2)*rij(2)
            temp3jyy = temp3jyy5 + temp3jyy7
            temp3iyy5 =  urr5*(scale3j + 2*rij(2)*e3(2,jbis))
            temp3iyy7 =  -urr7*scale3j*rij(2)*rij(2)
            temp3iyy = temp3iyy5 + temp3iyy7
c
c      alpha = beta = 3
c
            temp1jzz5 =  urr5*(scale1i + 2*rij(3)*e1(3,i))
            temp1jzz7 =  -urr7*scale1i*rij(3)*rij(3)
            temp1jzz = temp1jzz5 + temp1jzz7
            temp1izz5 =  urr5*(scale1j + 2*rij(3)*e1(3,jbis))
            temp1izz7 =  -urr7*scale1j*rij(3)*rij(3)
            temp1izz = temp1izz5 + temp1izz7

            temp2jzz5 =  urr5*(scale2i + 2*rij(3)*e2(3,i))
            temp2jzz7 =  -urr7*scale2i*rij(3)*rij(3)
            temp2jzz = temp2jzz5 + temp2jzz7
            temp2izz5 =  urr5*(scale2j + 2*rij(3)*e2(3,jbis))
            temp2izz7 =  -urr7*scale2j*rij(3)*rij(3)
            temp2izz = temp2izz5 + temp2izz7

            temp3jzz5 =  urr5*(scale3i + 2*rij(3)*e3(3,i))
            temp3jzz7 =  -urr7*scale3i*rij(3)*rij(3)
            temp3jzz = temp3jzz5 + temp3jzz7
            temp3izz5 =  urr5*(scale3j + 2*rij(3)*e3(3,jbis))
            temp3izz7 =  -urr7*scale3j*rij(3)*rij(3)
            temp3izz = temp3izz5 + temp3izz7
c
c      alpha = 2, beta = 3
c
            temp1jyz5 =  urr5*(rij(3)*e1(2,i)+rij(2)*e1(3,i))
            temp1jyz7 =  -urr7*scale1i*rij(2)*rij(3)
            temp1jyz = temp1jyz5 + temp1jyz7
            temp1iyz5 =  urr5*(rij(3)*e1(2,jbis)+rij(2)*e1(3,jbis))
            temp1iyz7 =  -urr7*scale1j*rij(2)*rij(3)
            temp1iyz = temp1iyz5 + temp1iyz7

            temp2jyz5 =  urr5*(rij(3)*e2(2,i)+rij(2)*e2(3,i))
            temp2jyz7 =  -urr7*scale2i*rij(2)*rij(3)
            temp2jyz = temp2jyz5 + temp2jyz7
            temp2iyz5 =  urr5*(rij(3)*e2(2,jbis)+rij(2)*e2(3,jbis))
            temp2iyz7 =  -urr7*scale2j*rij(2)*rij(3)
            temp2iyz = temp2iyz5 + temp2iyz7

            temp3jyz5 =  urr5*(rij(3)*e3(2,i)+rij(2)*e3(3,i))
            temp3jyz7 =  -urr7*scale3i*rij(2)*rij(3)
            temp3jyz = temp3jyz5 + temp3jyz7
            temp3iyz5 =  urr5*(rij(3)*e3(2,jbis)+rij(2)*e3(3,jbis))
            temp3iyz7 =  -urr7*scale3j*rij(2)*rij(3)
            temp3iyz = temp3iyz5 + temp3iyz7
c    
c      beta = 1
c
            depx =  
     $       d1(1,jbis)*temp1jxx + d1(2,jbis)*temp1jxy 
     $      + d1(3,jbis)*temp1jxz
     $      + d1(1,i)*temp1ixx + d1(2,i)*temp1ixy 
     $      + d1(3,i)*temp1ixz
     $      + d2(1,jbis)*temp2jxx + d2(2,jbis)*temp2jxy 
     $      + d2(3,jbis)*temp2jxz
     $      + d2(1,i)*temp2ixx + d2(2,i)*temp2ixy 
     $      + d2(3,i)*temp2ixz
     $      + d3(1,jbis)*temp3jxx + d3(2,jbis)*temp3jxy 
     $      + d3(3,jbis)*temp3jxz
     $      + d3(1,i)*temp3ixx + d3(2,i)*temp3ixy 
     $      + d3(3,i)*temp3ixz
            adTbcor(1,i) = adTbcor(1,i) + depx 
            adTbcor(1,jbis) = adTbcor(1,jbis) - depx 
c    
c      beta = 2
c
            depy =  
     $       d1(1,jbis)*temp1jxy + d1(2,jbis)*temp1jyy 
     $      + d1(3,jbis)*temp1jyz
     $      + d1(1,i)*temp1ixy + d1(2,i)*temp1iyy 
     $      + d1(3,i)*temp1iyz
     $      + d2(1,jbis)*temp2jxy + d2(2,jbis)*temp2jyy 
     $      + d2(3,jbis)*temp2jyz
     $      + d2(1,i)*temp2ixy + d2(2,i)*temp2iyy 
     $      + d2(3,i)*temp2iyz
     $      + d3(1,jbis)*temp3jxy + d3(2,jbis)*temp3jyy 
     $      + d3(3,jbis)*temp3jyz
     $      + d3(1,i)*temp3ixy + d3(2,i)*temp3iyy 
     $      + d3(3,i)*temp3iyz
            adTbcor(2,i) = adTbcor(2,i) + depy 
            adTbcor(2,jbis) = adTbcor(2,jbis) - depy 
c    
c      beta = 3
c
            depz =  
     $       d1(1,jbis)*temp1jxz + d1(2,jbis)*temp1jyz 
     $      + d1(3,jbis)*temp1jzz
     $      + d1(1,i)*temp1ixz + d1(2,i)*temp1iyz 
     $      + d1(3,i)*temp1izz
     $      + d2(1,jbis)*temp2jxz + d2(2,jbis)*temp2jyz 
     $      + d2(3,jbis)*temp2jzz
     $      + d2(1,i)*temp2ixz + d2(2,i)*temp2iyz 
     $      + d2(3,i)*temp2izz
     $      + d3(1,jbis)*temp3jxz + d3(2,jbis)*temp3jyz 
     $      + d3(3,jbis)*temp3jzz
     $      + d3(1,i)*temp3ixz + d3(2,i)*temp3iyz 
     $      + d3(3,i)*temp3izz
            adTbcor(3,i) = adTbcor(3,i) + depz 
            adTbcor(3,jbis) = adTbcor(3,jbis) - depz 

         end do !jj


         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = 1.0d0
            uscale(ip11(j,iglob)) = 1.0d0
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = 1.0d0
            uscale(ip12(j,iglob)) = 1.0d0
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = 1.0d0
            uscale(ip13(j,iglob)) = 1.0d0
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = 1.0d0
            uscale(ip14(j,iglob)) = 1.0d0
         end do
         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = 1.0d0
         end do
      end do

      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      adme_self(:,1,1:npoleloc) = term*a1(:,1:npoleloc)
      adme_self(:,2,1:npoleloc) = term*a2(:,1:npoleloc)

      ade  = - ade_cor  ! THIS -1 EVERYWHERE IS WEIRD
      adme = - adme_cor + adme_self ! But hey... 
      adte = - adte_cor !  ...it works.
      adtb = - adtbcor!

      deallocate(dscale, pscale, uscale)
      return
      end

c
      subroutine scalderfieldzmat3( a1, a2, d1, e1, d2, e2,
     $                             ade, adme, adte, adtb)

      use atmlst
      use atoms
      use bound
      use chgpot
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
      use potent
      use shunt
      use virial
      implicit none


      real*8, dimension(3,npolebloc), intent(in) :: a1, a2
      real*8, dimension(3,npolebloc), intent(in) :: e1, e2, d1, d2
      real*8, dimension(3,npolebloc), intent(out) :: adtb
      real*8, dimension(3,2,npolebloc), intent(out) :: adme, ade
      real*8, dimension(3,3,2,npolebloc), intent(out) :: adte
      
      real*8, allocatable, dimension(:,:,:) ::adE_cor,admE_cor,admE_self
      real*8, allocatable, dimension(:,:,:,:) ::adtE_cor




      integer :: ii, iipole, iglob,  i, jjpole, jglob,
     $           jjj, jbis, k, j

      real*8 :: pti, pdi, rrij3, rrij5, rrij7, rrij9,
     $         urr5, urr7, 
     $scale3, scale5, 
     $          scale7, scale9, dsc3, dsc5, dsc7, dsc9, psc3, psc5,
     $          psc7, psc9, usc5, usc7, scalmuir, scalmujr, 
     $quadterm, quadterm_bis, 
     $          ci, cj, bfac, invrij2, alsq2n, cutoff2, exp2a, ralpha,
     $          d, rij2, damp, expdamp, pgamma, alsq2, term

      real*8                :: bn(0:4) 
      real*8, dimension(2)  :: srr3, srr5, srr7, srr9
      real*8, dimension(3)  :: rij, thetajr, thetair, di, dj
      real*8, dimension(3,3):: qi, qj

      real*8, allocatable, dimension(:) :: dscale, pscale, uscale
      real*8, allocatable, dimension(:,:) :: adtbcor
      

      real*8:: a1xi,a1yi,a1zi,a1xj,a1yj,a1zj
      real*8:: a2xi,a2yi,a2zi,a2xj,a2yj,a2zj
      real*8,dimension(2):: term1,term2,term3,term4,term5,term6,term7
      real*8,dimension(2):: dexdxi,dexdxj,deydyi,deydyj,dezdzi,dezdzj
      real*8,dimension(2):: dexdyi,dexdyj,dexdzi,dexdzj,deydzi,deydzj
      real*8,dimension(2):: demxdx,demxdy,demxdz,demydy,demydz,demzdz
      real*8,dimension(2):: detxdx,detxdy,detxdz,detydy,detydz,detzdz
      real*8 :: depx,depy,depz,dempx,dempy,dempz,detx,dety,detz
      real*8 :: scale1i,scale1j,scale2i,scale2j
      real*8::temp1jxx5,temp1jxy5,temp1jxz5,temp1jyy5,temp1jyz5,
     $     temp1jzz5
      real*8::temp1jxx7,temp1jxy7,temp1jxz7,temp1jyy7,temp1jyz7,
     $     temp1jzz7
      real*8::temp1jxx,temp1jxy,temp1jxz,temp1jyy,temp1jyz,temp1jzz
      real*8::temp1ixx5,temp1ixy5,temp1ixz5,temp1iyy5,temp1iyz5,
     $     temp1izz5
      real*8::temp1ixx7,temp1ixy7,temp1ixz7,temp1iyy7,temp1iyz7,
     $     temp1izz7
      real*8::temp1ixx,temp1ixy,temp1ixz,temp1iyy,temp1iyz,temp1izz

      real*8::temp2jxx5,temp2jxy5,temp2jxz5,temp2jyy5,temp2jyz5,
     $     temp2jzz5
      real*8::temp2jxx7,temp2jxy7,temp2jxz7,temp2jyy7,temp2jyz7,
     $     temp2jzz7
      real*8::temp2jxx,temp2jxy,temp2jxz,temp2jyy,temp2jyz,temp2jzz
      real*8::temp2ixx5,temp2ixy5,temp2ixz5,temp2iyy5,temp2iyz5,
     $     temp2izz5
      real*8::temp2ixx7,temp2ixy7,temp2ixz7,temp2iyy7,temp2iyz7,
     $     temp2izz7
      real*8::temp2ixx,temp2ixy,temp2ixz,temp2iyy,temp2iyz,temp2izz












      real*8 :: f

      character(10) :: mode

 1000 format(' Warning, system moved too much since last neighbor list
     $   update, try lowering nlupdate')

      
      f = electric/dielec

      ade = 0d0
      adme = 0d0
      adte = 0d0
      adtb = 0d0

      mode = 'EWALD'
      if (use_polarshortreal) mode = 'SHORTEWALD'
      call switch (mode)
      cutoff2 = cut2

      allocate(dscale(n), pscale(n), uscale(n))
      allocate(adtbcor(3,npolebloc))
      allocate(adE_cor(3,2,npolebloc),admE_cor(3,2,npolebloc))
      allocate(admE_self(3,2,npolebloc),adtE_cor(3,3,2,npolebloc))
      adme_cor = 0d0
      adme_self = 0d0
      adte_cor = 0d0
      ade_cor = 0d0
      dscale = 1d0
      pscale = 1d0
      uscale = 1.0d0
      adtbcor = 0d0

      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = poleloc(iglob)
         if ((i.eq.0).or.(i.gt.nbloc)) then
            write(iout,1000)
            cycle
         end if
         pdi = pdamp(iipole)
         pti = thole(iipole)

         ci      = rpole(1,iipole)
         di(1)   = rpole(2,iipole)
         di(2)   = rpole(3,iipole)
         di(3)   = rpole(4,iipole)
         qi(1,1) = rpole(5, iipole)
         qi(2,1) = rpole(6, iipole)
         qi(3,1) = rpole(7, iipole)
         qi(1,2) = rpole(6, iipole)
         qi(2,2) = rpole(9, iipole)
         qi(3,2) = rpole(10,iipole)
         qi(1,3) = rpole(7, iipole)
         qi(2,3) = rpole(10,iipole)
         qi(3,3) = rpole(13,iipole)

         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = p2scale
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = p3scale
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = p4scale
            do k = 1, np11(iglob)
               if (i14(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i14(j,iglob)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = p5scale
         end do
         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = d1scale
            uscale(ip11(j,iglob)) = u1scale
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = d2scale
            uscale(ip12(j,iglob)) = u2scale
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = d3scale
            uscale(ip13(j,iglob)) = u3scale
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = d4scale
            uscale(ip14(j,iglob)) = u4scale
         end do

         do jjj = 1, nelst(ii)
            jjpole = elst(jjj,ii)
            jglob = ipole(jjpole)
            jbis = poleloc(jjpole)
            if (jbis.eq.0) then
               write(iout,1000)
               cycle
            end if

            !Distances
            rij(1) = x(jglob) - x(iglob)
            rij(2) = y(jglob) - y(iglob)
            rij(3) = z(jglob) - z(iglob)
            call image(rij(1), rij(2), rij(3))
            rij2 = dot_product(rij,rij)
            if (rij2 .gt. cutoff2) cycle
            d  = sqrt(rij2)
            invrij2 = 1d0/rij2
            rrij3 = invrij2/d
            rrij5 = 3d0*rrij3*invrij2
            rrij7 = 5d0*rrij5*invrij2
            rrij9 = 7d0*rrij7*invrij2

            !Multipoles
            cj      = rpole(1,jjpole)
            dj(1)   = rpole(2,jjpole)
            dj(2)   = rpole(3,jjpole)
            dj(3)   = rpole(4,jjpole)
            qj(1,1) = rpole(5, jjpole)
            qj(2,1) = rpole(6, jjpole)
            qj(3,1) = rpole(7, jjpole)
            qj(1,2) = rpole(6, jjpole)
            qj(2,2) = rpole(9, jjpole)
            qj(3,2) = rpole(10,jjpole)
            qj(1,3) = rpole(7, jjpole)
            qj(2,3) = rpole(10,jjpole)
            qj(3,3) = rpole(13,jjpole)

            scalmujr = dot_product(dj,rij)!sprod(3,dj,rij)
            scalmuir = dot_product(di,rij)!sprod(3,di,rij)
            thetajr(1) = dot_product(qj(1,:),rij)
            thetajr(2) = dot_product(qj(2,:),rij)
            thetajr(3) = dot_product(qj(3,:),rij)
            thetair(1) = dot_product(qi(1,:),rij)
            thetair(2) = dot_product(qi(2,:),rij)
            thetair(3) = dot_product(qi(3,:),rij)
            quadterm     = dot_product(thetajr,rij)!sprod(3,thetajr,rij)
            quadterm_bis = dot_product(thetair,rij)!sprod(3,thetair,rij)

            scale1j = dot_product(rij,e1(:,jbis))
            scale1i = dot_product(rij,e1(:,i))
            scale2j = dot_product(rij,e2(:,jbis))
            scale2i = dot_product(rij,e2(:,i))

            a1xi = a1(1,i)
            a1yi = a1(2,i)
            a1zi = a1(3,i)
            a1xj = a1(1,jbis)
            a1yj = a1(2,jbis)
            a1zj = a1(3,jbis)
            a2xi = a2(1,i)
            a2yi = a2(2,i)
            a2zi = a2(3,i)
            a2xj = a2(1,jbis)
            a2yj = a2(2,jbis)
            a2zj = a2(3,jbis)

            !Errfunc damping
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0d0 * aewald**2
            alsq2n = 0.0d0
            if (aewald .gt. 0.0d0)
     &        alsq2n = 1.0d0 / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 4 
              bfac = dble(j+j-1)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a)*invrij2
            end do

            !Scalings
            damp = pdi*pdamp(jjpole)
            scale3 = 1d0
            scale5 = 1d0
            scale7 = 1d0
            scale9 = 1d0
            if (damp.ne.0d0) then
               pgamma = min(pti,thole(jjpole))
               damp = -pgamma*(d/damp)**3
               if (damp .gt. -50d0) then
                  expdamp = exp(damp)
                  scale3 = scale3 * (1.0d0 - expdamp)
                  scale5 = scale5 * (1.0d0 - expdamp*(1.0d0 - damp))
                  scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                  scale9 = 1d0 - expdamp*(1d0 - damp + (18d0*damp**2 -
     $                     9d0*damp**3)/35d0)
               end if
            end if
            dsc3 = scale3*dscale(jglob)
            dsc5 = scale5*dscale(jglob)
            dsc7 = scale7*dscale(jglob)
            dsc9 = scale9*dscale(jglob)
            psc3 = scale3*pscale(jglob)
            psc5 = scale5*pscale(jglob)
            psc7 = scale7*pscale(jglob)
            psc9 = scale9*pscale(jglob)
            !mat
            usc5 = scale5*uscale(jglob)
            usc7 = scale7*uscale(jglob)
            srr3(1) = (1d0 - dsc3)*rrij3 - bn(1)
            srr3(2) = (1d0 - psc3)*rrij3 - bn(1)
            srr5(1) = (1d0 - dsc5)*rrij5 - bn(2)
            srr5(2) = (1d0 - psc5)*rrij5 - bn(2)
            srr7(1) = (1d0 - dsc7)*rrij7 - bn(3)
            srr7(2) = (1d0 - psc7)*rrij7 - bn(3)
            srr9(1) = (1d0 - dsc9)*rrij9 - bn(4)
            srr9(2) = (1d0 - psc9)*rrij9 - bn(4)
            !mat
            urr5 = (1d0 - usc5)*rrij5 - bn(2)
            urr7 = (1d0 - usc7)*rrij7 - bn(3)

            !Actual field calculations ! 



c dexdxj (alpha=beta=1)
            term1 = rij(1)*rij(1)*srr5 - srr3 
            term2 = 2d0*srr5*rij(1) 
            term3 = -srr7*rij(1)*rij(1) + srr5
            term4 = -4d0*rij(1)*srr7
            term5 = srr9*rij(1)*rij(1) - srr7

            dexdxj = cj*term1 + dj(1)*term2
     $      + scalmujr*term3 + 2d0*qj(1,1)*srr5
     $      + thetajr(1)*term4 + quadterm*term5

            dexdxi = ci*term1 - di(1)*term2
     $      - scalmuir*term3 + 2d0*qi(1,1)*srr5
     $      + thetair(1)*term4 + quadterm_bis*term5
            
            demxdx = term1
            detxdx = 2d0*srr5-srr7*rij(1)*rij(1)

c deydyj (alpha=beta=2)
            term1 = rij(2)*rij(2)*srr5 - srr3 
            term2 = 2d0*srr5*rij(2) 
            term3 = -srr7*rij(2)*rij(2) + srr5
            term4 = -4d0*rij(2)*srr7
            term5 = srr9*rij(2)*rij(2) - srr7

            deydyj = cj*term1 + dj(2)*term2
     $      + scalmujr*term3 + 2d0*qj(2,2)*srr5
     $      + thetajr(2)*term4 + quadterm*term5

            deydyi = ci*term1 - di(2)*term2
     $      - scalmuir*term3 + 2d0*qi(2,2)*srr5
     $      + thetair(2)*term4 + quadterm_bis*term5

            demydy = term1
            detydy = 2d0*srr5-srr7*rij(2)*rij(2)

c dexdyj (alpha=1 beta=2)
            term1 = rij(1)*rij(2)*srr5 
            term2 = srr5*rij(2) 
            term3 = srr5*rij(1) 
            term4 = -srr7*rij(1)*rij(2)
            term5 = -2d0*rij(2)*srr7
            term6 = -2d0*rij(1)*srr7
            term7 = srr9*rij(1)*rij(2)

            dexdyj = cj*term1 + dj(1)*term2 + dj(2)*term3
     $      + scalmujr*term4 + 2d0*qj(1,2)*srr5
     $      + thetajr(1)*term5 +  thetajr(2)*term6
     $      + quadterm*term7

            dexdyi = ci*term1 - di(1)*term2 - di(2)*term3
     $      - scalmuir*term4 + 2d0*qi(1,2)*srr5
     $      + thetair(1)*term5 +  thetair(2)*term6
     $      + quadterm_bis*term7

            demxdy = term1
            detxdy = -srr7*rij(1)*rij(2)


c dezdzj (alpha3=beta=3)
            term1 = rij(3)*rij(3)*srr5 - srr3 
            term2 = 2d0*srr5*rij(3) 
            term3 = -srr7*rij(3)*rij(3) + srr5
            term4 = -4d0*rij(3)*srr7
            term5 = srr9*rij(3)*rij(3) - srr7

            dezdzj = cj*term1 + dj(3)*term2
     $      + scalmujr*term3 + 2d0*qj(3,3)*srr5
     $      + thetajr(3)*term4 + quadterm*term5

            dezdzi = ci*term1 - di(3)*term2
     $      - scalmuir*term3 + 2d0*qi(3,3)*srr5
     $      + thetair(3)*term4 + quadterm_bis*term5

            demzdz = term1
            detzdz = 2d0*srr5-srr7*rij(3)*rij(3)

c dexdzj (alpha=1 beta=3)
            term1 = rij(1)*rij(3)*srr5 
            term2 = srr5*rij(3) 
            term3 = srr5*rij(1) 
            term4 = -srr7*rij(1)*rij(3)
            term5 = -2d0*rij(3)*srr7
            term6 = -2d0*rij(1)*srr7
            term7 = srr9*rij(1)*rij(3)

            dexdzj = cj*term1 + dj(1)*term2 + dj(3)*term3
     $      + scalmujr*term4 + 2d0*qj(1,3)*srr5
     $      + thetajr(1)*term5 +  thetajr(3)*term6
     $      + quadterm*term7

            dexdzi = ci*term1 - di(1)*term2 - di(3)*term3
     $      - scalmuir*term4 + 2d0*qi(1,3)*srr5
     $      + thetair(1)*term5 +  thetair(3)*term6
     $      + quadterm_bis*term7

            demxdz = term1
            detxdz = -srr7*rij(1)*rij(3)


c deydzj (alpha=2 beta=3)
            term1 = rij(2)*rij(3)*srr5 
            term2 = srr5*rij(3) 
            term3 = srr5*rij(2) 
            term4 = -srr7*rij(2)*rij(3)
            term5 = -2d0*rij(3)*srr7
            term6 = -2d0*rij(2)*srr7
            term7 = srr9*rij(2)*rij(3)

            deydzj = cj*term1 + dj(2)*term2 + dj(3)*term3
     $      + scalmujr*term4 + 2d0*qj(2,3)*srr5
     $      + thetajr(2)*term5 +  thetajr(3)*term6
     $      + quadterm*term7

            deydzi = ci*term1 - di(2)*term2 - di(3)*term3
     $      - scalmuir*term4 + 2d0*qi(2,3)*srr5
     $      + thetair(2)*term5 +  thetair(3)*term6
     $      + quadterm_bis*term7

            demydz = term1
            detydz = -srr7*rij(2)*rij(3)

            depx = -a1xi*dexdxj(1) + a1xj*dexdxi(1)
     $     -a1yi*dexdyj(1) + a1yj*dexdyi(1)
     $     -a1zi*dexdzj(1) + a1zj*dexdzi(1)
            depy = -a1xi*dexdyj(1) + a1xj*dexdyi(1)
     $     -a1yi*deydyj(1) + a1yj*deydyi(1)
     $     -a1zi*deydzj(1) + a1zj*deydzi(1)
            depz = -a1xi*dexdzj(1) + a1xj*dexdzi(1)
     $     -a1yi*deydzj(1) + a1yj*deydzi(1)
     $     -a1zi*dezdzj(1) + a1zj*dezdzi(1)

            ade_cor(1,1,i) = ade_cor(1,1,i) + depx
            ade_cor(2,1,i) = ade_cor(2,1,i) + depy
            ade_cor(3,1,i) = ade_cor(3,1,i) + depz

            ade_cor(1,1,jbis) = ade_cor(1,1,jbis) - depx
            ade_cor(2,1,jbis) = ade_cor(2,1,jbis) - depy
            ade_cor(3,1,jbis) = ade_cor(3,1,jbis) - depz

            depx = -a2xi*dexdxj(2) + a2xj*dexdxi(2)
     $     -a2yi*dexdyj(2) + a2yj*dexdyi(2)
     $     -a2zi*dexdzj(2) + a2zj*dexdzi(2)
            depy = -a2xi*dexdyj(2) + a2xj*dexdyi(2)
     $     -a2yi*deydyj(2) + a2yj*deydyi(2)
     $     -a2zi*deydzj(2) + a2zj*deydzi(2)
            depz = -a2xi*dexdzj(2) + a2xj*dexdzi(2)
     $     -a2yi*deydzj(2) + a2yj*deydzi(2)
     $     -a2zi*dezdzj(2) + a2zj*dezdzi(2)

            ade_cor(1,2,i) = ade_cor(1,2,i) + depx
            ade_cor(2,2,i) = ade_cor(2,2,i) + depy
            ade_cor(3,2,i) = ade_cor(3,2,i) + depz
cc
cc    increment virial
cc
c            vxx = rij(1) * 0.5*f*depx
c            vxy = rij(2) * 0.5*f*depx
c            vxz = rij(3) * 0.5*f*depx
c            vyy = rij(2) * 0.5*f*depy
c            vyz = rij(3) * 0.5*f*depy
c            vzz = rij(3) * 0.5*f*depz
c            vir(1,1) = vir(1,1) + vxx
c            vir(2,1) = vir(2,1) + vxy
c            vir(3,1) = vir(3,1) + vxz
c            vir(1,2) = vir(1,2) + vxy
c            vir(2,2) = vir(2,2) + vyy
c            vir(3,2) = vir(3,2) + vyz
c            vir(1,3) = vir(1,3) + vxz
c            vir(2,3) = vir(2,3) + vyz
c            vir(3,3) = vir(3,3) + vzz

            ade_cor(1,2,jbis) = ade_cor(1,2,jbis) - depx
            ade_cor(2,2,jbis) = ade_cor(2,2,jbis) - depy
            ade_cor(3,2,jbis) = ade_cor(3,2,jbis) - depz

            dempx = a1xj*demxdx(1) + a1yj*demxdy(1)
     $     + a1zj*demxdz(1)
            dempy = a1xj*demxdy(1) + a1yj*demydy(1)
     $     + a1zj*demydz(1)
            dempz = a1xj*demxdz(1) + a1yj*demydz(1)
     $     + a1zj*demzdz(1)

            adme_cor(1,1,i) = adme_cor(1,1,i) + dempx
            adme_cor(2,1,i) = adme_cor(2,1,i) + dempy
            adme_cor(3,1,i) = adme_cor(3,1,i) + dempz

            dempx = a1xi*demxdx(1) + a1yi*demxdy(1)
     $     + a1zi*demxdz(1)
            dempy = a1xi*demxdy(1) + a1yi*demydy(1)
     $     + a1zi*demydz(1)
            dempz = a1xi*demxdz(1) + a1yi*demydz(1)
     $     + a1zi*demzdz(1)

            adme_cor(1,1,jbis) = adme_cor(1,1,jbis) + dempx
            adme_cor(2,1,jbis) = adme_cor(2,1,jbis) + dempy
            adme_cor(3,1,jbis) = adme_cor(3,1,jbis) + dempz

            dempx = a2xj*demxdx(2) + a2yj*demxdy(2)
     $     + a2zj*demxdz(2)
            dempy = a2xj*demxdy(2) + a2yj*demydy(2)
     $     + a2zj*demydz(2)
            dempz = a2xj*demxdz(2) + a2yj*demydz(2)
     $     + a2zj*demzdz(2)

            adme_cor(1,2,i) = adme_cor(1,2,i) + dempx
            adme_cor(2,2,i) = adme_cor(2,2,i) + dempy
            adme_cor(3,2,i) = adme_cor(3,2,i) + dempz

            dempx = a2xi*demxdx(2) + a2yi*demxdy(2)
     $     + a2zi*demxdz(2)
            dempy = a2xi*demxdy(2) + a2yi*demydy(2)
     $     + a2zi*demydz(2)
            dempz = a2xi*demxdz(2) + a2yi*demydz(2)
     $     + a2zi*demzdz(2)

            adme_cor(1,2,jbis) = adme_cor(1,2,jbis) + dempx
            adme_cor(2,2,jbis) = adme_cor(2,2,jbis) + dempy
            adme_cor(3,2,jbis) = adme_cor(3,2,jbis) + dempz

            detx = -a1xj*detxdx(1)-a1yj*detxdy(1)
     $     -a1zj*detxdz(1) 
            dety = -a1xj*detxdy(1)-a1yj*detydy(1)
     $     -a1zj*detydz(1) 
            detz = -a1xj*detxdz(1)-a1yj*detydz(1)
     $     -a1zj*detzdz(1) 

            adtE_cor(1,:,1,i) = adtE_cor(1,:,1,i) + detx*rij(:)
            adtE_cor(2,:,1,i) = adtE_cor(2,:,1,i) + dety*rij(:)
            adtE_cor(3,:,1,i) = adtE_cor(3,:,1,i) + detz*rij(:)


            detx = a1xi*detxdx(1)+a1yi*detxdy(1)
     $      + a1zi*detxdz(1) 
            dety = a1xi*detxdy(1)+a1yi*detydy(1)
     $      + a1zi*detydz(1) 
            detz = a1xi*detxdz(1)+a1yi*detydz(1)
     $      + a1zi*detzdz(1) 

            adtE_cor(1,:,1,jbis) = adtE_cor(1,:,1,jbis) + detx*rij(:)
            adtE_cor(2,:,1,jbis) = adtE_cor(2,:,1,jbis) + dety*rij(:)
            adtE_cor(3,:,1,jbis) = adtE_cor(3,:,1,jbis) + detz*rij(:)

            detx = -a2xj*detxdx(2)-a2yj*detxdy(2)
     $     -a2zj*detxdz(2) 
            dety = -a2xj*detxdy(2)-a2yj*detydy(2)
     $     -a2zj*detydz(2) 
            detz = -a2xj*detxdz(2)-a2yj*detydz(2)
     $     -a2zj*detzdz(2) 

            adtE_cor(1,:,2,i) = adtE_cor(1,:,2,i) + detx*rij(:)
            adtE_cor(2,:,2,i) = adtE_cor(2,:,2,i) + dety*rij(:)
            adtE_cor(3,:,2,i) = adtE_cor(3,:,2,i) + detz*rij(:)

            detx = a2xi*detxdx(2)+a2yi*detxdy(2)
     $      + a2zi*detxdz(2) 
            dety = a2xi*detxdy(2)+a2yi*detydy(2)
     $      + a2zi*detydz(2) 
            detz = a2xi*detxdz(2)+a2yi*detydz(2)
     $      + a2zi*detzdz(2) 

            adtE_cor(1,:,2,jbis) = adtE_cor(1,:,2,jbis) + detx*rij(:)
            adtE_cor(2,:,2,jbis) = adtE_cor(2,:,2,jbis) + dety*rij(:)
            adtE_cor(3,:,2,jbis) = adtE_cor(3,:,2,jbis) + detz*rij(:)
c
c      also compute gradient of field of dipoles "e"
c
c      alpha = beta = 1
c
            temp1jxx5 =  urr5*(scale1i + 2*rij(1)*e1(1,i))
            temp1jxx7 =  -urr7*scale1i*rij(1)*rij(1)
            temp1jxx = temp1jxx5 + temp1jxx7
            temp1ixx5 =  urr5*(scale1j + 2*rij(1)*e1(1,jbis))
            temp1ixx7 =  -urr7*scale1j*rij(1)*rij(1)
            temp1ixx = temp1ixx5 + temp1ixx7

            temp2jxx5 =  urr5*(scale2i + 2*rij(1)*e2(1,i))
            temp2jxx7 =  -urr7*scale2i*rij(1)*rij(1)
            temp2jxx = temp2jxx5 + temp2jxx7
            temp2ixx5 =  urr5*(scale2j + 2*rij(1)*e2(1,jbis))
            temp2ixx7 =  -urr7*scale2j*rij(1)*rij(1)
            temp2ixx = temp2ixx5 + temp2ixx7

c      alpha = 2, beta = 1
c
            temp1jxy5 =  urr5*(rij(2)*e1(1,i)+rij(1)*e1(2,i))
            temp1jxy7 =  -urr7*scale1i*rij(2)*rij(1)
            temp1jxy = temp1jxy5 + temp1jxy7
            temp1ixy5 =  urr5*(rij(2)*e1(1,jbis)+rij(1)*e1(2,jbis))
            temp1ixy7 =  -urr7*scale1j*rij(2)*rij(1)
            temp1ixy = temp1ixy5 + temp1ixy7

            temp2jxy5 =  urr5*(rij(2)*e2(1,i)+rij(1)*e2(2,i))
            temp2jxy7 =  -urr7*scale2i*rij(2)*rij(1)
            temp2jxy = temp2jxy5 + temp2jxy7
            temp2ixy5 =  urr5*(rij(2)*e2(1,jbis)+rij(1)*e2(2,jbis))
            temp2ixy7 =  -urr7*scale2j*rij(2)*rij(1)
            temp2ixy = temp2ixy5 + temp2ixy7

c
c      alpha = 3, beta = 1
c
            temp1jxz5 =  urr5*(rij(3)*e1(1,i)+rij(1)*e1(3,i))
            temp1jxz7 =  -urr7*scale1i*rij(3)*rij(1)
            temp1jxz = temp1jxz5 + temp1jxz7
            temp1ixz5 =  urr5*(rij(3)*e1(1,jbis)+rij(1)*e1(3,jbis))
            temp1ixz7 =  -urr7*scale1j*rij(3)*rij(1)
            temp1ixz = temp1ixz5 + temp1ixz7

            temp2jxz5 =  urr5*(rij(3)*e2(1,i)+rij(1)*e2(3,i))
            temp2jxz7 =  -urr7*scale2i*rij(3)*rij(1)
            temp2jxz = temp2jxz5 + temp2jxz7
            temp2ixz5 =  urr5*(rij(3)*e2(1,jbis)+rij(1)*e2(3,jbis))
            temp2ixz7 =  -urr7*scale2j*rij(3)*rij(1)
            temp2ixz = temp2ixz5 + temp2ixz7

c
c      alpha = beta = 2
c
            temp1jyy5 =  urr5*(scale1i + 2*rij(2)*e1(2,i))
            temp1jyy7 =  -urr7*scale1i*rij(2)*rij(2)
            temp1jyy = temp1jyy5 + temp1jyy7
            temp1iyy5 =  urr5*(scale1j + 2*rij(2)*e1(2,jbis))
            temp1iyy7 =  -urr7*scale1j*rij(2)*rij(2)
            temp1iyy = temp1iyy5 + temp1iyy7

            temp2jyy5 =  urr5*(scale2i + 2*rij(2)*e2(2,i))
            temp2jyy7 =  -urr7*scale2i*rij(2)*rij(2)
            temp2jyy = temp2jyy5 + temp2jyy7
            temp2iyy5 =  urr5*(scale2j + 2*rij(2)*e2(2,jbis))
            temp2iyy7 =  -urr7*scale2j*rij(2)*rij(2)
            temp2iyy = temp2iyy5 + temp2iyy7

c
c      alpha = beta = 3
c
            temp1jzz5 =  urr5*(scale1i + 2*rij(3)*e1(3,i))
            temp1jzz7 =  -urr7*scale1i*rij(3)*rij(3)
            temp1jzz = temp1jzz5 + temp1jzz7
            temp1izz5 =  urr5*(scale1j + 2*rij(3)*e1(3,jbis))
            temp1izz7 =  -urr7*scale1j*rij(3)*rij(3)
            temp1izz = temp1izz5 + temp1izz7

            temp2jzz5 =  urr5*(scale2i + 2*rij(3)*e2(3,i))
            temp2jzz7 =  -urr7*scale2i*rij(3)*rij(3)
            temp2jzz = temp2jzz5 + temp2jzz7
            temp2izz5 =  urr5*(scale2j + 2*rij(3)*e2(3,jbis))
            temp2izz7 =  -urr7*scale2j*rij(3)*rij(3)
            temp2izz = temp2izz5 + temp2izz7

c
c      alpha = 2, beta = 3
c
            temp1jyz5 =  urr5*(rij(3)*e1(2,i)+rij(2)*e1(3,i))
            temp1jyz7 =  -urr7*scale1i*rij(2)*rij(3)
            temp1jyz = temp1jyz5 + temp1jyz7
            temp1iyz5 =  urr5*(rij(3)*e1(2,jbis)+rij(2)*e1(3,jbis))
            temp1iyz7 =  -urr7*scale1j*rij(2)*rij(3)
            temp1iyz = temp1iyz5 + temp1iyz7

            temp2jyz5 =  urr5*(rij(3)*e2(2,i)+rij(2)*e2(3,i))
            temp2jyz7 =  -urr7*scale2i*rij(2)*rij(3)
            temp2jyz = temp2jyz5 + temp2jyz7
            temp2iyz5 =  urr5*(rij(3)*e2(2,jbis)+rij(2)*e2(3,jbis))
            temp2iyz7 =  -urr7*scale2j*rij(2)*rij(3)
            temp2iyz = temp2iyz5 + temp2iyz7

c    
c      beta = 1
c
            depx =  
     $       d1(1,jbis)*temp1jxx + d1(2,jbis)*temp1jxy 
     $      + d1(3,jbis)*temp1jxz
     $      + d1(1,i)*temp1ixx + d1(2,i)*temp1ixy 
     $      + d1(3,i)*temp1ixz
     $      + d2(1,jbis)*temp2jxx + d2(2,jbis)*temp2jxy 
     $      + d2(3,jbis)*temp2jxz
     $      + d2(1,i)*temp2ixx + d2(2,i)*temp2ixy 
     $      + d2(3,i)*temp2ixz
            adTbcor(1,i) = adTbcor(1,i) + depx 
            adTbcor(1,jbis) = adTbcor(1,jbis) - depx 
c    
c      beta = 2
c
            depy =  
     $       d1(1,jbis)*temp1jxy + d1(2,jbis)*temp1jyy 
     $      + d1(3,jbis)*temp1jyz
     $      + d1(1,i)*temp1ixy + d1(2,i)*temp1iyy 
     $      + d1(3,i)*temp1iyz
     $      + d2(1,jbis)*temp2jxy + d2(2,jbis)*temp2jyy 
     $      + d2(3,jbis)*temp2jyz
     $      + d2(1,i)*temp2ixy + d2(2,i)*temp2iyy 
     $      + d2(3,i)*temp2iyz
            adTbcor(2,i) = adTbcor(2,i) + depy 
            adTbcor(2,jbis) = adTbcor(2,jbis) - depy 
c    
c      beta = 3
c
            depz =  
     $       d1(1,jbis)*temp1jxz + d1(2,jbis)*temp1jyz 
     $      + d1(3,jbis)*temp1jzz
     $      + d1(1,i)*temp1ixz + d1(2,i)*temp1iyz 
     $      + d1(3,i)*temp1izz
     $      + d2(1,jbis)*temp2jxz + d2(2,jbis)*temp2jyz 
     $      + d2(3,jbis)*temp2jzz
     $      + d2(1,i)*temp2ixz + d2(2,i)*temp2iyz 
     $      + d2(3,i)*temp2izz
            adTbcor(3,i) = adTbcor(3,i) + depz 
            adTbcor(3,jbis) = adTbcor(3,jbis) - depz 
cc
cc          increment virial          
cc
c            vxx = rij(1) * 0.5*f*depx
c            vxy = rij(2) * 0.5*f*depx
c            vxz = rij(3) * 0.5*f*depx
c            vyy = rij(2) * 0.5*f*depy
c            vyz = rij(3) * 0.5*f*depy
c            vzz = rij(3) * 0.5*f*depz
c            vir(1,1) = vir(1,1) + vxx
c            vir(2,1) = vir(2,1) + vxy
c            vir(3,1) = vir(3,1) + vxz
c            vir(1,2) = vir(1,2) + vxy
c            vir(2,2) = vir(2,2) + vyy
c            vir(3,2) = vir(3,2) + vyz
c            vir(1,3) = vir(1,3) + vxz
c            vir(2,3) = vir(2,3) + vyz
c            vir(3,3) = vir(3,3) + vzz
 

         end do !jj


         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = 1.0d0
            uscale(ip11(j,iglob)) = 1.0d0
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = 1.0d0
            uscale(ip12(j,iglob)) = 1.0d0
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = 1.0d0
            uscale(ip13(j,iglob)) = 1.0d0
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = 1.0d0
            uscale(ip14(j,iglob)) = 1.0d0
         end do
         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = 1.0d0
         end do
      end do

      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      adme_self(:,1,1:npoleloc) = term*a1(:,1:npoleloc)
      adme_self(:,2,1:npoleloc) = term*a2(:,1:npoleloc)

      ade  = - ade_cor  ! THIS -1 EVERYWHERE IS WEIRD
      adme = - adme_cor + adme_self ! But hey... 
      adte = - adte_cor !  ...it works.
      adtb = - adtbcor!



      deallocate(dscale, pscale, uscale)
      return
      end

!===================================================
!     sub scalderfieldzmatrec
!===================================================
! No scaling in E_recip.
! Given :
! - a pair of  vectors a and b, as well ad their 
! potential and subsequent derivatives fphiA, fphiB, 
! returns
! < a, dE_recip/dr> + < b, dE_recip/dr> = aderec
! < a, dE_recip/dmu > + < b, dE_recip/dmu > = admerec
! < a, dE_recip/dtheta > + < a, dE_recip/dtheta > = adterec
! 
! - three pairs Di, Ei, along with the potential and 
! successive derivatives arising fphiAi and fphiBi, 
! computes the reciprocal part of < A, T' B > 
! A and B should be given in cartesian coordinates
!
! Opti : uses no fft
! /!\  Finishes with an extra -1 factor !
!
      subroutine scalderfieldzmatrec6(a,abis,fphia,b,bbis,fphib, 
     $            d1,d1bis,fphid1,e1,e1bis,fphie1,d2,d2bis,fphid2, 
     $            e2,e2bis,fphie2,d3,d3bis,fphid3,e3,e3bis,fphie3, 
     $            aderec,admerec,adterec,adtb)
      use atmlst
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
      use potent
      use mpi
      implicit none



      real*8, intent(in), dimension(3,npolebloc) :: a, b, d1, d2, d3,
     $                                                e1, e2, e3
      real*8, intent(in), dimension(3,npolerecloc) :: abis,bbis,d1bis,
     $                                  d2bis,d3bis,e1bis,e2bis,e3bis
      real*8, intent(in), dimension(20,npolerecloc) :: fphia, fphib,
     $                  fphid1, fphid2, fphid3, fphie1, fphie2, fphie3
      real*8, intent(out), dimension(3,npolerecloc) :: aderec, adtb,
     $                                                      admerec
      real*8, intent(out), dimension(3,3,npolerecloc) :: adterec

      integer :: i, iipole, iglob, rankloc, j, k,
     $           l, ii, iloc
      integer, dimension(10) :: deriv1, deriv2, deriv3
      integer, dimension(3,10) :: derivi


      real*8, dimension(2) :: f1, f2, f3
      real*8, dimension(3) :: g1, g2, g3
      real*8, dimension(3,2) :: dm
      real*8, dimension(3,3,2) :: dt
      real*8, allocatable, dimension(:,:) :: cmp, fmp, a_frac, b_frac,
     $          fdervec, cdervec,
     $          fdervecdp, cdervecb, d1_frac, d2_frac, d3_frac,
     $          e1_frac, e2_frac, e3_frac, adtb_temp
      real*8, allocatable, dimension(:,:,:) :: frc, 
     $          frcdm, grc
      real*8, allocatable, dimension(:,:,:,:) :: frcdt
      

      rankloc = rank

      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
      derivi(1,:) = deriv1
      derivi(2,:) = deriv2
      derivi(3,:) = deriv3

      aderec = 0d0
      admerec = 0d0
      adterec = 0d0

      allocate(cmp(10,npolerecloc), fmp(10,npolerecloc))
      allocate(fdervec(20,npolerecloc), cdervec(20,npolerecloc),
     $fdervecdp(20,npolerecloc), 
     $          cdervecb(20,npolerecloc))
      allocate(frc(3,npolerecloc,2), grc(3,npolerecloc,3))
      allocate(a_frac(3,npolerecloc), b_frac(3,npolerecloc),
     $         adtb_temp(3,npolerecloc),
     $         d1_frac(3,npolerecloc), e1_frac(3,npolerecloc),
     $         d2_frac(3,npolerecloc), e2_frac(3,npolerecloc),
     $         d3_frac(3,npolerecloc), e3_frac(3,npolerecloc))
      allocate(frcdm(3,npolerecloc,2))
      allocate(frcdt(3,3,npolerecloc,2))

      cmp = 0d0
      fmp = 0d0
      frc = 0d0
      grc = 0d0
      adtb = 0d0
      
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         iloc = poleloc(iipole)
         cmp(1,i) =  rpole(1,iipole)
         cmp(2,i) =  rpole(2,iipole)
         cmp(3,i) =  rpole(3,iipole)
         cmp(4,i) =  rpole(4,iipole)
         cmp(5,i) =  rpole(5,iipole)
         cmp(6,i) =  rpole(9,iipole)
         cmp(7,i) =  rpole(13,iipole)
         cmp(8,i) =  rpole(6,iipole)  *2d0
         cmp(9,i) =  rpole(7,iipole)  *2d0
         cmp(10,i) = rpole(10,iipole)*2d0
         call cmp_to_fmp_site(cmp(1,i), fmp(1,i))
         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(abis(:,i), a_frac(:,i))
           call cart_to_frac_vec(bbis(:,i), b_frac(:,i))
         else
           call cart_to_frac_vec(a(:,iloc), a_frac(:,i))
           call cart_to_frac_vec(b(:,iloc), b_frac(:,i))
         end if
         call fphi_to_cphi_site(fphia(:,i), cdervec(:,i))
         call fphi_to_cphi_site(fphib(:,i), cdervecb(:,i))
         !mat
         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(d1bis(:,i), d1_frac(:,i))
           call cart_to_frac_vec(d2bis(:,i), d2_frac(:,i))
           call cart_to_frac_vec(d3bis(:,i), d3_frac(:,i))
           call cart_to_frac_vec(e1bis(:,i), e1_frac(:,i))
           call cart_to_frac_vec(e2bis(:,i), e2_frac(:,i))
           call cart_to_frac_vec(e3bis(:,i), e3_frac(:,i))
         else
           call cart_to_frac_vec(d1(:,iloc), d1_frac(:,i))
           call cart_to_frac_vec(d2(:,iloc), d2_frac(:,i))
           call cart_to_frac_vec(d3(:,iloc), d3_frac(:,i))
           call cart_to_frac_vec(e1(:,iloc), e1_frac(:,i))
           call cart_to_frac_vec(e2(:,iloc), e2_frac(:,i))
           call cart_to_frac_vec(e3(:,iloc), e3_frac(:,i))
         end if
      end do


      do i = 1, npolerecloc
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         !mat
         g1 = 0d0
         g2 = 0d0
         g3 = 0d0
         do k = 1, 10
            f1(1) = f1(1) + fmp(k,i)*fphia(deriv1(k),i)
            f2(1) = f2(1) + fmp(k,i)*fphia(deriv2(k),i)
            f3(1) = f3(1) + fmp(k,i)*fphia(deriv3(k),i)
            f1(2) = f1(2) + fmp(k,i)*fphib(deriv1(k),i)
            f2(2) = f2(2) + fmp(k,i)*fphib(deriv2(k),i)
            f3(2) = f3(2) + fmp(k,i)*fphib(deriv3(k),i)
         end do
         do k = 1, 3
            f1(1) = f1(1) + a_frac(k,i)*fphirec(deriv1(k+1),i)
            f2(1) = f2(1) + a_frac(k,i)*fphirec(deriv2(k+1),i)
            f3(1) = f3(1) + a_frac(k,i)*fphirec(deriv3(k+1),i)
            f1(2) = f1(2) + b_frac(k,i)*fphirec(deriv1(k+1),i)
            f2(2) = f2(2) + b_frac(k,i)*fphirec(deriv2(k+1),i)
            f3(2) = f3(2) + b_frac(k,i)*fphirec(deriv3(k+1),i)

            !mat
            g1(1) = g1(1) + d1_frac(k,i)*fphie1(deriv1(k+1),i) 
     $                    + e1_frac(k,i)*fphid1(deriv1(k+1),i)
            g2(1) = g2(1) + d1_frac(k,i)*fphie1(deriv2(k+1),i)
     $                    + e1_frac(k,i)*fphid1(deriv2(k+1),i)
            g3(1) = g3(1) + d1_frac(k,i)*fphie1(deriv3(k+1),i)
     $                    + e1_frac(k,i)*fphid1(deriv3(k+1),i)
            g1(2) = g1(2) + d2_frac(k,i)*fphie2(deriv1(k+1),i) 
     $                    + e2_frac(k,i)*fphid2(deriv1(k+1),i)
            g2(2) = g2(2) + d2_frac(k,i)*fphie2(deriv2(k+1),i)
     $                    + e2_frac(k,i)*fphid2(deriv2(k+1),i)
            g3(2) = g3(2) + d2_frac(k,i)*fphie2(deriv3(k+1),i)
     $                    + e2_frac(k,i)*fphid2(deriv3(k+1),i)
            g1(3) = g1(3) + d3_frac(k,i)*fphie3(deriv1(k+1),i) 
     $                    + e3_frac(k,i)*fphid3(deriv1(k+1),i)
            g2(3) = g2(3) + d3_frac(k,i)*fphie3(deriv2(k+1),i)
     $                    + e3_frac(k,i)*fphid3(deriv2(k+1),i)
            g3(3) = g3(3) + d3_frac(k,i)*fphie3(deriv3(k+1),i)
     $                    + e3_frac(k,i)*fphid3(deriv3(k+1),i)
         end do

         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3

         !mat
         g1 = dble(nfft1) * g1
         g2 = dble(nfft2) * g2
         g3 = dble(nfft3) * g3
         do k = 1,2
           frc(1,i,k)=recip(1,1)*f1(k)+recip(1,2)*f2(k)+recip(1,3)*f3(k)
           frc(2,i,k)=recip(2,1)*f1(k)+recip(2,2)*f2(k)+recip(2,3)*f3(k)
           frc(3,i,k)=recip(3,1)*f1(k)+recip(3,2)*f2(k)+recip(3,3)*f3(k)
         end do
         !mat
         do k = 1,3
           grc(1,i,k)=recip(1,1)*g1(k)+recip(1,2)*g2(k)+recip(1,3)*g3(k)
           grc(2,i,k)=recip(2,1)*g1(k)+recip(2,2)*g2(k)+recip(2,3)*g3(k)
           grc(3,i,k)=recip(3,1)*g1(k)+recip(3,2)*g2(k)+recip(3,3)*g3(k)
         end do
      end do

      frcdm = 0d0
      frcdt = 0d0
      do i = 1,npolerecloc
         dt = 0d0
         dm = 0d0
         do k = 1,3

            dm(k,1) = fphia(derivi(k,1), i)
            dm(k,2) = fphib(derivi(k,1), i)
            do l = 1,3
               dt(k,l,1) = cdervec(derivi(k,l+1), i)
               dt(k,l,2) = cdervecb(derivi(k,l+1), i)
            end do
         end do

         dm(1,:) = dble(nfft1)*dm(1,:)
         dm(2,:) = dble(nfft2)*dm(2,:)
         dm(3,:) = dble(nfft3)*dm(3,:)

         do k =1,3
            do j = 1,3
               do l = 1,2
                  frcdm(k,i,l) = frcdm(k,i,l) + recip(k,j)*dm(j,l)
                  frcdt(k,j,i,l) = dt(k,j,l)
               end do
            end do
         end do
      end do

      do ii = 1, npolerecloc
        do k = 1,3 
           !mat
           adTb(1,ii) = adTb(1,ii) + grc(1,ii,k)
           adTb(2,ii) = adTb(2,ii) + grc(2,ii,k)
           adTb(3,ii) = adTb(3,ii) + grc(3,ii,k)
           do l = 1,2
              adErec(k,ii) = aderec(k,ii) - frc(k,ii,l)
              admErec(k,ii) = admerec(k,ii) - frcdm(k,ii,l)
           end do
           do j = 1,3
              do l = 1,2
                 adtErec(k,j,ii) = adterec(k,j,ii) - frcdt(k,j,ii,l)
              end do
           end do
        end do
      end do
c
      deallocate(cmp, fmp, fdervec, cdervec)
      deallocate(a_frac, b_frac, frc, grc)
      deallocate(frcdt, frcdm)

      return
      end

!===================================================
!     sub scalderfieldzmatrec8
!===================================================
! No scaling in E_recip.
! Given :
! - a pair of  vectors a and b, as well ad their 
! potential and subsequent derivatives fphiA, fphiB, 
! returns
! < a, dE_recip/dr> + < b, dE_recip/dr> = aderec
! < a, dE_recip/dmu > + < b, dE_recip/dmu > = admerec
! < a, dE_recip/dtheta > + < a, dE_recip/dtheta > = adterec
! 
! - four pairs Di, Ei, along with the potential and 
! successive derivatives arising fphiAi and fphiBi, 
! computes the reciprocal part of < A, T' B > 
! A and B should be given in cartesian coordinates
!
! Opti : uses no fft
! /!\  Finishes with an extra -1 factor !
!
! 'bis' : only (3,n) vectors in; no nrhs
      subroutine scalderfieldzmatrec8(a,abis,fphia,b,bbis,fphib,
     $            d1,d1bis,fphid1,e1,e1bis,fphie1,d2,d2bis,fphid2, 
     $            e2,e2bis,fphie2,d3,d3bis,fphid3,e3,e3bis,fphie3, 
     $            d4,d4bis,fphid4,e4,e4bis,fphie4,
     $            aderec, admerec, adterec, adtb)
      use atmlst
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
      use potent
      use mpi
      implicit none



      real*8, intent(in), dimension(3,npolebloc) :: a, b, d1, d2, d3,
     $                                                d4,e1, e2, e3, e4
      real*8, intent(in), dimension(3,npolerecloc) :: abis,bbis,d1bis,
     $        d2bis,d3bis,d4bis,e1bis, e2bis, e3bis, e4bis
      real*8, intent(in), dimension(20,npolerecloc) :: fphia, fphib,
     $                  fphid1, fphid2, fphid3, fphid4,fphie1, fphie2, 
     $                  fphie3, fphie4
      real*8, intent(out), dimension(3,npolerecloc) :: aderec, adtb,
     $                                                      admerec
      real*8, intent(out), dimension(3,3,npolerecloc) :: adterec

      integer :: i, iipole, iglob, rankloc, j, k,
     $           l, ii, iloc
      integer, dimension(10) :: deriv1, deriv2, deriv3
      integer, dimension(3,10) :: derivi


      real*8, dimension(2) :: f1, f2, f3
      real*8, dimension(4) :: g1, g2, g3
      real*8, dimension(3,2) :: dm
      real*8, dimension(3,3,2) :: dt
      real*8, allocatable, dimension(:,:) :: cmp, fmp, a_frac, b_frac,
     $          fdervec, cdervec,
     $          fdervecdp, cdervecb, d1_frac, d2_frac, d3_frac,d4_frac,
     $          e1_frac, e2_frac, e3_frac, e4_frac, adtb_temp
      real*8, allocatable, dimension(:,:,:) :: frc, 
     $          frcdm, grc
      real*8, allocatable, dimension(:,:,:,:) :: frcdt
      

      rankloc = rank

      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
      derivi(1,:) = deriv1
      derivi(2,:) = deriv2
      derivi(3,:) = deriv3

      aderec = 0d0
      admerec = 0d0
      adterec = 0d0

      allocate(cmp(10,npolerecloc), fmp(10,npolerecloc))
      allocate(fdervec(20,npolerecloc), cdervec(20,npolerecloc),
     $ fdervecdp(20,npolerecloc), 
     $          cdervecb(20,npolerecloc))
      allocate(frc(3,npolerecloc,2), grc(3,npolerecloc,4))
      allocate(a_frac(3,npolerecloc), b_frac(3,npolerecloc),
     $         adtb_temp(3,npolerecloc),
     $         d1_frac(3,npolerecloc), e1_frac(3,npolerecloc),
     $         d2_frac(3,npolerecloc), e2_frac(3,npolerecloc),
     $         d3_frac(3,npolerecloc), e3_frac(3,npolerecloc), 
     $         d4_frac(3,npolerecloc), e4_frac(3,npolerecloc))
      allocate(frcdm(3,npolerecloc,2))
      allocate(frcdt(3,3,npolerecloc,2))

      cmp = 0d0
      fmp = 0d0
      frc = 0d0
      grc = 0d0
      adtb = 0d0
      
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         iloc = poleloc(iipole)
         cmp(1,i) =  rpole(1,iipole)
         cmp(2,i) =  rpole(2,iipole)
         cmp(3,i) =  rpole(3,iipole)
         cmp(4,i) =  rpole(4,iipole)
         cmp(5,i) =  rpole(5,iipole)
         cmp(6,i) =  rpole(9,iipole)
         cmp(7,i) =  rpole(13,iipole)
         cmp(8,i) =  rpole(6,iipole)  *2d0
         cmp(9,i) =  rpole(7,iipole)  *2d0
         cmp(10,i) = rpole(10,iipole)*2d0
         call cmp_to_fmp_site(cmp(1,i), fmp(1,i))
         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(abis(:,i), a_frac(:,i))
           call cart_to_frac_vec(bbis(:,i), b_frac(:,i))
         else
           call cart_to_frac_vec(a(:,iloc), a_frac(:,i))
           call cart_to_frac_vec(b(:,iloc), b_frac(:,i))
         end if
         call fphi_to_cphi_site(fphia(:,i), cdervec(:,i))
         call fphi_to_cphi_site(fphib(:,i), cdervecb(:,i))
         !mat
         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(d1bis(:,i), d1_frac(:,i))
           call cart_to_frac_vec(d2bis(:,i), d2_frac(:,i))
           call cart_to_frac_vec(d3bis(:,i), d3_frac(:,i))
           call cart_to_frac_vec(d4bis(:,i), d4_frac(:,i))
           call cart_to_frac_vec(e1bis(:,i), e1_frac(:,i))
           call cart_to_frac_vec(e2bis(:,i), e2_frac(:,i))
           call cart_to_frac_vec(e3bis(:,i), e3_frac(:,i))
           call cart_to_frac_vec(e4bis(:,i), e4_frac(:,i))
         else
           call cart_to_frac_vec(d1(:,iloc), d1_frac(:,i))
           call cart_to_frac_vec(d2(:,iloc), d2_frac(:,i))
           call cart_to_frac_vec(d3(:,iloc), d3_frac(:,i))
           call cart_to_frac_vec(d4(:,iloc), d4_frac(:,i))
           call cart_to_frac_vec(e1(:,iloc), e1_frac(:,i))
           call cart_to_frac_vec(e2(:,iloc), e2_frac(:,i))
           call cart_to_frac_vec(e3(:,iloc), e3_frac(:,i))
           call cart_to_frac_vec(e4(:,iloc), e4_frac(:,i))
         end if
      end do


      do i = 1, npolerecloc
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         !mat
         g1 = 0d0
         g2 = 0d0
         g3 = 0d0
         do k = 1, 10
            f1(1) = f1(1) + fmp(k,i)*fphia(deriv1(k),i)
            f2(1) = f2(1) + fmp(k,i)*fphia(deriv2(k),i)
            f3(1) = f3(1) + fmp(k,i)*fphia(deriv3(k),i)
            f1(2) = f1(2) + fmp(k,i)*fphib(deriv1(k),i)
            f2(2) = f2(2) + fmp(k,i)*fphib(deriv2(k),i)
            f3(2) = f3(2) + fmp(k,i)*fphib(deriv3(k),i)
         end do
         do k = 1, 3
            f1(1) = f1(1) + a_frac(k,i)*fphirec(deriv1(k+1),i)
            f2(1) = f2(1) + a_frac(k,i)*fphirec(deriv2(k+1),i)
            f3(1) = f3(1) + a_frac(k,i)*fphirec(deriv3(k+1),i)
            f1(2) = f1(2) + b_frac(k,i)*fphirec(deriv1(k+1),i)
            f2(2) = f2(2) + b_frac(k,i)*fphirec(deriv2(k+1),i)
            f3(2) = f3(2) + b_frac(k,i)*fphirec(deriv3(k+1),i)

            !mat
            g1(1) = g1(1) + d1_frac(k,i)*fphie1(deriv1(k+1),i) 
     $                    + e1_frac(k,i)*fphid1(deriv1(k+1),i)
            g2(1) = g2(1) + d1_frac(k,i)*fphie1(deriv2(k+1),i)
     $                    + e1_frac(k,i)*fphid1(deriv2(k+1),i)
            g3(1) = g3(1) + d1_frac(k,i)*fphie1(deriv3(k+1),i)
     $                    + e1_frac(k,i)*fphid1(deriv3(k+1),i)
            g1(2) = g1(2) + d2_frac(k,i)*fphie2(deriv1(k+1),i) 
     $                    + e2_frac(k,i)*fphid2(deriv1(k+1),i)
            g2(2) = g2(2) + d2_frac(k,i)*fphie2(deriv2(k+1),i)
     $                    + e2_frac(k,i)*fphid2(deriv2(k+1),i)
            g3(2) = g3(2) + d2_frac(k,i)*fphie2(deriv3(k+1),i)
     $                    + e2_frac(k,i)*fphid2(deriv3(k+1),i)
            g1(3) = g1(3) + d3_frac(k,i)*fphie3(deriv1(k+1),i) 
     $                    + e3_frac(k,i)*fphid3(deriv1(k+1),i)
            g2(3) = g2(3) + d3_frac(k,i)*fphie3(deriv2(k+1),i)
     $                    + e3_frac(k,i)*fphid3(deriv2(k+1),i)
            g3(3) = g3(3) + d3_frac(k,i)*fphie3(deriv3(k+1),i)
     $                    + e3_frac(k,i)*fphid3(deriv3(k+1),i)
            g1(4) = g1(4) + d4_frac(k,i)*fphie4(deriv1(k+1),i) 
     $                    + e4_frac(k,i)*fphid4(deriv1(k+1),i)
            g2(4) = g2(4) + d4_frac(k,i)*fphie4(deriv2(k+1),i)
     $                    + e4_frac(k,i)*fphid4(deriv2(k+1),i)
            g3(4) = g3(4) + d4_frac(k,i)*fphie4(deriv3(k+1),i)
     $                    + e4_frac(k,i)*fphid4(deriv3(k+1),i)
         end do

         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3

         !mat
         g1 = dble(nfft1) * g1
         g2 = dble(nfft2) * g2
         g3 = dble(nfft3) * g3
         do k = 1,2
           frc(1,i,k)=recip(1,1)*f1(k)+recip(1,2)*f2(k)+recip(1,3)*f3(k)
           frc(2,i,k)=recip(2,1)*f1(k)+recip(2,2)*f2(k)+recip(2,3)*f3(k)
           frc(3,i,k)=recip(3,1)*f1(k)+recip(3,2)*f2(k)+recip(3,3)*f3(k)
         end do
         !mat
         do k = 1,4
           grc(1,i,k)=recip(1,1)*g1(k)+recip(1,2)*g2(k)+recip(1,3)*g3(k)
           grc(2,i,k)=recip(2,1)*g1(k)+recip(2,2)*g2(k)+recip(2,3)*g3(k)
           grc(3,i,k)=recip(3,1)*g1(k)+recip(3,2)*g2(k)+recip(3,3)*g3(k)
         end do
      end do

      frcdm = 0d0
      frcdt = 0d0
      do i = 1,npolerecloc
         dt = 0d0
         dm = 0d0
         do k = 1,3

            dm(k,1) = fphia(derivi(k,1), i)
            dm(k,2) = fphib(derivi(k,1), i)
            do l = 1,3
               dt(k,l,1) = cdervec(derivi(k,l+1), i)
               dt(k,l,2) = cdervecb(derivi(k,l+1), i)
            end do
         end do

         dm(1,:) = dble(nfft1)*dm(1,:)
         dm(2,:) = dble(nfft2)*dm(2,:)
         dm(3,:) = dble(nfft3)*dm(3,:)

         do k =1,3
            do j = 1,3
               do l = 1,2
                  frcdm(k,i,l) = frcdm(k,i,l) + recip(k,j)*dm(j,l)
                  frcdt(k,j,i,l) = dt(k,j,l)
               end do
            end do
         end do
      end do

      do ii = 1, npolerecloc
        do k = 1,3 
           do l = 1,2
              adErec(k,ii) = aderec(k,ii) - frc(k,ii,l)
              admErec(k,ii) = admerec(k,ii) - frcdm(k,ii,l)
           end do
           do j = 1,3
              do l = 1,2
                 adtErec(k,j,ii) = adterec(k,j,ii) - frcdt(k,j,ii,l)
              end do
           end do
        end do
        do k = 1, 4
          !mat
          adTb(1,ii) = adTb(1,ii) + grc(1,ii,k)
          adTb(2,ii) = adTb(2,ii) + grc(2,ii,k)
          adTb(3,ii) = adTb(3,ii) + grc(3,ii,k)
        end do
      end do
c

      deallocate(cmp, fmp, fdervec, cdervec)
      deallocate(a_frac, b_frac, frc, grc)
      deallocate(frcdt, frcdm)

      return
      end

!===================================================
!     sub scalderfieldzmatrec3
!===================================================
! No scaling in E_recip.
! Given :
! - a pair of  vectors a and b, as well ad their 
! potential and subsequent derivatives fphiA, fphiB, 
! returns
! < a, dE_recip/dr> + < b, dE_recip/dr> = aderec
! < a, dE_recip/dmu > + < b, dE_recip/dmu > = admerec
! < a, dE_recip/dtheta > + < a, dE_recip/dtheta > = adterec
! 
! - two pairs Di, Ei, along with the potential and 
! successive derivatives arising fphiAi and fphiBi, 
! computes the reciprocal part of < A, T' B > 
! A and B should be given in cartesian coordinates
!
! Opti : uses no fft
! /!\  Finishes with an extra -1 factor !
!
      subroutine scalderfieldzmatrec3(a,abis,fphia,b,bbis,
     $            fphib,d1,d1bis,fphid1,e1,e1bis,fphie1,d2,d2bis,
     $            fphid2,e2,e2bis, fphie2, 
     $            aderec, admerec, adterec, adtb)
      use atmlst
      use bound
      use boxes
      use chgpot
      use domdec
      use sizes
      use energi
      use ewald
      use fft
      use math
      use mpole
      use pme
      use potent
      use mpi
      implicit none



      real*8, intent(in), dimension(3,npolebloc) :: a, b, d1, d2,
     $                                                e1, e2
      real*8, intent(in), dimension(3,npolerecloc) :: abis,bbis,
     $                 d1bis,d2bis,e1bis,e2bis
      real*8, intent(in), dimension(20,npolerecloc) :: fphia, fphib,
     $                  fphid1, fphid2, fphie1, fphie2
      real*8, intent(out), dimension(3,npolerecloc) :: aderec, adtb,
     $                                                      admerec
      real*8, intent(out), dimension(3,3,npolerecloc) :: adterec

      integer :: i, iipole, iglob, rankloc, j, k,
     $           l, ii, iloc
      integer, dimension(10) :: deriv1, deriv2, deriv3
      integer, dimension(3,10) :: derivi


      real*8, dimension(2) :: f1, f2, f3
      real*8, dimension(2) :: g1, g2, g3
      real*8, dimension(3,2) :: dm
      real*8, dimension(3,3,2) :: dt
      real*8, allocatable, dimension(:,:) :: cmp, fmp, a_frac, b_frac,
     $          fdervec, cdervec,
     $          fdervecdp, cdervecb, d1_frac, d2_frac, 
     $          e1_frac, e2_frac, adtb_temp
      real*8, allocatable, dimension(:,:,:) :: frc, 
     $          frcdm, grc
      real*8, allocatable, dimension(:,:,:,:) :: frcdt


      rankloc = rank

      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /


      derivi(1,:) = deriv1
      derivi(2,:) = deriv2
      derivi(3,:) = deriv3

      aderec = 0d0
      admerec = 0d0
      adterec = 0d0

      allocate(cmp(10,npolerecloc), fmp(10,npolerecloc))
      allocate(fdervec(20,npolerecloc), cdervec(20,npolerecloc),
     $fdervecdp(20,npolerecloc),cdervecb(20,npolerecloc))
      allocate(frc(3,npolerecloc,2), grc(3,npolerecloc,2))
      allocate(a_frac(3,npolerecloc), b_frac(3,npolerecloc),
     $         adtb_temp(3,npolerecloc),
     $         d1_frac(3,npolerecloc), e1_frac(3,npolerecloc),
     $         d2_frac(3,npolerecloc), e2_frac(3,npolerecloc))
      allocate(frcdm(3,npolerecloc,2))
      allocate(frcdt(3,3,npolerecloc,2))

      cmp = 0d0
      fmp = 0d0
      frc = 0d0
      grc = 0d0
      adtb = 0d0

      
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         iloc = poleloc(iipole)
         cmp(1,i) =  rpole(1,iipole)
         cmp(2,i) =  rpole(2,iipole)
         cmp(3,i) =  rpole(3,iipole)
         cmp(4,i) =  rpole(4,iipole)
         cmp(5,i) =  rpole(5,iipole)
         cmp(6,i) =  rpole(9,iipole)
         cmp(7,i) =  rpole(13,iipole)
         cmp(8,i) =  rpole(6,iipole)  *2d0
         cmp(9,i) =  rpole(7,iipole)  *2d0
         cmp(10,i) = rpole(10,iipole)*2d0
         call cmp_to_fmp_site(cmp(1,i), fmp(1,i))
         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(abis(:,i), a_frac(:,i))
           call cart_to_frac_vec(bbis(:,i), b_frac(:,i))
         else
           call cart_to_frac_vec(a(:,iloc), a_frac(:,i))
           call cart_to_frac_vec(b(:,iloc), b_frac(:,i))
         end if
         call fphi_to_cphi_site(fphia(:,i), cdervec(:,i))
         call fphi_to_cphi_site(fphib(:,i), cdervecb(:,i))
         !mat
         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(d1bis(:,i), d1_frac(:,i))
           call cart_to_frac_vec(d2bis(:,i), d2_frac(:,i))
           call cart_to_frac_vec(e1bis(:,i), e1_frac(:,i))
           call cart_to_frac_vec(e2bis(:,i), e2_frac(:,i))
         else
           call cart_to_frac_vec(d1(:,iloc), d1_frac(:,i))
           call cart_to_frac_vec(d2(:,iloc), d2_frac(:,i))
           call cart_to_frac_vec(e1(:,iloc), e1_frac(:,i))
           call cart_to_frac_vec(e2(:,iloc), e2_frac(:,i))
         end if
      end do


      do i = 1, npolerecloc
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         !mat
         g1 = 0d0
         g2 = 0d0
         g3 = 0d0
         do k = 1, 10
            f1(1) = f1(1) + fmp(k,i)*fphia(deriv1(k),i)
            f2(1) = f2(1) + fmp(k,i)*fphia(deriv2(k),i)
            f3(1) = f3(1) + fmp(k,i)*fphia(deriv3(k),i)
            f1(2) = f1(2) + fmp(k,i)*fphib(deriv1(k),i)
            f2(2) = f2(2) + fmp(k,i)*fphib(deriv2(k),i)
            f3(2) = f3(2) + fmp(k,i)*fphib(deriv3(k),i)
         end do
         do k = 1, 3
            f1(1) = f1(1) + a_frac(k,i)*fphirec(deriv1(k+1),i)
            f2(1) = f2(1) + a_frac(k,i)*fphirec(deriv2(k+1),i)
            f3(1) = f3(1) + a_frac(k,i)*fphirec(deriv3(k+1),i)
            f1(2) = f1(2) + b_frac(k,i)*fphirec(deriv1(k+1),i)
            f2(2) = f2(2) + b_frac(k,i)*fphirec(deriv2(k+1),i)
            f3(2) = f3(2) + b_frac(k,i)*fphirec(deriv3(k+1),i)

            !mat
            g1(1) = g1(1) + d1_frac(k,i)*fphie1(deriv1(k+1),i) 
     $                    + e1_frac(k,i)*fphid1(deriv1(k+1),i)
            g2(1) = g2(1) + d1_frac(k,i)*fphie1(deriv2(k+1),i)
     $                    + e1_frac(k,i)*fphid1(deriv2(k+1),i)
            g3(1) = g3(1) + d1_frac(k,i)*fphie1(deriv3(k+1),i)
     $                    + e1_frac(k,i)*fphid1(deriv3(k+1),i)
            g1(2) = g1(2) + d2_frac(k,i)*fphie2(deriv1(k+1),i) 
     $                    + e2_frac(k,i)*fphid2(deriv1(k+1),i)
            g2(2) = g2(2) + d2_frac(k,i)*fphie2(deriv2(k+1),i)
     $                    + e2_frac(k,i)*fphid2(deriv2(k+1),i)
            g3(2) = g3(2) + d2_frac(k,i)*fphie2(deriv3(k+1),i)
     $                    + e2_frac(k,i)*fphid2(deriv3(k+1),i)
         end do

         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3

         !mat
         g1 = dble(nfft1) * g1
         g2 = dble(nfft2) * g2
         g3 = dble(nfft3) * g3
         do k = 1,2
           frc(1,i,k)=recip(1,1)*f1(k)+recip(1,2)*f2(k)+recip(1,3)*f3(k)
           frc(2,i,k)=recip(2,1)*f1(k)+recip(2,2)*f2(k)+recip(2,3)*f3(k)
           frc(3,i,k)=recip(3,1)*f1(k)+recip(3,2)*f2(k)+recip(3,3)*f3(k)
         end do
         !mat
         do k = 1,2
           grc(1,i,k)=recip(1,1)*g1(k)+recip(1,2)*g2(k)+recip(1,3)*g3(k)
           grc(2,i,k)=recip(2,1)*g1(k)+recip(2,2)*g2(k)+recip(2,3)*g3(k)
           grc(3,i,k)=recip(3,1)*g1(k)+recip(3,2)*g2(k)+recip(3,3)*g3(k)
         end do
      end do

      frcdm = 0d0
      frcdt = 0d0
      do i = 1,npolerecloc
         dt = 0d0
         dm = 0d0
         do k = 1,3

            dm(k,1) = fphia(derivi(k,1), i)
            dm(k,2) = fphib(derivi(k,1), i)
            do l = 1,3
               dt(k,l,1) = cdervec(derivi(k,l+1), i)
               dt(k,l,2) = cdervecb(derivi(k,l+1), i)
            end do
         end do

         dm(1,:) = dble(nfft1)*dm(1,:)
         dm(2,:) = dble(nfft2)*dm(2,:)
         dm(3,:) = dble(nfft3)*dm(3,:)

         do k =1,3
            do j = 1,3
               do l = 1,2
                  frcdm(k,i,l) = frcdm(k,i,l) + recip(k,j)*dm(j,l)
                  frcdt(k,j,i,l) = dt(k,j,l)
               end do
            end do
         end do
      end do

      do ii = 1, npolerecloc

        do k = 1,3 
           do l = 1,2
              adErec(k,ii) = aderec(k,ii) - frc(k,ii,l)
              admErec(k,ii) = admerec(k,ii) - frcdm(k,ii,l)
           end do
           do j = 1,3
              do l = 1,2
                 adtErec(k,j,ii) = adterec(k,j,ii) - frcdt(k,j,ii,l)
              end do
           end do
        end do
        do k = 1,2
           !mat
           adTb(1,ii) = adTb(1,ii) + grc(1,ii,k)
           adTb(2,ii) = adTb(2,ii) + grc(2,ii,k)
           adTb(3,ii) = adTb(3,ii) + grc(3,ii,k)
         end do
      end do

      deallocate(cmp, fmp, fdervec, cdervec)
      deallocate(a_frac, b_frac, frc, grc)
      deallocate(frcdt, frcdm)

      return
      end

!===================================================
!     sub scalderfieldzmatrec1
!===================================================
! No scaling in E_recip.
! Given :
! - a pair of  vectors a and b, as well ad their 
! potential and subsequent derivatives fphiA, fphiB, 
! returns
! < a, dE_recip/dr> + < b, dE_recip/dr> = aderec
! < a, dE_recip/dmu > + < b, dE_recip/dmu > = admerec
! < a, dE_recip/dtheta > + < b, dE_recip/dtheta > = adterec
! 
! - ONE pair D, E, along with the potential and 
! successive derivatives arising fphiD and fphiE, 
! computes the reciprocal part of < D, T' E > 
! D and E should be given in cartesian coordinates
!
! Opti : uses no fft
! /!\  Finishes with an extra -1 factor !
!
      subroutine scalderfieldzmatrec1(a,abis,fphia,b,bbis,fphib,
     $            d1,d1bis,fphid1,e1,e1bis,fphie1, 
     $            aderec, admerec, adterec, adtb)
      use atmlst
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
      use potent
      use virial
      use mpi
      implicit none



      real*8, intent(in), dimension(3,npoleloc) :: a, b, d1,
     $                                                e1 
      real*8, intent(in), dimension(3,npolerecloc) :: abis, bbis, d1bis,
     $                                                e1bis
      real*8, intent(in), dimension(20,npolerecloc) :: fphia, fphib,
     $                  fphid1, fphie1
      real*8, intent(out), dimension(3,npolerecloc) :: aderec, adtb,
     $                                                      admerec
      real*8, intent(out), dimension(3,3,npolerecloc) :: adterec

      integer :: i, iipole, iglob, rankloc, j, k,
     $           l, ii, iloc
      integer, dimension(10) :: deriv1, deriv2, deriv3
      integer, dimension(3,10) :: derivi


      real*8, dimension(2) :: f1, f2, f3
      real*8 :: g1, g2, g3
      real*8, dimension(3,2) :: dm
      real*8, dimension(3,3,2) :: dt
      real*8, allocatable, dimension(:,:) :: cmp, fmp, a_frac, b_frac,
     $          fdervec, cdervec, admerec_temp, aderec_temp,
     $          fdervecdp, cdervecb, d1_frac, 
     $          e1_frac,  adtb_temp, grc
      real*8, allocatable, dimension(:,:,:) :: adtErec_temp, frc, 
     $          frcdm
      real*8, allocatable, dimension(:,:,:,:) :: frcdt
      real*8 :: vxx,vxy,vxz,vyy,vyz,vzz,f
      

      rankloc = rank
      f = electric/dielec

      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
      derivi(1,:) = deriv1
      derivi(2,:) = deriv2
      derivi(3,:) = deriv3

      aderec = 0d0
      admerec = 0d0
      adterec = 0d0

      allocate(cmp(10,npolerecloc), fmp(10,npolerecloc))
      allocate(fdervec(20,npolerecloc), cdervec(20,npolerecloc),
     $fdervecdp(20,npolerecloc), 
     $          cdervecb(20,npolerecloc))
      allocate(frc(3,npolerecloc,2), grc(3,npolerecloc))
      allocate(a_frac(3,npolerecloc), b_frac(3,npolerecloc),
     $         adtb_temp(3,nlocrec2),
     $         d1_frac(3,npolerecloc), e1_frac(3,npolerecloc))
      allocate(frcdm(3,npolerecloc,2))
      allocate(frcdt(3,3,npolerecloc,2))
      allocate(admerec_temp(3,nlocrec2), 
     $         adErec_temp(3,nlocrec2),
     $         adtErec_temp(3,3,nlocrec2))

      cmp = 0d0
      fmp = 0d0
      frc = 0d0
      grc = 0d0
      adtb = 0d0

      vxx = 0d0
      vxy = 0d0
      vxz = 0d0
      vyy = 0d0
      vyz = 0d0
      vzz = 0d0
      
      do i = 1, npolerecloc
         iipole = polerecglob(i)
         iglob = ipole(iipole)
         iloc = poleloc(iipole)
         cmp(1,i) =  rpole(1,iipole)
         cmp(2,i) =  rpole(2,iipole)
         cmp(3,i) =  rpole(3,iipole)
         cmp(4,i) =  rpole(4,iipole)
         cmp(5,i) =  rpole(5,iipole)
         cmp(6,i) =  rpole(9,iipole)
         cmp(7,i) =  rpole(13,iipole)
         cmp(8,i) =  rpole(6,iipole)  *2d0
         cmp(9,i) =  rpole(7,iipole)  *2d0
         cmp(10,i) = rpole(10,iipole)*2d0
         call cmp_to_fmp_site(cmp(1,i), fmp(1,i))
         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(abis(:,i), a_frac(:,i))
           call cart_to_frac_vec(bbis(:,i), b_frac(:,i))
         else
           call cart_to_frac_vec(a(:,iloc), a_frac(:,i))
           call cart_to_frac_vec(b(:,iloc), b_frac(:,i))
         end if
         call fphi_to_cphi_site(fphia(:,i), cdervec(:,i))
         call fphi_to_cphi_site(fphib(:,i), cdervecb(:,i))
         !mat
         if (repart(iglob).ne.rank) then
           call cart_to_frac_vec(d1bis(:,i), d1_frac(:,i))
           call cart_to_frac_vec(e1bis(:,i), e1_frac(:,i))
         else
           call cart_to_frac_vec(d1(:,iloc), d1_frac(:,i))
           call cart_to_frac_vec(e1(:,iloc), e1_frac(:,i))
         end if
c         vxx = vxx - fphirec(2,i)*(a_frac(1,i)+b_frac(1,i))
c         vxy = vxy - 0.5d0*(fphirec(2,i)*(a_frac(2,i)+b_frac(2,i))+
c     $    fphirec(3,i)*(a_frac(1,i)+b_frac(1,i)))
c         vxz = vxz - 0.5d0*(fphirec(2,i)*(a_frac(3,i)+b_frac(3,i))+
c     $    fphirec(4,i)*(a_frac(1,i)+b_frac(1,i)))
c         vyy = vyy - fphirec(3,i)*(a_frac(2,i)+b_frac(2,i))
c         vyz = vyz - 0.5d0*(fphirec(3,i)*(a_frac(3,i)+b_frac(3,i))+
c     $    fphirec(4,i)*(a_frac(2,i)+b_frac(2,i)))
c         vzz = vzz - fphirec(4,i)*(a_frac(3,i)+b_frac(3,i))
c         !mat
c         vxx = vxx - fphid1(2,i)*e1_frac(1,i)
c         vxy = vxy - 0.5d0*(fphid1(2,i)*e1_frac(2,i)+
c     $    fphid1(3,i)*e1_frac(1,i))
c         vxz = vxz - 0.5d0*(fphid1(2,i)*e1_frac(3,i)+
c     $    fphid1(4,i)*e1_frac(1,i))
c         vyy = vyy - fphid1(3,i)*e1_frac(2,i)
c         vyz = vyz - 0.5d0*(fphid1(3,i)*e1_frac(3,i)+
c     $    fphid1(4,i)*e1_frac(2,i))
c         vzz = vzz - fphid1(4,i)*e1_frac(3,i)
      end do
cc
cc    increment virial
cc
c      vir(1,1) = vir(1,1) -0.5*f*vxx
c      vir(2,1) = vir(2,1) -0.5*f*vxy
c      vir(3,1) = vir(3,1) -0.5*f*vxz
c      vir(1,2) = vir(1,2) -0.5*f*vxy
c      vir(2,2) = vir(2,2) -0.5*f*vyy
c      vir(3,2) = vir(3,2) -0.5*f*vyz
c      vir(1,3) = vir(1,3) -0.5*f*vxz
c      vir(2,3) = vir(2,3) -0.5*f*vyz
c      vir(3,3) = vir(3,3) -0.5*f*vzz


      do i = 1, npolerecloc
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         !mat
         g1 = 0d0
         g2 = 0d0
         g3 = 0d0
         do k = 1, 10
            f1(1) = f1(1) + fmp(k,i)*fphia(deriv1(k),i)
            f2(1) = f2(1) + fmp(k,i)*fphia(deriv2(k),i)
            f3(1) = f3(1) + fmp(k,i)*fphia(deriv3(k),i)
            f1(2) = f1(2) + fmp(k,i)*fphib(deriv1(k),i)
            f2(2) = f2(2) + fmp(k,i)*fphib(deriv2(k),i)
            f3(2) = f3(2) + fmp(k,i)*fphib(deriv3(k),i)
         end do
         do k = 1, 3
            f1(1) = f1(1) + a_frac(k,i)*fphirec(deriv1(k+1),i)
            f2(1) = f2(1) + a_frac(k,i)*fphirec(deriv2(k+1),i)
            f3(1) = f3(1) + a_frac(k,i)*fphirec(deriv3(k+1),i)
            f1(2) = f1(2) + b_frac(k,i)*fphirec(deriv1(k+1),i)
            f2(2) = f2(2) + b_frac(k,i)*fphirec(deriv2(k+1),i)
            f3(2) = f3(2) + b_frac(k,i)*fphirec(deriv3(k+1),i)

            !mat
            g1 = g1 + d1_frac(k,i)*fphie1(deriv1(k+1),i) 
     $              + e1_frac(k,i)*fphid1(deriv1(k+1),i)
            g2 = g2 + d1_frac(k,i)*fphie1(deriv2(k+1),i)
     $              + e1_frac(k,i)*fphid1(deriv2(k+1),i)
            g3 = g3 + d1_frac(k,i)*fphie1(deriv3(k+1),i)
     $              + e1_frac(k,i)*fphid1(deriv3(k+1),i)
         end do

         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3

         !mat
         g1 = dble(nfft1) * g1
         g2 = dble(nfft2) * g2
         g3 = dble(nfft3) * g3
         do k = 1,2
           frc(1,i,k)=recip(1,1)*f1(k)+recip(1,2)*f2(k)+recip(1,3)*f3(k)
           frc(2,i,k)=recip(2,1)*f1(k)+recip(2,2)*f2(k)+recip(2,3)*f3(k)
           frc(3,i,k)=recip(3,1)*f1(k)+recip(3,2)*f2(k)+recip(3,3)*f3(k)
         end do
         !mat
         grc(1,i)=recip(1,1)*g1+recip(1,2)*g2+recip(1,3)*g3
         grc(2,i)=recip(2,1)*g1+recip(2,2)*g2+recip(2,3)*g3
         grc(3,i)=recip(3,1)*g1+recip(3,2)*g2+recip(3,3)*g3
      end do

      frcdm = 0d0
      frcdt = 0d0
      do i = 1,npolerecloc
         dt = 0d0
         dm = 0d0
         do k = 1,3

            dm(k,1) = fphia(derivi(k,1), i)
            dm(k,2) = fphib(derivi(k,1), i)
            do l = 1,3
               dt(k,l,1) = cdervec(derivi(k,l+1), i)
               dt(k,l,2) = cdervecb(derivi(k,l+1), i)
            end do
         end do

         dm(1,:) = dble(nfft1)*dm(1,:)
         dm(2,:) = dble(nfft2)*dm(2,:)
         dm(3,:) = dble(nfft3)*dm(3,:)

         do k =1,3
            do j = 1,3
               do l = 1,2
                  frcdm(k,i,l) = frcdm(k,i,l) + recip(k,j)*dm(j,l)
                  frcdt(k,j,i,l) = dt(k,j,l)
               end do
            end do
         end do
      end do

      do ii = 1, npolerecloc

        !mat
        adTb(1,ii) = adTb(1,ii) + grc(1,ii)
        adTb(2,ii) = adTb(2,ii) + grc(2,ii)
        adTb(3,ii) = adTb(3,ii) + grc(3,ii)
        do k = 1,3 
           do l = 1,2
              adErec(k,ii) = aderec(k,ii) - frc(k,ii,l)
              admErec(k,ii) = admerec(k,ii) - frcdm(k,ii,l)
           end do
           do j = 1,3
              do l = 1,2
                 adtErec(k,j,ii) = adterec(k,j,ii) - frcdt(k,j,ii,l)
              end do
           end do
        end do
      end do

      deallocate(cmp, fmp, fdervec, cdervec)
      deallocate(a_frac, b_frac, frc, grc)
      deallocate(frcdt, frcdm)
      deallocate(adErec_temp,admErec_temp, adtErec_temp)
      deallocate(adtb_temp)

      return
      end

      subroutine scalderfieldzmat1( a1, a2, d1, e1, 
     $                             ade, adme, adte, adtb)

      use atmlst
      use atoms
      use bound
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
      use potent
      use shunt
      implicit none


      real*8, dimension(3,npolebloc), intent(in) :: a1, a2
      real*8, dimension(3,npolebloc), intent(in) :: e1, d1
      real*8, dimension(3,npolebloc), intent(out) :: adtb
      real*8, dimension(3,2,npolebloc), intent(out) :: adme, ade
      real*8, dimension(3,3,2,npolebloc), intent(out) :: adte
      
      real*8, allocatable, dimension(:,:,:) ::adE_cor,admE_cor,admE_self
      real*8, allocatable, dimension(:,:,:,:) ::adtE_cor




      integer :: ii, iipole, iglob, i, jjpole, jglob,
     $           jjj, jbis, k, j

      real*8 :: pti, pdi, rrij3, rrij5, rrij7, rrij9,
     $          urr5, urr7, 
     $scale3, scale5, 
     $          scale7, scale9, dsc3, dsc5, dsc7, dsc9, psc3, psc5,
     $          psc7, psc9, usc5, usc7, scalmuir, scalmujr, 
     $quadterm, quadterm_bis, 
     $          ci, cj, bfac, invrij2, alsq2n, cutoff2, exp2a, ralpha,
     $          d, rij2, damp, expdamp, pgamma, alsq2, term

      real*8                :: bn(0:4) 
      real*8, dimension(2)  :: srr3, srr5, srr7, srr9
      real*8, dimension(3)  :: rij, thetajr, thetair, di, dj
      real*8, dimension(3,3):: qi, qj

      real*8, allocatable, dimension(:) :: dscale, pscale, uscale
      real*8, allocatable, dimension(:,:) :: adtbcor
      

      real*8:: a1xi,a1yi,a1zi,a1xj,a1yj,a1zj
      real*8:: a2xi,a2yi,a2zi,a2xj,a2yj,a2zj
      real*8,dimension(2):: term1,term2,term3,term4,term5,term6,term7
      real*8,dimension(2):: dexdxi,dexdxj,deydyi,deydyj,dezdzi,dezdzj
      real*8,dimension(2):: dexdyi,dexdyj,dexdzi,dexdzj,deydzi,deydzj
      real*8,dimension(2):: demxdx,demxdy,demxdz,demydy,demydz,demzdz
      real*8,dimension(2):: detxdx,detxdy,detxdz,detydy,detydz,detzdz
      real*8 :: depx,depy,depz,dempx,dempy,dempz,detx,dety,detz
      real*8 :: scale1i,scale1j
      real*8::temp1jxx5,temp1jxy5,temp1jxz5,temp1jyy5,temp1jyz5,
     $     temp1jzz5
      real*8::temp1jxx7,temp1jxy7,temp1jxz7,temp1jyy7,temp1jyz7,
     $     temp1jzz7
      real*8::temp1jxx,temp1jxy,temp1jxz,temp1jyy,temp1jyz,temp1jzz
      real*8::temp1ixx5,temp1ixy5,temp1ixz5,temp1iyy5,temp1iyz5,
     $     temp1izz5
      real*8::temp1ixx7,temp1ixy7,temp1ixz7,temp1iyy7,temp1iyz7,
     $     temp1izz7
      real*8::temp1ixx,temp1ixy,temp1ixz,temp1iyy,temp1iyz,temp1izz

      character(10) :: mode

 1000 format(' Warning, system moved too much since last neighbor list
     $   update, try lowering nlupdate')

      

      ade = 0d0
      adme = 0d0
      adte = 0d0
      adtb = 0d0

      mode = 'EWALD'
      if (use_polarshortreal) mode = 'SHORTEWALD'

      call switch (mode)
      cutoff2 = cut2

      allocate(dscale(n), pscale(n), uscale(n))
      allocate(adtbcor(3,npolebloc))
      allocate(adE_cor(3,2,npolebloc),admE_cor(3,2,npolebloc))
      allocate(admE_self(3,2,npolebloc),adtE_cor(3,3,2,npolebloc))
      adme_cor = 0d0
      adme_self = 0d0
      adte_cor = 0d0
      ade_cor = 0d0
      dscale = 1d0
      pscale = 1d0
      uscale = 1.0d0
      adtbcor = 0d0

      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = poleloc(iglob)
         if ((i.eq.0).or.(i.gt.nbloc)) then
            write(iout,1000)
            cycle
         end if
         pdi = pdamp(iipole)
         pti = thole(iipole)

         ci      = rpole(1,iipole)
         di(1)   = rpole(2,iipole)
         di(2)   = rpole(3,iipole)
         di(3)   = rpole(4,iipole)
         qi(1,1) = rpole(5, iipole)
         qi(2,1) = rpole(6, iipole)
         qi(3,1) = rpole(7, iipole)
         qi(1,2) = rpole(6, iipole)
         qi(2,2) = rpole(9, iipole)
         qi(3,2) = rpole(10,iipole)
         qi(1,3) = rpole(7, iipole)
         qi(2,3) = rpole(10,iipole)
         qi(3,3) = rpole(13,iipole)

         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = p2scale
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = p3scale
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = p4scale
            do k = 1, np11(iglob)
               if (i14(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i14(j,iglob)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = p5scale
         end do
         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = d1scale
            uscale(ip11(j,iglob)) = u1scale
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = d2scale
            uscale(ip12(j,iglob)) = u2scale
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = d3scale
            uscale(ip13(j,iglob)) = u3scale
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = d4scale
            uscale(ip14(j,iglob)) = u4scale
         end do

         do jjj = 1, nelst(ii)
            jjpole = elst(jjj,ii)
            jglob = ipole(jjpole)
            jbis = poleloc(jjpole)
            if (jbis.eq.0) then
               write(iout,1000)
               cycle
            end if

            !Distances
            rij(1) = x(jglob) - x(iglob)
            rij(2) = y(jglob) - y(iglob)
            rij(3) = z(jglob) - z(iglob)
            call image(rij(1), rij(2), rij(3))
            rij2 = dot_product(rij,rij)
            if (rij2 .gt. cutoff2) cycle
            d  = sqrt(rij2)
            invrij2 = 1d0/rij2
            rrij3 = invrij2/d
            rrij5 = 3d0*rrij3*invrij2
            rrij7 = 5d0*rrij5*invrij2
            rrij9 = 7d0*rrij7*invrij2

            !Multipoles
            cj      = rpole(1,jjpole)
            dj(1)   = rpole(2,jjpole)
            dj(2)   = rpole(3,jjpole)
            dj(3)   = rpole(4,jjpole)
            qj(1,1) = rpole(5, jjpole)
            qj(2,1) = rpole(6, jjpole)
            qj(3,1) = rpole(7, jjpole)
            qj(1,2) = rpole(6, jjpole)
            qj(2,2) = rpole(9, jjpole)
            qj(3,2) = rpole(10,jjpole)
            qj(1,3) = rpole(7, jjpole)
            qj(2,3) = rpole(10,jjpole)
            qj(3,3) = rpole(13,jjpole)

            scalmujr = dot_product(dj,rij)!sprod(3,dj,rij)
            scalmuir = dot_product(di,rij)!sprod(3,di,rij)
            thetajr(1) = dot_product(qj(1,:),rij)
            thetajr(2) = dot_product(qj(2,:),rij)
            thetajr(3) = dot_product(qj(3,:),rij)
            thetair(1) = dot_product(qi(1,:),rij)
            thetair(2) = dot_product(qi(2,:),rij)
            thetair(3) = dot_product(qi(3,:),rij)
            quadterm     = dot_product(thetajr,rij)!sprod(3,thetajr,rij)
            quadterm_bis = dot_product(thetair,rij)!sprod(3,thetair,rij)

            scale1j = dot_product(rij,e1(:,jbis))
            scale1i = dot_product(rij,e1(:,i))

            a1xi = a1(1,i)
            a1yi = a1(2,i)
            a1zi = a1(3,i)
            a1xj = a1(1,jbis)
            a1yj = a1(2,jbis)
            a1zj = a1(3,jbis)
            a2xi = a2(1,i)
            a2yi = a2(2,i)
            a2zi = a2(3,i)
            a2xj = a2(1,jbis)
            a2yj = a2(2,jbis)
            a2zj = a2(3,jbis)

            !Errfunc damping
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0d0 * aewald**2
            alsq2n = 0.0d0
            if (aewald .gt. 0.0d0)
     &        alsq2n = 1.0d0 / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 4 
              bfac = dble(j+j-1)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a)*invrij2
            end do

            !Scalings
            damp = pdi*pdamp(jjpole)
            scale3 = 1d0
            scale5 = 1d0
            scale7 = 1d0
            scale9 = 1d0
            if (damp.ne.0d0) then
               pgamma = min(pti,thole(jjpole))
               damp = -pgamma*(d/damp)**3
               if (damp .gt. -50d0) then
                  expdamp = exp(damp)
                  scale3 = scale3 * (1.0d0 - expdamp)
                  scale5 = scale5 * (1.0d0 - expdamp*(1.0d0 - damp))
                  scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                  scale9 = 1d0 - expdamp*(1d0 - damp + (18d0*damp**2 -
     $                     9d0*damp**3)/35d0)
               end if
            end if
            dsc3 = scale3*dscale(jglob)
            dsc5 = scale5*dscale(jglob)
            dsc7 = scale7*dscale(jglob)
            dsc9 = scale9*dscale(jglob)
            psc3 = scale3*pscale(jglob)
            psc5 = scale5*pscale(jglob)
            psc7 = scale7*pscale(jglob)
            psc9 = scale9*pscale(jglob)
            !mat
            usc5 = scale5*uscale(jglob)
            usc7 = scale7*uscale(jglob)
            srr3(1) = (1d0 - dsc3)*rrij3 - bn(1)
            srr3(2) = (1d0 - psc3)*rrij3 - bn(1)
            srr5(1) = (1d0 - dsc5)*rrij5 - bn(2)
            srr5(2) = (1d0 - psc5)*rrij5 - bn(2)
            srr7(1) = (1d0 - dsc7)*rrij7 - bn(3)
            srr7(2) = (1d0 - psc7)*rrij7 - bn(3)
            srr9(1) = (1d0 - dsc9)*rrij9 - bn(4)
            srr9(2) = (1d0 - psc9)*rrij9 - bn(4)
            !mat
            urr5 = (1d0 - usc5)*rrij5 - bn(2)
            urr7 = (1d0 - usc7)*rrij7 - bn(3)

            !Actual field calculations ! 

c dexdxj (alpha=beta=1)
            term1 = rij(1)*rij(1)*srr5 - srr3 
            term2 = 2d0*srr5*rij(1) 
            term3 = -srr7*rij(1)*rij(1) + srr5
            term4 = -4d0*rij(1)*srr7
            term5 = srr9*rij(1)*rij(1) - srr7

            dexdxj = cj*term1 + dj(1)*term2
     $      + scalmujr*term3 + 2d0*qj(1,1)*srr5
     $      + thetajr(1)*term4 + quadterm*term5

            dexdxi = ci*term1 - di(1)*term2
     $      - scalmuir*term3 + 2d0*qi(1,1)*srr5
     $      + thetair(1)*term4 + quadterm_bis*term5
            
            demxdx = term1
            detxdx = 2d0*srr5-srr7*rij(1)*rij(1)

c deydyj (alpha=beta=2)
            term1 = rij(2)*rij(2)*srr5 - srr3 
            term2 = 2d0*srr5*rij(2) 
            term3 = -srr7*rij(2)*rij(2) + srr5
            term4 = -4d0*rij(2)*srr7
            term5 = srr9*rij(2)*rij(2) - srr7

            deydyj = cj*term1 + dj(2)*term2
     $      + scalmujr*term3 + 2d0*qj(2,2)*srr5
     $      + thetajr(2)*term4 + quadterm*term5

            deydyi = ci*term1 - di(2)*term2
     $      - scalmuir*term3 + 2d0*qi(2,2)*srr5
     $      + thetair(2)*term4 + quadterm_bis*term5

            demydy = term1
            detydy = 2d0*srr5-srr7*rij(2)*rij(2)

c dexdyj (alpha=1 beta=2)
            term1 = rij(1)*rij(2)*srr5 
            term2 = srr5*rij(2) 
            term3 = srr5*rij(1) 
            term4 = -srr7*rij(1)*rij(2)
            term5 = -2d0*rij(2)*srr7
            term6 = -2d0*rij(1)*srr7
            term7 = srr9*rij(1)*rij(2)

            dexdyj = cj*term1 + dj(1)*term2 + dj(2)*term3
     $      + scalmujr*term4 + 2d0*qj(1,2)*srr5
     $      + thetajr(1)*term5 +  thetajr(2)*term6
     $      + quadterm*term7

            dexdyi = ci*term1 - di(1)*term2 - di(2)*term3
     $      - scalmuir*term4 + 2d0*qi(1,2)*srr5
     $      + thetair(1)*term5 +  thetair(2)*term6
     $      + quadterm_bis*term7

            demxdy = term1
            detxdy = -srr7*rij(1)*rij(2)


c dezdzj (alpha3=beta=3)
            term1 = rij(3)*rij(3)*srr5 - srr3 
            term2 = 2d0*srr5*rij(3) 
            term3 = -srr7*rij(3)*rij(3) + srr5
            term4 = -4d0*rij(3)*srr7
            term5 = srr9*rij(3)*rij(3) - srr7

            dezdzj = cj*term1 + dj(3)*term2
     $      + scalmujr*term3 + 2d0*qj(3,3)*srr5
     $      + thetajr(3)*term4 + quadterm*term5

            dezdzi = ci*term1 - di(3)*term2
     $      - scalmuir*term3 + 2d0*qi(3,3)*srr5
     $      + thetair(3)*term4 + quadterm_bis*term5

            demzdz = term1
            detzdz = 2d0*srr5-srr7*rij(3)*rij(3)

c dexdzj (alpha=1 beta=3)
            term1 = rij(1)*rij(3)*srr5 
            term2 = srr5*rij(3) 
            term3 = srr5*rij(1) 
            term4 = -srr7*rij(1)*rij(3)
            term5 = -2d0*rij(3)*srr7
            term6 = -2d0*rij(1)*srr7
            term7 = srr9*rij(1)*rij(3)

            dexdzj = cj*term1 + dj(1)*term2 + dj(3)*term3
     $      + scalmujr*term4 + 2d0*qj(1,3)*srr5
     $      + thetajr(1)*term5 +  thetajr(3)*term6
     $      + quadterm*term7

            dexdzi = ci*term1 - di(1)*term2 - di(3)*term3
     $      - scalmuir*term4 + 2d0*qi(1,3)*srr5
     $      + thetair(1)*term5 +  thetair(3)*term6
     $      + quadterm_bis*term7

            demxdz = term1
            detxdz = -srr7*rij(1)*rij(3)


c deydzj (alpha=2 beta=3)
            term1 = rij(2)*rij(3)*srr5 
            term2 = srr5*rij(3) 
            term3 = srr5*rij(2) 
            term4 = -srr7*rij(2)*rij(3)
            term5 = -2d0*rij(3)*srr7
            term6 = -2d0*rij(2)*srr7
            term7 = srr9*rij(2)*rij(3)

            deydzj = cj*term1 + dj(2)*term2 + dj(3)*term3
     $      + scalmujr*term4 + 2d0*qj(2,3)*srr5
     $      + thetajr(2)*term5 +  thetajr(3)*term6
     $      + quadterm*term7

            deydzi = ci*term1 - di(2)*term2 - di(3)*term3
     $      - scalmuir*term4 + 2d0*qi(2,3)*srr5
     $      + thetair(2)*term5 +  thetair(3)*term6
     $      + quadterm_bis*term7

            demydz = term1
            detydz = -srr7*rij(2)*rij(3)

            depx = -a1xi*dexdxj(1) + a1xj*dexdxi(1)
     $     -a1yi*dexdyj(1) + a1yj*dexdyi(1)
     $     -a1zi*dexdzj(1) + a1zj*dexdzi(1)
            depy = -a1xi*dexdyj(1) + a1xj*dexdyi(1)
     $     -a1yi*deydyj(1) + a1yj*deydyi(1)
     $     -a1zi*deydzj(1) + a1zj*deydzi(1)
            depz = -a1xi*dexdzj(1) + a1xj*dexdzi(1)
     $     -a1yi*deydzj(1) + a1yj*deydzi(1)
     $     -a1zi*dezdzj(1) + a1zj*dezdzi(1)

            ade_cor(1,1,i) = ade_cor(1,1,i) + depx
            ade_cor(2,1,i) = ade_cor(2,1,i) + depy
            ade_cor(3,1,i) = ade_cor(3,1,i) + depz

            ade_cor(1,1,jbis) = ade_cor(1,1,jbis) - depx
            ade_cor(2,1,jbis) = ade_cor(2,1,jbis) - depy
            ade_cor(3,1,jbis) = ade_cor(3,1,jbis) - depz

            depx = -a2xi*dexdxj(2) + a2xj*dexdxi(2)
     $     -a2yi*dexdyj(2) + a2yj*dexdyi(2)
     $     -a2zi*dexdzj(2) + a2zj*dexdzi(2)
            depy = -a2xi*dexdyj(2) + a2xj*dexdyi(2)
     $     -a2yi*deydyj(2) + a2yj*deydyi(2)
     $     -a2zi*deydzj(2) + a2zj*deydzi(2)
            depz = -a2xi*dexdzj(2) + a2xj*dexdzi(2)
     $     -a2yi*deydzj(2) + a2yj*deydzi(2)
     $     -a2zi*dezdzj(2) + a2zj*dezdzi(2)

            ade_cor(1,2,i) = ade_cor(1,2,i) + depx
            ade_cor(2,2,i) = ade_cor(2,2,i) + depy
            ade_cor(3,2,i) = ade_cor(3,2,i) + depz

            ade_cor(1,2,jbis) = ade_cor(1,2,jbis) - depx
            ade_cor(2,2,jbis) = ade_cor(2,2,jbis) - depy
            ade_cor(3,2,jbis) = ade_cor(3,2,jbis) - depz


            dempx = a1xj*demxdx(1) + a1yj*demxdy(1)
     $     + a1zj*demxdz(1)
            dempy = a1xj*demxdy(1) + a1yj*demydy(1)
     $     + a1zj*demydz(1)
            dempz = a1xj*demxdz(1) + a1yj*demydz(1)
     $     + a1zj*demzdz(1)

            adme_cor(1,1,i) = adme_cor(1,1,i) + dempx
            adme_cor(2,1,i) = adme_cor(2,1,i) + dempy
            adme_cor(3,1,i) = adme_cor(3,1,i) + dempz

            dempx = a1xi*demxdx(1) + a1yi*demxdy(1)
     $     + a1zi*demxdz(1)
            dempy = a1xi*demxdy(1) + a1yi*demydy(1)
     $     + a1zi*demydz(1)
            dempz = a1xi*demxdz(1) + a1yi*demydz(1)
     $     + a1zi*demzdz(1)


            adme_cor(1,1,jbis) = adme_cor(1,1,jbis) + dempx
            adme_cor(2,1,jbis) = adme_cor(2,1,jbis) + dempy
            adme_cor(3,1,jbis) = adme_cor(3,1,jbis) + dempz

            dempx = a2xj*demxdx(2) + a2yj*demxdy(2)
     $     + a2zj*demxdz(2)
            dempy = a2xj*demxdy(2) + a2yj*demydy(2)
     $     + a2zj*demydz(2)
            dempz = a2xj*demxdz(2) + a2yj*demydz(2)
     $     + a2zj*demzdz(2)

            adme_cor(1,2,i) = adme_cor(1,2,i) + dempx
            adme_cor(2,2,i) = adme_cor(2,2,i) + dempy
            adme_cor(3,2,i) = adme_cor(3,2,i) + dempz

            dempx = a2xi*demxdx(2) + a2yi*demxdy(2)
     $     + a2zi*demxdz(2)
            dempy = a2xi*demxdy(2) + a2yi*demydy(2)
     $     + a2zi*demydz(2)
            dempz = a2xi*demxdz(2) + a2yi*demydz(2)
     $     + a2zi*demzdz(2)

            adme_cor(1,2,jbis) = adme_cor(1,2,jbis) + dempx
            adme_cor(2,2,jbis) = adme_cor(2,2,jbis) + dempy
            adme_cor(3,2,jbis) = adme_cor(3,2,jbis) + dempz

            detx = -a1xj*detxdx(1)-a1yj*detxdy(1)
     $     -a1zj*detxdz(1) 
            dety = -a1xj*detxdy(1)-a1yj*detydy(1)
     $     -a1zj*detydz(1) 
            detz = -a1xj*detxdz(1)-a1yj*detydz(1)
     $     -a1zj*detzdz(1) 

            adtE_cor(1,:,1,i) = adtE_cor(1,:,1,i) + detx*rij(:)
            adtE_cor(2,:,1,i) = adtE_cor(2,:,1,i) + dety*rij(:)
            adtE_cor(3,:,1,i) = adtE_cor(3,:,1,i) + detz*rij(:)


            detx = a1xi*detxdx(1)+a1yi*detxdy(1)
     $      + a1zi*detxdz(1) 
            dety = a1xi*detxdy(1)+a1yi*detydy(1)
     $      + a1zi*detydz(1) 
            detz = a1xi*detxdz(1)+a1yi*detydz(1)
     $      + a1zi*detzdz(1) 

            adtE_cor(1,:,1,jbis) = adtE_cor(1,:,1,jbis) + detx*rij(:)
            adtE_cor(2,:,1,jbis) = adtE_cor(2,:,1,jbis) + dety*rij(:)
            adtE_cor(3,:,1,jbis) = adtE_cor(3,:,1,jbis) + detz*rij(:)

            detx = -a2xj*detxdx(2)-a2yj*detxdy(2)
     $     -a2zj*detxdz(2) 
            dety = -a2xj*detxdy(2)-a2yj*detydy(2)
     $     -a2zj*detydz(2) 
            detz = -a2xj*detxdz(2)-a2yj*detydz(2)
     $     -a2zj*detzdz(2) 

            adtE_cor(1,:,2,i) = adtE_cor(1,:,2,i) + detx*rij(:)
            adtE_cor(2,:,2,i) = adtE_cor(2,:,2,i) + dety*rij(:)
            adtE_cor(3,:,2,i) = adtE_cor(3,:,2,i) + detz*rij(:)

            detx = a2xi*detxdx(2)+a2yi*detxdy(2)
     $      + a2zi*detxdz(2) 
            dety = a2xi*detxdy(2)+a2yi*detydy(2)
     $      + a2zi*detydz(2) 
            detz = a2xi*detxdz(2)+a2yi*detydz(2)
     $      + a2zi*detzdz(2) 

            adtE_cor(1,:,2,jbis) = adtE_cor(1,:,2,jbis) + detx*rij(:)
            adtE_cor(2,:,2,jbis) = adtE_cor(2,:,2,jbis) + dety*rij(:)
            adtE_cor(3,:,2,jbis) = adtE_cor(3,:,2,jbis) + detz*rij(:)
c
c      also compute gradient of field of dipoles "e"
c
c      alpha = beta = 1
c
            temp1jxx5 =  urr5*(scale1i + 2*rij(1)*e1(1,i))
            temp1jxx7 =  -urr7*scale1i*rij(1)*rij(1)
            temp1jxx = temp1jxx5 + temp1jxx7
            temp1ixx5 =  urr5*(scale1j + 2*rij(1)*e1(1,jbis))
            temp1ixx7 =  -urr7*scale1j*rij(1)*rij(1)
            temp1ixx = temp1ixx5 + temp1ixx7


c      alpha = 2, beta = 1
c
            temp1jxy5 =  urr5*(rij(2)*e1(1,i)+rij(1)*e1(2,i))
            temp1jxy7 =  -urr7*scale1i*rij(2)*rij(1)
            temp1jxy = temp1jxy5 + temp1jxy7
            temp1ixy5 =  urr5*(rij(2)*e1(1,jbis)+rij(1)*e1(2,jbis))
            temp1ixy7 =  -urr7*scale1j*rij(2)*rij(1)
            temp1ixy = temp1ixy5 + temp1ixy7

c
c      alpha = 3, beta = 1
c
            temp1jxz5 =  urr5*(rij(3)*e1(1,i)+rij(1)*e1(3,i))
            temp1jxz7 =  -urr7*scale1i*rij(3)*rij(1)
            temp1jxz = temp1jxz5 + temp1jxz7
            temp1ixz5 =  urr5*(rij(3)*e1(1,jbis)+rij(1)*e1(3,jbis))
            temp1ixz7 =  -urr7*scale1j*rij(3)*rij(1)
            temp1ixz = temp1ixz5 + temp1ixz7

c
c      alpha = beta = 2
c
            temp1jyy5 =  urr5*(scale1i + 2*rij(2)*e1(2,i))
            temp1jyy7 =  -urr7*scale1i*rij(2)*rij(2)
            temp1jyy = temp1jyy5 + temp1jyy7
            temp1iyy5 =  urr5*(scale1j + 2*rij(2)*e1(2,jbis))
            temp1iyy7 =  -urr7*scale1j*rij(2)*rij(2)
            temp1iyy = temp1iyy5 + temp1iyy7

c
c      alpha = beta = 3
c
            temp1jzz5 =  urr5*(scale1i + 2*rij(3)*e1(3,i))
            temp1jzz7 =  -urr7*scale1i*rij(3)*rij(3)
            temp1jzz = temp1jzz5 + temp1jzz7
            temp1izz5 =  urr5*(scale1j + 2*rij(3)*e1(3,jbis))
            temp1izz7 =  -urr7*scale1j*rij(3)*rij(3)
            temp1izz = temp1izz5 + temp1izz7


c
c      alpha = 2, beta = 3
c
            temp1jyz5 =  urr5*(rij(3)*e1(2,i)+rij(2)*e1(3,i))
            temp1jyz7 =  -urr7*scale1i*rij(2)*rij(3)
            temp1jyz = temp1jyz5 + temp1jyz7
            temp1iyz5 =  urr5*(rij(3)*e1(2,jbis)+rij(2)*e1(3,jbis))
            temp1iyz7 =  -urr7*scale1j*rij(2)*rij(3)
            temp1iyz = temp1iyz5 + temp1iyz7

c    
c      beta = 1
c
            depx =  
     $       d1(1,jbis)*temp1jxx + d1(2,jbis)*temp1jxy 
     $      + d1(3,jbis)*temp1jxz
     $      + d1(1,i)*temp1ixx + d1(2,i)*temp1ixy 
     $      + d1(3,i)*temp1ixz
            adTbcor(1,i) = adTbcor(1,i) + depx 
            adTbcor(1,jbis) = adTbcor(1,jbis) - depx 
c    
c      beta = 2
c
            depy =  
     $       d1(1,jbis)*temp1jxy + d1(2,jbis)*temp1jyy 
     $      + d1(3,jbis)*temp1jyz
     $      + d1(1,i)*temp1ixy + d1(2,i)*temp1iyy 
     $      + d1(3,i)*temp1iyz
            adTbcor(2,i) = adTbcor(2,i) + depy 
            adTbcor(2,jbis) = adTbcor(2,jbis) - depy 
c    
c      beta = 3
c
            depz =  
     $       d1(1,jbis)*temp1jxz + d1(2,jbis)*temp1jyz 
     $      + d1(3,jbis)*temp1jzz
     $      + d1(1,i)*temp1ixz + d1(2,i)*temp1iyz 
     $      + d1(3,i)*temp1izz
            adTbcor(3,i) = adTbcor(3,i) + depz 
            adTbcor(3,jbis) = adTbcor(3,jbis) - depz 

         end do !jj


         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = 1.0d0
            uscale(ip11(j,iglob)) = 1.0d0
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = 1.0d0
            uscale(ip12(j,iglob)) = 1.0d0
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = 1.0d0
            uscale(ip13(j,iglob)) = 1.0d0
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = 1.0d0
            uscale(ip14(j,iglob)) = 1.0d0
         end do
         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = 1.0d0
         end do
      end do

      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      adme_self(:,1,1:npoleloc) = term*a1(:,1:npoleloc)
      adme_self(:,2,1:npoleloc) = term*a2(:,1:npoleloc)

      ade  = - ade_cor  ! THIS -1 EVERYWHERE IS WEIRD
      adme = - adme_cor + adme_self ! But hey... 
      adte = - adte_cor !  ...it works.
      adtb = - adtbcor!

      deallocate(dscale, pscale, uscale)
      return
      end
c
      subroutine scalderfieldzmat8( a1, a2, d1, e1, d2, e2,d3,
     $                              e3, d4, e4, ade, adme, adte, adtb)
      use atmlst
      use atoms
      use bound
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
      use potent
      use shunt
      implicit none


      real*8, dimension(3,npolebloc), intent(in) :: a1, a2
      real*8, dimension(3,npolebloc), intent(in) :: e1, e2, e3,e4, d1,
     $                                               d2,d3,d4
      real*8, dimension(3,npolebloc), intent(out) :: adtb
      real*8, dimension(3,2,npolebloc), intent(out) :: adme, ade
      real*8, dimension(3,3,2,npolebloc), intent(out) :: adte
      
      real*8, allocatable, dimension(:,:,:) ::adE_cor,admE_cor,admE_self
      real*8, allocatable, dimension(:,:,:,:) ::adtE_cor




      integer :: ii, iipole, iglob, i, jjpole, jglob,
     $           jjj, jbis, k, j

      real*8 :: pti, pdi, rrij3, rrij5, rrij7, rrij9,
     $          urr5, urr7, 
     $scale3, scale5, 
     $          scale7, scale9, dsc3, dsc5, dsc7, dsc9, psc3, psc5,
     $          psc7, psc9, usc5, usc7, scalmuir, scalmujr, 
     $quadterm, quadterm_bis, 
     $          ci, cj, bfac, invrij2, alsq2n, cutoff2, exp2a, ralpha,
     $          d, rij2, damp, expdamp, pgamma, alsq2, term

      real*8                :: bn(0:4) 
      real*8, dimension(2)  :: srr3, srr5, srr7, srr9
      real*8, dimension(3)  :: rij, thetajr, thetair, di, dj
      real*8, dimension(3,3):: qi, qj

      real*8, allocatable, dimension(:) :: dscale, pscale, uscale
      real*8, allocatable, dimension(:,:) :: adtbcor
      

      real*8:: a1xi,a1yi,a1zi,a1xj,a1yj,a1zj
      real*8:: a2xi,a2yi,a2zi,a2xj,a2yj,a2zj
      real*8,dimension(2):: term1,term2,term3,term4,term5,term6,term7
      real*8,dimension(2):: dexdxi,dexdxj,deydyi,deydyj,dezdzi,dezdzj
      real*8,dimension(2):: dexdyi,dexdyj,dexdzi,dexdzj,deydzi,deydzj
      real*8,dimension(2):: demxdx,demxdy,demxdz,demydy,demydz,demzdz
      real*8,dimension(2):: detxdx,detxdy,detxdz,detydy,detydz,detzdz
      real*8 :: depx,depy,depz,dempx,dempy,dempz,detx,dety,detz
      real*8 :: scale1i,scale1j,scale2i,scale2j,scale3i,scale3j
      real*8 :: scale4i,scale4j
      real*8::temp1jxx5,temp1jxy5,temp1jxz5,temp1jyy5,temp1jyz5,
     $     temp1jzz5
      real*8::temp1jxx7,temp1jxy7,temp1jxz7,temp1jyy7,temp1jyz7,
     $     temp1jzz7
      real*8::temp1jxx,temp1jxy,temp1jxz,temp1jyy,temp1jyz,temp1jzz
      real*8::temp1ixx5,temp1ixy5,temp1ixz5,temp1iyy5,temp1iyz5,
     $     temp1izz5
      real*8::temp1ixx7,temp1ixy7,temp1ixz7,temp1iyy7,temp1iyz7,
     $     temp1izz7
      real*8::temp1ixx,temp1ixy,temp1ixz,temp1iyy,temp1iyz,temp1izz

      real*8::temp2jxx5,temp2jxy5,temp2jxz5,temp2jyy5,temp2jyz5,
     $     temp2jzz5
      real*8::temp2jxx7,temp2jxy7,temp2jxz7,temp2jyy7,temp2jyz7,
     $     temp2jzz7
      real*8::temp2jxx,temp2jxy,temp2jxz,temp2jyy,temp2jyz,temp2jzz
      real*8::temp2ixx5,temp2ixy5,temp2ixz5,temp2iyy5,temp2iyz5,
     $     temp2izz5
      real*8::temp2ixx7,temp2ixy7,temp2ixz7,temp2iyy7,temp2iyz7,
     $     temp2izz7
      real*8::temp2ixx,temp2ixy,temp2ixz,temp2iyy,temp2iyz,temp2izz

      real*8::temp3jxx5,temp3jxy5,temp3jxz5,temp3jyy5,temp3jyz5,
     $     temp3jzz5
      real*8::temp3jxx7,temp3jxy7,temp3jxz7,temp3jyy7,temp3jyz7,
     $     temp3jzz7
      real*8::temp3jxx,temp3jxy,temp3jxz,temp3jyy,temp3jyz,temp3jzz
      real*8::temp3ixx5,temp3ixy5,temp3ixz5,temp3iyy5,temp3iyz5,
     $     temp3izz5
      real*8::temp3ixx7,temp3ixy7,temp3ixz7,temp3iyy7,temp3iyz7,
     $     temp3izz7
      real*8::temp3ixx,temp3ixy,temp3ixz,temp3iyy,temp3iyz,temp3izz

      real*8::temp4jxx5,temp4jxy5,temp4jxz5,temp4jyy5,temp4jyz5,
     $     temp4jzz5
      real*8::temp4jxx7,temp4jxy7,temp4jxz7,temp4jyy7,temp4jyz7,
     $     temp4jzz7
      real*8::temp4jxx,temp4jxy,temp4jxz,temp4jyy,temp4jyz,temp4jzz
      real*8::temp4ixx5,temp4ixy5,temp4ixz5,temp4iyy5,temp4iyz5,
     $     temp4izz5
      real*8::temp4ixx7,temp4ixy7,temp4ixz7,temp4iyy7,temp4iyz7,
     $     temp4izz7
      real*8::temp4ixx,temp4ixy,temp4ixz,temp4iyy,temp4iyz,temp4izz

      character(10) :: mode

 1000 format(' Warning, system moved too much since last neighbor list
     $   update, try lowering nlupdate')

      

      ade = 0d0
      adme = 0d0
      adte = 0d0
      adtb = 0d0

      mode = 'EWALD'
      if (use_polarshortreal) mode = 'SHORTEWALD'
      call switch (mode)
      cutoff2 = cut2

      allocate(dscale(n), pscale(n), uscale(n))
      allocate(adtbcor(3,npolebloc))
      allocate(adE_cor(3,2,npolebloc),admE_cor(3,2,npolebloc))
      allocate(admE_self(3,2,npolebloc),adtE_cor(3,3,2,npolebloc))
      adme_cor = 0d0
      adme_self = 0d0
      adte_cor = 0d0
      ade_cor = 0d0
      dscale = 1d0
      pscale = 1d0
      uscale = 1.0d0
      adtbcor = 0d0

      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = poleloc(iglob)
         if ((i.eq.0).or.(i.gt.nbloc)) then
            write(iout,1000)
            cycle
         end if
         pdi = pdamp(iipole)
         pti = thole(iipole)

         ci      = rpole(1,iipole)
         di(1)   = rpole(2,iipole)
         di(2)   = rpole(3,iipole)
         di(3)   = rpole(4,iipole)
         qi(1,1) = rpole(5, iipole)
         qi(2,1) = rpole(6, iipole)
         qi(3,1) = rpole(7, iipole)
         qi(1,2) = rpole(6, iipole)
         qi(2,2) = rpole(9, iipole)
         qi(3,2) = rpole(10,iipole)
         qi(1,3) = rpole(7, iipole)
         qi(2,3) = rpole(10,iipole)
         qi(3,3) = rpole(13,iipole)

         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = p2scale
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = p3scale
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = p4scale
            do k = 1, np11(iglob)
               if (i14(j,iglob) .eq. ip11(k,iglob))
     &            pscale(i14(j,iglob)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = p5scale
         end do
         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = d1scale
            uscale(ip11(j,iglob)) = u1scale
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = d2scale
            uscale(ip12(j,iglob)) = u2scale
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = d3scale
            uscale(ip13(j,iglob)) = u3scale
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = d4scale
            uscale(ip14(j,iglob)) = u4scale
         end do

         do jjj = 1, nelst(ii)
            jjpole = elst(jjj,ii)
            jglob = ipole(jjpole)
            jbis = poleloc(jjpole)
            if (jbis.eq.0) then
               write(iout,1000)
               cycle
            end if

            !Distances
            rij(1) = x(jglob) - x(iglob)
            rij(2) = y(jglob) - y(iglob)
            rij(3) = z(jglob) - z(iglob)
            call image(rij(1), rij(2), rij(3))
            rij2 = dot_product(rij,rij)
            if (rij2 .gt. cutoff2) cycle
            d  = sqrt(rij2)
            invrij2 = 1d0/rij2
            rrij3 = invrij2/d
            rrij5 = 3d0*rrij3*invrij2
            rrij7 = 5d0*rrij5*invrij2
            rrij9 = 7d0*rrij7*invrij2

            !Multipoles
            cj      = rpole(1,jjpole)
            dj(1)   = rpole(2,jjpole)
            dj(2)   = rpole(3,jjpole)
            dj(3)   = rpole(4,jjpole)
            qj(1,1) = rpole(5, jjpole)
            qj(2,1) = rpole(6, jjpole)
            qj(3,1) = rpole(7, jjpole)
            qj(1,2) = rpole(6, jjpole)
            qj(2,2) = rpole(9, jjpole)
            qj(3,2) = rpole(10,jjpole)
            qj(1,3) = rpole(7, jjpole)
            qj(2,3) = rpole(10,jjpole)
            qj(3,3) = rpole(13,jjpole)

            scalmujr = dot_product(dj,rij)!sprod(3,dj,rij)
            scalmuir = dot_product(di,rij)!sprod(3,di,rij)
            thetajr(1) = dot_product(qj(1,:),rij)
            thetajr(2) = dot_product(qj(2,:),rij)
            thetajr(3) = dot_product(qj(3,:),rij)
            thetair(1) = dot_product(qi(1,:),rij)
            thetair(2) = dot_product(qi(2,:),rij)
            thetair(3) = dot_product(qi(3,:),rij)
            quadterm     = dot_product(thetajr,rij)!sprod(3,thetajr,rij)
            quadterm_bis = dot_product(thetair,rij)!sprod(3,thetair,rij)

            scale1j = dot_product(e1(:,jbis),rij)
            scale1i = dot_product(e1(:,i),rij)
            scale2j = dot_product(e2(:,jbis),rij)
            scale2i = dot_product(e2(:,i),rij)
            scale3j = dot_product(e3(:,jbis),rij)
            scale3i = dot_product(e3(:,i),rij)
            scale4j = dot_product(e4(:,jbis),rij)
            scale4i = dot_product(e4(:,i),rij)

            a1xi = a1(1,i)
            a1yi = a1(2,i)
            a1zi = a1(3,i)
            a1xj = a1(1,jbis)
            a1yj = a1(2,jbis)
            a1zj = a1(3,jbis)
            a2xi = a2(1,i)
            a2yi = a2(2,i)
            a2zi = a2(3,i)
            a2xj = a2(1,jbis)
            a2yj = a2(2,jbis)
            a2zj = a2(3,jbis)

            !Errfunc damping
            ralpha = aewald * d
            bn(0) = erfc(ralpha) / d
            alsq2 = 2.0d0 * aewald**2
            alsq2n = 0.0d0
            if (aewald .gt. 0.0d0)
     &        alsq2n = 1.0d0 / (sqrtpi*aewald)
            exp2a = exp(-ralpha**2)
            do j = 1, 4 
              bfac = dble(j+j-1)
              alsq2n = alsq2 * alsq2n
              bn(j) = (bfac*bn(j-1)+alsq2n*exp2a)*invrij2
            end do

            !Scalings
            damp = pdi*pdamp(jjpole)
            scale3 = 1d0
            scale5 = 1d0
            scale7 = 1d0
            scale9 = 1d0
            if (damp.ne.0d0) then
               pgamma = min(pti,thole(jjpole))
               damp = -pgamma*(d/damp)**3
               if (damp .gt. -50d0) then
                  expdamp = exp(damp)
                  scale3 = scale3 * (1.0d0 - expdamp)
                  scale5 = scale5 * (1.0d0 - expdamp*(1.0d0 - damp))
                  scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                  scale9 = 1d0 - expdamp*(1d0 - damp + (18d0*damp**2 -
     $                     9d0*damp**3)/35d0)
               end if
            end if
            dsc3 = scale3*dscale(jglob)
            dsc5 = scale5*dscale(jglob)
            dsc7 = scale7*dscale(jglob)
            dsc9 = scale9*dscale(jglob)
            psc3 = scale3*pscale(jglob)
            psc5 = scale5*pscale(jglob)
            psc7 = scale7*pscale(jglob)
            psc9 = scale9*pscale(jglob)
            !mat
            usc5 = scale5*uscale(jglob)
            usc7 = scale7*uscale(jglob)
            srr3(1) = (1d0 - dsc3)*rrij3 - bn(1)
            srr3(2) = (1d0 - psc3)*rrij3 - bn(1)
            srr5(1) = (1d0 - dsc5)*rrij5 - bn(2)
            srr5(2) = (1d0 - psc5)*rrij5 - bn(2)
            srr7(1) = (1d0 - dsc7)*rrij7 - bn(3)
            srr7(2) = (1d0 - psc7)*rrij7 - bn(3)
            srr9(1) = (1d0 - dsc9)*rrij9 - bn(4)
            srr9(2) = (1d0 - psc9)*rrij9 - bn(4)
            !mat
            urr5 = (1d0 - usc5)*rrij5 - bn(2)
            urr7 = (1d0 - usc7)*rrij7 - bn(3)

            !Actual field calculations ! 



c dexdxj (alpha=beta=1)
            term1 = rij(1)*rij(1)*srr5 - srr3 
            term2 = 2d0*srr5*rij(1) 
            term3 = -srr7*rij(1)*rij(1) + srr5
            term4 = -4d0*rij(1)*srr7
            term5 = srr9*rij(1)*rij(1) - srr7

            dexdxj = cj*term1 + dj(1)*term2
     $      + scalmujr*term3 + 2d0*qj(1,1)*srr5
     $      + thetajr(1)*term4 + quadterm*term5

            dexdxi = ci*term1 - di(1)*term2
     $      - scalmuir*term3 + 2d0*qi(1,1)*srr5
     $      + thetair(1)*term4 + quadterm_bis*term5
            
            demxdx = term1
            detxdx = 2d0*srr5-srr7*rij(1)*rij(1)

c deydyj (alpha=beta=2)
            term1 = rij(2)*rij(2)*srr5 - srr3 
            term2 = 2d0*srr5*rij(2) 
            term3 = -srr7*rij(2)*rij(2) + srr5
            term4 = -4d0*rij(2)*srr7
            term5 = srr9*rij(2)*rij(2) - srr7

            deydyj = cj*term1 + dj(2)*term2
     $      + scalmujr*term3 + 2d0*qj(2,2)*srr5
     $      + thetajr(2)*term4 + quadterm*term5

            deydyi = ci*term1 - di(2)*term2
     $      - scalmuir*term3 + 2d0*qi(2,2)*srr5
     $      + thetair(2)*term4 + quadterm_bis*term5

            demydy = term1
            detydy = 2d0*srr5-srr7*rij(2)*rij(2)

c dexdyj (alpha=1 beta=2)
            term1 = rij(1)*rij(2)*srr5 
            term2 = srr5*rij(2) 
            term3 = srr5*rij(1) 
            term4 = -srr7*rij(1)*rij(2)
            term5 = -2d0*rij(2)*srr7
            term6 = -2d0*rij(1)*srr7
            term7 = srr9*rij(1)*rij(2)

            dexdyj = cj*term1 + dj(1)*term2 + dj(2)*term3
     $      + scalmujr*term4 + 2d0*qj(1,2)*srr5
     $      + thetajr(1)*term5 +  thetajr(2)*term6
     $      + quadterm*term7

            dexdyi = ci*term1 - di(1)*term2 - di(2)*term3
     $      - scalmuir*term4 + 2d0*qi(1,2)*srr5
     $      + thetair(1)*term5 +  thetair(2)*term6
     $      + quadterm_bis*term7

            demxdy = term1
            detxdy = -srr7*rij(1)*rij(2)


c dezdzj (alpha3=beta=3)
            term1 = rij(3)*rij(3)*srr5 - srr3 
            term2 = 2d0*srr5*rij(3) 
            term3 = -srr7*rij(3)*rij(3) + srr5
            term4 = -4d0*rij(3)*srr7
            term5 = srr9*rij(3)*rij(3) - srr7

            dezdzj = cj*term1 + dj(3)*term2
     $      + scalmujr*term3 + 2d0*qj(3,3)*srr5
     $      + thetajr(3)*term4 + quadterm*term5
c            write(*,*) 'i = ',iglob,'j = ',jglob,dezdzj

            dezdzi = ci*term1 - di(3)*term2
     $      - scalmuir*term3 + 2d0*qi(3,3)*srr5
     $      + thetair(3)*term4 + quadterm_bis*term5

            demzdz = term1
            detzdz = 2d0*srr5-srr7*rij(3)*rij(3)

c dexdzj (alpha=1 beta=3)
            term1 = rij(1)*rij(3)*srr5 
            term2 = srr5*rij(3) 
            term3 = srr5*rij(1) 
            term4 = -srr7*rij(1)*rij(3)
            term5 = -2d0*rij(3)*srr7
            term6 = -2d0*rij(1)*srr7
            term7 = srr9*rij(1)*rij(3)

            dexdzj = cj*term1 + dj(1)*term2 + dj(3)*term3
     $      + scalmujr*term4 + 2d0*qj(1,3)*srr5
     $      + thetajr(1)*term5 +  thetajr(3)*term6
     $      + quadterm*term7

            dexdzi = ci*term1 - di(1)*term2 - di(3)*term3
     $      - scalmuir*term4 + 2d0*qi(1,3)*srr5
     $      + thetair(1)*term5 +  thetair(3)*term6
     $      + quadterm_bis*term7

            demxdz = term1
            detxdz = -srr7*rij(1)*rij(3)


c deydzj (alpha=2 beta=3)
            term1 = rij(2)*rij(3)*srr5 
            term2 = srr5*rij(3) 
            term3 = srr5*rij(2) 
            term4 = -srr7*rij(2)*rij(3)
            term5 = -2d0*rij(3)*srr7
            term6 = -2d0*rij(2)*srr7
            term7 = srr9*rij(2)*rij(3)

            deydzj = cj*term1 + dj(2)*term2 + dj(3)*term3
     $      + scalmujr*term4 + 2d0*qj(2,3)*srr5
     $      + thetajr(2)*term5 +  thetajr(3)*term6
     $      + quadterm*term7

            deydzi = ci*term1 - di(2)*term2 - di(3)*term3
     $      - scalmuir*term4 + 2d0*qi(2,3)*srr5
     $      + thetair(2)*term5 +  thetair(3)*term6
     $      + quadterm_bis*term7

            demydz = term1
            detydz = -srr7*rij(2)*rij(3)

            depx = -a1xi*dexdxj(1) + a1xj*dexdxi(1)
     $     -a1yi*dexdyj(1) + a1yj*dexdyi(1)
     $     -a1zi*dexdzj(1) + a1zj*dexdzi(1)
            depy = -a1xi*dexdyj(1) + a1xj*dexdyi(1)
     $     -a1yi*deydyj(1) + a1yj*deydyi(1)
     $     -a1zi*deydzj(1) + a1zj*deydzi(1)
            depz = -a1xi*dexdzj(1) + a1xj*dexdzi(1)
     $     -a1yi*deydzj(1) + a1yj*deydzi(1)
     $     -a1zi*dezdzj(1) + a1zj*dezdzi(1)

            ade_cor(1,1,i) = ade_cor(1,1,i) + depx
            ade_cor(2,1,i) = ade_cor(2,1,i) + depy
            ade_cor(3,1,i) = ade_cor(3,1,i) + depz

            ade_cor(1,1,jbis) = ade_cor(1,1,jbis) - depx
            ade_cor(2,1,jbis) = ade_cor(2,1,jbis) - depy
            ade_cor(3,1,jbis) = ade_cor(3,1,jbis) - depz

            depx = -a2xi*dexdxj(2) + a2xj*dexdxi(2)
     $     -a2yi*dexdyj(2) + a2yj*dexdyi(2)
     $     -a2zi*dexdzj(2) + a2zj*dexdzi(2)
            depy = -a2xi*dexdyj(2) + a2xj*dexdyi(2)
     $     -a2yi*deydyj(2) + a2yj*deydyi(2)
     $     -a2zi*deydzj(2) + a2zj*deydzi(2)
            depz = -a2xi*dexdzj(2) + a2xj*dexdzi(2)
     $     -a2yi*deydzj(2) + a2yj*deydzi(2)
     $     -a2zi*dezdzj(2) + a2zj*dezdzi(2)

            ade_cor(1,2,i) = ade_cor(1,2,i) + depx
            ade_cor(2,2,i) = ade_cor(2,2,i) + depy
            ade_cor(3,2,i) = ade_cor(3,2,i) + depz

            ade_cor(1,2,jbis) = ade_cor(1,2,jbis) - depx
            ade_cor(2,2,jbis) = ade_cor(2,2,jbis) - depy
            ade_cor(3,2,jbis) = ade_cor(3,2,jbis) - depz

            dempx = a1xj*demxdx(1) + a1yj*demxdy(1)
     $     + a1zj*demxdz(1)
            dempy = a1xj*demxdy(1) + a1yj*demydy(1)
     $     + a1zj*demydz(1)
            dempz = a1xj*demxdz(1) + a1yj*demydz(1)
     $     + a1zj*demzdz(1)

            adme_cor(1,1,i) = adme_cor(1,1,i) + dempx
            adme_cor(2,1,i) = adme_cor(2,1,i) + dempy
            adme_cor(3,1,i) = adme_cor(3,1,i) + dempz

            dempx = a1xi*demxdx(1) + a1yi*demxdy(1)
     $     + a1zi*demxdz(1)
            dempy = a1xi*demxdy(1) + a1yi*demydy(1)
     $     + a1zi*demydz(1)
            dempz = a1xi*demxdz(1) + a1yi*demydz(1)
     $     + a1zi*demzdz(1)

            adme_cor(1,1,jbis) = adme_cor(1,1,jbis) + dempx
            adme_cor(2,1,jbis) = adme_cor(2,1,jbis) + dempy
            adme_cor(3,1,jbis) = adme_cor(3,1,jbis) + dempz

            dempx = a2xj*demxdx(2) + a2yj*demxdy(2)
     $     + a2zj*demxdz(2)
            dempy = a2xj*demxdy(2) + a2yj*demydy(2)
     $     + a2zj*demydz(2)
            dempz = a2xj*demxdz(2) + a2yj*demydz(2)
     $     + a2zj*demzdz(2)

            adme_cor(1,2,i) = adme_cor(1,2,i) + dempx
            adme_cor(2,2,i) = adme_cor(2,2,i) + dempy
            adme_cor(3,2,i) = adme_cor(3,2,i) + dempz

            dempx = a2xi*demxdx(2) + a2yi*demxdy(2)
     $     + a2zi*demxdz(2)
            dempy = a2xi*demxdy(2) + a2yi*demydy(2)
     $     + a2zi*demydz(2)
            dempz = a2xi*demxdz(2) + a2yi*demydz(2)
     $     + a2zi*demzdz(2)

            adme_cor(1,2,jbis) = adme_cor(1,2,jbis) + dempx
            adme_cor(2,2,jbis) = adme_cor(2,2,jbis) + dempy
            adme_cor(3,2,jbis) = adme_cor(3,2,jbis) + dempz

            detx = -a1xj*detxdx(1)-a1yj*detxdy(1)
     $     -a1zj*detxdz(1) 
            dety = -a1xj*detxdy(1)-a1yj*detydy(1)
     $     -a1zj*detydz(1) 
            detz = -a1xj*detxdz(1)-a1yj*detydz(1)
     $     -a1zj*detzdz(1) 

            adtE_cor(1,:,1,i) = adtE_cor(1,:,1,i) + detx*rij(:)
            adtE_cor(2,:,1,i) = adtE_cor(2,:,1,i) + dety*rij(:)
            adtE_cor(3,:,1,i) = adtE_cor(3,:,1,i) + detz*rij(:)


            detx = a1xi*detxdx(1)+a1yi*detxdy(1)
     $      + a1zi*detxdz(1) 
            dety = a1xi*detxdy(1)+a1yi*detydy(1)
     $      + a1zi*detydz(1) 
            detz = a1xi*detxdz(1)+a1yi*detydz(1)
     $      + a1zi*detzdz(1) 

            adtE_cor(1,:,1,jbis) = adtE_cor(1,:,1,jbis) + detx*rij(:)
            adtE_cor(2,:,1,jbis) = adtE_cor(2,:,1,jbis) + dety*rij(:)
            adtE_cor(3,:,1,jbis) = adtE_cor(3,:,1,jbis) + detz*rij(:)

            detx = -a2xj*detxdx(2)-a2yj*detxdy(2)
     $     -a2zj*detxdz(2) 
            dety = -a2xj*detxdy(2)-a2yj*detydy(2)
     $     -a2zj*detydz(2) 
            detz = -a2xj*detxdz(2)-a2yj*detydz(2)
     $     -a2zj*detzdz(2) 

            adtE_cor(1,:,2,i) = adtE_cor(1,:,2,i) + detx*rij(:)
            adtE_cor(2,:,2,i) = adtE_cor(2,:,2,i) + dety*rij(:)
            adtE_cor(3,:,2,i) = adtE_cor(3,:,2,i) + detz*rij(:)

            detx = a2xi*detxdx(2)+a2yi*detxdy(2)
     $      + a2zi*detxdz(2) 
            dety = a2xi*detxdy(2)+a2yi*detydy(2)
     $      + a2zi*detydz(2) 
            detz = a2xi*detxdz(2)+a2yi*detydz(2)
     $      + a2zi*detzdz(2) 

            adtE_cor(1,:,2,jbis) = adtE_cor(1,:,2,jbis) + detx*rij(:)
            adtE_cor(2,:,2,jbis) = adtE_cor(2,:,2,jbis) + dety*rij(:)
            adtE_cor(3,:,2,jbis) = adtE_cor(3,:,2,jbis) + detz*rij(:)
c
c      also compute gradient of field of dipoles "e"
c
c      alpha = beta = 1
c
            temp1jxx5 =  urr5*(scale1i + 2*rij(1)*e1(1,i))
            temp1jxx7 =  -urr7*scale1i*rij(1)*rij(1)
            temp1jxx = temp1jxx5 + temp1jxx7
            temp1ixx5 =  urr5*(scale1j + 2*rij(1)*e1(1,jbis))
            temp1ixx7 =  -urr7*scale1j*rij(1)*rij(1)
            temp1ixx = temp1ixx5 + temp1ixx7

            temp2jxx5 =  urr5*(scale2i + 2*rij(1)*e2(1,i))
            temp2jxx7 =  -urr7*scale2i*rij(1)*rij(1)
            temp2jxx = temp2jxx5 + temp2jxx7
            temp2ixx5 =  urr5*(scale2j + 2*rij(1)*e2(1,jbis))
            temp2ixx7 =  -urr7*scale2j*rij(1)*rij(1)
            temp2ixx = temp2ixx5 + temp2ixx7

            temp3jxx5 =  urr5*(scale3i + 2*rij(1)*e3(1,i))
            temp3jxx7 =  -urr7*scale3i*rij(1)*rij(1)
            temp3jxx = temp3jxx5 + temp3jxx7
            temp3ixx5 =  urr5*(scale3j + 2*rij(1)*e3(1,jbis))
            temp3ixx7 =  -urr7*scale3j*rij(1)*rij(1)
            temp3ixx = temp3ixx5 + temp3ixx7

            temp4jxx5 =  urr5*(scale4i + 2*rij(1)*e4(1,i))
            temp4jxx7 =  -urr7*scale4i*rij(1)*rij(1)
            temp4jxx = temp4jxx5 + temp4jxx7
            temp4ixx5 =  urr5*(scale4j + 2*rij(1)*e4(1,jbis))
            temp4ixx7 =  -urr7*scale4j*rij(1)*rij(1)
            temp4ixx = temp4ixx5 + temp4ixx7

c
c      alpha = 2, beta = 1
c
            temp1jxy5 =  urr5*(rij(2)*e1(1,i)+rij(1)*e1(2,i))
            temp1jxy7 =  -urr7*scale1i*rij(2)*rij(1)
            temp1jxy = temp1jxy5 + temp1jxy7
            temp1ixy5 =  urr5*(rij(2)*e1(1,jbis)+rij(1)*e1(2,jbis))
            temp1ixy7 =  -urr7*scale1j*rij(2)*rij(1)
            temp1ixy = temp1ixy5 + temp1ixy7

            temp2jxy5 =  urr5*(rij(2)*e2(1,i)+rij(1)*e2(2,i))
            temp2jxy7 =  -urr7*scale2i*rij(2)*rij(1)
            temp2jxy = temp2jxy5 + temp2jxy7
            temp2ixy5 =  urr5*(rij(2)*e2(1,jbis)+rij(1)*e2(2,jbis))
            temp2ixy7 =  -urr7*scale2j*rij(2)*rij(1)
            temp2ixy = temp2ixy5 + temp2ixy7

            temp3jxy5 =  urr5*(rij(2)*e3(1,i)+rij(1)*e3(2,i))
            temp3jxy7 =  -urr7*scale3i*rij(2)*rij(1)
            temp3jxy = temp3jxy5 + temp3jxy7
            temp3ixy5 =  urr5*(rij(2)*e3(1,jbis)+rij(1)*e3(2,jbis))
            temp3ixy7 =  -urr7*scale3j*rij(2)*rij(1)
            temp3ixy = temp3ixy5 + temp3ixy7

            temp4jxy5 =  urr5*(rij(2)*e4(1,i)+rij(1)*e4(2,i))
            temp4jxy7 =  -urr7*scale4i*rij(2)*rij(1)
            temp4jxy = temp4jxy5 + temp4jxy7
            temp4ixy5 =  urr5*(rij(2)*e4(1,jbis)+rij(1)*e4(2,jbis))
            temp4ixy7 =  -urr7*scale4j*rij(2)*rij(1)
            temp4ixy = temp4ixy5 + temp4ixy7
            
c
c      alpha = 3, beta = 1
c
            temp1jxz5 =  urr5*(rij(3)*e1(1,i)+rij(1)*e1(3,i))
            temp1jxz7 =  -urr7*scale1i*rij(3)*rij(1)
            temp1jxz = temp1jxz5 + temp1jxz7
            temp1ixz5 =  urr5*(rij(3)*e1(1,jbis)+rij(1)*e1(3,jbis))
            temp1ixz7 =  -urr7*scale1j*rij(3)*rij(1)
            temp1ixz = temp1ixz5 + temp1ixz7

            temp2jxz5 =  urr5*(rij(3)*e2(1,i)+rij(1)*e2(3,i))
            temp2jxz7 =  -urr7*scale2i*rij(3)*rij(1)
            temp2jxz = temp2jxz5 + temp2jxz7
            temp2ixz5 =  urr5*(rij(3)*e2(1,jbis)+rij(1)*e2(3,jbis))
            temp2ixz7 =  -urr7*scale2j*rij(3)*rij(1)
            temp2ixz = temp2ixz5 + temp2ixz7

            temp3jxz5 =  urr5*(rij(3)*e3(1,i)+rij(1)*e3(3,i))
            temp3jxz7 =  -urr7*scale3i*rij(3)*rij(1)
            temp3jxz = temp3jxz5 + temp3jxz7
            temp3ixz5 =  urr5*(rij(3)*e3(1,jbis)+rij(1)*e3(3,jbis))
            temp3ixz7 =  -urr7*scale3j*rij(3)*rij(1)
            temp3ixz = temp3ixz5 + temp3ixz7


            temp4jxz5 =  urr5*(rij(3)*e4(1,i)+rij(1)*e4(3,i))
            temp4jxz7 =  -urr7*scale4i*rij(3)*rij(1)
            temp4jxz = temp4jxz5 + temp4jxz7
            temp4ixz5 =  urr5*(rij(3)*e4(1,jbis)+rij(1)*e4(3,jbis))
            temp4ixz7 =  -urr7*scale4j*rij(3)*rij(1)
            temp4ixz = temp4ixz5 + temp4ixz7

c
c      alpha = beta = 2
c
            temp1jyy5 =  urr5*(scale1i + 2*rij(2)*e1(2,i))
            temp1jyy7 =  -urr7*scale1i*rij(2)*rij(2)
            temp1jyy = temp1jyy5 + temp1jyy7
            temp1iyy5 =  urr5*(scale1j + 2*rij(2)*e1(2,jbis))
            temp1iyy7 =  -urr7*scale1j*rij(2)*rij(2)
            temp1iyy = temp1iyy5 + temp1iyy7

            temp2jyy5 =  urr5*(scale2i + 2*rij(2)*e2(2,i))
            temp2jyy7 =  -urr7*scale2i*rij(2)*rij(2)
            temp2jyy = temp2jyy5 + temp2jyy7
            temp2iyy5 =  urr5*(scale2j + 2*rij(2)*e2(2,jbis))
            temp2iyy7 =  -urr7*scale2j*rij(2)*rij(2)
            temp2iyy = temp2iyy5 + temp2iyy7

            temp3jyy5 =  urr5*(scale3i + 2*rij(2)*e3(2,i))
            temp3jyy7 =  -urr7*scale3i*rij(2)*rij(2)
            temp3jyy = temp3jyy5 + temp3jyy7
            temp3iyy5 =  urr5*(scale3j + 2*rij(2)*e3(2,jbis))
            temp3iyy7 =  -urr7*scale3j*rij(2)*rij(2)
            temp3iyy = temp3iyy5 + temp3iyy7

            temp4jyy5 =  urr5*(scale4i + 2*rij(2)*e4(2,i))
            temp4jyy7 =  -urr7*scale4i*rij(2)*rij(2)
            temp4jyy = temp4jyy5 + temp4jyy7
            temp4iyy5 =  urr5*(scale4j + 2*rij(2)*e4(2,jbis))
            temp4iyy7 =  -urr7*scale4j*rij(2)*rij(2)
            temp4iyy = temp4iyy5 + temp4iyy7

c
c      alpha = beta = 3
c
            temp1jzz5 =  urr5*(scale1i + 2*rij(3)*e1(3,i))
            temp1jzz7 =  -urr7*scale1i*rij(3)*rij(3)
            temp1jzz = temp1jzz5 + temp1jzz7
            temp1izz5 =  urr5*(scale1j + 2*rij(3)*e1(3,jbis))
            temp1izz7 =  -urr7*scale1j*rij(3)*rij(3)
            temp1izz = temp1izz5 + temp1izz7

            temp2jzz5 =  urr5*(scale2i + 2*rij(3)*e2(3,i))
            temp2jzz7 =  -urr7*scale2i*rij(3)*rij(3)
            temp2jzz = temp2jzz5 + temp2jzz7
            temp2izz5 =  urr5*(scale2j + 2*rij(3)*e2(3,jbis))
            temp2izz7 =  -urr7*scale2j*rij(3)*rij(3)
            temp2izz = temp2izz5 + temp2izz7

            temp3jzz5 =  urr5*(scale3i + 2*rij(3)*e3(3,i))
            temp3jzz7 =  -urr7*scale3i*rij(3)*rij(3)
            temp3jzz = temp3jzz5 + temp3jzz7
            temp3izz5 =  urr5*(scale3j + 2*rij(3)*e3(3,jbis))
            temp3izz7 =  -urr7*scale3j*rij(3)*rij(3)
            temp3izz = temp3izz5 + temp3izz7

            temp4jzz5 =  urr5*(scale4i + 2*rij(3)*e4(3,i))
            temp4jzz7 =  -urr7*scale4i*rij(3)*rij(3)
            temp4jzz = temp4jzz5 + temp4jzz7
            temp4izz5 =  urr5*(scale4j + 2*rij(3)*e4(3,jbis))
            temp4izz7 =  -urr7*scale4j*rij(3)*rij(3)
            temp4izz = temp4izz5 + temp4izz7
c
c      alpha = 2, beta = 3
c
            temp1jyz5 =  urr5*(rij(3)*e1(2,i)+rij(2)*e1(3,i))
            temp1jyz7 =  -urr7*scale1i*rij(2)*rij(3)
            temp1jyz = temp1jyz5 + temp1jyz7
            temp1iyz5 =  urr5*(rij(3)*e1(2,jbis)+rij(2)*e1(3,jbis))
            temp1iyz7 =  -urr7*scale1j*rij(2)*rij(3)
            temp1iyz = temp1iyz5 + temp1iyz7

            temp2jyz5 =  urr5*(rij(3)*e2(2,i)+rij(2)*e2(3,i))
            temp2jyz7 =  -urr7*scale2i*rij(2)*rij(3)
            temp2jyz = temp2jyz5 + temp2jyz7
            temp2iyz5 =  urr5*(rij(3)*e2(2,jbis)+rij(2)*e2(3,jbis))
            temp2iyz7 =  -urr7*scale2j*rij(2)*rij(3)
            temp2iyz = temp2iyz5 + temp2iyz7

            temp3jyz5 =  urr5*(rij(3)*e3(2,i)+rij(2)*e3(3,i))
            temp3jyz7 =  -urr7*scale3i*rij(2)*rij(3)
            temp3jyz = temp3jyz5 + temp3jyz7
            temp3iyz5 =  urr5*(rij(3)*e3(2,jbis)+rij(2)*e3(3,jbis))
            temp3iyz7 =  -urr7*scale3j*rij(2)*rij(3)
            temp3iyz = temp3iyz5 + temp3iyz7

            temp4jyz5 =  urr5*(rij(3)*e4(2,i)+rij(2)*e4(3,i))
            temp4jyz7 =  -urr7*scale4i*rij(2)*rij(3)
            temp4jyz = temp4jyz5 + temp4jyz7
            temp4iyz5 =  urr5*(rij(3)*e4(2,jbis)+rij(2)*e4(3,jbis))
            temp4iyz7 =  -urr7*scale4j*rij(2)*rij(3)
            temp4iyz = temp4iyz5 + temp4iyz7
c    
c      beta = 1
c
            depx =  
     $       d1(1,jbis)*temp1jxx + d1(2,jbis)*temp1jxy 
     $      + d1(3,jbis)*temp1jxz
     $      + d1(1,i)*temp1ixx + d1(2,i)*temp1ixy 
     $      + d1(3,i)*temp1ixz
     $      + d2(1,jbis)*temp2jxx + d2(2,jbis)*temp2jxy 
     $      + d2(3,jbis)*temp2jxz
     $      + d2(1,i)*temp2ixx + d2(2,i)*temp2ixy 
     $      + d2(3,i)*temp2ixz
     $      + d3(1,jbis)*temp3jxx + d3(2,jbis)*temp3jxy 
     $      + d3(3,jbis)*temp3jxz
     $      + d3(1,i)*temp3ixx + d3(2,i)*temp3ixy 
     $      + d3(3,i)*temp3ixz
     $      + d4(1,jbis)*temp4jxx + d4(2,jbis)*temp4jxy 
     $      + d4(3,jbis)*temp4jxz
     $      + d4(1,i)*temp4ixx + d4(2,i)*temp4ixy 
     $      + d4(3,i)*temp4ixz
            adTbcor(1,i) = adTbcor(1,i) + depx 
            adTbcor(1,jbis) = adTbcor(1,jbis) - depx 
c    
c      beta = 2
c
            depy =  
     $       d1(1,jbis)*temp1jxy + d1(2,jbis)*temp1jyy 
     $      + d1(3,jbis)*temp1jyz
     $      + d1(1,i)*temp1ixy + d1(2,i)*temp1iyy 
     $      + d1(3,i)*temp1iyz
     $      + d2(1,jbis)*temp2jxy + d2(2,jbis)*temp2jyy 
     $      + d2(3,jbis)*temp2jyz
     $      + d2(1,i)*temp2ixy + d2(2,i)*temp2iyy 
     $      + d2(3,i)*temp2iyz
     $      + d3(1,jbis)*temp3jxy + d3(2,jbis)*temp3jyy 
     $      + d3(3,jbis)*temp3jyz
     $      + d3(1,i)*temp3ixy + d3(2,i)*temp3iyy 
     $      + d3(3,i)*temp3iyz
     $      + d4(1,jbis)*temp4jxy + d4(2,jbis)*temp4jyy 
     $      + d4(3,jbis)*temp4jyz
     $      + d4(1,i)*temp4ixy + d4(2,i)*temp4iyy 
     $      + d4(3,i)*temp4iyz
            adTbcor(2,i) = adTbcor(2,i) + depy 
            adTbcor(2,jbis) = adTbcor(2,jbis) - depy 
c    
c      beta = 3
c
            depz =  
     $       d1(1,jbis)*temp1jxz + d1(2,jbis)*temp1jyz 
     $      + d1(3,jbis)*temp1jzz
     $      + d1(1,i)*temp1ixz + d1(2,i)*temp1iyz 
     $      + d1(3,i)*temp1izz
     $      + d2(1,jbis)*temp2jxz + d2(2,jbis)*temp2jyz 
     $      + d2(3,jbis)*temp2jzz
     $      + d2(1,i)*temp2ixz + d2(2,i)*temp2iyz 
     $      + d2(3,i)*temp2izz
     $      + d3(1,jbis)*temp3jxz + d3(2,jbis)*temp3jyz 
     $      + d3(3,jbis)*temp3jzz
     $      + d3(1,i)*temp3ixz + d3(2,i)*temp3iyz 
     $      + d3(3,i)*temp3izz
     $      + d4(1,jbis)*temp4jxz + d4(2,jbis)*temp4jyz 
     $      + d4(3,jbis)*temp4jzz
     $      + d4(1,i)*temp4ixz + d4(2,i)*temp4iyz 
     $      + d4(3,i)*temp4izz
            adTbcor(3,i) = adTbcor(3,i) + depz 
            adTbcor(3,jbis) = adTbcor(3,jbis) - depz 

         end do !jj


         do j = 1, np11(iglob)
            dscale(ip11(j,iglob)) = 1.0d0
            uscale(ip11(j,iglob)) = 1.0d0
         end do
         do j = 1, np12(iglob)
            dscale(ip12(j,iglob)) = 1.0d0
            uscale(ip12(j,iglob)) = 1.0d0
         end do
         do j = 1, np13(iglob)
            dscale(ip13(j,iglob)) = 1.0d0
            uscale(ip13(j,iglob)) = 1.0d0
         end do
         do j = 1, np14(iglob)
            dscale(ip14(j,iglob)) = 1.0d0
            uscale(ip14(j,iglob)) = 1.0d0
         end do
         do j = 1, n12(iglob)
            pscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            pscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            pscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            pscale(i15(j,iglob)) = 1.0d0
         end do
      end do

      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      adme_self(:,1,1:npoleloc) = term*a1(:,1:npoleloc)
      adme_self(:,2,1:npoleloc) = term*a2(:,1:npoleloc)

      ade  = - ade_cor  ! THIS -1 EVERYWHERE IS WEIRD
      adme = - adme_cor + adme_self ! But hey... 
      adte = - adte_cor !  ...it works.
      adtb = - adtbcor!



      deallocate(dscale, pscale, uscale)
      return
      end
