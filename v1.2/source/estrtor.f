c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine estrtor  --  stretch-torsion cross term energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "estrtor" calculates the stretch-torsion potential energy
c
c
      subroutine estrtor
      use atmlst
      use atoms
      use bond
      use bound
      use energi
      use group
      use inform
      use iounit
      use strtor
      use torpot
      use tors
      use usage
      implicit none
      integer i,k,istrtor,iistrtor
      integer ia,ib,ic,id
      real*8 e,dr
      real*8 rt2,ru2,rtru
      real*8 rba,rcb,rdc
      real*8 e1,e2,e3
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 sine,cosine
      real*8 sine2,cosine2
      real*8 sine3,cosine3
      real*8 phi1,phi2,phi3
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 fgrp
      logical proceed
c
      if (deb_Path) write(iout,*), 'estrtor '
c
c
c     zero out the stretch-torsion energy
c
      ebt = 0.0d0
c
c     calculate the stretch-torsion interaction energy term
c
      do istrtor = 1, nstrtorloc
         iistrtor = strtorglob(istrtor)
         i = ist(1,iistrtor)
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         if (use_group)  call groups (fgrp,ia,ib,ic,id,0,0)
         proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
            end if
            rba = sqrt(xba*xba + yba*yba + zba*zba)
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
            if (min(rba,rcb,rdc) .ne. 0.0d0) then
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
              rt2 = xt*xt + yt*yt + zt*zt
              rt2 = max(rt2,0.000001d0)
              ru2 = xu*xu + yu*yu + zu*zu
              ru2 = max(ru2,0.000001d0)
              rtru = sqrt(rt2 * ru2)
              rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
              cosine = (xt*xu + yt*yu + zt*zu) / rtru
              sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     compute multiple angle trigonometry and phase terms
c
              c1 = tors1(3,i)
              s1 = tors1(4,i)
              c2 = tors2(3,i)
              s2 = tors2(4,i)
              c3 = tors3(3,i)
              s3 = tors3(4,i)
              cosine2 = cosine*cosine - sine*sine
              sine2 = 2.0d0 * cosine * sine
              cosine3 = cosine*cosine2 - sine*sine2
              sine3 = cosine*sine2 + sine*cosine2
              phi1 = 1.0d0 + (cosine*c1 + sine*s1)
              phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
              phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
c
c     get the stretch-torsion values for the first bond
c
              v1 = kst(1,iistrtor)
              v2 = kst(2,iistrtor)
              v3 = kst(3,iistrtor)
              k = ist(2,iistrtor)
              dr = rba - bl(k)
              e1 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     get the stretch-torsion values for the second bond
c
              v1 = kst(4,iistrtor)
              v2 = kst(5,iistrtor)
              v3 = kst(6,iistrtor)
              k = ist(3,iistrtor)
              dr = rcb - bl(k)
              e2 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     get the stretch-torsion values for the third bond
c
              v1 = kst(7,iistrtor)
              v2 = kst(8,iistrtor)
              v3 = kst(9,iistrtor)
              k = ist(4,iistrtor)
              dr = rdc - bl(k)
              e3 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
              if (use_group) then
                 e1 = e1 * fgrp
                 e2 = e2 * fgrp
                 e3 = e3 * fgrp
              end if
c
c     increment the total stretch-torsion energy
c
              e = e1 + e2 + e3
              ebt = ebt + e
            end if
         end if
      end do
      return
      end
