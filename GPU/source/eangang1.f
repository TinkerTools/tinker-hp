c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangang1  --  angle-angle energy & derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangang1" calculates the angle-angle potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      subroutine eangang1
      use angang
      use angle
      use angpot
      use atmlst
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use group
      use math
      use usage
      use potent
      use tinheader
      use virial
      use mamd
      implicit none
      integer i,k,iangang,iangangloc
      integer ia,ib,ic,id,ie
      integer ialoc,ibloc,icloc,idloc,ieloc
      real(t_p) e,angle1
      real(t_p) dot,cosine
      real(t_p) dt1,deddt1
      real(t_p) dt2,deddt2
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xie,yie,zie
      real(t_p) xab,yab,zab
      real(t_p) xcb,ycb,zcb
      real(t_p) xdb,ydb,zdb
      real(t_p) xeb,yeb,zeb
      real(t_p) rab2,rcb2
      real(t_p) rdb2,reb2
      real(t_p) xp,yp,zp,rp
      real(t_p) xq,yq,zq,rq
      real(t_p) terma,termc
      real(t_p) termd,terme
      real(t_p) dedxia,dedyia,dedzia
      real(t_p) dedxib,dedyib,dedzib
      real(t_p) dedxic,dedyic,dedzic
      real(t_p) dedxid,dedyid,dedzid
      real(t_p) dedxie,dedyie,dedzie
      real(t_p) vxx,vyy,vzz
      real(t_p) vyx,vzx,vzy
      real(t_p) fgrp
      logical proceed
c
!$acc update host(deaa,vir,eaa)
c
c     zero out the angle-angle energy and first derivatives
c
      eaa = 0.0_ti_p
c
c     find the energy of each angle-angle interaction
c
      do iangangloc = 1, nangangloc
         iangang = angangglob(iangangloc)
         i = iaa(1,iangang)
         k = iaa(2,iangang)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(1,k)
         ie = iang(3,k)
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
         ieloc = loc(ie)
c
c     decide whether to compute the current interaction
c
         if (use_group)  call groups (fgrp,ia,ib,ic,id,ie,0)
         proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                               .or. use(id) .or. use(ie))
c
c     get the coordinates of the atoms in the angle
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
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
c
c     compute the values of the two bond angles
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            xeb = xie - xib
            yeb = yie - yib
            zeb = zie - zib
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
               call image (xdb,ydb,zdb)
               call image (xeb,yeb,zeb)
            end if
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            reb2 = xeb*xeb + yeb*yeb + zeb*zeb
            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            xq = yeb*zdb - zeb*ydb
            yq = zeb*xdb - xeb*zdb
            zq = xeb*ydb - yeb*xdb
            rp = sqrt(xp*xp + yp*yp + zp*zp)
            rq = sqrt(xq*xq + yq*yq + zq*zq)
            if (rp*rq .ne. 0.0_ti_p) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle1 = radian * acos(cosine)
               dt1 = angle1 - anat(i)
               dot = xdb*xeb + ydb*yeb + zdb*zeb
               cosine = dot / sqrt(rdb2*reb2)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle1 = radian * acos(cosine)
               dt2 = angle1 - anat(k)
c
c     get the energy and master chain rule terms for derivatives
c
               e = aaunit * kaa(iangang) * dt1 * dt2
               deddt1 = radian * e / dt1
               deddt2 = radian * e / dt2
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  deddt1 = deddt1 * fgrp
                  deddt2 = deddt2 * fgrp
               end if
c
c     find chain rule terms for the first bond angle deviation
c
               terma = -deddt1 / (rab2*rp)
               termc = deddt1 / (rcb2*rp)
               dedxia = terma * (yab*zp-zab*yp)
               dedyia = terma * (zab*xp-xab*zp)
               dedzia = terma * (xab*yp-yab*xp)
               dedxic = termc * (ycb*zp-zcb*yp)
               dedyic = termc * (zcb*xp-xcb*zp)
               dedzic = termc * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the second bond angle deviation
c
               termd = -deddt2 / (rdb2*rq)
               terme = deddt2 / (reb2*rq)
               dedxid = termd * (ydb*zq-zdb*yq)
               dedyid = termd * (zdb*xq-xdb*zq)
               dedzid = termd * (xdb*yq-ydb*xq)
               dedxie = terme * (yeb*zq-zeb*yq)
               dedyie = terme * (zeb*xq-xeb*zq)
               dedzie = terme * (xeb*yq-yeb*xq)
c
c     get the central atom derivative terms by difference
c
               dedxib = -dedxia - dedxic - dedxid - dedxie
               dedyib = -dedyia - dedyic - dedyid - dedyie
               dedzib = -dedzia - dedzic - dedzid - dedzie
c
c     increment the total angle-angle energy and derivatives
c
               eaa = eaa + e
               deaa(1,ibloc) = deaa(1,ibloc) + dedxib
               deaa(2,ibloc) = deaa(2,ibloc) + dedyib
               deaa(3,ibloc) = deaa(3,ibloc) + dedzib
c
               deaa(1,ialoc) = deaa(1,ialoc) + dedxia
               deaa(2,ialoc) = deaa(2,ialoc) + dedyia
               deaa(3,ialoc) = deaa(3,ialoc) + dedzia
c
               deaa(1,icloc) = deaa(1,icloc) + dedxic
               deaa(2,icloc) = deaa(2,icloc) + dedyic
               deaa(3,icloc) = deaa(3,icloc) + dedzic
c
               deaa(1,idloc) = deaa(1,idloc) + dedxid
               deaa(2,idloc) = deaa(2,idloc) + dedyid
               deaa(3,idloc) = deaa(3,idloc) + dedzid
c
               deaa(1,ieloc) = deaa(1,ieloc) + dedxie
               deaa(2,ieloc) = deaa(2,ieloc) + dedyie
               deaa(3,ieloc) = deaa(3,ieloc) + dedzie
c
c     aMD storage if waters are considered
c
               if (use_amd_wat1) then
               if (type(ialoc) == aMDwattype(1) .or. type(ibloc)
     $         == aMDwattype(1) .or. type(icloc) == aMDwattype(1)
     $         .or. type(idloc) == aMDwattype(1) .or. type(ieloc) ==
     $         aMDwattype(1)) then
                  eW1aMD = eW1aMD + e
                  deW1aMD(1,ialoc) = deW1aMD(1,ialoc) + dedxia
                  deW1aMD(2,ialoc) = deW1aMD(2,ialoc) + dedyia
                  deW1aMD(3,ialoc) = deW1aMD(3,ialoc) + dedzia
                  deW1aMD(1,ibloc) = deW1aMD(1,ibloc) + dedxib
                  deW1aMD(2,ibloc) = deW1aMD(2,ibloc) + dedyib
                  deW1aMD(3,ibloc) = deW1aMD(3,ibloc) + dedzib
                  deW1aMD(1,icloc) = deW1aMD(1,icloc) + dedxic
                  deW1aMD(2,icloc) = deW1aMD(2,icloc) + dedyic
                  deW1aMD(3,icloc) = deW1aMD(3,icloc) + dedzic
                  deW1aMD(1,idloc) = deW1aMD(1,idloc) + dedxid
                  deW1aMD(2,idloc) = deW1aMD(2,idloc) + dedyid
                  deW1aMD(3,idloc) = deW1aMD(3,idloc) + dedzid
                  deW1aMD(1,ieloc) = deW1aMD(1,ieloc) + dedxie
                  deW1aMD(2,ieloc) = deW1aMD(2,ieloc) + dedyie
                  deW1aMD(3,ieloc) = deW1aMD(3,ieloc) + dedzie
               end if
               end if
c
c     increment the internal virial tensor components
c
               vxx = xab*dedxia + xcb*dedxic + xdb*dedxid + xeb*dedxie
               vyx = yab*dedxia + ycb*dedxic + ydb*dedxid + yeb*dedxie
               vzx = zab*dedxia + zcb*dedxic + zdb*dedxid + zeb*dedxie
               vyy = yab*dedyia + ycb*dedyic + ydb*dedyid + yeb*dedyie
               vzy = zab*dedyia + zcb*dedyic + zdb*dedyid + zeb*dedyie
               vzz = zab*dedzia + zcb*dedzic + zdb*dedzid + zeb*dedzie
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
            end if
         end if
      end do
!$acc update device(deaa,vir,eaa)
      end
