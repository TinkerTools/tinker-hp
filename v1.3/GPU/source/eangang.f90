!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine eangang  --  angle-angle energy and analysis   ##
!     ##                                                            ##
!     ################################################################
!
!
!     "eangang" calculates the angle-angle potential energy;
!     also partitions the energy among the atoms
!
!
#include "tinker_macro.h"
subroutine eangang
   use angang
   use angle
   use angpot
   use atmlst
   use atmtyp
   use atoms
   use bound
   use energi
   use group
   use inform
   use iounit
   use math
   use tinheader
   use usage
   implicit none
   integer i,k,iangang,iangangloc
   integer ia,ib,ic,id,ie
   real(t_p) e,dt1,dt2
   real(t_p) angle1,dot,cosine
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
   real(t_p) fgrp
   logical proceed
   logical header,huge
!
!
!     zero out the angle-angle cross term energy
!
   eaa = 0.0_ti_p
!
!     find the energy of each angle-angle interaction
!
   do iangangloc = 1, nangangloc
      iangang = angangglob(iangangloc)
      i = iaa(1,iangang)
      k = iaa(2,iangang)
      ia = iang(1,i)
      ib = iang(2,i)
      ic = iang(3,i)
      id = iang(1,k)
      ie = iang(3,k)
!
!     decide whether to compute the current interaction
!
      if (use_group)  call groups (fgrp,ia,ib,ic,id,ie,0)
      proceed = (use(ia) .or. use(ib) .or. use(ic)&
         &.or. use(id) .or. use(ie))
!
!     get the coordinates of the atoms in the angle
!
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
!
!     compute the values of the two bond angles
!
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
         if (rab2*rcb2*rdb2*reb2 .ne. 0.0_ti_p) then
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
!
!     get the angle-angle interaction energy
!
            e = aaunit * kaa(iangang) * dt1 * dt2
!
!     scale the interaction based on its group membership
!
            if (use_group) e = e * fgrp
!
!     increment the total angle-angle energy
!
            eaa = eaa + e
!
!     print a message if the energy of this interaction is large
!
            huge = (e .gt. 5.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
10                format (/,' Individual Angle-Angle Interactions :',&
                     &//,' Type',10x,'Center',6x,'Angle1',&
                     &6x,'Angle2',4x,'dAngle1',&
                     &3x,'dAngle2',6x,'Energy',/)
               end if
               write (iout,20)  ib,name(ib),ia,ic,id,ie,dt1,dt2,e
20             format (' AngAng',4x,i7,'-',a3,4i6,2f10.4,f12.4)
            end if
         end if
      end if
   end do
end
