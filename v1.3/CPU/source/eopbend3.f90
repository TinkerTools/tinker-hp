!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine eopbend3  --  out-of-plane bending & analysis  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "eopbend3" computes the out-of-plane bend potential energy at
!     trigonal centers via a Wilson-Decius-Cross or Allinger angle;
!     also partitions the energy among the atoms
!
!
!> @brief 
!> computes the out-of-plane bend potential energy at
!> trigonal centers via a Wilson-Decius-Cross or Allinger angle;
!> also partitions the energy among the atoms
!> @param no params
subroutine eopbend3
   use action
   use analyz
   use angle
   use angpot
   use atmlst
   use atmtyp
   use atoms
   use bound
   use domdec
   use energi
   use group
   use inform
   use iounit
   use math
   use opbend
   use usage
   implicit none
   integer i,iopbend,iopbendloc
   integer ia,ib,ic,id,ibloc
   real*8 e,angle1,force
   real*8 cosine
   real*8 dt,dt2,dt3,dt4
   real*8 xia,yia,zia
   real*8 xib,yib,zib
   real*8 xic,yic,zic
   real*8 xid,yid,zid
   real*8 xab,yab,zab
   real*8 xcb,ycb,zcb
   real*8 xdb,ydb,zdb
   real*8 xad,yad,zad
   real*8 xcd,ycd,zcd
   real*8 rdb2,rad2,rcd2
   real*8 rab2,rcb2
   real*8 cc,ee,bkk2
   real*8 fgrp
   logical proceed
   logical header,huge
!
   if (deb_Path) write(iout,*), 'eopbend3 '
!
!
!
!     zero out the out-of-plane bend energy and partitioning
!
   neopb = 0
   eopb = 0.0d0
   aeopb = 0.0d0
   header = .true.
!
!     calculate the out-of-plane bending energy term
!
   do iopbendloc = 1, nopbendloc
      iopbend = opbendglob(iopbendloc)
      i = iopb(iopbend)
      ia = iang(1,i)
      ib = iang(2,i)
      ibloc = loc(ib)
      ic = iang(3,i)
      id = iang(4,i)
      force = opbk(iopbend)
!
!     decide whether to compute the current interaction
!
      if (use_group)  call groups (fgrp,ia,ib,ic,id,0,0)
      proceed = (use(ia) .or. use(ib) .or.&
      &use(ic) .or. use(id))
!
!     get the coordinates of the atoms at trigonal center
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
!
!     compute the out-of-plane bending angle
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
         xad = xia - xid
         yad = yia - yid
         zad = zia - zid
         xcd = xic - xid
         ycd = yic - yid
         zcd = zic - zid
         if (use_polymer) then
            call image (xab,yab,zab)
            call image (xcb,ycb,zcb)
            call image (xdb,ydb,zdb)
            call image (xad,yad,zad)
            call image (xcd,ycd,zcd)
         end if
!
!     W-D-C angle between A-B-C plane and B-D vector for D-B<AC
!
         if (opbtyp .eq. 'W-D-C') then
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            cc = rab2*rcb2 - (xab*xcb+yab*ycb+zab*zcb)**2
!
!     Allinger angle between A-C-D plane and D-B vector for D-B<AC
!
         else if (opbtyp .eq. 'ALLINGER') then
            rad2 = xad*xad + yad*yad + zad*zad
            rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
            cc = rad2*rcd2 - (xad*xcd+yad*ycd+zad*zcd)**2
         end if
!
!     find the out-of-plane angle bending energy
!
         ee = xdb*(yab*zcb-zab*ycb) + ydb*(zab*xcb-xab*zcb)&
         &+ zdb*(xab*ycb-yab*xcb)
         rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
         if (rdb2.ne.0.0d0 .and. cc.ne.0.0d0) then
            bkk2 = rdb2 - ee*ee/cc
            cosine = sqrt(bkk2/rdb2)
            cosine = min(1.0d0,max(-1.0d0,cosine))
            angle1 = radian * acos(cosine)
            dt = angle1
            dt2 = dt * dt
            dt3 = dt2 * dt
            dt4 = dt2 * dt2
            e = opbunit * force * dt2&
            &* (1.0d0+copb*dt+qopb*dt2+popb*dt3+sopb*dt4)
!
!     scale the interaction based on its group membership
!
            if (use_group)  e = e * fgrp
!
!     increment the total out-of-plane bending energy
!
            neopb = neopb + 1
            eopb = eopb + e
            aeopb(ibloc) = aeopb(ibloc) + e
!
!     print a message if the energy of this interaction is large
!
            huge = (e .gt. 2.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
10                format (/,' Individual Out-of-Plane Bend',&
                  &' Interactions :',&
                  &//,' Type',25x,'Atom Names',21x,'Angle',&
                  &6x,'Energy',/)
               end if
               write (iout,20)  id,name(id),ib,name(ib),ia,&
               &name(ia),ic,name(ic),angle1,e
20             format (' O-P-Bend',2x,4(i7,'-',a3),f11.4,f12.4)
            end if
         end if
      end if
   end do
   return
end
