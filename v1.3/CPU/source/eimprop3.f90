!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  subroutine eimprop3  --  imp. dihedral energy & analysis  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "eimprop3" calculates the improper dihedral potential
!     energy; also partitions the energy terms among the atoms
!
!
!> @brief 
!> calculates the improper dihedral potential
!> energy; also partitions the energy terms among the atoms
!> @param no params
subroutine eimprop3
   use action
   use analyz
   use atmlst
   use atmtyp
   use atoms
   use bound
   use domdec
   use energi
   use group
   use improp
   use inform
   use iounit
   use math
   use torpot
   use usage
   implicit none
   integer i,ia,ib,ic,id,ibloc,icloc
   integer iimprop
   real*8 e,dt
   real*8 ideal,force
   real*8 cosine,sine
   real*8 rcb,angle
   real*8 xt,yt,zt
   real*8 xu,yu,zu
   real*8 xtu,ytu,ztu
   real*8 rt2,ru2,rtru
   real*8 xia,yia,zia
   real*8 xib,yib,zib
   real*8 xic,yic,zic
   real*8 xid,yid,zid
   real*8 xba,yba,zba
   real*8 xcb,ycb,zcb
   real*8 xdc,ydc,zdc
   real*8 fgrp
   logical proceed
   logical header,huge
!
   if (deb_Path) write(iout,*), 'eimprop3 '
!
!
!     zero out improper dihedral energy and partitioning terms
!
   neid = 0
   eid = 0.0d0
   aeid = 0.0d0
   header = .true.
!
!     calculate the improper dihedral angle energy term
!
   do iimprop = 1, niproploc
      i = impropglob(iimprop)
      ia = iiprop(1,i)
      ib = iiprop(2,i)
      ibloc = loc(ib)
      ic = iiprop(3,i)
      icloc = loc(ic)
      id = iiprop(4,i)
!
!     decide whether to compute the current interaction
!
      if (use_group)  call groups (fgrp,ia,ib,ic,id,0,0)
      proceed = (use(ia) .or. use(ib) .or.&
      &use(ic) .or. use(id))
!
!     compute the value of the improper dihedral angle
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
         ru2 = xu*xu + yu*yu + zu*zu
         rtru = sqrt(rt2 * ru2)
         if (rtru .ne. 0.0d0) then
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            cosine = (xt*xu + yt*yu + zt*zu) / rtru
            sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
            cosine = min(1.0d0,max(-1.0d0,cosine))
            angle = radian * acos(cosine)
            if (sine .lt. 0.0d0)  angle = -angle
!
!     set the improper dihedral parameters for this angle
!
            ideal = vprop(i)
            force = kprop(i)
            if (abs(angle+ideal) .lt. abs(angle-ideal))&
            &ideal = -ideal
            dt = angle - ideal
            do while (dt .gt. 180.0d0)
               dt = dt - 360.0d0
            end do
            do while (dt .lt. -180.0d0)
               dt = dt + 360.0d0
            end do
!
!     calculate the improper dihedral energy
!
            e = idihunit * force * dt**2
!
!     scale the interaction based on its group membership
!
            if (use_group)  e = e * fgrp
!
!     increment the total improper dihedral energy
!
            neid = neid + 1
            eid = eid + e
            aeid(ibloc) = aeid(ibloc) + 0.5d0*e
            aeid(icloc) = aeid(icloc) + 0.5d0*e
!
!     print a message if the energy of this interaction is large
!
            huge = (e .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
10                format (/,' Individual Improper Dihedral',&
                  &' Interactions :',&
                  &//,' Type',25x,'Atom Names',21x,'Angle',&
                  &6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ic,&
               &name(ic),id,name(id),angle,e
20             format (' Improper',2x,4(i7,'-',a3),f11.4,f12.4)
            end if
         end if
      end if
   end do
   return
end
