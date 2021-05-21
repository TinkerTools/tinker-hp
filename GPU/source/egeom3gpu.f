c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egeom3  --  restraint energy terms & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egeom3" calculates the energy due to restraints on positions,
c     distances, angles and torsions as well as Gaussian basin and
c     droplet restraints; also partitions energy among the atoms
c
c
#include "tinker_precision.h"
      module egeom3gpu_inl
        contains
#include "image.f.inc"
#include "midpointimage.f.inc"
      end module

      subroutine egeom3gpu
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use domdec
      use energi
      use egeom3gpu_inl
      use group
      use inform
      use inter
      use iounit
      use kgeoms
      use math
      use molcul
      use tinheader ,only:ti_p,re_p
      use usage
      implicit none
      integer i,j,k
      integer iglob,kglob
      real(t_p) xk,yk,zk  
      integer ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer inpfix,indfix,inafix,intfix,ingfix,inchir
      logical docompute
      real(t_p) e,dt,dt2
      real(t_p) xr,yr,zr
      real(t_p) r,r2,r6,r12
      real(t_p) angle,target
      real(t_p) dot,force
      real(t_p) cosine,sine
      real(t_p) rab2,rcb2
      real(t_p) xt,yt,zt
      real(t_p) xu,yu,zu
      real(t_p) xtu,ytu,ztu
      real(t_p) rt2,ru2
      real(t_p) rtru,rcb
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xab,yab,zab
      real(t_p) xba,yba,zba
      real(t_p) xcb,ycb,zcb
      real(t_p) xdc,ydc,zdc
      real(t_p) xad,yad,zad
      real(t_p) xbd,ybd,zbd
      real(t_p) xcd,ycd,zcd
      real(t_p) df1,df2
      real(t_p) af1,af2
      real(t_p) tf1,tf2,t1,t2
      real(t_p) gf1,gf2
      real(t_p) weigh,size
      real(t_p) weigha,weighb
      real(t_p) xcm,ycm,zcm
      real(t_p) cf1,cf2,vol
      real(t_p) c1,c2,c3
      real(t_p) xi,yi,zi,ri
      real(t_p) a,b,buffer,term
      real(t_p) time0,time1
      logical proceed,intermol
      logical header,huge
      real(t_p) eg1
      integer neg1

      if (deb_Path) write(*,*) 'egeom3gpu'
c
c     zero out the restraint energy and partitioning terms
c
      neg = 0
      eg  = 0.0_re_p
c     aeg = 0.0_ti_p

!$acc enter data copyin(neg)
!$acc data present(mass,x,y,z,use)
!$acc&     present(neg, eg, einter)
c
c     compute the energy for position restraint terms
c
      if (npfix.ne.0) then

      header = .true.

!$acc parallel loop async
!$acc&         default(present)
      do inpfix = 1, npfixloc
         i = npfixglob(inpfix)
         ia = ipfix(i)
         ialoc = loc(ia)
         proceed = .true.
         if (proceed)  proceed = (use(ia))
         if (proceed) then
            xr = 0.0_ti_p
            yr = 0.0_ti_p
            zr = 0.0_ti_p
            if (kpfix(1,i) .ne. 0)  xr = x(ia) - xpfix(i)
            if (kpfix(2,i) .ne. 0)  yr = y(ia) - ypfix(i)
            if (kpfix(3,i) .ne. 0)  zr = z(ia) - zpfix(i)
            if (use_bounds)  call image_inl (xr,yr,zr)
            r = sqrt(xr*xr + yr*yr + zr*zr)
            force = pfix(1,i)
            dt = max(0.0_ti_p,r-pfix(2,i))
            dt2 = dt * dt
            e = force * dt2
            neg = neg + 1
            eg = eg+ e
!$acc atomic
            aeg(ialoc) = aeg(ialoc) + e
#ifndef _OPENACC
            huge = (e .gt. 10.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Atomic Position Restraint',
     &                       ' Terms :',
     &                    //,' Type',9x,'Atom Name',13x,'Target',
     &                       ' Position',7x,'Distance',6x,'Energy',/)
               end if
               dt = sqrt(dt2)
               if (kpfix(2,i).eq.0 .and. kpfix(3,i).eq.0) then
                  write (iout,20)  ia,name(ia),xpfix(i),dt,e
   20             format (' Position',2x,i7,'-',a3,4x,f10.4,
     &                       5x,'----',6x,'----',1x,f10.4,f12.4)
               else if (kpfix(1,i).eq.0 .and. kpfix(3,i).eq.0) then
                  write (iout,30)  ia,name(ia),ypfix(i),dt,e
   30             format (' Position',2x,i7,'-',a3,9x,'----',1x,
     &                       f10.4,5x,'----',1x,f10.4,f12.4)
               else if (kpfix(1,i).eq.0 .and. kpfix(2,i).eq.0) then
                  write (iout,40)  ia,name(ia),zpfix(i),dt,e
   40             format (' Position',2x,i7,'-',a3,9x,'----',
     &                       6x,'----',1x,2f10.4,f12.4)
               else
                 write (iout,50)  ia,name(ia),xpfix(i),ypfix(i),
     &                             zpfix(i),dt,e
   50             format (' Position',2x,i7,'-',a3,4x,4f10.4,f12.4)
               end if
            end if
#endif
         end if
      end do

      end if
      if (ndfix.ne.0) then
c
c     compute the energy for distance restraint terms
c
      header = .true.

!$acc parallel loop async
!$acc&         default(present)
      do indfix = 1, ndfixloc
         i = ndfixglob(indfix)
         ia = idfix(1,i)
         ialoc = loc(ia)
         ib = idfix(2,i)
         ibloc = loc(ib)
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib))
         if (proceed) then
            xr = x(ia) - x(ib)
            yr = y(ia) - y(ib)
            zr = z(ia) - z(ib)
            intermol = (molcule(ia) .ne. molcule(ib))
            if (use_bounds .and. intermol) 
     &         call image_inl(xr,yr,zr)
c            if (use_bounds)  call image (xr,yr,zr)
            r = sqrt(xr*xr + yr*yr + zr*zr)
            force = dfix(1,i)
            df1 = dfix(2,i)
            df2 = dfix(3,i)
            target = r
            if (r .lt. df1)  target = df1
            if (r .gt. df2)  target = df2
            dt = r - target
            dt2 = dt * dt
            e = force * dt2
            neg = neg + 1
            eg = eg + e
!$acc atomic update 
            aeg(ialoc) = aeg(ialoc) + 0.5_ti_p*e
!$acc atomic update 
            aeg(ibloc) = aeg(ibloc) + 0.5_ti_p*e
c           if (intermol) then
c              einter = einter + e
c           end if
#ifndef _OPENACC
            huge = (e .gt. 10.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,60)
   60             format (/,' Individual Interatomic Distance',
     &                       ' Restraint Terms :',
     &                    //,' Type',14x,'Atom Names',16x,'Ideal Range',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,70)  ia,name(ia),ib,name(ib),df1,df2,r,e
   70          format (' Distance',2x,2(i7,'-',a3),
     &                    7x,2f8.2,f10.4,f12.4)
            end if
#endif
         end if
      end do

      end if
      if (nafix.ne.0) then
c
c     compute the energy for angle restraint terms
c
!$acc parallel loop async
!$acc&         default(present)
      do inafix = 1, nafixloc
         i = nafixglob(inafix)
         ia = iafix(1,i)
         ialoc = loc(ia)
         ib = iafix(2,i)
         ibloc = loc(ib)
         ic = iafix(3,i)
         icloc = loc(ic)
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
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
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            if (rab2.ne.0.0_ti_p .and. rcb2.ne.0.0_ti_p) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle = radian * acos(cosine)
               force = afix(1,i)
               af1 = afix(2,i)
               af2 = afix(3,i)
               target = angle
               if (angle .lt. af1)  target = af1
               if (angle .gt. af2)  target = af2
               dt  = angle - target
               dt = dt / radian
               dt2 = dt * dt
               e   = force * dt2
               neg = neg + 1
               eg  = eg + e
!$acc atomic
               aeg(ibloc) = aeg(ibloc) + e
#ifndef _OPENACC
               huge = (e .gt. 10.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,80)
   80                format (/,' Individual Interatomic Angle',
     &                          ' Restraint Terms :',
     &                       //,' Type',14x,'Atom Numbers',14x,'Ideal',
     &                          ' Range',4x,'Actual',6x,'Energy',/)
                  end if
                  write (iout,90)  ia,ib,ic,af1,af2,angle,e
   90             format (' Angle',8x,3i6,8x,2f8.2,f10.4,f12.4)
               end if
#endif
            end if
         end if
      end do

      end if
      if (ntfix.ne.0) then
c
c     compute the energy for torsional restraint terms
c
      header = .true.

!$acc parallel loop gang vector async
!$acc&         default(present) reduction(+:eg,neg)
      do intfix = 1, ntfixloc
         i = ntfixglob(intfix)
         ia = itfix(1,i)
         ib = itfix(2,i)
         ibloc = loc(ib)
         ic = itfix(3,i)
         icloc = loc(ic)
         id = itfix(4,i)
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                            use(ic) .or. use(id))
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
            if (rtru .ne. 0.0_ti_p) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0_ti_p)  angle = -angle
               force = tfix(1,i)
               tf1 = tfix(2,i)
               tf2 = tfix(3,i)
               if (angle.gt.tf1 .and. angle.lt.tf2) then
                  target = angle
               else if (angle.gt.tf1 .and. tf1.gt.tf2) then
                  target = angle
               else if (angle.lt.tf2 .and. tf1.gt.tf2) then
                  target = angle
               else
                  t1 = angle - tf1
                  t2 = angle - tf2
                  if (t1 .gt. 180.0_ti_p) then
                     t1 = t1 - 360.0_ti_p
                  else if (t1 .lt. -180.0_ti_p) then
                     t1 = t1 + 360.0_ti_p
                  end if
                  if (t2 .gt. 180.0_ti_p) then
                     t2 = t2 - 360.0_ti_p
                  else if (t2 .lt. -180.0_ti_p) then
                     t2 = t2 + 360.0_ti_p
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     target = tf1
                  else
                     target = tf2
                  end if
               end if
               dt = angle - target
               if (dt .gt. 180.0_ti_p) then
                  dt = dt - 360.0_ti_p
               else if (dt .lt. -180.0_ti_p) then
                  dt = dt + 360.0_ti_p
               end if
               dt = dt / radian
               dt2 = dt * dt
               e = force * dt2
               neg = neg + 1
               eg = eg + e
!$acc atomic update        
               aeg(ibloc) = aeg(ibloc) + 0.5_ti_p*e
!$acc atomic update        
               aeg(icloc) = aeg(icloc) + 0.5_ti_p*e
c              if (molcule(ia).ne.molcule(ib) .or.
c    &             molcule(ia).ne.molcule(ic) .or.
c    &             molcule(ia).ne.molcule(id)) then
c                 einter = einter + e
c              end if
#ifndef _OPENACC
               huge = (e .gt. 10.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,100)
  100                format (/,' Individual Torsional Angle Restraint
     &',                          ' Terms :',
     &                       //,' Type',14x,'Atom Numbers',14x,'Ideal',
     &                          ' Range',4x,'Actual',6x,'Energy',/)
                  end if
                  write (iout,110)  ia,ib,ic,id,tf1,tf2,angle,e
  110             format (' Torsion',4x,4i6,4x,2f8.2,f10.4,f12.4)
               end if
#endif
            end if
         end if
      end do

      end if
      if (ngfix.ne.0) then
c
c     compute the energy for group distance restraint terms
c
      header = .true.

!$acc parallel loop gang vector_length(32) 
!$acc&         default(present) async
!$acc&         private(xcm,ycm,zcm)
      do ingfix = 1, ngfixloc
         i = ngfixglob(ingfix)
         ia = igfix(1,i)
         ib = igfix(2,i)
         xcm = 0.0_ti_p
         ycm = 0.0_ti_p
         zcm = 0.0_ti_p
!$acc loop vector
         do j = igrp(1,ia), igrp(2,ia)
            k = kgrp(j)
            weigh = mass(k)
            xcm = xcm + x(k)*weigh
            ycm = ycm + y(k)*weigh
            zcm = zcm + z(k)*weigh
         end do
         weigha = max(1.0_ti_p,grpmass(ia))
         xr = xcm / weigha
         yr = ycm / weigha
         zr = zcm / weigha
         xcm = 0.0_ti_p
         ycm = 0.0_ti_p
         zcm = 0.0_ti_p
!$acc loop vector
         do j = igrp(1,ib), igrp(2,ib)
            k = kgrp(j)
            weigh = mass(k)
            xcm = xcm + x(k)*weigh
            ycm = ycm + y(k)*weigh
            zcm = zcm + z(k)*weigh
         end do
         weighb = max(1.0_ti_p,grpmass(ib))
         xr = xr - xcm/weighb
         yr = yr - ycm/weighb
         zr = zr - zcm/weighb
         intermol = (molcule(kgrp(igrp(1,ia))) .ne.
     &               molcule(kgrp(igrp(1,ib))))
         if (use_bounds .and. intermol) call image_inl(xr,yr,zr)
c         if (use_bounds)  call image (xr,yr,zr)
         r = sqrt(xr*xr + yr*yr + zr*zr)
         force = gfix(1,i)
         gf1 = gfix(2,i)
         gf2 = gfix(3,i)
         target = r
         if (r .lt. gf1)  target = gf1
         if (r .gt. gf2)  target = gf2
         dt  = r - target
         dt2 = dt * dt
         e   = force * dt2
!$acc loop vector
         do j = 1,32
            if (j.eq.1) then
               neg = neg + 1
                eg =  eg + e
            end if
         end do
#ifndef _OPENACC
         size = real(igrp(2,ia - igrp(1,ia) + 1),t_p)
!$acc loop vector
         do j = igrp(1,ia), igrp(2,ia)
            k = kgrp(j)
!$acc atomic update 
            aeg(k) = aeg(k) + 0.5_ti_p*e/size
         end do
         size = real(igrp(2,ib - igrp(1,ib) + 1),t_p)
!$acc loop vector
         do j = igrp(1,ib), igrp(2,ib)
            k = kgrp(j)
!$acc atomic update 
            aeg(k) = aeg(k) + 0.5_ti_p*e/size
         end do
c        if (intermol) then
c           einter = einter + e
c        end if
         huge = (e .gt. 10.0_ti_p)
         if (debug .or. (verbose.and.huge)) then
            if (header) then
               header = .false.
               write (iout,120)
  120          format (/,' Individual Intergroup Distance',
     &                    ' Restraint Terms :',
     &                 //,' Type',13x,'Group Numbers',14x,'Ideal Range',
     &                    4x,'Actual',6x,'Energy',/)
            end if
            write (iout,130)  ia,ib,gf1,gf2,r,e
  130       format (' Distance',7x,2i7,10x,2f8.2,f10.4,f12.4)
         end if
#endif
      end do

      end if
      if (nchir.ne.0) then
c
c     compute the energy for chirality restraint terms
c
      header = .true.

!$acc parallel loop async gang vector
!$acc&         default(present) reduction(+:eg,neg)
      do inchir = 1, nchirloc
         i  = nchirglob(inchir)
         ia = ichir(1,i)
         ialoc = loc(ia)
         ib = ichir(2,i)
         ibloc = loc(ib)
         ic = ichir(3,i)
         icloc = loc(ic)
         id = ichir(4,i)
         idloc = loc(id)
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                            use(ic) .or. use(id))
         if (proceed) then
            xad = x(ia) - x(id)
            yad = y(ia) - y(id)
            zad = z(ia) - z(id)
            xbd = x(ib) - x(id)
            ybd = y(ib) - y(id)
            zbd = z(ib) - z(id)
            xcd = x(ic) - x(id)
            ycd = y(ic) - y(id)
            zcd = z(ic) - z(id)
            c1 = ybd*zcd - zbd*ycd
            c2 = ycd*zad - zcd*yad
            c3 = yad*zbd - zad*ybd
            vol   = xad*c1 + xbd*c2 + xcd*c3
            force = chir(1,i)
            cf1 = chir(2,i)
            cf2 = chir(3,i)
            target = vol
            if (vol .lt. min(cf1,cf2))  target = min(cf1,cf2)
            if (vol .gt. max(cf1,cf2))  target = max(cf1,cf2)
            dt  = vol - target
            dt2 = dt * dt
            e   = force * dt2
            eg  = eg + e
            neg = neg + 1
!$acc atomic update  
            aeg(ialoc) = aeg(ialoc) + 0.25_ti_p*e
!$acc atomic update  
            aeg(ibloc) = aeg(ibloc) + 0.25_ti_p*e
!$acc atomic update  
            aeg(icloc) = aeg(icloc) + 0.25_ti_p*e
!$acc atomic update  
            aeg(idloc) = aeg(idloc) + 0.25_ti_p*e
#ifndef _OPENACC
            huge = (e .gt. 10.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,140)
  140             format (/,' Individual Chirality Restraint Terms :'
     &,                    //,' Type',14x,'Atom Numbers',14x,'Ideal',
     &                       ' Range',4x,'Actual',6x,'Energy',/)
               end if
               write (iout,150)  ia,ib,ic,id,cf1,cf2,vol,e
  150          format (' Chiral',5x,4i6,4x,2f8.2,f10.4,f12.4)
            end if
#endif
         end if
      end do

      end if
c
c     compute the energy for a Gaussian basin restraint
c
      if (use_basin) then

         header = .true.
!$acc parallel loop default(present) async
         do i = 1, nbloc
            iglob = glob(i)
            xi = x(iglob)
            yi = y(iglob)
            zi = z(iglob)
!$acc loop vector
            do k = 1, nbloc
               kglob = glob(k)
               xk = x(kglob)
               yk = y(kglob)
               zk = z(kglob)
               proceed = .true.
               if (proceed)  proceed = (use(iglob) .or. use(kglob))
               if (proceed) then
                  if (kglob.le.iglob) cycle

                  xr = xi - xk
                  yr = yi - yk
                  zr = zi - zk
                  call midpoint_inl(xk,yk,zk,xr,yr,zr,docompute)
                  if (.not.(docompute)) cycle
c                  r2 = xr*xr + yr*yr + zr*zr
c                  term = -width * r2
c                  term = -width *(xr*xr + yr*yr + zr*zr)
                  term=-width *(xr**2+yr**2+zr**2)
                  e = 0.0_ti_p
                  if (term .gt. -50.0_ti_p)  e = depth * exp(term)
                  e = e - depth
                  neg = neg + 1
                  eg = eg + e
!$acc atomic update    
                  aeg(i) = aeg(i) + 0.5_ti_p*e
!$acc atomic update   
                  aeg(k) = aeg(k) + 0.5_ti_p*e
#ifndef _OPENACC
                  huge = (e .gt. 10.0_ti_p)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,160)
  160                   format (/,' Individual Gaussian Basin',
     &                             ' Restraint Terms :',
     &                          //,' Type',14x,'Atom Names',22x,'Ideal',
     &                             4x,'Actual',6x,'Energy',/)
                     end if
                     r = sqrt(r2)
                     write (iout,170)  i,name(iglob),k,name(kglob),
     $                 0.0_ti_p,r,e
  170                format (' Distance',2x,2(i7,'-',a3),
     &                          13x,2f10.4,f12.4)
                  end if
#endif
               end if
            end do
         end do

      end if
c
c     compute the energy for a spherical droplet restraint
c
      if (use_wall) then

         header = .true.
         buffer = 2.5_ti_p
         a = 2048.0_ti_p
         b = 64.0_ti_p
!$acc parallel loop default(present) async
         do i = 1, nloc
            iglob = glob(i)
            proceed = .true.
            if (proceed)  proceed = (use(iglob))
            if (proceed) then
               xi = x(iglob)
               yi = y(iglob)
               zi = z(iglob)
               ri = sqrt(xi**2 + yi**2 + zi**2)
               r = rwall + buffer - ri
               r2 = r * r
               r6 = r2 * r2 * r2
               r12 = r6 * r6
               e = a/r12 - b/r6
               neg = neg + 1
               eg = eg + e
               aeg(i) = aeg(i) + e
#ifndef _OPENACC
               huge = (e .gt. 10.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,180)
  180                format (/,' Individual Spherical Boundary',
     &                          ' Restraint Terms :',
     &                       //,' Type',14x,'Atom Name',30x,'Distance',
     &                          6x,'Energy',/)
                  end if
                  write (iout,190)  iglob,name(iglob),ri,e
  190             format (' Wall',11x,i7,'-',a3,29x,f10.4,f12.4)
               end if
#endif
            end if
         end do

      end if

!$acc update self(aeg) async
!$acc end data
!$acc exit data copyout(neg) async
!$acc wait
      end
