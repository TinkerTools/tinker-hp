c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eopdist3  --  out-of-plane distance & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eopdist3" computes the out-of-plane distance potential energy
c     at trigonal centers via the central atom height; also partitions
c     the energy among the atoms
c
c
#include "tinker_macro.h"
      subroutine eopdist3
      use action
      use analyz
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
      use opdist
      use tinheader ,only:ti_p,re_p
      use usage
      implicit none
      integer i,ia,ib,ic,id,iopdist
      integer ialoc
      real(t_p) e,force
      real(t_p) dt,dt2,dt3,dt4
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xid,yid,zid
      real(t_p) xad,yad,zad
      real(t_p) xbd,ybd,zbd
      real(t_p) xcd,ycd,zcd
      real(t_p) xt,yt,zt,rt2
      real(t_p) fgrp
      logical proceed
      logical header,huge
c
c
c     zero out the out-of-plane distance energy and partitioning
c
      neopd = 0
      eopd = 0.0_ti_p
      aeopd = 0.0_ti_p
      header = (rank.eq.0)
c
c     calculate the out-of-plane distance energy term
c
      do iopdist = 1, nopdistloc
         i = opdistglob(iopdist)
         ia = iopd(1,i)
         ialoc = loc(ia)
         ib = iopd(2,i)
         ic = iopd(3,i)
         id = iopd(4,i)
         force = opdk(i)
c
c     decide whether to compute the current interaction
c
         if (use_group)  call groups (fgrp,ia,ib,ic,id,0,0)
         proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the central and peripheral atoms
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
c
c     compute the out-of-plane distance for central atom
c
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xbd = xib - xid
            ybd = yib - yid
            zbd = zib - zid
            xcd = xic - xid
            ycd = yic - yid
            zcd = zic - zid
            if (use_polymer) then
               call image (xad,yad,zad)
               call image (xbd,ybd,zbd)
               call image (xcd,ycd,zcd)
            end if
            xt = ybd*zcd - zbd*ycd
            yt = zbd*xcd - xbd*zcd
            zt = xbd*ycd - ybd*xcd
            rt2 = xt*xt + yt*yt + zt*zt
            dt2 = (xt*xad + yt*yad + zt*zad)**2 / rt2
            dt = sqrt(dt2)
            dt3 = dt2 * dt
            dt4 = dt2 * dt2
c
c     find the out-of-plane distance energy
c
            e = opdunit * force * dt2
     &             * (1.0_ti_p+copd*dt+qopd*dt2+popd*dt3+sopd*dt4)

            if(use_group) then
              call groups(fgrp,ia,ib,ic,id,0,0)
              e = e*fgrp
            endif
c
c     scale the interaction based on its group membership
c
            if (use_group) e = e * fgrp
c
c     increment the total out-of-plane distance energy
c
            neopd = neopd + 1
            eopd = eopd + e
            aeopd(ialoc) = aeopd(ialoc) + e
c
c     print a message if the energy of this interaction is large
c
            huge = (e .gt. 2.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Out-of-Plane Distance',
     &                       ' Interactions :',
     &                    //,' Type',25x,'Atom Names',18x,'Distance',
     &                       6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ic,
     &                          name(ic),id,name(id),sqrt(dt2),e
   20          format (' O-P-Dist',2x,(i7,'-',a3),f11.4,f12.4)
            end if
         end if
      end do
      return
      end
