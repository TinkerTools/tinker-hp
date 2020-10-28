c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine eopdist  --  out-of-plane distance energy  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "eopdist" computes the out-of-plane distance potential
c     energy at trigonal centers via the central atom height
c
c
#include "tinker_precision.h"
      subroutine eopdist
      use angpot
      use atmlst
      use atoms
      use bound
      use energi
      use group
      use opdist
      use tinheader
      use usage
      implicit none
      integer i,ia,ib,ic,id,iopdist
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
      logical proceed
c
c
c     zero out the out-of-plane distance energy component
c
      eopd = 0.0_ti_p
c
c     calculate the out-of-plane distance energy term
c
      do iopdist = 1, nopdistloc
         i = opdistglob(iopdist)
         ia = iopd(1,i)
         ib = iopd(2,i)
         ic = iopd(3,i)
         id = iopd(4,i)
         force = opdk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the defining atoms
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
c
c     increment the total out-of-plane distance energy
c
            eopd = eopd + e
         end if
      end do
      return
      end
