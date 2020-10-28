c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine estrbnd3  --  stretch-bend energy & analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "estrbnd3" calculates the stretch-bend potential energy;
c     also partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      module estrbnd3_inl
        contains
#include "image.f.inc"
      end module

      subroutine estrbnd3
      use action
      use analyz
      use angle
      use angpot
      use atmlst
      use atmtyp
      use atoms
      use bond
      use bound
      use domdec
      use energi
      use estrbnd3_inl
      use group
      use inform
      use iounit
      use math
      use nvshmem
      use strbnd
      use tinheader ,only:ti_p,re_p
      use usage
      use sizes, only: tinkerdebug
      implicit none
      integer i,j,k,istrbnd,istrbndloc
      integer ia,ib,ic,ibloc
      integer ipe,ind
      real(t_p) e,dr1,dr2,dt
      real(t_p) angle1
      real(t_p) force1,force2
      real(t_p) dot,cosine
      real(t_p) xia,yia,zia
      real(t_p) xib,yib,zib
      real(t_p) xic,yic,zic
      real(t_p) xab,yab,zab
      real(t_p) xcb,ycb,zcb
      real(t_p) rab,rab2
      real(t_p) rcb,rcb2
      logical proceed
      logical header,huge
c
c
c     zero out the energy component and partitioning terms
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'estrbnd3'
      neba   = 0
      eba    = 0.0_ti_p
      aeba   = 0.0_ti_p
      if (rank.eq.0) then
         header = .true.
      else
         header=.false.
      end if
c
c     calculate the stretch-bend energy term
c
!$acc parallel loop default(present) present(eba) async
      do istrbndloc = 1, nstrbndloc
         istrbnd = strbndglob(istrbndloc)
         i       = isb(1,istrbnd)
#ifdef USE_NVSHMEM
         ipe     =     (i-1)/nangle_pe
         ind     = mod((i-1),nangle_pe) +1
         ia      = d_iang(ipe)%pel(1,ind)
         ib      = d_iang(ipe)%pel(2,ind)
         ic      = d_iang(ipe)%pel(3,ind)
#else
         ia      = iang(1,i)
         ib      = iang(2,i)
         ic      = iang(3,i)
#endif
         ibloc   = loc(ib)
         force1  = sbk(1,istrbnd)
         force2  = sbk(2,istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
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
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            if (use_polymer) then
               call image_inl (xab,yab,zab)
               call image_inl (xcb,ycb,zcb)
            end if
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            if (rab2.ne.0.0_ti_p .and. rcb2.ne.0.0_ti_p) then
               rab    = sqrt(rab2)
               rcb    = sqrt(rcb2)
               dot    = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0_ti_p,max(-1.0_ti_p,cosine))
               angle1 = radian * acos(cosine)
               dt     = angle1 - anat(i)
c
c     get the stretch-bend interaction energy
c
               j   = isb(2,istrbnd)
               k   = isb(3,istrbnd)
#ifdef USE_NVSHMEM_CUDA
               ipe = (j-1)/nbond_pe
               ind = mod((j-1),nbond_pe) +1
               dr1 = rab - d_bl(ipe)%pel(ind)
               ipe = (k-1)/nbond_pe
               ind = mod((k-1),nbond_pe) +1
               dr2 = rcb - d_bl(ipe)%pel(ind)
#else
               dr1 = rab - bl(j)
               dr2 = rcb - bl(k)
#endif
               e   = stbnunit * (force1*dr1+force2*dr2) * dt
c
c     increment the total stretch-bend energy
c
               neba = neba + 1
               eba  = eba + e
!$acc atomic
               aeba(ibloc) = aeba(ibloc) + e
#ifndef _OPENACC
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 2.0_ti_p)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Stretch-Bend',
     &                          ' Interactions :',
     &                       //,' Type',18x,'Atom Names',18x,'dSB 1',
     &                          5x,'dSB 2',6x,'Energy',/)
                  end if
                  write (iout,20)  ia,name(ia),ib,name(ib),
     &                             ic,name(ic),dr1*dt,dr2*dt,e
   20             format (' StrBend',3x,3(i7,'-',a3),2x,2f10.4,f12.4)
               end if
#endif
            end if
         end if
      end do
!$acc update self(aeba) async
      end
