c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ebond3  --  bond stretch energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ebond3" calculates the bond stretching energy; also
c     partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      subroutine ebond3gpu
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms    ,only: type
      use atomsMirror
      use bndpot
      use bond
      use bound
      use domdec
      use energi
      use group
      use inform
      use iounit
      use nvshmem
      use usage
      use timestat ,only: timer_enter,timer_exit,timer_ebond
     &             ,quiet_timers
      use tinheader
      implicit none
      integer i,ia,ib,ialoc,ibloc,ibond
      integer ipe,ind
      real(t_p) e,ideal,force
      real(t_p) expterm,bde
      real(t_p) dt,dt2
      real(t_p) xab,yab,zab,rab
      logical proceed
      logical header,huge
!$acc routine(image_acc) seq
c
c
c     zero out the bond energy and partitioning terms
c
      if(deb_Path) write(*,*) 'ebond3gpu'
      call timer_enter( timer_ebond )

      neb    = 0
      eb     = 0.0_re_p
c     aeb    = 0.0_ti_p
      if (rank.eq.0) then
         header = .true.
      else
         header=.false.
      end if
c
c     calculate the bond stretching energy term
c
!$acc parallel loop present(eb) async
#ifdef USE_NVSHMEM_CUDA
!$acc&         default(present) deviceptr(d_ibnd,d_bl,d_bk)
#else
!$acc&         present(x,y,z,use,loc,bndglob,ibnd,bl,bk,aeb)
#endif
      do ibond = 1, nbondloc
         i     = bndglob(ibond)
#ifdef USE_NVSHMEM_CUDA
         ipe   = (i-1)/nbond_pe
         ind   = mod((i-1),nbond_pe) +1
         ia    = d_ibnd(ipe)%pel(1,ind)
         ib    = d_ibnd(ipe)%pel(2,ind)
         ideal = d_bl  (ipe)%pel(ind)
         force = d_bk  (ipe)%pel(ind)
#else
         ia    = ibnd(1,i)
         ib    = ibnd(2,i)
         ideal = bl(i)
         force = bk(i)
#endif
         ialoc = loc(ia)
         ibloc = loc(ib)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            if (use_polymer)  call image_acc (xab,yab,zab)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt  = rab - ideal
c
c     harmonic potential uses Taylor expansion of Morse potential
c     through the fourth power of the bond length deviation
c
            if (bndtyp_i .eq. BND_HARMONIC) then
               dt2 = dt * dt
               e = bndunit * force * dt2 * (1.0_ti_p+cbnd*dt+qbnd*dt2)
c
c     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            else if (bndtyp_i .eq. BND_MORSE) then
               expterm = exp(-2.0_ti_p*dt)
               bde = 0.25_ti_p * bndunit * force
               e = bde * (1.0_ti_p-expterm)**2
            end if
c
c     increment the total bond energy and partition between atoms
c
            neb = neb + 1
            eb = eb + e
!$acc atomic update
            aeb(ialoc) = aeb(ialoc) + 0.5_ti_p*e
!$acc atomic update
            aeb(ibloc) = aeb(ibloc) + 0.5_ti_p*e
#ifndef _OPENACC
c
c     print a message if the energy of this interaction is large
c
            huge = (e .gt. 5.0_ti_p)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Bond Stretching',
     &                       ' Interactions :',
     &                    //,' Type',14x,'Atom Names',22x,'Ideal',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ideal,rab,e
   20          format (' Bond',6x,2(i7,'-',a3),13x,2f10.4,f12.4)
            end if
#endif
         end if
      end do
!$acc update host(aeb) async

      call timer_exit( timer_ebond )
      end
