
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal3  --  buffered 14-7 vdw energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehal3" calculates the buffered 14-7 van der Waals energy
c     and partitions the energy among the atoms
c
c
#include "tinker_precision.h"
      subroutine ehal3vec
      use analyz
      use atoms
      use domdec
      use energi
      use inform
      use iounit
      use tinheader ,only:ti_p,re_p
      use vdwpot

      implicit none
      integer i
      real(t_p) elrc,aelrc
c
c
c     choose the method for summing over pairwise interactions
c
      call ehal3cvec
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr.and.rank.eq.0) then
         call evcorr (elrc)
         ev = ev + elrc
         aelrc = elrc / real(n,t_p)
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do i = 1, nblocloop
            aev(i) = aev(i) + aelrc
         end do
         if (verbose .and. elrc.ne.0.0_ti_p) then
            write (iout,10)  elrc
   10       format (/,' Long Range vdw Correction :',9x,f12.4)
         end if
      end if
      return
      end

c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal3cvec  --  buffered 14-7 analysis via list  ##
c     ##                                                             ##
c     #################################################################
c
c     "ehal3cvec" calculates the buffered 14-7 van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine ehal3cvec
      use action
      use analyz
      use atmlst
      use atmtyp
      use atoms
      use bound
      use couple
      use domdec
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      use vec
      use vec_vdw
      use utilvec

      implicit none
      integer i,j,k,iglob
      integer ii,iv,it,iivdw
      integer kk,nn14,nn14loop
      integer nnvlst,nnvlst1,nvloop8,nvloop16
      integer iglobdefault
      integer kglobdefault

      real(t_p) xi,yi,zi
      real(r_p) e
      real(t_p) half,one
      !DIR$ ATTRIBUTES ALIGN:64::iglobvec,ivec,ivvec,rdnvec
      integer iglobvec(nvdwblocloop)
      integer ivec    (nvdwblocloop)
      integer ivvec   (nvdwblocloop)
      real(t_p)  rdnvec  (nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64::xred
      real(t_p)  xred(nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64::yred
      real(t_p)  yred(nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64::zred
      real(t_p)  zred(nvdwblocloop)

      character*10 mode
c
c     zero out the van der Waals energy and partitioning terms
c
      if(rank.eq.0.and.tinkerdebug) write(*,*) 'ehal3cvec'
      nev          = 0
      e            = 0.0_re_p
      half         = 0.5_ti_p
      one          = 1.0_ti_p
      iglobdefault = vlst(1,1)
      kglobdefault = vlst(1,1)

c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
!DIR$ ASSUME (mod(nvdwblocloop,16).eq.0)
!DIR$ NOFUSION
!DIR$ VECTOR ALIGNED
      do k = 1,nvdwblocloop
         if (k.le.nvdwbloc) then
            iglobvec (k) = ivdw (vdwglob (k))
         else
            iglobvec (k) = iglobdefault ! good value to point to
         endif
      enddo

!DIR ASSUME (mod(nvdwblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1 , nvdwblocloop
         ivec  (k) = loc  (iglobvec(k))
         ivvec (k) = ired (iglobvec(k))
         rdnvec(k) = kred (iglobvec(k))
      enddo
!DIR ASSUME (mod(nvdwblocloop,16).eq.0)
!DIR$ NOFUSION
c
c     apply any reduction factor to the atomic coordinates
c

!DIR ASSUME (mod(nvdwblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR G2S
      do k = 1 , nvdwblocloop
         xred(ivec(k)) =  rdnvec (k)          * x(iglobvec(k))
     &                  + (one - rdnvec(k)) * x(ivvec(k))

         yred(ivec(k)) =  rdnvec (k)          * y(iglobvec(k))
     &                  + (one - rdnvec(k)) * y(ivvec(k))

         zred(ivec(k)) =  rdnvec (k)          * z(iglobvec(k))
     &                  + (one - rdnvec(k)) * z(ivvec(k))
      enddo
c
c     find the van der Waals energy via neighbor list search
c
      MAINLOOP:
     &do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i     = loc(iglob)
         if (i.eq.0) then
           cycle MAINLOOP
         end if
         iv    = ired(iglob)
         it    = jvdw(iglob)
         xi    = xred(i)
         yi    = yred(i)
         zi    = zred(i)
c
c
c     fill vectors to compute interactions involving this atom
c
         nnvlst   = nvlst(ii)

         nvloop8  = (int(nnvlst / 8 ) + 1) * 8 ! First multiple of 8
         nvloop16 = (int(nnvlst / 16) + 1) * 16! First multiple of 16
         nvloop8  = merge(nnvlst,nvloop8 , mod(nnvlst,8 ).eq.0)
         nvloop16 = merge(nnvlst,nvloop16, mod(nnvlst,16).eq.0)

!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop16
            if (k.le.nnvlst) then
               kglobvec(k) = vlst(k,ii)
            else
               kglobvec(k) = kglobdefault !good value to point to by default
            endif
            kglobvec1  (k) = kglobdefault   !good value to point to by default
            kbisvec    (k) = loc(kglobvec (k))
         enddo
c
c     compute the energy contribution for this interaction
c

!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, nvloop16
            xposvec1(k) = xi - xred(kbisvec(k))
            yposvec1(k) = yi - yred(kbisvec(k))
            zposvec1(k) = zi - zred(kbisvec(k))
         enddo
         call image3dvec(xposvec1,yposvec1,zposvec1,nvloop16)

!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ NOFUSION
         do k = 1, nvloop8
            rik2vec1(k) =  xposvec1(k)**2 + yposvec1(k)**2
     &                   + zposvec1(k)**2
         enddo
c
c     check for an interaction distance less than the cutoff
c
         kk = 0
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop16
           if (rik2vec1(k).le.off2.and.k.le.nnvlst) then
             kk = kk + 1
             kglobvec1(kk) = kglobvec(k)
             rik2vec  (kk) = rik2vec1 (k)
           endif
         end do
         nnvlst1  = kk

         nvloop8  = (int(nnvlst1 / 8 ) + 1) * 8 ! First multiple of 8
         nvloop16 = (int(nnvlst1 / 16) + 1) * 16! First multiple of 16
         nvloop8  = merge(nnvlst1,nvloop8 , mod(nnvlst1,8 ).eq.0)
         nvloop16 = merge(nnvlst1,nvloop16, mod(nnvlst1,16).eq.0)

!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, nvloop16
            ktvec(k) = jvdw(kglobvec1(k))
         enddo

c
c     set interaction scaling coefficients for connected atoms
c
c        vscalevec = setscale2(iglob,kglobvec1,nvloop16,'vscale')
         call setscale(iglob,kglobvec1,nvloop16,'vscale',vscalevec)
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, nvloop8
            radminvec(k)   = radmin  (ktvec(k),it)
            radmin4vec(k)  = radmin4 (ktvec(k),it)
            epsilonvec(k)  = epsilon (ktvec(k),it)
            epsilon4vec(k) = epsilon4(ktvec(k),it)
         enddo

         nn14   = n14(iglob)
         nn14loop   = merge(nn14,(int(nn14 / 16) + 1) * 16,
     &                     mod(nn14,16).eq.0)
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
!DIR$ ASSUME (MOD(nn14loop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX=160
!DIR$ SIMD
         do k = 1, nvloop16
            mask1(k) = .false.
            do kk = 1, nn14loop
               if(kglobvec1(k) == i14(kk,iglob).and..not.mask1(k)) then
                  mask1(k) = .true.
               endif
            enddo
         enddo

!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX = 160
         do k = 1, nvloop8
            if (mask1(k)) then !mask epsilon and radmin
                rvvec(k) =  radmin4vec(k)
                epsvec(k) = epsilon4vec(k) * vscalevec(k)
            else
                rvvec(k) =  radminvec(k)
                epsvec(k) = epsilonvec(k)  * vscalevec(k)
            endif
         enddo

!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop8
            rv7vec (k) = rvvec(k)   ** 7.0_ti_p
            rikvec (k) = rik2vec(k) ** half
            rik3vec(k) = rik2vec(k) * rikvec(k)
            rik4vec(k) = rik3vec(k) * rikvec(k)
            rik5vec(k) = rik4vec(k) * rikvec(k)
            rik6vec(k) = rik5vec(k) * rikvec(k)
            rik7vec(k) = rik6vec(k) * rikvec(k)
         enddo

!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
         do k = 1, nvloop8
            invrhovec (k) = (rik7vec(k) + ghal * rv7vec(k)) ** ( - one )
            invtmpvec (k) = (rikvec (k) + dhal * rvvec (k)) ** ( - one )
            rv7orhovec(k) = rv7vec  (k)     * invrhovec(k)
            tau7vec   (k) = ((dhal + one) * invtmpvec(k)) ** 7.0_ti_p
         enddo

c
c     use energy switching if near the cutoff distance
c
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop8
            if(rik2vec(k) <= cut2) then ! mask energy switch
               tapervec(k) = one
            else
               tapervec(k) =  c5 * rik5vec(k) + c4 * rik4vec(k)
     &                      + c3 * rik3vec(k) + c2 * rik2vec(k)
     &                      + c1 * rikvec (k) + c0
            endif
            evec(k)  =  tapervec(k) * epsvec(k) * tau7vec(k) * rv7vec(k)
     &                * ((ghal + one) * rv7orhovec(k) - 2.0_ti_p)
c
c     increment the overall van der Waals energy components
c
            if(evec(k).ne.0.0_ti_p.and.k.le.nnvlst1) then
              nev = nev + 1
              e = e + evec(k)
            endif
         enddo
      enddo MAINLOOP

      ev = e  ! use ev outside the previous loop, because it prevents
              ! vectorization !

      return
      end
