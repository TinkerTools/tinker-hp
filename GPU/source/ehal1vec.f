c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1vec  --  buffered 14-7 energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c
c
#include "tinker_macro.h"
      subroutine ehal1vec

      use energi
      use virial
      use vdwpot
      implicit none
      real(t_p) elrc,vlrc
c
c
c     choose the method for summing over pairwise interactions
c
      call ehal1cvec
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr1 (elrc,vlrc)
         ev = ev + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1c  --  buffered 14-7 vdw derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1c" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c     using a pairwise neighbor list
c
c
      subroutine ehal1cvec
      use atmlst
      use atoms
      use bound
      use couple
      use deriv
      use domdec
      use energi
      use group
      use inter
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      use virial
      use vdw
      use vdwpot
      use vec
      use vec_vdw
      use utilvec

      implicit none
      integer i,k,kk
      integer iglob,iivdw
      integer ii,iv,it,ivloc
      integer nnvlst,nnvlst2,nvloop8,nvloop16
      integer iglobdefault
      integer kglobdefault
      integer ktdefault
      integer nn14,nn14loop
      real(t_p) vxx,vxy,vxz,vyy,vyz,vzz
      real(t_p) xi,yi,zi,redi,e
      real(t_p) half,one

!DIR$ ATTRIBUTES ALIGN:64::iglobvec
      integer iglobvec(nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64::ivec
      integer ivec(nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64::ivvec
      integer ivvec(nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64::xred
      real(t_p)  xred(nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64::yred
      real(t_p)  yred(nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64::zred
      real(t_p)  zred(nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: rdnvec
      real(t_p)  rdnvec(nvdwblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dev1vec
      real(t_p)  dev1vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dev2vec
      real(t_p)  dev2vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dev3vec
      real(t_p)  dev3vec(nblocloop)

      character*10 mode
c
 1000 format(' Warning, system moved too much since last neighbor list'
     $   'update, try lowering nlupdate VDW')
      if(rank.eq.0.and.tinkerdebug) write (*,*) 'ehal1cvec'
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0_re_p

!DIR$ ASSUME (mod(nblocloop,16).eq.0)
      do k = 1, nblocloop
         dev1vec(k) = 0.0_ti_p
         dev2vec(k) = 0.0_ti_p
         dev3vec(k) = 0.0_ti_p
      enddo

      half = 0.5_ti_p
      one  = 1.0_ti_p
      e    = 0.0_ti_p

c     set default values to point to for various indices
      iglobdefault = ivdw(vdwglob (nvdwbloc))
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
            iglobvec(k) = iglobdefault
         endif
      enddo
c
c     apply any reduction factor to the atomic coordinates
c
!DIR ASSUME (mod(nvdwblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1 , nvdwblocloop
         ivec    (k) = loc  (iglobvec(k))
         ivvec   (k) = ired (iglobvec(k))
         rdnvec  (k) = kred (iglobvec(k))
      enddo

!DIR ASSUME (mod(nvdwblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR G2S

      do k = 1 , nvdwblocloop
         xred(ivec(k)) =  rdnvec(k)           * x(iglobvec(k))
     &                  + (1.0_ti_p - rdnvec(k)) * x(ivvec   (k))
       
         yred(ivec(k)) =  rdnvec(k)           * y(iglobvec(k))
     &                  + (1.0_ti_p - rdnvec(k)) * y(ivvec   (k))
       
         zred(ivec(k)) =  rdnvec(k)           * z(iglobvec(k))
     &                  + (1.0_ti_p - rdnvec(k)) * z(ivvec   (k))
      enddo

c
c     find van der Waals energy and derivatives via neighbor list
c
c
      MAINLOOP:
     &do ii = 1, nvdwlocnl
         iivdw = vdwglobnl(ii)
         iglob = ivdw(iivdw)
         i = loc(iglob)
         if (i.eq.0) then
           write(iout,1000)
           cycle MAINLOOP
         end if
         iv   = ired(iglob)
         redi = merge (1.0_ti_p, kred(iglob),(iglob.eq.iv))
         ivloc = loc(iv)
         it = jvdw(iglob)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)

         nnvlst          = nvlst(ii)
         if(nnvlst.eq.0) cycle MAINLOOP
         nvloop8  = (int(nnvlst / 8 ) + 1) * 8 ! First multiple of 8
         nvloop16 = (int(nnvlst / 16) + 1) * 16! First multiple of 16

         nvloop8  = merge(nnvlst,nvloop8 , mod(nnvlst,8 ).eq.0)
         nvloop16 = merge(nnvlst,nvloop16, mod(nnvlst,16).eq.0)

         kglobdefault = vlst(nnvlst,ii)
         ktdefault    = 1
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ NOFUSION
         do k = 1, nvloop16
            if (k.le.nnvlst) then 
               kglobvec(k) = vlst(k,ii)
            else
               kglobvec(k) = kglobdefault
            endif
         enddo
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
!DIR$ NOFUSION
         do k = 1, nvloop16
            kbisvec  (k) = loc (kglobvec(k))
            kvvec    (k) = ired(kglobvec(k))
            kvlocvec (k) = loc (kvvec   (k))
            kbisvec2 (k) = kbisvec      (k)
            kvlocvec2(k) = kvlocvec     (k)
         enddo
c
c     compute the energy contribution for this interaction
c
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
         do k = 1, nvloop16
            xposvec1(k) = xi - xred(kbisvec(k))
            yposvec1(k) = yi - yred(kbisvec(k))
            zposvec1(k) = zi - zred(kbisvec(k))
         enddo

         call image3dvec(xposvec1,yposvec1,zposvec1,nvloop16)

c
c     decide whether to compute the current interaction
c     and check for an interaction distance less than the cutoff
c
         kk = 0
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
         do k = 1, nvloop16
            if (      kvlocvec(k) /= 0
     &          .and.kbisvec (k) <= nbloc
     &          .and.kvlocvec(k) <= nbloc
     &          .and.xposvec1(k)**2 + yposvec1(k)**2 + zposvec1(k)**2
     &               <=off2
     &          .and. k.le.nnvlst
     &         ) then
               kk = kk + 1
               kglobvec2(kk) = kglobvec(k)
               kbisvec2 (kk) = kbisvec (k)
               kvvec2   (kk) = kvvec   (k)
               kvlocvec2(kk) = kvlocvec(k)

               xposvec(kk) = xposvec1(k)
               yposvec(kk) = yposvec1(k)
               zposvec(kk) = zposvec1(k)
            endif
         end do
         nnvlst2 = kk
         if(nnvlst2.eq.0) cycle MAINLOOP ! No atoms selected

         nvloop8  = (int(nnvlst2 / 8 ) + 1) * 8 ! First multiple of 8
         nvloop16 = (int(nnvlst2 / 16) + 1) * 16! First multiple of 16

         nvloop8  = merge(nnvlst2,nvloop8 , mod(nnvlst2,8 ).eq.0)
         nvloop16 = merge(nnvlst2,nvloop16, mod(nnvlst2,16).eq.0)

!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
         do k = 1, nvloop8
            rik2vec(k) = xposvec(k)**2 + yposvec(k)**2 + zposvec(k)**2
            rikvec (k) = rik2vec (k) ** half
            invrikvec(k) = 1.0_ti_p/rikvec(k)
         enddo
c
c     compute the energy contribution for this interaction
c
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX=160
!DIR$ SIMD
         do k = 1, nvloop16
            if(k.le.nnvlst2) then
               ktvec2   (k) = jvdw(kglobvec2(k))
            else
               ktvec2 (k) =  ktdefault
            endif
        enddo
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!!DIR$ NOFUSION
!DIR$ SIMD
        do k = 1, nvloop16
            if (kglobvec2(k).ne.kvvec2(k)) then
               redkvec(k) = kred(kglobvec2(k))
            else
               redkvec(k) = 1.0_ti_p
            endif
         enddo
c
c     set interaction scaling coefficients for connected atoms
c
         call setscale(iglob,kglobvec2,nvloop16,"vscale",vscalevec)
c
c       get mask for radmin epsilon
c
         nn14   = n14(iglob)
         nn14loop   = merge(nn14,(int(nn14 / 16) + 1) * 16,
     &                     mod(nn14,16).eq.0)
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ ASSUME (MOD(nn14loop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX=80
!DIR$ SIMD
         do k = 1, nvloop8
            mask1(k) = .false.
            do kk = 1, nn14loop
               if(kglobvec2(k) == i14(kk,iglob).and..not.mask1(k)) then
                  mask1(k) = .true.
               endif
            enddo
         enddo
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop8
            if (mask1(k)) then !mask epsilon and radmin
               rvvec2(k)  = radmin4 (ktvec2(k),it)
               epsvec2(k) = epsilon4(ktvec2(k),it)
            else
               rvvec2(k)  = radmin  (ktvec2(k),it)
               epsvec2(k) = epsilon (ktvec2(k),it)
            endif
         enddo
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop8
            rv7vec (k) = rvvec2 (k)**7
            rik3vec(k) = rik2vec(k) * rikvec(k)
            rik4vec(k) = rik3vec(k) * rikvec(k)
            rik5vec(k) = rik4vec(k) * rikvec(k)
            rik6vec(k) = rik5vec(k) * rikvec(k)
            rik7vec(k) = rik6vec(k) * rikvec(k)
         enddo
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
         do k = 1, nvloop8
            rv7orhovec(k) =  rv7vec (k)
     &                    * (rik7vec(k) + ghal * rv7vec(k)) ** ( - one )
            tauvec (k)    =  (dhal + 1.0_ti_p)
     &                     * (rikvec(k) + dhal * rvvec2(k)) ** ( - one )
            tau7vec(k)    =  tauvec   (k)      ** 7.0_ti_p
         enddo
    
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop8
            dtauvec(k)    =  tauvec   (k)   / (dhal + 1.0_ti_p)
            gtauvec(k)    =  epsvec2  (k)   * tau7vec(k)
     &                     * rik6vec  (k)   * (ghal + 1.0_ti_p)
     &                     * rv7orhovec(k)  * rv7orhovec(k)
     &                     * vscalevec(k)
            evec(k)  =  epsvec2(k) * tau7vec(k) * rv7vec(k)
     &                * ((ghal + 1.0_ti_p) * rv7orhovec(k) - 2.0_ti_p)
     &                * vscalevec(k)
    
            devec(k) = - 7.0_ti_p * (dtauvec(k) * evec(k) + gtauvec(k))
         enddo
c
c     use energy switching if near the cutoff distance
c
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
         do k = 1, nvloop8
            if(rik2vec(k) > cut2) then ! mask energy switch
               tapervec(k) =  c5 * rik5vec(k) + c4 * rik4vec(k) 
     &                      + c3 * rik3vec(k) + c2 * rik2vec(k)
     &                      + c1 * rikvec(k) + c0           
               dtapervec(k)  =  5.0_ti_p * c5 * rik4vec(k)
     &                        + 4.0_ti_p * c4 * rik3vec(k)
     &                        + 3.0_ti_p * c3 * rik2vec(k)
     &                        + 2.0_ti_p * c2 * rikvec(k)
     &                        +         c1
            else
               tapervec(k)  = 1.0_ti_p
               dtapervec(k) = 0.0_ti_p
            endif
         enddo
c
c     increment the total van der Waals energy 
c
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop8
            if (k.le.nnvlst2) then
               devec(k) = (  evec (k) * dtapervec(k)
     &                     + devec(k) * tapervec (k)) * invrikvec(k)
               evec(k)  =    evec (k) * tapervec (k)
            else
               devec(k) = 0.0_ti_p
               evec (k) = 0.0_ti_p
            endif
         enddo
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop8
            dedxvec(k) = devec(k) * xposvec(k)
            dedyvec(k) = devec(k) * yposvec(k)
            dedzvec(k) = devec(k) * zposvec(k)
         enddo
c
c     increment the internal virial tensor components
c
         vxx = 0.0_ti_p
         vxy = 0.0_ti_p
         vxz = 0.0_ti_p
         vyy = 0.0_ti_p
         vyz = 0.0_ti_p
         vzz = 0.0_ti_p
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
         do k =1, nvloop8
            vxx = vxx + xposvec(k) * dedxvec(k)
            vxy = vxy + yposvec(k) * dedxvec(k)
            vxz = vxz + zposvec(k) * dedxvec(k)
            vyy = vyy + yposvec(k) * dedyvec(k)
            vyz = vyz + zposvec(k) * dedyvec(k)
            vzz = vzz + zposvec(k) * dedzvec(k)
         enddo
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz

!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop8
            e = e + evec(k)
            devxvec(k)  = redkvec(k) * dedxvec(k)
            devyvec(k)  = redkvec(k) * dedyvec(k)
            devzvec(k)  = redkvec(k) * dedzvec(k)

            devvxvec(k) = dedxvec(k) - devxvec(k)
            devvyvec(k) = dedyvec(k) - devyvec(k)
            devvzvec(k) = dedzvec(k) - devzvec(k)
         enddo
c
c     increment the van der Waals derivatives
c
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
         do k = 1, nvloop8
            dev1vec(kbisvec2(k)) = dev1vec(kbisvec2(k)) - devxvec(k)
            dev2vec(kbisvec2(k)) = dev2vec(kbisvec2(k)) - devyvec(k)
            dev3vec(kbisvec2(k)) = dev3vec(kbisvec2(k)) - devzvec(k)
         enddo
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
         do k = 1, nvloop8
            dev1vec(kvlocvec2(k)) = dev1vec(kvlocvec2(k)) - devvxvec(k)
            dev2vec(kvlocvec2(k)) = dev2vec(kvlocvec2(k)) - devvyvec(k)
            dev3vec(kvlocvec2(k)) = dev3vec(kvlocvec2(k)) - devvzvec(k)
         enddo
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, nvloop8
            dev1vec(i) =  dev1vec(i) + redi * dedxvec(k)
            dev2vec(i) =  dev2vec(i) + redi * dedyvec(k)
            dev3vec(i) =  dev3vec(i) + redi * dedzvec(k)
    
           dev1vec(ivloc)= dev1vec(ivloc) + (1.0_ti_p - redi)*dedxvec(k)
           dev2vec(ivloc)= dev2vec(ivloc) + (1.0_ti_p - redi)*dedyvec(k)
           dev3vec(ivloc)= dev3vec(ivloc) + (1.0_ti_p - redi)*dedzvec(k)
         enddo
      end do MAINLOOP
c
c     return calculated values in the original array
c
!DIR$ ASSUME (mod(nblocloop,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, nblocloop
         if(k.le.nbloc) then
            dev(1,k) = dev1vec(k)
            dev(2,k) = dev2vec(k)
            dev(3,k) = dev3vec(k)
         endif
      enddo
      ev = ev + e
      return
      end
