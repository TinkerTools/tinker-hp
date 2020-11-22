c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1  --  charge-charge energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine echarge1   
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      call echarge1c   
c
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1c  --  Ewald charge derivs via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1c" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a particle mesh Ewald summation and a neighbor list
c
c
      subroutine echarge1c   
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use deriv
      use domdec
      use energi
      use ewald
      use group
      use iounit
      use inter
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use virial
      use mpi
      use vec
      use vec_elec

      implicit none
      integer k
      integer ncloop8,ncloop16
      real*8  f,fs,e
      real*8  xi,yi,zi
      real*8  xd,yd,zd
      real*8  term
      real*8  time0,time1

!DIR$ ATTRIBUTES ALIGN:64:: iglobvec
      integer iglobvec(nionlocloop)
!DIR$ ATTRIBUTES ALIGN:64::ivec
      integer ivec(nionlocloop)
!DIR$ ATTRIBUTES ALIGN:64::pchgvec
      real*8 pchgvec(nionlocloop)
!DIR$ ATTRIBUTES ALIGN:64:: devec
      real*8  devec(nionlocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dec1vec
      real*8  dec1vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dec2vec
      real*8  dec2vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dec3vec
      real*8  dec3vec(nblocloop)

      if(rank.eq.0.and.tinkerdebug) write(*,*) 'echarge1cvec'

c     return if no charge
      if (nion .eq. 0)  return
c
c     zero out the Ewald summation energy and derivatives
c
      ec  = 0.0d0
c
c     compute the reciprocal space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.gt.ndir-1))
     &   then
         time0 = mpi_wtime()
         call ecrecip1
         time1 = mpi_wtime()
         timerec = timerec + time1 - time0
      end if
c
c     compute the real space part of the Ewald summation
c
      if ((.not.(use_pmecore)).or.(use_pmecore).and.(rank.le.ndir-1))
     &   then
         e = 0.0d0
         time0 = mpi_wtime()
         call ecreal1c
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
         do k = 1, nblocloop
            dec1vec(k) = 0.0d0
            dec2vec(k) = 0.0d0
            dec3vec(k) = 0.0d0
         enddo
c
c     set conversion factor, cutoff and switching coefficients
c
         f = electric / dielec
c
c     compute the Ewald self-energy term over all the atoms
c
         fs = -f * aewald / sqrtpi

!DIR$ ASSUME (mod(nionlocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1,nionlocloop
            if (k.le.nionloc) then 
               e = e + fs * pchg(chgglob(k))**2
            endif
         enddo

c
c     compute the cell dipole boundary correction term
c

         term = (2.0d0/3.0d0) * f * (pi/volbox)
         if (boundary .eq. 'VACUUM') then
            xd = 0.0d0
            yd = 0.0d0
            zd = 0.0d0
!DIR$ ASSUME (mod(nionlocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
            do k = 1,nionlocloop
               if (k.le.nionloc) then 
                  iglobvec(k) = iion(chgglob(k))
                  pchgvec (k) = pchg(chgglob(k))
               else
                  iglobvec(k) = iion(chgglob(nionloc))
                  pchgvec (k) = pchg(chgglob(nionloc))
               endif
               ivec       (k) = loc(iglobvec(k))
               xd = xd + pchgvec(k) * x (iglobvec(k))
               yd = yd + pchgvec(k) * y (iglobvec(k))
               zd = zd + pchgvec(k) * z (iglobvec(k))
            enddo

!DIR$ ASSUME (mod(nionlocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
            do k = 1,nionlocloop
               dec1vec(ivec(k)) =  dec1vec(ivec(k)) + xd
     &                           * 2.0d0 * term * pchgvec(k)
               dec2vec(ivec(k)) =  dec2vec(ivec(k)) + yd
     &                           * 2.0d0 * term * pchgvec(k)
               dec3vec(ivec(k)) =  dec3vec(ivec(k)) + zd
     &                           * 2.0d0 * term * pchgvec(k)
            enddo
            e = e + term * (xd**2 + yd**2 + zd**2 )
         endif
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, nblocloop
            if(k.le.nbloc) then
               dec(1,k) = dec(1,k) + dec1vec(k)
               dec(2,k) = dec(2,k) + dec2vec(k)
               dec(3,k) = dec(3,k) + dec3vec(k)
            endif
         enddo
         ec = ec + e
         time1 = mpi_wtime()
         timereal = timereal + time1 - time0
      endif
      return
      end
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine ecreal1c     --  Ewald real space derivs via list  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "ecreal1" calculates the real space portion of the charge-charge
c      interaction energy and first derivatives with respect to Cartesian
c      coordinates using a particle mesh Ewald summation and a neighbor list
c
c
      subroutine ecreal1c
      use atmlst
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use deriv
      use domdec
      use energi
      use ewald
      use group
      use iounit
      use inter
      use math
      use molcul
      use neigh
      use potent
      use shunt
      use timestat
      use usage
      use virial
      use mpi
      use vec
      use vec_elec
      use vec_charge
      use utilvec

      implicit none
      integer i,j
      integer ii,iichg,iglob,k,kk
      integer nnchg,nnchg1,nnchg2,countsel
      integer ncloop8,ncloop16
      integer kkchdefault,kglobdefault,kdefault
      real*8 f
      real*8 fs,fi,fact
      real*8 xi,yi,zi
      real*8 xd,yd,zd
      real*8 vxx, vxy, vxz, vyy, vyz, vzz
      real*8 term
      real*8 half,one
!DIR$ ATTRIBUTES ALIGN:64:: dec1vec
      real*8  dec1vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dec2vec
      real*8  dec2vec(nblocloop)
!DIR$ ATTRIBUTES ALIGN:64:: dec3vec
      real*8  dec3vec(nblocloop)

      character*6 mode
      if(rank.eq.0.and.tinkerdebug)write(*,*) 'ecreal1cvec'
1000  format(' Warning, system moved too much since last neighbor list'
     $   ' update, try lowering nlupdate')

!DIR$ ASSUME (mod(nblocloop,16).eq.0)
      do k = 1, nblocloop
         dec1vec(k) = 0.0d0
         dec2vec(k) = 0.0d0
         dec3vec(k) = 0.0d0
      enddo
      half = 0.5d0
      one  = 1.0d0
      fact = 2.0d0 * aewald / sqrtpi
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)

c     compute the real space Ewald energy and first derivatives
c
      MAINLOOP:
     &do ii = 1, nionlocnl
         iichg = chgglobnl(ii)
         iglob = iion(iichg)
         i = loc(iglob)
         if (i.eq.0) then
           write(iout,1000)
           cycle MAINLOOP
         end if
         countsel = 0
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         fi = f * pchg(iichg)
!DIR$ ASSUME_ALIGNED mask1:64
         nnchg    =  nelst(ii)
         if(nnchg.eq.0)      cycle MAINLOOP ! No atoms selected
         ncloop8  = (int(nnchg / 8 ) + 1) * 8 ! First multiple of 8
         ncloop16 = (int(nnchg / 16) + 1) * 16! First multiple of 16

         ncloop8  = merge(nnchg,ncloop8 , mod(nnchg,8 ).eq.0)
         ncloop16 = merge(nnchg,ncloop16, mod(nnchg,16).eq.0)

c     set the default value for indices

         kkchdefault  = elst(nnchg,ii)    ! Safe value to point to
         kglobdefault = iion(kkchdefault) ! Safe value to point to
         kdefault     = 0                 ! Exclusion value by default
         
!DIR$ ASSUME (mod(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, ncloop16
            if (k.le.nnchg) then
               kkchgvec(k) = elst(k,ii)
               kglobvec(k) = iion(kkchgvec(k))
               kvec    (k) = loc(kglobvec (k))
            else
               kkchgvec(k) = kkchdefault
               kglobvec(k) = kglobdefault
               kvec    (k) = kdefault
            endif
         enddo
c
c     compute the energy contribution for this interaction
c   
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, ncloop16
            xposvec1(k)  = xi - x(kglobvec(k))
            yposvec1(k)  = yi - y(kglobvec(k))
            zposvec1(k)  = zi - z(kglobvec(k))
         enddo

         call image3dvec(xposvec1,yposvec1,zposvec1,ncloop16)
c
c     find energy for interactions within real space cutoff
c
         kk = 0
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, ncloop16
c          if (r2vec(k).le.off2) then
           if (     kvec(k) /= 0
     &         .and.xposvec1(k)**2 + yposvec1(k)**2 + zposvec1(k)**2
     &              <=off2
     &         .and.k.le.nnchg) then

              kk = kk + 1

              kkchgvec1(kk) = kkchgvec(k)
              kglobvec1(kk) = kglobvec(k)
              kvec1    (kk) = kvec    (k)

              xposvec  (kk) = xposvec1 (k)
              yposvec  (kk) = yposvec1 (k)
              zposvec  (kk) = zposvec1 (k)
           endif
         enddo
         nnchg1 = kk
         if(nnchg1.eq.0) cycle MAINLOOP ! No atoms selected

         ncloop8  = (int(nnchg1 / 8 ) + 1) * 8 ! First multiple of 8
         ncloop16 = (int(nnchg1 / 16) + 1) * 16! First multiple of 16

         ncloop8  = merge(nnchg1,ncloop8 , mod(nnchg1,8 ).eq.0)
         ncloop16 = merge(nnchg1,ncloop16, mod(nnchg1,16).eq.0)

c
c     set exclusion coefficients for connected atoms
c
        call setscale(iglob,kglobvec1,ncloop16,'cscale',scalevec)
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, ncloop16
            fikvec(k)    =  fi *pchg(kkchgvec1(k))
         enddo
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
         do k = 1, ncloop16
            rvec(k)   = (xposvec(k)**2 + yposvec(k)**2 + zposvec(k)**2)
     &                   **   half
            rewvec(k) = aewald * rvec(k)
         enddo
         call vderfc(ncloop16,rewvec,erfvec)
!DIR$ ASSUME (MOD(ncloop8,8).eq.0)
         do k = 1, ncloop8
            term1vec(k)  = erfvec(k) + scalevec(k) - 1.0d0
            invrvec(k)   = rvec(k) ** ( - one )
            invrbvec(k)  = (rvec(k) + ebuffer) ** ( - one )
            invrb2vec(k) = invrbvec(k) ** 2
         enddo
c
c     form the chain rule terms for derivative expressions
c

!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, ncloop16
            if (k.le.nnchg1) then
               evec (k) =   fikvec(k) * invrbvec(k) * term1vec(k)
               devec(k) = - fikvec(k)
     &                     * (  term1vec(k) * invrb2vec(k)
     &                        + fact * exp( - rewvec(k)**2)
     &                          * invrvec(k)
     &                       ) * invrvec(k)
            else
               devec(k) = 0.0d0
               evec (k) = 0.0d0
            endif
         enddo
!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, ncloop16
            decxvec(k) = devec(k) * xposvec(k)
            decyvec(k) = devec(k) * yposvec(k)
            deczvec(k) = devec(k) * zposvec(k)
         enddo
         vxx = 0.0d0
         vxy = 0.0d0
         vxz = 0.0d0
         vyy = 0.0d0
         vyz = 0.0d0
         vzz = 0.0d0
c
c    increment the internal virial tensor components and the derivative expressions
c

!DIR$ ASSUME (MOD(ncloop16,16).eq.0)
!DIR$ VECTOR ALIGNED
         do k = 1, ncloop16
               vxx = vxx + xposvec(k) * decxvec(k)
               vxy = vxy + yposvec(k) * decxvec(k)
               vxz = vxz + zposvec(k) * decxvec(k)
               vyy = vyy + yposvec(k) * decyvec(k)
               vyz = vyz + zposvec(k) * decyvec(k)
               vzz = vzz + zposvec(k) * deczvec(k)

               dec1vec(i) = dec1vec(i) + decxvec(k)
               dec2vec(i) = dec2vec(i) + decyvec(k)
               dec3vec(i) = dec3vec(i) + deczvec(k)
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
c
c     increment the overall energy and the derivative expressions
c

!DIR$ ASSUME (MOD(ncloop8,8).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, ncloop8
            ec = ec + evec   (k)
            if (k.le.nnchg1) then
               dec1vec(kvec1(k)) = dec1vec(kvec1(k)) - decxvec(k)
               dec2vec(kvec1(k)) = dec2vec(kvec1(k)) - decyvec(k)
               dec3vec(kvec1(k)) = dec3vec(kvec1(k)) - deczvec(k)
            endif
         enddo
      enddo MAINLOOP
c
c     return calculated values in the original array
c
!DIR$ ASSUME (mod(nblocloop,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      do k = 1, nblocloop
         if(k.le.nbloc) then
            dec(1,k) = dec(1,k) + dec1vec(k)
            dec(2,k) = dec(2,k) + dec2vec(k)
            dec(3,k) = dec(3,k) + dec3vec(k)
         endif
      enddo
      return
      end
