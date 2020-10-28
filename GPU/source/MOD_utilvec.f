c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module utilvec  --   vectorized utility functions            ##
c     ##                                                               ##
c     ###################################################################
c
c
c
c
#include "tinker_precision.h"
      module utilvec
      use tinheader ,only: ti_p,re_p
      use cell
      use sizes
      use boxes
      use bound
      use domdec
      use sizes
      use atoms
      use vdwpot
      use chgpot
      use mplpot
      use polpot
      use polgrp
      use couple
      implicit none
      contains
c
c     ####################################################################
c     ##                                                                ##
c     ##  function imagevec3  --  compute the minimum image distance    ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "imagevec3" takes the components of pairwise distances between
c     two points in a periodic box and converts to the components
c     of the minimum image distances. Indice i designs x, y or z
c     direction
c
c     xcell    length of the a-axis of the complete replicated cell
c     ycell    length of the b-axis of the complete replicated cell
c     zcell    length of the c-axis of the complete replicated cell
c     xcell2   half the length of the a-axis of the replicated cell
c     ycell2   half the length of the b-axis of the replicated cell
c     zcell2   half the length of the c-axis of the replicated cell

c     n REALLY SHOULD be a multiple of 16 

      function imagevec3(pos,n,i)  result (imageout)
      integer, intent (in) ::i,n
      real(t_p),  intent (in) ::  pos(n)
      integer k
      real(t_p) numerator
      real(t_p) coordcell,coordcell2,mul
!DIR$ ATTRIBUTES ALIGN:64::imageout
      real(t_p) imageout(n)
c
      if (use_bounds) then
         SELECT CASE (i)
            CASE (1)
               coordcell  = xcell
               coordcell2 = xcell2
            CASE (2)
               coordcell  = ycell
               coordcell2 = ycell2
            CASE (3)
               coordcell  = zcell
               coordcell2 = zcell2
            CASE DEFAULT
               coordcell  = xcell
               coordcell2 = xcell2
         END SELECT
!DIR$ ASSUME (mod(n,16).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
         do k = 1, n
            imageout(k) =  pos(k)
     &                   - int( (abs(pos(k)) - coordcell2) / coordcell
     &                         + 1.0_ti_p
     &                        ) * sign (coordcell,pos(k))
         enddo
      endif
      end function imagevec3

c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine imagevec --  Use imagevec to compute minimal       ##
c     ##                          distance                              ##
c     ##                                                                ##
c     ####################################################################
c
      subroutine imagevec(pos,n,i)
      implicit none
      integer, intent (in) ::i,n
      real(t_p),  intent (inout) ::  pos(n)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      pos = imagevec3 (pos,n,i)
      return
      end subroutine imagevec

c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine image3dvec -- Use imagevec to compute minimal       ##
c     ##                          distance in all directions            ##
c     ##                                                                ##
c     ####################################################################
c
      subroutine image3dvec(posx,posy,posz,n)
      implicit none
      integer, intent (in) ::n
      real(t_p), intent (inout) ::  posx(n)
      real(t_p), intent (inout) ::  posy(n)
      real(t_p), intent (inout) ::  posz(n)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD
      posx = imagevec3 (posx,n,1)
!DIR$ SIMD
      posy = imagevec3 (posy,n,2)
!DIR$ SIMD
      posz = imagevec3 (posz,n,3)
      return
      end subroutine image3dvec

c     #################################################################
c     ##                                                             ##
c     ##  function setscale2  --  compute the scaling factors        ##
c     ##                                                             ##
c     #################################################################
c
c
c     "setscale2" takes the indice of the reference atom, the global indices
c     of the neighbours, the scaletype and returns the scaling factors

      function setscale2(indice,kvec,nmax,scaletype) result(coeffvec)

      implicit none
      integer, parameter :: maxind= 64
      integer    , intent (IN) :: indice,nmax
      integer    , intent (IN) :: kvec(nmax)
      character*6, intent (IN) :: scaletype

      integer i,idefault,j,k,ldefault
      integer j2,j3,j4,j5
      integer imax,indmax
      integer nn12,nn13,nn14,nn15,nnp11
!DIR$ ATTRIBUTES ALIGN:64::tmp11
      integer tmp11(maxind)
!DIR$ ATTRIBUTES ALIGN:64::tmp12
      integer tmp12(maxind)
!DIR$ ATTRIBUTES ALIGN:64::tmp13
      integer tmp13(maxind)
!DIR$ ATTRIBUTES ALIGN:64::tmp14
      integer tmp14(maxind)
!DIR$ ATTRIBUTES ALIGN:64::tmp15
      integer tmp15(maxind)

      real(t_p) scale2,scale3,scale4,scale41,scale5
!DIR$ ATTRIBUTES ALIGN:64::coeffvec
      real(t_p) coeffvec(nmax)
!DIR$ ATTRIBUTES ALIGN:64::scale41vec
      real(t_p) scale41vec(nmax)

cdeb  if(rank.eq.0.and.tinkerdebug) write(*,*) 'setscale2'

      idefault = nloop + 1 !default value is natoms + 1
                           !could not be a real indice
      ldefault = nloop + 2 !default value is natoms + 2
                           !p41scale could not be taken
c
c    get number of atoms  directly (1-2) or 1-3 or 1-4 bonded
c
      nnp11 = np11(indice)
      nn12  = n12 (indice)
      nn13  = n13 (indice)
      nn14  = n14 (indice)
      nn15  = n15 (indice)
     
c
c     imax is the biggest indice. Reduces the number of loops
c     to the minimum possible
      imax   = max(nnp11,nn12,nn13,nn14,nn15)
      indmax = (int(imax / 4 ) + 1) * 4 ! First multiple of 4 
      indmax = merge(imax,indmax , mod(imax,4 ).eq.0)
c
c  initialize scaling coefficients for all tomes and  for additional factor
c  for 1-4 intragroup polarization
c
!DIR$ ASSUME (MOD(nmax,16).eq.0)
!DIR$ VECTOR ALIGNED
      do k = 1, nmax
         coeffvec  (k) = 1.0_ti_p
         scale41vec(k) = 1.0_ti_p
      enddo
c     default value for indices in the list of bounded atoms
c     indice idefault cannot be of a real atom
c     indice ldefault cannot be of a real atom

!DIR$ ASSUME (MOD(indmax,4).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX=4
      do k = 1,indmax
            tmp11(k) = ldefault ! default value
            tmp12(k) = idefault ! default value
            tmp13(k) = idefault ! default value
            tmp14(k) = idefault ! default value
            tmp15(k) = idefault ! default value
      enddo
c
c     get list of indices 
c
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX=4
!DIR$ NOFUSION
      do k = 1,indmax
         if(k.le.nnp11) tmp11(k) = ip11(k,indice)
         if(k.le.nn12)  tmp12(k) = i12 (k,indice)
         if(k.le.nn13)  tmp13(k) = i13 (k,indice)
         if(k.le.nn14)  tmp14(k) = i14 (k,indice)
         if(k.le.nn15)  tmp15(k) = i15 (k,indice)
      enddo
c
c     select scaling coefficients
c
      select case (scaletype)
             case('cscale')
                            scale2  = c2scale
                            scale3  = c3scale
                            scale4  = c4scale
                            scale41 = 1.0_ti_p
                            scale5  = c5scale
             case('mscale')
                            scale2  = m2scale
                            scale3  = m3scale
                            scale4  = m4scale
                            scale41 = 1.0_ti_p
                            scale5  = m5scale
             case('pscale')
                            scale2  = p2scale
                            scale3  = p3scale
                            scale4  = p4scale
                            scale41 = p41scale
                            scale5  = p5scale
c
c                       Find where additional factor for 1-4 intragroup
c                       polarization should be taken into account
c
!DIR$ ASSUME (MOD(indmax,4).eq.0)
!DIR$ SIMD
                            do k = 1, nmax
                               do i = 1,indmax
!DIR$ ASSUME (MOD(nmax,16).eq.0)
                                  if(kvec(k) == tmp11(i)) then
                                     scale41vec(k) = scale41
                                     cycle
                                  endif
                               enddo
                            enddo
             case('vscale')
                            scale2  = v2scale
                            scale3  = v3scale
                            scale4  = v4scale
                            scale41 = 1.0_ti_p
                            scale5  = v5scale
             case default
                            scale2  = 1.0_ti_p
                            scale3  = 1.0_ti_p
                            scale4  = 1.0_ti_p
                            scale41 = 1.0_ti_p
                            scale5  = 1.0_ti_p
      endselect
c
c     for all atoms in the list, give the correct scale factor
c
!DIR$ ASSUME (MOD(indmax,4).eq.0)
!DIR$ ASSUME (MOD(nmax,16).eq.0)
!DIR$ SIMD
      do k = 1, nmax
         do i = 1, indmax
         j2 = tmp12(i)
         j3 = tmp13(i)
         j4 = tmp14(i)
         j5 = tmp15(i)
            if (kvec(k) == j2) coeffvec(k) = scale2
            if (kvec(k) == j3) coeffvec(k) = scale3
            if (kvec(k) == j4) coeffvec(k) = scale4 * scale41vec(k)
            if (kvec(k) == j5) coeffvec(k) = scale5
         enddo
      enddo
      return
      end function setscale2

c     #################################################################
c     ##                                                             ##
c     ##  function setscale2p  --  compute the scaling factors       ##
c     ##                                                             ##
c     #################################################################
c
c
c     "setscale2" takes the indice of the reference atom, the global indices
c     of the neighbours, the scaletype and returns the scaling factors

      function setscale2p(indice,kvec,nmax,scaletype) result(coeffvec)

      implicit none
      integer, parameter :: maxind = 64
      integer    , intent (IN) :: indice,nmax
      integer    , intent (IN) :: kvec(nmax)
      character*6, intent (IN) :: scaletype

      integer i,idefault,k
      integer j1,j2,j3,j4
      integer imax,indmax
      integer nnp11,nnp12,nnp13,nnp14
!DIR$ ATTRIBUTES ALIGN:64::tmp11
      integer tmp11(maxind)
!DIR$ ATTRIBUTES ALIGN:64::tmp12
      integer tmp12(maxind)
!DIR$ ATTRIBUTES ALIGN:64::tmp13
      integer tmp13(maxind)
!DIR$ ATTRIBUTES ALIGN:64::tmp14
      integer tmp14(maxind)

      real(t_p) scale1,scale2,scale3,scale4
!DIR$ ATTRIBUTES ALIGN:64::coeffvec
      real(t_p) coeffvec(nmax)
c
c     default value for indices in the list of bounded
c     indice idefault cannot be of a real atom

      idefault = nloop + 1 !default value is natoms + 1
                           !could not be a real indice
c
c    get number of atoms  directly (1-2) or 1-3 or 1-4 bonded
c
      nnp11 = np11(indice)
      nnp12 = np12(indice)
      nnp13 = np13(indice)
      nnp14 = np14(indice)
c
c     imax is the biggest indice. Reduces the number of loops
c     to the minimum possible
      imax   = max(nnp11,nnp12,nnp13,nnp14)
      indmax = (int(imax / 4 ) + 1) * 4 ! First multiple of 4 
      indmax = merge(imax,indmax , mod(imax,4 ).eq.0)
c
c  Initialize scaling coefficients
c
!DIR$ ASSUME_ALIGNED coeffvec:64
!DIR$ ASSUME (MOD(nmax,16).eq.0)
      do k = 1, nmax
         coeffvec(k) = 1.0_ti_p

      enddo
c     default value for indices in the list of bounded atoms
c     indice idefault cannot be of a real atom
c     indice ldefault cannot be of a real atom

!DIR$ ASSUME (MOD(indmax,4).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX=4
      do k = 1, indmax
         tmp11(k) = idefault ! default value
         tmp12(k) = idefault ! default value
         tmp13(k) = idefault ! default value
         tmp14(k) = idefault ! default value

      enddo
c
c     get list of indices 
c
!DIR$ ASSUME (MOD(indmax,4).eq.0)
!DIR$ VECTOR ALIGNED
!DIR$ LOOP COUNT MAX=4
!DIR$ NOFUSION
      do k = 1,indmax
         if(k.le.nnp11) tmp11(k) = ip11(k,indice)
         if(k.le.nnp12) tmp12(k) = ip12(k,indice)
         if(k.le.nnp13) tmp13(k) = ip13(k,indice)
         if(k.le.nnp14) tmp14(k) = ip14(k,indice)
      enddo
c
c     select scaling coefficients
c
      select case (scaletype)
             case('dscale')
                            scale1 = d1scale
                            scale2 = d2scale
                            scale3 = d3scale
                            scale4 = d4scale
             case('uscale')
                            scale1 = u1scale
                            scale2 = u2scale
                            scale3 = u3scale
                            scale4 = u4scale
             case default
                            scale1 = 1.0_ti_p
                            scale2 = 1.0_ti_p
                            scale3 = 1.0_ti_p
                            scale4 = 1.0_ti_p
      endselect
c
c     for all atoms in the list, give the correct scale factor
c
!DIR$ ASSUME (MOD(indmax,4).eq.0)
!DIR$ ASSUME (MOD(nmax,16).eq.0)
!DIR$ SIMD
      do k = 1, nmax
!DIR$ LOOP COUNT MAX=4
         do i = 1, indmax
            if (kvec(k) == tmp11(i)) coeffvec(k) = scale1
            if (kvec(k) == tmp12(i)) coeffvec(k) = scale2
            if (kvec(k) == tmp13(i)) coeffvec(k) = scale3
            if (kvec(k) == tmp14(i)) coeffvec(k) = scale4
         enddo
      enddo
      return
      end function setscale2p
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine setscale -- set interaction or exclusion scaling   ##
c     ##                         coefficients for connected atoms       ##
c     ##                                                                ##
c     ####################################################################
c

      subroutine setscale(indice,kvec,nmax,scaletype,coeffvec)
      implicit none
      integer, intent(IN) :: indice,nmax
      integer, intent(IN) :: kvec(nmax)
c!DIR$ ATTRIBUTES ALIGN:64::coeffvec
      real(t_p), intent (out)::  coeffvec(nmax)
      character*6 , intent(IN):: scaletype

      coeffvec = setscale2(indice,kvec,nmax,scaletype)
      return
      end subroutine setscale
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine setscalep -- set interaction or exclusion scaling  ##
c     ##                         coefficients for connected atoms       ##
c     ##                                                                ##
c     ####################################################################
c

      subroutine setscalep(indice,kvec,nmax,scaletype,coeffvec)
      implicit none
      integer, intent(IN) :: indice,nmax
      integer, intent(IN) :: kvec(nmax)
c!DIR$ ATTRIBUTES ALIGN:64::coeffvec
      real(t_p), intent (out)::  coeffvec(nmax)
      character*6 , intent(IN):: scaletype

      coeffvec = setscale2p(indice,kvec,nmax,scaletype)
      return
      end subroutine setscalep

      end module utilvec
