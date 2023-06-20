c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine bounds  --  check periodic boundary conditions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "bounds" finds the center of mass of each molecule and
c     translates any stray molecules back into the periodic box
c     atoms involved in a CV (through colvars or plumed) are dealt with
c     separately
c
c
      subroutine bounds
      use sizes
      use atmtyp
      use atoms
      use boxes
      use molcul
#ifdef COLVARS
      use colvars
#endif
#ifdef PLUMED
      use plumed
#endif
      implicit none
      integer i,j,k,l,lglob
      integer init,stop
      integer nlist
      integer, allocatable :: list(:)
c
      if (allocated(list)) deallocate(list)
      allocate (list(n))
c
#ifdef COLVARS
      if (use_colvars) then

      if (ncvatomsmol.gt.0) then
        nlist = ncvatomsmol
        do j = 1, ncvatomsmol
          list(j) = cvatomsmol(j)
        end do
        call boundslist(nlist,list)
      end if
c
c     then the other ones, molecule by molecule
c
        do i = 1, nmol
          nlist = 0
          list = 0
          init = imol(1,i)
          stop = imol(2,i)
          do j = init, stop
            k = kmol(j)
            do l = 1, ncvatomsmol
              lglob = cvatomsmol(l)
              if (lglob.eq.k) then
                goto 20
              end if
            end do
            nlist = nlist + 1
            list(nlist) = k
 20         continue
          end do
          call boundslist(nlist,list)
        end do
      else
#endif
#ifdef PLUMED
      if (lplumed) then
c
c      with PLUMED we don't know which atoms are in the cv so we wrap
c      the whole system at once
c
        nlist = n
        do i = 1, n
          list(i) = i
        end do
        call boundslist(nlist,list)
      else
#endif
c
c     wrap by molecules
c
        do i = 1, nmol
           nlist = 0
           list = 0
           init = imol(1,i)
           stop = imol(2,i)
           do j = init, stop
              k = kmol(j)
              nlist = nlist + 1
              list(nlist) = k
           end do
           call boundslist(nlist,list)
        end do
#ifdef PLUMED
      end if
#endif
#ifdef COLVARS
      end if
#endif
      deallocate (list)

      return
      end
c
c     "boundslist" finds the center of mass of a list of atoms
c     translates any stray atoms back into the periodic box
c
c
      subroutine boundslist(nlist,list)
      use sizes
      use atmtyp
      use atoms
      use boxes
      implicit none
      integer i,j,k
      integer nlist,list(nlist)
      real*8 weigh,weight
      real*8 xmid,ymid,zmid
      real*8 xfrac,yfrac,zfrac
      real*8 xcom,ycom,zcom
c
c
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      weight = 0.0d0
      do j = 1, nlist
         k = list(j)
         weigh = mass(k)
         weight = weight + weigh
         xmid = xmid + x(k)*weigh
         ymid = ymid + y(k)*weigh
         zmid = zmid + z(k)*weigh
      end do
      xmid = xmid / weight
      ymid = ymid / weight
      zmid = zmid / weight
c
c     get fractional coordinates of center of mass
c
      if (orthogonal .or. octahedron) then
         zfrac = zmid
         yfrac = ymid
         xfrac = xmid
      else if (monoclinic) then
         zfrac = zmid / beta_sin
         yfrac = ymid
         xfrac = xmid - zfrac*beta_cos
      else if (triclinic) then
         zfrac = zmid / gamma_term
         yfrac = (ymid - zfrac*beta_term) / gamma_sin
         xfrac = xmid - yfrac*gamma_cos - zfrac*beta_cos
      end if
c
c     translate center of mass into the periodic box
c
      do while (xfrac .gt. xbox2)
         xfrac = xfrac - xbox
      end do
      do while (xfrac .lt. -xbox2)
         xfrac = xfrac + xbox
      end do
      do while (yfrac .gt. ybox2)
         yfrac = yfrac - ybox
      end do
      do while (yfrac .lt. -ybox2)
         yfrac = yfrac + ybox
      end do
      do while (zfrac .gt. zbox2)
         zfrac = zfrac - zbox
      end do
      do while (zfrac .lt. -zbox2)
         zfrac = zfrac + zbox
      end do
c
c     truncated octahedron needs to have corners removed
c
      if (octahedron) then
         if (abs(xfrac)+abs(yfrac)+abs(zfrac) .gt. box34) then
            xfrac = xfrac - sign(xbox2,xfrac)
            yfrac = yfrac - sign(ybox2,yfrac)
            zfrac = zfrac - sign(zbox2,zfrac)
         end if
      end if
c
c     convert translated fraction center of mass to Cartesian
c
      if (orthogonal .or. octahedron) then
         xcom = xfrac
         ycom = yfrac
         zcom = zfrac
      else if (monoclinic) then
         xcom = xfrac + zfrac*beta_cos
         ycom = yfrac
         zcom = zfrac * beta_sin
      else if (triclinic) then
         xcom = xfrac + yfrac*gamma_cos + zfrac*beta_cos
         ycom = yfrac*gamma_sin + zfrac*beta_term
         zcom = zfrac * gamma_term
      end if
c
c     translate coordinates via offset from center of mass
c
      do j = 1, nlist
         k = list(j)
         x(k) = x(k) - xmid + xcom
         y(k) = y(k) - ymid + ycom
         z(k) = z(k) - zmid + zcom
      end do
      return
      end
