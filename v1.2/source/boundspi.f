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
      subroutine boundspi(polymer,ibead_beg,ibead_end)
      use sizes
      use atmtyp
      use atoms
      use boxes
      use molcul
      use beads
#ifdef COLVARS
      use colvars
#endif
#ifdef PLUMED
      use plumed
#endif
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer, intent(in) :: ibead_beg,ibead_end
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
        call boundslistpi(polymer,ibead_beg,ibead_end,nlist,list)
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
          call boundslistpi(polymer,ibead_beg,ibead_end,nlist,list)
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
        call boundslistpi(polymer,ibead_beg,ibead_end,nlist,list)
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
           call boundslistpi(polymer,ibead_beg,ibead_end,nlist,list)
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
      subroutine boundslistpi(polymer,ibead_beg,ibead_end,nlist,list)
      use sizes
      use atmtyp
      use atoms
      use boxes
      use beads
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer, intent(in) :: ibead_beg,ibead_end
      integer, intent(in) :: nlist,list(nlist)
      integer j,k,ibead
      real*8 weigh,weight
      real*8 xmid,ymid,zmid
      real*8 dx,dy,dz
      integer :: xshift,yshift,zshift
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
         xmid = xmid + polymer%eigpos(1,k,1)*weigh
         ymid = ymid + polymer%eigpos(2,k,1)*weigh
         zmid = zmid + polymer%eigpos(3,k,1)*weigh
      end do
      xmid = xmid / weight
      ymid = ymid / weight
      zmid = zmid / weight

      call compute_wrap_shifts(xmid,ymid,zmid,dx,dy,dz
     &                         ,xshift,yshift,zshift)
c
c     translate coordinates via offset from center of mass
c
      do j = 1, nlist
         k = list(j)
         polymer%eigpos(1,k,1) = polymer%eigpos(1,k,1) + dx
         polymer%eigpos(2,k,1) = polymer%eigpos(2,k,1) + dy
         polymer%eigpos(3,k,1) = polymer%eigpos(3,k,1) + dz
         pbcwrapindex(1,k) = pbcwrapindex(1,k) + xshift
         pbcwrapindex(2,k) = pbcwrapindex(2,k) + yshift
         pbcwrapindex(3,k) = pbcwrapindex(3,k) + zshift
      end do
      do ibead = ibead_beg,ibead_end 
        do j=1,nlist
          k=list(j)
          polymer%pos(1,k,ibead) = polymer%pos(1,k,ibead) + dx
          polymer%pos(2,k,ibead) = polymer%pos(2,k,ibead) + dy
          polymer%pos(3,k,ibead) = polymer%pos(3,k,ibead) + dz
        enddo
      enddo
      return
      end
