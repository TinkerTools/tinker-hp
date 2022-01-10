c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  subroutine krepulsion  --  repulsion parameters assignment       ##
c     ##                                                                   ##
c     #######################################################################
c
      subroutine krepulsion
      use atoms
      use atmtyp
      use cutoff
      use keys
      use repulsion
      use potent
      implicit none
      integer i,j,next,k
      logical header
      character*20 keyword
      character*120 record
      character*20 string

c
c     defaults for repulsion term (and parameters)
c                                                            
      use_repulsion = .false.
      use_replist = .false.
c
c     search keywords for repulsion commands
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:9) .eq. 'REPULSION ') then
            use_repulsion = .true.
            use_replist = .true.
         end if
      end do
c
      if (use_repulsion) then
c
c       allocate global arrays 
c   
        if (allocated(vdwrep)) deallocate (vdwrep)
        allocate (vdwrep(n))
        if (allocated(gorbrep)) deallocate (gorbrep)
        allocate (gorbrep(n))
c
c       fill the arrays
c
        do i = 1, n
          vdwrep(i)  = sibfarep(type(i))
          gorbrep(i) = gorb(type(i)) 
        end do
      end if
      end 
