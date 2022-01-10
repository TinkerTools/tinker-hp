c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #######################################################################
c     ##                                                                   ##
c     ##  subroutine kdispersion  -- dispersion parameters assignment      ##
c     ##                                                                   ##
c     #######################################################################
c
c     krepulsion does stuff to complete
c
      subroutine kdispersion
      use sizes
      use atoms
      use keys
      use potent
      use cutoff
      use dispersion
      implicit none
      integer i,next,k
      logical header
      character*20 keyword
      character*120 record
      character*20 string

c
c     defaults for dispersion term (and parameters)
c                                                            
      use_dispersion = .false.
      use_displist = .false.
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
         if (keyword(1:10) .eq. 'DISPERSION ') then
            use_dispersion = .true.
            use_displist = .true.
         end if
      end do
c
      if (use_dispersion) then
c
c       allocate global arrays 
c   
        if (allocated(vdwdisp1)) deallocate (vdwdisp1)
        allocate (vdwdisp1(n))
        if (allocated(vdwdisp2)) deallocate (vdwdisp2)
        allocate (vdwdisp2(n))
        if (allocated(vdwdisp3)) deallocate (vdwdisp3)
        allocate (vdwdisp3(n))
c
c       fill the arrays
c
        do i = 1, n
          vdwdisp1(i) = sibfadisp(1,type(i))
          vdwdisp2(i) = sibfadisp(2,type(i))
          vdwdisp3(i) = sibfadisp(3,type(i))
        end do
      end if
      end
