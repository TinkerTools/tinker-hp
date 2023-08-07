c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine prtdynbeads  --  output of MD restart information  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "prtdynbeads" writes out the information needed to restart a
c     molecular dynamics trajectory to an external disk file while
c     using PIMD
c
c
      subroutine prtdynbeads(polymer,ibead_beg,ibead_end)
      use atoms
      use atmtyp, only: mass
      use units, only: convert
      use beads
      use boxes
      use files
      use group
      use mdstuf
      use moldyn
      use titles
      use beads
      implicit none
      type(POLYMER_COMM_TYPE), intent(in) :: polymer
      integer, intent(in) :: ibead_beg,ibead_end
      integer i,idyn
      integer freeunit
      logical exist
      real*8 :: wt
      integer :: ibead
      character*2 atmc
      character*39 fstr
      character*3 numberbeads
      character*240 dynfile

      do ibead = ibead_beg,ibead_end 
c
c     update an existing restart file or open a new one
c
        write(numberbeads, '(i3.3)') ibead
        dynfile = filename(1:leng)//'_beads'//numberbeads//'.dyn'
        idyn = freeunit ()
        !write(*,*) dynfile
        inquire (file=dynfile,exist=exist)
        if (exist) then
           open (unit=idyn,file=dynfile,status='old')
           rewind (unit=idyn)
        else
           open (unit=idyn,file=dynfile,status='new')
        end if
c
c       save the number of atoms and the title string
c
        fstr = '('' Number of Atoms and Title :'')'
        write (idyn,fstr(1:32))
        atmc = 'i6'
        if (n .ge. 100000)  atmc = 'i7'
        if (n .ge. 1000000)  atmc = 'i8'
        if (ltitle .eq. 0) then
           fstr = '('//atmc//')'
           write (idyn,fstr(1:4))  n
        else
           fstr = '('//atmc//',2x,a)'
           write (idyn,fstr(1:9))  n,title(1:ltitle)
        end if
c
c       save the periodic box edge lengths and angles
c
        fstr = '('' Periodic Box Dimensions :'')'
        write (idyn,fstr(1:30))
        fstr = '(3d26.16)'
        write (idyn,fstr(1:9))  xbox,ybox,zbox
        write (idyn,fstr(1:9))  alpha,beta,gamma
c
c       save the atomic positions, velocities and accelerations
c
        fstr = '('' Current Atomic Positions :'')'
        write (idyn,fstr(1:31))
        fstr = '(3d26.16)'
        do i = 1, n
           write (idyn,fstr(1:9))  polymer%pos(:,i,ibead)
        end do
        fstr = '('' Current Atomic Velocities :'')'
        write (idyn,fstr(1:32))
        fstr = '(3d26.16)'
        do i = 1, n
           write (idyn,fstr(1:9)) polymer%vel(:,i,ibead) 
        end do
        fstr =  '('' Current Atomic Accelerations :'')'
        write (idyn,fstr(1:36))
        fstr = '(3d26.16)'
        do i = 1, n
           wt=convert/mass(i)
           write (idyn,fstr(1:9)) polymer%forces(:,i,ibead)*wt
        end do
        fstr =  '('' Alternate Atomic Accelerations :'')'
        write (idyn,fstr(1:38))
        fstr = '(3d26.16)'
        if(allocated(polymer%forces_slow)) then
          do i = 1, n
             wt=convert/mass(i)
             write (idyn,fstr(1:9)) polymer%forces_slow(:,i,ibead)*wt
          end do
        else
          do i = 1, n
             write (idyn,fstr(1:9)) 0.d0, 0.d0, 0.d0 
          end do
        endif
        fstr =  '('' pbc wrap index :'')'
        write (idyn,fstr(1:21))
        fstr = '(3i5)'
        do i = 1, n
           write (idyn,fstr(1:9))  pbcwrapindex(1,i),pbcwrapindex(2,i)
     &                                              ,pbcwrapindex(3,i)
        end do
c
c       close the dynamics trajectory restart file
c
        close (unit=idyn)
      
      enddo
      end
