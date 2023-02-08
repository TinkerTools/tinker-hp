c
c     Sorbonne Université
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     Umbrella Sampling control module
c

#include "tinker_macro.h"
      module USampling
      use kgeoms
      implicit none
      logical:: US_enable=.false.
      integer:: cpt_wh=1,step=1,step_save=-1
      integer:: USwrite=huge(USwrite)/2
      real(t_p),parameter:: UStime_period=0.5
      character(len=30)  :: USfile, Rdfile
      real(r_p):: timestep=0.0
      real(t_p),allocatable::Us_save(:)
      real(t_p),allocatable::Rd_save(:)

      contains

      subroutine init_USampling(nproc,arg)
      implicit none
      integer,intent(in) :: nproc
      character(len=240) :: arg(0:20)
      integer s_size, iwrite, iunit, freeunit
      real(r_p) dtdump

10    format( "=======================================================",
     &      /,"        WARNING  !!!  USampling disabled  !!!          ",
     &      /,"======================================================="
     &      )
11    format(/,"Tinker frame writing period is lower than Sampling's  ",
     &       /,"Increase it or consider reducing your sampling period !"
     &      )
12    format(7x,"Not supposed to work in parallel yet ")
20    format("---- U Sampling saving file for Fred ----",
     &      /, " ==== Restrain group ",
     &      /,
     &      /,6x,'time',8x,'distance',6x,'E Restrain',
     &      /)
21    format("---- U Sampling saving file for Fred ----",
     &      /, " ==== Restrain group ",
     &      /,
     &      /,6x,'time',8x,'distance',6x,'E Restrain',5x,'grpID'
     &      /)
30    format("---- U Sampling saving file for Fred ----",
     &      /, " ==== Restrain distance ",
     &      /,
     &      /,6x,'time',8x,'distance',6x,'E Restrain',
     &      /)
31    format("---- U Sampling saving file for Fred ----",
     &      /, " ==== Restrain distance ",
     &      /,
     &      /,6x,'time',8x,'distance',6x,'E Restrain',5x,'dstID'
     &      /)

      read (arg(3),*) timestep
      read (arg(4),*) dtdump

      cpt_wh   = 1
      step     = 0
      timestep = 1d-3*timestep
      USwrite  = nint(UStime_period/timestep)
      iwrite   = nint(dtdump/timestep)
      if (iwrite < USwrite) then 
         US_enable = .false.
         write(*,10)
         write(*,11)
      end if
      if (nproc.ne.1) then
         write(*,10)
         write(*,12)
         US_enable = .false.
      end if

      s_size   = iwrite/USwrite +1

      iunit   = freeunit()
      if (ngfix.gt.0) then
         USfile   = trim(arg(1))//'.US'
         call version(USfile,'new')
         open (unit=iunit,file=USfile,status='new')
         if (ngfix.gt.1) then
            write(iunit,21)
         else if (ngfix.eq.1) then
            write(iunit,20)
         end if
         close (iunit)

         if (allocated(Us_save)) deallocate(Us_save)
         allocate(Us_save(4*ngfix*s_size))
!$acc enter data create(Us_save)

         !print*,"USfile       ",USfile
         !print*,"US size save ",s_size
         !print*,"US write     ", USwrite
      end if

      if (ndfix.gt.0) then
         Rdfile   = trim(arg(1))//'.usd'
         call version(Rdfile,'new')
         open (unit=iunit,file=Rdfile,status='new')
         if (ndfix.gt.1) then
            write(iunit,31)
         else if (ndfix.eq.1) then
            write(iunit,30)
         end if
         close (iunit)
         if (allocated(Rd_save)) deallocate(Rd_save)
         allocate(Rd_save(4*ndfix*s_size))
!$acc enter data create(Rd_save)
      end if

      end subroutine

      subroutine prtUS
      implicit none
      integer iunit,i,j,loc,freeunit
      logical existf

      iunit   = freeunit()

30    format(F10.3,F16.6,F16.4)
35    format(F10.3,F16.6,F16.4,3x,I6,A1)
40    format(1x,"Updating US group in file",5X,A)
50    format(1x,"Updating US distance in file",5X,A)

      if (ngfix.gt.0) then

      inquire (file=USfile,exist=existf)
      if (existf) then
         open (unit=iunit,file=USfile,status='old',position='append')
      else
         print*, "There is a problem in US "
         print*, USfile ," should exist in this folder "
         print*, "Opening new file to proceed !!!"
         open (unit=iunit,file=USfile,status='new')
      end if

      write(*,40) USfile
!$acc update host(US_save)

      if (ngfix.gt.1) then
      do i = 1,cpt_wh-1
         do j = 1,ngfix
           loc = (i-1)*ngfix*4 + (j-1)*4
           write(iunit,35) US_save(loc+2),US_save(loc+3),
     &                     US_save(loc+4),US_save(loc+1)
           US_save(loc+1) = 0.0
           US_save(loc+2) = 0.0
           US_save(loc+3) = 0.0
           US_save(loc+4) = 0.0
         end do
      end do
      else
      do i = 1,cpt_wh-1
         loc = (i-1)*4
         write(iunit,30) US_save(loc+2),US_save(loc+3),US_save(loc+4)
         US_save(loc+2) = 0.0
         US_save(loc+3) = 0.0
         US_save(loc+4) = 0.0
      end do
      end if

      close(iunit)
!$acc update device(US_save)

      end if

      if (ndfix.gt.0) then

      inquire (file=Rdfile,exist=existf)
      if (existf) then
         open (unit=iunit,file=Rdfile,status='old',position='append')
      else
         print*, "There is a problem in usd "
         print*, Rdfile ," should exist in this folder "
         print*, "Opening new file to proceed !!!"
         open (unit=iunit,file=Rdfile,status='new')
      end if

      write(*,50) Rdfile
!$acc update host(Rd_save)

      if (ndfix.gt.1) then
      do i = 1,cpt_wh-1
         do j = 1,ndfix
           loc = (i-1)*ndfix*4 + (j-1)*4
           write(iunit,35) Rd_save(loc+2),Rd_save(loc+3),
     &                     Rd_save(loc+4),Rd_save(loc+1)
           Rd_save(loc+1) = 0.0
           Rd_save(loc+2) = 0.0
           Rd_save(loc+3) = 0.0
           Rd_save(loc+4) = 0.0
         end do
      end do
      else
      do i = 1,cpt_wh-1
         loc = (i-1)*4
         write(iunit,30) Rd_save(loc+2),Rd_save(loc+3),Rd_save(loc+4)
         Rd_save(loc+2) = 0.0
         Rd_save(loc+3) = 0.0
         Rd_save(loc+4) = 0.0
      end do
      end if

      close(iunit)
!$acc update device(Rd_save)

      end if
      cpt_wh = 1
      end subroutine

      end module

