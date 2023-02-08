c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine kangtor  --  angle-torsion parameter assignment  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "kangtor" assigns parameters for angle-torsion interactions
c     and processes new or changed parameter values
c
c
#include "tinker_macro.h"
      subroutine kangtor(init)
      use angtor
      use atmlst
      use atmtyp
      use atoms
      use couple
      use domdec
      use inform
      use iounit
      use keys
      use kantor
      use potent
      use tinMemory
      use tors
      use utils
      use utilgpu,only: rec_queue
#ifdef _OPENACC
     &           ,rec_stream
      use thrust
#endif
      implicit none
      integer i,j,k,l,m,nat
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer iitors,angtorcount,nangtorloc1
      integer size,next,isys,nang_cap
      real(t_p) at1,at2,at3
      real(t_p) at4,at5,at6
      logical header,swap
      character*4 pa,pb,pc,pd
      character*4 zeros
      integer(8),parameter:: blank=-1
      integer(8) pt
      character*20 keyword
      character*240 record
      character*240 string
      logical init
c
c
c     process keywords containing angle-torsion parameters
c
      if (init) then

        if (deb_Path) print*, 'kangtor init'
        header = .true.
        zeros = '0000'
        isys  = 0

        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:8) .eq. 'ANGTORS ') then
              ia = 0
              ib = 0
              ic = 0
              id = 0
              at1 = 0.0d0
              at2 = 0.0d0
              at3 = 0.0d0
              at4 = 0.0d0
              at5 = 0.0d0
              at6 = 0.0d0
              string = record(next:240)
              read (string,*,err=10,end=10)  ia,ib,ic,id,at1,at2,
     &                                       at3,at4,at5,at6
   10         continue
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    write (iout,20)
   20               format (/,' Additional Angle-Torsion Parameters :',
     &                      //,32x,'1st Angle',20x,'2nd Angle',
     &                      /,5x,'Atom Classes',7x,'1-Fold',3x,'2-Fold',
     &                        3x,'3-Fold',5x,'1-Fold',3x,'2-Fold',
     &                        3x,'3-Fold'/)
                 end if
                 write (iout,30)  ia,ib,ic,id,at1,at2,at3,at4,at5,at6
   30            format (2x,4i4,3x,3f9.3,2x,3f9.3)
              end if
              if (ib .lt. ic) then
                 call front_convert_base(0,ia,ib,ic,id,pt)
                 swap = .false.
              else if (ic .lt. ib) then
                 call front_convert_base(0,id,ic,ib,ia,pt)
                 swap = .true.
              else if (ia .le. id) then
                 call front_convert_base(0,ia,ib,ic,id,pt)
                 swap = .false.
              else if (id .lt. ia) then
                 call front_convert_base(0,id,ic,ib,ia,pt)
                 swap = .true.
              end if
              do j = 1, maxnat
                 if (kat(j).eq.blank .or. kat(j).eq.pt) then
                    kat(j) = pt
                    if (swap) then
                       atcon(1,j) = at4
                       atcon(2,j) = at5
                       atcon(3,j) = at6
                       atcon(4,j) = at1
                       atcon(5,j) = at2
                       atcon(6,j) = at3
                    else
                       atcon(1,j) = at1
                       atcon(2,j) = at2
                       atcon(3,j) = at3
                       atcon(4,j) = at4
                       atcon(5,j) = at5
                       atcon(6,j) = at6
                    end if
                    goto 50
                 end if
              end do
              write (iout,40)
   40         format (/,' KANGTOR  --  Too many Angle-Torsion',
     &                   ' Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        nat = maxnat
        do i = maxnat, 1, -1
           if (kat(i) .eq. blank)  nat = i - 1
        end do
c
c       assign the angle-torsion parameters for each torsion
c
        nangtor = 0
        if (nat .ne. 0) then
           do i = 1, ntors
              ia = itors(1,i)
              ib = itors(2,i)
              ic = itors(3,i)
              id = itors(4,i)
              ita = class(ia)
              itb = class(ib)
              itc = class(ic)
              itd = class(id)
              nbangtor(i) = nangtor
              if (itb .lt. itc) then
                 call front_convert_base(0,ita,itb,itc,itd,pt)
                 swap = .false.
              else if (itc .lt. itb) then
                 call front_convert_base(0,itd,itc,itb,ita,pt)
                 swap = .true.
              else if (ita .le. itd) then
                 call front_convert_base(0,ita,itb,itc,itd,pt)
                 swap = .false.
              else if (itd .lt. ita) then
                 call front_convert_base(0,itd,itc,itb,ita,pt)
                 swap = .true.
              end if
              do j = 1, nat
                 if (kat(j) .eq. pt) then
                    nangtor = nangtor + 1
                    if (swap) then
                       kant(1,nangtor) = atcon(4,j)
                       kant(2,nangtor) = atcon(5,j)
                       kant(3,nangtor) = atcon(6,j)
                       kant(4,nangtor) = atcon(1,j)
                       kant(5,nangtor) = atcon(2,j)
                       kant(6,nangtor) = atcon(3,j)
                    else
                       kant(1,nangtor) = atcon(1,j)
                       kant(2,nangtor) = atcon(2,j)
                       kant(3,nangtor) = atcon(3,j)
                       kant(4,nangtor) = atcon(4,j)
                       kant(5,nangtor) = atcon(5,j)
                       kant(6,nangtor) = atcon(6,j)
                    end if
                    ! Build kat_sys
                    if (.not.is_find8(kat_sys(1),isys,pt)) then
                       isys = isys+1
                       kat_sys(isys) = pt
                    end if
                    iat(1,nangtor) = i
                    m = 0
                    do k = 1, n12(ib)-1
                       do l = k+1, n12(ib)
                         m = m + 1
                         if ((i12(k,ib).eq.ia .and. i12(l,ib).eq.ic).or.
     &                    (i12(k,ib).eq.ic .and. i12(l,ib).eq.ia)) then
                            iat(2,nangtor) = anglist(m,ib)
                            goto 60
                         end if
                       end do
                    end do
   60               continue
                    m = 0
                    do k = 1, n12(ic)-1
                       do l = k+1, n12(ic)
                         m = m + 1
                         if ((i12(k,ic).eq.ib .and. i12(l,ic).eq.id).or.
     &                   (i12(k,ic).eq.id .and. i12(l,ic).eq.ib)) then
                            iat(3,nangtor) = anglist(m,ic)
                            goto 70
                         end if
                       end do
                    end do
                 end if
              end do
   70         continue
           end do
           kat_sys(0) = isys
        end if
c
c       turn off the angle-torsion potential if it is not used
c
        if (nangtor.eq.0)  use_angtor = .false.

        if (deb_Path) print*, 'nat',nat,isys
        if (use_angtor) then
           call upload_device_kangtor
        else
           call delete_data_kangtor
           return
        end if
      end if

      nat = kat_sys(0)
      call prmem_request(angtorglob,ntorsloc,async=.true.)

!$acc data present(nangtorloc)
!$acc serial async
      nangtorloc = 0
!$acc end serial
!$acc parallel loop async default(present)
      do i = 1, ntorsloc
         iitors = torsglob(i)
         ia = itors(1,iitors)
         ib = itors(2,iitors)
         ic = itors(3,iitors)
         id = itors(4,iitors)
         angtorcount = nbangtor(iitors)
         nangtorloc1 = 0
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         if (itb .lt. itc) then
            call front_convert_base(0,ita,itb,itc,itd,pt)
         else if (itc .lt. itb) then
            call front_convert_base(0,itd,itc,itb,ita,pt)
         else if (ita .le. itd) then
            call front_convert_base(0,ita,itb,itc,itd,pt)
         else if (itd .lt. ita) then
            call front_convert_base(0,itd,itc,itb,ita,pt)
         end if
!$acc loop seq
         do j = 1, nat
            if (kat_sys(j) .eq. pt) then
!$acc atomic capture
               nangtorloc  = nangtorloc + 1
               nang_cap    = nangtorloc
!$acc end atomic
               nangtorloc1 = nangtorloc1 + 1
               angtorglob(nang_cap) = angtorcount + nangtorloc1
               exit
            end if
         end do
      end do
!$acc update host(nangtorloc) async
!$acc end data

      end

      subroutine upload_device_kangtor
      use angtor
      use domdec,only: rank,hostcomm
      use inform,only: deb_Path
      use kantor
      use mpi   ,only: MPI_BARRIER
      implicit none
#ifdef _OPENACC
      integer ierr
 12   format(2x,'upload_device_kangtor')
      if(deb_Path) print 12
      call MPI_BARRIER(hostcomm,ierr)
!$acc enter data copyin(kat_sys,nangtorloc)
      ! Update if second call to routine
!$acc update device(kat_sys,nangtorloc)

!$acc update device(iat,nbangtor,kant)
#endif
      end subroutine

      subroutine delete_data_kangtor
      use angtor
      use inform,only: deb_Path
      use kantor
      use tinMemory
      implicit none
#ifdef _OPENACC
 12   format(2x,'delete_data_kangtor')
      if(deb_Path) print 12

!$acc exit data delete(nangtorloc,kat_sys)
#endif
      end subroutine
