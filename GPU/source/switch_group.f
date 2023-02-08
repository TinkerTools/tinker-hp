c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine switch_group_ene  --  group polarization energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "switch_group_ene" calculates the group induced dipole polarization energy:
c      if A and B are two groups and their interactions are scaled by alpha,
c      then the Epol(A-B) = alpha*(Epol(A+B) - Epol(A) - Epol(B)) st:
c      Epol(A+B) = alpha*Epol(A+B) + (1-alpha)*(Epol(A) + Epol(B))
c
c
#include "tinker_macro.h"
      subroutine switch_group(dograd)
      use atoms
      use domdec
      use energi
      use deriv
      use group
      use inform
      use iounit
      use mpole
      use mpi
      use potent
      use utils    ,only: set_to_zero1
      use tinheader,only: ti_p
      use utilgpu
      use virial
      implicit none
      logical, intent(in) :: dograd
      integer i,igroup,jgroup,iglob,iipole,iipolegroup,ilocgroup
      integer j,iglobgroup, iloc,igroup2,ibeggroup
      real(r_p) epoltot,ww
      real(r_p) :: epol_group(2)
      real(t_p) :: walpha, wbeta, wgamma
      integer, allocatable, save :: kgrp_bis(:)
      integer, save :: natgroups(2)
      logical skipscale, compute_group(2),do_polar,do_mpole
      logical, save :: is_bipartition
      logical, save :: first_in=.true.
      interface
       subroutine rot_mat_site(recBuild,nk,poleglobvec)
         logical,intent(in)::recBuild
         integer,intent(in)::nk
         integer,intent(in)::poleglobvec(nk)
       end subroutine
       function check_bipartition()
        logical :: check_bipartition
       end function
      end interface
c
c     return if more than 2 groups
c
      do_polar = use_group_polar .and. use_polar
      do_mpole = use_group_mpole .and. use_mpole
      if(.not. ANY([do_polar,do_mpole])) return

      if(first_in) then
        is_bipartition = check_bipartition()
      endif
      if (.not. is_bipartition) then
          if (rank.eq.0) then
            write(iout,*) 'Error: polar/mpole group scaling '
     $        ,' only possible with  a bipartition' 
     $        ,' (same wgrp for all groups but group 1)'
            call fatal
          end if
        return
      end if

      
c
c     return if no interactions between groups is scaled
c
      walpha = wgrp(1,1)
      wbeta = wgrp(2,2)
      wgamma = wgrp(1,2)
      if (wgamma.eq.1.0_ti_p .and. wbeta.eq.1.0_ti_p
     &      .and.walpha.eq.1.0_ti_p) return
      compute_group(1) = (wgamma-walpha).ne.0.0_ti_p
      compute_group(2) = (wgamma-wbeta).ne.0.0_ti_p

c      write(*,*) walpha.ne.1.0_ti_p, wbeta.ne.1.0_ti_p
c     & , wgamma.ne.1.0_ti_p, compute_group
      if(first_in) then
!$acc enter data create(epgroup,emgroup,vir_group) async
        first_in=.false.

        ! REGROUP ALL GROUPS > 1 WITH GROUP 0 (FOR BIPARTITION)
        ! AND REORDER kgrp INTO kgrp_bis
        allocate(kgrp_bis(n))
        natgroups(1) = igrp(2,0) - igrp(1,0) + 1
        do i=1,natgroups(1)
          kgrp_bis(i) = kgrp(igrp(1,0)+i-1)
        enddo
        if(ngrp>1) then
          do igroup=3,ngrp+1
            natgroup = igrp(2,igroup-1) - igrp(1,igroup-1) + 1
            do i=1,natgroup
              kgrp_bis(natgroups(1)+i) = kgrp(igrp(1,igroup-1)+i-1)
            enddo
            natgroups(1) = natgroups(1) + natgroup
          enddo
        endif
        natgroups(2) = igrp(2,1) - igrp(1,1) + 1
        do i=1,natgroups(2)
          kgrp_bis(natgroups(1)+i) = kgrp(igrp(1,1)+i-1)
        enddo
!$acc enter data copyin(kgrp_bis) async

      endif
c
c     save total unscaled polarization energy
c
!$acc serial async present(ep,em,wgrp)
      if(do_polar) then
        ep = real(wgrp(1,2),r_p)*ep
      endif
      if(do_mpole) then
        em = real(wgrp(1,2),r_p)*em
      endif
!$acc end serial

      if(dograd .and. (wgamma.ne.1.0_ti_p)) then
        if(do_polar) then
!$acc parallel loop collapse(2) async present(dep,wgrp)
          do i=1,nbloc; do j=1,3
            dep(j,i) = real(wgrp(1,2),r_p)*dep(j,i)
          enddo; enddo
!$acc parallel loop collapse(2) async present(deprec,wgrp)
          do i=1,nlocrec2; do j=1,3
            deprec(j,i) = real(wgrp(1,2),r_p)*deprec(j,i)
          enddo; enddo
        endif
        if(do_mpole) then
!$acc parallel loop collapse(2) async present(dem,wgrp)
          do i=1,nbloc; do j=1,3
            dem(j,i) = real(wgrp(1,2),r_p)*dem(j,i)
          enddo; enddo
!$acc parallel loop collapse(2) async present(demrec,wgrp)
          do i=1,nlocrec2; do j=1,3
            demrec(j,i) = real(wgrp(1,2),r_p)*demrec(j,i)
          enddo; enddo
        endif
      end if
c
c     loop on the groups and add their energy 
c
      do igroup = 1, 2
        if(.not.compute_group(igroup)) CYCLE
        !natgroup = igrp(2,igroup-1) - igrp(1,igroup-1) + 1
        natgroup = natgroups(igroup)
        if(natgroup<=0) cycle
        if(deb_path) 
     &     write(*,*) "compute group",igroup-1
        call prmem_request(globglobgroup,natgroup,async=.TRUE.)
        call prmem_request(loclocgroup,n,async=.TRUE.)
        call prmem_request(globgroup,natgroup,async=.TRUE.)
        call prmem_request(locgroup,natgroup,async=.TRUE.)
        call prmem_request(poleglobgroup,natgroup,async=.TRUE.)
        call prmem_request(globpolegroup,natgroup,async=.TRUE.)
        call prmem_request(polelocgroup,natgroup,async=.TRUE.)
        call prmem_request(ipolegroup,natgroup,async=.TRUE.)
        call prmem_request(pollistgroup,natgroup,async=.TRUE.)

        if (allocated(domlengroup)) deallocate (domlengroup)
        allocate (domlengroup(nproc))
        if (allocated(bufbeggroup)) deallocate (bufbeggroup)
        allocate (bufbeggroup(nproc))
        if (allocated(domlenpolegroup)) deallocate (domlenpolegroup)
        allocate (domlenpolegroup(nproc))
        if (allocated(bufbegpolegroup)) deallocate (bufbegpolegroup)
        allocate (bufbegpolegroup(nproc))

!$acc serial async present(vir_group)
        vir_group(1,1)=0.0_re_p
        vir_group(2,1)=0.0_re_p
        vir_group(3,1)=0.0_re_p
        vir_group(1,2)=0.0_re_p
        vir_group(2,2)=0.0_re_p
        vir_group(3,2)=0.0_re_p
        vir_group(1,3)=0.0_re_p
        vir_group(2,3)=0.0_re_p
        vir_group(3,3)=0.0_re_p
!$acc end serial

!$acc wait
        npolegroup = 0
        nlocatgroup = 0
        ibeggroup = (igroup-1)*natgroups(1)+1
!$acc parallel loop async present(kgrp,igrp,globglobgroup,loclocgroup
!$acc&  ,pollist,pollistgroup,ipolegroup,repart,globgroup,locgroup
!$acc&  ,globpolegroup,repart) 
!$acc&  copy(nlocatgroup,npolegroup)
        do i = 1, natgroup
          !iglob = kgrp(igrp(1,igroup-1)+i-1)
          iglob = kgrp_bis(ibeggroup+i-1)
          globglobgroup(i) = iglob
          loclocgroup(iglob) = i

          iipole = pollist(iglob)
          if (iipole/=0) then
!$acc atomic capture
            npolegroup = npolegroup + 1
            iipolegroup = npolegroup
!$acc end atomic
            pollistgroup(i) = iipolegroup
            ipolegroup(iipolegroup) = i
            globpolegroup(iipolegroup) = iipole
          endif

          if (repart(iglob).eq.rank) then
!$acc atomic capture
            nlocatgroup = nlocatgroup + 1
            ilocgroup = nlocatgroup
!$acc end atomic
            globgroup(ilocgroup) = i
            locgroup(i) = ilocgroup
          endif

        end do

        npolelocgroup = 0
!$acc wait

!$acc parallel loop async present(ipolegroup,polelocgroup,poleglobgroup
!$acc&  ,globglobgroup,repart)
!$acc&  copy(npolelocgroup)
        do i = 1, npolegroup
          iglobgroup = ipolegroup(i)
          iglob = globglobgroup(iglobgroup)
          if (repart(iglob).eq.rank) then
!$acc atomic capture
            npolelocgroup = npolelocgroup + 1
            iipolegroup = npolelocgroup
!$acc end atomic
            poleglobgroup(iipolegroup) = i
            polelocgroup(i) = iipolegroup
          end if
        end do
!$acc wait
        call orderbuffergroup

        if(do_mpole) call mpole_scaling_group
        if(do_polar) call polar_scaling1_group

        !if(do_polar .or. do_mpole) then
        !  call chkpolegpu_group
        !  call rot_mat_site(.false.,npolegroup,globpolegroup)
        !endif

        if(dograd) then
          if(do_mpole) call empole1_group
          if(do_polar) call epolar1_group
        else
          if(do_mpole) call empole3_group
          if(do_polar) call epolar3_group
        endif
c
c     add scaled energies and forces
c
          ww = real(wgrp(igroup,igroup),r_p) - real(wgrp(1,2),r_p)
!$acc serial async present(ep,epgroup,em,emgroup,wgrp)
        if(do_polar) then
          ep = ep + ww*epgroup
        endif
        if(do_mpole) then
          em = em + ww*emgroup
        endif
!$acc end serial

        if(dograd) then
c!$acc update host(demgroup,locgroup,loclocgroup,emgroup) async
c!$acc wait
c          i=locgroup(loclocgroup(1))
c          j=locgroup(loclocgroup(12))
c          write(*,'(f15.6,A,2f15.6)') emgroup,"  demgroup"
c     &       ,demgroup(1,i),demgroup(3,j)

!$acc parallel loop async collapse(2) 
!$acc& present(dep,depgroup,wgrp,globgroup,loc,globglobgroup
!$acc&        ,dem,demgroup) 
          do i = 1, nlocatgroup; do j=1,3
            iglobgroup = globgroup(i)
            iglob = globglobgroup(iglobgroup)
            iloc = loc(iglob)
            if(do_polar) then
              dep(j,iloc) = dep(j,iloc) 
     &          + ww*depgroup(j,i)
            endif
            if(do_mpole) then
              dem(j,iloc) = dem(j,iloc)
     &          + ww*demgroup(j,i)
            endif
          end do; enddo
          
          if(use_virial) then
!$acc serial async present(vir,vir_group)
            vir(1,1) = vir(1,1) + ww*vir_group(1,1)
            vir(2,1) = vir(2,1) + ww*vir_group(2,1)
            vir(3,1) = vir(3,1) + ww*vir_group(3,1)
            vir(1,2) = vir(1,2) + ww*vir_group(1,2)
            vir(2,2) = vir(2,2) + ww*vir_group(2,2)
            vir(3,2) = vir(3,2) + ww*vir_group(3,2)
            vir(1,3) = vir(1,3) + ww*vir_group(1,3)
            vir(2,3) = vir(2,3) + ww*vir_group(2,3)
            vir(3,3) = vir(3,3) + ww*vir_group(3,3)
!$acc end serial
          endif

        endif

      end do

      return
      end

      function check_bipartition() result(is_bipartition)
        use group
        use tinheader
        implicit none
        logical :: is_bipartition
        real(t_p) :: walpha,wgamma
        integer :: i,j

        is_bipartition=.true.
        if(ngrp==1) return

        wgamma=wgrp(1,2)
        walpha=wgrp(1,1)

        do i=2,ngrp
          if(wgrp(i+1,i+1) .ne. walpha) then
            ! if self weight is not the same as group 0 (except for group 1)
            ! then the system is not bipartitioned
            is_bipartition=.false.
            return
          endif
          do j=2,ngrp
            if((wgrp(i+1,j+1) .ne. walpha)
     &       .or. (wgrp(j+1,i+1) .ne. walpha)) then
              ! if weight between groups (except group 1) is not alpha
              ! then the system is not bipartitioned
              is_bipartition=.false.
              return
            endif
          enddo
          if((wgrp(i+1,1) .ne. walpha) 
     &       .or. (wgrp(1,i+1) .ne. walpha)) then
            ! if weight between groups and group 0 is not alpha
            ! then the system is not bipartitioned
            is_bipartition=.false.
            return
          endif
          if((wgrp(i+1,2) .ne. wgamma)
     &       .or. (wgrp(2,i+1) .ne. wgamma)) then
            ! if weight between groups and group 1
            ! is not the same as group 0
            ! then the system is not bipartitioned
            is_bipartition=.false.
            return
          endif
        end do

      end function check_bipartition
