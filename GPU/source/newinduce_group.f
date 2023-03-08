c 
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c      
c     evaluate induced dipole moments and the polarization energy
c     of a group of atom without PBC
c     using a (preconditioned) conjugate gradient algorithm 
c     or a JI/DIIS solver
c     useful to define the polarization to be retrieved from the total
c     one in "group-scaled" interactions
c      
c     literature reference:
c     "Scalable Evaluation of Polarization Energy and Associated Forces
c     in Polarizable Molecular Dynamics: II. Toward Massively Parallel
c     Computations Using Smooth Particle Mesh Ewald",L. Lagardere et al.,
c     J. Chem. Theory Comput., 2015, 11 (6), pp 2589–2599
c
#include "tinker_macro.h"
      module newinduce_group_inl
      use tintypes
      contains
#include "pair_efld.inc.f"
#include "pair_tmatxb.f.inc"
      end module

      subroutine newinduce_group
      use atoms     ,only: n,n_pe
      use atomsmirror
      use atmlst
      use domdec
      use deriv
      use ewald
      use iounit
      use inform    ,only: deb_Path,minmaxone
      use interfaces,only: inducepcg_pme2gpu,tmatxb_p
     &              , inducestepcg_pme2gpu, efld0_group
     &              , efld0_directgpu2, efld0_directgpu_p
     &              , commfieldfull, commdirdirfull
     &              , inducepcg_group
      use math
      use mpole
      use nvshmem
      use pme
      use polar
      use polar_temp,only: ef,mu
      use polpot
      use potent
      use units
      use uprior
      use utils
      use utilcomm ,buffermpi1=>buffermpi2d,buffermpi2=>buffermpi2d1
      use utilgpu
      use timestat
      use mpi
      use group
      implicit none
c
      integer i, j, k,iploc,iipolegroup
      integer, parameter :: nrhs=2
      integer iglob,iglobgroup
c
c     MPI
c
      integer iipole
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
c
c
      if (.not.use_polar) return
c
 1000 format(' illegal polalg in newinduce.')
  
      if (deb_Path) 
     &   write(*,'(3x,a)') 'newinduce_group'
c
c
c     allocate some memory and clear the arrays:
c
      call prmem_request(mu,3,nrhs,npolegroup,async=.true.)
      call prmem_request(ef,3,nrhs,npolegroup,async=.true.)
      call set_to_zero2(mu,ef,3*nrhs*npolegroup,def_queue)
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))
 
      !call orderbuffergroup
c
c     compute the electric fields:
c
!$acc wait
      if (deb_Path) 
     &   write(*,'(3x,a)') '>>>efld0_group'
      call efld0_group(nrhs,ef)
      if (deb_Path) 
     &   write(*,'(3x,a)') '<<<efld0_group'
c
      call commfieldfull(nrhs,ef)
c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
!$acc parallel loop async default(present) collapse(3)
      do i = 1, npolelocgroup; do k = 1, nrhs; do j = 1, 3
        iglobgroup = ipolegroup(poleglobgroup(i))
        iglob = globglobgroup(iglobgroup)
        iipole = pollist(iglob)
        mu(j,k,i) = polarity(iipole)*ef(j,k,i)
      enddo; enddo; enddo
c
      call commdirdirfull(nrhs,0,mu,reqrec,reqsend)
      call commdirdirfull(nrhs,1,mu,reqrec,reqsend)
      call commdirdirfull(nrhs,2,mu,reqrec,reqsend)
c
      
c
c     now, call the proper solver.
c
      if (polalg.eq.pcg_SId) then
        call inducepcg_group(nrhs,.true.,ef,mu)
      else
         if (rank.eq.0) write(iout,1000) 
         call fatal
      end if

c
c     move the computed dipoles in the module.
c
!$acc parallel loop async default(present) collapse(2)
      do i = 1, npolelocgroup; do j=1,3
        iipolegroup = poleglobgroup(i)
        uindgroup(j,iipolegroup) = mu(j,1,i)
        uinpgroup(j,iipolegroup) = mu(j,2,i)
      end do; end do


      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
      return
      end
c
      subroutine inducepcg_group(nrhs,precnd,ef,mu)
      use atmlst
      use domdec
      use ewald
      use inform     ,only: deb_Path,abort,minmaxone
      use interfaces ,only: commfieldfull, commdirdirfull, tmatxb_group
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent     ,only: use_pmecore
      use polar_temp ,only: res,h,pp,zr,diag
      use timestat
      use tinMemory  ,only: prmem_request
      use sizes
      use units
      use utils
      use utilcomm   ,only: skpPcomm
      use utilgpu
      use mpi
      use group
      implicit none
c
c     solves the polarizable dipoles linear equations by preconditioned
c     conjugate gradient. A diagonal preconditioner is used when precnd
c     is true, otherwise the standard conjugate gradient algorithm is
c     recovered by setting the preconditioner to one.
c
      integer  ,intent(in)   :: nrhs
      logical  ,intent(in)   :: precnd
      real(t_p),intent(in)   :: ef (:,:,:)
      real(t_p),intent(inout):: mu (:,:,:)

      integer i, it, j, k
      real(t_p) ggold(2), alphacg(2)
      real(r_p) gnorm(2), gg(2)
      real(r_p),target :: mbuf(4)
      real(r_p),pointer:: ggnew(:),ene(:)
      real(t_p),save:: ggold1,ggold2
      real(r_p),save:: ggnew1,ggnew2
      real(t_p),save:: ggdev1,ggdev2
      real(r_p),save:: gnorm1,gnorm2
      real(r_p),save:: gg1,gg2
      real(r_p),save:: ene1,ene2
      real(t_p),save:: alphacg1,alphacg2
      real(t_p) reg1,reg2,reg3
      real(t_p) zero, one
      real(r_p) zerom, pt5
      real(r_p) resnrm
      real(8) time1,time2
      logical tinker_isnan_m
      logical,save :: first_in=.true.
      parameter (zero=0.0_ti_p,zerom=0, pt5=0.50_re_p, one=1.0_ti_p)
c
c     MPI
c
      integer iglob, iipole, ierr
      integer iglobgroup
      integer req1, req2, req3, req4
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer,allocatable :: reqendrec(:),reqendsend(:)
      integer,allocatable :: req2endrec(:),req2endsend(:)
c
 1000 format(' cgiter shortreal converged after ',I3,' iterations.',/,
     $       ' final energy        = ',2D14.7,/,
     $       ' final residual norm = ',2D14.7)
 1010 format(' energy and residual norm at iteration ',I3,':',4D12.2)
 1020 format(' Conjugate gradient solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Conjugate gradient solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
 1040 format(' Using a diagonal preconditioner.')
 1050 format(' Induce PCG solver not converged after ',I6,' iterations'
     &    ,/,' Residual norm ', d10.5)

      if (deb_Path) 
     &   write(*,'(3x,a)') 'inducepcg_group'
c
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqendsend(nproc))
      allocate (reqendrec(nproc))
      allocate (req2endsend(nproc))
      allocate (req2endrec(nproc))
      
c
c     allocate some memory and setup the preconditioner:
c
      ggnew(1:2) => mbuf(1:2)
      ene(1:2)   => mbuf(3:4)

      call prmem_request(zr,3,nrhs,max(npolelocgroup,1),
     &     async=.true.)
      call prmem_request(res,3,nrhs,max(npolelocgroup,1),
     &     async=.true.)
      call prmem_request(h ,3,nrhs,max(npolegroup,1),
     &     async=.true.)
      call prmem_request(pp ,3,nrhs,max(npolegroup,1),
     &     async=.true.)
      call prmem_request(diag,max(npolelocgroup,1),async=.true.)

      if (first_in) then
!$acc enter data create(gg1,gg2)
         first_in=.false.
      end if
!$acc data present(res,h,pp,zr,diag,gg1,gg2)
!$acc&     present(mu,ef,poleglobgroup,polarity
!$acc&   ,ipolegroup,pollist,globglobgroup)

      if (precnd) then
!$acc parallel loop async default(present)
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
          iipole = pollist(iglob)
          diag(i) = polarity(iipole)
        end do
        if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
!$acc parallel loop async
        do i = 1, npolelocgroup
          diag(i) = one
        end do
      end if
c
c     initialize
c
      call set_to_zero2(res,zr,3*nrhs*npolelocgroup,def_queue)
      call set_to_zero2( pp, h,3*nrhs*npolegroup ,def_queue)
c
c     now, compute the initial direction
c
      ggold = zero
      ggold1 = zero
      ggold2 = zero
c
!$acc wait
      call tmatxb_group(nrhs,.true.,mu,h)

      call commfieldfull(nrhs,h)
c
!$acc parallel loop collapse(3) async
!$acc&         reduction(+:ggold1,ggold2)
      do k=1,npolelocgroup; do j=1,nrhs; do i=1,3
        res(i,j,k) = ef(i,j,k) - h(i,j,k)
        zr(i,j,k)  = res(i,j,k)*diag(k)
        pp(i,j,k)  = zr(i,j,k)
        if (btest(j,0)) then
          ggold1  = ggold1 + res(i,j,k)*zr(i,j,k)
        else
          ggold2  = ggold2 + res(i,j,k)*zr(i,j,k)
        end if
      end do; end do; end do
!$acc wait

      ggold(1) = ggold1
      ggold(2) = ggold2

      if(nproc>1) then
        call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_TPREC,
     $          MPI_SUM,COMM_TINKER,req1,ierr)
c  
        call commdirdirfull(nrhs,0,pp,reqrec,reqsend)
        call commdirdirfull(nrhs,1,pp,reqrec,reqsend)
        call commdirdirfull(nrhs,2,pp,reqrec,reqsend)
        call MPI_WAIT(req1,status,ierr)
      endif
c
c     now, start the main loop:
c
      do it = 1, politer
c
        call tmatxb_group(nrhs,.true.,pp,h)
        call commfieldfull(nrhs,h)
!$acc serial async
        gg1=zerom; gg2=zerom
!$acc end serial
c
!$acc parallel loop collapse(3) async(def_queue)
        do i = 1,npolelocgroup
           do j = 1,nrhs
              do k = 1, 3
                 if (btest(j,0)) then
                    gg1 = gg1 + pp(k,j,i)*h(k,j,i)
                 else
                    gg2 = gg2 + pp(k,j,i)*h(k,j,i)
                 end if
              end do
           end do
        end do
!$acc update host(gg1,gg2) async(def_queue)
!$acc wait
        gg(1)=gg1; gg(2)=gg2
        if(nproc>1) then
          call MPI_IALLREDUCE(MPI_IN_PLACE,gg(1),nrhs,MPI_RPREC,MPI_SUM,
     $      COMM_TINKER,req2,ierr)
          call MPI_WAIT(req2,status,ierr)
          gg1 = gg(1); gg2 = gg(2)
!$acc update device(gg1,gg2) async(def_queue)
        endif
        do k = 1, nrhs
          if (gg(k).eq.zero) goto 30
        end do
        alphacg  = ggold/gg
        ggnew  = zerom
        ene    = zerom
        ggnew1 = zerom; ggnew2 = zerom;
        ene1   = zerom; ene2   = zerom;
        ggdev1 = ggold1/gg1; ggdev2 = ggold2/gg2

!$acc parallel loop collapse(3) async(def_queue)
        do k=1,npolelocgroup; do j=1,nrhs; do i=1,3
          if (btest(j,0)) then
            mu (i,j,k) = mu (i,j,k) + ggdev1*pp(i,j,k)
            res(i,j,k) = res(i,j,k) - ggdev1*h (i,j,k)
            zr (i,j,k) = diag(k)* res(i,j,k)
            ggnew1     = ggnew1 + res(i,j,k)*zr(i,j,k)
            ene1    = ene1 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
          else
            mu (i,j,k) = mu (i,j,k) + ggdev2*pp(i,j,k)
            res(i,j,k) = res(i,j,k) - ggdev2*h (i,j,k)
            zr (i,j,k) = diag(k)* res(i,j,k)
            ggnew2     = ggnew2 + res(i,j,k)*zr(i,j,k)
            ene2    = ene2 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
          end if
        end do; end do; end do
!$acc wait(def_queue)
        ene2 = -pt5*ene2; ene1 = -pt5*ene1
        ggnew(1)=ggnew1; ggnew(2)=ggnew2
        ene(1)=  ene1;   ene(2)=  ene2        

        if(nproc>1) then
          call MPI_IALLREDUCE(MPI_IN_PLACE,mbuf(1),2*nrhs,MPI_RPREC,
     $          MPI_SUM,COMM_TINKER,req3,ierr)
           call MPI_WAIT(req3,status,ierr)
           ggnew1 = ggnew(1); ggnew2=ggnew(2)
        endif
        resnrm = zero
        do k = 1, nrhs
          gnorm(k) = sqrt(ggnew(k)/real(3*npolegroup,r_p))
          resnrm   = max(resnrm,gnorm(k))
        end do
        if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     $    it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)

        ggdev1 = ggnew(1)/ggold(1); ggdev2 = ggnew(2)/ggold(2)
        ggold(1:2) = ggnew(1:2)
        ggold1 = ggold(1); ggold2 = ggold(2)

        if ((it.eq.politer.and.resnrm.gt.poleps)
     &     .or.tinker_isnan_m(gnorm(1))) then
           ! Not converged Abort
           write(0,1050) it,max(gnorm(1),gnorm(2))
           abort = .true.
        endif
        if (resnrm.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     $      (ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
          goto 10
        end if

!$acc parallel loop collapse(3) async(def_queue)
        do k = 1,npolelocgroup; do j = 1,nrhs; do i = 1, 3
          if (btest(j,0)) then
            pp(i,j,k) = zr(i,j,k) + ggdev1*pp(i,j,k)
          else
            pp(i,j,k) = zr(i,j,k) + ggdev2*pp(i,j,k)
          end if
        end do; end do; end do
c
        call commdirdirfull(nrhs,0,pp,reqrec,reqsend)
        call commdirdirfull(nrhs,1,pp,reqrec,reqsend)
        call commdirdirfull(nrhs,2,pp,reqrec,reqsend)

      end do
 10   continue
c
      call commdirdirfull(nrhs,0,mu,reqendrec,reqendsend)
      call commdirdirfull(nrhs,1,mu,reqendrec,reqendsend)
      call commdirdirfull(nrhs,2,mu,reqendrec,reqendsend)
c
c     finalize and debug printing:
c
      if (polprt.ge.3) then
        write(iout,1020)
!$acc wait
!$acc update host(poleglobgroup,mu)
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
          write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
!$acc wait
!$acc update host(poleglobgroup,mu)
        do i = 1, npolelocgroup
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
          write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
        end do
      end if

  30  continue

!$acc end data
      return
      end


      subroutine efld0_group(nrhs,ef)
c
c     Compute the direct space contribution to the permanent electric field.
c      Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      !use couple  ,only: i12,i13,i14,i15,n12,n13,n14,n15
      use domdec  ,only: rank,loc,nbloc
      use ewald   ,only: aewald
      use newinduce_group_inl
      use inform  ,only: deb_Path
      use interfaces ,only: efld0_group_correct_scaling
      use math    ,only: sqrtpi
      use mpole   ,only: npolebloc,ipole,rpole,npolelocnl,pollist
      use polar   ,only: pdamp,thole
      use potent  , only : use_polarshortreal
      !use polgrp  ,only: ip11,ip12,ip13,ip14,np11,np12,np13,np14
      use shunt   ,only: cut2
      use utilgpu ,only: dir_queue,rec_queue,def_queue
#ifdef _OPENACC
     &                  ,dir_stream,stream_wait_async,
     &                   rec_stream,rec_event
#endif
      use timestat   ,only: timer_enter,timer_exit,timer_efld0_direct
      use tinheader  ,only: ti_p,re_p
      use group
      implicit none
      integer  , intent(in)    :: nrhs
      real(t_p), intent(inout) :: ef(:,:,:)

      integer :: i,iglob,iploc,kk,kkk
      integer :: ii,j,k,ksp,kd,iipole,kbis,kpole,kglob
      integer :: iipolegroup, iglobgroup,ilocgroup
      integer :: kglobgroup, klocgroup
      real(t_p) pdi,pti
      real(t_p) xi,yi,zi,d2
      real(t_p) pgamma,damp
      real(t_p) alsq2, alsq2n
      real(t_p), parameter ::  one=1.0_ti_p
      real(t_p) pscale,dscale
      real(t_p) d,bn1,bn2,sc3,sc5
      type(rpole_elt) ip,kp
      type(real3) fid,fip,fkd,fkp,pos
c
c
      if (deb_Path) then
!$acc wait
       write(0,'(3x,a)') 'efld0_group'
      endif
     
!$acc parallel loop gang vector_length(32)
!$acc&         present(ef)
!$acc&         present(poleglobgroup,ipolegroup,globglobgroup,
!$acc&  pollist ,locgroup,pdamp,polelocgroup,
!$acc&  thole,x,y,z,rpole)
!$acc&         private(ip)
!$acc&         async(def_queue)
      do ii = 1, npolelocgroup
        iipolegroup = poleglobgroup(ii)
        iglobgroup = ipolegroup(iipolegroup)
        iglob = globglobgroup(iglobgroup)
        iipole = pollist(iglob)

         pdi      = pdamp(iipole)
         pti      = thole(iipole)
         xi       = x    (iglob) 
         yi       = y    (iglob) 
         zi       = z    (iglob) 

         ip%c     = rpole(01,iipole)
         ip%dx    = rpole(02,iipole)
         ip%dy    = rpole(03,iipole)
         ip%dz    = rpole(04,iipole)
         ip%qxx   = rpole(05,iipole)
         ip%qxy   = rpole(06,iipole)
         ip%qxz   = rpole(07,iipole)
         ip%qyy   = rpole(09,iipole)
         ip%qyz   = rpole(10,iipole)
         ip%qzz   = rpole(13,iipole)
c
!$acc loop vector private(kp,fip,fid,fkp,fkd,pos)
        do kkk = iipolegroup+1, npolegroup
          kglobgroup = ipolegroup(kkk)
          kglob = globglobgroup(kglobgroup)
          kpole = pollist(kglob)
          kk = polelocgroup(kkk)
          klocgroup = locgroup(kglobgroup)

          pos%x  = x (kglob) - xi
          pos%y  = y (kglob) - yi
          pos%z  = z (kglob) - zi

          d2     = pos%x**2 + pos%y**2 + pos%z**2

          kp%c   = rpole( 1, kpole)
          kp%dx  = rpole( 2, kpole)
          kp%dy  = rpole( 3, kpole)
          kp%dz  = rpole( 4, kpole)
          kp%qxx = rpole( 5, kpole)
          kp%qxy = rpole( 6, kpole)
          kp%qxz = rpole( 7, kpole)
          kp%qyy = rpole( 9, kpole)
          kp%qyz = rpole(10, kpole)
          kp%qzz = rpole(13, kpole)

          damp   = pdi * pdamp(kpole)
          pgamma = min( pti,thole(kpole) )

          call efld0_couple(d2,pos,ip,kp,0._ti_p,0._ti_p,
     &            0._ti_p,damp,pgamma,1.0_ti_p,1.0_ti_p,
     &            fid,fip,fkd,fkp,d,bn1,bn2,sc3,sc5,.false.)

!$acc atomic update
            ef(1,1,ii) = ef(1,1,ii) + fid%x
!$acc atomic update
            ef(2,1,ii) = ef(2,1,ii) + fid%y
!$acc atomic update
            ef(3,1,ii) = ef(3,1,ii) + fid%z
!$acc atomic update
            ef(1,2,ii) = ef(1,2,ii) + fip%x
!$acc atomic update
            ef(2,2,ii) = ef(2,2,ii) + fip%y
!$acc atomic update
            ef(3,2,ii) = ef(3,2,ii) + fip%z
!$acc atomic update
            ef(1,1,kk)  = ef(1,1,kk ) + fkd%x
!$acc atomic update       
            ef(2,1,kk)  = ef(2,1,kk ) + fkd%y
!$acc atomic update       
            ef(3,1,kk)  = ef(3,1,kk ) + fkd%z
!$acc atomic update       
            ef(1,2,kk)  = ef(1,2,kk ) + fkp%x
!$acc atomic update       
            ef(2,2,kk)  = ef(2,2,kk ) + fkp%y
!$acc atomic update       
            ef(3,2,kk)  = ef(3,2,kk ) + fkp%z

        end do
      end do


c
      call efld0_group_correct_scaling(nrhs,ef)

      return
      end

       subroutine efld0_group_correct_scaling(nrhs,ef)
      use atmlst 
      use group
      use atoms    ,only: x,y,z
      use domdec   ,only: rank,loc,nbloc
      use ewald    ,only: aewald
      use newinduce_group_inl
      use inform   ,only: deb_Path
      use math     ,only: sqrtpi
      use mpole    ,only: npolebloc,ipole,rpole,poleloc,npolelocnl
      use neigh    ,only: nelst,elst
      use potent  , only : use_polarshortreal
      use polar    ,only: pdamp,thole
      use shunt    ,only: cut2
      use utilgpu  ,only: dir_queue,def_queue
      use timestat ,only: timer_enter,timer_exit,timer_efld0_direct
      use tinheader,only: ti_p
      implicit none

      ! shape(ef) = (/3,nrhs,npolegroup/)
      integer, intent(in) :: nrhs
      real(t_p),intent(inout):: ef(:,:,:)

      integer i,iglob,iploc,kk
      integer ii,j,k,ksp,kd,iipole,kbis,kpole,kglob
      integer iglobgroup, kglobgroup
      real(t_p) pdi,pti
      real(t_p) xi,yi,zi,d2
      real(t_p) thole1,pgamma,damp
      real(t_p) one
      real(t_p) pscale,dscale
      real(t_p) d,bn1,bn2,sc3,sc5
      type(real3) fid,fip,fkd,fkp,pos
      type(rpole_elt) ip,kp
      real(t_p) fgrp
      integer iga,igb

      parameter(one =1.0_ti_p)
c
      if (deb_Path) write(*,'(4x,a)') 'efld0_group_correct_scaling'


!$acc parallel loop gang vector
!$acc&         present(ef)
!$acc&         present(loclocgroup,locgroup,polelocgroup,pollistgroup,
!$acc&  pdamp,thole,x,y,z,rpole,nelst,elst,poleloc,
!$acc&  dpcorrect_ik_group,dpcorrect_scale_group)
!$acc&         private(fid,fip,fkd,fkp,pos,ip,kp)
!$acc&         async(def_queue)
      do ii = 1, n_dpscale_group
         iipole = dpcorrect_ik_group(2*(ii-1)+1)
         kpole  = dpcorrect_ik_group(2*(ii-1)+2)

         dscale = dpcorrect_scale_group(2*(ii-1)+1)
         pscale = dpcorrect_scale_group(2*(ii-1)+2)

         iglob  = ipole     (iipole)
         kglob  = ipole  (kpole)

         iglobgroup = loclocgroup(iglob)
         i = locgroup(iglobgroup)
         iploc = polelocgroup(pollistgroup(iglobgroup))

         kglobgroup = loclocgroup(kglob)
         k = locgroup(kglobgroup)
         kbis = polelocgroup(pollistgroup(kglobgroup))

         if (iploc.lt.1.or.iploc.gt.npolegroup.or.
     &        kbis.lt.1.or. kbis.gt.npolegroup) cycle

         pdi    = pdamp(iipole)
         pti    = thole(iipole)

         pos%x  = x(kglob) - x(iglob)
         pos%y  = y(kglob) - y(iglob)
         pos%z  = z(kglob) - z(iglob)
         !call image_inl(pos%x,pos%y,pos%z)
         d2     = pos%x**2 + pos%y**2 + pos%z**2

         ip%c   = rpole(01,iipole)
         ip%dx  = rpole(02,iipole)
         ip%dy  = rpole(03,iipole)
         ip%dz  = rpole(04,iipole)
         ip%qxx = rpole(05,iipole)
         ip%qxy = rpole(06,iipole)
         ip%qxz = rpole(07,iipole)
         ip%qyy = rpole(09,iipole)
         ip%qyz = rpole(10,iipole)
         ip%qzz = rpole(13,iipole)

         kp%c   = rpole(01, kpole)
         kp%dx  = rpole(02, kpole)
         kp%dy  = rpole(03, kpole)
         kp%dz  = rpole(04, kpole)
         kp%qxx = rpole(05, kpole)
         kp%qxy = rpole(06, kpole)
         kp%qxz = rpole(07, kpole)
         kp%qyy = rpole(09, kpole)
         kp%qyz = rpole(10, kpole)
         kp%qzz = rpole(13, kpole)

         thole1 = thole(kpole)
         damp   = pdi * pdamp(kpole)
         pgamma = min( pti,thole1 )

         call efld0_couple(d2,pos,ip,kp,0._ti_p,0._ti_p,
     &              0.0_ti_p,damp,pgamma,dscale,pscale,
     &              fid,fip,fkd,fkp,d,bn1,bn2,sc3,sc5,.true.)

         if (dscale.ne.0.0_ti_p) then
!$acc atomic update
            ef(1,1,iploc) = ef(1,1,iploc) + fid%x
!$acc atomic update
            ef(2,1,iploc) = ef(2,1,iploc) + fid%y
!$acc atomic update
            ef(3,1,iploc) = ef(3,1,iploc) + fid%z
!$acc atomic update
            ef(1,1,kbis)  = ef(1,1,kbis ) + fkd%x
!$acc atomic update       
            ef(2,1,kbis)  = ef(2,1,kbis ) + fkd%y
!$acc atomic update       
            ef(3,1,kbis)  = ef(3,1,kbis ) + fkd%z
         end if
         if (pscale.ne.0.0_ti_p) then
!$acc atomic update
            ef(1,2,iploc) = ef(1,2,iploc) + fip%x
!$acc atomic update
            ef(2,2,iploc) = ef(2,2,iploc) + fip%y
!$acc atomic update
            ef(3,2,iploc) = ef(3,2,iploc) + fip%z
!$acc atomic update       
            ef(1,2,kbis)  = ef(1,2,kbis ) + fkp%x
!$acc atomic update       
            ef(2,2,kbis)  = ef(2,2,kbis ) + fkp%y
!$acc atomic update       
            ef(3,2,kbis)  = ef(3,2,kbis ) + fkp%z
         end if
         end do
      end subroutine
c
      subroutine tmatxb_group(nrhs,dodiag,mu,efi)
c
c     Compute the direct space contribution to the electric field due to the current value
c     of the induced dipoles
c
      use atoms   , only : x, y, z
      use atmlst  , only : poleglobnl
      !use erf_mod
      use ewald   , only : aewald
      use inform  , only : deb_Path
      use interfaces ,only: tmatxb_correct_interactions_group
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl,pollist
      use neigh   , only : nelst,nelstc,elst,shortelst,nshortelstc
      use polar   , only : polarity, thole, pdamp,tinypol
      use shunt   , only : cut2
      use utils   , only : set_to_zero1
      use utilgpu , only : maxscaling1, def_queue, dir_queue,real3,real6
      use domdec  , only : rank, loc
      use polgrp  , only : np11, np12, np13, np14, ip11, ip12, ip13,ip14
      use polpot  , only : u1scale, u2scale, u3scale, u4scale
      use potent  , only : use_polarshortreal
      use tinheader ,only: ti_p
      use newinduce_group_inl
      use group
      implicit none
      integer, intent(in) ::  nrhs
      logical, intent(in) ::  dodiag
      real(t_p), intent(in) ::  mu(:,:,:)
      real(t_p), intent(inout) ::  efi(:,:,:)
      integer iglobgroup,kglobgroup
      integer :: iipolegroup
      integer  kk, kkk
      integer    :: i, j, k, irhs
      integer    :: start,finla,sized        ! Indexes associated to Matvec split
      integer    :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer    :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      integer,pointer :: lst(:,:),nlst(:)    ! pointers to neighbor list
      real(t_p)  :: d2                     ! Distance (d2 squared) between atoms
      real(t_p)  :: alsq2, alsq2n          ! ewald precomputed constants
      real(t_p)  :: sdamp,pgamma
      real(t_p)  :: tmp,ipolar
      type(real6):: dpui,dpuk  ! store coefficients of mu
      type(real3):: posi,dist  ! position of i pole and distance between i and k
      type(real3):: fid,fip    ! Contribution of local pair to efi(iploc)
      type(real3):: fkd,fkp    ! Contribution of local pair to efi(iploc)
      real(t_p) uscale,fgrp ! Scaling factor for interaction
      integer :: iga,igb
      
      if (deb_Path) print '(2X,A)','tmatxbgroup'

c
c     initialize the result vector
c
      call set_to_zero1(efi,3*nrhs*npolegroup,def_queue)
!$acc wait
      if (deb_Path) print '(2X,A)','tmatxbgroup start'
c
c     gather some parameters, then set up the damping factors.
c
c
!$acc parallel loop gang vector_length(32)
!$acc&         present(mu,efi,polarity,pollist,ipolegroup,polelocgroup,
!$acc&  loc,x,y,z,pdamp,thole,poleglobgroup,globglobgroup)
!$acc&         private(posi,dpui)
!$acc&         async(def_queue)
      do ii = 1, npolelocgroup
        iipolegroup =  poleglobgroup(ii)
        iglobgroup = ipolegroup(iipolegroup)
        iglob = globglobgroup(iglobgroup)
        iipole = pollist(iglob)

        posi%x  = x(iglob)
        posi%y  = y(iglob)
        posi%z  = z(iglob)
        dpui%x  = mu(1,1,ii)
        dpui%y  = mu(2,1,ii)
        dpui%z  = mu(3,1,ii)
        dpui%xx = mu(1,2,ii)
        dpui%yy = mu(2,2,ii)
        dpui%zz = mu(3,2,ii)
c
!$acc loop vector private(dist,dpuk,fid,fip,fkp,fkd)
        do kkk = iipolegroup+1, npolegroup
          kglobgroup = ipolegroup(kkk)
          kglob = globglobgroup(kglobgroup)
          kpole = pollist(kglob)
          kk = polelocgroup(kkk)
c
c     compute the distances and the scaling factors according to
c     Thole's model.
c
          dist%x = x(kglob) - posi%x
          dist%y = y(kglob) - posi%y
          dist%z = z(kglob) - posi%z
          !call image_inl(dist%x,dist%y,dist%z)
          d2 = dist%x*dist%x + dist%y*dist%y + dist%z*dist%z

          sdamp   = pdamp(iipole) * pdamp(kpole)
          pgamma  = min(thole(iipole),thole(kpole))

          dpuk%x   = mu(1,1,kk)
          dpuk%y   = mu(2,1,kk)
          dpuk%z   = mu(3,1,kk)
          dpuk%xx  = mu(1,2,kk)
          dpuk%yy  = mu(2,2,kk)
          dpuk%zz  = mu(3,2,kk)

          call tmatxb_couple(d2,dist,dpui,dpuk,
     &           sdamp,pgamma,0._ti_p,0._ti_p,0._ti_p,1.0_ti_p,
     &                         fid,fip,fkd,fkp,.false.)


            ! increment electric field for each atoms
!$acc atomic update
            efi(1,1,ii)    = efi(1,1,ii)    + fid%x
!$acc atomic update
            efi(2,1,ii)    = efi(2,1,ii)    + fid%y
!$acc atomic update
            efi(3,1,ii)    = efi(3,1,ii)    + fid%z
!$acc atomic update
            efi(1,2,ii)    = efi(1,2,ii)    + fip%x
!$acc atomic update
            efi(2,2,ii)    = efi(2,2,ii)    + fip%y
!$acc atomic update
            efi(3,2,ii)    = efi(3,2,ii)    + fip%z
!$acc atomic update
            efi(1,1,kk) = efi(1,1,kk) + fkd%x
!$acc atomic update
            efi(2,1,kk) = efi(2,1,kk) + fkd%y
!$acc atomic update
            efi(3,1,kk) = efi(3,1,kk) + fkd%z
!$acc atomic update
            efi(1,2,kk) = efi(1,2,kk) + fkp%x
!$acc atomic update
            efi(2,2,kk) = efi(2,2,kk) + fkp%y
!$acc atomic update
            efi(3,2,kk) = efi(3,2,kk) + fkp%z
        end do
      end do

      call tmatxb_correct_interactions_group(nrhs,mu,efi)

      if(dodiag) then
c
c     if dodiag is true, also compute the "self-induced" field,
c     i.e., the diagonal portion of the matrix/vector product.
c
!$acc  parallel loop collapse(3) async(def_queue)
!$acc&          default(present)
        do i = 1, npolelocgroup; do irhs = 1, nrhs; do j = 1, 3
          iglobgroup = ipolegroup(poleglobgroup(i))
          iglob = globglobgroup(iglobgroup)
          iipole = pollist(iglob)

          ipolar = polarity(iipole)
          !if no polarisability, take a negligeable value to allow convergence
          if (ipolar  == 0.0_ti_p) then
              tmp = mu(j,irhs,i)*(1.0_ti_p/tinypol)
          else
              tmp = mu(j,irhs,i)/ipolar
          end if

          efi(j,irhs,i) = efi(j,irhs,i) + tmp

        end do; enddo; enddo
      end if

      end

      subroutine tmatxb_correct_interactions_group(nrhs,mu,efi)
      use atoms   , only : x, y, z
      use atmlst  , only : poleglobnl
      !use erf_mod
      use ewald   , only : aewald
      use math    , only : sqrtpi
      use mpole   , only : ipole,poleloc,npolelocnl,npolebloc
      use polar   , only : polarity, thole, pdamp
      use shunt   , only : cut2
      use utilgpu , only : def_queue, dir_queue,real3,real6
      use utilcomm, only : no_commdir
      use domdec  , only : loc
      use polpot  , only : n_uscale,ucorrect_ik,ucorrect_scale
      use tinheader ,only: ti_p
      use newinduce_group_inl
      use group
      use inform
      implicit none
      integer, intent(in) :: nrhs
      real(t_p) ,intent(in) :: mu (:,:,:)
      real(t_p) ,intent(out):: efi(:,:,:)

      integer :: iglobgroup,kglobgroup,kbis
      integer    :: i, j, k, irhs
      integer    :: ii, iipole, iglob, iploc ! Indexes associated to local atom
      integer    :: kpole, kpoleloc, kglob   ! Indexes associated to neighbor atom
      real(t_p)  :: d2                     ! Distance (d2 squared) between atoms
      real(t_p)  :: alsq2, alsq2n          ! ewald precomputed constants
      real(t_p)  :: sdamp,pgamma
      type(real6):: dpui,dpuk  ! store coefficients of mu
      type(real3):: posi,dist  ! position of i pole and distance between i and k
      type(real3):: fid,fip    ! Contribution of local pair to efi(iploc)
      type(real3):: fkd,fkp    ! Contribution of local pair to efi(iploc)
      integer :: iga,igb
      real(t_p) uscale,fgrp ! Scaling factor for interaction

      if (n_uscale_group.eq.0) return

      if (deb_Path) print '(2X,A)','tmatxb_correct_interactions_group'

      ! Scaling corrections loop
!$acc parallel loop default(present)
!$acc&         private(dpui,dpuk,dist,fid,fip,fkd,fkp)
!$acc&         async(def_queue)
      do ii = 1,n_uscale_group
         iipole   = ucorrect_ik_group(2*(ii-1)+1)
         kpole    = ucorrect_ik_group(2*(ii-1)+2)
         uscale   = ucorrect_scale_group(ii)

         iglob  = ipole     (iipole)
         kglob  = ipole  (kpole)

         iglobgroup = loclocgroup(iglob)
         i = locgroup(iglobgroup)
         iploc = polelocgroup(pollistgroup(iglobgroup))

         kglobgroup = loclocgroup(kglob)
         k = locgroup(kglobgroup)
         kbis = polelocgroup(pollistgroup(kglobgroup))
         !skip atom if it is not polarizable
         !FIXME do we need to test on polarity
         if (polarity(iipole) == 0) cycle
         if (iploc.eq.0.or.iploc.gt.npolegroup) cycle
         if (kbis.eq.0.or.kbis.gt.npolegroup) cycle

         dpui%x   = mu(1,1,iploc)
         dpui%y   = mu(2,1,iploc)
         dpui%z   = mu(3,1,iploc)
         dpui%xx  = mu(1,2,iploc)
         dpui%yy  = mu(2,2,iploc)
         dpui%zz  = mu(3,2,iploc)

         dpuk%x   = mu(1,1,kbis)
         dpuk%y   = mu(2,1,kbis)
         dpuk%z   = mu(3,1,kbis)
         dpuk%xx  = mu(1,2,kbis)
         dpuk%yy  = mu(2,2,kbis)
         dpuk%zz  = mu(3,2,kbis)

         dist%x   = x(kglob) - x(iglob)
         dist%y   = y(kglob) - y(iglob)
         dist%z   = z(kglob) - z(iglob)
         !call image_inl(dist%x,dist%y,dist%z)
         d2  = dist%x*dist%x + dist%y*dist%y + dist%z*dist%z

         sdamp    = pdamp(iipole) * pdamp(kpole)
         pgamma   = min(thole(iipole),thole(kpole))

         call tmatxb_couple(d2,dist,dpui,dpuk,
     &        sdamp,pgamma,0._ti_p,0._ti_p,0._ti_p,uscale,
     &                      fid,fip,fkd,fkp,.true.)

         ! increment electric field for each atoms
!$acc atomic update
         efi(1,1,iploc)    = efi(1,1,iploc)    + fid%x
!$acc atomic update
         efi(2,1,iploc)    = efi(2,1,iploc)    + fid%y
!$acc atomic update
         efi(3,1,iploc)    = efi(3,1,iploc)    + fid%z
!$acc atomic update
         efi(1,2,iploc)    = efi(1,2,iploc)    + fip%x
!$acc atomic update
         efi(2,2,iploc)    = efi(2,2,iploc)    + fip%y
!$acc atomic update
         efi(3,2,iploc)    = efi(3,2,iploc)    + fip%z
!$acc atomic update
         efi(1,1,kbis) = efi(1,1,kbis) + fkd%x
!$acc atomic update
         efi(2,1,kbis) = efi(2,1,kbis) + fkd%y
!$acc atomic update
         efi(3,1,kbis) = efi(3,1,kbis) + fkd%z
!$acc atomic update
         efi(1,2,kbis) = efi(1,2,kbis) + fkp%x
!$acc atomic update
         efi(2,2,kbis) = efi(2,2,kbis) + fkp%y
!$acc atomic update
         efi(3,2,kbis) = efi(3,2,kbis) + fkp%z
      end do
      end

c
c    subroutine orderbuffergroup : get sizes and indexes of all the domains
c
      
      subroutine orderbuffergroup
      use atoms
      use domdec
      use group
      use mpole
      use potent
      use mpi
      use inform
      implicit none
      integer i,iproc,tag,status(MPI_STATUS_SIZE)
      integer ierr,commloc,idomlen,ibeg
      integer count
      integer, allocatable :: reqrec(:), reqsend(:)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return

      if (deb_Path) print '(2X,A)','>>>orderbuffergroup'

      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
c     deal with atoms first
c
      domlengroup(rank+1) = nlocatgroup

      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(domlengroup(iproc+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlengroup(rank+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
      bufbeggroup(rank+1) = 1
      count = domlengroup(rank+1)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          if (domlengroup(iproc+1).eq.0) then
            bufbeggroup(iproc+1) = 1
          else
            bufbeggroup(iproc+1) = count + 1
          end if
          count = count + domlengroup(iproc+1)
        end if
      end do
c
c     get the indexes of the atoms in the group frame
c
!$acc host_data use_device(globgroup)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_IRECV(globgroup(bufbeggroup(iproc+1)),
     $     domlengroup(iproc+1),MPI_INT,iproc,tag,
     $     COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*iproc + rank + 1
          call MPI_ISEND(globgroup,domlengroup(rank+1),MPI_INT,
     $     iproc,tag,COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
!$acc end host_data
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
        end if
      end do

      do iproc = 1, nproc
        if (iproc.ne.(rank+1)) then
          idomlen = domlengroup(iproc)
          ibeg    = bufbeggroup(iproc)
!$acc parallel loop async present(globgroup,locgroup)
          do i = 1, idomlen
            locgroup(globgroup(ibeg+i-1)) =
     $        ibeg+i-1
          end do
        end if
      end do
c
c     then deal with multipoles
c
      domlenpolegroup(rank+1) = npolelocgroup
c
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*rank + iproc + 1
        call MPI_IRECV(domlenpolegroup(iproc+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
        tag = nproc*iproc + rank + 1
        call MPI_ISEND(domlenpolegroup(rank+1),1,MPI_INT,iproc,tag,
     $   COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
c
      bufbegpolegroup(rank+1) = 1
      count = domlenpolegroup(rank+1)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          if (domlenpolegroup(iproc+1).eq.0) then
            bufbegpolegroup(iproc+1) = 1
          else
            bufbegpolegroup(iproc+1) = count + 1
          end if
          count = count + domlenpolegroup(iproc+1)
        end if
      end do
c
c     get the indexes of the mulipoles in the group frame
c
!$acc host_data use_device(poleglobgroup)
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*rank + iproc + 1
          call MPI_IRECV(poleglobgroup(bufbegpolegroup(iproc+1)),
     $     domlenpolegroup(iproc+1),MPI_INT,iproc,tag,
     $     COMM_TINKER,reqrec(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          tag = nproc*iproc + rank + 1
          call MPI_ISEND(poleglobgroup,domlenpolegroup(rank+1),MPI_INT,
     $     iproc,tag,COMM_TINKER,reqsend(iproc+1),ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqrec(iproc+1),status,ierr)
        end if
      end do
      do iproc = 0, nproc-1
        if (iproc.ne.rank) then
          call MPI_WAIT(reqsend(iproc+1),status,ierr)
        end if
      end do
!$acc end host_data

      do iproc = 1, nproc
        if (iproc.ne.(rank+1)) then
          idomlen = domlenpolegroup(iproc)
          ibeg    = bufbegpolegroup(iproc)
!$acc parallel loop async present(poleglobgroup,polelocgroup)
          do i = 1, idomlen
            polelocgroup(poleglobgroup(ibeg+i-1)) =
     $        ibeg+i-1
          end do
        end if
      end do
c
      deallocate (reqrec)
      deallocate (reqsend)
!$acc wait
      if (deb_Path) print '(2X,A)','<<<orderbuffergroup'
      return
      end
c
c     subroutine commfieldfull : communicate some direct fields (Newton s third law)
c
      subroutine commfieldfull(nrhs,ef)
      use atoms
      use domdec
      use group
      use mpole
      use potent
      use mpi
      use utilgpu ,only: prmem_request,dir_queue
      use utilcomm,only: buff_field
      implicit none
      integer, intent(in) :: nrhs
      real(t_p), intent(inout) ::  ef(:,:,:)
      integer i,j,k,l,tag,status(MPI_STATUS_SIZE),c
      integer ierr,commloc
      integer, allocatable :: reqrec(:), reqsend(:)
c
      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return

      if (use_pmecore) then
        commloc = comm_dir
      else
        commloc = COMM_TINKER
      end if
c
      call prmem_request(buff_field,3,nrhs,
     &        max(npolelocgroup,1),nproc,async=.true.)
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
c
!$acc wait
!$acc host_data use_device(buff_field)
      c=0
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          c=c+1
          tag = nproc*rank + i
          call MPI_IRECV(buff_field(1,1,1,c)
     $     ,3*nrhs*npolelocgroup,MPI_TPREC,
     $     i-1,0,commloc,reqrec(i),ierr)
        end if
      end do
!$acc end host_data
c
!$acc host_data use_device(ef)
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*(i-1) + rank + 1
          call MPI_ISEND(ef(1,1,bufbegpolegroup(i)),
     $     3*nrhs*domlenpolegroup(i),MPI_TPREC,i-1,
     $    0,commloc,reqsend(i),ierr)
        end if
      end do
!$acc end host_data
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*rank + i
          call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*(i-1) + rank + 1
          call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
!$acc parallel loop collapse(3) present(ef,buff_field) async
        do j = 1, npolelocgroup; do k = 1, nrhs; do l=1,3
!$acc loop seq
          do c = 1, nproc-1
            ef(l,k,j) = ef(l,k,j) + buff_field(l,k,j,c)
          enddo
        end do; enddo; enddo
c
      deallocate (reqrec)
      deallocate (reqsend)
      return
      end
c
c    subroutine commdirdirfull : deal with communications of direct dipoles
c
c    rule determines what to do:
c        - 0: start reception of direct dipoles
c        - 1: start sending of direct dipoles 
c        - 2: wait for the communications to be done
c
      subroutine commdirdirfull(nrhs,rule,mu,reqrec,reqsend)
      use domdec
      use group
      use iounit
      use mpole
      use potent
      use mpi
      use inform  , only : deb_Path
      implicit none
      integer, intent(in) :: nrhs,rule
      integer, intent(inout) :: reqrec(nproc),reqsend(nproc)
      real(t_p), intent(inout) ::  mu(:,:,:)
      integer ierr,status(MPI_STATUS_SIZE),tag,i
 1000 format(' illegal rule in commdirdirfull.')
 41   format(7x,'>> ',A20,   3x,'recv')
 42   format(7x,   3x,A20,' >>','send')
 43   format(7x,'<< ',A20,   3x,'wait')

      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return

      if (rule.eq.0) then
        if (deb_Path) write(*,41) 'commdirdirfull       '
c
c     MPI : begin reception
c
!$acc host_data use_device(mu)
        do i = 1, nproc
          if (i.ne.(rank+1)) then
            tag = nproc*rank + i
            call MPI_IRECV(mu(1,1,bufbegpolegroup(i)),
     $       3*nrhs*domlenpolegroup(i),MPI_TPREC,i-1,tag,
     $       COMM_TINKER,reqrec(i),ierr)
          end if
        end do
!$acc end host_data

      else if (rule.eq.1) then
         if (deb_Path) write(*,42) 'commdirdirfull       '
c
!$acc host_data use_device(mu)
        do i = 1, nproc
          if (i.ne.(rank+1)) then
            tag = nproc*(i-1) + rank + 1
            call MPI_ISEND(mu,3*nrhs*npolelocgroup,
     $         MPI_TPREC,i-1,tag,COMM_TINKER,
     $         reqsend(i),ierr)
          end if
        end do
!$acc end host_data
        if (tinkerdebug.gt.0) call MPI_BARRIER(hostcomm,ierr)

      else if (rule.eq.2) then

        do i = 1, nproc
          if (i.ne.(rank+1)) then
            call MPI_WAIT(reqrec(i),status,ierr)
          end if
        end do
        do i = 1, nproc
          if (i.ne.(rank+1)) then
            call MPI_WAIT(reqsend(i),status,ierr)
          end if
        end do

      else
        if (rank.eq.0) write(iout,1000) 
        call fatal 
      end if

      return
      end
c
c     subroutine commforcesgroup : deal with communications of "polarization group" forces
c
      subroutine commforcesgroup
      use sizes
      use atoms
      use deriv
      use domdec
      use group
      use potent
      use mpi
      use tinMemory ,only: prmem_requestm
      use utilcomm  ,only: buffMpi_f1
      use inform  , only : deb_Path
      implicit none
      integer i,j,k,tag,ierr,sz,c
      integer  , allocatable:: reqsend(:),reqrec(:)
      mdyn_rtyp, pointer:: buffer(:,:,:)
      integer status(MPI_STATUS_SIZE)

      if (nproc.eq.1.or.((use_pmecore).and.(rank.gt.ndir-1))) return
      if (deb_Path) write(*,*) '   >> commforcesgroup'
c
c     allocate some arrays
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))

      sz = 3*max(1,nlocatgroup)*nproc
      call prmem_requestm(buffMpi_f1,sz,async=.false.)
      buffer(1:3,1:max(1,nlocatgroup),1:nproc) => buffMpi_f1(1:sz)
!$acc enter data attach(buffer)

c
c     MPI : begin reception in buffer
c
!$acc host_data use_device(buffer)
      c=0
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          c=c+1
          tag = nproc*rank + i
          call MPI_IRECV(buffer(1,1,c),3*nlocatgroup,
     $     MPI_MDTYP,i-1,tag,COMM_TINKER,reqrec(i),ierr)
        end if
      end do
!$acc end host_data
c
c     MPI : begin sending
c
!$acc host_data use_device(depgroup)
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          tag = nproc*(i-1) + rank + 1
          call MPI_ISEND(depgroup(1,bufbeggroup(i)),
     $     3*domlengroup(i),MPI_MDTYP,i-1,tag,COMM_TINKER,
     $     reqsend(i),ierr)
        end if
      end do
!$acc end host_data
c
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          call MPI_WAIT(reqrec(i),status,ierr)
        end if
      end do
      do i = 1, nproc
        if (i.ne.(rank+1)) then
          call MPI_WAIT(reqsend(i),status,ierr)
        end if
      end do
c
c     MPI : move in global arrays
c
!$acc parallel loop collapse(2) async present(depgroup,buffer)
      do j=1,nlocatgroup; do k=1,3
!$acc loop seq
        do c = 1, nproc-1
          depgroup(k,j) = depgroup(k,j) +  buffer(k,j,c)
        end do
      enddo; enddo
      deallocate (reqsend)
      deallocate (reqrec)

!$acc exit data detach(buffer) async
      nullify(buffer)
      if (deb_Path) write(*,*) '   << commforcesgroup'
      
      return
      end
