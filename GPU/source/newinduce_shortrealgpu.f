 
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c      
c     evaluate induced dipole moments and the polarization energy
c     using either a (preconditioned) conjugate gradient algorithm or
c     Jacobi iterations coupled with DIIS extrapolation.
c      
c     literature reference:
c     "Scalable Evaluation of Polarization Energy and Associated Forces
c     in Polarizable Molecular Dynamics: II. Toward Massively Parallel
c     Computations Using Smooth Particle Mesh Ewald",L. Lagardere et al.,
c     J. Chem. Theory Comput., 2015, 11 (6), pp 2589–2599
c
#include "tinker_precision.h"
      subroutine newinduce_shortrealgpu
      use atoms ,only: n,n_pe
      use atmlst
      use domdec
      use ewald
      use iounit
      use interfaces,only:inducepcg_shortrealgpu
     &              ,inducestepcg_pme2gpu
     &              ,tmatxb_p1,efld0_directgpu2
     &              ,efld0_directgpu_p1
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use polar_temp,only:mu,ef,murec
      use timestat
      use tinMemory
      use units
      use uprior
      use utils
      use utilgpu
      use mpi
      implicit none
c
      integer i, j, k, nrhs
c
c     MPI
c
      integer iipole
      integer, dimension(nproc):: reqsend,reqrec
     &       , req2send,req2rec
      integer npoleloc_e
c
      parameter (nrhs=2)
      real(t_p)  udsum, upsum
c
      external tmatxb_shortreal
c
      if (.not.use_polar) return

#ifdef USE_NVSHMEM_CUDA
      if (use_pmecore.and.rank.ge.ndir) then
         call prmem_request(mu,3,nrhs,max(npolebloc,1),async=.true.)
         goto 200
      end if
#else
      if (use_pmecore.and.rank.ge.ndir) return
#endif

c
 1000 format(' illegal polalg in newinduce.')
c1010 format(' time for the ',a,F14.5)
c1020 format(' total elapsed time in newinduce: ',F14.5)
 
c
c     allocate some memory and clear the arrays:
c
      call prmem_request(mu,3,nrhs,max(npolebloc,1),async=.true.)
      call prmem_request(ef,3,nrhs,max(npolebloc,1),async=.true.)
c
      def_queue = dir_queue

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) call start_dir_stream_cover
#endif

      call set_to_zero2(mu,ef,3*nrhs*npolebloc,def_queue)
c
c     compute the electric fields:
c
      call timer_enter( timer_realdip )
      call efld0_directgpu_p1(nrhs,ef)
      call timer_exit( timer_realdip,quiet_timers )
c
      call commfieldshort(nrhs,ef)
c
      call commdirdirshort(nrhs,0,mu,reqrec,reqsend)
c
c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
 200  continue
      if (use_pred.and. nualt.eq.maxualt) then
         call pred_InitField(mu)
      end if
#ifdef USE_NVSHMEM_CUDA
      if (use_pmecore.and.rank.ge.ndir) goto 300
#endif

      if (.not.(use_pred.and.nualt.eq.maxualt)) then
      if (polgsf.eq.0) then
!$acc parallel loop collapse(3) async(def_queue) default(present)
        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
          iipole = poleglob(i)
              mu(j,k,i) = polarity(iipole)*ef(j,k,i)
            end do
          end do
        end do
      else
!$acc parallel loop collapse(2) async(def_queue) default(present)
        do i = 1, npoleloc
          do j = 1, 3
            iipole = poleglob(i)
            mu(j,1,i) = uind(j,iipole)
            mu(j,2,i) = uinp(j,iipole)
          end do
        end do
      end if
      end if
c
      call commdirdirshort(nrhs,1,mu,reqrec,reqsend)
      call commdirdirshort(nrhs,2,mu,reqrec,reqsend)
c
c     now, call the proper solver.
c
      if (polalg.eq.pcg_SId) then
         call inducepcg_shortrealgpu(tmatxb_p1,nrhs,.true.,ef,mu)
      else if (polalg.eq.jacobi_SId) then
         call inducejac_shortrealgpu(tmatxb_p1,nrhs,.true.,ef,mu)
#if 0
      else if (polalg.eq.step_pcg_SId.or.polalg.eq.step_pcg_short_SId)
     &   then
         call inducestepcg_pme2gpu(tmatxb_p1,nrhs,.true.,ef,mu,murec)
#endif
      else
         if (rank.eq.0) write(iout,1000) 
         call fatal
      end if
c
!$acc parallel loop collapse(2) async(def_queue) default(present)
      do i = 1, npolebloc
         do j = 1, 3
            iipole = poleglob(i)
            uind(j,iipole) = mu(j,1,i)
            uinp(j,iipole) = mu(j,2,i)
         end do
      end do

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif
c
c     update the lists of previous induced dipole values
c
  300 continue
      if (use_pred) call pred_SetAlt

#ifdef USE_NVSHMEM_CUDA
      if (use_pmecore.and.rank.ge.ndir) return
#endif

      end
c
      module inducepcg_s_mod
      logical   :: first_in=.true.
      integer   :: precnd1, ips_call=0
      integer   :: lastIter=0,sameIter
      logical   :: exitlast=.false.
      real(t_p) :: ggold1,ggold2
      real(r_p) :: gg1,gg2,ggnew1,ggnew2,ene1,ene2
      real(t_p),target :: gbuff(4)
      real(r_p),target :: gbuff1(6)
      real(r_p),pointer:: ggnew(:),ene(:),gnorm(:)
      !real(t_p)    ggdev1,ggdev2
      !real(t_p),save:: alphacg1,alphacg2
      parameter(precnd1=0)
      end module
      subroutine inducepcg_shortrealgpu(matvec,nrhs,precnd,ef,mu)
      use atmlst
      use domdec
      use ewald
      use inform     ,only: deb_Path,abort,minmaxone,app_id,dynamic_a
      use inducepcg_s_mod
      use interfaces ,only: tmatxb_pmegpu
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
      procedure(tmatxb_pmegpu) :: matvec
c
      integer   i, it, j, k, npoleloc_e
      real(t_p) zero, one, reg1, reg2
      real(r_p) zerom, pt5, resnrm
      parameter(zero=0.0, zerom=0, pt5=0.5, one=1.0)
c
c     MPI
c
      integer iglob, iipole, ierr, commloc
      integer req1, req2, req3, req4
      integer status(MPI_STATUS_SIZE)
      integer,dimension(nproc) ::reqrec,reqsend,req2rec,req2send
     &       , reqendrec,reqendsend,req2endrec,req2endsend
      logical tinker_isnan_m,skip_chk
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
 1050 format(' Short Induce PCG solver not converged after',I6,' steps '
     &    ,/,' Residual norm ', d10.5)
c
      if (deb_Path) 
     &   write(*,'(3x,a)') 'inducepcg_shortrealgpu'

      ips_call   = ips_call + 1
      commloc    = merge(comm_dir,COMM_TINKER,use_pmecore)
      npoleloc_e = merge(npolebloc,npoleloc,(skpPcomm))
      gbuff1     =  0.0_re_p
      def_queue  = dir_queue

      if (ips_call.eq.1) then
!$acc enter data create(ggold1,ggold2,gg1,gg2,ggnew1,ggnew2
!$acc&                 ,ene1,ene2,gbuff,gbuff1)
         !call Tinker_shellEnv("PRECND",precnd1,0)
         ggnew(1:2) => gbuff1(1:2)
         ene  (1:2) => gbuff1(3:4)
         gnorm(1:2) => gbuff1(5:6)
      end if

      if (precnd1) call polarEingenVal
c
c     allocate some memory and setup the preconditioner:
c
      call prmem_request(zr,3,nrhs,max(npoleloc_e,1),
     &     async=.true.)
      call prmem_request(res,3,nrhs,max(npoleloc_e,1),
     &     async=.true.)
      call prmem_request(h ,3,nrhs,max(npolebloc,1),
     &     async=.true.)
      call prmem_request(pp ,3,nrhs,max(npolebloc,1),
     &     async=.true.)
      call prmem_request(diag,npoleloc_e,async=.true.)

      if (precnd) then
!$acc parallel loop async(def_queue) default(present)
         do i = 1, npoleloc_e
            iipole  = poleglob(i)
            reg1    = polarity(iipole)
            diag(i) = merge(tinypol,reg1,reg1.eq.0.0_ti_p)
         end do
         if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
!$acc parallel loop async(def_queue) default(present)
         do i = 1, npoleloc_e
            diag(i) = one
         end do
      end if
c
c     initialize
c
      call set_to_zero2(res,zr,3*nrhs*npoleloc_e,def_queue)
      call set_to_zero2( pp, h,3*nrhs*npolebloc ,def_queue)
c
c     now, compute the initial direction
c
!$acc serial async(def_queue) present(ggold1,ggold2)
      ggold1 = zero
      ggold2 = zero
!$acc end serial
c
      if (precnd1) then
         call projectorOpe(mu,mu,0)
         call invertDefaultQmat(ef,mu,1.0_ti_p)
      end if
c
      call timer_enter( timer_realdip )
      call matvec(nrhs,.true.,mu,h)
      call timer_exit ( timer_realdip,quiet_timers )
      call commfieldshort(nrhs,h)
c
      !We do exit this sub-routine if (nproc.eq.1)
c
      call timer_enter ( timer_other )
!$acc parallel loop collapse(3) async(def_queue) default(present)
!$acc&         present(ggold1,ggold2) reduction(+:ggold1,ggold2)
      do k=1,npoleloc; do j=1,nrhs; do i=1,3
         res(i,j,k) = ef(i,j,k) - h(i,j,k)
         zr(i,j,k)  = res(i,j,k)*diag(k)
         pp(i,j,k)  = zr(i,j,k)
         if (btest(j,0)) then
            ggold1  = ggold1 + res(i,j,k)*zr(i,j,k)
         else
            ggold2  = ggold2 + res(i,j,k)*zr(i,j,k)
         end if
      end do; end do; end do
c
      if (skpPcomm) then
!$acc parallel loop collapse(3) async(def_queue) default(present)
      do k = npoleloc+1,npolebloc; do j = 1,nrhs; do i = 1,3
         reg1  = ef(i,j,k) - h(i,j,k)
         reg2  = reg1*diag(k)
         res(i,j,k) = reg1
          zr(i,j,k) = reg2
          pp(i,j,k) = reg2
      end do; end do; end do
      end if
      if (precnd1) call projectorOpe(pp,pp,0)

      call timer_exit( timer_other,quiet_timers )

      if (nproc.ne.1) then
!$acc serial async present(gbuff,ggold1,ggold2)
         gbuff(1) = ggold1
         gbuff(2) = ggold2
!$acc end serial
!$acc wait(def_queue)
!$acc host_data use_device(gbuff)
         call MPI_ALLREDUCE(MPI_IN_PLACE,gbuff,nrhs,MPI_TPREC
     $                     ,MPI_SUM,commloc,ierr)
         !call MPI_WAIT(req1,status,ierr)
!$acc end host_data
!$acc serial async present(ggold1,ggold2,gbuff)
         ggold1 = gbuff(1)
         ggold2 = gbuff(2)
!$acc end serial
      end if
c
c     now, start the main loop:
c
      do it = 1, politer
c
        if (.not.skpPcomm) then
        ! Direct space electrical field comm to neighbors
        call commdirdirshort(nrhs,0,pp,reqrec,reqsend)
        call commdirdirshort(nrhs,1,pp,reqrec,reqsend)
        call commdirdirshort(nrhs,2,pp,reqrec,reqsend)
        end if
c
        call timer_enter( timer_realdip )
        call matvec(nrhs,.true.,pp,h)
        call timer_exit ( timer_realdip,quiet_timers )

        call commfieldshort(nrhs,h)
!$acc serial async(def_queue)
!$acc&       present(gg1,gg2,ggnew1,ggnew2,ene1,ene2)
        gg1    = zerom; gg2    = zerom;
        ggnew1 = zerom; ggnew2 = zerom;
        ene1   = zerom; ene2   = zerom;
!$acc end serial

        call timer_enter( timer_other )
!$acc parallel loop collapse(3) async(def_queue)
!$acc&         present(gg1,gg2)
        do i = 1,npoleloc; do j = 1,nrhs; do k = 1, 3
           if (btest(j,0)) then
              gg1 = gg1 + pp(k,j,i)*h(k,j,i)
           else
              gg2 = gg2 + pp(k,j,i)*h(k,j,i)
           end if
        end do; end do; end do

        if (nproc.ne.1) then
!$acc wait(def_queue)
!$acc host_data use_device(gg1,gg2)
           call MPI_ALLREDUCE(MPI_IN_PLACE,gg1,1,MPI_RPREC,
     $          MPI_SUM,commloc,ierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE,gg2,1,MPI_RPREC,
     $          MPI_SUM,commloc,ierr)
!$acc end host_data
        end if

        if (deb_Path) then
!$acc wait(def_queue)
!$acc update host(gg1,gg2)
          if (gg1.eq.zero) goto 30
          if (gg2.eq.zero) goto 30
        end if

!$acc parallel loop collapse(3) async(def_queue)
!$acc&         default(present) 
!$acc&         present(gg1,gg2,ggold1,ggold2,ggnew1,ggnew2,ene1,ene2)
!$acc&         reduction(+:ggnew1,ggnew2,ene1,ene2)
        do k=1,npoleloc; do j=1,nrhs; do i=1,3
           if (btest(j,0)) then
              mu (i,j,k) = mu (i,j,k) + (ggold1/gg1)*pp(i,j,k)
              res(i,j,k) = res(i,j,k) - (ggold1/gg1)*h (i,j,k)
              zr (i,j,k) = diag(k)* res(i,j,k)
              ggnew1     = ggnew1 + res(i,j,k)*zr(i,j,k)
              ene1    = ene1 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
            else
              mu (i,j,k) = mu (i,j,k) + (ggold2/gg2)*pp(i,j,k)
              res(i,j,k) = res(i,j,k) - (ggold2/gg2)*h (i,j,k)
              zr (i,j,k) = diag(k)* res(i,j,k)
              ggnew2     = ggnew2 + res(i,j,k)*zr(i,j,k)
              ene2    = ene2 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
           end if
        end do; end do; end do
        if (skpPcomm) then
!$acc parallel loop collapse(3) async(def_queue)
!$acc&         default(present) present(gg1,gg2,ggold1,ggold2)
        do k=npoleloc+1,npoleloc_e; do j=1,nrhs; do i=1,3
           if (btest(j,0)) then
              mu (i,j,k) = mu (i,j,k) + (ggold1/gg1)*pp(i,j,k)
              res(i,j,k) = res(i,j,k) - (ggold1/gg1)*h (i,j,k)
              zr (i,j,k) = diag(k)* res(i,j,k)
            else
              mu (i,j,k) = mu (i,j,k) + (ggold2/gg2)*pp(i,j,k)
              res(i,j,k) = res(i,j,k) - (ggold2/gg2)*h (i,j,k)
              zr (i,j,k) = diag(k)* res(i,j,k)
           end if
        end do; end do; end do
        end if

        if (nproc.ne.1) then
!$acc serial async(def_queue) present(ene1,ene2,ggnew1,ggnew2,gbuff1)
           gbuff1(1)=ggnew1; gbuff1(2)=ggnew2
           gbuff1(3)=  ene1; gbuff1(4)=  ene2
!$acc end serial
!$acc wait(def_queue)
!$acc host_data use_device(gbuff1)
           call MPI_ALLREDUCE(MPI_IN_PLACE,gbuff1(1),2*nrhs,MPI_RPREC
     $                       ,MPI_SUM,commloc,ierr)
           !call MPI_WAIT(req3,status,ierr)
!$acc end host_data
!$acc serial async(def_queue) present(ggnew1,ggnew2,ene1,ene2,gbuff1)
           ggnew1 = gbuff1(1); ggnew2= gbuff1(2)
           ene1   = gbuff1(3); ene2  = gbuff1(4)
!$acc end serial
        end if

        !ggdev1 = ggnew(1)/ggold(1); ggdev2 = ggnew(2)/ggold(2)
        ! Compute Next direction
!$acc host_data use_device(pp,zr)
!$acc parallel loop collapse(3) async(def_queue)
!$acc&         present(ggold1,ggold2,ggnew1,ggnew2) deviceptr(pp,zr)
        do k = 1,npoleloc_e; do j = 1,nrhs; do i = 1, 3
           if (btest(j,0)) then
             pp(i,j,k) = zr(i,j,k) + (ggnew1/ggold1)*pp(i,j,k)
           else
             pp(i,j,k) = zr(i,j,k) + (ggnew2/ggold2)*pp(i,j,k)
           end if
        end do; end do; end do
!$acc end host_data
        if (precnd1) call projectorOpe(pp,pp,0)

!$acc serial async(def_queue)
!$acc&       present(ggold1,ggold2,ggnew1,ggnew2,ene1,ene2,gbuff1)
        ggold1 = ggnew1
        ggold2 = ggnew2
        if (nproc.eq.1) then
           gbuff1(1)= ggnew1; gbuff1(2)= ggnew2
           gbuff1(3)= ene1  ; gbuff1(4)= ene2  
        end if
!$acc end serial
        call timer_exit( timer_other,quiet_timers )

         ! Skip Iteration check
         skip_chk = it.lt.2.or.(exitlast.and.(it.lt.lastIter))
         ! Skip last Iteration check
         if (exitlast.and.it.eq.lastIter.and.mod(sameIter,3).ne.0)
     &      goto 10

        if (skip_chk) then
           resnrm = 2*poleps
        else
!$acc update host(gbuff1) async(def_queue)
!$acc wait(def_queue)
           do k = 1, nrhs
             gnorm(k) = sqrt(ggnew(k)/real(3*npolar,r_p))
           end do
           resnrm = max(gnorm(1),gnorm(2))
        end if

        ene = -pt5*ene
        if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     $    it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)
c
        if ((it.eq.politer.and.resnrm.gt.poleps)
     &     .or.tinker_isnan_m(gnorm(1))) then
!$acc update host(gbuff)
           ! Not converged Abort
           write(0,1050) it,max(gnorm(1),gnorm(2))
           write(0,'(6f12.7)') gbuff1(1:4),gbuff(1:2)
           call minmaxone(ef,size(ef),'ef')
           call minmaxone(mu,size(mu),'mu')
           abort = .true.
           goto 10  !Exit Loop
        end if
        if (resnrm.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     $      (ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
          goto 10   !Exit Loop
        end if
      end do
 10   continue

      ! Save and Check Convergence Iteration Index
      if (app_id.eq.dynamic_a) then
      sameIter = merge(0,sameIter+1,lastIter.ne.it)
      lastIter = it
      exitlast = polprt.eq.0.and.sameIter.gt.25
     &           .and.mod(ips_call,100).ne.0
      if (deb_Path) print*, 'pcg_conv iter',lastIter,sameIter,exitlast
      end if

      ! Solution field comm to direct neighbors
      call commdirdirshort(nrhs,0,mu,reqendrec,reqendsend) ! Reception chanels
      call commdirdirshort(nrhs,1,mu,reqendrec,reqendsend) ! Send solution
      call commdirdirshort(nrhs,2,mu,reqendrec,reqendsend) ! Wait comm
c
c     finalize and debug printing:
c
      if (polprt.ge.3) then
        write(iout,1020)
!$acc wait
!$acc update host(poleglob,mu)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = ipole(iipole)
          write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
!$acc wait
!$acc update host(poleglob,mu)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = ipole(iipole)
          write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
        end do
      end if

  30  continue
      end
c
      subroutine inducejac_shortrealgpu(matvec,nrhs,dodiis,ef,mu)
      use atmlst
      use domdec
      use ewald
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent
      use timestat
      use tinheader ,only:ti_p,re_p
      use mpi
      implicit none
c
c     solves the polarizable dipoles linear equations by jacobi
c     iterations coupled to the DIIS extrapolation.
c
      integer nrhs, info
      real(t_p)  ef(3,nrhs,*), mu(3,nrhs,*)
      real(t_p), allocatable :: munew(:,:,:), h(:,:,:)
      real(t_p), allocatable :: xdiis(:,:), ediis(:,:), bmat(:,:),
     $   bloc(:,:), cex(:)
      integer i, j, k, it, ind, ndismx, nmat, lenb
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real(t_p)  zero, one, rnorm(2), rr, xx(1)

      real(t_p) time1,time2
      save    zero, one, xx
      external matvec
c
c     MPI
c
      integer iglob, iipole, ierr
      integer reqnorm,reqdiis(2*ndismx+1)
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
c
 1000 format(' itsolv shortreal converged after ',I3,' iterations.',/,
     $       ' final residual norm = ',3D14.7)
 1010 format(' residual norm at iteration ',I3,':',3D12.2)
 1020 format(' Jacobi/DIIS solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Jacobi/DIIS solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
c
c
      zero  = 0.0_ti_p
      one   = 1.0_ti_p
      xx(1) = 0.0_ti_p

      allocate (munew(3,nrhs,max(1,npolebloc)))
      munew = 0_ti_p
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
c
      allocate (h(3,nrhs,max(1,npolebloc)))
      h = 0_ti_p
      if (dodiis) then
        nmat = 1
        lenb = ndismx + 1
        allocate (xdiis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (ediis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (bmat(lenb,lenb))
        allocate (bloc(lenb,lenb))
        allocate (cex(lenb))
        bmat = 0_ti_p
      end if
c
c     main loop:
c
      do it = 1, politer
        h = 0_ti_p
        rnorm = 0_ti_p
c
        call timer_enter( timer_realdip )
        call matvec(nrhs,.false.,mu,h)
        call timer_exit( timer_realdip,quiet_timers )
        call commfieldshort(nrhs,h)
c
        call commdirdirshort(nrhs,0,mu,reqrec,reqsend)
c
c     jacobi step:
c
        do i = 1, npoleloc
          iipole = poleglob(i)
          do k = 1, nrhs
            do j = 1, 3
              munew(j,k,i) = polarity(iipole)*(ef(j,k,i) - h(j,k,i))
            end do
          end do
        end do
c
        do k = 1, nrhs
          rnorm(k) = sum((munew(1:3,k,1:npoleloc)-
     $      mu(1:3,k,1:npoleloc))**2)
        end do
c
        call MPI_IALLREDUCE(MPI_IN_PLACE,rnorm(1),nrhs,MPI_TPREC,
     $    MPI_SUM,COMM_TINKER,reqnorm,ierr)
        if (dodiis) then
          ind = 0
          do i = 1, npoleloc
            do k = 1, nrhs
              do j = 1, 3
                ind = ind + 1
                xdiis(ind,nmat) = munew(j,k,i)
                ediis(ind,nmat) = munew(j,k,i) - mu(j,k,i)
              end do
            end do
          end do
c
c         Compute Pulay's Matrix and extrapolate
c
          call diis(ndismx,3*nrhs*npoleloc,xdiis,ediis,
     $      bmat,nmat,reqdiis,COMM_TINKER)
c
          do i = 1, 2*nmat-3
            call MPI_WAIT(reqdiis(i),status,ierr)
          end do
          bloc = bmat
          cex = 0_ti_p
          cex(1) = one
#ifdef _OPENACC
          call fatal_device("inducejac_shortrealgpu")
#else
          call M_gesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#endif
          munew = 0_ti_p
          call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
        end if
c
        call commdirdirshort(nrhs,1,munew,reqrec,reqsend)
        call commdirdirshort(nrhs,2,mu,reqrec,reqsend)
        mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)
c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/dble(3*npolar))
          rr = max(rnorm(k),rr)
        end do
        if (rr.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1000) it, (rnorm(k), k = 1, nrhs)
          goto 10
        end if
        if (polprt.ge.1.and.rank.eq.0)
     $     write(6,1010) it, (rnorm(k), k = 1, nrhs)
      end do

 10   continue
      if (polprt.ge.3) then
        write(iout,1020)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = ipole(iipole)
          write(iout,1030) iglob, (mu(j,1,i), j = 1, 3)
        end do
      end if
      if (polprt.ge.4) then
        write(iout,1021)
        do i = 1, npoleloc
          iipole = poleglob(i)
          iglob = poleglob(iipole)
          write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
        end do
      end if
c
c     free the memory.
c
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (munew)
      deallocate (h)
      if (dodiis) then
        deallocate (xdiis)
        deallocate (ediis)
        deallocate (bmat)
        deallocate (bloc)
        deallocate (cex)
      end if
      end
