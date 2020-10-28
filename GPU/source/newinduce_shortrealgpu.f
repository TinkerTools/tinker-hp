 
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
      use interfaces,only:inducepcg_shortrealgpu,tmatxb_p,
     &                    tmatxb_pmevec,efld0_directgpu2,
     &                    efld0_directgpu_p
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      use timestat
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
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
c
      parameter (nrhs=2)
      real(t_p)  udsum, upsum

      real(t_p), allocatable :: ef(:,:,:), mu(:,:,:)
c
      external tmatxb_shortreal
c
      if (.not.use_polar) return

#ifdef USE_NVSHMEM_CUDA
      if (use_pmecore.and.rank.ge.ndir) then
         allocate (mu(3,nrhs,max(1,npolebloc)))
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
      allocate (mu(3,nrhs,max(1,npolebloc)))
      allocate (ef(3,nrhs,max(1,npolebloc)))
c
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))

      def_queue = dir_queue

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
         ! Start stream overlapping
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

!$acc enter data create(mu,ef) async(def_queue)

      call set_to_zero2(mu,ef,3*nrhs*npolebloc,def_queue)
c
c     compute the electric fields:
c
      call timer_enter( timer_realdip )
      call efld0_directgpu_p(nrhs,ef)
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

#ifdef USE_NVSHMEM_CUDA
!$acc data present(poleglob,uind,uinp)
#else
!$acc data present(poleglob,udshortalt,upshortalt,uind,uinp)
#endif
      if (.not.(use_pred.and.nualt.eq.maxualt)) then
      if (polgsf.eq.0) then
!$acc parallel loop collapse(3) async(def_queue)
        do i = 1, npoleloc
          do k = 1, nrhs
            do j = 1, 3
          iipole = poleglob(i)
              mu(j,k,i) = polarity(iipole)*ef(j,k,i)
            end do
          end do
        end do
      else
!$acc parallel loop collapse(2) async(def_queue)
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
      if (polalg.eq.1) then
         call inducepcg_shortrealgpu(tmatxb_p,nrhs,.true.,ef,mu)
      else if (polalg.eq.2) then
         call inducejac_shortrealgpu(tmatxb_p,nrhs,.true.,ef,mu)
      else
         if (rank.eq.0) write(iout,1000) 
         call fatal
      end if
c
!$acc parallel loop collapse(2) async(def_queue)
      do i = 1, npolebloc
         do j = 1, 3
            iipole = poleglob(i)
            uind(j,iipole) = mu(j,1,i)
            uinp(j,iipole) = mu(j,2,i)
         end do
      end do
!$acc end data

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
         ! End stream overlapping
         call stream_wait_async(dir_stream,rec_stream,rec_event)
      end if
#endif
c
c     update the lists of previous induced dipole values
c
  300 continue
      if (use_pred) call pred_SetAlt

#ifdef USE_NVSHMEM_CUDA
      if (use_pmecore.and.rank.ge.ndir) then
         deallocate(mu)
         return
      end if
#endif

      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
!$acc exit data delete(mu,ef) async(def_queue)
      deallocate (ef)
      deallocate (mu)
      end
c
      subroutine inducepcg_shortrealgpu(matvec,nrhs,precnd,ef,mu)
      use atmlst
      use domdec
      use ewald
      use inform     ,only: deb_Path,abort
      use interfaces ,only: tmatxb_pmegpu
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent     ,only: use_pmecore
      use timestat
      use units
      use utils
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

      integer i, it, j, k
      real(t_p) ggold(2), alphacg(2)
      real(r_p) gnorm(2), ggnew(2), ene(2), gg(2)
      real(t_p),save:: ggold1,ggold2
      real(r_p),save:: ggnew1,ggnew2
      real(t_p),save:: ggdev1,ggdev2
      real(r_p),save:: gnorm1,gnorm2
      real(r_p),save:: gg1,gg2
      real(r_p),save:: ene1,ene2
      real(t_p),save:: alphacg1,alphacg2
      real(t_p) zero, one
      real(r_p) zerom, pt5
      real(r_p) resnrm
      real(8) time1,time2
      real(t_p),allocatable :: res(:,:,:),h(:,:,:),pp(:,:,:),
     $          zr(:,:,:),diag(:)
      parameter (zero=0.0_ti_p,zerom=0, pt5=0.50_re_p, one=1.0_ti_p)
c
c     MPI
c
      integer iglob, iipole, ierr, commloc
      integer req1, req2, req3, req4
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer,allocatable :: reqendrec(:),reqendsend(:)
      integer,allocatable :: req2endrec(:),req2endsend(:)
      logical,save :: first_in=.true.
      logical tinker_isnan_m
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

      if (use_pmecore) then
         commloc = comm_dir
      else
         commloc = COMM_TINKER
      end if

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
      allocate (res(3,nrhs,max(1,npoleloc)))
      allocate (h(3,nrhs,max(1,npolebloc)))
      allocate (pp(3,nrhs,max(1,npolebloc)))
      allocate (zr(3,nrhs,max(1,npoleloc)))
      allocate (diag(npoleloc))

      if (first_in) then
!$acc enter data create(gg1,gg2) async(def_queue)
         first_in=.false.
      end if

!!$acc enter data create(res,h,pp,zr,diag) async(def_queue)
!$acc data create(res,h,pp,zr,diag)
!$acc&     present(mu,ef,poleglob,polarity,ipole)
!$acc&     async(def_queue)

      if (precnd) then
!$acc parallel loop async(def_queue)
         do i = 1, npoleloc
            iipole  = poleglob(i)
            diag(i) = polarity(iipole)
         end do
         if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
!$acc parallel loop async(dir_queue)
         do i = 1, npoleloc
            diag(i) = one
         end do
      end if
c
c     initialize
c
      call set_to_zero2(res,zr,3*nrhs* npoleloc,def_queue)
      call set_to_zero2( pp, h,3*nrhs*npolebloc,def_queue)
c
c     now, compute the initial direction
c
      ggold  = zero
      ggold1 = zero
      ggold2 = zero
c
      call timer_enter( timer_realdip )
      call matvec(nrhs,.true.,mu,h)
      call timer_exit ( timer_realdip,quiet_timers )
      call commfieldshort(nrhs,h)
c
      !We do exit this sub-routine if (nproc.eq.1)
      call commdirdirshort(nrhs,0,pp,reqrec,reqsend)
c
      call timer_enter ( timer_other )
      !FIXME This implicit reduction is undetected by PGI-20.4 Compiler
      !  Why ? ( Probably a Bug ? )
!$acc parallel loop collapse(3) async(def_queue)
!$acc&         reduction(+:ggold1,ggold2)
      do k=1,npoleloc
         do j=1,nrhs
            do i=1,3
               res(i,j,k) = ef(i,j,k) - h(i,j,k)
               zr(i,j,k)  = res(i,j,k)*diag(k)
               pp(i,j,k)  = zr(i,j,k)
               if (btest(j,0)) then
                  ggold1  = ggold1 + res(i,j,k)*zr(i,j,k)
               else
                  ggold2  = ggold2 + res(i,j,k)*zr(i,j,k)
               end if
            end do
         end do
      end do
!$acc wait(def_queue)

      ggold(1) = ggold1
      ggold(2) = ggold2
      call timer_exit ( timer_other,quiet_timers )

      if (nproc.ne.1) then
!$acc wait
         call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_TPREC,
     $        MPI_SUM,commloc,req1,ierr)
c
c     MPI : begin sending
c
         call commdirdirshort(nrhs,1,pp,reqrec,reqsend)
         call commdirdirshort(nrhs,2,mu,reqrec,reqsend)
         call MPI_WAIT(req1,status,ierr)
         ggold1 = ggold(1)
         ggold2 = ggold(2)
      end if
c
c     now, start the main loop:
c
      do it = 1, politer
c
        call timer_enter( timer_realdip )
        call matvec(nrhs,.true.,pp,h)
        call timer_exit ( timer_realdip,quiet_timers )

        call commfieldshort(nrhs,h)
c
c     MPI : begin reception
c
        call commdirdirshort(nrhs,0,pp,reqrec,reqsend)
!$acc serial async(def_queue) present(gg1,gg2)
        gg1=zerom; gg2=zerom
!$acc end serial

        call timer_enter( timer_other )
!$acc parallel loop collapse(3) present(gg1,gg2) async(def_queue)
        do i = 1,npoleloc
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

        if (nproc.ne.1) then
           call MPI_IALLREDUCE(MPI_IN_PLACE,gg(1),nrhs,MPI_RPREC,
     $          MPI_SUM,commloc,req2,ierr)
           call MPI_WAIT(req2,status,ierr)
           gg1 = gg(1); gg2 = gg(2)
!$acc update device(gg1,gg2) async(def_queue)
        end if
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
         do k=1,npoleloc
            do j=1,nrhs
               do i=1,3
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
               end do
            end do
         end do
c Implicit wait seems to be removed with PGI-20
!$acc wait(def_queue)
         ene2 = -pt5*ene2; ene1 = -pt5*ene1
         ggnew(1)=ggnew1; ggnew(2)=ggnew2
           ene(1)=  ene1;   ene(2)=  ene2

        if (nproc.ne.1) then
           call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_RPREC,
     $          MPI_SUM,commloc,req3,ierr)
           call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_RPREC,
     $          MPI_SUM,commloc,req4,ierr)
           call MPI_WAIT(req3,status,ierr)
           call MPI_WAIT(req4,status,ierr)
           ggnew1 = ggnew(1); ggnew2=ggnew(2)
        end if
        resnrm = zero
        do k = 1, nrhs
          gnorm(k) = sqrt(ggnew(k)/real(3*npolar,r_p))
          resnrm   = max(resnrm,gnorm(k))
        end do
        if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     $    it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)

        ggdev1 = ggnew(1)/ggold(1); ggdev2 = ggnew(2)/ggold(2)
!$acc parallel loop collapse(3) async(def_queue)
        do k = 1,npoleloc
           do j = 1,nrhs
              do i = 1, 3
               if (btest(j,0)) then
                 pp(i,j,k) = zr(i,j,k) + ggdev1*pp(i,j,k)
               else
                 pp(i,j,k) = zr(i,j,k) + ggdev2*pp(i,j,k)
               end if
              end do
           end do
        end do
        call timer_exit( timer_other,quiet_timers )
c
        call commdirdirshort(nrhs,1,pp,reqrec,reqsend)
c
        call commdirdirshort(nrhs,2,pp,reqrec,reqsend)
c
        ggold = ggnew
        ggold1 = ggold(1); ggold2 = ggold(2)
        if (resnrm.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     $      (ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
          goto 10
        end if
        if ((it.eq.politer.and.resnrm.gt.poleps)
     &     .or.tinker_isnan_m(resnrm)) then
           ! Not converged Exit
           write(0,1050) it,resnrm
           abort = .true.
        end if
      end do
 10   continue
c
c     MPI : begin reception
c
      call commdirdirshort(nrhs,0,mu,reqendrec,reqendsend)
c
c     MPI : begin sending
c
      call commdirdirshort(nrhs,1,mu,reqendrec,reqendsend)
c
c
      call commdirdirshort(nrhs,2,mu,reqendrec,reqendsend)
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
!$acc end data
!!$acc exit data delete(res,h,pp,zr,diag) async(def_queue)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqendsend)
      deallocate (reqendrec)
      deallocate (req2endsend)
      deallocate (req2endrec)
      deallocate (res)
      deallocate (h)
      deallocate (pp)
      deallocate (zr)
      deallocate (diag)
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
          call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
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
