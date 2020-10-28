c      
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
      subroutine newinduce_pme2gpu
      use atoms     ,only: n,n_pe
      use atmlst
      use domdec
      use ewald
      use iounit
      use inform    ,only: deb_Path
      use interfaces,only:inducepcg_pme2gpu,tmatxb_p,
     &                    tmatxb_pmevec,efld0_directgpu2,
     &                    efld0_directgpu_p
      use math
      use mpole
      use nvshmem
      use pme
      use polar
      use polpot
      use potent
      use units
      use uprior
      use utils
      use utilcomm ,buffermpi1=>buffermpi2d,buffermpi2=>buffermpi2d1
      use utilgpu
      use timestat
      use mpi
      implicit none
c
c     without separate cores for reciprocal part
c
      integer i, j, k, nrhs, proc
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
      integer npl
#endif
c
c     MPI
c
      integer iglob, ierr, iipole
      integer, allocatable :: reqrecdirsend(:),reqrecdirrec(:)
      integer, allocatable :: reqsend(:),reqrec(:)
      integer, allocatable :: req2send(:),req2rec(:)
      integer tag,status(MPI_STATUS_SIZE)
      parameter (nrhs=2)
c
      real(t_p),parameter :: zero=0.0_ti_p
      real(8) wtime0, wtime1, wtime2
      real(t_p) term, xx(1), comp
      real(t_p),save, allocatable :: ef(:,:,:), mu(:,:,:),
     &                          murec(:,:,:)
      real(t_p),save, allocatable :: cphi(:,:)
c
c     real(t_p), allocatable ::buffermpi1(:,:),buffermpi2(:,:),
c    &                buffermpimu1(:,:,:),buffermpimu2(:,:,:)
c
      if (deb_Path) 
     &   write(*,'(2x,a)') 'newinduce_pme2gpu'
      if (.not.use_polar) return
c
 1000 format(' illegal polalg in newinduce.')
 1010 format(' time for the ',a,F14.5)
 1020 format(' total elapsed time in newinduce: ',F14.5)
 
c
c     allocate some memory and clear the arrays:
c
      call timer_enter( timer_other )
c     allocate (mu   (3,nrhs,max(1,npolebloc)))
c     allocate (murec(3,nrhs,max(1,npolerecloc)))
c     allocate (ef   (3,nrhs,max(1,npolebloc)))
c
      call prmem_request(mu,3,nrhs,max(npolebloc,1),async=.true.)
      call prmem_request(murec,3,nrhs,max(npolerecloc,1),
     &     async=.true.)
      call prmem_request(ef,3,nrhs,max(npolebloc,1),async=.true.)
      call prmem_request(cphi,10,npoleloc,async=.true.)

      if (nproc.ge.1) then
      call prmem_request(buffermpimu1,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(buffermpimu2,3,nrhs,max(npolerecloc,1),
     &     async=.true.)
      call prmem_request(buffermpi1,10,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(buffermpi2,10,max(npolerecloc,1),
     &     async=.true.)
      end if
c
c     allocate (cphi(10,max(npoleloc,1)))

      if (.not.use_mpole) then
         j = max(npolerecloc,1)
         call prmem_request(cphirec,10,j,async=.true.)
         call prmem_request(fphirec,20,j,async=.true.)
      end if

      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))

!$acc data present(ef,mu,murec,cphi)
!$acc&     present(cphirec,fphirec,rpole,uind,uinp,polarity,
#ifndef USE_NVSHMEM_CUDA
!$acc&   upalt,udalt,
#endif
!$acc&   bpred,poleloc,poleglob)

      ! set arrays to zero
      call set_to_zero2(mu,ef,3*nrhs*npolebloc,rec_queue)
      call set_to_zero1(fphirec,20*npolerecloc,rec_queue)
      if (nproc.eq.1) then
         call set_to_zero1(murec,3*nrhs*npolerecloc,rec_queue)
         call set_to_zero1(cphi   ,10*npoleloc   ,rec_queue)
         call set_to_zero1(cphirec,10*npolerecloc,rec_queue)
      else
         call set_to_zero1(buffermpimu1,3*nrhs*npoleloc ,rec_queue)
         call set_to_zero2(murec,buffermpimu2,3*nrhs*npolerecloc,
     &                    rec_queue)
         call set_to_zero2(cphi   ,buffermpi1,10*npoleloc   ,rec_queue)
         call set_to_zero2(cphirec,buffermpi2,10*npolerecloc,rec_queue)
      end if
      call timer_exit (timer_other,quiet_timers)
c
c     compute the electric fields:
c
      wtime0 = mpi_wtime()
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)
#endif
c
c     compute the direct space contribution (fields)
c
      call timer_enter( timer_realdip )
      call efld0_directgpu_p(nrhs,ef)
      call timer_exit( timer_realdip,quiet_timers )
c
c     compute the reciprocal space contribution (fields)
c
      call timer_enter( timer_recdip )
      call efld0_recipgpu(cphi)
      call timer_exit( timer_recdip,quiet_timers )
c
      call commrecdirfieldsgpu(0,cphirec,cphi,buffermpi1,buffermpi2,
     &     reqrecdirrec,reqrecdirsend)
      call commrecdirfieldsgpu(1,cphirec,cphi,buffermpi1,buffermpi2,
     &     reqrecdirrec,reqrecdirsend)
      call commrecdirfieldsgpu(2,cphirec,cphi,buffermpi1,buffermpi2,
     &     reqrecdirrec,reqrecdirsend)
c
      call commfieldgpu(nrhs,ef)
      call commdirdirgpu(nrhs,0,mu,reqrec,reqsend)
c
c     add direct and reciprocal fields
c
      term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
      def_queue = rec_queue
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
         call stream_wait_async(dir_stream,rec_stream,dir_event)
      end if
#endif

!$acc parallel loop collapse(3) async(def_queue)
      do i = 1,npoleloc
         do j = 1,nrhs
            do k = 1,3
               iipole = poleglob(i)
               comp = term*rpole(k+1,iipole) - cphi(k+1,i)
               ef(k,j,i) = ef(k,j,i) + comp
            end do
         end do
      end do
      wtime1 = mpi_wtime()
c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
      if (use_pred .and. nualt.eq.maxualt) then
         call pred_InitField(mu)
      else if (polgsf.eq.0) then
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
!$acc parallel loop collapse(3) async(def_queue)
         do i = 1, npoleloc
            do j = 1,nrhs
               do k = 1, 3
                  iipole = poleglob(i)
                  if (btest(j,0)) then !j==1
                     mu(k,j,i) = uind(k,iipole)
                  else  !j/=1
                     mu(k,j,i) = uinp(k,iipole)
                  end if
               end do
            end do
         end do
      end if
      call commdirdirgpu(nrhs,1,mu,reqrec,reqsend)
      call commdirdirgpu(nrhs,2,mu,reqrec,reqsend)
c
      call commrecdirdipgpu(nrhs,0,murec,mu,buffermpimu1,buffermpimu2,
     &     req2rec,req2send)
      call commrecdirdipgpu(nrhs,1,murec,mu,buffermpimu1,buffermpimu2,
     &     req2rec,req2send)
      call commrecdirdipgpu(nrhs,2,murec,mu,buffermpimu1,buffermpimu2,
     &     req2rec,req2send)
c
c     now, call the proper solver.
c
      if (polalg.eq.1) then
         call inducepcg_pme2gpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
      else if (polalg.eq.2) then ! FIXME a porter
!$acc update host(ef,mu,murec)
         call inducejac_pme2gpu(tmatxb_pmevec,nrhs,.true.,ef,mu,murec)
!$acc update device(ef,mu,murec)
      else
         if (rank.eq.0) write(iout,1000)
         call fatal
      end if
      wtime2 = mpi_wtime()

      if (polprt.ge.1.and.rank.eq.0.and.tinkertime) then
         if (polprt.ge.2) then
            write(iout,1010) 'fields:  ', wtime1-wtime0
            write(iout,1010) 'dipoles: ', wtime2-wtime1
         end if
         write(iout,1020) wtime2 - wtime0
      end if
c
c     move the computed dipoles in the common block.
c
      call timer_enter( timer_other )
      def_queue = rec_queue
!$acc parallel loop collapse(2) async(def_queue)
      do i = 1, npolebloc
         do j = 1, 3
            iipole         = poleglob(i)
            uind(j,iipole) = mu(j,1,i)
            uinp(j,iipole) = mu(j,2,i)
         end do
      end do
!$acc parallel loop collapse(2) async(def_queue)
      do i = 1, npolerecloc
         do j = 1, 3
            iipole = polerecglob(i)
            iglob  = ipole(iipole)
            if (repart(iglob).ne.rank) then
               uind(j,iipole) = murec(j,1,i)
               uinp(j,iipole) = murec(j,2,i)
            else
               uind(j,iipole) = mu(j,1,poleloc(iipole))
               uinp(j,iipole) = mu(j,2,poleloc(iipole))
            end if
         end do
      end do
      call timer_exit( timer_other,quiet_timers )
c
c     update the lists of previous induced dipole values
c
      if (use_pred) call pred_SetAlt
!$acc end data

      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
      deallocate (reqrecdirrec)
      deallocate (reqrecdirsend)
c     deallocate (buffermpi1)
c     deallocate (buffermpi2)
c     deallocate (buffermpimu1)
c     deallocate (buffermpimu2)
c     deallocate (ef)
c     deallocate (fphi)
c     deallocate (mu)
c     deallocate (murec)
c     deallocate (cphi)
c     deallocate (cphirec)
      return
      end

      subroutine pred_InitField(mu)
      use atmlst,only: poleglob
      use atoms ,only: n,n_pe
      use domdec,only: rank,ndir
      use inform,only: deb_Path
      use mpi
      use mpole ,only: npoleloc,npolebloc
      use nvshmem
      use potent,only: use_pmecore
     &          ,use_polarshortreal
      use tinheader,only: ti_p
      use timestat
      use uprior
      use utilgpu  ,only: def_queue
      implicit none
      integer,parameter:: nrhs=2,d3=3
      real(t_p),intent(inout)::mu(d3,nrhs,npolebloc)
      integer i,j,k,iipole,ierr
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) udsum,upsum

#ifndef USE_NVSHMEM_CUDA
      if (use_pmecore.and.rank.ge.ndir) return
#endif
      if(deb_Path) print '(3x,A)',"pred_InitField"
      call timer_enter( timer_ulspred )
      call ulspredgpu

#ifdef USE_NVSHMEM_CUDA

      if (use_polarshortreal) then
         udalt_p1 => d_udshortalt
         upalt_p1 => d_upshortalt
      else
         udalt_p1 => d_udalt
         upalt_p1 => d_upalt
      end if

!$acc parallel loop gang vector async(def_queue)
!$acc&         present(bpred) deviceptr(udalt_p1,upalt_p1,d_altbuf)
      do i = 1, n_pe
!$acc loop seq
         do j = 1, 3
            udsum = 0.0_ti_p
            upsum = 0.0_ti_p
            do k = 1, nualt - 1
              udsum = udsum + bpred(k)*udalt_p1(mype)%pel(k,j,i)
              upsum = upsum + bpred(k)*upalt_p1(mype)%pel(k,j,i)
            end do
            d_altbuf(mype)%pel(j,1,i) = udsum
            d_altbuf(mype)%pel(j,2,i) = upsum
         end do
      end do

      ! Synchronise to avoid overwriting
      if(npes.ne.1) call mpi_barrier(COMM_NVSHMEM,ierr)
      if(.not.use_pmecore.or.(use_pmecore.and.rank.lt.ndir)) then

!$acc parallel loop async(def_queue) collapse(3)
!$acc&         present(mu)
!$acc&         deviceptr(d_altbuf)
      do i = 1, npoleloc
         do j = 1, nrhs
            do k = 1, d3
               iipole = poleglob(i)
               ipe    =    (iipole-1)/n_pe
               ind    = mod(iipole-1 ,n_pe)+1
               mu (k,j,i) = d_altbuf(ipe)%pel(k,j,ind)
            end do
         end do
      end do

      end if
#else
      if (use_polarshortreal) then
         udalt_p0 => udshortalt
         upalt_p0 => upshortalt
      else
         udalt_p0 => udalt
         upalt_p0 => upalt
      end if

!$acc parallel loop gang vector async(def_queue)
!$acc&         present(poleglob,bpred,udalt_p0,upalt_p0,mu)
      do i = 1, npoleloc
         iipole = poleglob(i)
!$acc loop seq
         do j = 1, 3
            udsum = 0.0_ti_p
            upsum = 0.0_ti_p
            do k = 1, nualt - 1
              udsum = udsum + bpred(k)*udalt_p0(k,j,iipole)
              upsum = upsum + bpred(k)*upalt_p0(k,j,iipole)
            end do
            mu(j,1,i) = udsum
            mu(j,2,i) = upsum
         end do
      end do
#endif

 50   continue
      call timer_exit( timer_ulspred )
      end subroutine

      subroutine pred_SetAlt
      use atmlst ,only: poleglob
      use atoms  ,only: n_pe
      use domdec ,only: nproc,rank,ndir
      use inform ,only: deb_Path
      use mpole
      use nvshmem
      use polar
      use potent ,only: use_pmecore,use_polarshortreal
      use timestat
      use tinheader ,only: ti_p
      use uprior
      use utilgpu   ,only: def_queue
      implicit none
      integer i,j,k,iipole,ierr
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif

      ! TODO use shortnualt instead
      if (deb_Path) print '(3x,A)',"pred_SetAlt"
      call timer_enter( timer_ulspred )

      nualt = min(nualt+1,maxualt)

#ifdef USE_NVSHMEM_CUDA
      if (use_polarshortreal) then
         udalt_p1 => d_udshortalt
         upalt_p1 => d_upshortalt
      else
         udalt_p1 => d_udalt
         upalt_p1 => d_upalt
      end if

      if (nproc.eq.npes.and..not.use_pmecore) then

!$acc parallel loop collapse(2) async(def_queue)
!$acc&     present(poleglob,uind,uinp) deviceptr(udalt_p1,upalt_p1)
      do i = 1, npolebloc
         do j = 1, 3
            iipole = poleglob(i)
            ipe    = (iipole-1)/n_pe
            ind    = mod((iipole-1),n_pe) +1
            if (mype.eq.ipe) then
!$acc loop seq
            do k = nualt, 2, -1
              udalt_p1(ipe)%pel(k,j,ind) = udalt_p1(ipe)%pel(k-1,j,ind)
              upalt_p1(ipe)%pel(k,j,ind) = upalt_p1(ipe)%pel(k-1,j,ind)
            end do
            udalt_p1(ipe)%pel(1,j,ind) = uind(j,iipole)
            upalt_p1(ipe)%pel(1,j,ind) = uinp(j,iipole)
            end if
         end do
      end do

      else

!$acc parallel loop collapse(2) async(def_queue)
!$acc&         present(poleglob) 
!$acc&         deviceptr(udalt_p1,upalt_p1)
      do i = 1, n_pe
         do j = 1, 3
!$acc loop seq
            do k = nualt, 2, -1
              udalt_p1(mype)%pel(k,j,i) = udalt_p1(mype)%pel(k-1,j,i)
              upalt_p1(mype)%pel(k,j,i) = upalt_p1(mype)%pel(k-1,j,i)
            end do
         end do
      end do
      if (npes.ne.1) call mpi_barrier(COMM_NVSHMEM,ierr)

      !TODO Find a way to reduce RMA in this kernel
!$acc parallel loop collapse(2) async(def_queue)
!$acc&         present(poleglob,uind,uinp)
!$acc&         deviceptr(udalt_p1,upalt_p1)
      do i = 1, npolebloc
         do j = 1, 3
            iipole = poleglob(i)
            ipe    = (iipole-1)/n_pe
            ind    = mod((iipole-1),n_pe) +1
            udalt_p1(ipe)%pel(1,j,ind) = uind(j,iipole)
            upalt_p1(ipe)%pel(1,j,ind) = uinp(j,iipole)
         end do
      end do

      end if
#else
      ! PME-CORE exeption
      if (.not.use_pmecore.or.use_pmecore.and.rank.lt.ndir) then

      if (use_polarshortreal) then
         udalt_p0 => udshortalt
         upalt_p0 => upshortalt
      else
         udalt_p0 => udalt
         upalt_p0 => upalt
      end if

!$acc parallel loop collapse(2) async(def_queue)
!$acc&     present(poleglob,uind,uinp,upalt_p0,udalt_p0)
      do i = 1, npolebloc
         do j = 1, 3
            iipole = poleglob(i)
!$acc loop seq
            do k = nualt, 2, -1
               udalt_p0(k,j,iipole) = udalt_p0(k-1,j,iipole)
               upalt_p0(k,j,iipole) = upalt_p0(k-1,j,iipole)
            end do
            udalt_p0(1,j,iipole) = uind(j,iipole)
            upalt_p0(1,j,iipole) = uinp(j,iipole)
          end do
      end do

      end if
#endif
      call timer_exit( timer_ulspred )
      end subroutine
c
c
c
      subroutine inducepcg_pme2gpu(matvec,nrhs,precnd,ef,mu,murec)
      use atmlst
      use domdec
      use ewald
      use iounit
      use inform    ,only: deb_Path,abort
      use interfaces,only: tmatxb_pmegpu
      use math
      use mpole
      use polar
      use polpot
      use sizes
      use timestat
      use units
      use utils
      use utilcomm ,buffermpi1=>buffermpi3d,buffermpi2=>buffermpi3d1
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
      real(t_p),intent(in)   :: ef   (:,:,:)
      real(t_p),intent(inout):: mu   (:,:,:)
      real(t_p),intent(inout):: murec(:,:,:)
      procedure(tmatxb_pmegpu) :: matvec

      integer i, it, j, k
      real(t_p) ggold(2), alphacg(2), gg(2)
      real(r_p) gnorm(2), ggnew(2), ene(2)
      real(t_p),save:: ggold1=0.0,ggold2=0.0
      real(r_p),save:: ggnew1=0.0,ggnew2=0.0
      real(t_p),save:: ggdev1=0.0,ggdev2=0.0
      real(r_p),save:: gnorm1=0.0,gnorm2=0.0
      real(t_p),save:: gg1=0.0,gg2=0.0
      real(r_p),save:: ene1=0.0,ene2=0.0
      real(t_p),save:: alphacg1=0.0,alphacg2=0.0
      real(t_p)  zero, one, term
      real(r_p) zerom, pt5
      real(r_p) resnrm
      real(t_p),save,allocatable :: res(:,:,:), h(:,:,:),
     &                pp(:,:,:), zr(:,:,:), diag(:)
      parameter(zero=0.0,zerom=0,pt5=0.5,one=1.0)
c
c     MPI
c
      real(t_p),save,allocatable :: dipfield(:,:,:),
     &                       dipfieldbis(:,:,:)
      integer iglob, iipole, ierr, tag, proc, sizeloc, sizebloc
      integer ::place, place_tmp=0, place1=0
      integer req1, req2, req3, req4, req5
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer,allocatable :: reqendrec(:),reqendsend(:)
      integer,allocatable :: req2endrec(:),req2endsend(:)
      logical,save::first_in=.true.
      logical tinker_isnan_m
c
 1000 format(' cgiter converged after ',I3,' iterations.',/,
     &       ' final energy        = ',2D14.7,/,
     &       ' final residual norm = ',2D14.7)
 1010 format(' energy and residual norm at iteration ',I3,':',4D12.2)
 1020 format(' Conjugate gradient solver: induced dipoles',/,
     &  ' ipole       mux         muy         muz')
 1021 format(' Conjugate gradient solver: induced p-dipoles',/,
     &  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
 1040 format(' Using a diagonal preconditioner.')
 1050 format(' Induce PCG solver not converged after ',I6,' iterations'
     &    ,/,' Residual norm ', d10.5)

      if (deb_Path) 
     &   write(*,'(3x,a)') 'inducepcg_pme2gpu'
c

      call timer_enter ( timer_other )
      if (nproc.gt.1) then
      call prmem_request(buffermpi1,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(buffermpimu1,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(buffermpi2,3,nrhs,max(npolerecloc,1),
     &     async=.true.)
      call prmem_request(buffermpimu2,3,nrhs,max(npolerecloc,1),
     &     async=.true.)
      end if

      call prmem_request(dipfield,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(dipfieldbis,3,nrhs,max(npolerecloc,1),
     &     async=.true.)
      call prmem_request(zr,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(res,3,nrhs,max(npoleloc,1),
     &     async=.true.)
      call prmem_request(h ,3,nrhs,max(npolebloc,1),
     &     async=.true.)
      call prmem_request(pp ,3,nrhs,max(npolebloc,1),
     &     async=.true.)
      call prmem_request(diag,npoleloc,async=.true.)
c
      allocate (reqrecdirrec(nproc))
      allocate (reqrecdirsend(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
      allocate (reqendsend(nproc))
      allocate (reqendrec(nproc))
      allocate (req2endsend(nproc))
      allocate (req2endrec(nproc))
c     allocate (dipfield(3,nrhs,max(1,npoleloc)))
c     allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
c
c     allocate some memory and setup the preconditioner:
c
c     allocate (zr (3,nrhs,max(1,npoleloc)))
c     allocate (res(3,nrhs,max(1,npoleloc)))
c     allocate (h  (3,nrhs,max(1,npolebloc)))
c     allocate (pp (3,nrhs,max(1,npolebloc)))
c     allocate (diag(npoleloc))
!!$acc enter data create(res,pp,dipfield,dipfieldbis) async
!!$acc data create(zr,h,diag) async

      sizeloc  = 3*nrhs*npoleloc
      sizebloc = 3*nrhs*npolebloc
      !TODO This wait might be useless
!$acc wait
      if (first_in) then
!$acc enter data create(gg1,gg2) async
         first_in=.false.
      end if


      if (precnd) then
!$acc parallel loop present(poleglob,polarity) async
         do i = 1, npoleloc
            iipole = poleglob(i)
            diag(i) = polarity(iipole)
         end do
         if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
!$acc parallel loop async
         do i = 1, npoleloc
            diag(i) = one
         end do
      end if
c
c     initialize to zero
c
      if (nproc.eq.1) then
      call set_to_zero3(res,zr,dipfield,3*nrhs*npoleloc,rec_queue)
      call set_to_zero1(dipfieldbis,3*nrhs*npolerecloc,rec_queue)
      call set_to_zero1(pp,3*nrhs*npolebloc  ,rec_queue)
      else
      call set_to_zero5(buffermpi1,buffermpimu1,res,zr,dipfield,
     &                     3*nrhs*npoleloc   ,rec_queue)
      call set_to_zero3(buffermpimu2,buffermpi2,dipfieldbis,
     &                     3*nrhs*npolerecloc,rec_queue)
      call set_to_zero1(pp,3*nrhs*npolebloc  ,rec_queue)
      end if
      call timer_exit( timer_other,quiet_timers )
c     now, compute the initial direction
c
      ggold  = 0.0_ti_p
      ggold1 = 0.0_ti_p
      ggold2 = 0.0_ti_p

#ifdef _OPENACC
      ! Start async overlapping
      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)
#endif
c
      call timer_enter( timer_realdip )
      call matvec(nrhs,.true.,mu,h)
      call timer_exit ( timer_realdip,quiet_timers )

      call timer_enter( timer_recdip )
      call tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
      call timer_exit ( timer_recdip,quiet_timers )

      call commfieldgpu(nrhs,h)
      call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $     buffermpi2,reqrecdirrec,reqrecdirsend)
      call commdirdirgpu(nrhs,0,pp,reqrec,reqsend)
c
      call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $     buffermpi2,reqrecdirrec,reqrecdirsend)
      call commrecdirdipgpu(nrhs,0,murec,pp,buffermpimu1,
     $     buffermpimu2,req2rec,req2send)
      call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $     buffermpi2,reqrecdirrec,reqrecdirsend)
      term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi

      call timer_enter( timer_other )
#ifdef _OPENACC
      ! Start async overlapping
      if (dir_queue.ne.rec_queue) then
         call stream_wait_async(dir_stream,rec_stream,dir_event)
      end if
#endif
!$acc parallel loop collapse(3) present(ef) async
      do k=1,npoleloc
         do j=1,nrhs
            do i=1,3
               h(i,j,k)   = h(i,j,k)  + dipfield(i,j,k) - term*mu(i,j,k)
               res(i,j,k) = ef(i,j,k) - h(i,j,k)
               zr(i,j,k)  = res(i,j,k)*diag(k)
               pp(i,j,k)  = zr(i,j,k)
               if (btest(j,0)) then  ! test de 1
                  ggold1  = ggold1 + res(i,j,k)*zr(i,j,k)
               else
                  ggold2  = ggold2 + res(i,j,k)*zr(i,j,k)
               end if
            end do
         end do
      end do
!$acc wait

      ggold(1) = ggold1
      ggold(2) = ggold2

      call MPI_IALLREDUCE(MPI_IN_PLACE,ggold(1),nrhs,MPI_TPREC,MPI_SUM,
     &     COMM_TINKER,req1,ierr)
      call timer_exit ( timer_other )

c
c     MPI : begin sending
c
      call commdirdirgpu(nrhs,1,pp,reqrec,reqsend)
      call commdirdirgpu(nrhs,2,pp,reqrec,reqsend)
      call commrecdirdipgpu(nrhs,1,murec,pp,buffermpimu1,
     &     buffermpimu2,req2rec,req2send)
      call commrecdirdipgpu(nrhs,2,murec,pp,buffermpimu1,
     &     buffermpimu2,req2rec,req2send)
      call MPI_WAIT(req1,status,ierr)

      ggold1 = ggold(1)
      ggold2 = ggold(2)

#ifdef _OPENACC
      !start Async overlapping with matvec
      if (dir_queue.ne.rec_queue)
     &   call stream_wait_async(rec_stream,dir_stream,rec_event)
#endif

!!$acc update host(pp,murec,zr,res,h,dipfield,dipfieldbis,diag)
c
c     now, start the main loop:
c
      do it = 1,politer
c
         call timer_enter( timer_realdip )
         call matvec(nrhs,.true.,pp,h)
         call timer_exit ( timer_realdip,quiet_timers )

         call timer_enter( timer_recdip )
         call tmatxbrecipgpu(pp,murec,nrhs,dipfield,dipfieldbis)
         call timer_exit ( timer_recdip,quiet_timers )

!$acc serial async present(gg1,gg2)
         gg1=zero; gg2=zero
!$acc end serial
c
         call commfieldgpu(nrhs,h)
c
c     Begin the reception of the reciprocal fields
c
         call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,reqrecdirrec,reqrecdirsend)
c
c     MPI : begin reception
c
         call commdirdirgpu(nrhs,0,pp,reqrec,reqsend)
         call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,reqrecdirrec,reqrecdirsend)
c
c     Begin reception of mu for PME
c
         call commrecdirdipgpu(nrhs,0,murec,pp,buffermpimu1,
     &        buffermpimu2,req2rec,req2send)
c
c       Wait for the reciprocal fields
c
         call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,reqrecdirrec,reqrecdirsend)

         term = (4.0_ti_p/3.0_ti_p)*aewald**3 / sqrtpi

         call timer_enter( timer_other )
         if (dir_queue.ne.rec_queue) then
!$acc wait(dir_queue)
         end if
!$acc parallel loop collapse(3) present(gg1,gg2) async
!$acc&         default(present)
         do k=1,npoleloc
            do j=1,nrhs
               do i=1,3
                  h(i,j,k) = h(i,j,k) + dipfield(i,j,k) - term*pp(i,j,k)
                  if (btest(j,0)) then
                     gg1 = gg1 + pp(i,j,k)*h(i,j,k)
                  else
                     gg2 = gg2 + pp(i,j,k)*h(i,j,k)
                  end if
               end do
            end do
         end do

         if (nproc.gt.1) then
!$acc wait
!$acc host_data use_device(gg1,gg2)
            call MPI_ALLREDUCE(MPI_IN_PLACE,gg1,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,gg2,1,MPI_TPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
!$acc end host_data
            !call MPI_WAIT(req2,status,ierr)
            !call MPI_WAIT(req5,status,ierr)
         end if

         ! This check might be unecessary
         ! the only for residu this to reach zero is for the field itself to be equal to zero
         ! TODO Find a proper place for that check
         if (deb_Path) then
!$acc update host(gg1,gg2) async
!$acc wait
            if (gg1.eq.zerom) then
               print*, 'Residu equals zero in polarisation solver'
               goto 50
            end if
            if (gg2.eq.zerom) goto 50
         end if

         ! reduce division operation in next kernel
         !alphacg1  = ggold1/gg1
         !alphacg2  = ggold2/gg2

         ggnew1 = zerom; ggnew2 = zerom
         ene1 = zerom; ene2 = zerom

!$acc parallel loop collapse(3) async present(gg1,gg2,mu,ef)
         do k=1,npoleloc
            do j=1,nrhs
               do i=1,3
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
               end do
            end do
         end do
c Implicit wait seems to be removed with PGI-20
!$acc wait
         ene2 = -pt5*ene2; ene1 = -pt5*ene1
         ggnew(1)=ggnew1; ggnew(2)=ggnew2
           ene(1)=  ene1;   ene(2)=  ene2

         call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew(1),nrhs,MPI_RPREC,
     &        MPI_SUM,COMM_TINKER,req3,ierr)
         call MPI_IALLREDUCE(MPI_IN_PLACE,ene(1),nrhs,MPI_RPREC,MPI_SUM,
     &        COMM_TINKER,req4,ierr)
         call MPI_WAIT(req3,status,ierr)
         call MPI_WAIT(req4,status,ierr)
         resnrm = 0.0_re_p

c        do k = 1, nrhs
c           pp(:,k,1:npoleloc) = zr(:,k,1:npoleloc)+ggnew(k)/ggold(k)*
c    &                         pp(:,k,1:npoleloc)
c        end do

         ggdev1 = ggnew(1)/ggold(1); ggdev2 = ggnew(2)/ggold(2)
!$acc parallel loop collapse(3) async
         do k=1,npoleloc
            do j=1,nrhs
               do i=1,3
                  if (btest(j,0)) then
                     pp(i,j,k) = zr(i,j,k)+ ggdev1*pp(i,j,k)
                  else
                     pp(i,j,k) = zr(i,j,k)+ ggdev2*pp(i,j,k)
                  end if
               end do
            end do
         end do

         call timer_exit( timer_other )
c
         call commdirdirgpu(nrhs,1,pp,reqrec,reqsend)
c
         do k = 1, nrhs
            gnorm(k) = sqrt(ggnew(k)/real(3*npolar,r_p))
            resnrm   = max(resnrm,gnorm(k))
         end do

         if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     &      it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)
c
         ggold  = ggnew
         ggold1 = ggnew(1)
         ggold2 = ggnew(2)
         call commdirdirgpu(nrhs,2,pp,reqrec,reqsend)
         call commrecdirdipgpu(nrhs,1,murec,pp,buffermpimu1,
     &        buffermpimu2,req2rec,req2send)
         call commrecdirdipgpu(nrhs,2,murec,pp,buffermpimu1,
     &        buffermpimu2,req2rec,req2send)
c
         if (resnrm.lt.poleps) then
            if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     &         (ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
            goto 10
         end if
         if ((it.eq.politer.and.resnrm.gt.poleps)
     &      .or.tinker_isnan_m(resnrm)) then
            ! Not converged Exit
            if (rank.eq.0) write(0,1050) it,resnrm
            abort = .true.
         end if
c
      end do
 10   continue

c
c     Begin reception of mu for PME
c
      call commdirdirgpu(nrhs,0,mu,reqendrec,reqendsend)
      call commrecdirdipgpu(nrhs,0,murec,mu,buffermpimu1,
     &     buffermpimu2,req2endrec,req2endsend)
c
c     MPI : begin direct space communication
c
      call commdirdirgpu(nrhs,1,mu,reqendrec,reqendsend)
      call commdirdirgpu(nrhs,2,mu,reqendrec,reqendsend)
c
c     Exchange data between direct and reciproqual space
c
      call commrecdirdipgpu(nrhs,1,murec,mu,buffermpimu1,
     &     buffermpimu2,req2endrec,req2endsend)
      call commrecdirdipgpu(nrhs,2,murec,mu,buffermpimu1,
     &     buffermpimu2,req2endrec,req2endsend)
c
c     finalize and debug printing:
c
      if (polprt.ge.3) then
!$acc wait
!$acc update host(mu)
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
            iglob = ipole(iipole)
            write(iout,1030) iglob, (mu(j,2,i), j = 1, 3)
         end do
      end if
 50   continue

!!$acc end data
!!$acc exit data delete(res,pp,dipfield,dipfieldbis) async
      deallocate (reqrecdirsend)
      deallocate (reqrecdirrec)
      deallocate (reqsend)
      deallocate (reqrec)
      deallocate (req2send)
      deallocate (req2rec)
      deallocate (reqendsend)
      deallocate (reqendrec)
      deallocate (req2endsend)
      deallocate (req2endrec)
c     deallocate (buffermpi1)
c     deallocate (buffermpimu1)
c     deallocate (buffermpimu2)
c     deallocate (buffermpi2)
c     deallocate (dipfield)
c     deallocate (dipfieldbis)
c     deallocate (res)
c     deallocate (h)
c     deallocate (pp)
c     deallocate (zr)
c     deallocate (diag)
      end
c
      subroutine inducejac_pme2gpu(matvec,nrhs,dodiis,ef,mu,murec)
      use atoms
      use atmlst
      use domdec
      use ewald
      use inform ,only: deb_Path
      use iounit
      use math
      use mpole
      use polar
      use polpot
      use potent
      use timestat
      use mpi
      use utilgpu
      implicit none
c
c     solves the polarizable dipoles linear equations by jacobi
c     iterations coupled to the DIIS extrapolation.
c
      integer nrhs, info, proc
      real(t_p)  ef(3,nrhs,*), mu(3,nrhs,*), murec(3,nrhs,*)
      real(t_p), allocatable :: munew(:,:,:), h(:,:,:)
      real(t_p), allocatable :: xdiis(:,:), ediis(:,:),
     $       bmat(:,:), bloc(:,:), cex(:)
      integer i, j, k, it, ind, ndismx, nmat, lenb, tag
      logical dodiis
      parameter (ndismx=25)
      integer ipiv(ndismx+1)
      real(t_p)  zero, one, rnorm(2), rr, xx(1)
      real(t_p) term
      real(t_p) time0,time1,time2,time3
      real(t_p), allocatable :: dipfield(:,:,:),dipfieldbis(:,:,:)
c
c     MPI
c
      real(t_p), allocatable::buffermpi1(:,:,:),buffermpi2(:,:,:)
      real(t_p), allocatable::buffermpimu1(:,:,:),buffermpimu2(:,:,:)
      integer iglob, iipole, ierr
      integer reqnorm,reqdiis(2*ndismx+1)
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      external matvec
      parameter(zero=0.0_ti_p,one=1.0_ti_p,xx=0.0_ti_p)
c
 1000 format(' itsolv converged after ',I3,' iterations.',/,
     $       ' final residual norm = ',3D14.7)
 1010 format(' residual norm at iteration ',I3,':',3D12.2)
 1020 format(' Jacobi/DIIS solver: induced dipoles',/,
     $  ' ipole       mux         muy         muz')
 1021 format(' Jacobi/DIIS solver: induced p-dipoles',/,
     $  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
      if (deb_Path) 
     &   write(*,'(3x,a)') 'inducejac_pme2gpu'
      
      call openacc_abort("induce jac unported use polar-alg 1")
c
c
      allocate (buffermpi1(3,nrhs,max(npoleloc,1)))
      buffermpi1 = zero
      allocate (buffermpimu1(3,nrhs,max(npoleloc,1)))
      buffermpimu1 = zero
      allocate (buffermpi2(3,nrhs,max(npolerecloc,1)))
      buffermpi2 = zero
      allocate (buffermpimu2(3,nrhs,max(npolerecloc,1)))
      buffermpimu2 = zero
c
      allocate (munew(3,nrhs,max(1,npolebloc)))
      munew = zero
      allocate (dipfield(3,nrhs,max(1,npoleloc)))
      dipfield = zero
      allocate (dipfieldbis(3,nrhs,max(1,npolerecloc)))
      dipfieldbis = zero
      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqrec(nproc))
      allocate (reqsend(nproc))
      allocate (req2rec(nproc))
      allocate (req2send(nproc))
c
      allocate (h(3,nrhs,max(1,npolebloc)))
      h = zero
      if (dodiis) then
        nmat = 1
        lenb = ndismx + 1
        allocate (xdiis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (ediis(3*nrhs*max(npoleloc,1),ndismx))
        allocate (bmat(lenb,lenb))
        allocate (bloc(lenb,lenb))
        allocate (cex(lenb))
        bmat = zero
      end if
c
c     main loop:
c
      do it = 1, politer
        h = zero
        rnorm = zero
c
        call timer_enter( timer_recdip )
        call tmatxbrecip(mu,murec,nrhs,dipfield,dipfieldbis)
        call timer_exit ( timer_recdip,quiet_timers )

        call timer_enter( timer_realdip )
        call matvec(nrhs,.false.,mu,h)
        call timer_exit ( timer_realdip,quiet_timers )

        call commfield(nrhs,h)
c
        call commrecdirsolv(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c
        call commdirdirgpu(nrhs,0,mu,reqrec,reqsend)
        call commrecdirsolv(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c
        call commrecdirdip(nrhs,0,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
c
c     jacobi step:
c
        call commrecdirsolv(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $    buffermpi2,reqrecdirrec,reqrecdirsend)
c        
        h(1:3,1:nrhs,1:npoleloc) = h(1:3,1:nrhs,1:npoleloc) + 
     $   dipfield(1:3,1:nrhs,1:npoleloc)

        term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
        h(1:3,1:nrhs,1:npoleloc) = h(1:3,1:nrhs,1:npoleloc) 
     $   - term*mu(1:3,1:nrhs,1:npoleloc)
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
          rnorm(k) = sum((munew(1:3,k,1:npoleloc) -
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
          cex = zero
          cex(1) = one
#ifdef SINGLE
          call sgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#elif defined(MIXED)
          call sgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#else
          call dgesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
#endif
          munew = zero
          call extrap(3*nrhs*npoleloc,nmat-1,xdiis,cex,munew)
        end if
c
        call commdirdirgpu(nrhs,1,munew,reqrec,reqsend)
        call commrecdirdip(nrhs,1,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
        call commdirdirgpu(nrhs,2,mu,reqrec,reqsend)
        call commrecdirdip(nrhs,2,murec,munew,buffermpimu1,
     $  buffermpimu2,req2rec,req2send)
        mu(:,:,1:npoleloc) = munew(:,:,1:npoleloc)

c
c     compute the norm of the increment.
c
        call MPI_WAIT(reqnorm,status,ierr)
        rr = zero
        do k = 1, nrhs
         rnorm(k) = sqrt(rnorm(k)/real(3*npolar,t_p))
          rr = max(rnorm(k),rr)
        end do
        if (rr.lt.poleps) then
          if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1000) it, (rnorm(k), k = 1, nrhs)
          goto 10
        end if
        if (polprt.ge.1.and.rank.eq.0)
     $       write(6,1010) it, (rnorm(k), k = 1, nrhs)
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
      deallocate (reqrecdirsend)
      deallocate (reqrecdirrec)
      deallocate (buffermpi1)
      deallocate (buffermpimu2)
      deallocate (buffermpi2)
      deallocate (munew)
      deallocate (dipfield)
      deallocate (dipfieldbis)
      deallocate (h)
      if (dodiis) then
        deallocate (xdiis)
        deallocate (ediis)
        deallocate (bmat)
        deallocate (bloc)
        deallocate (cex)
      end if
      return
      end
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ulspred  --  induced dipole prediction coeffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ulspred" uses standard extrapolation or a least squares fit
c     to set coefficients of an induced dipole predictor polynomial
c
c     literature references:
c
c     J. Kolafa, "Time-Reversible Always Stable Predictor-Corrector
c     Method for Molecular Dynamics of Polarizable Molecules", Journal
c     of Computational Chemistry, 25, 335-342 (2004)
c
c     W. Wang and R. D. Skeel, "Fast Evaluation of Polarizable Forces",
c     Journal of Chemical Physics, 123, 164107 (2005)
c
c
      subroutine ulspredgpu
      use domdec
      use inform ,only: deb_Path
      use mpole
      use sizes
      use uprior
      use utils
      use utilgpu
      implicit none
      integer i,j,k,m
      integer iwh
      integer:: nblock=1
      integer:: array(maxualt)
      real(t_p) coeff,udk,upk
      real(t_p) amax,apmax
      real(t_p) cmk,cpmk
      real(t_p) zero,atm
      real(t_p) b(maxualt)
      real(t_p) bp(maxualt)
      real(t_p) a(maxualt*(maxualt+1)/2)
      real(t_p) ap(maxualt*(maxualt+1)/2)
      real(t_p) c(maxualt,maxualt)
      real(t_p) cp(maxualt,maxualt)
      real(t_p),allocatable:: rc(:,:,:),rcp(:,:,:)  
      parameter(zero=0.0_ti_p)


      if (deb_Path) write(*,'(4x,a)') 'ulspredgpu'
c
!$acc data create(b,bp,a,ap,c,cp,rc,rcp,array)
!$acc&     present(gear,aspc,bpred,bpredp,bpreds,bpredps)
!$acc&     async(def_queue)
c
c     set the Gear predictor binomial coefficients
c
      if (polpred .eq. 'GEAR') then
!$acc parallel loop vector_length(32) async(def_queue)
         do i = 1, nualt
            coeff      = gear(i)
            bpred(i)   = coeff
            bpredp(i)  = coeff
            bpreds(i)  = coeff
            bpredps(i) = coeff
         end do
c
c     set always stable predictor-corrector (ASPC) coefficients
c
      else if (polpred .eq. 'ASPC') then
!$acc parallel loop vector_length(32) async(def_queue)
         do i = 1, nualt
            coeff      = aspc(i)
            bpred(i)   = coeff
            bpredp(i)  = coeff
            bpreds(i)  = coeff
            bpredps(i) = coeff
         end do
c
c     derive normal equations corresponding to least squares fit
c
      else
         if (npole.gt.cores_SMP) nblock=cores_SMP
         allocate( rc(nblock,maxualt,maxualt))
         allocate(rcp(nblock,maxualt,maxualt))
!$acc data create(rc,rcp) async(def_queue)
c
!$acc parallel loop collapse(2)
!$acc&         vector_length(32) 
!$acc&         async(def_queue)
         do k = 1, nualt
            do m = 1, nualt
                 if (m.lt.k) cycle
                  b(k) = zero
                 bp(k) = zero
                c(m,k) = zero
               cp(m,k) = zero
            end do
         end do
         call set_to_zero2(rc,rcp,nblock*maxualt*maxualt,def_queue)
!$acc parallel loop collapse(2) async(def_queue)
         do i = 1, npole
            do j = 1, 3
               iwh = 1+mod(3*(i-1)+j-1,nblock)
!$acc loop seq
               do k = 1, nualt
                  udk = udalt(k,j,i)
                  upk = upalt(k,j,i)
!$acc loop seq
                  do m = k, nualt
                     atm =  udk*udalt(m,j,i)
!$acc atomic update
                      rc(iwh,m,k) = rc(iwh,m,k)  + atm
                     atm =  upk*upalt(m,j,i)
!$acc atomic update
                     rcp(iwh,m,k) = rcp(iwh,m,k) + atm
                  end do
               end do
            end do
         end do

!$acc parallel loop collapse(2) async(def_queue)
         do k= 1, nualt
            do m= 1, nualt
               if (m.lt.k) cycle
               cmk = zero
               cpmk= zero
!$acc loop vector
               do i= 1,nblock
                  cmk  = cmk  +  rc(i,m,k)
                  cpmk = cpmk + rcp(i,m,k)
               end do
                c(m,k) = cmk
               cp(m,k) = cpmk
            end do
         end do

c!$acc parallel loop collapse(3) async(def_queue)
c         do i = 1, npole
c            do j = 1, 3
c               do k = 1, nualt
c                  udk = udalt(k,j,i)
c                  upk = upalt(k,j,i)
c!$acc loop seq
c                  do m = k, nualt
c                     atm =  udk*udalt(m,j,i)
c!$acc atomic update
c                      c(m,k) =  c(m,k) + atm
c                     atm =  upk*upalt(m,j,i)
c!$acc atomic update
c                     cp(m,k) = cp(m,k) + atm
c                  end do
c               end do
c            end do
c         end do
         !i = 0
!$acc kernels async(def_queue)
         array = (/0,0,6,11,15,18,20/)
!$acc loop collapse(2)
         do k = 2, nualt
            do m = 2, nualt
               if (m.lt.k) cycle
               b(k-1)  = c(k,1)
               bp(k-1) = cp(k,1)
               i = array(k)+m-k+1
               !i = i + 1
               a(i)  = c(m,k)
               ap(i) = cp(m,k)
            end do
         end do
c
c     check for nonzero coefficients and solve normal equations
c
         k = nualt - 1
         amax  = zero
         apmax = zero
!$acc loop
         do i = 1, k*(k+1)/2
            amax  = max(amax,a(i))
            apmax = max(apmax,ap(i))
         end do
!$acc end kernels
         if ( amax .ne. zero)  call cholesky (k,a,b)
         if (apmax .ne. zero)  call cholesky (k,ap,bp)
c
c     transfer the final solution to the coefficient vector
c
!$acc kernels async(def_queue)
         do i = 1, nualt-1
            bpred  (i) = b (i)
            bpredp (i) = bp(i)
            bpreds (i) = b (i)
            bpredps(i) = bp(i)
         end do
         bpred  (nualt) = zero
         bpredp (nualt) = zero
         bpreds (nualt) = zero
         bpredps(nualt) = zero
!$acc end kernels
!$acc exit data delete(rc,rcp) async(def_queue)

!$acc end data

         deallocate(rc)
         deallocate(rcp)
      end if

!$acc end data
      return
      end
c
!===================================================
!     sub diagvecvec
!===================================================
! Performs product of vector a with polarisabilities
!
      subroutine diagvecgpu(nrhs, A, B)
      use atmlst
      use mpole
      use polar
      implicit none

      integer, intent(in) :: nrhs
      real(t_p), dimension(3,nrhs,npolebloc) :: A,B
      integer :: i,iipole, irhs, j

!$acc parallel loop collapse(3) present(A,B)
      do i = 1, npolebloc
         do irhs = 1, nrhs
            do j = 1,3
               iipole = poleglob(i)
               B(j,irhs,i) = A(j,irhs,i)*polarity(iipole)
            end do
         end do
      end do

      return
      end
