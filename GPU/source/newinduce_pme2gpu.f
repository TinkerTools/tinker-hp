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
#include "tinker_macro.h"
      subroutine newinduce_pme2gpu
      use atoms     ,only: n,n_pe
      use atomsmirror
      use atmlst
      use domdec
      use deriv
      use ewald
      use iounit
      use inform    ,only: deb_Path,minmaxone
      use interfaces,only: inducepcg_pme2gpu,tmatxb_p
     &                   , inducestepcg_pme2gpu
     &                   , efld0_directgpu2, efld0_directgpu_p
      use math
      use mdstate   ,only: track_mds,ms_back_p
      use moldyn    ,only: step_c,stepint
      use mpole
      use mutant    ,only: nmut,elambda
      use nvshmem
      use pme
      use polar
      use polar_temp,only: ef,mu,murec,cphi
      use polpot
      use potent
      use units
      use uprior
      use utils
      use utilcomm  ,buffermpi1=>buffermpi2d,buffermpi2=>buffermpi2d1
      use utilgpu
      use timestat
      use tinheader ,only: i_init
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
c
      if (.not.use_polar) return
c
 1000 format(' illegal polalg(',I3,') in newinduce.')
 1010 format(' time for the ',a,F14.5)
 1020 format(' total elapsed time in newinduce: ',F14.5)
 
      call timer_enter( timer_other )
      if (deb_Path) write(*,'(2x,a)') 'newinduce_pme2gpu'
c
c     allocate some memory and clear the arrays:
c
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
      j = max(npolerecloc,1)
      call prmem_request(cphirec,10,j,async=.true.)
      call prmem_request(fphirec,20,j,async=.true.)

      allocate (reqrecdirsend(nproc))
      allocate (reqrecdirrec(nproc))
      allocate (reqsend(nproc))
      allocate (reqrec(nproc))
      allocate (req2send(nproc))
      allocate (req2rec(nproc))

      ! set arrays to zero
      call set_to_zero2(mu,ef,3*nrhs*npolebloc,rec_queue)
      if (nproc.eq.1) then
         call set_to_zero1(murec,3*nrhs*npolerecloc,rec_queue)
         call set_to_zero1(cphi   ,10*npoleloc   ,rec_queue)
      else
         call set_to_zero1(buffermpimu1,3*nrhs*npoleloc ,rec_queue)
         call set_to_zero2(murec,buffermpimu2,3*nrhs*npolerecloc,
     &                    rec_queue)
         call set_to_zero2(cphi   ,buffermpi1,10*npoleloc   ,rec_queue)
         call set_to_zero1(buffermpi2,10*npolerecloc,rec_queue)
      end if
      call timer_exit (timer_other,quiet_timers)
c
c     compute the electric fields:
c
      wtime0 = mpi_wtime()
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
      call efld0_recipgpu(cphi,ef)
      call timer_exit( timer_recdip,quiet_timers )
c
      if (aewald.ge.1d-6) then
      call commrecdirfieldsgpu(0,cphirec,cphi,buffermpi1,buffermpi2,
     &     reqrecdirrec,reqrecdirsend)
      call commrecdirfieldsgpu(1,cphirec,cphi,buffermpi1,buffermpi2,
     &     reqrecdirrec,reqrecdirsend)
      call commrecdirfieldsgpu(2,cphirec,cphi,buffermpi1,buffermpi2,
     &     reqrecdirrec,reqrecdirsend)
      end if
c
      call commfieldgpu(nrhs,ef)
#ifdef _OPENACC
      ! Prevent overlap between commfield and guess construct
      if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif
      call commdirdirgpu(nrhs,0,mu,reqrec,reqsend)
c
c     add direct and reciprocal fields
c
      term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
      def_queue = rec_queue

!$acc parallel loop collapse(3) async(def_queue) default(present)
      do i = 1,npoleloc; do j = 1,nrhs; do k = 1,3
         iipole = poleglob(i)
         comp = term*rpole(k+1,iipole) - cphi(k+1,i)
         ef(k,j,i) = ef(k,j,i) + comp
      end do; end do; end do
      wtime1 = mpi_wtime()
c
c     Reset predictor if MD State Tracking is enabled 
c
      if (track_mds.and.mod(step_c,ms_back_p).eq.1.and.stepint
     &   .eq.i_init) call pred_reset
c
c     guess the dipoles:
c
c     predicted values for always stable predictor-corrector method
c
      if (use_pred .and. nualt.eq.maxualt) then
         call pred_InitField(mu)
      else if (polgsf.eq.0) then
!$acc parallel loop collapse(3) async(def_queue) default(present)
         do i = 1,npoleloc; do k = 1,nrhs; do j = 1,3
            iipole = poleglob(i)
            mu(j,k,i) = polarity(iipole)*ef(j,k,i)
         end do; end do; end do
      else
!$acc parallel loop collapse(3) async(def_queue) default(present)
         do i = 1,npoleloc; do j = 1,nrhs; do k = 1,3
            iipole = poleglob(i)
            if (btest(j,0)) then !j==1
               mu(k,j,i) = uind(k,iipole)
            else  !j/=1
               mu(k,j,i) = uinp(k,iipole)
            end if
         end do; end do; end do
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
      if (polalg.eq.pcg_SId) then
         call inducepcg_pme2gpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
      else if (polalg.eq.jacobi_SId) then ! FIXME a porter
         call inducejac_pme2gpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
#if !TINKERHP_REL_BUILD
      else if (polalg.eq.step_pcg_SId) then
         call inducestepcg_pme2gpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
      else if (polalg.eq.step_pcg_short_SId) then
         if (use_mpolelong) then
            call inducepcg_pme2gpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
         else
            call inducestepcg_pme2gpu(tmatxb_p,nrhs,.true.,ef,mu,murec)
         end if
#endif
      else
         if (rank.eq.0) write(iout,1000) polalg
         call fatal
      end if
      wtime2 = mpi_wtime()

      if (polprt.ge.1.and.rank.eq.0) then
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
      if (nmut.gt.0.and.elambda.eq.0.0) then
!$acc parallel loop collapse(3) async(def_queue) default(present)
         do i = 1, npolebloc; do j = 1, nrhs; do k = 1,3
            if (polarity(poleglob(i)).eq.0.0) mu(k,j,i) = 0.0
         end do; end do; end do
      end if
!$acc parallel loop collapse(2) async(def_queue) default(present)
      do i = 1, npolebloc; do j = 1, 3
         iipole         = poleglob(i)
         uind(j,iipole) = mu(j,1,i)
         uinp(j,iipole) = mu(j,2,i)
      end do; end do
!$acc parallel loop collapse(2) async(def_queue) default(present)
      do i = 1, npolerecloc; do j = 1, 3
         iipole = polerecglob(i)
         iglob  = ipole(iipole)
         if (repart(iglob).ne.rank) then
            uind(j,iipole) = murec(j,1,i)
            uinp(j,iipole) = murec(j,2,i)
         else
            uind(j,iipole) = mu(j,1,poleloc(iipole))
            uinp(j,iipole) = mu(j,2,poleloc(iipole))
         end if
      end do; end do
      call timer_exit( timer_other,quiet_timers )
c
c     update the lists of previous induced dipole values
c
      if (use_pred) call pred_SetAlt

      deallocate (reqrec)
      deallocate (reqsend)
      deallocate (req2rec)
      deallocate (req2send)
      deallocate (reqrecdirrec)
      deallocate (reqrecdirsend)
      end

      subroutine pred_reset
      use uprior
      nualt       = 0
      lshalt      = maxualt+1
      lalt        = maxualt+1
      end subroutine

      subroutine pred_zero
      use uprior
      use utilgpu
      implicit none
      integer(8) siz
      real(t_p) ,parameter:: val=0.0
      siz  = size(upalt)
      call mem_set(upalt,val,siz,rec_stream)
      call mem_set(udalt,val,siz,rec_stream)
      if (allocated(upshortalt)) then
         siz  = size(upshortalt)
         call mem_set(upshortalt,val,siz,rec_stream)
         call mem_set(udshortalt,val,siz,rec_stream)
      end if
      end subroutine

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
      use utilcomm ,only: skpPcomm
      use utilgpu  ,only: def_queue
      implicit none
      integer,parameter:: nrhs=2,d3=3
      real(t_p),intent(inout)::mu(d3,nrhs,npolebloc)
      integer i,j,k,k1,iipole,ierr
      integer npoleloc_e
      integer calt
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif
      real(t_p) udsum,upsum

#ifndef USE_NVSHMEM_CUDA
      if (use_pmecore.and.rank.ge.ndir) return
#endif
      if(deb_Path) print '(3x,A,2I0)',"pred_InitField",lshalt,lalt
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
         calt = lshalt-1
         !npoleloc_e= merge(npolebloc,npoleloc,(skpPcomm))
         npoleloc_e= npoleloc
      else
         udalt_p0 => udalt
         upalt_p0 => upalt
         calt = lalt-1
         npoleloc_e= npoleloc
      end if
c     write(*,'(I3,$)') (mod(calt+i-1,maxualt)+1,i=1,nualt-1)
c     write(*,*)

!$acc parallel loop collapse(2) async(def_queue)
!$acc&         present(poleglob,bpred,udalt_p0,upalt_p0,mu)
      do i = 1, npoleloc_e
         do j = 1, 3
            iipole = poleglob(i)
            udsum  = 0.0_ti_p
            upsum  = 0.0_ti_p
!$acc loop seq
            do k = 1, nualt-1
                k1   = mod(calt+k-1,maxualt) + 1
               udsum = udsum + bpred(k)*udalt_p0(j,iipole,k1)
               upsum = upsum + bpred(k)*upalt_p0(j,iipole,k1)
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
      use atoms  ,only: n_pe,n
      use domdec ,only: nproc,rank,ndir
      use inform ,only: deb_Path
      use mpole
      use mutant ,only: elambda
      use nvshmem
      use polar
      use potent ,only: use_pmecore,use_polarshortreal,use_lambdadyn
      use timestat
      use tinheader ,only: ti_p
      use uprior
      use utilgpu   ,only: def_queue
      implicit none
      integer i,j,k,iipole,ierr
      integer calt
#ifdef USE_NVSHMEM_CUDA
      integer ipe,ind
#endif

      if (use_lambdadyn.and.elambda.gt.0.0) return

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
         lshalt = lshalt - 1
         if (lshalt.eq.0) lshalt=maxualt
         calt   = lshalt
      else
         udalt_p0 => udalt
         upalt_p0 => upalt
         lalt   = lalt - 1
         if (lalt.eq.0) lalt=maxualt
         calt   = lalt
      end if

!$acc parallel loop collapse(2) async(def_queue)
!$acc&     present(poleglob,uind,uinp,upalt_p0,udalt_p0)
      do i = 1, n
         do j = 1, 3
            !iipole = poleglob(i)
            udalt_p0(j,i,calt) = uind(j,i)
            upalt_p0(j,i,calt) = uinp(j,i)
          end do
      end do

      end if
#endif
      call timer_exit( timer_ulspred )
      end subroutine
c
c
c
      module inducepcg_mod
      integer  :: indpcg_call=0
      integer  :: precnd1
      integer  :: lastIter=0,sameIter=0
      logical  :: exitlast=.false.
      real(t_p):: ggold1,ggold2
      real(r_p):: gg1,gg2,ggnew1,ggnew2,ene1,ene2
      real(t_p),target :: gbuff(4)
      real(r_p),target :: gbuff1(6)
      real(r_p),pointer:: gnorm(:), ggnew(:), ene(:)
      !real(t_p):: ggdev1=0.0,ggdev2=0.0
      !real(t_p):: alphacg1=0.0,alphacg2=0.0
      end module
      subroutine inducepcg_pme2gpu(matvec,nrhs,precnd,ef,mu,murec)
      use atmlst
      use domdec
      use ewald
      use iounit
      use inform    ,only: deb_Path,abort,minmaxone
     &              ,app_id,dynamic_a,pimd_a
      use inducepcg_mod
      use interfaces,only: tmatxb_pmegpu,pcg_a,pcg_aRec,pcg_b
     &              ,pcg_newDirection
      use math
      use mdstate   ,only: track_mds
      use mpole
      use polar
      use polpot
      use sizes
      use timestat
      use units
      use utils
      use utilcomm ,buffermpi1=>buffermpi3d,buffermpi2=>buffermpi3d1
      use utilgpu
      use polar_temp,only: res,h,pp,zr,diag
     &              , dipfield,dipfieldbis
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
c
      integer   i, it, j, k
      real(t_p) zero, one, reg1, term
      real(r_p) zerom, pt5, resnrm
      parameter(zero=0.0,zerom=0,pt5=0.5,one=1.0)
c
c     MPI
c
      integer iglob, iipole, ierr, tag, proc, sizeloc, sizebloc
      integer req1, req2, req3, req4, req5
      integer status(MPI_STATUS_SIZE)
      integer,allocatable :: reqrecdirrec(:),reqrecdirsend(:)
      integer,allocatable :: reqrec(:),reqsend(:)
      integer,allocatable :: req2rec(:),req2send(:)
      integer,allocatable :: reqendrec(:),reqendsend(:)
      integer,allocatable :: req2endrec(:),req2endsend(:)
      logical tinker_isnan_m,skip_chk
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

      call timer_enter ( timer_other )
      indpcg_call = indpcg_call+1
c
c     Allocate some memory
c
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
      gbuff1     = 0.0_re_p
      ggnew(1:2) => gbuff1(1:2)
        ene(1:2) => gbuff1(3:4)
      gnorm(1:2) => gbuff1(5:6)

      sizeloc  = 3*nrhs*npoleloc
      sizebloc = 3*nrhs*npolebloc
c
c     setup the preconditioner
c
      if (precnd) then
!$acc parallel loop default(present) async(rec_queue)
         do i = 1, npoleloc
            iipole  = poleglob(i)
            reg1    = polarity(iipole)
            diag(i) = merge(tinypol,reg1,reg1.eq.0.0_ti_p)
         end do
         if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
!$acc parallel loop default(present) async(rec_queue)
         do i = 1, npoleloc
            diag(i) = one
         end do
      end if
c
      if (indpcg_call.eq.1) then
!$acc enter data create(gg1,gg2,ggnew1,ene1,ggnew2,ene2
!$acc&          ,ggold1,ggold2,gbuff,gbuff1)
        !Set up Deflation Preconditioner
         call Tinker_shellEnv("PRECND",precnd1,0)
         !if (btest(precnd1,0)) call polarEingenVal
      end if
      !if (btest(precnd1,0)) call restartDeflat
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
!$acc serial async(rec_queue) present(ggold1,ggold2)
      ggold1 = 0.0_ti_p
      ggold2 = 0.0_ti_p
!$acc end serial

#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) call start_dir_stream_cover
#endif
      !if (btest(precnd1,0)) then
      !   call ProjectorPTransGeneric(mu,mu)
      !   call ApplyQxV(ef,mu,1.0_ti_p)
      !end if
c
      call timer_enter( timer_realdip )
      call matvec(nrhs,.true.,mu,h)
      call timer_exit ( timer_realdip,quiet_timers )
c
      call timer_enter( timer_recdip )
      call tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
      call timer_exit ( timer_recdip,quiet_timers )

      call commfieldgpu(nrhs,h)
      if (nproc.gt.1) then
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
      end if
      term = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi

      call timer_enter( timer_other )

      if (btest(precnd1,0)) then

!$acc parallel loop collapse(3) default(present) async(rec_queue)
      do k=1,npoleloc; do j=1,nrhs; do i=1,3
         h(i,j,k)   = h(i,j,k)  + dipfield(i,j,k) - term*mu(i,j,k)
         res(i,j,k) = ef(i,j,k) - h(i,j,k)
      end do; end do; end do
      call adaptedDeflat2(res,pp,diag)
!$acc parallel loop collapse(3) default(present) async(rec_queue)
!$acc&        present(ggold1,ggold2) reduction(+:ggold1,ggold2)
      do k=1,npoleloc; do j=1,nrhs; do i=1,3
         if (btest(j,0)) then; ggold1= ggold1 + res(i,j,k)*pp(i,j,k)
         else;                 ggold2= ggold2 + res(i,j,k)*pp(i,j,k)
         end if
      end do; end do; end do

      else

!$acc parallel loop collapse(3) default(present) async(rec_queue)
!$acc&              present(ggold1,ggold2) reduction(+:ggold1,ggold2)
      do k=1,npoleloc; do j=1,nrhs; do i=1,3
         h(i,j,k)   = h(i,j,k)  + dipfield(i,j,k) - term*mu(i,j,k)
         res(i,j,k) = ef(i,j,k) - h(i,j,k)
         zr(i,j,k)  = res(i,j,k)*diag(k)
         pp(i,j,k)  = zr(i,j,k)
         if (btest(j,0)) then  ! test de 1
            ggold1  = ggold1 + res(i,j,k)*zr(i,j,k)
         else
            ggold2  = ggold2 + res(i,j,k)*zr(i,j,k)
         end if
      end do; end do; end do

      end if

      if (nproc.gt.1) then
!$acc serial async present(gbuff,ggold1,ggold2)
         gbuff(1) = ggold1
         gbuff(2) = ggold2
!$acc end serial
!$acc wait(rec_queue)
!$acc host_data use_device(gbuff)
         call MPI_ALLREDUCE(MPI_IN_PLACE,gbuff,nrhs,MPI_TPREC
     &                     ,MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
!$acc serial async present(ggold1,ggold2,gbuff)
         ggold1 = gbuff(1)
         ggold2 = gbuff(2)
!$acc end serial
         !if (nproc.gt.1) call MPI_WAIT(req1,status,ierr)
      end if
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

c     -------------------------
c     now, start the main loop:
c     -------------------------
      do it = 1,politer
c
         call timer_enter( timer_realdip )
         call matvec(nrhs,.true.,pp,h)
         call timer_exit ( timer_realdip,quiet_timers )

         call timer_enter( timer_recdip )
         call tmatxbrecipgpu(pp,murec,nrhs,dipfield,dipfieldbis)
         call timer_exit ( timer_recdip,quiet_timers )

!$acc serial async present(gg1,gg2,ggnew1,ggnew2,ene1,ene2)
         gg1=zero; gg2=zero
         ggnew1 = zerom; ggnew2 = zerom
         ene1 = zerom; ene2 = zerom
!$acc end serial
c
         call commfieldgpu(nrhs,h)
c
c     Begin the reception of the reciprocal fields
c
         if (nproc.gt.1) then
         call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,reqrecdirrec,reqrecdirsend)
c
c     MPI : begin reception
c
         call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,reqrecdirrec,reqrecdirsend)
c
c       Wait for the reciprocal fields
c
         call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,reqrecdirrec,reqrecdirsend)
         end if

         term = (4.0_ti_p/3.0_ti_p)*aewald**3 / sqrtpi

         call timer_enter( timer_other )
         if (dir_queue.ne.rec_queue) then
!$acc wait(dir_queue)
            def_queue = rec_queue
         end if

         call pcg_aRec(6*npoleloc,term,pp,dipfield,h,gg1,gg2)
c!$acc host_data use_device(h,pp,dipfield)
c!$acc parallel loop collapse(3) present(gg1,gg2) async(rec_queue)
c!$acc&         deviceptr(h,dipfield,pp)
c         do k=1,npoleloc; do j=1,nrhs; do i=1,3
c            h(i,j,k) = h(i,j,k) + dipfield(i,j,k) - term*pp(i,j,k)
c            if (btest(j,0)) then
c               gg1 = gg1 + pp(i,j,k)*h(i,j,k)
c            else
c               gg2 = gg2 + pp(i,j,k)*h(i,j,k)
c            end if
c         end do; end do; end do
c!$acc end host_data

         if (nproc.gt.1) then
!$acc host_data use_device(gg1,gg2)
!$acc wait(rec_queue)
            call MPI_ALLREDUCE(MPI_IN_PLACE,gg1,1,MPI_RPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,gg2,1,MPI_RPREC,MPI_SUM,
     $           COMM_TINKER,ierr)
!$acc end host_data
         end if

         ! This check might be unecessary
         ! the only for residu this to reach zero is for the field itself to be equal to zero
         ! TODO Find a proper place for that check
         if (deb_Path) then
!$acc update host(gg1,gg2) async(rec_queue)
!$acc wait
            if (gg1.eq.zerom) then
               print*, 'Residu equals zero in polarisation solver'
               call timer_exit(timer_other)
               goto 50
            end if
            if (gg2.eq.zerom) goto 50
         end if

         ! reduce division operation in next kernel
         !alphacg1  = ggold1/gg1
         !alphacg2  = ggold2/gg2

         ggnew1 = zerom; ggnew2 = zerom
         ene1 = zerom; ene2 = zerom

         if (btest(precnd1,0)) then
         !   call loadKrylovBase(res,h,ggdev1,ggdev2)

!$acc parallel loop collapse(3) async(rec_queue)
!$acc&         present(gg1,gg2,ggold1,ggold2,ene1,ene2)
!$acc&         default(present)
         do k=1,npoleloc; do j=1,nrhs; do i=1,3
            if (btest(j,0)) then
               mu (i,j,k) = mu (i,j,k) + (ggold1/gg1)*pp(i,j,k)
               res(i,j,k) = res(i,j,k) - (ggold1/gg1)*h (i,j,k)
               ene1    = ene1 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
            else
               mu (i,j,k) = mu (i,j,k) + (ggold2/gg2)*pp(i,j,k)
               res(i,j,k) = res(i,j,k) - (ggold2/gg2)*h (i,j,k)
               ene2    = ene2 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
            end if
         end do; end do; end do

         call adaptedDeflat2(res,zr,diag)

!$acc parallel loop collapse(3) async(rec_queue)
!$acc&         present(gg1,gg2,ggnew1,ggnew2,ene1,ene2)
!$acc&         default(present)
         do k=1,npoleloc; do j=1,nrhs; do i=1,3
            if (btest(j,0)) then; ggnew1= ggnew1 + res(i,j,k)*zr(i,j,k)
            else; ggnew2 = ggnew2 + res(i,j,k)*zr(i,j,k); end if
         end do; end do; end do

         else

        call pcg_b(npoleloc,npoleloc,gg1,ggold1,gg2,ggold2
     &            ,ggnew1,ene1,ggnew2,ene2
     &            ,mu,pp,res,h,zr,diag,ef)
c!$acc host_data use_device(h,pp,dipfield,mu,res,ef,zr,diag)
c!$acc parallel loop collapse(3) async(rec_queue)
c!$acc&         present(gg1,gg2,ggold1,ggold2,ggnew1,ggnew2,ene1,ene2)
c!$acc&         deviceptr(mu,pp,res,h,zr,diag,ef)
c!$acc&         reduction(+:ggnew1,ggnew2,ene1,ene2)
c         do k=1,npoleloc; do j=1,nrhs; do i=1,3
c            if (btest(j,0)) then
c               mu (i,j,k) = mu (i,j,k) + (ggold1/gg1)*pp(i,j,k)
c               res(i,j,k) = res(i,j,k) - (ggold1/gg1)*h (i,j,k)
c               zr (i,j,k) = diag(k)* res(i,j,k)
c               ggnew1     = ggnew1 + res(i,j,k)*zr(i,j,k)
c               ene1    = ene1 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
c            else
c               mu (i,j,k) = mu (i,j,k) + (ggold2/gg2)*pp(i,j,k)
c               res(i,j,k) = res(i,j,k) - (ggold2/gg2)*h (i,j,k)
c               zr (i,j,k) = diag(k)* res(i,j,k)
c               ggnew2     = ggnew2 + res(i,j,k)*zr(i,j,k)
c               ene2    = ene2 + (mu(i,j,k)*(res(i,j,k)+ef(i,j,k)))
c            end if
c         end do; end do; end do
c!$acc end host_data

         end if

         if (nproc.gt.1) then
c Implicit wait seems to be removed with PGI-20
!$acc serial async(rec_queue) present(ene1,ene2,ggnew1,ggnew2,gbuff1)
         gbuff1(1)=ggnew1; gbuff1(2)=ggnew2
         gbuff1(3)=  ene1; gbuff1(4)=  ene2
!$acc end serial
!$acc wait(rec_queue)
!$acc host_data use_device(gbuff1)
         call MPI_ALLREDUCE(MPI_IN_PLACE,gbuff1,2*nrhs,MPI_RPREC,
     &        MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
         !call MPI_WAIT(req3,status,ierr)
         !call MPI_WAIT(req4,status,ierr)
!$acc serial async(rec_queue) present(ggnew1,ggnew2,ene1,ene2,gbuff1)
         ggnew1 = gbuff1(1); ggnew2 = gbuff1(2)
         ene1   = gbuff1(3); ene2   = gbuff1(4)
!$acc end serial
         end if

        !Compute Next direction
         call pcg_newDirection(6*npoleloc,ggnew1,ggold1,ggnew2,ggold2
     &                      ,zr,pp)
c
         call timer_exit( timer_other )
c
         ! Skip Iteration check
         skip_chk = (it.lt.2.or.(exitlast.and.(it.lt.lastIter)))
     &              .and.polprt.eq.0
         ! Skip last Iteration check
         if (exitlast.and.it.eq.lastIter.and.mod(sameIter,3).ne.0
     &      .and..not.track_mds)
     &      goto 10
c
         call commdirdirgpu(nrhs,1,pp,reqrec,reqsend)
         call commdirdirgpu(nrhs,0,pp,reqrec,reqsend)
 
         ! Begin reception of mu for PME
         call commrecdirdipgpu(nrhs,0,murec,pp,buffermpimu1,
     &        buffermpimu2,req2rec,req2send)
c
!$acc serial async(rec_queue)
!$acc&       present(ggold1,ggold2,ggnew1,ggnew2,ene1,ene2,gbuff1)
         ggold1 = ggnew1
         ggold2 = ggnew2
         if (nproc.eq.1) then
            gbuff1(1)= ggnew1; gbuff1(2)= ggnew2
            gbuff1(3)= ene1  ; gbuff1(4)= ene2  
         end if
!$acc end serial

         if (skip_chk) then
            resnrm = 2*poleps
         else
!$acc update host(gbuff1) async(rec_queue)
!$acc wait(rec_queue)
            do k = 1, nrhs
               gnorm(k) = sqrt(ggnew(k)/real(3*npolar,r_p))
            end do
            resnrm = max(gnorm(1),gnorm(2))
c
            ene = -pt5*ene
            if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     &         it, (ene(k)*coulomb, gnorm(k), k = 1, nrhs)
         end if
c
         call commdirdirgpu(nrhs,2,pp,reqrec,reqsend)
         call commrecdirdipgpu(nrhs,1,murec,pp,buffermpimu1,
     &        buffermpimu2,req2rec,req2send)
         call commrecdirdipgpu(nrhs,2,murec,pp,buffermpimu1,
     &        buffermpimu2,req2rec,req2send)
c
         if ((it.eq.politer.and.resnrm.gt.poleps)
     &      .or.tinker_isnan_m(gnorm(1))) then
            ! Not converged Exit
            if (rank.eq.0) write(0,1050) it,max(gnorm(1),gnorm(2))
            abort = .true.
            goto 10
         end if
c
         if (resnrm.lt.poleps) then
            if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     &         (ene(k)*coulomb, k = 1, nrhs), (gnorm(k), k = 1, nrhs)
            goto 10
         end if
      end do
 10   continue

      ! Save and Check Convergence Iteration Index
      if (app_id.eq.dynamic_a.or.app_id.eq.pimd_a) then
      sameIter = merge(0,sameIter+1,lastIter.ne.it)
      lastIter = it
      exitlast = polprt.eq.0.and.sameIter.gt.25
     &           .and.mod(indpcg_call,100).ne.0
      if (deb_Path) print*, 'pcg_conv iter',lastIter,sameIter,exitlast
      end if

c
c     Begin reception of mu for PME
c
      call commrecdirdipgpu(nrhs,0,murec,mu,buffermpimu1,
     &     buffermpimu2,req2endrec,req2endsend)
c
c     MPI : begin direct space communication
c
      call commdirdirgpu(nrhs,0,mu,reqendrec,reqendsend)
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
      !if (btest(precnd1,0)) then
      !   call FinalizeKrylovBase
      !end if

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
#ifdef _OPENACC
          call fatal_device("inducejac_pme2gpu")
#else
          call M_gesv(nmat,1,bloc,ndismx+1,ipiv,cex,nmat,info)
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
#if 1
      module step_cg
#if TINKER_DOUBLE_PREC
      integer,parameter:: cgstep_max=20
#else
      integer,parameter:: cgstep_max=15
#endif
      integer s,sshort,slong
      integer,private:: step_save=0
      real(t_p),allocatable,target:: sb0(:)
      real(r_p),allocatable,target:: sb1(:)
      real(t_p),allocatable,target:: stepMatV(:)
      real(r_p),allocatable,target:: stepMatVr(:)
      real(t_p),pointer:: MatBs(:,:,:),MatWs(:,:,:)
     &         ,MatW_s(:,:,:),Vas(:,:)
      real(r_p),pointer:: MatBs_r(:,:,:),MatWs_r(:,:,:)
     &         ,MatW_sr(:,:,:),Vas_r(:,:)
      real(t_p),pointer:: MRes_s(:,:,:,:),MRes_s0(:,:,:)
     &         ,MQs(:,:,:,:),MTs(:,:,:,:)
      real(r_p),pointer:: MPs(:,:,:,:)
      integer ipiv(cgstep_max)

      ! One macro-iteration optimisation data set
      integer MaxIt1InARow,short_m,normal_m
      integer indStepCall(2),OneItInARow(2),MoreInARow(2)
      logical optOneIt(2)
      parameter(MaxIt1InARow=20,short_m=2,normal_m=1)

      ! Shift attributes
      logical ::shift_enabled=.true.
      real(t_p),dimension(0:cgstep_max):: sh_the
     &         ,sh_sig,sh_gam

      contains

      subroutine man_stepWorkSpace(nrhs)
      implicit none
      integer,intent(in)::nrhs
      integer lsize,s2

      if (step_save.ne.s) then
        s2    = s*s
        lsize = nrhs*s*(3*s+1)
        step_save = s
        if (associated(MatBs)) then
!$acc exit data detach(MatBs,MatWs,MatW_s,Vas
!$acc&         ,MatBs_r,MatWs_r,MatW_sr,Vas_r)
        end if
        MatBs (1:s,1:s,1:nrhs) => stepMatV(1:nrhs*s2)
        MatWs (1:s,1:s,1:nrhs) => stepMatV(  nrhs*s2+1:2*nrhs*s2)
        Vas   (1:s,1:nrhs)     => stepMatV(2*nrhs*s2+1:nrhs*s*(2*s+1))
        MatW_s(1:s,1:s,1:nrhs) => stepMatV(nrhs*s*(2*s+1)+1:lsize)
        MatBs_r (1:s,1:s,1:nrhs) => stepMatVr(1:nrhs*s2)
        MatWs_r (1:s,1:s,1:nrhs) => stepMatVr(  nrhs*s2+1:2*nrhs*s2)
        Vas_r   (1:s,1:nrhs)   => stepMatVr(2*nrhs*s2+1:nrhs*s*(2*s+1))
        MatW_sr (1:s,1:s,1:nrhs) => stepMatVr(nrhs*s*(2*s+1)+1:lsize)
!$acc enter data attach(MatBs,MatWs,MatW_s,Vas
!$acc&          ,MatBs_r,MatWs_r,MatW_sr,Vas_r)
      end if
      end subroutine
c     subroutine select_step(iter)
c     implicit none
c     integer,intent(in):: iter
c     if (iter.eq.2) then
c        s = 4
c     else
c        s = cgstep_max
c     end if
c     end subroutine

      ! Shift methods
      subroutine init_shift_attr
      implicit none
      sh_the = 0
      sh_gam = 1
      sh_sig = 0
      end subroutine

      subroutine shift_field(elecf,vec_1,vec_2,sz,s)
      implicit none
      integer,intent(in):: s,sz
      real(t_p),intent(in):: vec_1(*),vec_2(*)
      real(t_p),intent(inout):: elecf(*)

      integer i, offset,start,quo
      real(t_p) theta0_gam,inv_gamma0,theta,sig_gam

      start = sz/2
      if (s.eq.1) then
         theta0_gam = sh_the(0)/sh_gam(0)
         inv_gamma0 = 1.0/sh_gam(0)
!$acc parallel loop async default(present)
         do i = 1,sz
            quo = (i-1)/3
            if (btest(quo,0)) then
               offset = start
               quo    = (((quo-1)/2)*3) + mod(i-1,3) + 1
            else
               offset = 0
               quo    = (( quo   /2)*3) + mod(i-1,3) + 1
            end if
           elecf(i) = inv_gamma0*elecf(i) - theta0_gam*vec_1(offset+quo)
         end do
      else
         inv_gamma0 = 1.0/sh_gam(s-1)
         theta      =  sh_the(s-1)/sh_gam(s-1)
         sig_gam    = -sh_sig(s-2)/sh_gam(s-1) 
!$acc parallel loop async default(present)
         do i = 1,sz
            quo = (i-1)/3
            if (btest(quo,0)) then
               offset = start
               quo    = (((quo-1)/2)*3) + mod(i-1,3) + 1
            else
               offset = 0
               quo    = (( quo   /2)*3) + mod(i-1,3) + 1
            end if
            elecf(i) = inv_gamma0*elecf(i) -   theta*vec_1(offset+quo)
     &                                     + sig_gam*vec_2(offset+quo)
         end do
      end if

      end subroutine

      subroutine norm_field(elecf,sz,s,normd,normp)
      implicit none
      integer,intent(in):: s,sz
      real(r_p) normd,normp
      real(t_p),intent(inout):: elecf(*)
      integer i,offset,start,quo

      normd=0; normp=0;
      start = sz/2
!$acc parallel loop async default(present)
      do i = 1,sz
         quo = (i-1)/3
         if (btest(quo,0)) then
            offset = start
            quo    = (((quo-1)/2)*3) + mod(i-1,3) + 1
            normd  = normd + elecf(i)*elecf(i)
         else
            offset = 0
            quo    = ((quo/2)*3) + mod(i-1,3) + 1
            normp  = normp + elecf(i)*elecf(i)
         end if
      end do
!$acc wait
      normp=sqrt(normp); normd=sqrt(normd)
!$acc parallel loop async default(present)
      do i = 1,sz
         quo = (i-1)/3
         if (btest(quo,0)) then
            offset = start
            quo    = (((quo-1)/2)*3) + mod(i-1,3) + 1
            elecf(i) = elecf(i)/normd
         else
            offset = 0
            quo    = ((quo/2)*3) + mod(i-1,3) + 1
            elecf(i) = elecf(i)/normp
         end if
      end do
      end subroutine

      subroutine orth_field(elecf,vec_1,vec_2,sz,s)
      implicit none
      integer  ,intent(in)   :: s,sz
      real(t_p),intent(in)   :: vec_1(*),vec_2(*)
      real(t_p),intent(inout):: elecf(*)

      integer i, offset,start,quo
      real(t_p) theta0_gam,inv_gamma0,theta,sig_gam
      real(r_p),save:: scl0,scl1,l2nrm0,l2nrm1
      real(t_p)     :: scald,scalp,l2normd,l2normp

      start = sz/2
      ! Scalar Product
      scl0=0; scl1=0
!$acc parallel loop async default(present)
      do i = 1,sz
         quo = (i-1)/3
         if (btest(quo,0)) then
            offset = start
            quo    = (((quo-1)/2)*3) + mod(i-1,3) + 1
            scl0   = scl0 + elecf(i)*vec_1(offset+quo)
         else
            offset = 0
            quo    = ((quo/2)*3) + mod(i-1,3) + 1
            scl1   = scl1 + elecf(i)*vec_1(offset+quo)
         end if
      end do
!$acc wait
      ! Orthogonalization (GS)
      scald=scl0; scalp=scl1;
      if (s.eq.1) l2nrm0=0
      if (s.eq.1) l2nrm1=0
      l2normd = l2nrm0
      l2normp = l2nrm1
!$acc parallel loop async default(present)
      do i = 1,sz
         quo = (i-1)/3
         if (btest(quo,0)) then
            offset = start
            quo    = (((quo-1)/2)*3) + mod(i-1,3) + 1
            elecf(i) = elecf(i) - scald*vec_1(offset+quo)
     &                        - l2normd*vec_2(offset+quo)
         else
            offset = 0
            quo    = ((quo/2)*3) + mod(i-1,3) + 1
            elecf(i) = elecf(i) - scalp*vec_1(offset+quo)
     &                        - l2normp*vec_2(offset+quo)
         end if
      end do
      call norm_field(elecf,sz,s,l2nrm0,l2nrm1)
      if (s.eq.1) then
         ! Scalar Product
         scl0=0; scl1=0
!$acc parallel loop async default(present)
         do i = 1,sz
            quo = (i-1)/3
            if (btest(quo,0)) then
               offset = start
               quo    = (((quo-1)/2)*3) + mod(i-1,3) + 1
               scl0   = scl0 + elecf(i)*vec_1(offset+quo)
            else
               offset = 0
               quo    = ((quo/2)*3) + mod(i-1,3) + 1
               scl1   = scl1 + elecf(i)*vec_1(offset+quo)
            end if
         end do
!$acc wait
         if (abs(scl0).gt.1d-6.or.abs(scl1).gt.1d-6) then
            print*, s,'ortho check',scl0,scl1
         end if
      else
         ! Scalar Product
         scl0=0; scl1=0
!$acc parallel loop async default(present)
         do i = 1,sz
            quo = (i-1)/3
            if (btest(quo,0)) then
               offset = start
               quo    = (((quo-1)/2)*3) + mod(i-1,3) + 1
               scl0   = scl0 + elecf(i)*vec_2(offset+quo)
            else
               offset = 0
               quo    = ((quo/2)*3) + mod(i-1,3) + 1
               scl1   = scl1 + elecf(i)*vec_2(offset+quo)
            end if
         end do
!$acc wait
         if (abs(scl0).gt.1d-6.or.abs(scl1).gt.1d-6) then
            print*, s,'ortho check',scl0,scl1
         end if
      end if

      end subroutine

      ! Utility methods
      subroutine amove(n,a,b,queue)
      implicit none
      integer n,queue
      real(t_p) a(*),b(*)
      integer i
!$acc parallel loop async(queue) present(a,b)
      do i = 1,n
         b(i) = a(i)
      end do
      end subroutine
      subroutine amover(n,a,b,queue)
      implicit none
      integer n,queue
      real(r_p) a(*),b(*)
      integer i
!$acc parallel loop async(queue) present(a,b)
      do i = 1,n
         b(i) = a(i)
      end do
      end subroutine
      subroutine amovec(n,a,b,queue)
      implicit none
      integer n,queue
      real(r_p) a(*)
      real(t_p) b(*)
      integer i
!$acc parallel loop async(queue) present(a,b)
      do i = 1,n
         b(i) = a(i)
      end do
      end subroutine
      subroutine amovec1(n,a,b,queue)
      implicit none
      integer n,queue
      real(t_p) a(*)
      real(r_p) b(*)
      integer i
!$acc parallel loop async(queue) present(a,b)
      do i = 1,n
         b(i) = a(i)
      end do
      end subroutine
      subroutine prod_scal(n,a,b,res)
      implicit none
      integer n
      real(t_p) a(*),b(*)
      real(r_p) res
      integer i
!$acc routine vector
!$acc loop vector reduction(res)
      do i = 1,n
         res = res + a(i)*b(i)
      end do
      end subroutine
      subroutine prod_scal1(n,a,b,res)
      implicit none
      integer n
      real(t_p) a(*)
      real(r_p) b(*)
      real(r_p) res
      integer i
!$acc routine vector
!$acc loop vector reduction(res)
      do i = 1,n
         res = res + a(i)*b(i)
      end do
      end subroutine
      subroutine prod_scal2(n,a,b,res)
      implicit none
      integer n
      real(r_p) a(*)
      real(t_p) b(*)
      real(r_p) res
      integer i
!$acc routine vector
!$acc loop vector reduction(res)
      do i = 1,n
         res = res + a(i)*b(i)
      end do
      end subroutine
      subroutine prod_scal_n(n,a,b,res)
      implicit none
      integer n
      real(t_p) a(*),b(*)
      real(r_p) res
      integer i
!$acc routine vector
!$acc loop vector reduction(res)
      do i = 1,n
         res = res - a(i)*b(i)
      end do
      end subroutine
      subroutine prod_scal_n1(n,a,b,res)
      implicit none
      integer n
      real(t_p) a(*)
      real(r_p) b(*)
      real(r_p) res
      integer i
!$acc routine vector
!$acc loop vector reduction(res)
      do i = 1,n
         res = res - a(i)*b(i)
      end do
      end subroutine
      end module

      subroutine load_shift_attr(theta,sigma,k)
      use domdec ,only: rank,nproc,COMM_TINKER
      use step_cg
      use sizes  ,only: tinkerdebug
      use mpi
      implicit none
      integer,intent(in)::k
      real(r_p) theta(k),sigma(k)
      integer i,j

      do i = 0,min(k-1,cgstep_max)
         sh_the(i) = theta(i+1)
c        sh_sig(i) = sigma(i+1)
      end do

c     sh_the(1) = theta(1)

      if (tinkerdebug.gt.0) then
         if(rank.eq.0) 
     &      print*, ' ***** Step PCG base Shift Enabled ***** '

         do i = 1,nproc
            write(0,'(A,I0,$)')'theta ',rank
            write(0,'(F8.4,$)')(sh_the(j),j=0,cgstep_max)
            write(0,*) ''
            call MPI_BARRIER(COMM_TINKER,j)
         end do
      end if
      !print*,'sigma',(real(sh_sig(i),4),i=0,cgstep_max)
      end subroutine

      subroutine inducestepcg_pme2gpu(matvec,nrhs,precnd,ef,mu,murec)
      use atmlst
      use domdec
      use ewald
      use iounit
      use inform    ,only: deb_Path,abort,minmaxone
      use interfaces,only: tmatxb_pmegpu
#ifdef _OPENACC
     &              ,initcuSolver,cuGESVm
#endif
      use math
      use mpole
      use polar
      use polpot
      use sizes
      use step_cg
      use timestat
      use units
      use utils
      use utilcomm ,buffermpi1=>buffermpi3d,buffermpi2=>buffermpi3d1
      use utilgpu
      use polar_temp,only: res,h,pp,zr,diag
     &              , dipfield,dipfieldbis
      use potent
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

      integer i, it, j, k, l, ll, lc, its, it2, mode
      integer sdec,ndec,ndec1,stdec,ledec
      real(r_p) gnorm(2), ggnew(2), ene(2), gg0(2)
      real(r_p),save:: ggold1=0.0,ggold2=0.0
      real(r_p),save:: ggnew1=0.0,ggnew2=0.0
      real(r_p),save:: gnorm1=0.0,gnorm2=0.0
      real(r_p),save:: ene1=0.0,ene2=0.0
      real(t_p) zero, one
      real(r_p) zerom
      real(r_p) resnrm
      real(r_p) tempo
      real(t_p) term,prcdt,rt
      real(r_p) rest(cgstep_max),row(cgstep_max),rtr
      parameter( zero=0.0,zerom=0,one=1.0 )
c
c     MPI
c
      integer iglob, iipole, ierr, tag, proc, sizeloc, sizebloc
      integer req1, req2, req3, req4, req5
      integer status(MPI_STATUS_SIZE)
      integer,dimension(nproc) :: reqrecdirrec,reqrecdirsend
     &       ,reqrec,reqsend,req2rec,req2send
     &       ,reqendrec,reqendsend,req2endrec,req2endsend
      logical,save::first_in=.true.
      logical tinker_isnan_m
      integer info 
      logical d_MQR !(dinstinct space between MRes_s & MQ)
c
 1000 format(' stepcgiter converged after ',I3,' iterations.',/,
     &       ' final residual norm = ',2D14.7)
 1010 format(' residual norm at iteration ',I3,':',2D12.3,1x,2D12.3)
 1020 format(' Conjugate gradient solver: induced dipoles',/,
     &  ' ipole       mux         muy         muz')
 1021 format(' Conjugate gradient solver: induced p-dipoles',/,
     &  ' ipole       mux         muy         muz')
 1030 format(i6,2x,f10.7,2x,f10.7,2x,f10.7)
 1040 format(' Using a diagonal preconditioner.')
 1050 format(' Induce step PCG solver not converged after '
     &      ,I6,' iterations',/,' Residual norm ', 2d10.5)

      if (deb_Path) 
     &   write(*,'(3x,a)') 'inducestepcg_pme2gpu'
      call timer_enter ( timer_stepcg )
c
c     Allocate some memory
c
      if (nproc.gt.1) then
      call prmem_request(buffermpi1,3,nrhs,max(npoleloc,1))
      call prmem_request(buffermpimu1,3,nrhs,max(npoleloc,1))
      call prmem_request(buffermpi2,3,nrhs,max(npolerecloc,1))
      call prmem_request(buffermpimu2,3,nrhs,max(npolerecloc,1))
      end if

      call prmem_request(dipfield,3,nrhs,max(npoleloc,1))
      call prmem_request(dipfieldbis,3,nrhs,max(npolerecloc,1))
      call prmem_request(zr,3,nrhs,max(npoleloc,1))
      call prmem_request(res,3,nrhs,max(npolebloc,1))
      call prmem_request(h ,3,nrhs,max(npolebloc,1))
c     call prmem_request(pp ,3,nrhs,max(npolebloc,1))
      call prmem_request(diag,npoleloc)
c
      if (first_in) then
#ifdef _OPENACC
         call initcuSolver(rec_stream)
#endif
!$acc enter data create(ipiv)
         call Tinker_shellEnv('SSDIM',sshort,4)
         call Tinker_shellEnv('SDIM',slong,5)
         sshort = min(sshort,cgstep_max)
         slong  = min(slong ,cgstep_max)
c
         if (shift_enabled) then
            call init_shift_attr
            call polarEingenVal
         end if

         indStepCall(:) = 0
         OneItInARow(:) = 0
         MoreInARow(:)  = 0
         optOneIt(:)    = .false.

         first_in = .false.
      end if
c
      mode     = merge(short_m,normal_m,use_polarshortreal)
      indStepCall(mode) = indStepCall(mode) + 1

      d_MQR    = merge(.true.,.false.,precnd.or.shift_enabled)
      ! step WorkSpace buffer segmentation
      s        = merge( sshort,slong,use_polarshortreal )
c
      req1     = merge(3*npoleloc*nrhs*(2*s+1),3*npoleloc*nrhs*(s+1)
     &                ,d_MQR)
      req2     = 3*npoleloc*nrhs*s
      req3     = nrhs*s*(3*s+1)
      sizeloc  = 3*nrhs*npoleloc
      sizebloc = 3*nrhs*npolebloc
      sdec     = 128
      ndec1    = 3*npoleloc/sdec
      ndec     = merge(ndec1+1,ndec1,mod(3*npoleloc,sdec).gt.0)
c
      call prmem_request(sb0,req1)
      call prmem_requestm(sb1,req2)
      call prmem_request(stepMatV ,req3)
      call prmem_requestm(stepMatVr,req3)

      ! WorkSpace buffer segmentation
      if (d_MQR) then
      req4 = req2+sizeloc
      MTs (1:3,1:npoleloc,1:nrhs,1:2*s+1) => sb0(1:req1)
      MRes_s(1:3,1:npoleloc,1:nrhs,1:s) => sb0(sizeloc+1:req4)
      MQs   (1:3,1:npoleloc,1:nrhs,1:s) => sb0(req4+1:req1)
      d_MQR=.true.
      else
      MTs (1:3,1:npoleloc,1:nrhs,1:s+1) => sb0(1:req1)
      MRes_s(1:3,1:npoleloc,1:nrhs,1:s) => sb0(1:req2)
      MQs   (1:3,1:npoleloc,1:nrhs,1:s) => sb0(sizeloc+1:req1)
      d_MQR=.false.
      end if
      MRes_s0(1:3,1:npoleloc,1:nrhs   ) => sb0(1:sizeloc)
      MPs   (1:3,1:npoleloc,1:nrhs,1:s) => sb1(1:req2)
!$acc enter data attach(MRes_s,Mres_s0,MQs,MTs,MPs)

      call man_stepWorkSpace(nrhs)
c

!$acc data present(diag,zr,res,h,MRes_s,MQs,MPs,MTs
!$acc&            ,ef,mu,dipfield)
c
c     setup the preconditioner
c
      if (precnd) then
!$acc parallel loop present(poleglob,polarity) async(rec_queue)
         do i = 1, npoleloc
            iipole  = poleglob(i)
            diag(i) = polarity(poleglob(i))
         end do
         if (polprt.ge.2.and.rank.eq.0) write (iout,1040)
      else
!$acc parallel loop async(rec_queue)
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
      !call set_to_zero1(pp,3*nrhs*npolebloc  ,rec_queue)
      else
      call set_to_zero5(buffermpi1,buffermpimu1,res,zr,dipfield,
     &                     3*nrhs*npoleloc   ,rec_queue)
      call set_to_zero3(buffermpimu2,buffermpi2,dipfieldbis,
     &                     3*nrhs*npolerecloc,rec_queue)
      !call set_to_zero1(pp,3*nrhs*npolebloc  ,rec_queue)
      end if
      call timer_exit( timer_stepcg,quiet_timers )
c
#ifdef _OPENACC
      if (use_polarshortreal.and.(dir_queue.ne.rec_queue))
     &   call start_dir_stream_cover
#endif

      ! compute the initial residu
      call timer_enter( timer_realdip )
      call matvec(nrhs,.true.,mu,h)
      call timer_exit ( timer_realdip,quiet_timers )
c
      if (.not.use_polarshortreal) then
      call timer_enter( timer_recdip )
      call tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
      call timer_exit ( timer_recdip,quiet_timers )
      end if

      if (use_polarshortreal) then
         call commfieldshort(nrhs,h)
      else
         call commfieldgpu(nrhs,h)

         call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $        buffermpi2,reqrecdirrec,reqrecdirsend)
         call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $        buffermpi2,reqrecdirrec,reqrecdirsend)
         call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $        buffermpi2,reqrecdirrec,reqrecdirsend)
      end if

      call timer_enter( timer_stepcg )

      term   = (4.0_ti_p/3.0_ti_p) * aewald**3 / sqrtpi
      ggold1 = 0.0_re_p; ggold2 = 0.0_re_p
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif
c
!$acc parallel loop collapse(3) async(rec_queue)
      do k=1,npoleloc; do j=1,nrhs; do i=1,3
         if (.not.use_polarshortreal)
     &   h(i,j,k) = h(i,j,k)  + dipfield(i,j,k) - term*mu(i,j,k)
         rt       = (ef(i,j,k) - h(i,j,k))
         !pp(i,j,k)  = zr(i,j,k)
         !if (btest(j,0)) then  ! test de 1
         !   ggold1  = ggold1 + diag(k)*rt**2
         !else
         !   ggold2  = ggold2 + diag(k)*rt**2
         !end if
         res(i,j,k) = rt
      end do; end do; end do
!!$acc wait(rec_queue)

c     gg0(1) = ggold1; gg0(2) = ggold2
c     if (nproc.gt.1) then
c        call MPI_IALLREDUCE(MPI_IN_PLACE,gg0(1),nrhs,MPI_RPREC
c    &                      ,MPI_SUM,COMM_TINKER,req1,ierr)
c     end if
      call timer_exit ( timer_stepcg,quiet_timers )
c
c     now, start the main loop:
c
      do it = 1,politer
         if(d_MQR) then
!$acc parallel loop collapse(3) async(rec_queue)
            do j=1,nrhs; do k=1,npoleloc; do i=1,3
               rt  = res(i,j,k)
               Mres_s0(i,k,j) = rt
               res(i,j,k)     = rt*diag(k)
            end do; end do; end do
         end if
c        call norm_field(res,6*npoleloc,its,ene1,ene2)
c
         ! Store first residu in tall matrix and apply preconditionner
!$acc parallel loop collapse(3) async(rec_queue)
         do j=1,nrhs; do k=1,npoleloc; do i=1,3
            Mres_s(i,k,j,1) = res(i,j,k)
         end do; end do; end do
c        call minmaxone(Mres_s(1,1,1,1),6*npoleloc,'res0')

#ifdef _OPENACC
         if (dir_queue.ne.rec_queue) call start_dir_stream_cover
#endif
         !------------------------
         !  Matrix-Power Kernels |
         !------------------------
         do its = 1,s

         ! get neighbor contribution of field
         if (use_polarshortreal) then
            if (.not.skpPcomm) then
            call commdirdirshort(nrhs,1,res,reqrec,reqsend)
            call commdirdirshort(nrhs,0,res,reqrec,reqsend)
            call commdirdirshort(nrhs,2,res,reqrec,reqsend)
            end if
         else
            call commdirdirgpu(nrhs,0,res,reqrec,reqsend)
            call commdirdirgpu(nrhs,1,res,reqrec,reqsend)
            call commdirdirgpu(nrhs,2,res,reqrec,reqsend)

            call commrecdirdipgpu(nrhs,0,murec,res,buffermpimu1,
     &           buffermpimu2,req2rec,req2send)
            call commrecdirdipgpu(nrhs,1,murec,res,buffermpimu1,
     &           buffermpimu2,req2rec,req2send)
            call commrecdirdipgpu(nrhs,2,murec,res,buffermpimu1,
     &           buffermpimu2,req2rec,req2send)
         end if
c
         ! ------- Matrix vector --------
         call timer_enter( timer_realdip )
         call matvec(nrhs,.true.,res,h)
         call timer_exit ( timer_realdip,quiet_timers )

         if (.not.use_polarshortreal) then
         call timer_enter( timer_recdip )
         call tmatxbrecipgpu(res,murec,nrhs,dipfield,dipfieldbis)
         call timer_exit ( timer_recdip,quiet_timers )
         end if
         ! ------------------------------

         ! Add local missing interactions
         if (use_polarshortreal) then
            call commfieldshort(nrhs,h)
         else
            call commfieldgpu(nrhs,h)
 
         !Communicate the reciprocal fields
         call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,reqrecdirrec,reqrecdirsend)
         call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,reqrecdirrec,reqrecdirsend)
         call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     &        buffermpi2,reqrecdirrec,reqrecdirsend)
         end if

         call timer_enter( timer_stepcg )
         term = (4.0_ti_p/3.0_ti_p)*aewald**3 / sqrtpi

         ! Store result in tall matrix
         if (use_polarshortreal) then
!$acc parallel loop collapse(3) async(dir_queue)
            do j=1,nrhs; do k=1,npoleloc; do i=1,3
               rt        = h(i,j,k)
               h(i,j,k)  =  rt*diag(k)
               MQs(i,k,j,its) = rt
            end do; end do; end do

            if (its.ne.s) then  !Unecessary on last
            if (shift_enabled) then
            call shift_field(h,Mres_s(1:,1,1,its),MRes_s(1:,1,1,its-1)
     &                        ,3*nrhs*npoleloc,its)
            end if
!$acc parallel loop collapse(3) async(dir_queue)
            do j=1,nrhs; do k=1,npoleloc; do i=1,3
               rt          = h(i,j,k)
               res(i,j,k)  = rt
               if (d_MQR) MRes_s(i,k,j,its+1)=rt
            end do; end do; end do
            end if
         else
!$acc parallel loop collapse(3) async(rec_queue)
            do j=1,nrhs; do k=1,npoleloc; do i=1,3
               rt = h(i,j,k) + dipfield(i,j,k) - term*res(i,j,k)
               h(i,j,k)       = rt*diag(k)
               MQs(i,k,j,its) = rt
            end do; end do; end do
            if (its.ne.s) then
            if (shift_enabled) then
            call shift_field(h,MRes_s(1:,1,1,its),MRes_s(1:,1,1,its-1)
     &                        ,3*nrhs*npoleloc,its)
            end if
!$acc parallel loop collapse(3) async(rec_queue)
            do j=1,nrhs; do k=1,npoleloc; do i=1,3
               rt             = h(i,j,k)
               res(i,j,k)     = rt
               if (d_MQR) MRes_s(i,k,j,its+1) = rt
            end do; end do; end do
            end if
         end if
         call timer_exit ( timer_stepcg,quiet_timers )
c        if (its.ne.s)
c    &   call minmaxone(Mres_s(1,1,1,its+1),6*npoleloc,'res')

         end do  ! Matrix Power loop
#ifdef _OPENACC
         if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif

         call timer_enter( timer_stepcg )
         ! ReSet step Matrix & Vector
!$acc parallel loop async(rec_queue) present(stepMatV,stepMatVr)
         do i = 1,nrhs*s*(2*s+1)
            stepMatV(i) = 0.0
           stepMatVr(i) = 0.0
         end do

         if (it.eq.1) then
            call amovec1(3*nrhs*npoleloc*s,MRes_s,MPs,rec_queue)
            !call minmaxone(MPs(1,1,1,1),6*npoleloc*s,'MPs')
         else
            ! Block dot-products
            ! Compute  C_i = -Q_i^T * P_{i-1} in MatBs
!$acc parallel loop gang collapse(4) async(rec_queue) present(MatBs_r)
            do j = 1,nrhs; do ll = 1,s; do lc = 1,s; do i = 1,ndec;
               tempo = 0
               stdec = (i-1)*sdec+1
               ledec = min(sdec,3*npoleloc-stdec+1)
               call prod_scal_n1(ledec,MQs(stdec:,1,j,ll)
     &                         ,MPs(stdec:,1,j,lc),tempo)
!$acc loop vector
               do k = 1,32
                  if (k.eq.1) then
!$acc atomic
                    MatBs_r(lc,ll,j) = MatBs_r(lc,ll,j) + tempo
                  end if
               end do
            end do; end do; end do; end do;
            if (nproc.gt.1) then
!$acc host_data use_device(stepMatVr)
!$acc wait(rec_queue)
               call MPI_ALLREDUCE(MPI_IN_PLACE,stepMatVr,nrhs*(s**2)
     $                           ,MPI_RPREC,MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
            end if

            ! Find scalars B
            ! Resolve W*B = C_i
            do j = 1, nrhs
#ifdef _OPENACC
!$acc host_data use_device(MatW_sr,ipiv,MatBs_r)
               call cuGESVm(s,s,MatW_sr(1,1,j),s,ipiv,MatBs_r(1,1,j),s
     &                     ,rec_stream)
!$acc end host_data
#else
               call M_gesvm(s,s,MatW_sr(1:,1,j),s,ipiv,MatBs_r(1:,1,j)
     &                     ,s,info)
#endif
            end do
            call amovec(nrhs*s*s,stepMatVr,stepMatV,rec_queue)

#if 0
            if (rank.eq.0) then
!$acc wait
!$acc update host(matBs_r)
            write(*,'(A)') 'matBs'
 16         format(D14.6,$)
            do j=1,s
            write(*,16) (matBs_r(j,i,1),i=1,s),(matBs_r(j,i,2),i=1,s)
            write(*,*)
            end do
            end if
#endif

           !  Next s Block Direction
           !  P_i = R + P_{i-1}*B
!$acc parallel loop gang vector collapse(3) private(rest,row)
!$acc&         async(rec_queue) present(MatBs)
            do j = 1,nrhs; do k = 1,npoleloc; do i = 1,3
!$acc loop seq
               do its = 1,s
                  row(its) = MPs(i,k,j,its)
                 rest(its) = 0
               end do
               do its = 1,s
                  do lc = 1,s
                  rest(its) = rest(its)+ row(lc)*MatBs_r(lc,its,j)
                  end do
                  rest(its) = rest(its) + MRes_s(i,k,j,its)
               end do
!$acc loop seq
               do its = 1,s
                  MPs(i,k,j,its) = rest(its)
               end do
            end do; end do; end do
            !call minmaxone(MPs(1,1,1,1),6*npoleloc*s,'MPs')
            !call minmaxone(MRes_s(1,1,1,1),6*npoleloc*s,'MRes_s')

         end if

         ! Second Block dot-Products
         ! W_i = Q_i^T * P_i
!$acc parallel loop gang collapse(4) async(rec_queue) present(MatWs_r)
         do j = 1,nrhs; do ll = 1,s; do lc = 1,s; do i = 1,ndec
            tempo = 0
            stdec = (i-1)*sdec+1
            ledec = min(sdec,3*npoleloc-stdec+1)
            call prod_scal1(ledec,MQs(stdec:,1,j,ll)
     &                    ,MPs(stdec:,1,j,lc),tempo)
!$acc loop vector
            do k = 1,32
               if (k.eq.1) then
!$acc atomic
                  MatWs_r(lc,ll,j) = MatWs_r(lc,ll,j) + tempo
               end if
            end do
         end do; end do; end do; end do

         ! Compute
         ! g_i = P_i^T * res_i  (use container for vector a)
!$acc parallel loop gang collapse(3) async(rec_queue) present(Vas_r)
         do j = 1,nrhs; do ll = 1,s; do i = 1,ndec
            tempo = 0
            stdec = (i-1)*sdec+1
            ledec = min(sdec,3*npoleloc-stdec+1)
            call prod_scal2(ledec,MPs(stdec:,1,j,ll)
     &                    ,MRes_s0(stdec:,1,j),tempo)
!$acc loop vector
            do k = 1,32
               if (k.eq.1) then
!$acc atomic
                  Vas_r(ll,j) = Vas_r(ll,j) + tempo
               end if
            end do
         end do; end do; end do

         if (nproc.gt.1) then
!$acc host_data use_device(stepMatVr)
!$acc wait(rec_queue)
            call MPI_ALLREDUCE(MPI_IN_PLACE,stepMatVr(nrhs*s**2+1)
     $          ,nrhs*s*(s+1),MPI_RPREC,MPI_SUM,COMM_TINKER,ierr)
!$acc end host_data
         end if

         ! Save Old W_{i-1}
!$acc parallel loop async(rec_queue) default(present)
         do i = 1,nrhs*s**2
            MatW_sr(i,1,1) = MatWs_r(i,1,1)
         end do

         ! Find alpha scalar
         ! Solve W_i *(a_i) = g_i
         do j = 1,nrhs
#ifdef _OPENACC
!$acc host_data use_device(MatWs_r,ipiv,Vas_r)
           call cuGESVm(s,1,MatWs_r(1,1,j),s,ipiv,Vas_r(1,j)
     &                 ,s,rec_stream)
!$acc end host_data
#else
           call M_gesvm(s,1,MatWs_r(1:,1,j),s,ipiv,Vas_r(1:,j),s,info)
#endif
         end do
         call amovec(size(stepMatVr),stepMatVr,stepMatV,rec_queue)

#if 0
         if (rank.Eq.0) then
!$acc wait
!$acc update host(stepMatVr)
            if (s.eq.5) then
               write(*,'(A,137X,A)') 'matWs','alpha'
            else
               write(*,'(A,193X,A)') 'matWs','alpha'
            end if
 18         format(D14.6,$)
            do j = 1,s
            write(*,18) (matW_sr(j,i,1),i=1,s),(matW_sr(j,i,2),i=1,s)
     &                 ,(Vas_r(j,i),i=1,nrhs)
            write(*,*)
            end do
         end if
#endif

         ! mu = mu + P_i*a_i
!$acc parallel loop gang vector collapse(3) private(row)
!$acc&         async(rec_queue) present(Vas_r)
         do k = 1,npoleloc; do j = 1,nrhs; do i = 1,3
!$acc loop seq
            do its = 1,s
               row(its) = MPs(i,k,j,its)
            end do
            rtr = 0
            do its = 1,s
               rtr = rtr + row(its)*Vas_r(its,j)
            end do
            mu(i,j,k) = mu(i,j,k) + rtr
         end do; end do; end do
c
         call timer_exit ( timer_stepcg,quiet_timers )

         ! Update New residu
         ! Compute Matrix Vector Product
         ! -----------------------------

         ! get neighbor contribution of field
         if (use_polarshortreal) then
#ifdef _OPENACC
            if (dir_queue.ne.rec_queue) call start_dir_stream_cover
#endif
            if (.not.skpPcomm) then
            call commdirdirshort(nrhs,0,mu,reqendrec,reqendsend)
            call commdirdirshort(nrhs,1,mu,reqendrec,reqendsend)
            call commdirdirshort(nrhs,2,mu,reqendrec,reqendsend)
            end if
         else
            call commdirdirgpu(nrhs,0,mu,reqendrec,reqendsend)
            call commdirdirgpu(nrhs,1,mu,reqendrec,reqendsend)
            call commdirdirgpu(nrhs,2,mu,reqendrec,reqendsend)

            ! Begin Communication of mu for PME
            call commrecdirdipgpu(nrhs,0,murec,mu,buffermpimu1,
     &           buffermpimu2,req2rec,req2send)
            call commrecdirdipgpu(nrhs,1,murec,mu,buffermpimu1,
     &           buffermpimu2,req2rec,req2send)
            call commrecdirdipgpu(nrhs,2,murec,mu,buffermpimu1,
     &           buffermpimu2,req2rec,req2send)
         end if

         ! Check one iteration optimisation and exit solver
         if (optOneIt(mode)) goto 10

         ! Real space
         call timer_enter( timer_realdip )
         call matvec(nrhs,.true.,mu,h)
         call timer_exit ( timer_realdip,quiet_timers )

         ! Reciprocal space
         if (.not.use_polarshortreal) then
         call timer_enter( timer_recdip )
         call tmatxbrecipgpu(mu,murec,nrhs,dipfield,dipfieldbis)
         call timer_exit ( timer_recdip,quiet_timers )
         end if

         ! get missing real space contribution
         if (use_polarshortreal) then
            call commfieldshort(nrhs,h)
#ifdef _OPENACC
            if (dir_queue.ne.rec_queue) call end_dir_stream_cover
#endif
         else
            call commfieldgpu(nrhs,h)

         ! get missing reciprocal space contribution
         call commrecdirsolvgpu(nrhs,0,dipfieldbis,dipfield,buffermpi1,
     $        buffermpi2,reqrecdirrec,reqrecdirsend)
         call commrecdirsolvgpu(nrhs,1,dipfieldbis,dipfield,buffermpi1,
     $        buffermpi2,reqrecdirrec,reqrecdirsend)
         call commrecdirsolvgpu(nrhs,2,dipfieldbis,dipfield,buffermpi1,
     $        buffermpi2,reqrecdirrec,reqrecdirsend)
         end if

         call timer_enter( timer_stepcg )
         ggnew1 = 0; ggnew2 = 0;
!$acc parallel loop collapse(3) async(rec_queue)
!$acc&         default(present)
         do k=1,npoleloc; do j=1,nrhs; do i=1,3
            if (.not.use_polarshortreal)
     &      h(i,j,k) = h(i,j,k) + dipfield(i,j,k) - term*mu(i,j,k)
            rt       = (ef(i,j,k) - h(i,j,k))
            if (btest(j,0)) then
               ggnew1 = ggnew1 + diag(k)*rt**2
            else
               ggnew2 = ggnew2 + diag(k)*rt**2
            end if
            res(i,j,k) = rt
         end do; end do; end do

!$acc wait(rec_queue)
         ggnew(1) = ggnew1; ggnew(2)=ggnew2

c        if (nproc.gt.1.and.it.eq.1) call MPI_WAIT(req1,status,ierr)

c        if (it.eq.1.and.polprt.lt.2) then
c           call timer_exit( timer_stepcg,quiet_timers )
c           cycle ! Avoid Convergence check at first iteration
c        end if

         ! ----------  Compute norm & Convergence verification  ----------
         if (nproc.gt.1) then
            call MPI_IALLREDUCE(MPI_IN_PLACE,ggnew,nrhs,MPI_RPREC
     &          ,MPI_SUM,COMM_TINKER,req3,ierr)
            call MPI_WAIT(req3,status,ierr)
            ggnew1=ggnew(1); ggnew2=ggnew(2)
         end if

         resnrm = -0.0_re_p
         do k = 1,nrhs
            gnorm(k) = sqrt(ggnew(k)/real(3*npolar,r_p))
            resnrm   = max(resnrm,gnorm(k))
         end do

         if (polprt.ge.2.and.rank.eq.0) write(iout,1010)
     &      it, (gnorm(k), k=1,nrhs), (ggnew(k), k=1,nrhs)
c
         if ((it.eq.politer.and.resnrm.gt.poleps)
     &          .or.tinker_isnan_m(gnorm(1))
     &          .or.tinker_isnan_m(gnorm(2))) then
            ! Not converged Exit
!$acc wait
            if (rank.eq.0) write(0,1050) it,gnorm(1),gnorm(2)
            abort = .true.
         end if
c
         if (resnrm.lt.poleps) then
            if (polprt.ge.1.and.rank.eq.0) write(iout,1000) it,
     &         (gnorm(k), k = 1, nrhs)
            call timer_exit( timer_stepcg,quiet_timers )
            goto 10
         end if
         ! ---------------------------------------------------------------
c
         call timer_exit( timer_stepcg,quiet_timers )

      end do
 10   continue

      if (it.eq.1) then
         OneItInARow(mode) = OneItInARow(mode) + 1
         MoreInARow(mode)  = 0
      else
         OneItInARow(mode) = 0
         MoreInARow(mode)  = MoreInARow(mode) + 1
         optOneIt(mode)    = .false.
      end if

      if (MoreInARow(mode).eq.20) then
         MoreInARow(mode)  = 0
         if (use_polarshortreal) then
            sshort= min(cgstep_max,sshort+1)
         else
            slong = min(cgstep_max, slong+1)
         end if
         if (rank.eq.0.and.tinkerdebug.gt.0)
     &      print*,"Increase krylov subspace dimension",sshort,slong 
      end if

      if (OneItInARow(mode).eq.MaxIt1InARow-1) then
         optOneIt(mode)=.false.
      else if (OneItInARow(mode).eq.MaxIt1InARow) then
         optOneIt(mode)= .true.
         OneItInARow(mode) = 0
      end if

 
      ! Debug printing
      ! <<<<<<<<<<<<<<
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
      ! >>>>>>>>>>>>>>

 50   continue
!$acc end data
!$acc exit data detach(MRes_s,MQs,MTs,MPs) async(rec_queue)
      end subroutine
#endif
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
!$acc data present(gear,aspc,bpred,bpredp,bpreds,bpredps)
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
!$acc data create(b,bp,a,ap,c,cp,rc,rcp,array) async(def_queue)
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

!$acc parallel loop collapse(3) present(A,B) async
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
