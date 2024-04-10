c
c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c
c     "bspline_fill_site" finds B-spline coefficients and derivatives
c     for PME i-th atomic sites along the fractional coordinate axes
c
c
#define TINKER_CUF
#include "tinker_macro.h"
#include "tinker_cudart.h"
      module pmestuffcu
        use math   ,only: pi
        use pme    ,only: maxorder
        use tinheader
        use utilcu ,only: PME_BLOCK_DIM,PME_GRID_BDIM,PME_FPHI_DIM1
     &             , octahedron,use_virial
        use utilgpu,only: RED_BUFF_SIZE

        integer,parameter ::  level=4
        integer,parameter ::  level1=2
#if defined(BSORDER_VALUE) && BSORDER_VALUE>5
        integer,parameter :: bsorder=BSORDER_VALUE
#else
#   if BSORDER_VALUE<5
#      undef BSORDER_VALUE
#      define BSORDER_VALUE 5
#   endif
        integer,parameter :: bsorder=5
#endif
        ! PME pointer variable
        real(t_p),constant :: recip_c(3,3)
#if (defined(SINGLE)||defined(MIXED))
        real(t_p),constant :: pme_eps
#else
        real(t_p),parameter :: pme_eps=1d-8
#endif
        integer  ,device, pointer:: igrid_t(:,:)
     &           ,kstart1_t(:),kend1_t(:),jstart1_t(:)
     &           ,jend1_t(:),istart1_t(:),iend1_t(:)
     &           ,prec_send_t(:)
        real(t_p),device, pointer:: x_t(:),y_t(:),z_t(:)
     &           ,thetai1_t(:,:,:),thetai2_t(:,:,:),thetai3_t(:,:,:)
     &           ,qgrid2in_t(:,:,:,:,:),qgridin_t(:,:,:,:,:)

        contains

        subroutine copy_recip_cu
        use boxes, only:recip
        implicit none
        recip_c = recip
        !call print_recip_c<<<1,1>>>
        end subroutine

        subroutine setcu_pme_eps(pme_eps_)
        implicit none
        real(t_p),intent(in)::pme_eps_
#if (defined(SINGLE)||defined(MIXED))
        pme_eps = pme_eps_
#endif
        end subroutine

        logical function pme_check_cuda_bsorder(tinker_bsorder)
     &          result(test)
        implicit none
        integer,intent(in)::tinker_bsorder
        if (tinker_bsorder.eq.bsorder) then
           test=.true.
        else
           test=.false.
        end if
        end function

        attributes(global) subroutine print_recip_c
        implicit none
        print*,'  recip_on_device'
        print*,recip_c(1,1),recip_c(2,1),recip_c(3,1)
        print*,recip_c(1,2),recip_c(2,2),recip_c(3,2)
        print*,recip_c(1,3),recip_c(2,3),recip_c(3,3)
        end

c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine bsplgen  --  B-spline coefficients for an atom  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "bsplgen" gets B-spline depth 4 coefficients and derivatives for
c     a single PME atomic site along a particular direction
c
c
        attributes(device) subroutine ibsplgen (w,thetai,temp)
        implicit none
        integer i,j,k
        integer isite
        real(t_p) w,denom
        real(t_p) thetai(4,bsorder)
        real(t_p) temp(bsorder*bsorder)
c
c       set B-spline depth for partial charges or multipoles
c

c
c       initialization to get to 2nd order recursion
c
        temp(2+bsorder*1) = w
        temp(2+bsorder*0) = 1.0_ti_p - w
c
c       perform one pass to get to 3rd order recursion
c
        temp(3+bsorder*2) = 0.5_ti_p * w * temp(2+bsorder*1)
        temp(3+bsorder*1) = 0.5_ti_p * ((1.0_ti_p+w)*temp(2+bsorder*0)+
     &                                  (2.0_ti_p-w)*temp(2+bsorder*1))
        temp(3+bsorder*0) = 0.5_ti_p * (1.0_ti_p-w) * temp(2+bsorder*0)
c
c       compute standard B-spline recursion to desired order
c
        do i = 4, bsorder
           k = i - 1
           denom = 1.0_ti_p / real(k,t_p)
           temp(i+bsorder*(i-1)) = denom * w * temp(k+bsorder*(k-1))
           do j = 1, i-2
              temp(i+bsorder*(i-j-1)) = denom * 
     &                      ((w+real(j,t_p))*temp(k+bsorder*(i-j-1-1))
     &                    +(real(i-j,t_p)-w)*temp(k+bsorder*(i-j-1)))
           end do
           temp(i+bsorder*0) = denom * (1.0_ti_p-w) * temp(k+bsorder*0)
        end do
c
c       get coefficients for the B-spline first derivative
c
        k = bsorder - 1
        temp(k+bsorder*(bsorder-1)) = temp(k+bsorder*(bsorder-1-1))
        do i = bsorder-1, 2, -1
           temp(k+bsorder*(i-1)) = temp(k+bsorder*(i-1-1)) 
     &                           - temp(k+bsorder*(i-1))
        end do
        temp(k+bsorder*0) = -temp(k+bsorder*0)
        ! if (level.eq.4) then
c
c       get coefficients for the B-spline second derivative
c
        k = bsorder - 2
        temp(k+bsorder*(bsorder-2)) = temp(k+bsorder*(bsorder-3))
        do i = bsorder-2, 2, -1
           temp(k+bsorder*(i-1)) = temp(k+bsorder*(i-2))
     &                           - temp(k+bsorder*(i-1))
        end do
        temp(k+bsorder*0) = -temp(k+bsorder*0)
        temp(k+bsorder*(bsorder-1)) = temp(k+bsorder*(bsorder-2))
        do i = bsorder-1, 2, -1
           temp(k+bsorder*(i-1)) = temp(k+bsorder*(i-1-1))
     &                           - temp(k+bsorder*(i-1))
        end do
        temp(k+bsorder*0) = -temp(k+bsorder*0)
c
c       get coefficients for the B-spline third derivative
c
        k = bsorder - 3
        temp(k+bsorder*(bsorder-3)) = temp(k+bsorder*(bsorder-4))
        do i = bsorder-3, 2, -1
           temp(k+bsorder*(i-1)) = temp(k+bsorder*(i-2))
     &                           - temp(k+bsorder*(i-1))
        end do
        temp(k+bsorder*0) = -temp(k+bsorder*0)
        temp(k+bsorder*(bsorder-2)) = temp(k+bsorder*(bsorder-3))
        do i = bsorder-2, 2, -1
           temp(k+bsorder*(i-1)) = temp(k+bsorder*(i-2))
     &                           - temp(k+bsorder*(i-1))
        end do
        temp(k+bsorder*0) = -temp(k+bsorder*0)
        temp(k+bsorder*(bsorder-1)) = temp(k+bsorder*(bsorder-2))
        do i = bsorder-1, 2, -1
           temp(k+bsorder*(i-1)) = temp(k+bsorder*(i-2))
     &                           - temp(k+bsorder*(i-1))
        end do
        temp(k+bsorder*0) = -temp(k+bsorder*0)
        ! end if
c
c       copy coefficients from temporary to permanent storage
c
        do i = 1, bsorder
           do j = 1, level
              thetai(j,i) = temp(bsorder-j+1+bsorder*(i-1))
           end do
        end do
        end

c       "bsplgen" gets B-spline depth 2 coefficients and derivatives for
c       a single PME atomic site along a particular direction
c        ( Optimized for partial charge )
        attributes(device) subroutine ibsplgen_2 (w,thetai,temp)
        implicit none
        integer i,j,k
        integer isite
        real(t_p) w,denom
        real(t_p) thetai(level1*bsorder)
        real(t_p) temp(bsorder*bsorder)
c
c       set B-spline depth for partial charges or multipoles
c

c
c       initialization to get to 2nd order recursion
c
        temp(2+bsorder*1) = w
        temp(2+bsorder*0) = 1.0_ti_p - w
c
c       perform one pass to get to 3rd order recursion
c
        temp(3+bsorder*2) = 0.5_ti_p * w * temp(2+bsorder*1)
        temp(3+bsorder*1) = 0.5_ti_p * ((1.0_ti_p+w)*temp(2+bsorder*0)+
     &                                  (2.0_ti_p-w)*temp(2+bsorder*1))
        temp(3+bsorder*0) = 0.5_ti_p * (1.0_ti_p-w) * temp(2+bsorder*0)
c
c       compute standard B-spline recursion to desired order
c
        do i = 4, bsorder
           k = i - 1
           denom = 1.0_ti_p / real(k,t_p)
           temp(i+bsorder*(i-1)) = denom * w * temp(k+bsorder*(k-1))
           do j = 1, i-2
              temp(i+bsorder*(i-j-1)) = denom * 
     &                      ((w+real(j,t_p))*temp(k+bsorder*(i-j-1-1))
     &                    +(real(i-j,t_p)-w)*temp(k+bsorder*(i-j-1)))
           end do
           temp(i+bsorder*0) = denom * (1.0_ti_p-w) * temp(k+bsorder*0)
        end do
c
c       get coefficients for the B-spline first derivative
c
        k = bsorder - 1
        temp(k+bsorder*(bsorder-1)) = temp(k+bsorder*(bsorder-1-1))
        do i = bsorder-1, 2, -1
           temp(k+bsorder*(i-1)) = temp(k+bsorder*(i-1-1)) 
     &                           - temp(k+bsorder*(i-1))
        end do
        temp(k+bsorder*0) = -temp(k+bsorder*0)
c
c       copy coefficients from temporary to permanent storage
c
        do i = 1, bsorder; do j = 1, level1
           thetai(j+(i-1)*level1) = temp(bsorder-j+1+bsorder*(i-1))
        end do; end do
        end

c       "bsplgen" gets B-spline depth 2 coefficients and derivatives for
c       a single PME atomic site along a particular direction
c        ( Optimized for partial charge grid interpolation on the fly )
        attributes(device) subroutine ibsplgen_21 (w,thetai,temp)
        implicit none
        real(t_p),intent(in):: w
        real(t_p),intent(out):: thetai(bsorder)
        real(t_p),intent(out):: temp(bsorder*bsorder)
c
c       initialization to get to 2nd order recursion
c
        temp(2+bsorder*1) = w
        temp(2+bsorder*0) = 1.0_ti_p - w
c
c       perform one pass to get to 3rd order recursion
c
        temp(3+bsorder*2) = 0.5_ti_p * w * temp(2+bsorder*1)
        temp(3+bsorder*1) = 0.5_ti_p * ((1.0_ti_p+w)*temp(2+bsorder*0)+
     &                                  (2.0_ti_p-w)*temp(2+bsorder*1))
        temp(3+bsorder*0) = 0.5_ti_p * (1.0_ti_p-w) * temp(2+bsorder*0)
c
c       compute standard B-spline recursion to desired order
c
        block
        integer i,j,k
        real(t_p) denom
        do i = 4, bsorder
           k = i - 1
           denom = 1.0_ti_p / real(k,t_p)
           temp(i+bsorder*(i-1)) = denom * w * temp(k+bsorder*(k-1))
           do j = 1, i-2
              temp(i+bsorder*(i-j-1)) = denom * 
     &                      ((w+real(j,t_p))*temp(k+bsorder*(i-j-1-1))
     &                    +(real(i-j,t_p)-w)*temp(k+bsorder*(i-j-1)))
           end do
           temp(i+bsorder*0) = denom * (1.0_ti_p-w) * temp(k+bsorder*0)
        end do
        end block
c
c       get coefficients for the B-spline first derivative
c
        block; integer k,i
        k = bsorder - 1
        temp(k+bsorder*(bsorder-1)) = temp(k+bsorder*(bsorder-1-1))
        do i = bsorder-1, 2, -1
           temp(k+bsorder*(i-1)) = temp(k+bsorder*(i-1-1)) 
     &                           - temp(k+bsorder*(i-1))
        end do
        temp(k+bsorder*0) = -temp(k+bsorder*0)
        end block
        end

        attributes(device) subroutine get_thetai_21(thetai,temp)
        implicit none
        real(t_p),intent(out):: thetai(bsorder)
        real(t_p),intent(out):: temp(bsorder*bsorder)
c
c       copy coefficients from temporary to permanent storage
c
        block; integer i
        do i = 1, bsorder
           thetai(i) = temp(bsorder+bsorder*(i-1))
        end do
        !thetai(1) = temp(bsorder+bsorder*(1-1))
        !thetai(2) = temp(bsorder+bsorder*(2-1))
        !thetai(3) = temp(bsorder+bsorder*(3-1))
        !thetai(4) = temp(bsorder+bsorder*(4-1))
        !thetai(5) = temp(bsorder+bsorder*(5-1))
c#if BSORDER_VALUE == 6
c        thetai(6) = temp(bsorder+bsorder*(6-1))
c#endif
c#if BSORDER_VALUE == 7
c        thetai(7) = temp(bsorder+bsorder*(7-1))
c#endif
c#if BSORDER_VALUE == 8
c        thetai(8) = temp(bsorder+bsorder*(8-1))
c#endif
c#if BSORDER_VALUE > 8
c#   warning(" bsorder value is not take under consideration")
c#endif

        end block
        end
c
c
c       "grid_mpole_site" places the i-th fractional atomic multipole onto
c       the particle mesh Ewald grid
c       CUDA Fortran written kernel
c       Work Only in Sequential
c
c
        attributes(global) subroutine grid_mpole_sitecu_core_1p
     &            (polerecglob,ipole,igrid,thetai1,thetai2,thetai3
     &            ,kstat,ked,jstat,jed,istat,ied
     &            ,nfft1,nfft2,nfft3,npolerecloc
     &            ,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff
     &            ,order,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &            ,fmpvec,qgrid2loc,fir)
        implicit none
        integer istart,iend,jstart,jend,kstart,kend
        integer i,j,k,impi
        integer mk,mj,mi,m
        integer iproc,proc,rankloc
        integer ii,jj,kk
        integer iipole,iglob,isite,iatm
        integer offsetx,offsety
        integer offsetz
        real(t_p) term
        real(t_p) term0,term1,term2
        real(t_p) rstat
        integer igrid1,igrid2,igrid3

        integer  ,intent(in),value:: order,nlpts,grdoff,twonlpts_1
     &           ,twonlpts_12,nlptsit,kstat,ked,jstat,jed,istat,ied
     &           ,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &           ,nfft1,nfft2,nfft3,npolerecloc
        logical  ,value,intent(in):: fir
        integer  ,device,intent(in):: polerecglob(*),igrid(3,*)
     &           ,ipole(*)
        real(t_p),device,intent(in):: thetai1(4*order,*)
     &           ,thetai2(4*order,*),thetai3(4*order,*)
        real(t_p),device,dimension(10,*),intent(in) :: fmpvec
        real(t_p),device:: qgrid2loc(2,n1mpimax,n2mpimax,n3mpimax,
     &                      nrec_send+1)
        real(t_p),shared::theta1(4*bsorder),theta2(4*bsorder)
     &           ,theta3(4*bsorder)
        real(t_p),shared::  fmp(10)
        real(t_p) vut(6)

        do impi = blockIdx%x, npolerecloc, gridDim%x
           isite     = (polerecglob(impi))
           igrid1    = igrid(1,isite)
           igrid2    = igrid(2,isite)
           igrid3    = igrid(3,isite)
           offsetx   = 1 - ( igrid1 + grdoff - nlpts )
           offsety   = 1 - ( igrid2 + grdoff - nlpts )
           offsetz   = 1 - ( igrid3 + grdoff - nlpts )

           ! Load data into shared memory
           if (threadIdx%x<11) then
              fmp(threadIdx%x) = fmpvec(threadIdx%x,isite)
           else
              do i=threadIdx%x-10,4*order,blockDim%x-10
                 theta1(i) = thetai1(i,impi)
                 theta2(i) = thetai2(i,impi)
                 theta3(i) = thetai3(i,impi)
              end do
           end if
           if (blockDim%x>warpsize) call syncthreads
c
c       Three dimensional loop on the grid collapse by hand
c
           do kk = threadIdx%x-1, nlptsit, blockDim%x

              k  = igrid3 + grdoff + kk/twonlpts_12 - nlpts
              j  = igrid2 + grdoff
     &           + mod(kk/twonlpts_1,twonlpts_1) - nlpts
              i  = igrid1 + grdoff
     &           + mod(kk,twonlpts_1)-nlpts
              mk     = (k + offsetz -1)*4
              if (k .lt. 1) k = k + nfft3
              mj     = (j + offsety -1)*4
              if (j .lt. 1) j = j + nfft2
              mi     = (i + offsetx -1)*4
              if (i .lt. 1) i = i + nfft1

              term0 = fmp(01)*theta2(1+mj)*theta3(1+mk)
     &              + fmp(03)*theta2(2+mj)*theta3(1+mk)
     &              + fmp(04)*theta2(1+mj)*theta3(2+mk)
     &              + fmp(06)*theta2(3+mj)*theta3(1+mk)
     &              + fmp(07)*theta2(1+mj)*theta3(3+mk)
     &              + fmp(10)*theta2(2+mj)*theta3(2+mk)
              term1 = fmp(02)*theta2(1+mj)*theta3(1+mk)
     &              + fmp(08)*theta2(2+mj)*theta3(1+mk)
     &              + fmp(09)*theta2(1+mj)*theta3(2+mk)
              term2 = fmp(05)*theta2(1+mj)*theta3(1+mk)

              term  = term0*theta1(1+mi) + term1*theta1(2+mi)
     &              + term2*theta1(3+mi)

              rstat = AtomicAdd(
     &        qgrid2loc(1,i-istat+1,j-jstat+1,k-kstat+1,1)
     &        , term)

           end do
        end do
        end
c
c
c       "grid_pchg_site" places the i-th fractional atomic charge onto
c       the particle mesh Ewald grid
c       CUDA Fortran written kernel
c       Work Only in Sequential
c
c
        attributes(global) subroutine grid_put_site_kcu_o
     &            (kind_id,atom_id,igrid,pchg,thetai1,thetai2,thetai3
     &            ,x_sp,y_sp,z_sp
     &            ,kstat,ked,jstat,jed,istat,ied
     &            ,nfft1,nfft2,nfft3,nionrecloc
     &            ,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff
     &            ,order,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &            ,qgrid2loc,fir)
        implicit none
        integer  ,intent(in),value:: order,nlpts,grdoff,twonlpts_1
     &           ,twonlpts_12,nlptsit,kstat,ked,jstat,jed,istat,ied
     &           ,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &           ,nfft1,nfft2,nfft3,nionrecloc
        logical  ,value,intent(in):: fir
        integer  ,device,intent(in):: kind_id(*),igrid(3,*)
     &           ,atom_id(*)
        real(t_p),device,intent(in):: pchg(*),thetai1(2*order,*)
     &           ,thetai2(2*order,*),thetai3(2*order,*)
     &           ,x_sp(*),y_sp(*),z_sp(*)
        real(t_p),device:: qgrid2loc(2,n1mpimax,n2mpimax,n3mpimax,
     &                      nrec_send+1)

        real(t_p),shared::theta1(level1*bsorder),theta2(level1*bsorder)
     &           ,theta3(level1*bsorder)
        integer   impi,kk,i,isite,igrid1,igrid2,igrid3
     &           ,offsetx,offsety,offsetz
        real(t_p) q

        do impi = blockIdx%x, nionrecloc, gridDim%x
           isite  = atom_id(impi)
           igrid1 = igrid(1,isite)
           igrid2 = igrid(2,isite)
           igrid3 = igrid(3,isite)
           ! Load data into shared memory
           do i=threadIdx%x,2*order,blockDim%x
              theta1(i) = thetai1(i,impi)
              theta2(i) = thetai2(i,impi)
              theta3(i) = thetai3(i,impi)
           end do
           if (blockDim%x>warpsize) call syncthreads
           q         = pchg(kind_id(impi))
           offsetx   = 1 - ( igrid1 + grdoff - nlpts )
           offsety   = 1 - ( igrid2 + grdoff - nlpts )
           offsetz   = 1 - ( igrid3 + grdoff - nlpts )
c
c       Three dimensional loop on the grid collapse by hand
c
           do kk = threadIdx%x-1, nlptsit, blockDim%x; block
           integer k,j,i,mk,mj,mi
           real(t_p) term,rstat

              k  = igrid3 + grdoff + kk/twonlpts_12 - nlpts
              j  = igrid2 + grdoff
     &           + mod(kk/twonlpts_1,twonlpts_1) - nlpts
              i  = igrid1 + grdoff
     &           + mod(kk,twonlpts_1)-nlpts
              mk     = (k + offsetz -1)*2
              if (k .lt. 1) k = k + nfft3
              mj     = (j + offsety -1)*2
              if (j .lt. 1) j = j + nfft2
              mi     = (i + offsetx -1)*2
              if (i .lt. 1) i = i + nfft1

              term  = q*theta2(1+mj)*theta1(1+mi)*theta3(1+mk)

              rstat = AtomicAdd(
     &        qgrid2loc(1,i-istat+1,j-jstat+1,k-kstat+1,1),term )

           end block; end do
        end do
        end

        attributes(global) subroutine grid_put_site_kcu1
     &            (kind_id,atom_id,igrid,attr,thetai1,thetai2,thetai3
     &            ,x_sp,y_sp,z_sp
     &            ,kstat,ked,jstat,jed,istat,ied
     &            ,nfft1,nfft2,nfft3,nionrecloc
     &            ,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff
     &            ,order,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &            ,qgrid2loc,fir)
        implicit none
        integer  ,intent(in),value:: order,nlpts,grdoff,twonlpts_1
     &           ,twonlpts_12,nlptsit,kstat,ked,jstat,jed,istat,ied
     &           ,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &           ,nfft1,nfft2,nfft3,nionrecloc
        logical  ,value,intent(in):: fir
        integer  ,device,intent(in):: kind_id(*),igrid(3,*)
     &           ,atom_id(*)
        real(t_p),device,intent(in):: attr(*),thetai1(2*order,*)
     &           ,thetai2(2*order,*),thetai3(2*order,*)
     &           ,x_sp(*),y_sp(*),z_sp(*)
        real(t_p),device:: qgrid2loc(2,n1mpimax,n2mpimax,n3mpimax,
     &                      nrec_send+1)

        real(t_p),shared::temp(bsorder*bsorder*PME_GRID_BDIM)
        real(t_p)::theta1(bsorder),theta2(bsorder)
     &            ,theta3(bsorder)

        integer   impi,kk,i,isite,igrid1,igrid2,igrid3
     &           ,offsetx,offsety,offsetz,offset
        real(t_p) q

        offset    = (threadIdx%x-1)*bsorder**2+1

        do impi = (blockIdx%x-1)*blockDim%x + threadIdx%x, nionrecloc
     &          , blockDim%x*gridDim%x
           q      = attr(kind_id(impi))
c
c          get the b-spline coefficients for the i-th atomic site
c          it faster to recompute theta*
c
           block
           integer ifr
           real(t_p) xi,yi,zi,w,fr
           xi     = x_sp(impi)
           yi     = y_sp(impi)
           zi     = z_sp(impi)

           w      = xi*recip_c(1,1) + yi*recip_c(2,1) + zi*recip_c(3,1)
           fr     = nfft1 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           igrid1 = ifr - bsorder
           w      = fr - real(ifr,t_p)
           call ibsplgen_21(w,theta1,temp(offset))
           ! FIXME Why is this call using so much registers
           call get_thetai_21(theta1,temp(offset))

           w      = xi*recip_c(1,2) + yi*recip_c(2,2) + zi*recip_c(3,2)
           fr     = nfft2 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           igrid2 = ifr - bsorder
           w      = fr - real(ifr,t_p)
           call ibsplgen_21(w,theta2,temp(offset))
           call get_thetai_21(theta2,temp(offset))

           w      = xi*recip_c(1,3) + yi*recip_c(2,3) + zi*recip_c(3,3)
           fr     = nfft3 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           igrid3 = ifr - bsorder
           w      = fr - real(ifr,t_p)
           call ibsplgen_21(w,theta3,temp(offset))
           call get_thetai_21(theta3,temp(offset))
           end block

           offsetx = 1 - ( igrid1 + grdoff - nlpts )
           offsety = 1 - ( igrid2 + grdoff - nlpts )
           offsetz = 1 - ( igrid3 + grdoff - nlpts )
c
c          Three dimensional loop on the grid collapse by hand
c
           block; integer k1,k
           do k1 = 0, bsorder-1;
              k  =  1 - offsetz + k1
              if (k.lt.1) k = k + nfft3

              block; integer j1,j
              do j1 = 0, bsorder-1
                 j  =  1 - offsety + j1
                 if (j.lt.1) j = j + nfft2

                 block; integer i1,i
                 real(t_p) term,rstat
                 do i1 = 0, bsorder-1
                    i  =  1 - offsetx + i1 
                    if (i.lt.1) i = i + nfft1

                 term  = q*theta2(1+j1)*theta1(1+i1)*theta3(1+k1)
                    rstat = AtomicAdd(
     &              qgrid2loc(1,i-istat+1,j-jstat+1,k-kstat+1,1),term)
                 end do; end block
              end do; end block
           end do; end block 
        end do
        end
c
c
c       "grid_uind_site" places the fractional induced dipoles onto the
c       particle mesh Ewald grid
c       CUDA Fortran written kernel
c       Work Only in Sequential
c
c
        attributes(global) subroutine grid_uind_sitecu_core_1p
     &            (polerecglob,ipole,igrid,thetai1,thetai2,thetai3
     &            ,kstat,ked,jstat,jed,istat,ied
     &            ,nfft1,nfft2,nfft3,npolerecloc
     &            ,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff
     &            ,order,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &            ,fuindvec,fuinpvec,qgrid2loc,fir)
        implicit none
        integer istart,iend,jstart,jend,kstart,kend
        integer i,j,k,impi
        integer mk,mj,mi
        integer iproc,proc,rankloc
        integer ii,jj,kk
        integer iipole,iglob,isite,iatm
        integer offsetx,offsety
        integer offsetz
        real(t_p) v0,u0,t0
        real(t_p) v1,u1,t1
        real(t_p) term
        real(t_p) term01,term11
        real(t_p) term02,term12
        real(t_p) rstat
        integer igrid1,igrid2,igrid3

        integer  ,intent(in),value:: order,nlpts,grdoff,twonlpts_1
     &           ,twonlpts_12,nlptsit,kstat,ked,jstat,jed,istat,ied
     &           ,n1mpimax,n2mpimax,n3mpimax,nrec_send
     &           ,nfft1,nfft2,nfft3,npolerecloc
        logical  ,value,intent(in):: fir
        integer  ,device,intent(in):: polerecglob(*),igrid(3,*)
     &           ,ipole(*)
        real(t_p),device,intent(in):: thetai1(4*order,*)
     &           ,thetai2(4*order,*),thetai3(4*order,*)
        real(t_p),device,dimension(3,*),intent(in) :: fuindvec,fuinpvec
        real(t_p),device:: qgrid2loc(2,n1mpimax,n2mpimax,n3mpimax,
     &                      nrec_send+1)
        real(t_p),shared::theta1(4*bsorder),theta2(4*bsorder)
     &           ,theta3(4*bsorder)
        real(t_p),shared::  fuinp(3),fuind(3)

        do impi = blockIdx%x, npolerecloc, gridDim%x
           isite     = (polerecglob(impi))
           offsetx   = 1 - (igrid(1,isite) + grdoff - nlpts)
           offsety   = 1 - (igrid(2,isite) + grdoff - nlpts)
           offsetz   = 1 - (igrid(3,isite) + grdoff - nlpts)

           ! Load data into shared memory
           if (threadIdx%x<4) then
              fuind(threadIdx%x) = fuindvec(threadIdx%x,isite)
              fuinp(threadIdx%x) = fuinpvec(threadIdx%x,isite)
           else
              do i=threadIdx%x-3,4*order,blockDim%x-3
                 theta1(i) = thetai1(i,impi)
                 theta2(i) = thetai2(i,impi)
                 theta3(i) = thetai3(i,impi)
              end do
           end if
           if (blockDim%x>warpsize) call syncthreads
c
c       Three dimensional loop on the grid collapse by hand
c
           do kk = threadIdx%x-1, nlptsit, blockDim%x

              k  = igrid(3,isite) + grdoff + kk/twonlpts_12 - nlpts
              j  = igrid(2,isite) + grdoff 
     &           + mod(kk/twonlpts_1,twonlpts_1) - nlpts
              i  = igrid(1,isite) + grdoff 
     &           + mod(kk,twonlpts_1)-nlpts
              mk     = k + offsetz -1
              if (k .lt. 1) k = k + nfft3
              mj     = j + offsety -1
              if (j .lt. 1) j = j + nfft2
              mi     = i + offsetx -1
              if (i .lt. 1) i = i + nfft1

              term01 = fuind(2)*theta2(2+mj*4)*theta3(1+mk*4)
     &               + fuind(3)*theta2(1+mj*4)*theta3(2+mk*4)
              term11 = fuind(1)*theta2(1+mj*4)*theta3(1+mk*4)
              term02 = fuinp(2)*theta2(2+mj*4)*theta3(1+mk*4)
     &               + fuinp(3)*theta2(1+mj*4)*theta3(2+mk*4)
              term12 = fuinp(1)*theta2(1+mj*4)*theta3(1+mk*4)
              term01 = term01*theta1(1+mi*4)
     &               + term11*theta1(2+mi*4)
              term02 = term02*theta1(1+mi*4)
     &               + term12*theta1(2+mi*4)
c
c             if (((k.ge.kstat).and.(k.le.ked)).and.
c    &            ((j.ge.jstat).and.(j.le.jed)).and.
c    &            ((i.ge.istat).and.(i.le.ied)))    then
                 rstat = AtomicAdd(
     &           qgrid2loc(1,i-istat+1,j-jstat+1,k-kstat+1,1)
     &         , term01 )
                 rstat = AtomicAdd(
     &           qgrid2loc(2,i-istat+1,j-jstat+1,k-kstat+1,1)
     &         , term02 )
c                cycle
c             end if

c             do iproc = 1, nrec_send
c                proc   = prec_send_t(iproc)
c                kstart = kstart1_t  (proc+1)
c                kend   = kend1_t    (proc+1)
c                jstart = jstart1_t  (proc+1)
c                jend   = jend1_t    (proc+1)
c                istart = istart1_t  (proc+1)
c                iend   = iend1_t    (proc+1)
c                if (((k.ge.kstart).and.(k.le.kend)).and.
c    &               ((j.ge.jstart).and.(j.le.jend)).and.
c    &               ((i.ge.istart).and.(i.le.iend))) then
c                rstat = AtomicAdd(
c    &             qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
c    &           , term01 )
c                rstat = AtomicAdd(
c    &             qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)
c    &           , term02 )
c                   exit
c                end if
c             end do
           end do
        end do
        end
c
c     "fphi_uind_sitegpu2" extracts the induced dipole potential at the i-th site from
c     the particle mesh Ewald grid second version
c
c
        attributes(global)  subroutine fphi_mpole_core
     &            (kstat,ked,jstat,jed,istat,ied,nfft1,nfft2,nfft3
     &            ,npolerecloc,nlocrec,order,nrec_send,nproc
     &            ,polerecglob,ipole,igrid
     &            ,fphirec)
        implicit none
        integer  ,intent(in),value:: npolerecloc,nlocrec,order
     &           ,nfft1,nfft2,nfft3,nrec_send,nproc
        integer  ,intent(in),value:: kstat,ked,jstat,jed,istat,ied
        integer  ,intent(in),device::polerecglob(*),ipole(*),igrid(3,*)
        real(t_p),device::fphirec(20,npolerecloc)

        integer istart,iend,jstart,jend,kstart,kend
        integer i,j,k,impi
        integer iproc,proc
        integer isite,iatm,iipole
        integer i0,j0,k0
        integer it1,it2,it3
        integer igrd0,jgrd0,kgrd0
        real(t_p) v0,v1,v2,v3
        real(t_p) u0,u1,u2,u3
        real(t_p) t0,t1,t2,t3,tq
        real(t_p) tu00,tu10,tu01,tu20,tu11
        real(t_p) tu02,tu21,tu12,tu30,tu03
        real(t_p) tuv000,tuv100,tuv010,tuv001
        real(t_p) tuv200,tuv020,tuv002,tuv110
        real(t_p) tuv101,tuv011,tuv300,tuv030
        real(t_p) tuv003,tuv210,tuv201,tuv120
        real(t_p) tuv021,tuv102,tuv012,tuv111
        real(t_p) xi,yi,zi
        real(t_p) w,fr
        integer ifr
        real(t_p),dimension(4,bsorder):: theta1,theta2,theta3
        real(t_p),shared::temp(bsorder*bsorder*PME_BLOCK_DIM)
        real(t_p),parameter::eps=1.0d-4

        do impi = threadIdx%x+(blockIdx%x-1)*blockDim%x,npolerecloc,
     &            blockDim%x*gridDim%x
           iipole = polerecglob(impi)
           iatm   = ipole(iipole)
#if 1
c
c       get the b-spline coefficients for the i-th atomic site
c       it faster to recompute theta*
c
           xi     = x_t(iatm)
           yi     = y_t(iatm)
           zi     = z_t(iatm)
           w      = xi*recip_c(1,1) + yi*recip_c(2,1) + zi*recip_c(3,1)
           fr     = nfft1 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           igrd0  = ifr - bsorder
           call ibsplgen (w,theta1,temp((threadIdx%x-1)*bsorder**2+1))
           w      = xi*recip_c(1,2) + yi*recip_c(2,2) + zi*recip_c(3,2)
           fr     = nfft2 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           jgrd0  = ifr - bsorder
           call ibsplgen (w,theta2,temp((threadIdx%x-1)*bsorder**2+1))
           w      = xi*recip_c(1,3) + yi*recip_c(2,3) + zi*recip_c(3,3)
           fr     = nfft3 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           kgrd0  = ifr - bsorder
           call ibsplgen (w,theta3,temp((threadIdx%x-1)*bsorder**2+1))
#else
           igrd0  = igrid(1,iatm)
           jgrd0  = igrid(2,iatm)
           kgrd0  = igrid(3,iatm)
#endif
           tuv000 = 0.0_ti_p
           tuv001 = 0.0_ti_p
           tuv010 = 0.0_ti_p
           tuv100 = 0.0_ti_p
           tuv200 = 0.0_ti_p
           tuv020 = 0.0_ti_p
           tuv002 = 0.0_ti_p
           tuv110 = 0.0_ti_p
           tuv101 = 0.0_ti_p
           tuv011 = 0.0_ti_p
           tuv300 = 0.0_ti_p
           tuv030 = 0.0_ti_p
           tuv003 = 0.0_ti_p
           tuv210 = 0.0_ti_p
           tuv201 = 0.0_ti_p
           tuv120 = 0.0_ti_p
           tuv021 = 0.0_ti_p
           tuv102 = 0.0_ti_p
           tuv012 = 0.0_ti_p
           tuv111 = 0.0_ti_p
           k0     = kgrd0
           do it3 = 1, bsorder
              k0   = k0 + 1
c             k   = k0 + 1 + (nfft3-isign(nfft3,k0))/2
              k   = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
#if 1
              v0   = theta3(1,it3)
              v1   = theta3(2,it3)
              v2   = theta3(3,it3)
              v3   = theta3(4,it3)
#else
              v0   = thetai3_t(1,it3,impi)
              v1   = thetai3_t(2,it3,impi)
              v2   = thetai3_t(3,it3,impi)
              v3   = thetai3_t(4,it3,impi)
#endif
              tu00 = 0.0_ti_p
              tu10 = 0.0_ti_p
              tu01 = 0.0_ti_p
              tu20 = 0.0_ti_p
              tu11 = 0.0_ti_p
              tu02 = 0.0_ti_p
              tu30 = 0.0_ti_p
              tu21 = 0.0_ti_p
              tu12 = 0.0_ti_p
              tu03 = 0.0_ti_p
              j0   = jgrd0
              do it2 = 1, bsorder
                 j0 = j0 + 1
c                j  = j0 + 1 + (nfft2-isign(nfft2,j0))/2
                 j  = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
#if 1
                 u0 = theta2(1,it2)
                 u1 = theta2(2,it2)
                 u2 = theta2(3,it2)
                 u3 = theta2(4,it2)
#else
                 u0 = thetai2_t(1,it2,impi)
                 u1 = thetai2_t(2,it2,impi)
                 u2 = thetai2_t(3,it2,impi)
                 u3 = thetai2_t(4,it2,impi)
#endif
                 t0 = 0.0_ti_p
                 t1 = 0.0_ti_p
                 t2 = 0.0_ti_p
                 t3 = 0.0_ti_p
                 i0 = igrd0
                 do it1 = 1, bsorder
                    i0 = i0 + 1
c                   i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                    i = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)

                    tq     = 0.0_ti_p
                    if (((k.ge.kstat).and.(k.le.ked)).and.
     &                  ((j.ge.jstat).and.(j.le.jed)).and.
     &                  ((i.ge.istat).and.(i.le.ied))) then
                    tq = qgridin_t(1,i-istat+1,j-jstat+1,
     &                               k-kstat+1,1)
                      goto 10
                    end if
                    do iproc = 1, nrec_send
                      proc   = prec_send_t(iproc)
                      kstart = kstart1_t(proc+1)
                      kend   = kend1_t  (proc+1)
                      jstart = jstart1_t(proc+1)
                      jend   = jend1_t  (proc+1)
                      istart = istart1_t(proc+1)
                      iend   = iend1_t  (proc+1)
                      if (((k.ge.kstart).and.(k.le.kend)).and.
     $                    ((j.ge.jstart).and.(j.le.jend)).and.
     $                    ((i.ge.istart).and.(i.le.iend))) then
                      tq  = qgridin_t(1,i-istart+1,j-jstart+1
     $                      ,k-kstart+1,iproc+1)
                        goto 10
                      end if
                    end do
                    cycle
c       
 10                 continue
#if 1
                    t0 = t0 + tq*theta1(1,it1)
                    t1 = t1 + tq*theta1(2,it1)
                    t2 = t2 + tq*theta1(3,it1)
                    t3 = t3 + tq*theta1(4,it1)
#else
                    t0 = t0 + tq*thetai1_t(1,it1,impi)
                    t1 = t1 + tq*thetai1_t(2,it1,impi)
                    t2 = t2 + tq*thetai1_t(3,it1,impi)
                    t3 = t3 + tq*thetai1_t(4,it1,impi)
#endif
                 end do
                 tu00 = tu00 + t0*u0
                 tu10 = tu10 + t1*u0
                 tu01 = tu01 + t0*u1
                 tu20 = tu20 + t2*u0
                 tu11 = tu11 + t1*u1
                 tu02 = tu02 + t0*u2
                 tu30 = tu30 + t3*u0
                 tu21 = tu21 + t2*u1
                 tu12 = tu12 + t1*u2
                 tu03 = tu03 + t0*u3
              end do
              tuv000 = tuv000 + tu00*v0
              tuv100 = tuv100 + tu10*v0
              tuv010 = tuv010 + tu01*v0
              tuv001 = tuv001 + tu00*v1
              tuv200 = tuv200 + tu20*v0
              tuv020 = tuv020 + tu02*v0
              tuv002 = tuv002 + tu00*v2
              tuv110 = tuv110 + tu11*v0
              tuv101 = tuv101 + tu10*v1
              tuv011 = tuv011 + tu01*v1
              tuv300 = tuv300 + tu30*v0
              tuv030 = tuv030 + tu03*v0
              tuv003 = tuv003 + tu00*v3
              tuv210 = tuv210 + tu21*v0
              tuv201 = tuv201 + tu20*v1
              tuv120 = tuv120 + tu12*v0
              tuv021 = tuv021 + tu02*v1
              tuv102 = tuv102 + tu10*v2
              tuv012 = tuv012 + tu01*v2
              tuv111 = tuv111 + tu11*v1
           end do
           fphirec(1,impi) = tuv000
           fphirec(2,impi) = tuv100
           fphirec(3,impi) = tuv010
           fphirec(4,impi) = tuv001
           fphirec(5,impi) = tuv200
           fphirec(6,impi) = tuv020
           fphirec(7,impi) = tuv002
           fphirec(8,impi) = tuv110
           fphirec(9,impi) = tuv101
           fphirec(10,impi) = tuv011
           fphirec(11,impi) = tuv300
           fphirec(12,impi) = tuv030
           fphirec(13,impi) = tuv003
           fphirec(14,impi) = tuv210
           fphirec(15,impi) = tuv201
           fphirec(16,impi) = tuv120
           fphirec(17,impi) = tuv021
           fphirec(18,impi) = tuv102
           fphirec(19,impi) = tuv012
           fphirec(20,impi) = tuv111
        end do
        end


        attributes(global) subroutine fphi_uind_sitecu2_core
     &             (kstat,ked,jstat,jed,istat,ied
     &             ,npolerecloc,nlocrec,order,n1mpimax,n2mpimax,n3mpimax
     &             ,nrec_send,nproc,prec_send,psize,nfft1,nfft2,nfft3
     &             ,kstart1,kend1,jstart1,jend1,istart1,iend1
     &             ,qgrid_in,polerecglob,ipole,igrid,x,y,z,recip
     &             ,thetai1,thetai2,thetai3
     &             ,fdip_phi1,fdip_phi2,fir)
        implicit none

        integer  ,intent(in),value:: npolerecloc,nlocrec,order
     &           ,n1mpimax,n2mpimax,n3mpimax,nrec_send,psize
     &           ,nfft1,nfft2,nfft3,nproc
        integer  ,intent(in),value:: kstat,ked,jstat,jed,istat,ied
        integer  ,intent(in),device::polerecglob(*),ipole(*),igrid(3,*)
     &           ,kstart1(psize),kend1(psize),jstart1(psize)
     &           ,jend1(psize),istart1(psize),iend1(psize)
     &           ,prec_send(nproc)
        real(t_p),intent(in),device::x(*),y(*),z(*),recip(3,3)
     &           ,thetai1(4,order,nlocrec),thetai2(4,order,nlocrec)
     &           ,thetai3(4,order,nlocrec)
     &           ,qgrid_in(2,n1mpimax,n2mpimax,n3mpimax,nrec_send)
        real(t_p),device,dimension(10,npolerecloc)::fdip_phi1,fdip_phi2
        logical,value,intent(in):: fir

        integer istart,iend,jstart,jend,kstart,kend
        integer i,j,k,impi,ifr
        integer iproc,proc
        integer isite,iatm,iipole
        integer i0,j0,k0
        integer it1,it2,it3
        integer igrd0,jgrd0,kgrd0
        integer ind
        logical notfind
        real(t_p) v0,v1,v2,v3
        real(t_p) u0,u1,u2,u3
        real(t_p) t0,t1,t2,t3
        real(t_p) t0_1,t0_2,t1_1,t1_2
        real(t_p) t2_1,t2_2,tq_1,tq_2
        real(t_p) tu00,tu10,tu01,tu20,tu11
        real(t_p) tu02,tu30,tu21,tu12,tu03
        real(t_p) tu00_1,tu01_1,tu10_1
        real(t_p) tu00_2,tu01_2,tu10_2
        real(t_p) tu20_1,tu11_1,tu02_1
        real(t_p) tu20_2,tu11_2,tu02_2
        real(t_p) tuv100_1,tuv010_1,tuv001_1
        real(t_p) tuv100_2,tuv010_2,tuv001_2
        real(t_p) tuv200_1,tuv020_1,tuv002_1
        real(t_p) tuv110_1,tuv101_1,tuv011_1
        real(t_p) tuv200_2,tuv020_2,tuv002_2
        real(t_p) tuv110_2,tuv101_2,tuv011_2
        real(t_p) tuv000,tuv100,tuv010,tuv001
        real(t_p) tuv200,tuv020,tuv002,tuv110
        real(t_p) tuv101,tuv011,tuv300,tuv030
        real(t_p) tuv003,tuv210,tuv201,tuv120
        real(t_p) tuv021,tuv102,tuv012,tuv111
        real(t_p) tuv1(9),tuv2(9)
        real(t_p) xi,yi,zi
        real(t_p) w,fr
        real(t_p),dimension(4,bsorder):: theta1,theta2,theta3
        real(t_p),shared::temp(bsorder*bsorder*PME_BLOCK_DIM)

c
        do impi = threadIdx%x+(blockIdx%x-1)*blockDim%x,npolerecloc,
     &            blockDim%x*gridDim%x
           iipole = polerecglob(impi)
           iatm   = ipole(iipole)
#if 1
c
c       get the b-spline coefficients for the i-th atomic site
c       it faster to recompute theta*
c
           xi     = x_t(iatm)
           yi     = y_t(iatm)
           zi     = z_t(iatm)
           w      = xi*recip_c(1,1) + yi*recip_c(2,1) + zi*recip_c(3,1)
           fr     = nfft1 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           igrd0  = ifr - bsorder
           call ibsplgen (w,theta1,temp((threadIdx%x-1)*bsorder**2+1))
           w      = xi*recip_c(1,2) + yi*recip_c(2,2) + zi*recip_c(3,2)
           fr     = nfft2 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           jgrd0  = ifr - bsorder
           call ibsplgen (w,theta2,temp((threadIdx%x-1)*bsorder**2+1))
           w      = xi*recip_c(1,3) + yi*recip_c(2,3) + zi*recip_c(3,3)
           fr     = nfft3 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           kgrd0  = ifr - bsorder
           call ibsplgen (w,theta3,temp((threadIdx%x-1)*bsorder**2+1))
#else
           igrd0  = igrid(1,iatm)
           jgrd0  = igrid(2,iatm)
           kgrd0  = igrid(3,iatm)
#endif
              tuv1(01) = 0.0_ti_p; tuv2(01) = 0.0_ti_p
              tuv1(02) = 0.0_ti_p; tuv2(02) = 0.0_ti_p
              tuv1(03) = 0.0_ti_p; tuv2(03) = 0.0_ti_p
              tuv1(04) = 0.0_ti_p; tuv2(04) = 0.0_ti_p
              tuv1(05) = 0.0_ti_p; tuv2(05) = 0.0_ti_p
              tuv1(06) = 0.0_ti_p; tuv2(06) = 0.0_ti_p
              tuv1(07) = 0.0_ti_p; tuv2(07) = 0.0_ti_p
              tuv1(08) = 0.0_ti_p; tuv2(08) = 0.0_ti_p
              tuv1(09) = 0.0_ti_p; tuv2(09) = 0.0_ti_p
c          k0       = kgrd0
           do it3 = 1, bsorder
c             k0     = k0 + 1
c             k      = k0 + 1 + (nfft3-isign(nfft3,kgrd0+it3))/2
              k      = kgrd0 +it3 +1 + 
     &                 ishft(nfft3-isign(nfft3,kgrd0+it3),-1)
#if 1
              v0     = theta3(1,it3)
              v1     = theta3(2,it3)
              v2     = theta3(3,it3)
#else
              v0     = thetai3_t(1,it3,impi)
              v1     = thetai3_t(2,it3,impi)
              v2     = thetai3_t(3,it3,impi)
#endif  
              tu00_1 = 0.0_ti_p
              tu01_1 = 0.0_ti_p
              tu10_1 = 0.0_ti_p
              tu20_1 = 0.0_ti_p
              tu11_1 = 0.0_ti_p
              tu02_1 = 0.0_ti_p
              tu00_2 = 0.0_ti_p
              tu01_2 = 0.0_ti_p
              tu10_2 = 0.0_ti_p
              tu20_2 = 0.0_ti_p
              tu11_2 = 0.0_ti_p
              tu02_2 = 0.0_ti_p
c             j0 = jgrd0
              do it2 = 1, bsorder
c                j0 = j0 + 1
c                j    = jgrd0 +it2 + 1 + (nfft2-isign(nfft2,jgrd0+it2))/2
                 j    = jgrd0 +it2 + 1 +
     &                        ishft(nfft2-isign(nfft2,jgrd0+it2),-1)
#if 1
                 u0   = theta2(1,it2)
                 u1   = theta2(2,it2)
                 u2   = theta2(3,it2)
#else
                 u0   = thetai2_t(1,it2,impi)
                 u1   = thetai2_t(2,it2,impi)
                 u2   = thetai2_t(3,it2,impi)
#endif  
                 t0_1 = 0.0_ti_p
                 t1_1 = 0.0_ti_p
                 t2_1 = 0.0_ti_p
                 t0_2 = 0.0_ti_p
                 t1_2 = 0.0_ti_p
                 t2_2 = 0.0_ti_p
                 i0   = igrd0
                 do it1 = 1, bsorder
c                   i0  = i0 + 1
c                   i    = igrd0 +it1 +1 + (nfft1-isign(nfft1,igrd0 +it1))/2
                    i    = igrd0 +it1 +1 +
     &                           ishft(nfft1-isign(nfft1,igrd0+it1),-1)
c
                    tq_1 = 0_ti_p
                    tq_2 = 0_ti_p
                    if (((k.ge.kstat).and.(k.le.ked)).and.
     &                  ((j.ge.jstat).and.(j.le.jed)).and.
     &                  ((i.ge.istat).and.(i.le.ied))) then
c                 tq_1 = qgrid2in_p(ind)
c                 tq_2 = qgrid2in_p(ind+1)
                    tq_1 = qgrid_in(1,i-istat+1,j-jstat+1,k-kstat+1,1)
                    tq_2 = qgrid_in(2,i-istat+1,j-jstat+1,k-kstat+1,1)
                    else 
                      iproc = 1
                      notfind  = .true.
                      do while (notfind.and.(iproc.le.nrec_send))
                        proc   = prec_send_t(iproc)
                        kstart = kstart1_t (proc+1)
                        kend   = kend1_t   (proc+1)
                        jstart = jstart1_t (proc+1)
                        jend   = jend1_t   (proc+1)
                        istart = istart1_t (proc+1)
                        iend   = iend1_t   (proc+1)
                      if (((k.ge.kstart).and.(k.le.kend)).and.
     &                    ((j.ge.jstart).and.(j.le.jend)).and.
     &                    ((i.ge.istart).and.(i.le.iend))) then
                  tq_1 = qgrid_in(1,i-istart+1,j-jstart+1,k-kstart+1,
     &                       iproc+1)
                  tq_2 = qgrid_in(2,i-istart+1,j-jstart+1,k-kstart+1,
     &                       iproc+1)
                  notfind = .false.
                      end if
                      iproc = iproc + 1
                      end do
                    end if
c
#if 1
                    t0   = theta1(1,it1)
                    t1   = theta1(2,it1)
                    t2   = theta1(3,it1)
#else
                    t0   = thetai1_t(1,it1,impi)
                    t1   = thetai1_t(2,it1,impi)
                    t2   = thetai1_t(3,it1,impi)
#endif
                    t0_1 = t0_1 + tq_1*t0
                    t1_1 = t1_1 + tq_1*t1
                    t2_1 = t2_1 + tq_1*t2
                    t0_2 = t0_2 + tq_2*t0
                    t1_2 = t1_2 + tq_2*t1
                    t2_2 = t2_2 + tq_2*t2
                 end do
                 tu00_1 = tu00_1 + t0_1*u0
                 tu10_1 = tu10_1 + t1_1*u0
                 tu01_1 = tu01_1 + t0_1*u1
                 tu20_1 = tu20_1 + t2_1*u0
                 tu11_1 = tu11_1 + t1_1*u1
                 tu02_1 = tu02_1 + t0_1*u2
                 tu00_2 = tu00_2 + t0_2*u0
                 tu10_2 = tu10_2 + t1_2*u0
                 tu01_2 = tu01_2 + t0_2*u1
                 tu20_2 = tu20_2 + t2_2*u0
                 tu11_2 = tu11_2 + t1_2*u1
                 tu02_2 = tu02_2 + t0_2*u2
              end do
              tuv1(1) = tuv1(1) + tu10_1*v0
              tuv1(2) = tuv1(2) + tu01_1*v0
              tuv1(3) = tuv1(3) + tu00_1*v1
              tuv1(4) = tuv1(4) + tu20_1*v0
              tuv1(5) = tuv1(5) + tu02_1*v0
              tuv1(6) = tuv1(6) + tu00_1*v2
              tuv1(7) = tuv1(7) + tu11_1*v0
              tuv1(8) = tuv1(8) + tu10_1*v1
              tuv1(9) = tuv1(9) + tu01_1*v1
              tuv2(1) = tuv2(1) + tu10_2*v0
              tuv2(2) = tuv2(2) + tu01_2*v0
              tuv2(3) = tuv2(3) + tu00_2*v1
              tuv2(4) = tuv2(4) + tu20_2*v0
              tuv2(5) = tuv2(5) + tu02_2*v0
              tuv2(6) = tuv2(6) + tu00_2*v2
              tuv2(7) = tuv2(7) + tu11_2*v0
              tuv2(8) = tuv2(8) + tu10_2*v1
              tuv2(9) = tuv2(9) + tu01_2*v1
           end do
c          do i = 2,10
c             fdip_phi1(i,impi) = tuv1(i-1)
c             fdip_phi2(i,impi) = tuv2(i-1)
c          end do
           fdip_phi1(02,impi) = tuv1(1)
           fdip_phi1(03,impi) = tuv1(2)
           fdip_phi1(04,impi) = tuv1(3)
           fdip_phi1(05,impi) = tuv1(4)
           fdip_phi1(06,impi) = tuv1(5)
           fdip_phi1(07,impi) = tuv1(6)
           fdip_phi1(08,impi) = tuv1(7)
           fdip_phi1(09,impi) = tuv1(8)
           fdip_phi1(10,impi) = tuv1(9)
           fdip_phi2(02,impi) = tuv2(1)
           fdip_phi2(03,impi) = tuv2(2)
           fdip_phi2(04,impi) = tuv2(3)
           fdip_phi2(05,impi) = tuv2(4)
           fdip_phi2(06,impi) = tuv2(5)
           fdip_phi2(07,impi) = tuv2(6)
           fdip_phi2(08,impi) = tuv2(7)
           fdip_phi2(09,impi) = tuv2(8)
           fdip_phi2(10,impi) = tuv2(9)
        end do
        end

        attributes(global) subroutine fphi_uind_sitecu2_core_1p
     &             (kstat,ked,jstat,jed,istat,ied
     &             ,npolerecloc,nlocrec,order,n1mpimax,n2mpimax,n3mpimax
     &             ,nfft1,nfft2,nfft3
     &             ,qgrid_in,polerecglob,ipole,igrid,x,y,z,recip
     &             ,thetai1,thetai2,thetai3
     &             ,fdip_phi1,fdip_phi2,fir)
        implicit none

        integer  ,intent(in),value:: npolerecloc,nlocrec,order
     &           ,n1mpimax,n2mpimax,n3mpimax
     &           ,nfft1,nfft2,nfft3
        integer  ,intent(in),value:: kstat,ked,jstat,jed,istat,ied
        integer  ,intent(in),device::polerecglob(*),ipole(*),igrid(3,*)
        real(t_p),intent(in),device::x(*),y(*),z(*),recip(3,3)
     &           ,thetai1(4,order,nlocrec),thetai2(4,order,nlocrec)
     &           ,thetai3(4,order,nlocrec)
     &           ,qgrid_in(2,n1mpimax,n2mpimax,n3mpimax)
        real(t_p),device,dimension(10,npolerecloc)::fdip_phi1,fdip_phi2
        logical,value,intent(in):: fir

        integer istart,iend,jstart,jend,kstart,kend
        integer i,j,k,impi,ifr
        integer iproc,proc
        integer isite,iatm
        integer i0,j0,k0
        integer it1,it2,it3
        integer igrd0,jgrd0,kgrd0
        integer ind
        logical notfind
        real(t_p) v0,v1,v2,v3
        real(t_p) u0,u1,u2,u3
        real(t_p) t0,t1,t2,t3
        real(t_p) t0_1,t0_2,t1_1,t1_2
        real(t_p) t2_1,t2_2,tq_1,tq_2
        real(t_p) tu00,tu10,tu01,tu20,tu11
        real(t_p) tu02,tu30,tu21,tu12,tu03
        real(t_p) tu00_1,tu01_1,tu10_1
        real(t_p) tu00_2,tu01_2,tu10_2
        real(t_p) tu20_1,tu11_1,tu02_1
        real(t_p) tu20_2,tu11_2,tu02_2
        real(t_p) tuv100_1,tuv010_1,tuv001_1
        real(t_p) tuv100_2,tuv010_2,tuv001_2
        real(t_p) tuv200_1,tuv020_1,tuv002_1
        real(t_p) tuv110_1,tuv101_1,tuv011_1
        real(t_p) tuv200_2,tuv020_2,tuv002_2
        real(t_p) tuv110_2,tuv101_2,tuv011_2
        real(t_p) tuv000,tuv100,tuv010,tuv001
        real(t_p) tuv200,tuv020,tuv002,tuv110
        real(t_p) tuv101,tuv011,tuv300,tuv030
        real(t_p) tuv003,tuv210,tuv201,tuv120
        real(t_p) tuv021,tuv102,tuv012,tuv111
        real(t_p) tuv1(9),tuv2(9)
        real(t_p) xi,yi,zi
        real(t_p) w,fr
        real(t_p),dimension(4,bsorder):: theta1,theta2,theta3
        real(t_p),shared::temp(bsorder*bsorder*PME_BLOCK_DIM)

c
        do impi = threadIdx%x+(blockIdx%x-1)*blockDim%x,npolerecloc,
     &            blockDim%x*gridDim%x
           !iipole = polerecglob(impi)  !no need with one MPI process
           iatm   = ipole(impi)
#if 1
c
c       get the b-spline coefficients for the i-th atomic site
c       it faster to recompute theta*
c
           xi     = x_t(iatm)
           yi     = y_t(iatm)
           zi     = z_t(iatm)
           w      = xi*recip_c(1,1) + yi*recip_c(2,1) + zi*recip_c(3,1)
           fr     = nfft1 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           igrd0  = ifr - bsorder
           call ibsplgen (w,theta1,temp((threadIdx%x-1)*bsorder**2+1))
           w      = xi*recip_c(1,2) + yi*recip_c(2,2) + zi*recip_c(3,2)
           fr     = nfft2 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           jgrd0  = ifr - bsorder
           call ibsplgen (w,theta2,temp((threadIdx%x-1)*bsorder**2+1))
           w      = xi*recip_c(1,3) + yi*recip_c(2,3) + zi*recip_c(3,3)
           fr     = nfft3 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           kgrd0  = ifr - bsorder
           call ibsplgen (w,theta3,temp((threadIdx%x-1)*bsorder**2+1))
#else
           igrd0  = igrid(1,iatm)
           jgrd0  = igrid(2,iatm)
           kgrd0  = igrid(3,iatm)
#endif
              tuv1(01) = 0.0_ti_p; tuv2(01) = 0.0_ti_p
              tuv1(02) = 0.0_ti_p; tuv2(02) = 0.0_ti_p
              tuv1(03) = 0.0_ti_p; tuv2(03) = 0.0_ti_p
              tuv1(04) = 0.0_ti_p; tuv2(04) = 0.0_ti_p
              tuv1(05) = 0.0_ti_p; tuv2(05) = 0.0_ti_p
              tuv1(06) = 0.0_ti_p; tuv2(06) = 0.0_ti_p
              tuv1(07) = 0.0_ti_p; tuv2(07) = 0.0_ti_p
              tuv1(08) = 0.0_ti_p; tuv2(08) = 0.0_ti_p
              tuv1(09) = 0.0_ti_p; tuv2(09) = 0.0_ti_p
c          k0       = kgrd0
           do it3 = 1, bsorder
c             k0     = k0 + 1
c             k      = k0 + 1 + (nfft3-isign(nfft3,kgrd0+it3))/2
              k      = kgrd0 +it3 +1 + 
     &                 ishft(nfft3-isign(nfft3,kgrd0+it3),-1)
#if 1
              v0     = theta3(1,it3)
              v1     = theta3(2,it3)
              v2     = theta3(3,it3)
#else
              v0     = thetai3_t(1,it3,impi)
              v1     = thetai3_t(2,it3,impi)
              v2     = thetai3_t(3,it3,impi)
#endif  
              tu00_1 = 0.0_ti_p
              tu01_1 = 0.0_ti_p
              tu10_1 = 0.0_ti_p
              tu20_1 = 0.0_ti_p
              tu11_1 = 0.0_ti_p
              tu02_1 = 0.0_ti_p
              tu00_2 = 0.0_ti_p
              tu01_2 = 0.0_ti_p
              tu10_2 = 0.0_ti_p
              tu20_2 = 0.0_ti_p
              tu11_2 = 0.0_ti_p
              tu02_2 = 0.0_ti_p
c             j0 = jgrd0
              do it2 = 1, bsorder
c                j0 = j0 + 1
c                j    = jgrd0 +it2 + 1 + (nfft2-isign(nfft2,jgrd0+it2))/2
                 j    = jgrd0 +it2 + 1 +
     &                        ishft(nfft2-isign(nfft2,jgrd0+it2),-1)
#if 1
                 u0   = theta2(1,it2)
                 u1   = theta2(2,it2)
                 u2   = theta2(3,it2)
#else
                 u0   = thetai2_t(1,it2,impi)
                 u1   = thetai2_t(2,it2,impi)
                 u2   = thetai2_t(3,it2,impi)
#endif  
                 t0_1 = 0.0_ti_p
                 t1_1 = 0.0_ti_p
                 t2_1 = 0.0_ti_p
                 t0_2 = 0.0_ti_p
                 t1_2 = 0.0_ti_p
                 t2_2 = 0.0_ti_p
                 i0   = igrd0
                 do it1 = 1, bsorder
c                   i0  = i0 + 1
c                   i    = igrd0 +it1 +1 + (nfft1-isign(nfft1,igrd0 +it1))/2
                    i    = igrd0 +it1 +1 +
     &                           ishft(nfft1-isign(nfft1,igrd0+it1),-1)
c
                    tq_1 = qgrid_in(1,i-istat+1,j-jstat+1,k-kstat+1)
                    tq_2 = qgrid_in(2,i-istat+1,j-jstat+1,k-kstat+1)
c
#if 1
                    t0   = theta1(1,it1)
                    t1   = theta1(2,it1)
                    t2   = theta1(3,it1)
#else
                    t0   = thetai1_t(1,it1,impi)
                    t1   = thetai1_t(2,it1,impi)
                    t2   = thetai1_t(3,it1,impi)
#endif
                    t0_1 = t0_1 + tq_1*t0
                    t1_1 = t1_1 + tq_1*t1
                    t2_1 = t2_1 + tq_1*t2
                    t0_2 = t0_2 + tq_2*t0
                    t1_2 = t1_2 + tq_2*t1
                    t2_2 = t2_2 + tq_2*t2
                 end do
                 tu00_1 = tu00_1 + t0_1*u0
                 tu10_1 = tu10_1 + t1_1*u0
                 tu01_1 = tu01_1 + t0_1*u1
                 tu20_1 = tu20_1 + t2_1*u0
                 tu11_1 = tu11_1 + t1_1*u1
                 tu02_1 = tu02_1 + t0_1*u2
                 tu00_2 = tu00_2 + t0_2*u0
                 tu10_2 = tu10_2 + t1_2*u0
                 tu01_2 = tu01_2 + t0_2*u1
                 tu20_2 = tu20_2 + t2_2*u0
                 tu11_2 = tu11_2 + t1_2*u1
                 tu02_2 = tu02_2 + t0_2*u2
              end do
              tuv1(1) = tuv1(1) + tu10_1*v0
              tuv1(2) = tuv1(2) + tu01_1*v0
              tuv1(3) = tuv1(3) + tu00_1*v1
              tuv1(4) = tuv1(4) + tu20_1*v0
              tuv1(5) = tuv1(5) + tu02_1*v0
              tuv1(6) = tuv1(6) + tu00_1*v2
              tuv1(7) = tuv1(7) + tu11_1*v0
              tuv1(8) = tuv1(8) + tu10_1*v1
              tuv1(9) = tuv1(9) + tu01_1*v1
              tuv2(1) = tuv2(1) + tu10_2*v0
              tuv2(2) = tuv2(2) + tu01_2*v0
              tuv2(3) = tuv2(3) + tu00_2*v1
              tuv2(4) = tuv2(4) + tu20_2*v0
              tuv2(5) = tuv2(5) + tu02_2*v0
              tuv2(6) = tuv2(6) + tu00_2*v2
              tuv2(7) = tuv2(7) + tu11_2*v0
              tuv2(8) = tuv2(8) + tu10_2*v1
              tuv2(9) = tuv2(9) + tu01_2*v1
           end do
           fdip_phi1(02,impi) = tuv1(1)
           fdip_phi1(03,impi) = tuv1(2)
           fdip_phi1(04,impi) = tuv1(3)
           fdip_phi1(05,impi) = tuv1(4)
           fdip_phi1(06,impi) = tuv1(5)
           fdip_phi1(07,impi) = tuv1(6)
           fdip_phi1(08,impi) = tuv1(7)
           fdip_phi1(09,impi) = tuv1(8)
           fdip_phi1(10,impi) = tuv1(9)
           fdip_phi2(02,impi) = tuv2(1)
           fdip_phi2(03,impi) = tuv2(2)
           fdip_phi2(04,impi) = tuv2(3)
           fdip_phi2(05,impi) = tuv2(4)
           fdip_phi2(06,impi) = tuv2(5)
           fdip_phi2(07,impi) = tuv2(6)
           fdip_phi2(08,impi) = tuv2(7)
           fdip_phi2(09,impi) = tuv2(8)
           fdip_phi2(10,impi) = tuv2(9)
        end do
        end


        attributes(global) subroutine fphi_uind_sitecu1_core
     &             (kstat,ked,jstat,jed,istat,ied,nfft1,nfft2,nfft3
     &             ,npolerecloc,nlocrec,order,nrec_send,nproc
     &             ,polerecglob,ipole,igrid
     &             ,fdip_phi1,fdip_phi2,fdip_sum_phi)
        implicit none
        integer  ,intent(in),value:: npolerecloc,nlocrec,order
     &           ,nfft1,nfft2,nfft3,nrec_send,nproc
        integer  ,intent(in),value:: kstat,ked,jstat,jed,istat,ied
        integer  ,intent(in),device::polerecglob(*),ipole(*),igrid(3,*)
c       real(t_p),intent(in),device::thetai1(4,order,nlocrec)
c    &           ,thetai2(4,order,nlocrec),thetai3(4,order,nlocrec)
        real(t_p),device,dimension(10,npolerecloc)::fdip_phi1,fdip_phi2
        real(t_p),device:: fdip_sum_phi(20,npolerecloc)
        
 
        integer istart,iend,jstart,jend,kstart,kend
        integer i,j,k,impi
        integer iproc,proc
        integer isite,iatm,iipole
        integer i0,j0,k0
        integer it1,it2,it3
        integer igrd0,jgrd0,kgrd0
        real(t_p) v0,v1,v2,v3
        real(t_p) u0,u1,u2,u3
        real(t_p) t0,t1,t2,t3
        real(t_p) t0_1,t0_2,t1_1,t1_2
        real(t_p) t2_1,t2_2,tq_1,tq_2
        real(t_p) tu00,tu10,tu01,tu20,tu11
        real(t_p) tu02,tu30,tu21,tu12,tu03
        real(t_p) tu00_1,tu01_1,tu10_1
        real(t_p) tu00_2,tu01_2,tu10_2
        real(t_p) tu20_1,tu11_1,tu02_1
        real(t_p) tu20_2,tu11_2,tu02_2
        real(t_p) tuv100_1,tuv010_1,tuv001_1
        real(t_p) tuv100_2,tuv010_2,tuv001_2
        real(t_p) tuv200_1,tuv020_1,tuv002_1
        real(t_p) tuv110_1,tuv101_1,tuv011_1
        real(t_p) tuv200_2,tuv020_2,tuv002_2
        real(t_p) tuv110_2,tuv101_2,tuv011_2
        real(t_p) tuv000,tuv100,tuv010,tuv001
        real(t_p) tuv200,tuv020,tuv002,tuv110
        real(t_p) tuv101,tuv011,tuv300,tuv030
        real(t_p) tuv003,tuv210,tuv201,tuv120
        real(t_p) tuv021,tuv102,tuv012,tuv111
        real(t_p) xi,yi,zi
        real(t_p) w,fr
        integer ifr
        real(t_p),dimension(4,bsorder):: theta1,theta2,theta3
        real(t_p),shared::temp(bsorder*bsorder*PME_BLOCK_DIM)
c       
        do impi = threadIdx%x+(blockIdx%x-1)*blockDim%x,npolerecloc,
     &            blockDim%x*gridDim%x
           iipole = polerecglob(impi)
           iatm   = ipole(iipole)
#if 1
c
c       get the b-spline coefficients for the i-th atomic site
c       it faster to recompute theta*
c
           xi     = x_t(iatm)
           yi     = y_t(iatm)
           zi     = z_t(iatm)
           w      = xi*recip_c(1,1) + yi*recip_c(2,1) + zi*recip_c(3,1)
           fr     = nfft1 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           igrd0  = ifr - bsorder
           call ibsplgen (w,theta1,temp((threadIdx%x-1)*bsorder**2+1))
           w      = xi*recip_c(1,2) + yi*recip_c(2,2) + zi*recip_c(3,2)
           fr     = nfft2 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           jgrd0  = ifr - bsorder
           call ibsplgen (w,theta2,temp((threadIdx%x-1)*bsorder**2+1))
           w      = xi*recip_c(1,3) + yi*recip_c(2,3) + zi*recip_c(3,3)
           fr     = nfft3 * (w-anint(w)+0.5_ti_p)
           ifr    = int(fr-pme_eps)
           w      = fr - real(ifr,t_p)
           kgrd0  = ifr - bsorder
           call ibsplgen (w,theta3,temp((threadIdx%x-1)*bsorder**2+1))
#else
           igrd0  = igrid(1,iatm)
           jgrd0  = igrid(2,iatm)
           kgrd0  = igrid(3,iatm)
#endif
           tuv100_1 = 0.0_ti_p
           tuv010_1 = 0.0_ti_p
           tuv001_1 = 0.0_ti_p
           tuv200_1 = 0.0_ti_p
           tuv020_1 = 0.0_ti_p
           tuv002_1 = 0.0_ti_p
           tuv110_1 = 0.0_ti_p
           tuv101_1 = 0.0_ti_p
           tuv011_1 = 0.0_ti_p
           tuv100_2 = 0.0_ti_p
           tuv010_2 = 0.0_ti_p
           tuv001_2 = 0.0_ti_p
           tuv200_2 = 0.0_ti_p
           tuv020_2 = 0.0_ti_p
           tuv002_2 = 0.0_ti_p
           tuv110_2 = 0.0_ti_p
           tuv101_2 = 0.0_ti_p
           tuv011_2 = 0.0_ti_p
           tuv000   = 0.0_ti_p
           tuv001   = 0.0_ti_p
           tuv010   = 0.0_ti_p
           tuv100   = 0.0_ti_p
           tuv200   = 0.0_ti_p
           tuv020   = 0.0_ti_p
           tuv002   = 0.0_ti_p
           tuv110   = 0.0_ti_p
           tuv101   = 0.0_ti_p
           tuv011   = 0.0_ti_p
           tuv300   = 0.0_ti_p
           tuv030   = 0.0_ti_p
           tuv003   = 0.0_ti_p
           tuv210   = 0.0_ti_p
           tuv201   = 0.0_ti_p
           tuv120   = 0.0_ti_p
           tuv021   = 0.0_ti_p
           tuv102   = 0.0_ti_p
           tuv012   = 0.0_ti_p
           tuv111   = 0.0_ti_p
           k0       = kgrd0
!$acc      loop seq
           do it3 = 1, bsorder
              k0 = k0 + 1
c             k      = k0 + 1 + (nfft3-isign(nfft3,k0))/2
              k      = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
#if 1
              v0     = theta3(1,it3)
              v1     = theta3(2,it3)
              v2     = theta3(3,it3)
              v3     = theta3(4,it3)
#else
              v0     = thetai3_t(1,it3,impi)
              v1     = thetai3_t(2,it3,impi)
              v2     = thetai3_t(3,it3,impi)
              v3     = thetai3_t(4,it3,impi)
#endif
              tu00_1 = 0.0_ti_p
              tu01_1 = 0.0_ti_p
              tu10_1 = 0.0_ti_p
              tu20_1 = 0.0_ti_p
              tu11_1 = 0.0_ti_p
              tu02_1 = 0.0_ti_p
              tu00_2 = 0.0_ti_p
              tu01_2 = 0.0_ti_p
              tu10_2 = 0.0_ti_p
              tu20_2 = 0.0_ti_p
              tu11_2 = 0.0_ti_p
              tu02_2 = 0.0_ti_p
              tu00   = 0.0_ti_p
              tu10   = 0.0_ti_p
              tu01   = 0.0_ti_p
              tu20   = 0.0_ti_p
              tu11   = 0.0_ti_p
              tu02   = 0.0_ti_p
              tu30   = 0.0_ti_p
              tu21   = 0.0_ti_p
              tu12   = 0.0_ti_p
              tu03   = 0.0_ti_p
              j0 = jgrd0
!$acc      loop seq
              do it2 = 1, bsorder
                 j0 = j0 + 1
c                j    = j0 + 1 + (nfft2-isign(nfft2,j0))/2
                 j    = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
#if 1
                 u0   = theta2(1,it2)
                 u1   = theta2(2,it2)
                 u2   = theta2(3,it2)
                 u3   = theta2(4,it2)
#else
                 u0   = thetai2_t(1,it2,impi)
                 u1   = thetai2_t(2,it2,impi)
                 u2   = thetai2_t(3,it2,impi)
                 u3   = thetai2_t(4,it2,impi)
#endif
                 t0_1 = 0.0_ti_p
                 t1_1 = 0.0_ti_p
                 t2_1 = 0.0_ti_p
                 t0_2 = 0.0_ti_p
                 t1_2 = 0.0_ti_p
                 t2_2 = 0.0_ti_p
                 t3   = 0.0_ti_p
                 i0   = igrd0
!$acc      loop seq
                 do it1 = 1, bsorder
                    i0 = i0 + 1
c                   i    = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                    i    = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)
c        
                    tq_1 = 0_ti_p
                    tq_2 = 0_ti_p
                    if (((k.ge.kstat).and.(k.le.ked)).and.
     &                  ((j.ge.jstat).and.(j.le.jed)).and.
     &                  ((i.ge.istat).and.(i.le.ied))) then
                  tq_1 = qgrid2in_t(1,i-istat+1,j-jstat+1,k-kstat+1,1)
                  tq_2 = qgrid2in_t(2,i-istat+1,j-jstat+1,k-kstat+1,1)
                      goto 10
                    end if
!$acc      loop seq
                    do iproc = 1, nrec_send
                      proc   = prec_send_t(iproc)
                      kstart = kstart1_t(proc+1)
                      kend   = kend1_t  (proc+1)
                      jstart = jstart1_t(proc+1)
                      jend   = jend1_t  (proc+1)
                      istart = istart1_t(proc+1)
                      iend   = iend1_t  (proc+1)
                      if (((k.ge.kstart).and.(k.le.kend)).and.
     &                    ((j.ge.jstart).and.(j.le.jend)).and.
     &                    ((i.ge.istart).and.(i.le.iend))) then
                  tq_1 = qgrid2in_t(1,i-istart+1,j-jstart+1,k-kstart+1,
     &                       iproc+1)
                  tq_2 = qgrid2in_t(2,i-istart+1,j-jstart+1,k-kstart+1,
     &                       iproc+1)
                        goto 10
                      end if
                    end do
                    cycle
 10                 continue
c       
#if 1
                    t0_1 = t0_1 +  tq_1      *theta1(1,it1)
                    t1_1 = t1_1 +  tq_1      *theta1(2,it1)
                    t2_1 = t2_1 +  tq_1      *theta1(3,it1)
                    t0_2 = t0_2 +        tq_2*theta1(1,it1)
                    t1_2 = t1_2 +        tq_2*theta1(2,it1)
                    t2_2 = t2_2 +        tq_2*theta1(3,it1)
                    t3   = t3   + (tq_1+tq_2)*theta1(4,it1)
#else
                    t0_1 = t0_1 +  tq_1      *thetai1_t(1,it1,impi)
                    t1_1 = t1_1 +  tq_1      *thetai1_t(2,it1,impi)
                    t2_1 = t2_1 +  tq_1      *thetai1_t(3,it1,impi)
                    t0_2 = t0_2 +        tq_2*thetai1_t(1,it1,impi)
                    t1_2 = t1_2 +        tq_2*thetai1_t(2,it1,impi)
                    t2_2 = t2_2 +        tq_2*thetai1_t(3,it1,impi)
                    t3   = t3   + (tq_1+tq_2)*thetai1_t(4,it1,impi)
#endif
                 end do
                 tu00_1 = tu00_1 + t0_1*u0
                 tu10_1 = tu10_1 + t1_1*u0
                 tu01_1 = tu01_1 + t0_1*u1
                 tu20_1 = tu20_1 + t2_1*u0
                 tu11_1 = tu11_1 + t1_1*u1
                 tu02_1 = tu02_1 + t0_1*u2
                 tu00_2 = tu00_2 + t0_2*u0
                 tu10_2 = tu10_2 + t1_2*u0
                 tu01_2 = tu01_2 + t0_2*u1
                 tu20_2 = tu20_2 + t2_2*u0
                 tu11_2 = tu11_2 + t1_2*u1
                 tu02_2 = tu02_2 + t0_2*u2
                 t0     = t0_1 + t0_2
                 t1     = t1_1 + t1_2
                 t2     = t2_1 + t2_2
                 tu00   = tu00 + t0*u0
                 tu10   = tu10 + t1*u0
                 tu01   = tu01 + t0*u1
                 tu20   = tu20 + t2*u0
                 tu11   = tu11 + t1*u1
                 tu02   = tu02 + t0*u2
                 tu30   = tu30 + t3*u0
                 tu21   = tu21 + t2*u1
                 tu12   = tu12 + t1*u2
                 tu03   = tu03 + t0*u3
              end do
              tuv100_1 = tuv100_1 + tu10_1*v0
              tuv010_1 = tuv010_1 + tu01_1*v0
              tuv001_1 = tuv001_1 + tu00_1*v1
              tuv200_1 = tuv200_1 + tu20_1*v0
              tuv020_1 = tuv020_1 + tu02_1*v0
              tuv002_1 = tuv002_1 + tu00_1*v2
              tuv110_1 = tuv110_1 + tu11_1*v0
              tuv101_1 = tuv101_1 + tu10_1*v1
              tuv011_1 = tuv011_1 + tu01_1*v1
              tuv100_2 = tuv100_2 + tu10_2*v0
              tuv010_2 = tuv010_2 + tu01_2*v0
              tuv001_2 = tuv001_2 + tu00_2*v1
              tuv200_2 = tuv200_2 + tu20_2*v0
              tuv020_2 = tuv020_2 + tu02_2*v0
              tuv002_2 = tuv002_2 + tu00_2*v2
              tuv110_2 = tuv110_2 + tu11_2*v0
              tuv101_2 = tuv101_2 + tu10_2*v1
              tuv011_2 = tuv011_2 + tu01_2*v1
              tuv000   = tuv000 + tu00*v0
              tuv100   = tuv100 + tu10*v0
              tuv010   = tuv010 + tu01*v0
              tuv001   = tuv001 + tu00*v1
              tuv200   = tuv200 + tu20*v0
              tuv020   = tuv020 + tu02*v0
              tuv002   = tuv002 + tu00*v2
              tuv110   = tuv110 + tu11*v0
              tuv101   = tuv101 + tu10*v1
              tuv011   = tuv011 + tu01*v1
              tuv300   = tuv300 + tu30*v0
              tuv030   = tuv030 + tu03*v0
              tuv003   = tuv003 + tu00*v3
              tuv210   = tuv210 + tu21*v0
              tuv201   = tuv201 + tu20*v1
              tuv120   = tuv120 + tu12*v0
              tuv021   = tuv021 + tu02*v1
              tuv102   = tuv102 + tu10*v2
              tuv012   = tuv012 + tu01*v2
              tuv111   = tuv111 + tu11*v1
           end do
           fdip_phi1   ( 2,impi) = tuv100_1
           fdip_phi1   ( 3,impi) = tuv010_1
           fdip_phi1   ( 4,impi) = tuv001_1
           fdip_phi1   ( 5,impi) = tuv200_1
           fdip_phi1   ( 6,impi) = tuv020_1
           fdip_phi1   ( 7,impi) = tuv002_1
           fdip_phi1   ( 8,impi) = tuv110_1
           fdip_phi1   ( 9,impi) = tuv101_1
           fdip_phi1   (10,impi) = tuv011_1
           fdip_phi2   ( 2,impi) = tuv100_2
           fdip_phi2   ( 3,impi) = tuv010_2
           fdip_phi2   ( 4,impi) = tuv001_2
           fdip_phi2   ( 5,impi) = tuv200_2
           fdip_phi2   ( 6,impi) = tuv020_2
           fdip_phi2   ( 7,impi) = tuv002_2
           fdip_phi2   ( 8,impi) = tuv110_2
           fdip_phi2   ( 9,impi) = tuv101_2
           fdip_phi2   (10,impi) = tuv011_2
           fdip_sum_phi( 1,impi) = tuv000
           fdip_sum_phi( 2,impi) = tuv100
           fdip_sum_phi( 3,impi) = tuv010
           fdip_sum_phi( 4,impi) = tuv001
           fdip_sum_phi( 5,impi) = tuv200
           fdip_sum_phi( 6,impi) = tuv020
           fdip_sum_phi( 7,impi) = tuv002
           fdip_sum_phi( 8,impi) = tuv110
           fdip_sum_phi( 9,impi) = tuv101
           fdip_sum_phi(10,impi) = tuv011
           fdip_sum_phi(11,impi) = tuv300
           fdip_sum_phi(12,impi) = tuv030
           fdip_sum_phi(13,impi) = tuv003
           fdip_sum_phi(14,impi) = tuv210
           fdip_sum_phi(15,impi) = tuv201
           fdip_sum_phi(16,impi) = tuv120
           fdip_sum_phi(17,impi) = tuv021
           fdip_sum_phi(18,impi) = tuv102
           fdip_sum_phi(19,impi) = tuv012
           fdip_sum_phi(20,impi) = tuv111
        end do
      end subroutine

      attributes(global) subroutine fphi_chg_site_kcu
     &                  (kind_id,atom_id,igrid
     &                  ,kstat,ked,jstat,jed,istat,ied
     &                  ,nrec_send,nfft1,nfft2,nfft3,nionrecloc,n
     &                  ,fphirec)
      implicit none
      integer  ,value,intent(in)::kstat,ked,jstat,jed,istat,ied
     &         ,nrec_send,nfft1,nfft2,nfft3,nionrecloc,n
      !real(t_p),value,intent(in)::
      integer  ,device,intent(in)::kind_id(nionrecloc),atom_id(n)
     &         ,igrid(3,n)
      real(t_p),device:: fphirec(4,*)

      integer   i,j,k,impi,rankloc,iglob,iatm,iichg,i0,j0,k0
     &         ,it1,it2,it3,igrd0,jgrd0,kgrd0,proc
     &         ,iproc,istart,iend,jstart,jend,kstart,kend
      real(t_p) v0,v1,u0,u1,t0,t1,tq,tu00,tu10,tu01
     &         ,tuv000,tuv100,tuv010,tuv001
#if 1
      integer ifr
      real(t_p),xi,yi,zi,w,fr
      real(t_p),dimension(2*bsorder):: theta1,theta2,theta3
      real(t_p),shared::temp(bsorder*bsorder*PME_BLOCK_DIM)
#endif
      do impi = (blockIdx%x-1)*blockDim%x + threadIdx%x, nionrecloc,
     &           blockDim%x*gridDim%x

         iichg  = kind_id(impi)
         iglob  = atom_id(iichg)
#if 1
c
c       get the b-spline coefficients for the i-th atomic site
c       it faster to recompute theta*
c
         xi     = x_t(iglob)
         yi     = y_t(iglob)
         zi     = z_t(iglob)
         w      = xi*recip_c(1,1) + yi*recip_c(2,1) + zi*recip_c(3,1)
         fr     = nfft1 * (w-anint(w)+0.5_ti_p)
         ifr    = int(fr-pme_eps)
         w      = fr - real(ifr,t_p)
         igrd0  = ifr - bsorder
         call ibsplgen_2 (w,theta1,temp((threadIdx%x-1)*bsorder**2+1))
         w      = xi*recip_c(1,2) + yi*recip_c(2,2) + zi*recip_c(3,2)
         fr     = nfft2 * (w-anint(w)+0.5_ti_p)
         ifr    = int(fr-pme_eps)
         w      = fr - real(ifr,t_p)
         jgrd0  = ifr - bsorder
         call ibsplgen_2 (w,theta2,temp((threadIdx%x-1)*bsorder**2+1))
         w      = xi*recip_c(1,3) + yi*recip_c(2,3) + zi*recip_c(3,3)
         fr     = nfft3 * (w-anint(w)+0.5_ti_p)
         ifr    = int(fr-pme_eps)
         w      = fr - real(ifr,t_p)
         kgrd0  = ifr - bsorder
         call ibsplgen_2 (w,theta3,temp((threadIdx%x-1)*bsorder**2+1))
#else
         igrd0  = igrid(1,iglob)
         jgrd0  = igrid(2,iglob)
         kgrd0  = igrid(3,iglob)
#endif
         tuv000 = 0.0_ti_p
         tuv001 = 0.0_ti_p
         tuv010 = 0.0_ti_p
         tuv100 = 0.0_ti_p
         k0     = kgrd0
         do it3 = 1, bsorder
            k0   = k0 + 1
            k    = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
            v0   = theta3(1+(it3-1)*level1)
            v1   = theta3(2+(it3-1)*level1)
            tu00 = 0.0_ti_p
            tu10 = 0.0_ti_p
            tu01 = 0.0_ti_p
            j0   = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j  = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
               u0 = theta2(1+(it2-1)*level1)
               u1 = theta2(2+(it2-1)*level1)
               t0 = 0.0_ti_p
               t1 = 0.0_ti_p
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)
c
                  tq     = 0.0_ti_p
                  if (((k.ge.kstat).and.(k.le.ked)).and.
     $                ((j.ge.jstat).and.(j.le.jed)).and.
     $                ((i.ge.istat).and.(i.le.ied))) then
                     tq = qgridin_t(1,i-istat+1,j-jstat+1,k-kstat+1,1)
                     goto 10
                  end if
                  do iproc = 1, nrec_send
                    proc   = prec_send_t(iproc)
                    kstart = kstart1_t  (proc+1)
                    kend   = kend1_t    (proc+1)
                    jstart = jstart1_t  (proc+1)
                    jend   = jend1_t    (proc+1)
                    istart = istart1_t  (proc+1)
                    iend   = iend1_t    (proc+1)
                    if (((k.ge.kstart).and.(k.le.kend)).and.
     $                  ((j.ge.jstart).and.(j.le.jend)).and.
     $                  ((i.ge.istart).and.(i.le.iend))) then
                       tq= qgridin_t(1,i-istart+1,j-jstart+1,k-kstart+1
     $                             ,iproc+1)
                       goto 10
                    end if
                  end do
                  cycle
c
 10               continue
                  t0 = t0 + tq*theta1(1+(it1-1)*level1)
                  t1 = t1 + tq*theta1(2+(it1-1)*level1)
               end do
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
            end do
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
         end do
         fphirec(1,impi) = tuv000
         fphirec(2,impi) = tuv100
         fphirec(3,impi) = tuv010
         fphirec(4,impi) = tuv001
      end do
      end subroutine

      attributes(global) subroutine grid_calc_frc_kcu
     &                 (kind_id,atom_id,locrec,igrid
     &                 ,attrb,thetai1,thetai2,thetai3
     &                 ,derec
     &                 ,kstat,ked,jstat,jed,istat,ied
     &                 ,nrec_send,nionrecloc,n,nfft1,nfft2,nfft3
     &                 ,f,scal,dn1,dn2,dn3)
      implicit none
      integer,value,intent(in)::kstat,ked,jstat,jed,istat,ied,nrec_send
     &       ,nfft1,nfft2,nfft3,nionrecloc,n
      real(t_p),value,intent(in):: dn1,dn2,dn3,f,scal
      integer,device,intent(in)::kind_id(nionrecloc),atom_id(n)
     &       ,locrec(n),igrid(3,n)
      real(t_p),device,intent(in):: attrb(n),thetai1(2,bsorder,*)
     &         ,thetai2(2,bsorder,*),thetai3(2,bsorder,*)
      real(r_p),device,intent(out):: derec(3,*)

      integer iichg,iglob,iloc,isite,rankloc
     &       ,igrd0,jgrd0,kgrd0,i,j,k,i0,j0,k0,it1,it2,it3
     &       ,istart,iend,jstart,jend,kstart,kend
     &       ,iproc,proc
      real(t_p) fi,de1,de2,de3,t1,t2,t3,dt1,dt2,dt3,term
#if 1
      integer ifr
      real(t_p),xi,yi,zi,w,fr
      real(t_p),dimension(2*bsorder):: theta1,theta2,theta3
      real(t_p),shared::temp(bsorder*bsorder*PME_BLOCK_DIM)
#endif

      do isite = (blockIdx%x-1)*blockDim%x + threadIdx%x, nionrecloc,
     &            blockDim%x*gridDim%x
        iichg = kind_id (isite)
        iglob = atom_id   (iichg)
        iloc  = locrec (iglob)
#if 1
c
c       get the b-spline coefficients for the i-th atomic site
c       it faster to recompute theta*
c
        xi     = x_t(iglob)
        yi     = y_t(iglob)
        zi     = z_t(iglob)
        w      = xi*recip_c(1,1) + yi*recip_c(2,1) + zi*recip_c(3,1)
        fr     = nfft1 * (w-anint(w)+0.5_ti_p)
        ifr    = int(fr-pme_eps)
        w      = fr - real(ifr,t_p)
        igrd0  = ifr - bsorder
        call ibsplgen_2 (w,theta1,temp((threadIdx%x-1)*bsorder**2+1))
        w      = xi*recip_c(1,2) + yi*recip_c(2,2) + zi*recip_c(3,2)
        fr     = nfft2 * (w-anint(w)+0.5_ti_p)
        ifr    = int(fr-pme_eps)
        w      = fr - real(ifr,t_p)
        jgrd0  = ifr - bsorder
        call ibsplgen_2 (w,theta2,temp((threadIdx%x-1)*bsorder**2+1))
        w      = xi*recip_c(1,3) + yi*recip_c(2,3) + zi*recip_c(3,3)
        fr     = nfft3 * (w-anint(w)+0.5_ti_p)
        ifr    = int(fr-pme_eps)
        w      = fr - real(ifr,t_p)
        kgrd0  = ifr - bsorder
        call ibsplgen_2 (w,theta3,temp((threadIdx%x-1)*bsorder**2+1))
#else
        igrd0  = igrid(1,iglob)
        jgrd0  = igrid(2,iglob)
        kgrd0  = igrid(3,iglob)
#endif
        fi     = f * attrb(iichg)
        de1    = 0.0_ti_p
        de2    = 0.0_ti_p
        de3    = 0.0_ti_p
        k0     = kgrd0
        do it3 = 1, bsorder
           k0  = k0 + 1
           k   = k0 + 1 + (nfft3-sign(nfft3,k0))/2
#if 1
           t3  =       theta3(1+(it3-1)*level1)
           dt3 = dn3 * theta3(2+(it3-1)*level1)
#else
           t3  =       thetai3(1,it3,iglob)
           dt3 = dn3 * thetai3(2,it3,iglob)
#endif
           j0  = jgrd0
           do it2 = 1, bsorder
              j0  = j0 + 1
              j   = j0 + 1 + (nfft2-sign(nfft2,j0))/2
#if 1
              t2  =       theta2(1+(it2-1)*level1)
              dt2 = dn2 * theta2(2+(it2-1)*level1)
#else
              t2  =       thetai2(1,it2,iglob)
              dt2 = dn2 * thetai2(2,it2,iglob)
#endif
              i0  = igrd0
              do it1 = 1, bsorder
                 i0  = i0 + 1
                 i   = i0 + 1 + (nfft1-sign(nfft1,i0))/2
#if 1
                 t1  =       theta1(1+(it1-1)*level1)
                 dt1 = dn1 * theta1(2+(it1-1)*level1)
#else
                 t1  =       thetai1(1,it1,iglob)
                 dt1 = dn1 * thetai1(2,it1,iglob)
#endif
c
                 if (((k.ge.kstat).and.(k.le.ked)).and.
     $               ((j.ge.jstat).and.(j.le.jed)).and.
     $               ((i.ge.istat).and.(i.le.ied)))  then
                    term = qgridin_t(1,i-istat+1,j-jstat+1
     $                                ,k-kstat+1,1)
                    goto 100
                 end if
c
                 do iproc = 1, nrec_send
                    proc   = prec_send_t(iproc)
                    kstart =   kstart1_t(proc+1)
                    kend   =     kend1_t(proc+1)
                    jstart =   jstart1_t(proc+1)
                    jend   =     jend1_t(proc+1)
                    istart =   istart1_t(proc+1)
                    iend   =     iend1_t(proc+1)
                    if (((k.ge.kstart).and.(k.le.kend)).and.
     $                  ((j.ge.jstart).and.(j.le.jend)).and.
     $                  ((i.ge.istart).and.(i.le.iend))) then
                       term=qgridin_t(1,i-istart+1,j-jstart+1
     $                                 ,k-kstart+1,iproc+1)
                       goto 100
                    end if
                  end do
 100              continue
c
                 de1 = de1 + term*dt1*t2*t3
                 de2 = de2 + term*dt2*t1*t3
                 de3 = de3 + term*dt3*t1*t2
              end do
           end do
        end do
        derec(1,iloc) =derec(1,iloc)+fi*scal*(recip_c(1,1)*de1
     &                   +recip_c(1,2)*de2+recip_c(1,3)*de3)
        derec(2,iloc) =derec(2,iloc)+fi*scal*(recip_c(2,1)*de1
     &                   +recip_c(2,2)*de2+recip_c(2,3)*de3)
        derec(3,iloc) =derec(3,iloc)+fi*scal*(recip_c(3,1)*de1
     &                   +recip_c(3,2)*de2+recip_c(3,3)*de3)
      end do
      end subroutine

      ! Scalar sum for reciprocal space (point charge)
      attributes(global) subroutine grid_calc_frc_kcu1
     &                 (kind_id,atom_id,locrec,igrid
     &                 ,attrb,thetai1,thetai2,thetai3
     &                 ,derec
     &                 ,kstat,ked,jstat,jed,istat,ied
     &                 ,nrec_send,nionrecloc,n,nfft1,nfft2,nfft3
     &                 ,f,scal,dn1,dn2,dn3)
      implicit none
      integer,value,intent(in)::kstat,ked,jstat,jed,istat,ied,nrec_send
     &       ,nfft1,nfft2,nfft3,nionrecloc,n
      real(t_p),value,intent(in):: f,dn1,dn2,dn3,scal
      integer,device,intent(in)::kind_id(nionrecloc),atom_id(n)
     &       ,locrec(n),igrid(3,n)
      real(t_p),device,intent(in):: attrb(n),thetai1(2,bsorder,*)
     &         ,thetai2(2,bsorder,*),thetai3(2,bsorder,*)
      real(r_p),device,intent(out):: derec(3,*)

      integer iichg,iglob,iloc,isite,rankloc
     &       ,igrd0,jgrd0,kgrd0,i,j,k,i0,j0,k0,it1,it2,it3
     &       ,istart,iend,jstart,jend,kstart,kend
     &       ,iproc,proc
      real(t_p) fi,de1,de2,de3,t1,t2,t3,dt1,dt2,dt3,term
#if 1
      integer ifr
      real(t_p),xi,yi,zi,w,fr
      real(t_p),dimension(2*bsorder):: theta1,theta2,theta3
      real(t_p),shared::temp(bsorder*bsorder*PME_BLOCK_DIM)
#endif

      do isite = (blockIdx%x-1)*blockDim%x + threadIdx%x, nionrecloc,
     &            blockDim%x*gridDim%x
        iichg = kind_id (isite)
        iglob = atom_id (iichg)
#if 1
c
c       get the b-spline coefficients for the i-th atomic site
c       it faster to recompute theta*
c
        xi     = x_t(iglob)
        yi     = y_t(iglob)
        zi     = z_t(iglob)
        w      = xi*recip_c(1,1) + yi*recip_c(2,1) + zi*recip_c(3,1)
        fr     = nfft1 * (w-anint(w)+0.5_ti_p)
        ifr    = int(fr-pme_eps)
        w      = fr - real(ifr,t_p)
        igrd0  = ifr - bsorder
        call ibsplgen_2 (w,theta1,temp((threadIdx%x-1)*bsorder**2+1))
        w      = xi*recip_c(1,2) + yi*recip_c(2,2) + zi*recip_c(3,2)
        fr     = nfft2 * (w-anint(w)+0.5_ti_p)
        ifr    = int(fr-pme_eps)
        w      = fr - real(ifr,t_p)
        jgrd0  = ifr - bsorder
        call ibsplgen_2 (w,theta2,temp((threadIdx%x-1)*bsorder**2+1))
        w      = xi*recip_c(1,3) + yi*recip_c(2,3) + zi*recip_c(3,3)
        fr     = nfft3 * (w-anint(w)+0.5_ti_p)
        ifr    = int(fr-pme_eps)
        w      = fr - real(ifr,t_p)
        kgrd0  = ifr - bsorder
        call ibsplgen_2 (w,theta3,temp((threadIdx%x-1)*bsorder**2+1))
#else
        igrd0  = igrid(1,iglob)
        jgrd0  = igrid(2,iglob)
        kgrd0  = igrid(3,iglob)
#endif
        fi     = f * attrb(iichg)
        de1    = 0.0_ti_p
        de2    = 0.0_ti_p
        de3    = 0.0_ti_p
        k0     = kgrd0
        do it3 = 1, bsorder
           k0  = k0 + 1
           k   = k0 + 1 + (nfft3-sign(nfft3,k0))/2
#if 1
           t3  =       theta3(1+(it3-1)*level1)
           dt3 = dn3 * theta3(2+(it3-1)*level1)
#else
           t3  =       thetai3(1,it3,isite)
           dt3 = dn3 * thetai3(2,it3,isite)
#endif
           j0  = jgrd0
           do it2 = 1, bsorder
              j0  = j0 + 1
              j   = j0 + 1 + (nfft2-sign(nfft2,j0))/2
#if 1
              t2  =       theta2(1+(it2-1)*level1)
              dt2 = dn2 * theta2(2+(it2-1)*level1)
#else
              t2  =       thetai2(1,it2,isite)
              dt2 = dn2 * thetai2(2,it2,isite)
#endif
              i0  = igrd0
              do it1 = 1, bsorder
                 i0  = i0 + 1
                 i   = i0 + 1 + (nfft1-sign(nfft1,i0))/2
#if 1
                 t1  =       theta1(1+(it1-1)*level1)
                 dt1 = dn1 * theta1(2+(it1-1)*level1)
#else
                 t1  =       thetai1(1,it1,isite)
                 dt1 = dn1 * thetai1(2,it1,isite)
#endif
                 term = qgridin_t(1,i-istat+1,j-jstat+1
     $                             ,k-kstat+1,1)
c
                 de1 = de1 + term*dt1*t2*t3
                 de2 = de2 + term*dt2*t1*t3
                 de3 = de3 + term*dt3*t1*t2
              end do
           end do
        end do
        derec(1,iglob)=derec(1,iglob)+fi*scal*(recip_c(1,1)*de1
     &                    +recip_c(1,2)*de2+recip_c(1,3)*de3)
        derec(2,iglob)=derec(2,iglob)+fi*scal*(recip_c(2,1)*de1
     &                    +recip_c(2,2)*de2+recip_c(2,3)*de3)
        derec(3,iglob)=derec(3,iglob)+fi*scal*(recip_c(3,1)*de1
     &                    +recip_c(3,2)*de2+recip_c(3,3)*de3)
      end do
      end subroutine


      ! Scalar sum for reciprocal space (charges)
      attributes(global) subroutine pme_conv_kcu
     &                   (bsmod1,bsmod2,bsmod3
     &                   ,kst2,jst2,ist2,qsz1,qsz2,qsz12,qsz
     &                   ,nff,nf1,nf2,nf3,nfft1,nfft2,nfft3
     &                   ,f,pterm,volterm,xbox,calc_e,use_bounds
     &                   ,qgrid,e_buff,v_buff)
      implicit none
      integer  ,value:: nf1,nf2,nf3,nff,nfft1,nfft2,nfft3
     &         ,kst2,jst2,ist2,qsz1,qsz2,qsz12,qsz
      logical  ,value:: calc_e,use_bounds
      real(t_p),value:: pterm,f,xbox,volterm
      real(t_p),device:: bsmod1(*),bsmod2(*),bsmod3(*)
     &         ,qgrid(2*qsz)
      real(r_p),device:: e_buff(  RED_BUFF_SIZE)
      real(t_p),device:: v_buff(6*RED_BUFF_SIZE)

      integer   ithread,nthread,it
     &         ,k,k1,k2,k3,m1,m2,m3
      real(t_p) vterm,term,e0,stat
     &         ,expterm
     &         ,denom,hsq,struc2
     &         ,h1,h2,h3,r1,r2,r3
      
      ithread = (blockIdx%x-1)*blockDim%x + threadIdx%x
           it = iand(ithread-1,RED_BUFF_SIZE-1) + 1
c
c     use scalar sum to get reciprocal space energy and virial
c
      do k  = ithread-1,qsz-1, blockDim%x*gridDim%x
         k1 = mod(k,qsz1)      + ist2
         k2 = mod(k/qsz1,qsz2) + jst2
         k3 = k/qsz12          + kst2
         m1 = k1 - 1 - merge(nfft1,0,(k1.gt.nf1))
         m2 = k2 - 1 - merge(nfft2,0,(k2.gt.nf2))
         m3 = k3 - 1 - merge(nfft3,0,(k3.gt.nf3))
         if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
         r1 = real(m1,t_p)
         r2 = real(m2,t_p)
         r3 = real(m3,t_p)
         h1 = recip_c(1,1)*r1 + recip_c(1,2)*r2 + recip_c(1,3)*r3
         h2 = recip_c(2,1)*r1 + recip_c(2,2)*r2 + recip_c(2,3)*r3
         h3 = recip_c(3,1)*r1 + recip_c(3,2)*r2 + recip_c(3,3)*r3
         hsq     = h1*h1 + h2*h2 + h3*h3
         term    = -pterm * hsq
         expterm = 0.0_ti_p
         if (term.gt.-50.0_ti_p) then
            denom   = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2).ne.0)  expterm = 0.0_ti_p
            end if
            struc2 = qgrid(1+2*k)**2 + qgrid(2+2*k)**2
            if (calc_e.or.use_virial) then; block;
            real(r_p) rstat
              e0   = f * expterm * struc2
            vterm  = (2.0_ti_p/hsq) * (1.0_ti_p-term) * e0
            rstat= atomicAdd( e_buff(it),real(e0,r_p))
            stat = atomicAdd( v_buff(0*RED_BUFF_SIZE+it), h1*h1*vterm
     &                      - e0 )
            stat = atomicAdd( v_buff(1*RED_BUFF_SIZE+it), h1*h2*vterm )
            stat = atomicAdd( v_buff(2*RED_BUFF_SIZE+it), h1*h3*vterm )
            stat = atomicAdd( v_buff(3*RED_BUFF_SIZE+it), h2*h2*vterm
     &                      - e0 )
            stat = atomicAdd( v_buff(4*RED_BUFF_SIZE+it), h3*h2*vterm )
            stat = atomicAdd( v_buff(5*RED_BUFF_SIZE+it), h3*h3*vterm
     &                      - e0 )
            end block; end if
         end if
         qgrid(1+2*k) = expterm * qgrid(1+2*k)
         qgrid(2+2*k) = expterm * qgrid(2+2*k)
 10      continue
      end do

      !account for zeroth grid point for nonperiodic system
      if (.not.use_bounds .and. ithread.eq.1
     &    .and.ist2.eq.1.and.jst2.eq.1.and.kst2.eq.1) then
         expterm  = 0.5_ti_p * pi / xbox
         struc2   = qgrid(1)**2 + qgrid(2)**2
          if (calc_e) then
             e0   = f*expterm*struc2
           stat   = atomicAdd( e_buff(it),real(e0,r_p) )
          end if
         qgrid(1) = expterm*qgrid(1)
         qgrid(2) = expterm*qgrid(2)
      end if
      end subroutine

      ! Scalar sum for reciprocal space (dispersion)
      attributes(global) subroutine pme_convd_kcu
     &          (bsmod1,bsmod2,bsmod3,kst2,jst2,ist2,qsz1,qsz2,qsz12,qsz
     &          ,nff,nf1,nf2,nf3,nfft1,nfft2,nfft3,bfac,fac1,fac2,fac3
     &          ,denom0,pterm,xbox,calc_e,use_bounds
     &          ,qgrid,e_buff,v_buff)
      implicit none
      integer  ,value:: nf1,nf2,nf3,nff,nfft1,nfft2,nfft3
     &         ,kst2,jst2,ist2,qsz1,qsz2,qsz12,qsz
      logical  ,value:: calc_e,use_bounds
      real(t_p),value:: pterm,f,xbox,volterm,bfac,fac1,fac2,fac3,denom0
      real(t_p),device:: bsmod1(*),bsmod2(*),bsmod3(*)
     &         ,qgrid(2*qsz)
      real(r_p),device:: e_buff(  RED_BUFF_SIZE)
      real(t_p),device:: v_buff(6*RED_BUFF_SIZE)

      integer   ithread,nthread,it
     &         ,k,k1,k2,k3,m1,m2,m3
      real(t_p) vterm,term,e0,stat
     &         ,expterm,erfcterm
     &         ,denom,hsq,h,b,hhh,struc2
     &         ,h1,h2,h3,r1,r2,r3

      ithread = (blockIdx%x-1)*blockDim%x + threadIdx%x
           it = iand(ithread-1,RED_BUFF_SIZE-1) + 1

      if ((ist2.eq.1).and.(jst2.eq.1).and.(kst2.eq.1).and.
     &    ithread.eq.1) then
         qgrid(1) = 0.0
         qgrid(2) = 0.0
      end if
c
c     use scalar sum to get reciprocal space energy and virial
c
      do k  = ithread-1,qsz-1, blockDim%x*gridDim%x
         k1   = mod(k,qsz1)      + ist2
         k2   = mod(k/qsz1,qsz2) + jst2
         k3   = k/qsz12          + kst2
         m1   = k1 - 1 - merge(nfft1,0,(k1.gt.nf1))
         m2   = k2 - 1 - merge(nfft2,0,(k2.gt.nf2))
         m3   = k3 - 1 - merge(nfft3,0,(k3.gt.nf3))
         if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
         r1   = real(m1,t_p)
         r2   = real(m2,t_p)
         r3   = real(m3,t_p)
         h1   = recip_c(1,1)*r1 + recip_c(1,2)*r2 + recip_c(1,3)*r3
         h2   = recip_c(2,1)*r1 + recip_c(2,2)*r2 + recip_c(2,3)*r3
         h3   = recip_c(3,1)*r1 + recip_c(3,2)*r2 + recip_c(3,3)*r3
         h    = sqrt(hsq)
         b    = h*bfac
         hhh  = h*hsq
         term = -b * b
         expterm    = 0.0_ti_p
        denom = denom0*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
         if (term.gt.-50.0_ti_p) then
             expterm = exp(term)
            erfcterm = erfc(b)
            if (.not. use_bounds) then
               expterm =  expterm * (1.0-cos(pi*xbox*h))
              erfcterm = erfcterm * (1.0-cos(pi*xbox*h))
            else if (octahedron) then
               if (mod(m1+m2+m3,2).ne.0) then
                  expterm = 0.0_ti_p
                 erfcterm = 0.0_ti_p
               end if
            end if
             term = fac1*erfcterm*hhh + expterm*(fac2 + fac3*hsq)
            struc2 = qgrid(1+2*k)**2 + qgrid(2+2*k)**2
            if (calc_e.or.use_virial) then; block;
            real(r_p) rstat
              e0 = -(term / denom) * struc2
            vterm= 3.0*(fac1*erfcterm*h + fac3*expterm)*struc2/denom
            rstat= atomicAdd( e_buff(it),real(e0,r_p))
            stat = atomicAdd( v_buff(0*RED_BUFF_SIZE+it), h1*h1*vterm
     &                      - e0 )
            stat = atomicAdd( v_buff(1*RED_BUFF_SIZE+it), h1*h2*vterm )
            stat = atomicAdd( v_buff(2*RED_BUFF_SIZE+it), h1*h3*vterm )
            stat = atomicAdd( v_buff(3*RED_BUFF_SIZE+it), h2*h2*vterm
     &                      - e0 )
            stat = atomicAdd( v_buff(4*RED_BUFF_SIZE+it), h3*h2*vterm )
            stat = atomicAdd( v_buff(5*RED_BUFF_SIZE+it), h3*h3*vterm
     &                      - e0 )
            end block; end if
         end if
         qgrid(1+2*k) = expterm * qgrid(1+2*k)
         qgrid(2+2*k) = expterm * qgrid(2+2*k)
 10      continue
      end do
      end subroutine

      end module
