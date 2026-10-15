!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ################################################################
!     ##                                                            ##
!     ##  program testgrad  --  derivative test; Cartesian version  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "testgrad" computes and compares the analytical and numerical
!     gradient vectors of the potential energy function with respect
!     to Cartesian coordinates
!
!
#include "tinker_macro.h"
program testgrad
   use config, only: raw_cli_t, capture_raw_cli
   use parallel, only: get_domain_comm
   use iounit, only: iout
   use mdevice, only: bind_device
   use mpi
   use status, only: failed, status_t
   implicit none
   integer ierr,world_rank
   type(raw_cli_t) raw_cli
   type(status_t) status_value
   call bind_device(status_value)
   if (failed(status_value)) stop 1
   call capture_raw_cli(raw_cli)
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,world_rank,ierr)
   iout = 6
   call testgrad_bis(raw_cli,status_value)
   if (failed(status_value)) then
      if (world_rank.eq.0) write(iout,*) trim(status_value%message)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      stop 1
   end if
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)

   contains
#include "convert.inc.f90"

subroutine testgrad_bis(raw_cli,status_value)
   use atomsMirror
   use bound, only: use_bounds
   use boxes, only: volbox
   use cell
   use config, only: config_t, max_raw_cli_args, parse_common_cli, raw_cli_length, raw_cli_t
   use deriv
   use domdec
   use parallel, only: get_domain_comm
   use energi
   use sizes
   use inform
   use initialization, only: initialize, testgrad_application
   use iounit
   use mdevice, only: channel_t, finalize_device
   use mpole
   use parallel, only: parallel_t, release, assign_direct
   use status, only: failed, status_t
   use structure, only: structure_t
   use tinheader ,only: ti_p,re_p
   use usage
   use virial, only: dedv, vir, virnum, use_virial
   use mpi
   implicit none
   type(raw_cli_t), intent(in) :: raw_cli
   type(status_t), intent(out) :: status_value
   integer i,j,next,iglob,ierr
   integer iloc,value_count,value_index
   real(r_p) etot,f,f0,eps,eps0,old,energy
   real(r_p), parameter :: virial_eps=1.0e-3_re_p
   real(r_p) eb0,ea0,eba0,eub0,eaa0,eopb0
   real(r_p) eopd0,eid0,eit0,et0,ept0,eat0,ebt0
   real(r_p) ett0,ev0,er0,edsp0,ect0,ec0,em0,ep0
   real(r_p) eg0,ex0
   real(r_p) virial_analytic,virial_numeric,virial_residual
   real(r_p) virial_relative_residual
   real(r_p) virial_analytic_tensor(3,3),virial_tensor_norm
   real(r_p) grad_diff_norm,grad_scale,grad_relative_error
   real(r_p) totnorm,ntotnorm,rms,nrms
   real(r_p), allocatable :: denorm(:)
   real(r_p), allocatable :: ndenorm(:)
   real(r_p), allocatable :: detot(:,:)
   real(r_p), allocatable :: ndetot(:,:)
   real(r_p), allocatable :: ndeb(:,:)
   real(r_p), allocatable :: ndea(:,:)
   real(r_p), allocatable :: ndeba(:,:)
   real(r_p), allocatable :: ndeub(:,:)
   real(r_p), allocatable :: ndeaa(:,:)
   real(r_p), allocatable :: ndeopb(:,:)
   real(r_p), allocatable :: ndeopd(:,:)
   real(r_p), allocatable :: ndeid(:,:)
   real(r_p), allocatable :: ndeit(:,:)
   real(r_p), allocatable :: ndet(:,:)
   real(r_p), allocatable :: ndept(:,:)
   real(r_p), allocatable :: ndeat(:,:)
   real(r_p), allocatable :: ndebt(:,:)
   real(r_p), allocatable :: ndett(:,:)
   real(r_p), allocatable :: ndev(:,:)
   real(r_p), allocatable :: nder(:,:)
   real(r_p), allocatable :: ndedsp(:,:)
   real(r_p), allocatable :: ndect(:,:)
   real(r_p), allocatable :: ndem(:,:)
   real(r_p), allocatable :: ndec(:,:)
   real(r_p), allocatable :: ndep(:,:)
   real(r_p), allocatable :: ndeg(:,:)
   real(r_p), allocatable :: ndex(:,:)
   logical exist,query
   logical doanalyt,donumer,dofull
   logical virnum_save,use_virial_save
   character*1 answer
   character*1 axis(3)
   character*240 record
   character*240 string
   character(raw_cli_length) :: values(max_raw_cli_args)
   character(raw_cli_length) :: xyz_path
   type(config_t) configuration
   type(structure_t) molecular_structure
   type(parallel_t) platform
   type(channel_t) reciprocal_work
   type(status_t) cleanup_status
   data axis  / 'X','Y','Z' /
   interface
      subroutine dedvcalc(delta_)
         use tinheader, only: re_p
         real(r_p), optional, intent(in) :: delta_
      end subroutine dedvcalc
   end interface
!
!     set up the structure and mechanics calculation
!
   call initialize(testgrad_application, raw_cli, configuration, molecular_structure, &
      &platform, reciprocal_work, status_value)
   if (failed(status_value)) then
      if (ranktot.eq.0) write(iout,*) trim(status_value%message)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      __TINKER_FATAL__
   end if
!
   dotstgrad = .true.
!
   call nblist(0)
   call allocstep
!
!     parse the gradient controls after Configuration consumes the structure
!     argument and any option/value pairs
!
   call parse_common_cli(raw_cli,xyz_path,values,value_count,status_value)
   if (failed(status_value)) then
      if (ranktot.eq.0) write(iout,*) trim(status_value%message)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      __TINKER_FATAL__
   end if
   value_index = 1
!
!     decide whether to do an analytical gradient calculation
!
   doanalyt = .true.
   if (value_index.gt.value_count) then
      if(ranktot.eq.0) write (iout,10)
10    format (/,' You Need To Enter if You want to Compute the Analytical Gradient Vector Y or N')
      call MPI_BARRIER(get_domain_comm(),ierr)
      __TINKER_FATAL__ 
   end if
   answer = values(value_index)
   value_index = value_index + 1
   call upcase (answer)
   if (answer .eq. 'N')  doanalyt = .false.
!
   !     decide whether to do a numerical gradient calculation
!
   donumer = .true.
   if (value_index.gt.value_count) then
      if (ranktot.eq.0) write (iout,30)
30    format (/,' You Need To Enter if You want to Compute the Numerical Gradient Vector: Y or N')
      call MPI_BARRIER(get_domain_comm(),ierr)
      __TINKER_FATAL__
   end if
   answer = values(value_index)
   value_index = value_index + 1
   call upcase (answer)
   if (answer .eq. 'N')  donumer = .false.
!
!     get the stepsize for numerical gradient calculation
!
   if (donumer) then
      eps = -1.0_ti_p
      eps0 = 1.0d-5
      query = .true.
      if (value_index.gt.value_count) then
         if (ranktot.eq.0) write (iout,60)  eps0
60       format (/,' You Need To Enter a Numerical Stepsize')
         __TINKER_FATAL__
      else
         read (values(value_index),*,err=50,end=50)  eps
         value_index = value_index + 1
      end if
50    continue
      if (eps .le. 0.0_ti_p)  eps = eps0
   end if
!
!     decide whether to output results by gradient component
!
   dofull = .true.
   if (n .gt. 100) then
      dofull = .false.
      if (value_index.gt.value_count) then
         if (ranktot.eq.0) write (iout,80)
80       format (/,' You Need To Enter If You Want The Output Breakdown by Gradient Component [Y|N]')
         call MPI_BARRIER(get_domain_comm(),ierr)
         __TINKER_FATAL__
      end if
      answer = values(value_index)
      value_index = value_index + 1
      call upcase (answer)
      if (answer .eq. 'Y')  dofull = .true.
   end if
!
!     perform dynamic allocation of some local arrays
!
   allocate (detot(3,nbloc))
   detot = 0_ti_p
!$acc enter data create(etot)
!$acc enter data copyin(detot)
   call assign_direct
   call reassignpme(.false.)
   call reinitnl(0)
   call mechanicstep(0)
   call nblist(0)
   !
   !     compute the analytical gradient components
   !
   if (doanalyt) then
      virnum_save = virnum
      use_virial_save = use_virial
      virnum = .true.
      use_virial = .true.
      call zero_forces
      call gradient (etot,detot)
!$acc wait
!$acc update self(etot)
!$acc update self(vir)
!$acc wait
!
!    MPI : communicate analytical forces and torques
!
      ! Diagnostic consumer; the exchange remains owned by deriv.
      call commforceststgrad(detot)
      call commforcesrec(detot)
      if (nproc.eq.1) then
         if (ftot_l) then
            call get_ftot(detot,nbloc)
         end if
!$acc update self(detot(1:3,1:nbloc)) async(reciprocal_work%queue)
!$acc wait(reciprocal_work%queue)
      end if
!
!    MPI : get total energies from gradient routine
!
      call reduceen(etot)
      if (use_bounds) then
         virial_analytic_tensor = vir
         virial_tensor_norm = sqrt(sum(virial_analytic_tensor**2))
         virial_analytic = (vir(1,1)+vir(2,2)+vir(3,3)) /&
            &(3.0_re_p*volbox)
      end if
      virnum = virnum_save
      use_virial = use_virial_save
   end if
!
!     print the total potential energy of the system
!
   if (doanalyt.and.rank.eq.0) then
      if (digits .ge. 8) then
         write (iout,100)  etot
100      format (/,' Total Potential Energy :',8x,f20.8,' Kcal/mole')
      else if (digits .ge. 6) then
         write (iout,110)  etot
110      format (/,' Total Potential Energy :',8x,f18.6,' Kcal/mole')
      else
         write (iout,120)  etot
120      format (/,' Total Potential Energy :',8x,f16.4,' Kcal/mole')
      end if
!
!     print the energy breakdown over individual components
!
      write (iout,130)
130   format (/,' Potential Energy Breakdown by Individual',&
         &' Components :')
!         if (digits .ge. 8) then
!            write (iout,140)
!  140       format (/,'  Energy',7x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
!     &              /,'  Terms',8x,'EAA',13x,'EOPB',12x,'EOPD',
!     &                 12x,'EID',
!     &              /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
!     &              /,15x,'ETT',13x,'EV',14x,'EM',14x,'EP')
!            write (iout,150)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
!     &                        eit,et,ept,ebt,ett,ev,em,ep
!  150       format (/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,
!     &                 /,6x,4f16.8,/,6x,4f16.8)
!         else if (digits .ge. 6) then
!            write (iout,160)
!  160       format (/,'  Energy',7x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
!     &              /,'  Terms',8x,'EAA',13x,'EOPB',12x,'EOPD',
!     &                 12x,'EID',
!     &              /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
!     &              /,15x,'ETT',13x,'EV',14x,'EM',14x,'EP')
!            write (iout,170)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
!     &                        eit,et,ept,ebt,ett,ev,em,ep
!  170       format (/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,
!     &                 /,6x,4f14.6)
!         else
      write (iout,180)
180   format (/,'  Energy',7x,'EB',14x,'EA',14x,'EBA',13x,'EUB',&
         &11x,'EX',&
         &/,'  Terms',8x,'EAA',13x,'EOPB',12x,'EOPD',&
         &12x,'EID',11x,'EC',&
         &/,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',&
         &11x,'EREP',/,15x,'ETT',13x,'EV',14x,'EM',14x,'EP',&
         &12x,'EXDISP',&
         &/,15x,'EAT')
      write (iout,190)  eb,ea,eba,eub,ex,eaa,eopb,eopd,eid,ec,&
         &eit,et,ept,ebt,enr2en(er),ett,ev,em,ep,enr2en(edsp),&
         &enr2en(ect),&
         &eat
190   format (/,6x,5f15.4,/,6x,5f15.4,/,6x,5f15.4,/,6x,5f15.4,&
         &/,6x,2f15.4)
!         end if
   end if
!
!     print a header for the gradients of individual potentials
!
   if (dofull.and.rank.eq.0) then
      write (iout,200)
200   format (/,' Cartesian Gradient Breakdown by Individual',&
         &' Components :')
      if (digits .ge. 8) then
         write (iout,210)
210      format (/,2x,'Atom',9x,'d EB',11x,'d EA',11x,'d EBA',&
            &10x,'d EUB',8x,'d EX',&
            &/,2x,'Axis',9x,'d EAA',9x,'d EOPB',10x,'d EOPD',&
            &9x,'d EID',8x,'d EC',&
            &/,2x,'Type',9x,'d EIT',10x,'d ET',11x,'d EPT',&
            &10x,'d EBT',8x,'d EREP',&
            &/,15x,'d ETT',10x,'d EV',11x,'d EM',11x,'d EP',9x,&
            &'d EDISP',/,15x,'d ECT',10x,'d EAT')
      else if (digits .ge. 6) then
         write (iout,220)
220      format (/,2x,'Atom',9x,'d EB',11x,'d EA',11x,'d EBA',&
            &10x,'d EUB',8x,' d EX',&
            &/,2x,'Axis',9x,'d EAA',9x,'d EOPB',10x,'d EOPD',&
            &9x,'d EID',8x,'d EC',&
            &/,2x,'Type',9x,'d EIT',10x,'d ET',11x,'d EPT',&
            &10x,'d EBT',8x,'d ER'&
            &/,15x,'d ETT',10x,'d EV',11x,'d EM',11x,'d EP',9x,&
            &'d EDISP',/,15x,'d ECT',10x,'d EAT')
      else
         write (iout,230)
230      format (/,2x,'Atom',9x,'d EB',11x,'d EA',11x,'d EBA',&
            &10x,'d EUB',8x,'d EX',&
            &/,2x,'Axis',9x,'d EAA',10x,'d EOPB',9x,'d EOPD',&
            &9x,'d EID',8x,'d EC',&
            &/,2x,'Type',9x,'d EIT',10x,'d ET',11x,'d EPT',&
            &10x,'d EBT',8x,'d ER',&
            &/,15x,'d ETT',10x,'d EV',11x,'d EM',11x,'d EP',9x,&
            &'d EDISP',/,15x,'d ECT',10x,'d EAT')
      end if
   end if
!
!     perform dynamic allocation of some local arrays
!
   allocate (ndetot(3,nloc))
   allocate (ndeb(3,nloc))
   allocate (ndea(3,nloc))
   allocate (ndeba(3,nloc))
   allocate (ndeub(3,nloc))
   allocate (ndeaa(3,nloc))
   allocate (ndeopb(3,nloc))
   allocate (ndeopd(3,nloc))
   allocate (ndeid(3,nloc))
   allocate (ndeit(3,nloc))
   allocate (ndet(3,nloc))
   allocate (ndept(3,nloc))
   allocate (ndeat(3,nloc))
   allocate (ndebt(3,nloc))
   allocate (ndett(3,nloc))
   allocate (ndev(3,nloc))
   allocate (nder(3,nloc))
   allocate (ndedsp(3,nloc))
   allocate (ndect(3,nloc))
   allocate (ndem(3,nloc))
   allocate (ndec(3,nloc))
   allocate (ndep(3,nloc))
   allocate (ndeg(3,nloc))
   allocate (ndex(3,nloc))
   ndetot = 0_re_p
   ndeb = 0_re_p
   ndea = 0_re_p
   ndeba = 0_re_p
   ndeub = 0_re_p
   ndeaa = 0_re_p
   ndeopb = 0_re_p
   ndeopd = 0_re_p
   ndeid = 0_re_p
   ndet = 0_re_p
   ndept = 0_re_p
   ndebt = 0_re_p
   ndett = 0_re_p
   ndev = 0_re_p
   ndem = 0_re_p
   ndep = 0_re_p
   ndeg = 0_re_p
   ndex = 0_re_p
   ndeg = 0_re_p
   ndex = 0_re_p
   ndep = 0.0_re_p
   ndeg = 0.0_re_p
   ndex = 0.0_re_p
!
!     get the Cartesian component two-sided numerical gradients
!
   do i = 1, n
      if (repart(i).eq.rank) iloc = loc(i)
      call MPI_BARRIER(get_domain_comm(),ierr)
      if (donumer .and. use(i)) then
         do j = 1, 3
            if (j .eq. 1) then
               old = x(i)
               x(i) = x(i) - 0.5_re_p*eps
            else if (j .eq. 2) then
               old = y(i)
               y(i) = y(i) - 0.5_re_p*eps
            else if (j .eq. 3) then
               old = z(i)
               z(i) = z(i) - 0.5_re_p*eps
            end if
!$acc update device(x(i),y(i),z(i))
            call reCast_position
!$acc wait
            f0 = energy ()
            call allreduceen(f0)
            eb0 = eb
            ea0 = ea
            eba0 = eba
            eub0 = eub
            eaa0 = eaa
            eopb0 = eopb
            eopd0 = eopd
            eid0 = eid
            eit0 = eit
            et0 = et
            ept0 = ept
            eat0 = eat
            ebt0 = ebt
            ett0 = ett
            ev0 = ev
            er0 = er
            edsp0 = edsp
            ect0 = ect
            ec0 = ec
            em0 = em
            ep0 = ep
            eg0 = eg
            ex0 = ex
            eg0 = eg
            if (j .eq. 1) then
               x(i) = x(i) + eps
            else if (j .eq. 2) then
               y(i) = y(i) + eps
            else if (j .eq. 3) then
               z(i) = z(i) + eps
            end if
!$acc update device(x(i),y(i),z(i))
            call reCast_position
!$acc wait
            f = energy ()
            call allreduceen(f)
            if (j .eq. 1) then
               x(i) = old
            else if (j .eq. 2) then
               y(i) = old
            else if (j .eq. 3) then
               z(i) = old
            end if
!$acc update device(x(i),y(i),z(i))
            call reCast_position
!$acc wait
!
!   si i dans domenloc alors nde et affichage
!
            if (repart(i).eq.rank) then
               ndetot(j,iloc) = (f - f0) / eps
               ndeb(j,iloc) = (eb - eb0) / eps
               ndea(j,iloc) = (ea - ea0) / eps
               ndeba(j,iloc) = (eba - eba0) / eps
               ndeub(j,iloc) = (eub - eub0) / eps
               ndeaa(j,iloc) = (eaa - eaa0) / eps
               ndeopb(j,iloc) = (eopb - eopb0) / eps
               ndeopd(j,iloc) = (eopd - eopd0) / eps
               ndeid(j,iloc) = (eid - eid0) / eps
               ndeit(j,iloc) = (eit - eit0) / eps
               ndet(j,iloc) = (et - et0) / eps
               ndept(j,iloc) = (ept - ept0) / eps
               ndeat(j,iloc) = (eat - eat0) / eps
               ndebt(j,iloc) = (ebt - ebt0) / eps
               ndett(j,iloc) = (ett - ett0) / eps
               ndev(j,iloc) = (ev - ev0) / eps
               nder(j,iloc) = (er - er0) / eps
               ndedsp(j,iloc) = (edsp - edsp0) / eps
               ndect(j,iloc) = (ect - ect0) / eps
               ndem(j,iloc) = (em - em0) / eps
               ndec(j,iloc) = (ec - ec0) / eps
               ndep(j,iloc) = (ep - ep0) / eps
               ndeg(j,iloc) = (eg - eg0) / eps
               ndex(j,iloc) = (ex - ex0) / eps
            end if
         end do
      end if
!
!     print analytical gradients of each energy term for each atom
!
      if (dofull .and. use(i)) then
         do j = 1, 3
            if (doanalyt.and.repart(i).eq.rank) then
!                  if (digits .ge. 8) then
!                     write (iout,240) iglob,deb(j,i),dea(j,i),deba(j,i),
!     &                                deub(j,i),axis(j),deaa(j,i),
!     &                                deopb(j,i),deopd(j,i),deid(j,i),
!     &                                deit(j,i),det(j,i),dept(j,i),
!     &                                debt(j,i),dett(j,i),dev(j,i),
!     &                                dem(j,i),dep(j,i)
!  240                format (/,i6,4f16.8,/,5x,a1,4f16.8,
!     &                          /,' Anlyt',4f16.8,/,6x,4f16.8,
!     &                          /,6x,4f16.8,/,6x,4f16.8)
!                  else if (digits .ge. 6) then
!                     write (iout,250) iglob,deb(j,i),dea(j,i),deba(j,i),
!     &                                deub(j,i),deaa(j,i),axis(j),
!     &                                deopb(j,i),deopd(j,i),deid(j,i),
!     &                                deit(j,i),det(j,i),dept(j,i),
!     &                                debt(j,i),dett(j,i),dev(j,i),
!     &                                dem(j,i),dep(j,i)
!  250                format (/,i6,5f14.6,/,5x,a1,5f14.6,/,' Anlyt',
!     &                          5f14.6,/,6x,5f14.6,/,6x,4f14.6)
!                  else
               if (associated(deb).and.associated(dea).and.&
                  &associated(deba).and.associated(deub).and.&
                  &associated(dex).and.associated(deaa).and.&
                  &associated(deopb).and.associated(deopd).and.&
                  &associated(deid).and.associated(dec).and.&
                  &associated(deit).and.associated(det).and.&
                  &associated(dept).and.associated(debt).and.&
                  &associated(der).and.associated(dett).and.&
                  &associated(dev).and.associated(dem).and.&
                  &associated(dep).and.associated(dedsp).and.&
                  &associated(dect).and.associated(deat)) then
                  write (iout,260) i,deb(j,iloc),dea(j,iloc),&
                     &deba(j,iloc),deub(j,iloc),dex(j,iloc),&
                     &axis(j),&
                     &deaa(j,iloc),deopb(j,iloc),deopd(j,iloc),&
                     &deid(j,iloc),dec(j,iloc),deit(j,iloc),&
                     &det(j,iloc),&
                     &dept(j,iloc),debt(j,iloc),der(j,iloc),&
                     &dett(j,iloc),dev(j,iloc),dem(j,iloc),&
                     &dep(j,iloc),dedsp(j,iloc),dect(j,iloc),&
                     &deat(j,iloc)
               end if
260            format (/,i6,5f15.4,/,5x,a1,5f15.4,/,' Anlyt',&
                  &5f15.4,/,6x,5f15.4,/,6x,2f15.4)
!                  end if
            end if
!
!     print numerical gradients of each energy term for each atom
!
            if (donumer.and.repart(i).eq.rank) then
!                  if (digits .ge. 8) then
!                     write (iout,270) iglob,ndeb(j,i),ndea(j,i),
!     &                                ndeba(j,i),ndeub(j,i),ndeaa(j,i),
!     &                                ndeopb(j,i),axis(j),
!     &                                ndeaa(j,i),ndeopb(j,i),
!     &                                ndeopd(j,i),ndeid(j,i),ndeit(j,i), c     &                                ndet(j,i),ndept(j,i),ndebt(j,i),
!     &                                ndett(j,i),ndev(j,i),
!     &                                ndem(j,i),ndep(j,i)
!  270                format (/,i6,4f16.8,/,5x,a1,4f16.8,
!     &                          /,' Numer',4f16.8,/,6x,4f16.8,
!     &                          /,6x,4f16.8,/,6x,4f16.8)
!                  else if (digits .ge. 6) then
!                     write (iout,280) iglob,ndeb(j,i),ndea(j,i),
!     &                                ndeba(j,i),ndeub(j,i),ndeaa(j,i),
!     &                                ndeopb(j,i),axis(j),
!     &                                ndeaa(j,i),ndeopb(j,i),
!     &                                ndeopd(j,i),ndeid(j,i),ndeit(j,i),
!     &                                ndet(j,i),ndept(j,i),ndebt(j,i),
!     &                                ndett(j,i),ndev(j,i),
!     &                                ndem(j,i),ndep(j,i)
!  280                format (/,i6,5f14.6,/,5x,a1,5f14.6,/,' Numer',
!     &                          5f14.6,/,6x,5f14.6,/,6x,4f14.6)
!                  else
               write (iout,290) i,ndeb(j,iloc),ndea(j,iloc),&
                  &ndeba(j,iloc),ndeub(j,iloc),&
                  &ndex(j,iloc),axis(j),&
                  &ndeaa(j,iloc),ndeopb(j,iloc),&
                  &ndeopd(j,iloc),ndeid(j,iloc),&
                  &ndec(j,iloc),ndeit(j,iloc),&
                  &ndet(j,iloc),&
                  &ndept(j,iloc),ndebt(j,iloc),&
                  &nder(j,iloc),ndett(j,iloc),&
                  &ndev(j,iloc),ndem(j,iloc),ndep(j,iloc),&
                  &ndedsp(j,iloc),ndect(j,iloc),&
                  &ndeat(j,iloc)

290            format (/,i6,5f15.4,/,5x,a1,5f15.4,/,' Numer',&
                  &5f15.4,/,6x,5f15.4,/,6x,2f15.4)
!                  end if
            end if
         end do
      end if
   end do
   call MPI_BARRIER(get_domain_comm(),ierr)
!
!     perform deallocation of some local arrays
!
   deallocate (ndeb)
   deallocate (ndea)
   deallocate (ndeba)
   deallocate (ndeub)
   deallocate (ndeaa)
   deallocate (ndeopb)
   deallocate (ndeopd)
   deallocate (ndeid)
   deallocate (ndeit)
   deallocate (ndet)
   deallocate (ndept)
   deallocate (ndeat)
   deallocate (ndebt)
   deallocate (ndett)
   deallocate (ndev)
   deallocate (nder)
   deallocate (ndedsp)
   deallocate (ndect)
   deallocate (ndem)
   deallocate (ndec)
   deallocate (ndep)
   deallocate (ndeg)
   deallocate (ndex)
!
!
!     perform dynamic allocation of some local arrays
!
   allocate (denorm(nloc))
   allocate (ndenorm(nloc))
   denorm = 0_re_p
   ndenorm = 0_re_p
!
   call MPI_BARRIER(get_domain_comm(),ierr)
!
!     print the total gradient components for each atom
!
   if (rank.eq.0) then
      write (iout,300)
300   format (/,' Cartesian Gradient Breakdown over Individual Atoms :')
      if (digits .ge. 8) then
         write (iout,310)
310      format (/,2x,'Type',4x,'Atom',10x,'dE/dX',11x,'dE/dY',&
            &11x,'dE/dZ',11x,'Norm',/)
      else if (digits .ge. 6) then
         write (iout,320)
320      format (/,2x,'Type',6x,'Atom',11x,'dE/dX',9x,'dE/dY',&
            &9x,'dE/dZ',11x,'Norm',/)
      else
         write (iout,330)
330      format (/,2x,'Type',6x,'Atom',14x,'dE/dX',7x,'dE/dY',&
            &7x,'dE/dZ',10x,'Norm',/)
      end if
   end if
   totnorm = 0.0_re_p
   ntotnorm = 0.0_re_p
   grad_diff_norm = 0.0_re_p
   grad_scale = 0.0_re_p
   do i = 1, nloc
      iglob = glob(i)
      if (doanalyt .and. use(iglob)) then
         denorm(i) = detot(1,i)**2 + detot(2,i)**2 +&
            &detot(3,i)**2
         totnorm = totnorm + denorm(i)
         denorm(i) = sqrt(denorm(i))
!            if (digits .ge. 8) then
!               write (iout,340)  iglob,(detot(j,i),j=1,3),denorm(i)
!  340          format (' Anlyt',i8,1x,3f16.8,f16.8)
!            else if (digits .ge. 6) then
!               write (iout,350)  iglob,(detot(j,i),j=1,3),denorm(i)
!  350          format (' Anlyt',2x,i8,3x,3f14.6,2x,f14.6)
!            else
         write (iout,360)  iglob,(detot(j,i),j=1,3),denorm(i)
360      format (' Anlyt',2x,i8,7x,3f12.4,2x,f12.4)
!            end if
      end if
      if (donumer .and. use(iglob)) then
         ndenorm(i) = ndetot(1,i)**2 + ndetot(2,i)**2 +&
            &ndetot(3,i)**2
         ntotnorm = ntotnorm + ndenorm(i)
         if (doanalyt) then
            grad_diff_norm = grad_diff_norm + sum((detot(:,i)-ndetot(:,i))**2)
            grad_scale = grad_scale + max(sum(detot(:,i)**2),sum(ndetot(:,i)**2))
         end if
         ndenorm(i) = sqrt(ndenorm(i))
!            if (digits .ge. 8) then
!               write (iout,370)  iglob,(ndetot(j,i),j=1,3),ndenorm(i)
!  370          format (' Numer',i8,1x,3f16.8,f16.8)
!            else if (digits .ge. 6) then
!               write (iout,380)  iglob,(ndetot(j,i),j=1,3),ndenorm(i)
!  380          format (' Numer',2x,i8,3x,3f14.6,2x,f14.6)
!            else
         write (iout,390)  iglob,(ndetot(j,i),j=1,3),ndenorm(i)
390      format (' Numer',2x,i8,7x,3f12.4,2x,f12.4)
!            end if
      end if
   end do
   call MPI_BARRIER(get_domain_comm(),ierr)
   if (rank.eq.0) then
      call MPI_REDUCE(MPI_IN_PLACE,totnorm,1,MPI_RPREC,MPI_SUM,0,&
         &get_domain_comm(),ierr)
   else
      call MPI_REDUCE(totnorm,totnorm,1,MPI_RPREC,MPI_SUM,0,&
         &get_domain_comm(),ierr)
   end if
   if (rank.eq.0) then
      call MPI_REDUCE(MPI_IN_PLACE,ntotnorm,1,MPI_RPREC,MPI_SUM,&
         &0,get_domain_comm(),ierr)
   else
      call MPI_REDUCE(ntotnorm,ntotnorm,1,MPI_RPREC,MPI_SUM,0,&
         &get_domain_comm(),ierr)
   end if
   if (doanalyt.and.donumer) then
      if (rank.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,grad_diff_norm,1,MPI_RPREC,MPI_SUM,0,&
            &get_domain_comm(),ierr)
         call MPI_REDUCE(MPI_IN_PLACE,grad_scale,1,MPI_RPREC,MPI_SUM,0,&
            &get_domain_comm(),ierr)
      else
         call MPI_REDUCE(grad_diff_norm,grad_diff_norm,1,MPI_RPREC,MPI_SUM,0,&
            &get_domain_comm(),ierr)
         call MPI_REDUCE(grad_scale,grad_scale,1,MPI_RPREC,MPI_SUM,0,&
            &get_domain_comm(),ierr)
      end if
   end if
!
!     perform deallocation of some local arrays
!
!$acc exit data delete(etot)
!$acc exit data delete(detot)
   deallocate (detot)
   deallocate (ndetot)
   deallocate (denorm)
   deallocate (ndenorm)
!
!
!     print the total norm for the analytical gradient
!
   if (rank.eq.0) then
      write (iout,400)
400   format (/,' Total Gradient Norm and RMS Gradient per Atom :')
      if (doanalyt) then
         totnorm = sqrt(totnorm)
!           if (digits .ge. 8) then
!              write (iout,410)  totnorm
!  410         format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f20.8)
!           else if (digits .ge. 6) then
!              write (iout,420)  totnorm
!  420         format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f18.6)
!           else
         write (iout,430)  totnorm
430      format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f16.4)
!           end if
      end if
   end if
!
!     print the total norm for the numerical gradient
!
   if (donumer) then
      ntotnorm = sqrt(ntotnorm)
!         if (digits .ge. 8) then
!            write (iout,440)  ntotnorm
!  440       format (' Numer',6x,'Total Gradient Norm Value',6x,f20.8)
!         else if (digits .ge. 6) then
!            write (iout,450)  ntotnorm
!  450       format (' Numer',6x,'Total Gradient Norm Value',6x,f18.6)
!         else
      if (rank.eq.0) write (iout,460)  ntotnorm
460   format (' Numer',6x,'Total Gradient Norm Value',6x,f16.4)
!         end if
   end if
!
!     print the rms per atom norm for the analytical gradient
!
   if (doanalyt) then
      rms = totnorm / sqrt(real(nuse,t_p))
!         if (digits .ge. 8) then
!            write (iout,470)  rms
!  470       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
!     &                 4x,f20.8)
!         else if (digits .ge. 6) then
!            write (iout,480)  rms
!  480       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
!     &                 4x,f18.6)
!         else
      if (rank.eq.0) write (iout,490)  rms
490   format (/,' Anlyt',6x,'RMS Gradient over All Atoms',&
         &4x,f16.4)
!         end if
   end if
!
!     print the rms per atom norm for the numerical gradient
!
   if (donumer) then
      nrms = ntotnorm / sqrt(real(nuse,t_p))
!         if (digits .ge. 8) then
!            write (iout,500)  nrms
!  500       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f20.8)
!         else if (digits .ge. 6) then
!            write (iout,510)  nrms
!  510       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f18.6)
!         else
      if (rank.eq.0) write (iout,520)  nrms
520   format (' Numer',6x,'RMS Gradient over All Atoms',4x,f16.4)
!         end if
   end if
   if (doanalyt.and.donumer.and.rank.eq.0) then
      grad_relative_error = sqrt(grad_diff_norm) /&
         &sqrt(max(grad_scale,tiny(1.0_re_p)))
      write (iout,525) sqrt(grad_diff_norm),grad_relative_error
525   format (/,' Gradient L2 Difference :',8x,f16.8,/,&
         &6x,'Relative Difference',6x,f16.8)
   end if
!
!     perform any final tasks before program exit
!
   if (doanalyt) then
      if (use_bounds) then
         virnum = .true.
         use_virial = .true.
         call dedvcalc(virial_eps)
         virial_numeric = dedv
         virial_residual = virial_analytic - virial_numeric
         virial_relative_residual = abs(virial_residual) /&
            &max(max(abs(virial_analytic),abs(virial_numeric)),tiny(1.0_re_p))
      end if
      virnum = virnum_save
      use_virial = use_virial_save
   end if
   if (doanalyt.and.use_bounds.and.rank.eq.0) then
      write (iout,125) virial_analytic,virial_numeric,virial_residual,&
         &virial_relative_residual
125   format (/,' Isotropic Virial dE/dV Check :',/&
         &6x,'Anlyt',6x,f16.8,/,6x,'Numer',6x,f16.8,/,6x,&
         &'Residual',3x,f16.8,/,6x,'Relative Residual',f16.8)
      write (iout,126) ((virial_analytic_tensor(i,j),j=1,3),i=1,3),&
         &virial_tensor_norm
126   format (/,' Analytical Virial Tensor :',/,&
         &6x,'XX',6x,f16.8,6x,'XY',6x,f16.8,6x,'XZ',6x,f16.8,/,&
         &6x,'YX',6x,f16.8,6x,'YY',6x,f16.8,6x,'YZ',6x,f16.8,/,&
         &6x,'ZX',6x,f16.8,6x,'ZY',6x,f16.8,6x,'ZZ',6x,f16.8,/,&
         &6x,'Frobenius Norm',2x,f16.8)
   end if
   call final(status_value)
   call finalize_device(cleanup_status)
   call release(platform, cleanup_status)
   if (failed(cleanup_status)) status_value = cleanup_status
end subroutine testgrad_bis

end program testgrad
