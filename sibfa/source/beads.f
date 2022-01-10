c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module beads   --  pimd variables                              ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     nbeads  number of global replicas used for PIMD simulations
c     nbeadsloc  number of local replicas used for PIMD simulations
c     nproctot number of process for beads parallelism
c     nproc number of process for gradient parallelism
c     ncomm number of communicators for beads parallelism
c     globbead : array to switch from local to global beads
c     locbead : array to switch from global to local beads
c     ibeadsloc : number of the current local beads
c     
c  dynamical variables replicated for PIMD simulations
c
c     pospi : array of positions
c     velpi : array of velocities
c     api : array of accelerations
c     masspi : array of PIMD masses (may differ for original masses, used for sampling)
c     
c     array used for domain decomposition within a bead:
c      glob
c      loc
c      repart
c      repartrec
c      domlen
c      domlenrec number of reciprocal atoms in the reciprocal domains
c      domlenpole
c      domlenpolerec
c      globrec
c      locrec
c      globrec1
c      locrec1
c      bufbegrec
c      bufbegpole
c      bufbeg
c     buflen1,buflen2,buf1,buf2,bufbeg1,bufbeg2 explicit direct-reciprocal atomic correspondance, 
c      nloc
c      nbloc
c      nlocrec
c      nblocrec local + neighbors reciprocal number of atoms
c      nlocnl local nl number of atoms
c      nblocrecdir local + neighbors direct+reciprocal number of atoms
c
c      molecule 
c      nmolelocpi,molculeglobpi
c
c      VDW :
c        nvdwblocpi,vdwglobpi,nvdwlocnlpi,vdwglobnlpi
c        nvlstpi,vlstpi
c
c      BONDS:
c        nbondlocpi,bndglobpi
c
c      STRETCH-BEND:
c        nstrbndlocpi,strbndglobpi
c
c      ANGLE-ANGLE:
c        nanganglocpi,angangglobpi
c
c      OP-BENDING:
c        nopbendlocpi,opbendglobpi
c
c      OP-DIST:
c        nopdistlocpi,opdistglobpi
c
c      IMPROP:
c        niproplocpi,impropglobpi
c
c      IMPTOR:
c        nitorslocpi,imptorglobpi 
c
c      TORSION:
c        ntorslocpi,torsglobpi 
c
c      PITORSION:
c        npitorslocpi,pitorsglobpi 
c
c      STRETCH-TORSION:
c        nstrtorlocpi,strtorglobpi 
c
c      TORSION-TORSION:
c        ntortorlocpi,tortorglobpi 
c
c      ANGLE:
c        nanglelocpi,angleglobpi 
c
c      CHARGE:
c        nionreclocpi,chgrecglobpi 
c        nionlocpi,chgglobpi
c        nionlocnlpi,chgglobnlpi
c        nelstpi,elstpi
c
c      MULTIPOLE:
c        npolereclocpi,polerecglobpi 
c        npolelocpi,poleglobpi,polelocpi
c        npolelocnlpi,poleglobnlpi
c        
c      POLARIZATION:
c        udaltpi,upaltpi,nualtpi
c        uindpi,uinppi
c        
c
c      STATS:
c        etot_sumpi,etot2_sumpi,eint_sumpi,eint2_sumpi
c        epot_sumpi,epot2_sumpi,ekin_sumpi,ekin2_sumpi
c        temp_sumpi,temp2_sumpi,pres_sumpi,pres2_sumpi
c        dens_sumpi,dens2_sumpi
c
c      TIME:
c        timestep
c
      module beads
      implicit none
      integer :: ibeadsloc
      integer :: nbeads,nbeadsloc,nprocbeads,ncomm
      integer, allocatable :: globbead(:),locbead(:)
      real*8 :: temppi,temppi_cl
      real*8, allocatable :: epotpi(:),etotpi(:)
      real*8, allocatable :: eksumpi(:),ekinpi(:,:,:)
      real*8 :: epotpi_loc,etotpi_loc
      real*8 :: eksumpi_loc,ekinpi_loc(3,3)
      real*8, allocatable :: dedv_pi(:)
      real*8 :: prespi
      
      real*8, allocatable :: eigmat(:,:)
      real*8, allocatable :: eigmattr(:,:)
      real*8, allocatable ::  omkpi(:)
c
c     parallelism
c
      integer, allocatable :: nlocpi(:),nblocpi(:),nlocrecpi(:)
      integer, allocatable :: nblocrecpi(:),nlocnlpi(:)
      integer, allocatable :: nblocrecdirpi(:)
      integer, allocatable :: globpi(:,:),locpi(:,:),repartpi(:,:)
      integer, allocatable :: ineignlpi(:,:),locnlpi(:,:) 
      integer, allocatable :: repartrecpi(:,:), domlenpi(:,:)
      integer, allocatable :: domlenrecpi(:,:),domlenpolepi(:,:)
      integer, allocatable :: domlenpolerecpi(:,:),globrecpi(:,:)
      integer, allocatable :: locrecpi(:,:),globrec1pi(:,:)
      integer, allocatable :: locrec1pi(:,:)
      integer, allocatable :: bufbegrecpi(:,:)
      integer, allocatable :: bufbegpolepi(:,:),bufbegpi(:,:)
      integer, allocatable :: buflen1pi(:,:),buflen2pi(:,:)
      integer, allocatable :: buf1pi(:,:),buf2pi(:,:)
      integer, allocatable :: bufbeg1pi(:,:),bufbeg2pi(:,:)
      integer, allocatable :: nmolelocpi(:),molculeglobpi(:,:)
      real*8, allocatable :: pospi(:,:,:),velpi(:,:,:),masspi(:)
      real*8, allocatable :: api(:,:,:)
      real*8, allocatable :: Ekcentroid
c
c     BOND-STRETCHING
c
      integer, allocatable :: nbondlocpi(:),bndglobpi(:,:)
c
c     STRETCH-BENDING
c
      integer, allocatable :: nstrbndlocpi(:),strbndglobpi(:,:)
c
c     UREY-BRADLEY
c
      integer, allocatable :: nureylocpi(:),ureyglobpi(:,:)
c
c     ANGLE-ANGLE
c
      integer, allocatable :: nanganglocpi(:),angangglobpi(:,:)
c
c     OP-BENDING
c
      integer, allocatable :: nopbendlocpi(:),opbendglobpi(:,:)
c
c     OP-DIST
c
      integer, allocatable :: nopdistlocpi(:),opdistglobpi(:,:)
c
c     IMPROP
c
      integer, allocatable :: niproplocpi(:),impropglobpi(:,:)
c
c     IMPTOR
c
      integer, allocatable :: nitorslocpi(:),imptorglobpi(:,:)
c
c     TORSION
c
      integer, allocatable :: ntorslocpi(:),torsglobpi(:,:)
c
c     PITORSION
c
      integer, allocatable :: npitorslocpi(:),pitorsglobpi(:,:)
c
c     STRETCH-TORSION
c
      integer, allocatable :: nstrtorlocpi(:),strtorglobpi(:,:)
c
c     TORSION-TORSION
c
      integer, allocatable :: ntortorlocpi(:),tortorglobpi(:,:)
c
c     ANGLE
c
      integer, allocatable :: nanglelocpi(:),angleglobpi(:,:)
c
c     CHARGE
c
      integer, allocatable :: nionreclocpi(:),chgrecglobpi(:,:)
      integer, allocatable :: nionlocpi(:),chgglobpi(:,:)
      integer, allocatable :: nionlocnlpi(:),chgglobnlpi(:,:)
      integer, allocatable :: nelstpi(:,:),elstpi(:,:,:)
c
c     MULTIPOLE
c
      integer, allocatable :: npolereclocpi(:),polerecglobpi(:,:)
      integer, allocatable :: npolelocpi(:),poleglobpi(:,:)
      integer, allocatable :: npoleblocpi(:)
      integer, allocatable :: polelocpi(:,:)
      integer, allocatable :: npolelocnlpi(:),poleglobnlpi(:,:)
c 
c      POLARIZATION
c
      integer, allocatable :: nualtpi(:)
      real*8, allocatable :: udaltpi(:,:,:,:),upaltpi(:,:,:,:)
      real*8, allocatable :: uindpi(:,:,:),uinppi(:,:,:)
c
c     VDW
c
      integer, allocatable :: nvdwblocpi(:),vdwglobpi(:,:)
      integer, allocatable :: nvdwlocnlpi(:),vdwglobnlpi(:,:)
      integer, allocatable :: nvlstpi(:,:),vlstpi(:,:,:)
c
c     STAT
c
      real*8, allocatable :: etot_sumpi(:),etot2_sumpi(:)
      real*8, allocatable :: eint_sumpi(:),eint2_sumpi(:)
      real*8, allocatable :: epot_sumpi(:),epot2_sumpi(:)
      real*8, allocatable :: ekin_sumpi(:),ekin2_sumpi(:)
      real*8, allocatable :: temp_sumpi(:),temp2_sumpi(:)
      real*8, allocatable :: pres_sumpi(:),pres2_sumpi(:)
      real*8, allocatable :: dens_sumpi(:),dens2_sumpi(:)
      real*8 :: ekprim_ave, epotpi_ave, temppi_ave
c
c     TIME
c
      real*8, allocatable :: timesteppi(:)
      save

      contains
      subroutine allocpi(init,istep)
      use angle
      use atoms
      use bath
      use bitor
      use domdec
      use molcul
      use neigh
      use sizes
      use pitors
      use potent
      use tors
      use uprior
      use units
      use mdstuf
      implicit none
      logical init
      integer nblocrecdirmax,nlocnlmax,istep,modnl
      integer ibead,i,ierr
      real*8, allocatable :: WORK(:)
c
      if (init) then
        if(ranktot.eq.0) then
           Ekcentroid=nfree*nbeads*kelvin*gasconst        
        
          if (allocated(eigmat)) deallocate (eigmat)
          allocate (eigmat(nbeads,nbeads))
          if (allocated(eigmattr)) deallocate (eigmattr)
          allocate (eigmattr(nbeads,nbeads))
          if (allocated(omkpi)) deallocate (omkpi)
          allocate (omkpi(nbeads))
          
          allocate(WORK(3*nbeads))
          eigmat=0
          DO i=1,nbeads-1
            eigmat(i,i)=2
            eigmat(i+1,i)=-1
            eigmat(i,i+1)=-1
          ENDDO
          eigmat(1,nbeads)=-1
          eigmat(nbeads,1)=-1
          eigmat(nbeads,nbeads)=2
          call DSYEV('V','U',nbeads,eigMat,nbeads, 
     $      omkpi,WORK,3*nbeads,ierr)		
          omkpi(1)=0
          omkpi(:)=sqrt(omkpi)*(nbeads*boltzmann*kelvin/hbar)
          eigmattr=transpose(eigmat)
          deallocate(WORK)
        endif
      
      
        if (allocated(nlocpi)) deallocate (nlocpi)
        allocate (nlocpi(nbeadsloc))
        if (allocated(nblocpi)) deallocate (nblocpi)
        allocate (nblocpi(nbeadsloc))
        if (allocated(nlocrecpi)) deallocate (nlocrecpi)
        allocate (nlocrecpi(nbeadsloc))
        if (allocated(nblocrecpi)) deallocate (nblocrecpi)
        allocate (nblocrecpi(nbeadsloc))
        if (allocated(nlocnlpi)) deallocate (nlocnlpi)
        allocate (nlocnlpi(nbeadsloc))
        if (allocated(nblocrecdirpi)) deallocate (nblocrecdirpi)
        allocate (nblocrecdirpi(nbeadsloc))
        if (allocated(globpi)) deallocate (globpi)
        allocate (globpi(n,nbeadsloc))
        if (allocated(locpi)) deallocate (locpi)
        allocate (locpi(n,nbeadsloc))
        if (allocated(repartpi)) deallocate (repartpi)
        allocate (repartpi(n,nbeadsloc))
        if (allocated(repartrecpi)) deallocate (repartrecpi)
        allocate (repartrecpi(n,nbeadsloc))
        if (allocated(domlenpi)) deallocate (domlenpi)
        allocate (domlenpi(nproc,nbeadsloc))
        if (allocated(domlenrecpi)) deallocate (domlenrecpi)
        allocate (domlenrecpi(nproc,nbeadsloc))
        if (allocated(domlenpolepi)) deallocate (domlenpolepi)
        allocate (domlenpolepi(nproc,nbeadsloc))
        if (allocated(domlenpolerecpi)) deallocate (domlenpolerecpi)
        allocate (domlenpolerecpi(nproc,nbeadsloc))
        if (allocated(globrecpi)) deallocate (globrecpi)
        allocate (globrecpi(n,nbeadsloc))
        if (allocated(locrecpi)) deallocate (locrecpi)
        allocate (locrecpi(n,nbeadsloc))
        if (allocated(globrec1pi)) deallocate (globrec1pi)
        allocate (globrec1pi(n,nbeadsloc))
        if (allocated(locrec1pi)) deallocate (locrec1pi)
        allocate (locrec1pi(n,nbeadsloc))
        if (allocated(bufbegrecpi)) deallocate (bufbegrecpi)
        allocate (bufbegrecpi(nproc,nbeadsloc))
        if (allocated(bufbegpolepi)) deallocate (bufbegpolepi)
        allocate (bufbegpolepi(nproc,nbeadsloc))
        if (allocated(bufbegpi)) deallocate (bufbegpi)
        allocate (bufbegpi(nproc,nbeadsloc))
        if (allocated(buflen1pi)) deallocate (buflen1pi)
        allocate (buflen1pi(nproc,nbeadsloc))
        if (allocated(buflen2pi)) deallocate (buflen2pi)
        allocate (buflen2pi(nproc,nbeadsloc))
        if (allocated(bufbeg1pi)) deallocate (bufbeg1pi)
        if (allocated (buf1pi)) deallocate(buf1pi)
        allocate (buf1pi(n,nbeadsloc))
        if (allocated (buf2pi)) deallocate(buf2pi)
        allocate (buf2pi(n,nbeadsloc))
        allocate (bufbeg1pi(nproc,nbeadsloc))
        if (allocated(bufbeg2pi)) deallocate (bufbeg2pi)
        allocate (bufbeg2pi(nproc,nbeadsloc))
        if (allocated(pospi)) deallocate (pospi)
        allocate (pospi(3,n,nbeadsloc))
        if (allocated(velpi)) deallocate (velpi)
        allocate (velpi(3,n,nbeadsloc))
        if (allocated(api)) deallocate (api)
        allocate (api(3,n,nbeadsloc))
        if (allocated(masspi)) deallocate (masspi)
        allocate (masspi(n))
        if (allocated(nmolelocpi)) deallocate (nmolelocpi)
        allocate (nmolelocpi(nbeadsloc))
        if (allocated(molculeglobpi)) deallocate (molculeglobpi)
        allocate (molculeglobpi(nmol,nbeadsloc))
        if (allocated(ineignlpi)) deallocate (ineignlpi)
        allocate (ineignlpi(n,nbeadsloc))
        if (allocated(locnlpi)) deallocate (locnlpi)
        allocate (locnlpi(n,nbeadsloc))


        if (allocated(epotpi)) deallocate (epotpi)
        allocate (epotpi(nbeadsloc))
        if (allocated(etotpi)) deallocate (etotpi)
        allocate (etotpi(nbeadsloc))
        if (allocated(eksumpi)) deallocate (eksumpi)
        allocate (eksumpi(nbeadsloc))
        if (allocated(dedv_pi)) deallocate (dedv_pi)
        allocate (dedv_pi(nbeadsloc))
        if (allocated(ekinpi)) deallocate (ekinpi)
        allocate (ekinpi(3,3,nbeadsloc))

        if (use_vdw) then
          if (allocated(nvdwblocpi)) deallocate (nvdwblocpi)
          allocate (nvdwblocpi(nbeadsloc))
          if (allocated(vdwglobpi)) deallocate (vdwglobpi)
          allocate (vdwglobpi(n,nbeadsloc))
          if (allocated(vdwglobnlpi)) deallocate (vdwglobnlpi)
          allocate (vdwglobnlpi(n,nbeadsloc))
          if (allocated(nvdwlocnlpi)) deallocate (nvdwlocnlpi) 
          allocate (nvdwlocnlpi(nbeadsloc))
        end if

        if (use_bond) then
         if (allocated(nbondlocpi)) deallocate(nbondlocpi)
         allocate (nbondlocpi(nbeadsloc))
         if (allocated(bndglobpi)) deallocate(bndglobpi)
         allocate (bndglobpi(8*n,nbeadsloc))
        end if

        if (use_strbnd) then
          if (allocated(nstrbndlocpi)) deallocate(nstrbndlocpi)
          allocate (nstrbndlocpi(nbeadsloc))
          if (allocated(strbndglobpi)) deallocate(strbndglobpi)
          allocate (strbndglobpi(nangle,nbeadsloc))
        end if
        
        if (use_urey) then
          if (allocated(nureylocpi)) deallocate(nureylocpi)
          allocate (nureylocpi(nbeadsloc))
          if (allocated(ureyglobpi)) deallocate(ureyglobpi)
          allocate (ureyglobpi(nangle,nbeadsloc))
        end if

        if (use_angang) then
          if (allocated(nanganglocpi)) deallocate(nanganglocpi)
          allocate (nanganglocpi(nbeadsloc))
          if (allocated(angangglobpi)) deallocate(angangglobpi)
          allocate (angangglobpi(ntors,nbeadsloc))
        end if

        if (use_opbend) then
          if (allocated(nopbendlocpi)) deallocate(nopbendlocpi)
          allocate (nopbendlocpi(nbeadsloc))
          if (allocated(opbendglobpi)) deallocate(opbendglobpi)
          allocate (opbendglobpi(nangle,nbeadsloc))
        end if

        if (use_opdist) then
          if (allocated(nopdistlocpi)) deallocate(nopdistlocpi)
          allocate (nopdistlocpi(nbeadsloc))
          if (allocated(opdistglobpi)) deallocate(opdistglobpi)
          allocate (opdistglobpi(n,nbeadsloc))
        end if

        if (use_improp) then
          if (allocated(niproplocpi)) deallocate(niproplocpi)
          allocate (niproplocpi(nbeadsloc))
          if (allocated(impropglobpi)) deallocate(impropglobpi)
          allocate (impropglobpi(6*n,nbeadsloc))
        end if

        if (use_imptor) then
          if (allocated(nitorslocpi)) deallocate(nitorslocpi)
          allocate (nitorslocpi(nbeadsloc))
          if (allocated(imptorglobpi)) deallocate(imptorglobpi)
          allocate (imptorglobpi(6*n,nbeadsloc))
        end if

        if (use_tors) then
          if (allocated(ntorslocpi)) deallocate(ntorslocpi)
          allocate (ntorslocpi(nbeadsloc))
          if (allocated(torsglobpi)) deallocate(torsglobpi)
          allocate (torsglobpi(6*n,nbeadsloc))
        end if

        if (use_pitors) then 
          if (allocated(npitorslocpi)) deallocate(npitorslocpi)
          allocate (npitorslocpi(nbeadsloc))
          if (allocated(pitorsglobpi)) deallocate(pitorsglobpi)
          allocate (pitorsglobpi(ntors,nbeadsloc))
        end if

        if (use_strtor) then
          if (allocated(nstrtorlocpi)) deallocate(nstrtorlocpi)
          allocate (nstrtorlocpi(nbeadsloc))
          if (allocated(strtorglobpi)) deallocate(strtorglobpi)
          allocate (strtorglobpi(ntors,nbeadsloc))
        end if

        if (use_tortor) then
          if (allocated(ntortorlocpi)) deallocate(ntortorlocpi)
          allocate (ntortorlocpi(nbeadsloc))
          if (allocated(tortorglobpi)) deallocate(tortorglobpi)
          allocate (tortorglobpi(nbitor,nbeadsloc))
        end if

        if (use_angle) then
          if (allocated(nanglelocpi)) deallocate(nanglelocpi)
          allocate (nanglelocpi(nbeadsloc))
          if (allocated(angleglobpi)) deallocate(angleglobpi)
          allocate (angleglobpi(4*n,nbeadsloc))
        end if

        if (use_charge) then
          if (allocated(nionreclocpi)) deallocate(nionreclocpi)
          allocate (nionreclocpi(nbeadsloc))
          if (allocated(chgrecglobpi)) deallocate(chgrecglobpi)
          allocate (chgrecglobpi(n,nbeadsloc))
          if (allocated(nionlocpi)) deallocate(nionlocpi)
          allocate (nionlocpi(nbeadsloc))
          if (allocated(chgglobpi)) deallocate(chgglobpi)
          allocate (chgglobpi(n,nbeadsloc))
          if (allocated(nionlocnlpi)) deallocate(nionlocnlpi)
          allocate (nionlocnlpi(nbeadsloc))
          if (allocated(chgglobnlpi)) deallocate(chgglobnlpi)
          allocate (chgglobnlpi(n,nbeadsloc))
        end if

        if (use_polar) then
          if (allocated(nualtpi)) deallocate(nualtpi)
          allocate (nualtpi(nbeadsloc))
          if (allocated(udaltpi)) deallocate(udaltpi)
          allocate (udaltpi(maxualt,3,n,nbeadsloc))
          if (allocated(upaltpi)) deallocate(upaltpi)
          allocate (upaltpi(maxualt,3,n,nbeadsloc))
          if (allocated(uindpi)) deallocate(uindpi)
          allocate (uindpi(3,n,nbeadsloc))
          if (allocated(uinppi)) deallocate(uinppi)
          allocate (uinppi(3,n,nbeadsloc))
        end if

        if (use_mpole) then
          if (allocated(npolereclocpi)) deallocate(npolereclocpi)
          allocate (npolereclocpi(nbeadsloc))
          if (allocated(polerecglobpi)) deallocate(polerecglobpi)
          allocate (polerecglobpi(n,nbeadsloc))
          if (allocated(poleglobpi)) deallocate(poleglobpi)
          allocate (poleglobpi(n,nbeadsloc))
          if (allocated(polelocpi)) deallocate(polelocpi)
          allocate (polelocpi(n,nbeadsloc))
          if (allocated(npolelocpi)) deallocate(npolelocpi)
          allocate (npolelocpi(nbeadsloc))
          if (allocated(npoleblocpi)) deallocate(npoleblocpi)
          allocate (npoleblocpi(nbeadsloc))
          if (allocated(npolelocnlpi)) deallocate(npolelocnlpi)
          allocate (npolelocnlpi(nbeadsloc))
          if (allocated(poleglobnlpi)) deallocate(poleglobnlpi)
          allocate (poleglobnlpi(n,nbeadsloc))
        end if

 
        if (allocated(etot_sumpi)) deallocate(etot_sumpi)
        allocate (etot_sumpi(nbeadsloc))
        if (allocated(etot2_sumpi)) deallocate(etot2_sumpi)
        allocate (etot2_sumpi(nbeadsloc))
        if (allocated(eint_sumpi)) deallocate(eint_sumpi)
        allocate (eint_sumpi(nbeadsloc))
        if (allocated(eint2_sumpi)) deallocate(eint2_sumpi)
        allocate (eint2_sumpi(nbeadsloc))
        if (allocated(epot_sumpi)) deallocate(epot_sumpi)
        allocate (epot_sumpi(nbeadsloc))
        if (allocated(epot2_sumpi)) deallocate(epot2_sumpi)
        allocate (epot2_sumpi(nbeadsloc))
        if (allocated(ekin_sumpi)) deallocate(ekin_sumpi)
        allocate (ekin_sumpi(nbeadsloc))
        if (allocated(ekin2_sumpi)) deallocate(ekin2_sumpi)
        allocate (ekin2_sumpi(nbeadsloc))
        if (allocated(temp_sumpi)) deallocate(temp_sumpi)
        allocate (temp_sumpi(nbeadsloc))
        if (allocated(temp2_sumpi)) deallocate(temp2_sumpi)
        allocate (temp2_sumpi(nbeadsloc))
        if (allocated(pres_sumpi)) deallocate(pres_sumpi)
        allocate (pres_sumpi(nbeadsloc))
        if (allocated(pres2_sumpi)) deallocate(pres2_sumpi)
        allocate (pres2_sumpi(nbeadsloc))
        if (allocated(dens_sumpi)) deallocate(dens_sumpi)
        allocate (dens_sumpi(nbeadsloc))
        if (allocated(dens2_sumpi)) deallocate(dens2_sumpi)
        allocate (dens2_sumpi(nbeadsloc))

        if (allocated(timesteppi)) deallocate(timesteppi)
        allocate (timesteppi(nbeadsloc))
        do ibead = 1, nbeadsloc
          nlocnlpi(ibead) = nlocnl
        end do
      end if


c
      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return

      nlocnlmax = maxval(nlocnlpi)
      if (use_vdw) then
        if (allocated(nvlstpi)) deallocate (nvlstpi)
        allocate (nvlstpi(nlocnlmax,nbeadsloc))
        if (allocated(vlstpi)) deallocate (vlstpi)
        allocate (vlstpi(maxvlst,nlocnlmax,nbeadsloc))
      end if
      if (use_mpole.or.use_charge) then
        if (allocated(nelstpi)) deallocate (nelstpi)
        allocate (nelstpi(nlocnlmax,nbeadsloc))
        if (allocated(elstpi)) deallocate (elstpi)
        allocate (elstpi(maxelst,nlocnlmax,nbeadsloc))
      end if

      end subroutine

      subroutine initbead(ibead,istep)
      use angang
      use angle
      use atoms
      use atmlst
      use atmtyp
      use bond
      use charge
      use domdec
      use energi
      use improp
      use imptor
      use moldyn
      use mpole
      use neigh
      use opbend
      use opdist
      use pitors
      use polar
      use potent
      use sizes
      use stat
      use strbnd
      use strtor
      use timestat
      use tors
      use tortor
      use uprior
      use urey
      use vdw
      use virial
      implicit none
      integer ibead,i,istep,modnl

c
c     positions, speed, mass
c
      pospi(1,:,ibead) = x
      pospi(2,:,ibead) = y
      pospi(3,:,ibead) = z
      velpi(:,:,ibead) = v
      api(:,:,ibead) = a

      epotpi(ibead)=epotpi_loc
      eksumpi(ibead)=eksumpi_loc
      etotpi(ibead)=epotpi_loc+eksumpi_loc
      ekinpi(:,:,ibead)=ekinpi_loc
      dedv_pi(ibead) = dedv

      nlocpi(ibead) = nloc
      globpi(:,ibead) = glob
c
c     STAT
c
      etot_sumpi(ibead) = etot_sum
      etot2_sumpi(ibead)= etot2_sum
      eint_sumpi(ibead) = eint_sum
      eint2_sumpi(ibead)= eint2_sum
      epot_sumpi(ibead) = epot_sum
      epot2_sumpi(ibead)= epot2_sum
      ekin_sumpi(ibead) = ekin_sum
      ekin2_sumpi(ibead)= ekin2_sum
      temp_sumpi(ibead) = temp_sum
      temp2_sumpi(ibead)= temp2_sum
      pres_sumpi(ibead) = pres_sum
      pres2_sumpi(ibead)= pres2_sum
      dens_sumpi(ibead) = dens_sum
      dens2_sumpi(ibead)= dens2_sum

      if(nbeadsloc.eq.1) return

      modnl = mod(istep,ineigup)
c
c     parallelism
c
c      nlocpi(ibead) = nloc
      nblocpi(ibead) = nbloc
      nlocrecpi(ibead) = nlocrec
      nblocrecpi(ibead) = nblocrec
      nlocnlpi(ibead) = nlocnl
      nblocrecdirpi(ibead) = nblocrecdir
c      globpi(:,ibead) = glob
      locpi(:,ibead) = loc
      ineignlpi(:,ibead) = ineignl
      locnlpi(:,ibead) = locnl
      repartpi(:,ibead) = repart
      repartrecpi(:,ibead) = repartrec
      domlenpi(:,ibead) = domlen
      domlenrecpi(:,ibead) = domlenrec
      domlenpolepi(:,ibead) = domlenpole
      domlenpolerecpi(:,ibead) = domlenpolerec
      globrecpi(:,ibead) = globrec
      locrecpi(:,ibead) = locrec
      globrec1pi(:,ibead) = globrec1
      locrec1pi(:,ibead) = locrec1
      bufbegrecpi(:,ibead) = bufbegrec
      bufbegpolepi(:,ibead) = bufbegpole
      bufbegpi(:,ibead) = bufbeg
      buflen1pi(:,ibead) = buflen1
      buf1pi(1:nblocrecdir,ibead) = buf1
      buflen2pi(:,ibead) = buflen2
      buf2pi(1:nblocrecdir,ibead) = buf2
      bufbeg1pi(:,ibead) =  bufbeg1
      bufbeg2pi(:,ibead) = bufbeg2

c      epotpi(ibead)=epotpi_loc
c      eksumpi(ibead)=eksumpi_loc
c      etotpi(ibead)=epotpi_loc+eksumpi_loc
c      ekinpi(:,:,ibead)=ekinpi_loc
c      dedv_pi(ibead) = dedv
c
c     positions, speed, mass
c
c      pospi(1,:,ibead) = x
c      pospi(2,:,ibead) = y
c      pospi(3,:,ibead) = z
c      velpi(:,:,ibead) = v
c      api(:,:,ibead) = a
      masspi = mass
c
c     VDW 
c
      if (use_vdw) then
        nvdwblocpi(ibead) = nvdwbloc
        vdwglobpi(:,ibead) = vdwglob
        vdwglobnlpi(:,ibead) = vdwglobnl
        nvdwlocnlpi(ibead) = nvdwlocnl
c        if (modnl.eq.0) then
c          nvlstpi(1:nlocnl,ibead) = nvlst
c          vlstpi(:,1:nlocnl,ibead) = vlst
c        end if
      end if
c
c     BONDS
c
      if (use_bond) then
        bndglobpi(:,ibead) = bndglob
        nbondlocpi(ibead) = nbondloc
      end if
c
c     STRETCH-BEND
c
      if (use_strbnd) then
        strbndglobpi(:,ibead) = strbndglob
        nstrbndlocpi(ibead) = nstrbndloc
      end if
c
c     UREY-BRADLEY
c
      if (use_urey) then
        ureyglobpi(:,ibead) = ureyglob
        nureylocpi(ibead) = nureyloc
      end if
c
c     ANGlE-ANGLE
c
      if (use_angang) then
        angangglobpi(:,ibead) = angangglob
        nanganglocpi(ibead) = nangangloc
      end if
c
c     OP-BENDING
c
      if (use_opbend) then
        opbendglobpi(:,ibead) = opbendglob
        nopbendlocpi(ibead) = nopbendloc
      end if
c
c     OP-DIST
c
      if (use_opdist) then
        opdistglobpi(:,ibead) = opdistglob
        nopdistlocpi(ibead) = nopdistloc
      end if
c
c     IMPROP
c
      if (use_improp) then
        impropglobpi(:,ibead) = impropglob
        niproplocpi(ibead) = niproploc
      end if
c
c     IMPTOR
c
      if (use_imptor) then
        imptorglobpi(:,ibead) = imptorglob
        nitorslocpi(ibead) = nitorsloc
      end if
c
c     TORSION
c
      if (use_tors) then
        torsglobpi(:,ibead) = torsglob
        ntorslocpi(ibead) = ntorsloc
      end if
c
c     PITORSION
c
      if (use_pitors) then
        pitorsglobpi(:,ibead) = pitorsglob
        npitorslocpi(ibead) = npitorsloc
      end if
c
c     STRETCH-TORSION
c
      if (use_strtor) then
        strtorglobpi(:,ibead) = strtorglob
        nstrtorlocpi(ibead) = nstrtorloc
      end if
c
c     TORSION-TORSION
c
      if (use_tortor) then
        tortorglobpi(:,ibead) = tortorglob
        ntortorlocpi(ibead) = ntortorloc
      end if
c
c     ANGLE
c
      if (use_angle) then
        angleglobpi(:,ibead) = angleglob
        nanglelocpi(ibead) = nangleloc
      end if
c
c     CHARGE
c
      if (use_charge) then
        chgrecglobpi(:,ibead) = chgrecglob
        nionreclocpi(ibead) = nionrecloc
        chgglobpi(:,ibead) = chgglob
        nionlocpi(ibead) = nionloc
        chgglobnlpi(:,ibead) = chgglobnl
        nionlocnlpi(ibead) = nionlocnl
      end if
c      if ((use_charge.or.use_mpole).and.(modnl.eq.0)) then
c        nelstpi(1:nlocnl,ibead) = nelst
c        elstpi(:,1:nlocnl,ibead) = elst
c      end if
c
c     MULTIPOLE
c
      if (use_mpole) then
        polerecglobpi(:,ibead) = polerecglob
        npolereclocpi(ibead) = npolerecloc
        poleglobpi(:,ibead) = poleglob
        polelocpi(:,ibead) = poleloc
        npolelocpi(ibead) = npoleloc
        npoleblocpi(ibead) = npolebloc
        poleglobnlpi(:,ibead) = poleglobnl
        npolelocnlpi(ibead) = npolelocnl
      end if
c
c     POLARIZATION
c
      if (use_polar) then
        nualtpi(ibead) = nualt
        udaltpi(:,:,:,ibead) = udalt
        upaltpi(:,:,:,ibead) = upalt
        uindpi(:,:,ibead) = uind
        uinppi(:,:,ibead) = uinp
      end if

c
c     STAT
c
c      etot_sumpi(ibead) = etot_sum
c      etot2_sumpi(ibead)= etot2_sum
c      eint_sumpi(ibead) = eint_sum
c      eint2_sumpi(ibead)= eint2_sum
c      epot_sumpi(ibead) = epot_sum
c      epot2_sumpi(ibead)= epot2_sum
c      ekin_sumpi(ibead) = ekin_sum
c      ekin2_sumpi(ibead)= ekin2_sum
c      temp_sumpi(ibead) = temp_sum
c      temp2_sumpi(ibead)= temp2_sum
c      pres_sumpi(ibead) = pres_sum
c      pres2_sumpi(ibead)= pres2_sum
c      dens_sumpi(ibead) = dens_sum
c      dens2_sumpi(ibead)= dens2_sum
c
c     TIME
c
      timesteppi(ibead) = timestep

      end subroutine
c
      subroutine savebeadnl(ibead,istep)
      use domdec
      use neigh
      use potent
      implicit none
      integer ibead,istep,modnl

      if (nbeadsloc.eq.1) return

      modnl = mod(istep,ineigup)
      if (modnl.ne.0) return

      if (use_vdw) then
        nvlstpi(1:nlocnl,ibead) = nvlst
        vlstpi(:,1:nlocnl,ibead) = vlst
      end if
      
      if ((use_charge).or.(use_mpole)) then
        nelstpi(1:nlocnl,ibead) = nelst
        elstpi(:,1:nlocnl,ibead) = elst
      end if
      return
      end
      
c
      subroutine pushbead(ibead,istep)
      use angang
      use angle
      use atoms
      use atmlst
      use atmtyp
      use bond
      use charge
      use domdec
      use energi
      use improp
      use imptor
      use moldyn
      use mpole
      use neigh
      use opbend
      use opdist
      use pitors
      use polar
      use potent
      use sizes
      use stat
      use strbnd
      use strtor
      use timestat
      use tors
      use tortor
      use uprior
      use urey
      use vdw
      use virial
      implicit none
      integer ibead,istep,modnl

      epotpi_loc=epotpi(ibead)
      eksumpi_loc=eksumpi(ibead)
      ekinpi_loc=ekinpi(:,:,ibead)
      dedv=dedv_pi(ibead)

      x = pospi(1,:,ibead)
      y = pospi(2,:,ibead)
      z = pospi(3,:,ibead) 
      v = velpi(:,:,ibead)
      a= api(:,:,ibead) 

      nloc = nlocpi(ibead)
      glob = globpi(:,ibead)
c
c     STAT
c
      etot_sum  = etot_sumpi(ibead)
      etot2_sum = etot2_sumpi(ibead)
      eint_sum  = eint_sumpi(ibead)
      eint2_sum = eint2_sumpi(ibead)
      epot_sum  = epot_sumpi(ibead)
      epot2_sum = epot2_sumpi(ibead)
      ekin_sum  = ekin_sumpi(ibead)
      ekin2_sum = ekin2_sumpi(ibead)
      temp_sum  = temp_sumpi(ibead)
      temp2_sum = temp2_sumpi(ibead)
      pres_sum  = pres_sumpi(ibead)
      pres2_sum = pres2_sumpi(ibead)
      dens_sum  = dens_sumpi(ibead)
      dens2_sum = dens2_sumpi(ibead)
c
      if (nbeadsloc.eq.1) return

      modnl = mod(istep,ineigup)
c
c     parallelism
c
c      nloc = nlocpi(ibead)
      nbloc = nblocpi(ibead)
      nlocrec = nlocrecpi(ibead)
      nblocrec = nblocrecpi(ibead)
      nlocnl = nlocnlpi(ibead)
      nblocrecdir = nblocrecdirpi(ibead)
c      glob = globpi(:,ibead)
      loc = locpi(:,ibead)
      ineignl = ineignlpi(:,ibead)
      locnl = locnlpi(:,ibead)
      repart = repartpi(:,ibead)
      repartrec = repartrecpi(:,ibead)
      domlen = domlenpi(:,ibead)
      domlenrec = domlenrecpi(:,ibead)
      domlenpole = domlenpolepi(:,ibead)
      domlenpolerec = domlenpolerecpi(:,ibead)
      globrec = globrecpi(:,ibead)
      locrec = locrecpi(:,ibead)
      globrec1 = globrec1pi(:,ibead)
      locrec1 = locrec1pi(:,ibead)
      bufbegrec = bufbegrecpi(:,ibead)
      bufbegpole = bufbegpolepi(:,ibead)
      bufbeg = bufbegpi(:,ibead)
      buflen1 = buflen1pi(:,ibead)
      buf1 = buf1pi(:,ibead)
      buflen2 = buflen2pi(:,ibead)
      buf2 = buf2pi(:,ibead)
      bufbeg1 = bufbeg1pi(:,ibead)
      bufbeg2 = bufbeg2pi(:,ibead)

c      epotpi_loc=epotpi(ibead)
c      eksumpi_loc=eksumpi(ibead)
c      ekinpi_loc=ekinpi(:,:,ibead)
c      dedv=dedv_pi(ibead)
c
c     positions, speed, mass
c
c      x = pospi(1,:,ibead)
c      y = pospi(2,:,ibead)
c      z = pospi(3,:,ibead) 
c      v = velpi(:,:,ibead)
c      a = api(:,:,ibead)
      mass = masspi
c
c     VDW
c
      if (use_vdw) then
        nvdwbloc = nvdwblocpi(ibead)
        vdwglob = vdwglobpi(:,ibead)
        vdwglobnl = vdwglobnlpi(:,ibead)
        nvdwlocnl = nvdwlocnlpi(ibead)
c        if (modnl.eq.0) then
          nvlst = nvlstpi(1:nlocnl,ibead)
          vlst = vlstpi(:,1:nlocnl,ibead)
c        end if
      end if
c
c     BOND
c
      if (use_bond) then
        nbondloc = nbondlocpi(ibead)
        bndglob = bndglobpi(:,ibead)
      end if
c
c     STRETCH-BEND
c
      if (use_strbnd) then
        nstrbndloc = nstrbndlocpi(ibead)
        strbndglob = strbndglobpi(:,ibead)
      end if
c
c     UREY-BRADLEY
c
      if (use_urey) then
        nureyloc = nureylocpi(ibead)
        ureyglob = ureyglobpi(:,ibead)
      end if
c
c     ANGLE-ANGLE
c
      if (use_angang) then
        nangangloc = nanganglocpi(ibead)
        angangglob = angangglobpi(:,ibead)
      end if
c
c     OP-BENDING
c
      if (use_opbend) then
        nopbendloc = nopbendlocpi(ibead)
        opbendglob = opbendglobpi(:,ibead)
      end if
c
c     OP-DIST
c
      if (use_opdist) then
        nopdistloc = nopdistlocpi(ibead)
        opdistglob = opdistglobpi(:,ibead)
      end if
c
c     IMPROP
c
      if (use_improp) then
        niproploc = niproplocpi(ibead)
        impropglob = impropglobpi(:,ibead)
      end if
c
c     IMPTOR
c
      if (use_imptor) then
        nitorsloc = nitorslocpi(ibead)
        imptorglob = imptorglobpi(:,ibead)
      end if
c
c     TORSION
c
      if (use_tors) then
        ntorsloc = ntorslocpi(ibead)
        torsglob = torsglobpi(:,ibead)
      end if
c
c     PITORSION
c
      if (use_pitors) then
        npitorsloc = npitorslocpi(ibead)
        pitorsglob = pitorsglobpi(:,ibead)
      end if
c
c     STRETCH-TORSION
c
      if (use_strtor) then
        nstrtorloc = nstrtorlocpi(ibead)
        strtorglob = strtorglobpi(:,ibead)
      end if
c
c     TORSION-TORSION
c
      if (use_tortor) then
        ntortorloc = ntortorlocpi(ibead)
        tortorglob = tortorglobpi(:,ibead)
      end if
c
c     ANGLE
c
      if (use_angle) then
        nangleloc = nanglelocpi(ibead)
        angleglob = angleglobpi(:,ibead)
      end if
c
c     CHARGE
c
      if (use_charge) then
        nionrecloc = nionreclocpi(ibead)
        chgrecglob = chgrecglobpi(:,ibead)
        nionloc = nionlocpi(ibead)
        chgglob = chgglobpi(:,ibead)
        chgglobnl = chgglobnlpi(:,ibead)
        nionlocnl = nionlocnlpi(ibead)
      end if
      if (use_charge.or.use_mpole) then
        nelst = nelstpi(1:nlocnl,ibead)
        elst = elstpi(:,1:nlocnl,ibead)
      end if
c
c     MULTIPOLE
c
      if (use_mpole) then
        npolerecloc = npolereclocpi(ibead)
        polerecglob = polerecglobpi(:,ibead)
        npoleloc = npolelocpi(ibead)
        npolebloc = npoleblocpi(ibead)
        poleglob = poleglobpi(:,ibead)
        poleloc = polelocpi(:,ibead)
        poleglobnl = poleglobnlpi(:,ibead)
        npolelocnl = npolelocnlpi(ibead)
      end if
c
c     POLARIZATION
c
      if (use_polar) then
        nualt = nualtpi(ibead)
        udalt = udaltpi(:,:,:,ibead)
        upalt = upaltpi(:,:,:,ibead)
        uind = uindpi(:,:,ibead)
        uinp = uinppi(:,:,ibead)
      end if
c
c     STAT
c
c      etot_sum  = etot_sumpi(ibead)
c      etot2_sum = etot2_sumpi(ibead)
c      eint_sum  = eint_sumpi(ibead)
c      eint2_sum = eint2_sumpi(ibead)
c      epot_sum  = epot_sumpi(ibead)
c      epot2_sum = epot2_sumpi(ibead)
c      ekin_sum  = ekin_sumpi(ibead)
c      ekin2_sum = ekin2_sumpi(ibead)
c      temp_sum  = temp_sumpi(ibead)
c      temp2_sum = temp2_sumpi(ibead)
c      pres_sum  = pres_sumpi(ibead)
c      pres2_sum = pres2_sumpi(ibead)
c      dens_sum  = dens_sumpi(ibead)
c      dens2_sum = dens2_sumpi(ibead)
c
c     TIME
c 
      timestep = timesteppi(ibead)

      end subroutine


      subroutine compute_observables_pi(pos,forces,dedv_mean
     &                        ,ek0,ekprim,ekvir,presvir)
        use atoms
        use units
        use atmtyp
        use math
        use boxes
        use mdstuf
        use bath
        IMPLICIT NONE
        real*8, intent(in) :: pos(:,:,:),forces(:,:,:),dedv_mean,ek0
        real*8, intent(out) :: ekprim,ekvir,presvir
        real*8, allocatable :: centroid(:,:)
        real*8 :: omp,omp2!,hbar
        integer :: ibead,i,j

        !hbar=(planck*1.d11*avogadro)/(2*pi)

        allocate(centroid(3,n))
        centroid(:,:)=0.d0
        DO ibead=1,nbeads
          DO i=1,n
        !  if (atomic(i).eq.0) cycle
            DO j=1,3
              centroid(j,i)=centroid(j,i)+pos(j,i,ibead)
            ENDDO
           ENDDO
        ENDDO  
        centroid(:,:)=centroid(:,:)/REAL(nbeads)

        omp=nbeads*boltzmann*kelvin/hbar        
        omp2=omp*omp

c       COMPUTE PRIMITIVE KINETIC ENERGY
        ekprim=0.d0
        DO ibead=1,nbeads-1
          DO i=1,n
        !    if (atomic(i).eq.0) cycle
            DO j=1,3
              ekprim = ekprim - 0.5*mass(i)*omp2
     &          *(pos(j,i,ibead+1)-pos(j,i,ibead))**2
            ENDDO
          ENDDO
        ENDDO  
        DO i=1,n
         ! if (atomic(i).eq.0) cycle
          DO j=1,3
            ekprim = ekprim - 0.5*mass(i)*omp2
     &          *(pos(j,i,nbeads)-pos(j,i,1))**2
          ENDDO
        ENDDO  
        ekprim = (ekprim/nbeads
     &          + 0.5*nbeads*nfree*boltzmann*kelvin)/convert

c       COMPUTE VIRIAL KINETIC ENERGY
        ekvir=0.d0
        DO ibead=1,nbeads
          DO i=1,n
         !  if (atomic(i).eq.0) cycle
            DO j=1,3
              ekvir=ekvir+(pos(j,i,ibead)-centroid(j,i))
     &                      *forces(j,i,ibead)
            ENDDO
           ENDDO
        ENDDO

        presvir = prescon*( -dedv_mean + (ek0
     &               - ekvir)/(3*nbeads*volbox) )

        ekvir=0.5d0*(nfree*boltzmann*kelvin/convert-ekvir/nbeads)


      end subroutine
      
      subroutine compute_observables_pi_4site(pos,forces,dedv_mean,
     &                      ek0,ekprim,ekvir,presvir)
        use atoms
        use units
        use atmtyp
        use math
        use boxes
        use mdstuf
        use bath
        IMPLICIT NONE
        real*8, intent(in) :: pos(:,:,:),forces(:,:,:),dedv_mean,ek0
        real*8, intent(out) :: ekprim,ekvir,presvir
        real*8, allocatable :: centroid(:,:)
        real*8 :: omp,omp2!,hbar
        integer :: ibead,i,j

        !hbar=(planck*1.d11*avogadro)/(2*pi)

        allocate(centroid(3,n))
        centroid(:,:)=0.d0
        DO ibead=1,nbeads
          DO i=1,n
          if (atomic(i).eq.0) cycle
            DO j=1,3
              centroid(j,i)=centroid(j,i)+pos(j,i,ibead)
            ENDDO
           ENDDO
        ENDDO  
        centroid(:,:)=centroid(:,:)/REAL(nbeads)

        omp=nbeads*boltzmann*kelvin/hbar        
        omp2=omp*omp

c       COMPUTE PRIMITIVE KINETIC ENERGY
        ekprim=0.d0
        DO ibead=1,nbeads-1
          DO i=1,n
            if (atomic(i).eq.0) cycle
            DO j=1,3
              ekprim = ekprim - 0.5*mass(i)*omp2
     &          *(pos(j,i,ibead+1)-pos(j,i,ibead))**2
            ENDDO
          ENDDO
        ENDDO  
        DO i=1,n;
          if (atomic(i).eq.0) cycle
          DO j=1,3
            ekprim = ekprim - 0.5*mass(i)*omp2
     &          *(pos(j,i,nbeads)-pos(j,i,1))**2
          ENDDO
        ENDDO  
        ekprim = (ekprim/nbeads
     &          + 0.5*nbeads*nfree*boltzmann*kelvin)/convert

c       COMPUTE VIRIAL KINETIC ENERGY
        ekvir=0.d0
        DO ibead=1,nbeads
          DO i=1,n
           if (atomic(i).eq.0) cycle
            DO j=1,3
              ekvir=ekvir+(pos(j,i,ibead)-centroid(j,i))
     &                      *forces(j,i,ibead)
            ENDDO
           ENDDO
        ENDDO

        presvir = prescon*( -dedv_mean + (ek0
     &               - ekvir)/(3*nbeads*volbox) )

        ekvir=0.5d0*(nfree*boltzmann*kelvin/convert-ekvir/nbeads)


      end subroutine

      end
