c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kmlpo  --  ML potential parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kmlpot" assigns Machine learning potential parameters 
c
c
#include "tinker_precision.h"
      subroutine kmlpot(init,istep)
      use ani
      use atoms
      use bound
      use boxes
      use chunks
      use cutoff
      use domdec
      use ewald
      use fft
      use group
      use inform
      use iounit
      use keys
      use mpi
      use pme
      use potent
      use sizes     ,only: tinkerdebug
      use tinheader
      use tinMemory
      use usage
      implicit none
      integer istep
      logical init
      integer i,j,k
      integer next,gnumml,nloc_cap
      integer iproc,iglob,iloc,idomlen,ibufbeg
      integer nfixed
      integer, allocatable :: fixed(:)
      real(t_p) distcut,d,mcell
      character(:), allocatable :: value
      character*25 word
      character*25 keyword
      character*240 text
      character*240 record
      character*240 record_raw
      character*240 string

      MLpot = "ANI2X"
      use_ani_only = .FALSE.
      use_ml_embedding=.FALSE.
      mlpotscale = 1.0_ti_p
      ml_embedding_mode=1
c
c     parse the line to extract any possible keyword
c
      do i = 1, nkey
         record_raw = keyline(i)
         next = 1
         record = record_raw
         call upcase (record)
         call gettext (record,keyword,next)

         if (keyword(1:8) .eq. 'ANI1CCX ') then
            call getword (record,word,next)
            value = trim(word)
            if (value .eq. 'ONLY') then
               call potoff; use_ani_only=.true.;
            end if
            use_mlpot    = .true.
            MLpot = "ANI1CCX"
            if (value .eq. 'NONE')  use_mlpot = .false.
         else if (keyword(1:6) .eq. 'ANI1X ') then
            call getword (record,word,next)
            value = trim(word)
            if (value .eq. 'ONLY') then
               call potoff; use_ani_only=.true.
            end if
            use_mlpot    = .true.
            MLpot = "ANI1X"
            if (value .eq. 'NONE')  use_mlpot = .false.
         else if (keyword(1:6) .eq. 'ANI2X ') then
            call getword (record,word,next)
            value = trim(word)
            if (value .eq. 'ONLY') then
               call potoff; use_ani_only=.true.
            end if
            use_mlpot    = .true.
            MLpot = "ANI2X"
            if (value .eq. 'NONE')  use_mlpot = .false.
         else if (keyword(1:6) .eq. 'MLPOT ') then
            call getword (record,word,next)
            value = trim(word)
            use_mlpot    = .true.
            if (value .eq. 'ONLY') then
              call potoff; use_ani_only=.true.
            elseif (value .eq. 'EMBEDDING') then
              use_ml_embedding=.true.
            elseif (value .eq. 'NONE') then
              use_mlpot = .false.
            end if
          else if (keyword(1:9) .eq. 'ML-MODEL ') then
            call getword (record,MLpot,next)
            if (trim(MLpot)=="ANI_GENERIC") then
              call getword (record_raw,model_file,next)
            elseif (trim(MLpot)=="DEEPMD") then
              call getword (record_raw,model_file,next)
            endif
         end if
      end do

      if (.not.use_mlpot) return

      if(use_ml_embedding .AND. (.not. use_group)) then
        if(ranktot==0) then
          write(0,*) 'ML embedding requires the group 1 to be defined'
          call fatal  
        endif
      endif

      if (use_ml_embedding) then
         if (rank.eq.0) write(*,*) '--- TinkerHP : USING ML EMBEDDING'
         call fetchkeya('INTEGRATOR',record,'')
         if (index(record,'RESPA1').gt.0) ml_embedding_mode=2
      end if

      call ddpme3d
c
c     search keywords for MLP commands
c
      call prmem_request(iaml,n)
      call prmem_request(laml,n)

      naml    = 0
      laml(:) = .false.
      iaml(:) = 0

      if(use_ml_embedding) then
        naml = igrp(2,1) - igrp(1,1) + 1
        if((naml<=0 .or. igrp(1,i)==0) .and. ranktot==0) then
          write(0,*) 'ML embedding requires the group 1 to be defined'
          call fatal
        endif
        do i=1,naml
          iaml(i)  = kgrp(igrp(1,1)+i-1)
          laml(iaml(i)) = .true.
        enddo

        if(ranktot==0) then
          write(*,*) 'ML embedding attached to group 1'
        endif

      else
        allocate(fixed(n))
        fixed(:) = 0 
        nfixed   = 0

        do i = 1, nkey
           next   = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           string = record(next:240)
           if (keyword(1:6) .eq. 'MLATOM ') then
              read (string,*,err=20,end=20)  (fixed(j),j=nfixed+1,n)
   20         continue
              do while (fixed(nfixed+1) .ne. 0)
                 nfixed = nfixed + 1
                 fixed(nfixed) = max(-n,min(n,fixed(nfixed)))
              end do
           end if
        end do
c
c     set ml atoms to those in fixed
c
        i = 1
        do while (fixed(i) .ne. 0)
           if (fixed(i) .gt. 0) then
              j = fixed(i)
              if (use(j)) then
                 naml        = naml + 1
                 laml(j)     = .true.
                 iaml(naml)  = j
              end if
              i = i + 1
           else
              do j = abs(fixed(i)), abs(fixed(i+1))
                 if (use(j)) then
                    naml       = naml + 1
                    laml(j)    = .true.
                    iaml(naml) = j 
                 end if
              end do
              i = i + 2
           end if
        end do
c
c     Default case
c
        if (i.eq.1) then
           do j = 1,n
              laml(j) = .true.
              iaml(j) = j
           end do
           naml = n
        else
           call sort (naml,iaml)
        end if

      endif

!$acc update device(iaml,laml)

      if (init) then
!$acc enter data copyin(namloc)
      end if

      call kmlpot_reset(istep)

      
      end

      subroutine kmlpot_reset(istep)
      use ani
      use atoms
      use domdec
      use neigh
      use sizes     ,only: tinkerdebug
      use tinMemory
#ifdef _OPENACC
      use thrust
      use utilgpu  ,only: rec_stream
#endif
      implicit none
      integer,intent(in)::istep
      integer i,iglob,naml_cap
      real(t_p) d,distcut
!$acc routine(distprocpart1)

      distcut = mlpotcut+lbuffer

      if (mod(istep,ineigup).ne.0.or.istep.lt.0) return
      call prmem_request(amloc,nbloc,async=.true.)

!$acc serial async present(namloc)
      namloc = 0
!$acc end serial

!$acc parallel loop async default(present)
!$acc&         present(namloc)
      do i = 1,nbloc
         iglob = glob(i)
         if (i.le.nloc) then
            if (laml(iglob)) then
!$acc atomic capture
                 namloc   = namloc + 1
                 naml_cap = namloc
!$acc end atomic
               amloc(naml_cap) = iglob
            end if
         else
            if (laml(iglob)) then
               call distprocpart1(iglob,rank,d,.true.,x,y,z)
               if (d.le.(distcut)) then
!$acc atomic capture
                 namloc   = namloc + 1
                 naml_cap = namloc
!$acc end atomic
                 amloc(naml_cap) = iglob
               end if
            end if
         end if
      end do
!$acc update host(namloc) async

#ifdef _OPENACC
      if (nproc.eq.1) then
!$acc wait
!$acc host_data use_device(amloc)
         call thrust_sort(amloc,namloc,rec_stream)
!$acc end host_data
      end if
#endif

!$acc wait
      if (rank.eq.0.and.tinkerdebug.gt.0) then 
         print*, 'kmlpot_reset',naml,namloc,istep,rank
      end if

      end subroutine
