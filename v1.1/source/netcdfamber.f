c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "netcdfamber" writes out a set of Cartesian coordinates
c     to an external disk file in the amber netcdf format
c
      subroutine netcdfamber(init,istep,dt)
      use atoms
      use atmtyp
      use bound
      use boxes
      use couple
      use files
      use inform
      use netcdf
      use titles
      implicit none
      integer i,j,k,ixyz
      integer size,crdsiz
      integer ncid,status,nf90_create_
      integer dimframe,dimspatial,dimatoms,dimcells,dimcella
      integer dimidchar(1),dimidcellspatial(1),dimidcellangular(1)
      integer dimidframe(1),dimidcoordinates(3)
      integer dimidcelll(2),dimidcella(2)
      integer varid,istep
      real*8 crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 xyzfile
      real, allocatable :: temp(:,:)
      real*8 dt
      real tempcell(3)
      real pico
      integer starttime(1),starttab(3),starttab1(2),counttab(3)
      integer counttab1(2)
      logical init,exist
      character*240 netcdffile
c
      netcdffile = filename(1:leng)//'.nc'
      if (init) then
c  
c     create the netcdf dataset
c
        status = nf90_create(netcdffile,nf90_clobber,ncid)
c
c     create the netcdf dimensions
c
        status = nf90_def_dim(ncid,"frame",NF90_unlimited,dimframe)
        status = nf90_def_dim(ncid,"spatial",3,dimspatial)
        status = nf90_def_dim(ncid,"atom",n,dimatoms)
        status = nf90_def_dim(ncid,"cell_spatial",3,dimcells)
        status = nf90_def_dim(ncid,"cell_angular",3,dimcella)
        dimidchar(1) = dimspatial
        dimidcellspatial(1) = dimcells
        dimidcellangular(1) = dimcella
        dimidframe(1) = dimframe
        dimidcoordinates(1) = dimspatial
        dimidcoordinates(2) = dimatoms
        dimidcoordinates(3) = dimframe
        dimidcelll(1) = dimcells
        dimidcelll(2) = dimframe
        dimidcella(1) = dimcella
        dimidcella(2) = dimframe
c
c     create the netcdf variables
c
        status = nf90_def_var(ncid,"spatial",NF90_CHAR,dimidchar,varid)

        status = nf90_def_var(ncid,"cell_spatial",NF90_CHAR,
     $     dimidcellspatial,varid)

        status = nf90_def_var(ncid,"cell_angular",NF90_CHAR,
     $     dimidcellangular,varid)

        status = nf90_def_var(ncid,"time",NF90_FLOAT,dimidframe,varid)
        status = nf90_put_att(ncid,varid,"units","picosecond")

        status = nf90_def_var(ncid,"coordinates",NF90_FLOAT,
     $     dimidcoordinates,varid)
        status = nf90_put_att(ncid,varid,"units","angstrom")

        status = nf90_def_var(ncid,"cell_lengths",NF90_FLOAT,
     $     dimidcelll,varid)
        status = nf90_put_att(ncid,varid,"units","angstrom")

        status = nf90_def_var(ncid,"cell_angles",NF90_FLOAT,
     $     dimidcella,varid)
        status = nf90_put_att(ncid,varid,"units","degree")

        status = nf90_put_att(ncid,NF90_GLOBAL,"title","amber output")
        status = nf90_put_att(ncid,NF90_GLOBAL,"application","tinker")
        status = nf90_put_att(ncid,NF90_GLOBAL,"program","tinker")
        status = nf90_put_att(ncid,NF90_GLOBAL,"programVersion","8.2.1")
        status = nf90_put_att(ncid,NF90_GLOBAL,"Conventions","AMBER")
        status = nf90_put_att(ncid,NF90_GLOBAL,"ConventionVersion",
     $   "1.0")
         
        status = nf90_enddef(ncid)
      else
c
c     open the netcdf dataset
c
        status = nf90_open(netcdffile,nf90_write,ncid)
      endif
c
c     append the netcdf file
c
      status = nf90_inq_varid(ncid,"time",varid)
      pico = real(istep) * real(dt)
      starttime(1) = istep
      status = nf90_put_var(ncid,varid,pico,starttime)

      status = nf90_inq_varid(ncid,"cell_lengths",varid)
      tempcell(1) = real(xbox)
      tempcell(2) = real(ybox)
      tempcell(3) = real(zbox)
      starttab1(1) = 1
      starttab1(2) = istep 
      counttab1(1) = 3
      counttab1(2) = 1
      status = nf90_put_var(ncid,varid,tempcell,starttab1,counttab1)

      status = nf90_inq_varid(ncid,"coordinates",varid)
      allocate (temp(3,n))
      do i = 1, n
        temp(1,i) = real(x(i))
        temp(2,i) = real(y(i))
        temp(3,i) = real(z(i))
      end do
      starttab(1) = 1 
      starttab(2) = 1 
      starttab(3) = istep 
      counttab(1) = 3
      counttab(2) = n
      counttab(3) = 1
      status = nf90_put_var(ncid,varid,temp,starttab,counttab)
      deallocate (temp)

      status = nf90_inq_varid(ncid,"cell_angles",varid)
      tempcell(1) = real(alpha)
      tempcell(2) = real(beta)
      tempcell(3) = real(gamma) 
      starttab1(1) = 1
      starttab1(2) = istep 
      counttab1(1) = 3
      counttab1(2) = 1
      status = nf90_put_var(ncid,varid,tempcell,starttab1,counttab1)
c
c     close the netcdf dataset
c
      status = nf90_close(ncid)
      return
      end
