c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readprm  --  input of force field parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readprm" processes the potential energy parameter file
c     in order to define the default force field parameters
c
c
#include "tinker_precision.h"
      subroutine readprm
      use fields
      use iounit
      use kanang
      use kangs
      use kantor
      use katoms
      use kbonds
      use kcflux
      use kchrge
      use kcpen
      use kctrn
      use kdsp
      use kiprop
      use kitors
      use khbond
      use kmulti
      use kopbnd
      use kopdst
      use kpitor
      use kpolr
      use krepl
      use kstbnd
      use ksttor
      use ktorsn
      use ktrtor
      use kurybr
      use kvdwpr
      use kvdws
      use merck
      use params
      use tinheader
      use utils
      implicit none
      integer i,j,iprm
      integer ia,ib,ic,id,ie
      integer if,ig,ih,ii
      integer size,next
      integer length,trimtext
      integer nb,nb5,nb4,nb3,nel,nbm,nbm4
      integer na,na5,na4,na3,nap,naf,naps
      integer nups,nuq
      integer nsb,nu,nopb,nopd
      integer ndi,nti,nt,nt5,nt4
      integer npt,nbt,nat,ntt,nd,nd5
      integer nd4,nd3,nvp,nhb,nmp
      integer npi,npi5,npi4
      integer cls,atn,lig
      integer nx,ny,nxy
      integer bt,at,sbt,tt
      integer ncfa,ncfb
      integer(8) res
      integer ft(6),pg(maxvalue)
      real(r_p) wght
      real(t_p) rd,ep,rdn
      real(t_p) an1,an2,an3
      real(t_p) ba1,ba2
      real(t_p) aa1,aa2,aa3
      real(t_p) bt1,bt2,bt3
      real(t_p) bt4,bt5,bt6
      real(t_p) bt7,bt8,bt9
      real(t_p) at1,at2,at3
      real(t_p) at4,at5,at6
      real(t_p) an,pr,ds,dk
      real(t_p) vd,cg
      real(t_p) fc,bd,dl
      real(t_p) pt,pol,thl
      real(t_p) abc,cba
      real(t_p) gi,alphi
      real(t_p) nni,factor
      real(t_p) vt(6),st(6)
      real(t_p) pl(13)
      real(t_p) emtp1,emtp2,emtp3
      real(t_p) tx(maxtgrd2)
      real(t_p) ty(maxtgrd2)
      real(t_p) tf(maxtgrd2)
      real(t_p) spr,apr,epr
      real(t_p) ctrn,atrn
      real(t_p) cfb,cfb1,cfb2
      real(t_p) cfa1,cfa2
      real(t_p) cdp,adp,ddp
      real(t_p) pal,pel

      logical header,swap
      character*1 da1
      character*4 pa,pb,pc
      character*4 pd,pe
      character*8 axt
      character*20 keyword
      character*20 text
      character*240 record
      character*240 string
c
c
c     initialize the counters for some parameter types
c
      nvp = 0
      nhb = 0
      nb = 0
      nb5 = 0
      nb4 = 0
      nb3 = 0
      nbm = 0
      nbm4 = 0
      nel = 0
      na = 0
      na5 = 0
      na4 = 0
      na3 = 0
      nap = 0
      naf = 0
      nsb = 0
      nu = 0
      naps = 0
      nups = 0
      nuq = 0
      nopb = 0
      nopd = 0
      ndi = 0
      nti = 0
      nt = 0
      nt5 = 0
      nt4 = 0
      npt = 0
      nbt = 0
      nat = 0
      ntt = 0
      nd = 0
      nd5 = 0
      nd4 = 0
      nd3 = 0
      nmp = 0
      ncfb = 0
      ncfa = 0
      npi = 0
      npi5 = 0
      npi4 = 0
c
c     number of characters in an atom number text string
c
      size = 4
c
c     set blank line header before echoed comment lines
c
      header = .true.
c
c     process each line of the parameter file, first
c     extract the keyword at the start of each line
c
      iprm = 0
      do while (iprm .lt. nprm)
         iprm = iprm + 1
         record = prmline(iprm)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
c
c     check for a force field modification keyword
c
         call prmkey (record)
c
c     comment line to be echoed to the output
c
         if (keyword(1:5) .eq. 'ECHO ') then
            string = record(next:240)
            length = trimtext (string)
            if (header) then
               header = .false.
               write (iout,10)
   10          format ()
            end if
            if (length .eq. 0) then
               write (iout,20)
   20          format ()
            else
               write (iout,30)  string(1:length)
   30          format (a)
            end if
c
c     atom type definitions and parameters
c
         else if (keyword(1:5) .eq. 'ATOM ') then
            ia = 0
            cls = 0
            atn = 0
            wght = 0.0_ti_p
            lig = 0
            call getnumb (record,ia,next)
            call getnumb (record,cls,next)
            if (cls .eq. 0)  cls = ia
            atmcls(ia) = cls
            if (ia .ge. maxtyp) then
               write (iout,40)
   40          format (/,' READPRM  --  Too many Atom Types;',
     &                    ' Increase MAXTYP')
               call fatal
            else if (cls .ge. maxclass) then
               write (iout,50)
   50          format (/,' READPRM  --  Too many Atom Classes;',
     &                    ' Increase MAXCLASS')
               call fatal
            end if
            if (ia .ne. 0) then
               call gettext (record,symbol(ia),next)
               call getstring (record,describe(ia),next)
               string = record(next:240)
               read (string,*,err=60,end=60)  atn,wght,lig
   60          continue
               atmnum(ia) = atn
               weight(ia) = wght
               ligand(ia) = lig
            end if
c
c     van der Waals parameters for individual atom types
c
         else if (keyword(1:4) .eq. 'VDW ') then
            ia = 0
            rd = 0.0_ti_p
            ep = 0.0_ti_p
            rdn = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,rd,ep,rdn
   70       continue
            if (ia .ne. 0) then
               rad(ia) = rd
               eps(ia) = ep
               reduct(ia) = rdn
            end if
c
c     van der Waals 1-4 parameters for individual atom types
c
         else if (keyword(1:6) .eq. 'VDW14 ') then
            ia = 0
            rd = 0.0_ti_p
            ep = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=80,end=80)  ia,rd,ep
   80       continue
            if (ia .ne. 0) then
               rad4(ia) = rd
               eps4(ia) = ep
            end if
c
c     van der Waals parameters for specific atom pairs
c
         else if (keyword(1:6) .eq. 'VDWPR ') then
            ia = 0
            ib = 0
            rd = 0.0_ti_p
            ep = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=90,end=90)  ia,ib,rd,ep
   90       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nvp = nvp + 1
            if (ia .le. ib) then
               kvpr(nvp) = pa//pb
            else
               kvpr(nvp) = pb//pa
            end if
            radpr(nvp) = rd
            epspr(nvp) = ep
c
c     van der Waals parameters for hydrogen bonding pairs
c
         else if (keyword(1:6) .eq. 'HBOND ') then
            ia = 0
            ib = 0
            rd = 0.0_ti_p
            ep = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=100,end=100)  ia,ib,rd,ep
  100       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nhb = nhb + 1
            if (ia .le. ib) then
               khb(nhb) = pa//pb
            else
               khb(nhb) = pb//pa
            end if
            radhb(nhb) = rd
            epshb(nhb) = ep
c
c     bond stretching parameters
c
         else if (keyword(1:5) .eq. 'BOND ') then
            ia = 0
            ib = 0
            fc = 0.0_ti_p
            bd = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=110,end=110)  ia,ib,fc,bd
  110       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb = nb + 1
            if (ia .le. ib) then
               kb(nb) = pa//pb
            else
               kb(nb) = pb//pa
            end if
            bcon(nb) = fc
            blen(nb) = bd
        else if (keyword(1:6) .eq. 'MORSE ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            ep = 2.0d0
            string = record(next:240)
            read (string,*)  ia,ib,fc,bd,ep
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nbm = nbm + 1
            if (ia .le. ib) then
               kbm(nbm) = pa//pb
            else
               kbm(nbm) = pb//pa
            end if
            bmor(1,nbm) = fc
            bmor(2,nbm) = bd
            bmor(3,nbm) = ep
         else if (keyword(1:7) .eq. 'MORSE4 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            ep = 2.0d0
            string = record(next:240)
            read (string,*)  ia,ib,fc,bd,ep
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nbm4 = nbm4 + 1
            if (ia .le. ib) then
               kbm4(nbm4) = pa//pb
            else
               kbm4(nbm4) = pb//pa
            end if
            bmor4(1,nbm4) = fc
            bmor4(2,nbm4) = bd
            bmor4(3,nbm4) = ep
c
c     bond stretching parameters for 5-membered rings
c
         else if (keyword(1:6) .eq. 'BOND5 ') then
            ia = 0
            ib = 0
            fc = 0.0_ti_p
            bd = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=120,end=120)  ia,ib,fc,bd
  120       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb5 = nb5 + 1
            if (ia .le. ib) then
               kb5(nb5) = pa//pb
            else
               kb5(nb5) = pb//pa
            end if
            bcon5(nb5) = fc
            blen5(nb5) = bd
c
c     bond stretching parameters for 4-membered rings
c
         else if (keyword(1:6) .eq. 'BOND4 ') then
            ia = 0
            ib = 0
            fc = 0.0_ti_p
            bd = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=130,end=130)  ia,ib,fc,bd
  130       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb4 = nb4 + 1
            if (ia .le. ib) then
               kb4(nb4) = pa//pb
            else
               kb4(nb4) = pb//pa
            end if
            bcon4(nb4) = fc
            blen4(nb4) = bd
c
c     bond stretching parameters for 3-membered rings
c
         else if (keyword(1:6) .eq. 'BOND3 ') then
            ia = 0
            ib = 0
            fc = 0.0_ti_p
            bd = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=140,end=140)  ia,ib,fc,bd
  140       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb3 = nb3 + 1
            if (ia .le. ib) then
               kb3(nb3) = pa//pb
            else
               kb3(nb3) = pb//pa
            end if
            bcon3(nb3) = fc
            blen3(nb3) = bd
c
c     electronegativity bond length correction parameters
c
         else if (keyword(1:9) .eq. 'ELECTNEG ') then
            ia = 0
            ib = 0
            ic = 0
            dl = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=150,end=150)  ia,ib,ic,dl
  150       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            nel = nel + 1
            if (ia .le. ic) then
               kel(nel) = pa//pb//pc
            else
               kel(nel) = pc//pb//pa
            end if
            dlen(nel) = dl
c
c     bond angle bending parameters
c
         else if (keyword(1:6) .eq. 'ANGLE ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            an1 = 0.0_ti_p
            an2 = 0.0_ti_p
            an3 = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=160,end=160)  ia,ib,ic,fc,an1,an2,an3
  160       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na = na + 1
            if (ia .le. ic) then
               ka(na) = pa//pb//pc
            else
               ka(na) = pc//pb//pa
            end if
            acon(na) = fc
            ang(1,na) = an1
            ang(2,na) = an2
            ang(3,na) = an3
c
c     angle bending parameters for 5-membered rings
c
         else if (keyword(1:7) .eq. 'ANGLE5 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            an1 = 0.0_ti_p
            an2 = 0.0_ti_p
            an3 = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=170,end=170)  ia,ib,ic,fc,an1,an2,an3
  170       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na5 = na5 + 1
            if (ia .le. ic) then
               ka5(na5) = pa//pb//pc
            else
               ka5(na5) = pc//pb//pa
            end if
            acon5(na5) = fc
            ang5(1,na5) = an1
            ang5(2,na5) = an2
            ang5(3,na5) = an3
c
c     angle bending parameters for 4-membered rings
c
         else if (keyword(1:7) .eq. 'ANGLE4 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            an1 = 0.0_ti_p
            an2 = 0.0_ti_p
            an3 = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=180,end=180)  ia,ib,ic,fc,an1,an2,an3
  180       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na4 = na4 + 1
            if (ia .le. ic) then
               ka4(na4) = pa//pb//pc
            else
               ka4(na4) = pc//pb//pa
            end if
            acon4(na4) = fc
            ang4(1,na4) = an1
            ang4(2,na4) = an2
            ang4(3,na4) = an3
c
c     angle bending parameters for 3-membered rings
c
         else if (keyword(1:7) .eq. 'ANGLE3 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            an1 = 0.0_ti_p
            an2 = 0.0_ti_p
            an3 = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=190,end=190)  ia,ib,ic,fc,an1,an2,an3
  190       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na3 = na3 + 1
            if (ia .le. ic) then
               ka3(na3) = pa//pb//pc
            else
               ka3(na3) = pc//pb//pa
            end if
            acon3(na3) = fc
            ang3(1,na3) = an1
            ang3(2,na3) = an2
            ang3(3,na3) = an3
c
c     Fourier bond angle bending parameters
c
         else if (keyword(1:7) .eq. 'ANGLEF ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            an = 0.0_ti_p
            pr = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=200,end=200)  ia,ib,ic,fc,an,pr
  200       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            naf = naf + 1
            if (ia .le. ic) then
               kaf(naf) = pa//pb//pc
            else
               kaf(naf) = pc//pb//pa
            end if
            aconf(naf) = fc
            angf(1,naf) = an
            angf(2,naf) = pr
c
c     Partridge-Schwenke bond angle bending parameters
c
         else if (keyword(1:8) .eq. 'ANGLEPS ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an = 0.0d0
            pr = 0.0d0
            string = record(next:240)
            ! ia ib ic  theta_e  r_e  beta
            read (string,*)  ia,ib,ic,fc,an,pr
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            naps = naps + 1
            if (ia .le. ic) then
               kaps(naps) = pa//pb//pc
            else
               kaps(naps) = pc//pb//pa
            end if
            angps(1,naps) = fc 
            angps(2,naps) = an
            angps(3,naps) = pr
c
c     in-plane projected angle bending parameters
c
         else if (keyword(1:7) .eq. 'ANGLEP ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            string = record(next:240)
            read (string,*,err=800,end=800)  ia,ib,ic,fc,an1,an2
  800       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            nap = nap + 1
            if (ia .le. ic) then
               kap(nap) = pa//pb//pc
            else
               kap(nap) = pc//pb//pa
            end if
            aconp(nap) = fc
            angp(1,nap) = an1
            angp(2,nap) = an2
c
c     stretch-bend parameters
c
         else if (keyword(1:7) .eq. 'STRBND ') then
            ia = 0
            ib = 0
            ic = 0
            ba1 = 0.0_ti_p
            ba2 = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=210,end=210)  ia,ib,ic,ba1,ba2
  210       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
c           call numeral (ic,pc,size)
            nsb = nsb + 1
            if (ia .le. ic) then
               call front_convert_base(ia,ib,ic,ksb(nsb))
c              ksb(nsb) = pa//pb//pc
               stbn(1,nsb) = ba1
               stbn(2,nsb) = ba2
            else
               call front_convert_base(ic,ib,ia,ksb(nsb))
c              ksb(nsb) = pc//pb//pa
               stbn(1,nsb) = ba2
               stbn(2,nsb) = ba1
            end if
c
c     Urey-Bradley parameters
c
         else if (keyword(1:9) .eq. 'UREYBRAD ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            ds = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=220,end=220)  ia,ib,ic,fc,ds
  220       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
c           call numeral (ic,pc,size)
            nu = nu + 1
            if (ia .le. ic) then
               call front_convert_base(ia,ib,ic,ku(nu))
c              ku(nu) = pa//pb//pc
            else
               call front_convert_base(ic,ib,ia,ku(nu))
c              ku(nu) = pc//pb//pa
            end if
            ucon(nu) = fc
            dst13(nu) = ds
c
c     Partridge-Schwenke-type Angle repulsion (HH potential)
c
         else if (keyword(1:7) .eq. 'ANGREP ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            ds = 0.0d0
            string = record(next:240)
            read (string,*,err=221,end=221)  ia,ib,ic,fc,ds
  221       continue
c            call numeral (ia,pa,size)
c            call numeral (ib,pb,size)
c            call numeral (ic,pc,size)
            nups = nups + 1
            if (ia .le. ic) then
               call front_convert_base(ia,ib,ic,kups(nups))
c              ku(nu) = pa//pb//pc
            else
               call front_convert_base(ic,ib,ia,kups(nups))
c              ku(nu) = pc//pb//pa
            end if
            uconps(nups) = fc
            dst13ps(nups) = ds
c
c         Quartic Urey term
c
         else if (keyword(1:12) .eq. 'UREYQUARTIC ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            ds = 0.0d0
            string = record(next:240)
            read (string,*,err=222,end=222)  ia,ib,ic,fc,ds
  222       continue
c            call numeral (ia,pa,size)
c            call numeral (ib,pb,size)
c            call numeral (ic,pc,size)
            nuq = nuq + 1
            if (ia .le. ic) then
               call front_convert_base(ia,ib,ic,kuq(nuq))
c              ku(nu) = pa//pb//pc
            else
               call front_convert_base(ic,ib,ia,kuq(nuq))
c              ku(nu) = pc//pb//pa
            end if
            uconq(nuq) = fc
            dst13q(nuq) = ds
c
c     angle-angle parameters
c
         else if (keyword(1:7) .eq. 'ANGANG ') then
            ia = 0
            aa1 = 0.0_ti_p
            aa2 = 0.0_ti_p
            aa3 = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=230,end=230)  ia,aa1,aa2,aa3
  230       continue
            if (ia .ne. 0) then
               anan(1,ia) = aa1
               anan(2,ia) = aa2
               anan(3,ia) = aa3
            end if
c
c     out-of-plane bend parameters
c
         else if (keyword(1:7) .eq. 'OPBEND ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fc = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=240,end=240)  ia,ib,ic,id,fc
  240       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
c           call numeral (ic,pc,size)
c           call numeral (id,pd,size)
            nopb = nopb + 1
c           kopb(nopb) = pa//pb//pc//pd
            call front_convert_base(0,ia,ib,ic,id,kopb(nopb))
            opbn(nopb) = fc
c
c     out-of-plane distance parameters
c
         else if (keyword(1:7) .eq. 'OPDIST ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fc = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=250,end=250)  ia,ib,ic,id,fc
  250       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nopd = nopd + 1
            kopd(nopd) = pa//pb//pc//pd
            opds(nopd) = fc
c
c     improper dihedral parameters
c
         else if (keyword(1:9) .eq. 'IMPROPER ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            dk = 0.0_ti_p
            vd = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=260,end=260)  ia,ib,ic,id,dk,vd
  260       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
c           call numeral (ic,pc,size)
c           call numeral (id,pd,size)
            ndi = ndi + 1
c           kdi(ndi) = pa//pb//pc//pd
            call front_convert_base(0,ia,ib,ic,id,kdi(ndi))
            dcon(ndi) = dk
            tdi(ndi) = vd
c
c     improper torsional parameters
c
         else if (keyword(1:8) .eq. 'IMPTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0_ti_p
               st(i) = 0.0_ti_p
               ft(i) = 0
            end do
            string = record(next:240)
            read (string,*,err=270,end=270)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  270       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
c           call numeral (ic,pc,size)
c           call numeral (id,pd,size)
            nti = nti + 1
c           kti(nti) = pa//pb//pc//pd
            call torphase (ft,vt,st)
            call front_convert_base(0,ia,ib,ic,id,kti(nti))
            ti1(1,nti) = vt(1)
            ti1(2,nti) = st(1)
            ti2(1,nti) = vt(2)
            ti2(2,nti) = st(2)
            ti3(1,nti) = vt(3)
            ti3(2,nti) = st(3)
c
c     torsional parameters
c
         else if (keyword(1:8) .eq. 'TORSION ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0_ti_p
               st(i) = 0.0_ti_p
               ft(i) = 0
            end do
            string = record(next:240)
            read (string,*,err=280,end=280)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  280       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt = nt + 1
            if (ib .lt. ic) then
               kt(nt) = pa//pb//pc//pd
            else if (ic .lt. ib) then
               kt(nt) = pd//pc//pb//pa
            else if (ia .le. id) then
               kt(nt) = pa//pb//pc//pd
            else if (id .lt. ia) then
               kt(nt) = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            t1(1,nt) = vt(1)
            t1(2,nt) = st(1)
            t2(1,nt) = vt(2)
            t2(2,nt) = st(2)
            t3(1,nt) = vt(3)
            t3(2,nt) = st(3)
            t4(1,nt) = vt(4)
            t4(2,nt) = st(4)
            t5(1,nt) = vt(5)
            t5(2,nt) = st(5)
            t6(1,nt) = vt(6)
            t6(2,nt) = st(6)
c
c     torsional parameters for 5-membered rings
c
         else if (keyword(1:9) .eq. 'TORSION5 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0_ti_p
               st(i) = 0.0_ti_p
               ft(i) = 0
            end do
            string = record(next:240)
            read (string,*,err=290,end=290)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  290       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt5 = nt5 + 1
            if (ib .lt. ic) then
               kt5(nt5) = pa//pb//pc//pd
            else if (ic .lt. ib) then
               kt5(nt5) = pd//pc//pb//pa
            else if (ia .le. id) then
               kt5(nt5) = pa//pb//pc//pd
            else if (id .lt. ia) then
               kt5(nt5) = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            t15(1,nt5) = vt(1)
            t15(2,nt5) = st(1)
            t25(1,nt5) = vt(2)
            t25(2,nt5) = st(2)
            t35(1,nt5) = vt(3)
            t35(2,nt5) = st(3)
            t45(1,nt5) = vt(4)
            t45(2,nt5) = st(4)
            t55(1,nt5) = vt(5)
            t55(2,nt5) = st(5)
            t65(1,nt5) = vt(6)
            t65(2,nt5) = st(6)
c
c     torsional parameters for 4-membered rings
c
         else if (keyword(1:9) .eq. 'TORSION4 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0_ti_p
               st(i) = 0.0_ti_p
               ft(i) = 0
            end do
            string = record(next:240)
            read (string,*,err=300,end=300)  ia,ib,ic,id,
     &                                       (vt(i),st(i),ft(i),i=1,6)
  300       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt4 = nt4 + 1
            if (ib .lt. ic) then
               kt4(nt4) = pa//pb//pc//pd
            else if (ic .lt. ib) then
               kt4(nt4) = pd//pc//pb//pa
            else if (ia .le. id) then
               kt4(nt4) = pa//pb//pc//pd
            else if (id .lt. ia) then
               kt4(nt4) = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            t14(1,nt4) = vt(1)
            t14(2,nt4) = st(1)
            t24(1,nt4) = vt(2)
            t24(2,nt4) = st(2)
            t34(1,nt4) = vt(3)
            t34(2,nt4) = st(3)
            t44(1,nt4) = vt(4)
            t44(2,nt4) = st(4)
            t54(1,nt4) = vt(5)
            t54(2,nt4) = st(5)
            t64(1,nt4) = vt(6)
            t64(2,nt4) = st(6)
c
c     pi-orbital torsion parameters
c
         else if (keyword(1:7) .eq. 'PITORS ') then
            ia = 0
            ib = 0
            pt = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=310,end=310)  ia,ib,pt
  310       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
            npt = npt + 1
            if (ia .le. ib) then
               call front_convert_base(0,ia,ib,kpt(npt))
c              kpt(npt) = pa//pb
            else
               call front_convert_base(0,ib,ia,kpt(npt))
c              kpt(npt) = pb//pa
            end if
            ptcon(npt) = pt
c
c     stretch-torsion parameters
c
         else if (keyword(1:8) .eq. 'STRTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            bt1 = 0.0_ti_p
            bt2 = 0.0_ti_p
            bt3 = 0.0_ti_p
            bt4 = 0.0_ti_p
            bt5 = 0.0_ti_p
            bt6 = 0.0_ti_p
            bt7 = 0.0_ti_p
            bt8 = 0.0_ti_p
            bt9 = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=320,end=320)  ia,ib,ic,id,bt1,bt2,bt3,
     &                                       bt4,bt5,bt6,bt7,bt8,bt9
  320       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
c           call numeral (ic,pc,size)
c           call numeral (id,pd,size)
            nbt = nbt + 1
            if (ib .lt. ic) then
               call front_convert_base(0,ia,ib,ic,id,kbt(nbt))
               swap = .false.
            else if (ic .lt. ib) then
               call front_convert_base(0,id,ic,ib,ia,kbt(nbt))
c              kbt(nbt) = pd//pc//pb//pa
               swap = .true.
            else if (ia .le. id) then
               call front_convert_base(0,ia,ib,ic,id,kbt(nbt))
c              kbt(nbt) = pa//pb//pc//pd
               swap = .false.
            else if (id .lt. ia) then
               call front_convert_base(0,id,ic,ib,ia,kbt(nbt))
c              kbt(nbt) = pd//pc//pb//pa
               swap = .true.
            end if
            btcon(4,nbt) = bt4
            btcon(5,nbt) = bt5
            btcon(6,nbt) = bt6
            if (swap) then
               btcon(1,nbt) = bt7
               btcon(2,nbt) = bt8
               btcon(3,nbt) = bt9
               btcon(7,nbt) = bt1
               btcon(8,nbt) = bt2
               btcon(9,nbt) = bt3
            else
               btcon(1,nbt) = bt1
               btcon(2,nbt) = bt2
               btcon(3,nbt) = bt3
               btcon(7,nbt) = bt7
               btcon(8,nbt) = bt8
               btcon(9,nbt) = bt9
            end if
c
c     angle-torsion parameters
c
         else if (keyword(1:8) .eq. 'ANGTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            at1 = 0.0d0
            at2 = 0.0d0
            at3 = 0.0d0
            at4 = 0.0d0
            at5 = 0.0d0
            at6 = 0.0d0
            string = record(next:240)
            read (string,*,err=330,end=330)  ia,ib,ic,id,at1,at2,
     &                                       at3,at4,at5,at6
  330       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
c           call numeral (ic,pc,size)
c           call numeral (id,pd,size)
            nat = nat + 1
            if (ib .lt. ic) then
               call front_convert_base(0,ia,ib,ic,id,kat(nat))
               !kat(nat) = pa//pb//pc//pd
               swap = .false.
            else if (ic .lt. ib) then
               call front_convert_base(0,id,ic,ib,ia,kat(nat))
               !kat(nat) = pd//pc//pb//pa
               swap = .true.
            else if (ia .le. id) then
               call front_convert_base(0,ia,ib,ic,id,kat(nat))
               !kat(nat) = pa//pb//pc//pd
               swap = .false.
            else if (id .lt. ia) then
               call front_convert_base(0,id,ic,ib,ia,kat(nat))
               !kat(nat) = pd//pc//pb//pa
               swap = .true.
            end if
            if (swap) then
               atcon(1,nat) = at4
               atcon(2,nat) = at5
               atcon(3,nat) = at6
               atcon(4,nat) = at1
               atcon(5,nat) = at2
               atcon(6,nat) = at3
            else
               atcon(1,nat) = at1
               atcon(2,nat) = at2
               atcon(3,nat) = at3
               atcon(4,nat) = at4
               atcon(5,nat) = at5
               atcon(6,nat) = at6
            end if
c
c     torsion-torsion parameters
c
         else if (keyword(1:8) .eq. 'TORTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            ie = 0
            nx = 0
            ny = 0
            nxy = 0
            do i = 1, maxtgrd2
               tx(i) = 0.0_ti_p
               ty(i) = 0.0_ti_p
               tf(i) = 0.0_ti_p
            end do
            string = record(next:240)
            read (string,*,err=340,end=340)  ia,ib,ic,id,ie,nx,ny
            nxy = nx * ny
            do i = 1, nxy
               iprm = iprm + 1
               record = prmline(iprm)
               read (record,*,err=340,end=340)  tx(i),ty(i),tf(i)
            end do
  340       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
c           call numeral (ic,pc,size)
c           call numeral (id,pd,size)
c           call numeral (ie,pe,size)
            ntt = ntt + 1
c           ktt(ntt) = pa//pb//pc//pd//pe
            call front_convert_base(ia,ib,ic,id,ie,ktt(ntt))
            nx = nxy
            call sort9 (nx,tx)
            ny = nxy
            call sort9 (ny,ty)
            tnx(ntt) = nx
            tny(ntt) = ny
            do i = 1, nx
               ttx(i,ntt) = tx(i)
            end do
            do i = 1, ny
               tty(i,ntt) = ty(i)
            end do
            do i = 1, nxy
               tbf(i,ntt) = tf(i)
            end do
c
c     atomic partial charge parameters
c
         else if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=350,end=350)  ia,cg
  350       continue
            if (ia .ne. 0)  chg(ia) = cg
cc
cc     bond dipole moment parameters
cc
c         else if (keyword(1:7) .eq. 'DIPOLE ') then
c            ia = 0
c            ib = 0
c            dp = 0.0_ti_p
c            ps = 0.5_ti_p
c            string = record(next:240)
c            read (string,*,err=360,end=360)  ia,ib,dp,ps
c  360       continue
c            call numeral (ia,pa,size)
c            call numeral (ib,pb,size)
c            nd = nd + 1
c            if (ia .le. ib) then
c               kd(nd) = pa//pb
c            else
c               kd(nd) = pb//pa
c            end if
c            dpl(nd) = dp
c            pos(nd) = ps
cc
cc     bond dipole moment parameters for 5-membered rings
cc
c         else if (keyword(1:8) .eq. 'DIPOLE5 ') then
c            ia = 0
c            ib = 0
c            dp = 0.0_ti_p
c            ps = 0.5_ti_p
c            string = record(next:240)
c            read (string,*,err=360,end=360)  ia,ib,dp,ps
c  360       continue
c            call numeral (ia,pa,size)
c            call numeral (ib,pb,size)
c            nd5 = nd5 + 1
c            if (ia .le. ib) then
c               kd5(nd5) = pa//pb
c            else
c               kd5(nd5) = pb//pa
c            end if
c            dpl5(nd5) = dp
c            pos5(nd5) = ps
cc
cc     bond dipole moment parameters for 4-membered rings
cc
c         else if (keyword(1:8) .eq. 'DIPOLE4 ') then
c            ia = 0
c            ib = 0
c            dp = 0.0_ti_p
c            ps = 0.5_ti_p
c            string = record(next:240)
c            read (string,*,err=370,end=370)  ia,ib,dp,ps
c  370       continue
c            call numeral (ia,pa,size)
c            call numeral (ib,pb,size)
c            nd4 = nd4 + 1
c            if (ia .le. ib) then
c               kd4(nd4) = pa//pb
c            else
c               kd4(nd4) = pb//pa
c            end if
c            dpl4(nd4) = dp
c            pos4(nd4) = ps
cc
cc     bond dipole moment parameters for 3-membered rings
cc
c         else if (keyword(1:8) .eq. 'DIPOLE3 ') then
c            ia = 0
c            ib = 0
c            dp = 0.0_ti_p
c            ps = 0.5_ti_p
c            string = record(next:240)
c            read (string,*,err=380,end=380)  ia,ib,dp,ps
c  380       continue
c            call numeral (ia,pa,size)
c            call numeral (ib,pb,size)
c            nd3 = nd3 + 1
c            if (ia .le. ib) then
c               kd3(nd3) = pa//pb
c            else
c               kd3(nd3) = pb//pa
c            end if
c            dpl3(nd3) = dp
c            pos3(nd3) = ps
c
c     atomic multipole moment parameters
c
         else if (keyword(1:10) .eq. 'MULTIPOLE ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            axt = 'Z-then-X'
            do i = 1, 13
               pl(i) = 0.0_ti_p
            end do
            string = record(next:240)
            read (string,*,err=390,end=390)  ia,ib,ic,id,pl(1)
            goto 420
  390       continue
            id = 0
            read (string,*,err=400,end=400)  ia,ib,ic,pl(1)
            goto 420
  400       continue
            ic = 0
            read (string,*,err=410,end=410)  ia,ib,pl(1)
            goto 420
  410       continue
            ib = 0
            read (string,*,err=430,end=430)  ia,pl(1)
  420       continue
            iprm = iprm + 1
            record = prmline(iprm)
            read (record,*,err=430,end=430)  pl(2),pl(3),pl(4)
            iprm = iprm + 1
            record = prmline(iprm)
            read (record,*,err=430,end=430)  pl(5)
            iprm = iprm + 1
            record = prmline(iprm)
            read (record,*,err=430,end=430)  pl(8),pl(9)
            iprm = iprm + 1
            record = prmline(iprm)
            read (record,*,err=430,end=430)  pl(11),pl(12),pl(13)
  430       continue
            if (ib .eq. 0)  axt = 'None'
            if (ib.ne.0 .and. ic.eq.0)  axt = 'Z-Only'
            if (ib.lt.0 .or. ic.lt.0)  axt = 'Bisector'
            if (ic.lt.0 .and. id.lt.0)  axt = 'Z-Bisect'
            if (max(ib,ic,id) .lt. 0)  axt = '3-Fold'
            ib = abs(ib)
            ic = abs(ic)
            id = abs(id)
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nmp = nmp + 1
            kmp(nmp) = pa//pb//pc//pd
            mpaxis(nmp) = axt
            multip(1,nmp) = pl(1)
            multip(2,nmp) = pl(2)
            multip(3,nmp) = pl(3)
            multip(4,nmp) = pl(4)
            multip(5,nmp) = pl(5)
            multip(6,nmp) = pl(8)
            multip(7,nmp) = pl(11)
            multip(8,nmp) = pl(8)
            multip(9,nmp) = pl(9)
            multip(10,nmp) = pl(12)
            multip(11,nmp) = pl(11)
            multip(12,nmp) = pl(12)
            multip(13,nmp) = pl(13)
c
c     atomic emtp formula parameters
c
         else if (keyword(1:8) .eq. 'SIBFACP ') then
            ia = 0
            emtp1 = 0.0_ti_p
            emtp2 = 0.0_ti_p
            emtp3 = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=720,end=720)  ia,emtp1,emtp2,emtp3
  720       continue
            if (ia.ne.0) then
              sibfacp(1,ia) = emtp1
              sibfacp(2,ia) = emtp2
              sibfacp(3,ia) = emtp3
            end if
c
c     atomic dipole polarizability parameters
c
         else if (keyword(1:9) .eq. 'POLARIZE ') then
            ia = 0
            pol = 0.0_ti_p
            thl = 0.0_ti_p
            ddp = 0.0_ti_p
            do i = 1, maxvalue
               pg(i) = 0
            end do
            string = record(next:240)
            next   = 1
            call getnumb (string,ia,next)
            call gettext (string,text,next)
            read (text,*,err=440,end=440)  pol
            call gettext (string,text,next)
            i = 1
            call getnumb (text,pg(1),i)
            if (pg(1) .eq. 0) then
               read (text,*,err=440,end=440)  thl
               call gettext (string,text,next)
               i = 1
               call getnumb (text,pg(1),i)
               string = string(next:240)
               if (pg(1) .eq. 0) then
                  read (text,*,err=440,end=440)  ddp
                  read (string,*,err=440,end=440)  (pg(i),i=1,maxvalue)
               else
                  read (string,*,err=440,end=440)  (pg(i),i=2,maxvalue)
               end if
            else
               string = string(next:240)
               read (string,*,err=440,end=440)  (pg(i),i=2,maxvalue)
            end if
  440       continue
            if (ia .ne. 0) then
               polr(ia) = pol
               athl(ia) = thl
               ddir(ia) = ddp
               do i = 1, maxvalue
                  pgrp(i,ia) = pg(i)
               end do
            end if
cc
cc     conjugated pisystem atom parameters
cc
c         else if (keyword(1:7) .eq. 'PIATOM ') then
c            ia = 0
c            el = 0.0_ti_p
c            iz = 0.0_ti_p
c            rp = 0.0_ti_p
c            string = record(next:240)
c            read (string,*,err=450,end=450)  ia,el,iz,rp
c  450       continue
c            if (ia .ne. 0) then
c               electron(ia) = el
c               ionize(ia) = iz
c               repulse(ia) = rp
c            end if
cc
cc     conjugated pisystem bond parameters
cc
c         else if (keyword(1:7) .eq. 'PIBOND ') then
c            ia = 0
c            ib = 0
c            ss = 0.0_ti_p
c            ts = 0.0_ti_p
c            string = record(next:240)
c            read (string,*,err=460,end=460)  ia,ib,ss,ts
c  460       continue
c            call numeral (ia,pa,size)
c            call numeral (ib,pb,size)
c            npi = npi + 1
c            if (ia .le. ib) then
c               kpi(npi) = pa//pb
c            else
c               kpi(npi) = pb//pa
c            end if
c            sslope(npi) = ss
c            tslope(npi) = ts
cc
cc     conjugated pisystem bond parameters for 5-membered rings
cc
c         else if (keyword(1:8) .eq. 'PIBOND5 ') then
c            ia = 0
c            ib = 0
c            ss = 0.0_ti_p
c            ts = 0.0_ti_p
c            string = record(next:240)
c            read (string,*,err=470,end=470)  ia,ib,ss,ts
c  470       continue
c            call numeral (ia,pa,size)
c            call numeral (ib,pb,size)
c            npi5 = npi5 + 1
c            if (ia .le. ib) then
c               kpi5(npi5) = pa//pb
c            else
c               kpi5(npi5) = pb//pa
c            end if
c            sslope5(npi5) = ss
cc            tslope5(npi5) = ts
ccc
cc     conjugated pisystem bond parameters for 4-membered rings
cc
c         else if (keyword(1:8) .eq. 'PIBOND4 ') then
c            ia = 0
c            ib = 0
c            ss = 0.0_ti_p
c            ts = 0.0_ti_p
c            string = record(next:240)
c            read (string,*,err=480,end=480)  ia,ib,ss,ts
c  480       continue
c            call numeral (ia,pa,size)
c            call numeral (ib,pb,size)
c            npi4 = npi4 + 1
c            if (ia .le. ib) then
c               kpi4(npi4) = pa//pb
c            else
c               kpi4(npi4) = pb//pa
c            end if
c            sslope4(npi4) = ss
c            tslope4(npi4) = ts
c
c     metal ligand field splitting parameters
c
         else if (keyword(1:6) .eq. 'METAL ') then
            string = record(next:240)
            read (string,*,err=490,end=490)  ia
  490       continue
c
c     biopolymer atom type conversion definitions
c
         else if (keyword(1:8) .eq. 'BIOTYPE ') then
            ia = 0
            ib = 0
            string = record(next:240)
            read (string,*,err=500,end=500)  ia
            call getword (record,string,next)
            call getstring (record,string,next)
            string = record(next:240)
            read (string,*,err=500,end=500)  ib
  500       continue
            if (ia .ge. maxbio) then
               write (iout,40)
c 510          format (/,' READPRM  --  Too many Biopolymer Types;',
c    &                    ' Increase MAXBIO')
               call fatal
            end if
            if (ia .ne. 0)  biotyp(ia) = ib
c
c     MMFF van der Waals parameters
c
         else if (keyword(1:8) .eq. 'MMFFVDW ') then
            ia = 0
            rd = 0.0_ti_p
            ep = 0.0_ti_p
            rdn = 0.0_ti_p
            da1 = 'C'
            string = record(next:240)
            read (string,*,err=520,end=520)  ia,rd,alphi,nni,gi,da1
  520       continue
            if (ia .ne. 0) then
               rad(ia) = rd
               g(ia) = gi
               alph(ia) = alphi
               nn(ia) = nni
               da(ia) = da1
            end if
c
c     MMFF bond stretching parameters
c
         else if (keyword(1:9) .eq. 'MMFFBOND ') then
            ia = 0
            ib = 0
            fc = 0.0_ti_p
            bd = 0.0_ti_p
            bt = 2
            string = record(next:240)
            read (string,*,err=530,end=530)  ia,ib,fc,bd,bt
  530       continue
            nb = nb + 1
            if (bt .eq. 0) then
               mmff_kb(ia,ib) = fc
               mmff_kb(ib,ia) = fc
               mmff_b0(ia,ib) = bd
               mmff_b0(ib,ia) = bd
            else if (bt .eq. 1) then
               mmff_kb1(ia,ib) = fc
               mmff_kb1(ib,ia) = fc
               mmff_b1(ia,ib) = bd
               mmff_b1(ib,ia) = bd
            end if
c
c     MMFF bond stretching empirical rule parameters
c
         else if (keyword(1:11) .eq. 'MMFFBONDER ') then
            ia = 0
            ib = 0
            fc = 0.0_ti_p
            bd = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=540,end=540)  ia,ib,fc,bd
  540       continue
            r0ref(ia,ib) = fc
            r0ref(ib,ia) = fc
            kbref(ia,ib) = bd
            kbref(ib,ia) = bd
c
c     MMFF bond angle bending parameters
c
         else if (keyword(1:10) .eq. 'MMFFANGLE ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0_ti_p
            an1 = 0.0_ti_p
            at = 3
            string = record(next:240)
            read (string,*,err=550,end=550)  ia,ib,ic,fc,an1,at
  550       continue
            na = na + 1
            if (an1 .ne. 0.0_ti_p) then
               if (at .eq. 0) then
                  mmff_ka(ia,ib,ic) = fc
                  mmff_ka(ic,ib,ia) = fc
                  mmff_ang0(ia,ib,ic) = an1
                  mmff_ang0(ic,ib,ia) = an1
               else if (at .eq. 1) then
                  mmff_ka1(ia,ib,ic) = fc
                  mmff_ka1(ic,ib,ia) = fc
                  mmff_ang1(ia,ib,ic) = an1
                  mmff_ang1(ic,ib,ia) = an1
               else if (at .eq. 2) then
                  mmff_ka2(ia,ib,ic) = fc
                  mmff_ka2(ic,ib,ia) = fc
                  mmff_ang2(ia,ib,ic) = an1
                  mmff_ang2(ic,ib,ia) = an1
               else if (at .eq. 3) then
                  mmff_ka3(ia,ib,ic) = fc
                  mmff_ka3(ic,ib,ia) = fc
                  mmff_ang3(ia,ib,ic) = an1
                  mmff_ang3(ic,ib,ia) = an1
               else if (at .eq. 4) then
                  mmff_ka4(ia,ib,ic) = fc
                  mmff_ka4(ic,ib,ia) = fc
                  mmff_ang4(ia,ib,ic) = an1
                  mmff_ang4(ic,ib,ia) = an1
               else if (at .eq. 5) then
                  mmff_ka5(ia,ib,ic) = fc
                  mmff_ka5(ic,ib,ia) = fc
                  mmff_ang5(ia,ib,ic) = an1
                  mmff_ang5(ic,ib,ia) = an1
               else if (at .eq. 6) then
                  mmff_ka6(ia,ib,ic) = fc
                  mmff_ka6(ic,ib,ia) = fc
                  mmff_ang6(ia,ib,ic) = an1
                  mmff_ang6(ic,ib,ia) = an1
               else if (at .eq. 7) then
                  mmff_ka7(ia,ib,ic) = fc
                  mmff_ka7(ic,ib,ia) = fc
                  mmff_ang7(ia,ib,ic) = an1
                  mmff_ang7(ic,ib,ia) = an1
               else if (at .eq. 8) then
                  mmff_ka8(ia,ib,ic) = fc
                  mmff_ka8(ic,ib,ia) = fc
                  mmff_ang8(ia,ib,ic) = an1
                  mmff_ang8(ic,ib,ia) = an1
               end if
            end if
c
c     MMFF stretch-bend parameters
c
         else if (keyword(1:11) .eq. 'MMFFSTRBND ') then
            ia = 0
            ib = 0
            ic = 0
            abc = 0.0_ti_p
            cba = 0.0_ti_p
            sbt = 4
            string = record(next:240)
            read (string,*,err=560,end=560)  ia,ib,ic,abc,cba,sbt
  560       continue
            if (ia .ne. 0) then
               if (sbt .eq. 0) then
                  stbn_abc(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc(ic,ib,ia) = cba
                  stbn_cba(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba(ic,ib,ia) = abc
               else if (sbt .eq. 1) then
                  stbn_abc1(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc1(ic,ib,ia) = cba
                  stbn_cba1(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba1(ic,ib,ia) = abc
               else if (sbt .eq. 2) then
                  stbn_abc2(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc2(ic,ib,ia) = cba
                  stbn_cba2(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba2(ic,ib,ia) = abc
               else if (sbt .eq. 3) then
                  stbn_abc3(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc3(ic,ib,ia) = cba
                  stbn_cba3(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba3(ic,ib,ia) = abc
               else if (sbt .eq. 4) then
                  stbn_abc4(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc4(ic,ib,ia) = cba
                  stbn_cba4(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba4(ic,ib,ia) = abc
               else if (sbt .eq. 5) then
                  stbn_abc5(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc5(ic,ib,ia) = cba
                  stbn_cba5(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba5(ic,ib,ia) = abc
               else if (sbt .eq. 6) then
                  stbn_abc6(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc6(ic,ib,ia) = cba
                  stbn_cba6(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba6(ic,ib,ia) = abc
               else if (sbt .eq. 7) then
                  stbn_abc7(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc7(ic,ib,ia) = cba
                  stbn_cba7(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba7(ic,ib,ia) = abc
               else if (sbt .eq. 8) then
                  stbn_abc8(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc8(ic,ib,ia) = cba
                  stbn_cba8(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba8(ic,ib,ia) = abc
               else if (sbt .eq. 9) then
                  stbn_abc9(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc9(ic,ib,ia) = cba
                  stbn_cba9(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba9(ic,ib,ia) = abc
               else if (sbt .eq. 10) then
                  stbn_abc10(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc10(ic,ib,ia) = cba
                  stbn_cba10(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba10(ic,ib,ia) = abc
               else if (sbt .eq. 11) then
                  stbn_abc11(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc11(ic,ib,ia) = cba
                  stbn_cba11(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba11(ic,ib,ia) = abc
               end if
            end if
c
c     MMFF out-of-plane bend parameters
c
         else if (keyword(1:11) .eq. 'MMFFOPBEND ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fc = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=570,end=570)  ia,ib,ic,id,fc
  570       continue
c           call numeral (ia,pa,size)
c           call numeral (ib,pb,size)
c           call numeral (ic,pc,size)
c           call numeral (id,pd,size)
            nopb = nopb + 1
            if (ic .le. id) then
c              kopb(nopb) = pa//pb//pc//pd
               call front_convert_base(0,ia,ib,ic,id,kopb(nopb))
            else
c              kopb(nopb) = pa//pb//pd//pc
               call front_convert_base(0,ia,ib,id,ic,kopb(nopb))
            end if
            opbn(nopb) = fc
c           if (ic.gt.0 .or. id.gt.0) then
c              nopb = nopb + 1
c              if (ib .le. id) then
c                 kopb(nopb) = pc//pb//pb//pd
c              else
c                 kopb(nopb) = pc//pb//pd//pb
c              end if
c              opbn(nopb) = fc
c              nopb = nopb + 1
c              if (ia .le. ic) then
c                 kopb(nopb) = pd//pb//pa//pc
c              else
c                 kopb(nopb) = pd//pb//pc//pa
c              end if
c              opbn(nopb) = fc
c           end if
c
c     MMFF torsional parameters
c
         else if (keyword(1:12) .eq. 'MMFFTORSION ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0_ti_p
               st(i) = 0.0_ti_p
               ft(i) = 0
            end do
            tt = 3
            string = record(next:240)
            read (string,*,err=580,end=580)  ia,ib,ic,id,(vt(j),
     &                                       st(j),ft(j),j=1,3),tt
  580       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt = nt + 1
            if (tt .eq. 0) then
               if (ib .lt. ic) then
                  kt(nt) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt(nt) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt(nt) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt(nt) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t1(1,nt) = vt(1)
               t1(2,nt) = st(1)
               t2(1,nt) = vt(2)
               t2(2,nt) = st(2)
               t3(1,nt) = vt(3)
               t3(2,nt) = st(3)
            else if (tt .eq. 1) then
               if (ib .lt. ic) then
                  kt_1(nt) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt_1(nt) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt_1(nt) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt_1(nt) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t1_1(1,nt) = vt(1)
               t1_1(2,nt) = st(1)
               t2_1(1,nt) = vt(2)
               t2_1(2,nt) = st(2)
               t3_1(1,nt) = vt(3)
               t3_1(2,nt) = st(3)
            else if (tt .eq. 2) then
               if (ib .lt. ic) then
                  kt_2(nt) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt_2(nt) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt_2(nt) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt_2(nt) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t1_2(1,nt) = vt(1)
               t1_2(2,nt) = st(1)
               t2_2(1,nt) = vt(2)
               t2_2(2,nt) = st(2)
               t3_2(1,nt) = vt(3)
               t3_2(2,nt) = st(3)
            else if (tt .eq. 4) then
               nt4 = nt4 + 1
               if (ib .lt. ic) then
                  kt4(nt4) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt4(nt4) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt4(nt4) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt4(nt4) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t14(1,nt4) = vt(1)
               t14(2,nt4) = st(1)
               t24(1,nt4) = vt(2)
               t24(2,nt4) = st(2)
               t34(1,nt4) = vt(3)
               t34(2,nt4) = st(3)
            else if (tt .eq. 5) then
               nt5 = nt5 + 1
               if (ib .lt. ic) then
                  kt5(nt5) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt5(nt5) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt5(nt5) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt5(nt5) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t15(1,nt5) = vt(1)
               t15(2,nt5) = st(1)
               t25(1,nt5) = vt(2)
               t25(2,nt5) = st(2)
               t35(1,nt5) = vt(3)
               t35(2,nt5) = st(3)
            end if
c
c     MMFF bond charge increment parameters
c
         else if (keyword(1:8) .eq. 'MMFFBCI ') then
            ia = 0
            ib = 0
            cg = 1000.0_ti_p
            bt = 2
            string = record(next:240)
            read (string,*,err=590,end=590)  ia,ib,cg,bt
  590       continue
            if (ia .ne. 0) then
               if (bt .eq. 0) then
                  bci(ia,ib) = cg
                  bci(ib,ia) = -cg
               else if (bt .eq. 1) then
                  bci_1(ia,ib) = cg
                  bci_1(ib,ia) = -cg
               end if
            end if
c
c     MMFF partial bond charge increment parameters
c
         else if (keyword(1:9) .eq. 'MMFFPBCI ') then
            ia = 0
            string = record(next:240)
            read (string,*,err=600,end=600)  ia,cg,factor
  600       continue
            if (ia .ne. 0) then
               pbci(ia) = cg
               fcadj(ia) = factor
            end if
c
c     MMFF atom class equivalency parameters
c
         else if (keyword(1:10) .eq. 'MMFFEQUIV ') then
            string = record(next:240)
            ia = 1000
            ib = 1000
            ic = 1000
            id = 1000
            ie = 1000
            if = 0
            read (string,*,err=610,end=610)  ia,ib,ic,id,ie,if
  610       continue
            eqclass(if,1) = ia
            eqclass(if,2) = ib
            eqclass(if,3) = ic
            eqclass(if,4) = id
            eqclass(if,5) = ie
c
c     MMFF default stretch-bend parameters
c
         else if (keyword(1:12) .eq. 'MMFFDEFSTBN ') then
            string = record(next:240)
            ia = 1000
            ib = 1000
            ic = 1000
            abc = 0.0_ti_p
            cba = 0.0_ti_p
            read (string,*,err=620,end=620)  ia,ib,ic,abc,cba
  620       continue
            defstbn_abc(ia,ib,ic) = abc
            defstbn_cba(ia,ib,ic) = cba
            defstbn_abc(ic,ib,ia) = cba
            defstbn_cba(ic,ib,ia) = abc
c
c     MMFF covalent radius and electronegativity parameters
c
         else if (keyword(1:11) .eq. 'MMFFCOVRAD ') then
            ia = 0
            fc = 0.0_ti_p
            bd = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=630,end=630)  ia,fc,bd
  630       continue
            rad0(ia) = fc
            paulel(ia) = bd
c
c     MMFF property parameters
c
         else if (keyword(1:9) .eq. 'MMFFPROP ') then
            string = record(next:240)
            ia = 1000
            ib = 1000
            ic = 1000
            id = 1000
            ie = 1000
            if = 1000
            ig = 1000
            ih = 1000
            ii = 1000
            read (string,*,err=640,end=640)  ia,ib,ic,id,ie,
     &                                       if,ig,ih,ii
  640       continue
            crd(ia) = ic
            val(ia) = id
            pilp(ia) = ie
            mltb(ia) = if
            arom(ia) = ig
            lin(ia) = ih
            sbmb(ia) = ii
c
c     MMFF aromatic ion parameters
c
         else if (keyword(1:9) .eq. 'MMFFAROM ') then
            string = record(next:240)
            read (string,*,err=650,end=650)  ia,ib,ic,id,ie,if
  650       continue
            if (ie.eq.0 .and. id.eq.0) then
               mmffarom(ia,if) = ic
            else if (id .eq. 1) then
               mmffaromc(ia,if) = ic
            else if (ie .eq. 1) then
               mmffaroma(ia,if) = ic
            end if
c
c     Pauli repulsion parameters
c
         else if (keyword(1:10) .eq. 'REPULSION ') then
            ia = 0
            spr = 0.0_ti_p
            apr = 0.0_ti_p
            epr = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=660,end=660)  ia,spr,apr,epr
  660       continue
            if (ia .ne. 0) then
               prsiz(ia) = spr
               prdmp(ia) = apr
               prele(ia) = -abs(epr)
            end if
c
c     charge transfer parameters
c
         else if (keyword(1:7) .eq. 'CHGTRN ') then
            ia = 0
            ctrn = 0.0_ti_p
            atrn = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=670,end=670)  ia,ctrn,atrn
  670       continue
            if (ia .ne. 0) then
               ctchg(ia) = ctrn
               ctdmp(ia) = atrn
            end if
c
c     bond charge flux parameters
c
         else if (keyword(1:9) .eq. 'BNDCFLUX ') then
            ia = 0
            ib = 0
            cfb = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=680,end=680)  ia,ib,cfb
  680       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            ncfb = ncfb + 1
            if (ia .lt. ib) then
               kcfb(ncfb) = pa//pb
               cflb(ncfb) = cfb
            else if (ib .lt. ia) then
               kcfb(ncfb) = pb//pa
               cflb(ncfb) = -cfb
            else
               kcfb(ncfb) = pa//pb
               cflb(ncfb) = 0.0d0
            end if
c
c     angle charge flux parameters
c
         else if (keyword(1:9) .eq. 'ANGCFLUX ') then
            ia = 0
            ib = 0
            ic = 0
            cfa1 = 0.0_ti_p
            cfa2 = 0.0_ti_p
            cfb1 = 0.0_ti_p
            cfb2 = 0.0_ti_p
            string = record(next:240)
            read (string,*,err=690,end=690)  ia,ib,ic,cfa1,cfa2,
     &                                       cfb1,cfb2
  690       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            ncfa = ncfa + 1
            if (ia .le. ic) then
               kcfa(ncfa) = pa//pb//pc
            else
               kcfa(ncfa) = pc//pb//pa
            end if
            cfla(1,ncfa) = cfa1
            cfla(2,ncfa) = cfa2
            cflab(1,ncfa) = cfb1
            cflab(2,ncfa) = cfb2
c
c     charge penetration parameters
c
         else if (keyword(1:7) .eq. 'CHGPEN ') then
            ia = 0
            pel = 0.0d0
            pal = 0.0d0
            string = record(next:240)
            read (string,*,err=700,end=700)  ia,pel,pal
  700       continue
            if (ia .ne. 0) then
               cpele(ia) = abs(pel)
               cpalp(ia) = pal
            end if
c
c     damped dispersion parameters
c
         else if (keyword(1:11) .eq. 'DISPERSION ') then
            ia = 0
            cdp = 0.0d0
            adp = 0.0d0
            string = record(next:240)
            read (string,*,err=710,end=710)  ia,cdp,adp
  710       continue
            if (ia .ne. 0) then
               dspsix(ia) = cdp
               dspdmp(ia) = adp
            end if
         end if
      end do
c     call sleep(1)
      return
      end
