#
#   Makefile for Tinker-HP
#
RunF77 = mpif90
MKLDIR = /opt/intel/mkl/
#FFTDECOMPDIR = $(HOME)/codes/tinker/pmedev/2decomp_fft/
FFTDECOMPDIR = ../2decomp_fft/
FFTDECOMP = -I$(FFTDECOMPDIR)/include -L$(FFTDECOMPDIR)/lib -l2decomp_fft 
BLAS   = -I$(MKLDIR)/include  -L$(MKLDIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread #-lm
##FFLAGS = -O3 -xHost -openmp -no-ipo -no-prec-div -assume cc_omp -inline -heap-arrays#-vec-report=5#-align array64byte
#FFLAGS = -O3 -xHost -openmp -no-ipo -no-prec-div -openmp -inline -heap-arrays #-vec-report=5#-align array64byte
#FFLAGS = -O3 -xHost -openmp -no-ipo -no-prec-div -assume cc_omp -check uninit -check pointers -check bounds -traceback -inline -g -mcmodel=medium -shared-intel 
#FFLAGS2 = -fno-omit-frame-pointer -recursive -g -mcmodel=medium -shared-intel
FFLAGS = -O3 -fopenmp -ffast-math
#FFLAGS = -fopenmp -fbacktrace -fcheck=bounds
LIBS   = libtinker.a 
MODSAM = MOD_sizes.o\
  MOD_action.o MOD_analyz.o MOD_angang.o MOD_angle.o\
  MOD_angpot.o MOD_argue.o  MOD_ascii.o  MOD_atmlst.o\
  MOD_atoms.o  MOD_atmtyp.o MOD_bath.o   MOD_bitor.o\
  MOD_bndpot.o MOD_bond.o   MOD_bound.o  MOD_boxes.o\
  MOD_cell.o   MOD_charge.o MOD_chgpot.o MOD_chunks.o\
  MOD_couple.o MOD_cutoff.o MOD_deriv.o  MOD_divcon.o\
  MOD_domdec.o MOD_energi.o MOD_ewald.o  MOD_fft.o\
  MOD_fields.o MOD_files.o  MOD_freeze.o MOD_group.o\
  MOD_improp.o MOD_imptor.o MOD_inform.o MOD_inter.o\
  MOD_iounit.o MOD_kanang.o MOD_kangs.o  MOD_katoms.o\
  MOD_kbonds.o MOD_kchrge.o MOD_keys.o   MOD_khbond.o\
  MOD_kiprop.o MOD_kgeoms.o MOD_kitors.o MOD_kmulti.o\
  MOD_kopbnd.o MOD_kopdst.o MOD_kpitor.o MOD_kpolr.o\
  MOD_kstbnd.o MOD_ksttor.o MOD_ktorsn.o MOD_ktrtor.o\
  MOD_kurybr.o MOD_kvdwpr.o MOD_kvdws.o  MOD_langevin.o\
  MOD_linmin.o
MODSNZ = MOD_math.o MOD_merck.o  MOD_mpole.o    MOD_mdstuf.o\
  MOD_minima.o MOD_molcul.o   MOD_moldyn.o MOD_mplpot.o\
  MOD_mutant.o MOD_neigh.o    MOD_opbend.o MOD_opdist.o\
  MOD_output.o MOD_params.o   MOD_pitors.o MOD_pme.o\
  MOD_polar.o  MOD_polgrp.o   MOD_polpot.o\
  MOD_potent.o MOD_precis.o   MOD_ptable.o MOD_resdue.o\
  MOD_ring.o   MOD_scales.o   MOD_shunt.o  MOD_strbnd.o\
  MOD_strtor.o MOD_timestat.o MOD_titles.o MOD_torpot.o\
  MOD_tors.o   MOD_tortor.o   MOD_units.o  MOD_uprior.o\
  MOD_urey.o   MOD_urypot.o   MOD_usage.o  MOD_virial.o\
  MOD_vdw.o    MOD_vdwpot.o\

OBJSAM = active.o analysis.o analyze.o angles.o\
   attach.o basefile.o beeman.o bicubic.o\
   baoab.o baoabrespa.o bbk.o\
   bitors.o bonds.o bounds.o calendar.o\
   chkpole.o chkring.o chkxyz.o cholesky.o cluster.o\
   command.o commstuff.o control.o cspline.o\
   cutoffs.o diis.o drivermpi.o dynamic.o\
   dcinduce_pme.o dcinduce_pme2.o\
   eangang.o eangang1.o eangang3.o\
   eangle.o eangle1.o eangle3.o\
   ebond.o ebond1.o ebond3.o\
   echarge.o echarge1.o echarge3.o\
   egeom.o egeom1.o egeom3.o\
   ehal.o ehal1.o ehal3.o\
   eimprop.o eimprop1.o eimprop3.o\
   eimptor.o eimptor1.o eimptor3.o\
   elj.o elj1.o elj3.o\
   empole.o empole0.o empole1.o empole3.o\
   energy.o eopbend.o eopbend1.o eopbend3.o\
   eopdist.o eopdist1.o eopdist3.o\
   epitors.o epitors1.o epitors3.o\
   epolar.o epolar1.o epolar3.o\
   erf.o estrbnd.o estrbnd1.o estrbnd3.o\
   estrtor.o estrtor1.o estrtor3.o\
   etors.o etors1.o etors3.o\
   etortor.o etortor1.o etortor3.o\
   eurey.o eurey1.o eurey3.o\
   evcorr.o extra.o extra1.o extra3.o\
   fatal.o fft_mpi.o field.o\
   final.o freeunit.o geometry.o getkey.o\
   getnumb.o getprm.o getstring.o gettext.o\
   getword.o getxyz.o gradient.o\
   hybrid.o image.o\
   initatom.o initial.o initprm.o initres.o\
   invert.o kangang.o kangle.o\
   katom.o kbond.o kcharge.o kewald.o kgeom.o\
   kimprop.o kimptor.o kinetic.o kmpole.o\
   kopbend.o kopdist.o kpitors.o kpolar.o\
   kstrbnd.o kstrtor.o ktors.o ktortor.o\
   kurey.o kvdw.o lattice.o lbfgs.o\
   linalg.o maxwell.o mdinit.o mdrest.o\
   mdsave.o mdstat.o mechanic.o minimize.o\
   molecule.o mutate.o\

OBJSNZ = nblist.o newinduce_pme.o newinduce_pme2.o\
   nextarg.o nexttext.o nspline.o number.o\
   numeral.o openend.o optsave.o pmestuff.o\
   precise.o pressure.o prime.o prmkey.o promo.o\
   prtdyn.o prtxyz.o random.o readdyn.o\
   readprm.o readxyz.o respa.o\
   rings.o rotpole.o search.o sort.o\
   rattle.o shakeup.o\
   suffix.o switch.o temper.o testgrad.o\
   torphase.o torque.o torsions.o trimtext.o\
   unitcell.o verlet.o version.o

MODS   = $(MODSAM) $(MODSNZ)
OBJS   = $(OBJSAM) $(OBJSNZ)

all: libtinker analyze dynamic testgrad minimize

everything: libtinker analyze dynamic testgrad

analyze: $(MODS) $(OBJS)
	./library.make
	$(RunF77) $(FFLAGS) $(FFLAGS2) -o analyze.x analyze.o $(LIBS) $(BLAS) $(FFTDECOMP)
	strip analyze.x
	mkdir -p ../bin
	mv analyze.x ../bin/analyze
#
dynamic: $(MODS) $(OBJS)
	./library.make
	$(RunF77) $(FFLAGS) $(FFLAGS2) -o dynamic.x dynamic.o $(LIBS) $(BLAS) $(FFTDECOMP)
	strip dynamic.x
	mkdir -p ../bin
	mv dynamic.x ../bin/dynamic
#
testgrad: $(MODS) $(OBJS)
	./library.make
	$(RunF77) $(FFLAGS) $(FFLAGS2) -o testgrad.x testgrad.o $(LIBS) $(BLAS) $(FFTDECOMP)
	strip testgrad.x
	mkdir -p ../bin
	mv testgrad.x ../bin/testgrad
#
minimize: $(MODS) $(OBJS)
	./library.make
	$(RunF77) $(FFLAGS) $(FFLAGS2) -o minimize.x minimize.o $(LIBS) $(BLAS) $(FFTDECOMP)
	strip minimize.x
	mkdir -p ../bin
	mv minimize.x ../bin/minimize
#
libtinker: $(MODS) $(OBJS)
	./library.make
#
%.o: %.f
	$(RunF77) $(FFLAGS) -c $*.f 
%.o: %.F
	$(RunF77) $(FFLAGS) -c $*.F
%.o: %.f90
	$(RunF77) $(FFLAGS) -c $*.f90 $(FFTDECOMP)
#
clean:
	rm -fr $(MODS) $(OBJS) *.exe *.mod
	rm -f libtinker.a 
#
tags: tags-mark
	-'rm' -f tags
	ctags *.f
	-'rm' -f tags-mark

tags-mark:
	touch tags-mark

debug:
	cat *.f > xxx.f
	ftnchek xxx.f
	rm xxx.f
