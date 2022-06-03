ifeq ($(arch),$(filter $(arch),device gpu))
echargecu.o: MOD_utilcu.o pair_charge.f.inc echargecu.tpl.f
echgljcu.o: MOD_utilcu.o MOD_utilgpu.o MOD_vdw.o echgljcu.tpl.f pair_lj.inc.f pair_charge.f.inc
eljcu.o: MOD_utilcu.o MOD_utilgpu.o pair_lj.inc.f eljcu.tpl.f
eliaison1cu.o: MOD_utilcu.o MOD_utilgpu.o MOD_vdw.o
ehal1cu.o: MOD_utilcu.o MOD_utilgpu.o MOD_vdw.o ehalcu.tpl.f pair_ehal.inc.f
nblistcu.o: MOD_utilcu.o MOD_utilgpu.o
tmatxb_pmecu.o: MOD_utilcu.o MOD_utilgpu.o pair_tmatxb.f.inc
epolar1cu.o: MOD_utilcu.o MOD_utilgpu.o pair_polar.f.inc
empole1cu.o: MOD_utilcu.o MOD_utilgpu.o empolecu.tpl.f pair_mpole1.f.inc
pmestuffcu.o: MOD_utilcu.o MOD_utilgpu.o

cu_CholeskySolver.o: utils.h
cu_nblist.o: utils.h image.h
cu_tmatxb_pme.o: utils.h image.h
cu_mpole1.o: utils.h image.h

elj1gpu.o: eljcu.o
elj3gpu.o: eljcu.o
eliaison1gpu.o: MOD_utilcu.o eliaison1cu.o
ehal1gpu.o: ehal1cu.o
ehal3gpu.o: ehal1cu.o
echarge1gpu.o: echargecu.o echgljcu.o
echarge3gpu.o: echargecu.o
empole1gpu.o: empole1cu.o
empole3gpu.o: empole1cu.o
epolar1gpu.o: epolar1cu.o
epolar3gpu.o: epolar1cu.o
kscalfactor.o: MOD_utilcu.o
kvdw.o: ehal1cu.o
lattice.o: pmestuffcu.o MOD_utilcu.o
nblistgpu.o: nblistcu.o
nblist_build.gpu.o: nblistcu.o
efld0_directgpu.o: tmatxb_pmecu.o
tmatxb_pmegpu.o: MOD_utilcomm.o tmatxb_pmecu.o
pmestuffgpu.o: pmestuffcu.o
endif

MOD_subAtoms.o: MOD_atoms.o MOD_utilgpu.o
MOD_subMemory.o: MOD_memory.o MOD_atoms.o MOD_vdw.o
MOD_subDeriv.o: MOD_deriv.o MOD_atoms.o MOD_domdec.o MOD_neigh.o MOD_utilgpu.o
MOD_subInform.o: MOD_inform.o MOD_neigh.o MOD_usage.o MOD_utilgpu.o MOD_vdw.o
MOD_precompute_polegpu.o: MOD_atoms.o MOD_atmlst.o MOD_couple.o MOD_chgpot.o MOD_domdec.o MOD_inform.o MOD_mpole.o MOD_mplpot.o MOD_neigh.o MOD_polar.o MOD_utilgpu.o MOD_sizes.o MOD_vdwpot.o
active.o:
analysis.o: MOD_vdw.o
angles.o: MOD_utilgpu.o
attach.o: MOD_utilgpu.o
basefile.o:
beeman.o:
bicubic.o:
baoab.o:
baoabrespa.o: MOD_utilgpu.o
baoabrespa1.o: MOD_utilgpu.o
bbk.o:
bitors.o: MOD_utilgpu.o
bonds.o: MOD_neigh.o MOD_vdw.o MOD_utilgpu.o
bounds.o: MOD_neigh.o
calendar.o:
chkpole.o:
chkpolegpu.o: MOD_utilgpu.o
chkring.o:
chkxyz.o:
cholesky.o: MOD_utilgpu.o
cluster.o: MOD_neigh.o MOD_utilgpu.o
colvars.o: MOD_atoms.o MOD_atmtyp.o MOD_bath.o MOD_boxes.o MOD_colvars.o MOD_deriv.o MOD_domdec.o MOD_energi.o MOD_files.o MOD_iounit.o MOD_mutant.o MOD_potent.o
command.o:
control.o:
cspline.o:
cutoffs.o: MOD_neigh.o
diis.o:
domdecstuff.o:
dcdio.o:
dcinduce_pme.o: MOD_neigh.o MOD_pme.o
dcinduce_pme2.o: MOD_neigh.o MOD_pme.o
dcinduce_pmegpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
dcinduce_pme2gpu.o: MOD_pme.o MOD_utilgpu.o
dcinduce_shortreal.o: MOD_neigh.o MOD_pme.o
dcinduce_shortrealgpu.o: MOD_neigh.o MOD_pme.o
domdecstuff.o: MOD_neigh.o
eamd1.o: MOD_mamd.o MOD_utilgpu.o
eangang.o:
eangang1.o:
eangang3.o:
eangle.o:
eangle1.o:
eangle1gpu.o:
eangle3.o:
eangle3gpu.o :
eangtor.o:
eangtor1.o:
eangtor3.o:
ebond.o:
ebond1.o:
ebond1gpu.o:
ebond3.o:
ebond3gpu.o:
eliaison1gpu.o:
echarge.o: MOD_neigh.o MOD_pme.o
echarge1.o: MOD_neigh.o MOD_pme.o
echarge1gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o MOD_vdw.o
echarge3.o: MOD_neigh.o MOD_pme.o
echarge3gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
efld0_direct.o: MOD_neigh.o
efld0_directgpu.o: MOD_neigh.o MOD_utilgpu.o
egeom.o:
egeom1.o:
egeom1gpu.o: MOD_utilgpu.o
egeom3.o:
egeom3gpu.o:
ehal1.o: MOD_neigh.o MOD_vdw.o
ehal3.o: MOD_neigh.o MOD_vdw.o
ehal1gpu.o: MOD_neigh.o MOD_utilgpu.o MOD_vdw.o
ehal3gpu.o: MOD_neigh.o MOD_vdw.o
eimprop.o:
eimprop1.o:
eimprop1gpu.o:
eimprop3.o:
eimptor.o:
eimptor1.o:
eimptor1gpu.o:
eimptor3.o:
elj.o: MOD_neigh.o MOD_vdw.o
elj1.o: MOD_neigh.o MOD_vdw.o
elj3.o: MOD_neigh.o MOD_vdw.o
elj1gpu.o: MOD_neigh.o MOD_utilgpu.o MOD_vdw.o
elj3gpu.o: MOD_neigh.o MOD_utilgpu.o MOD_vdw.o
empole0.o: MOD_neigh.o MOD_pme.o
empole1.o: MOD_neigh.o MOD_pme.o
empole1gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
empole3.o: MOD_neigh.o MOD_pme.o
empole3gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
energy.o: MOD_utilgpu.o
epolar.o: MOD_neigh.o MOD_pme.o
epolar1.o: MOD_neigh.o MOD_pme.o
epolar1gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
epolar3.o: MOD_neigh.o MOD_pme.o
epolar3gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
esmd1.o:
estrbnd.o:
estrbnd1.o:
estrbnd1gpu.o:
estrbnd3.o:
estrtor.o:
estrtor1.o:
estrtor1gpu.o:
estrtor3.o:
etors.o:
etors1.o:
etors1gpu.o:
etors3.o:
etortor.o:
etortor1.o:
etortor1gpu.o:
etortor3.o:
eurey.o:
eurey1.o:
eurey1gpu.o :
eurey3.o:
eurey3gpu.o:
evcorr.o: MOD_vdw.o MOD_utilgpu.o
extra.o:
extra1.o:
extra3.o:
fatal.o: MOD_pme.o
fft_mpi.o:
field.o:
final.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o MOD_vdw.o
freeunit.o:
geometry.o:
getkey.o:
getnumb.o:
getprm.o:
getstring.o:
gettext.o:
getword.o:
getxyz.o:
gradient.o: MOD_utilgpu.o
image.o:
initatom.o:
initial.o: MOD_neigh.o MOD_utilgpu.o
initprm.o:
initres.o:
invert.o:
kamd.o: MOD_mamd.o
kangtor.o: MOD_utilgpu.o
kcharge.o: MOD_neigh.o MOD_utilgpu.o
kchgflx.o:
kchgtrn.o: MOD_pme.o
kdisp.o: MOD_pme.o
kewald.o: MOD_pme.o
kgeom.o: MOD_vdw.o MOD_utilgpu.o
kimprop.o: MOD_utilgpu.o
kimptor.o:
kmpole.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
kopbend.o: MOD_utilgpu.o
kpitors.o: MOD_utilgpu.o
kpolar.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
kscalfactor.o: MOD_neigh.o MOD_utilgpu.o MOD_vdw.o
kstrbnd.o: MOD_utilgpu.o
kstrtor.o: MOD_utilgpu.o
ktortor.o: MOD_utilgpu.o
kurey.o: MOD_utilgpu.o
kvdw.o: MOD_neigh.o MOD_utilgpu.o MOD_vdw.o
lattice.o: MOD_pme.o
linalg.o: MOD_utilgpu.o
mdinit.o: MOD_neigh.o MOD_utilgpu.o
mechanic.o: MOD_neigh.o MOD_utilgpu.o MOD_vdw.o
minimize.o: MOD_utilgpu.o
mpistuff.o: MOD_neigh.o MOD_pme.o MOD_utilcomm.o MOD_utilgpu.o
nblist.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o MOD_vdw.o
nblistgpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o MOD_vdw.o
nblist_vdw.gpu.o: MOD_neigh.o MOD_utilgpu.o
newinduce_pmegpu.o: MOD_pme.o MOD_utilgpu.o
newinduce_pme2gpu.o: MOD_pme.o MOD_utilgpu.o
newinduce_shortreal.o: MOD_neigh.o MOD_pme.o MOD_utilcomm.o
newinduce_shortrealgpu.o: MOD_pme.o MOD_utilcomm.o MOD_utilgpu.o
orthogonalize.o: MOD_orthogonalize.o MOD_utilgpu.o
pmestuff.o: MOD_pme.o MOD_utilgpu.o
pmestuffgpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
prmkey.o: MOD_neigh.o MOD_vdw.o
respa.o:MOD_utilgpu.o
respa1.o: MOD_utilgpu.o
rotpolegpu.o: MOD_pme.o MOD_utilgpu.o
scalders.o: MOD_neigh.o
tcgstuff.o: MOD_neigh.o
tmatxb_pme.o: MOD_neigh.o
tmatxb_pmegpu.o: MOD_neigh.o MOD_utilgpu.o
torphase.o:
torque.o:
torquegpu.o: MOD_utilgpu.o
torsions.o: MOD_utilgpu.o
trimtext.o:
scalders.o:
unitcell.o:
verlet.o: MOD_utilgpu.o
version.o:
