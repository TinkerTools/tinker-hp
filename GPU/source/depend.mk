MOD_beads.o: MOD_dcdio.o
MOD_couple.o: MOD_sizes.o
MOD_commstuffpi.o: MOD_beads.o MOD_qtb.o
MOD_fft.o: MOD_sizes.o
MOD_fields.o: MOD_sizes.o
MOD_files.o: MOD_iounit.o
MOD_group.o: MOD_sizes.o
MOD_kanang.o: MOD_sizes.o
MOD_katoms.o: MOD_sizes.o
MOD_kchrge.o: MOD_sizes.o
MOD_kcpen.o: MOD_sizes.o
MOD_krepl.o: MOD_sizes.o
MOD_params.o: MOD_sizes.o
MOD_qtb.o: MOD_atmtyp.o MOD_atoms.o MOD_bath.o MOD_domdec.o MOD_langevin.o\
           MOD_math.o MOD_moldyn.o MOD_sizes.o MOD_usage.o MOD_units.o MOD_random.o\
           MOD_utilgpu.o MOD_subAtoms.o MOD_mdstuf.o MOD_spectra.o MOD_beads.o
MOD_spectra.o: MOD_atmtyp.o MOD_atoms.o MOD_bath.o MOD_domdec.o MOD_langevin.o\
           MOD_math.o MOD_moldyn.o MOD_sizes.o MOD_usage.o MOD_units.o MOD_random.o\
           MOD_utilgpu.o MOD_subAtoms.o MOD_mdstuf.o MOD_boxes.o MOD_charge.o MOD_mpole.o MOD_atmlst.o\
           MOD_utils.o MOD_uprior.o
MOD_utilbaoab.o:
MOD_utilbaoabpi.o: MOD_beads.o MOD_deriv.o MOD_langevin.o MOD_utils.o MOD_cell.o \
					 MOD_utilgpu.o MOD_couple.o MOD_commstuffpi.o MOD_spectra.o MOD_qtb.o
MOD_utilgpu.o: MOD_couple.o MOD_cell.o MOD_polgrp.o MOD_polpot.o MOD_sizes.o MOD_vdwpot.o
MOD_utils.o: MOD_sizes.o
MOD_utilvec.o: MOD_polgrp.o MOD_polpot.o MOD_vdwpot.o
ifeq ($(arch),$(filter $(arch),device gpu))
echargecu.o: MOD_utilcu.o MOD_utilgpu.o pair_charge.f.inc echargecu.tpl.f
echgljcu.o: MOD_utilcu.o MOD_utilgpu.o MOD_vdw.o echgljcu.tpl.f pair_lj.inc.f pair_charge.f.inc
echgtrncu.o: MOD_utilcu.o MOD_utilgpu.o echgtrncu.tpl.f pair_chgtrn.inc.f
efld0_cpencu.o: MOD_utilcu.o MOD_utilgpu.o pair_efld_cp.inc.f
eljcu.o: MOD_utilcu.o MOD_utilgpu.o pair_lj.inc.f eljcu.tpl.f
eliaison1cu.o: MOD_utilcu.o MOD_utilgpu.o MOD_vdw.o
ehal1cu.o: MOD_utilcu.o MOD_utilgpu.o MOD_vdw.o ehalcu.tpl.f pair_ehal.inc.f
nblistcu.o: MOD_utilcu.o MOD_utilgpu.o
tmatxb_pmecu.o: MOD_utilcu.o MOD_utilgpu.o pair_tmatxb.f.inc
tmatxb_pme_cpen.cu.o: MOD_utilcu.o MOD_utilgpu.o pair_tmatxb.f.inc
empole1cu.o: MOD_utilcu.o MOD_utilgpu.o empolecu.tpl.f pair_mpole1.f.inc
empole_cpencu.o: MOD_utilcu.o MOD_utilgpu.o empole_cpencu.tpl.f pair_mpole1.f.inc
epolar1cu.o: MOD_utilcu.o MOD_utilgpu.o pair_polar.inc.f
epolar_cpencu.o: MOD_utilcu.o MOD_utilgpu.o pair_polar_cpen.inc.f epolar_cpencu.tpl.f
pmestuffcu.o: MOD_utilcu.o MOD_utilgpu.o

cu_CholeskySolver.o: utils.h
cu_nblist.o: utils.h image.h
cu_tmatxb_pme.o: utils.h image.h damping.h
cu_mpole1.o: utils.h image.h

cluster.o: MOD_utilcu.o
elj1gpu.o: eljcu.o
elj3gpu.o: eljcu.o
eliaison1gpu.o: MOD_utilcu.o eliaison1cu.o
ehal1gpu.o: ehal1cu.o
ehal3gpu.o: ehal1cu.o
echarge1gpu.o: echargecu.o echgljcu.o
echarge3gpu.o: echargecu.o
echgtrn1gpu.o: echgtrncu.o
echgtrn3gpu.o: echgtrncu.o
empole1gpu.o: empole1cu.o
empole3gpu.o: empole1cu.o
epolar1gpu.o: epolar1cu.o
epolar3gpu.o: epolar1cu.o
kscalfactor.o: MOD_utilcu.o
kvdw.o: ehal1cu.o
lattice.o: pmestuffcu.o MOD_utilcu.o
nblistgpu.o: nblistcu.o
nblist_build.gpu.o: nblistcu.o
efld0_directgpu.o: tmatxb_pmecu.o efld0_cpencu.o
tmatxb_pmegpu.o: MOD_utilcomm.o tmatxb_pmecu.o tmatxb_pme_cpen.cu.o
pmestuffgpu.o: pmestuffcu.o
endif
MOD_vec.o: MOD_sizes.o MOD_couple.o MOD_polgrp.o
MOD_subAtoms.o: MOD_atoms.o MOD_utilgpu.o
MOD_subMemory.o: MOD_memory.o MOD_atoms.o MOD_vdw.o
MOD_subDeriv.o: MOD_deriv.o MOD_atoms.o MOD_domdec.o MOD_neigh.o MOD_utilgpu.o MOD_beads.o
MOD_subInform.o: MOD_inform.o MOD_neigh.o MOD_usage.o MOD_utilgpu.o MOD_vdw.o
MOD_precompute_polegpu.o: MOD_atoms.o MOD_atmlst.o MOD_couple.o MOD_chgpot.o MOD_domdec.o MOD_inform.o MOD_mpole.o MOD_mplpot.o MOD_neigh.o MOD_polar.o MOD_utilgpu.o MOD_sizes.o MOD_vdwpot.o
active.o:
alterchg.o: image.f.inc MOD_atoms.o MOD_atmlst.o MOD_atmlst.o MOD_inform.o MOD_mplpot.o MOD_potent.o
alterchg.gpu.o: image.f.inc MOD_atoms.o MOD_atmlst.o MOD_atmlst.o MOD_inform.o MOD_mplpot.o MOD_potent.o MOD_utilcomm.o
analysis.o: MOD_vdw.o
angles.o: MOD_utilgpu.o
ani.o: MOD_ani.o MOD_atoms.o MOD_bath.o MOD_domdec.o
attach.o: MOD_utilgpu.o
basefile.o:
beeman.o:
bicubic.o:
baoab.o:
baoab_util.o: MOD_atoms.o MOD_atmtyp.o MOD_bath.o MOD_domdec.o MOD_freeze.o MOD_moldyn.o MOD_usage.o
baoabpi.o: MOD_utilbaoabpi.o
baoabpiston.o: MOD_utilgpu.o
baoabrespa.o: MOD_ani.o MOD_utilgpu.o MOD_utilbaoab.o
baoabrespapi.o: MOD_utilbaoabpi.o
baoabrespa1.o: MOD_ani.o MOD_utilgpu.o MOD_utilbaoab.o
bbk.o:
beads.o: MOD_beads.o MOD_atoms.o MOD_bitor.o MOD_bond.o MOD_boxes.o MOD_charge.o MOD_cutoff.o MOD_dcdio.o MOD_domdec.o MOD_energi.o MOD_freeze.o MOD_inform.o MOD_improp.o MOD_imptor.o\
					MOD_math.o MOD_mdstuf.o MOD_molcul.o MOD_moldyn.o MOD_mpole.o MOD_neigh.o\
					MOD_opbend.o MOD_opdist.o MOD_pitors.o MOD_potent.o MOD_polar.o MOD_strbnd.o MOD_tors.o\
					MOD_strtor.o MOD_timestat.o MOD_tortor.o MOD_units.o MOD_uprior.o MOD_urey.o MOD_usage.o MOD_utilgpu.o MOD_utils.o MOD_vdw.o MOD_virial.o
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
cutoffs.o: MOD_ani.o MOD_neigh.o
dcflux.o:
dcflux.gpu.o: MOD_cflux.o
diis.o:
domdecstuff.o: MOD_ani.o
dcdio.o: MOD_dcdio.o MOD_subAtoms.o MOD_boxes.o MOD_files.o MOD_inform.o  MOD_iounit.o\
					MOD_domdec.o MOD_atoms.o
dcinduce_pme.o: MOD_neigh.o MOD_pme.o
dcinduce_pme2.o: MOD_neigh.o MOD_pme.o
dcinduce_pmegpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
dcinduce_pme2gpu.o: MOD_pme.o MOD_utilgpu.o
dcinduce_shortreal.o: MOD_neigh.o MOD_pme.o
dcinduce_shortrealgpu.o: MOD_neigh.o MOD_pme.o
dynamic.o: MOD_ani.o
domdecstuff.o: MOD_neigh.o
eamd1.o: MOD_mamd.o MOD_utilgpu.o
eangang.o:
eangang1.o:
eangang3.o:
eangle.o:
eangle1.o: ker_angle.inc.f
eangle1gpu.o: ker_angle.inc.f
eangle3.o:
eangle3gpu.o: ker_angle.inc.f
eangtor.o: ker_angtor.inc.f
eangtor1.o: ker_angtor.inc.f
eangtor3.o: ker_angtor.inc.f
ebond.o: ker_bond.inc.f
ebond1.o:
ebond1gpu.o: ker_bond.inc.f
ebond3.o:
ebond3gpu.o: ker_bond.inc.f
eliaison1gpu.o:
echarge.o: MOD_neigh.o MOD_pme.o
echarge1.o: MOD_neigh.o MOD_pme.o
echarge1gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o MOD_vdw.o
echarge3.o: MOD_neigh.o MOD_pme.o
echarge3gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
efld0_direct.o: MOD_neigh.o
efld0_directgpu.o: MOD_neigh.o MOD_utilgpu.o pair_efld.inc.f
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
energy.o: MOD_ani.o MOD_utilgpu.o
eopbend.o: ker_opbend.inc.f
eopbend1gpu.o: ker_opbend.inc.f
eopbend3gpu.o: ker_opbend.inc.f
epolar.o: MOD_neigh.o MOD_pme.o
epolar1.o: MOD_neigh.o MOD_pme.o
epolar1gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
epolar3.o: MOD_neigh.o MOD_pme.o
epolar3gpu.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
esmd1.o:
estrbnd.o: ker_strbnd.inc.f
estrbnd1.o:
estrbnd1gpu.o: ker_strbnd.inc.f
estrbnd3.o: ker_strbnd.inc.f
estrtor.o:
estrtor1.o:
estrtor1gpu.o:
estrtor3.o:
etors.o: ker_tors.inc.f
etors1.o:
etors1gpu.o: ker_tors.inc.f
etors3.o: ker_tors.inc.f
etortor.o:
etortor1.o:
etortor1gpu.o:
etortor3.o:
eurey.o: ker_urey.inc.f
eurey1.o:
eurey1gpu.o: ker_urey.inc.f
eurey3.o:
eurey3gpu.o: ker_urey.inc.f
evcorr.o: MOD_vdw.o MOD_utilgpu.o
extra.o:
extra1.o:
extra3.o:
fatal.o: MOD_pme.o
fft_mpi.o:
field.o:
final.o: MOD_ani.o MOD_neigh.o MOD_pme.o MOD_utilgpu.o MOD_vdw.o
freeunit.o:
geometry.o:
getkey.o:
getnumb.o:
getprm.o:
getstring.o:
gettext.o:
getword.o:
getxyz.o:
gradient.o: MOD_ani.o MOD_utilgpu.o
image.o:
initatom.o:
initial.o: MOD_beads.o MOD_neigh.o MOD_utilgpu.o
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
kmlpot.o: MOD_ani.o
kopbend.o: MOD_utilgpu.o
kpitors.o: MOD_utilgpu.o
kpolar.o: MOD_neigh.o MOD_pme.o MOD_utilgpu.o
kscalfactor.o: MOD_neigh.o MOD_utilgpu.o MOD_vdw.o
kstrbnd.o: MOD_utilgpu.o
kstrtor.o: MOD_utilgpu.o
ktortor.o: MOD_utilgpu.o
kurey.o: MOD_utilgpu.o MOD_urey.o MOD_urypot.o MOD_kurybr.o
kvdw.o: MOD_neigh.o MOD_utilgpu.o MOD_vdw.o
lattice.o: MOD_pme.o
linalg.o: MOD_utilgpu.o
maxwell.o: MOD_atoms.o MOD_domdec.o MOD_langevin.o MOD_tinheader.o MOD_random.o MOD_units.o
mdinit.o: MOD_beads.o MOD_neigh.o MOD_qtb.o MOD_utilgpu.o
mdinitbead.o:  MOD_beads.o
mdsavebeads.o: MOD_beads.o
mdstatpi.o: MOD_beads.o
mechanic.o: MOD_neigh.o MOD_sizes.o MOD_utilgpu.o MOD_vdw.o
minimize.o: MOD_utilgpu.o
molecule.o: MOD_sizes.o
mpistuff.o: MOD_neigh.o MOD_pme.o MOD_sizes.o MOD_utilcomm.o MOD_utilgpu.o
mpiUtils.gpu.o: MOD_domdec.o MOD_utilcomm.o
mutate.o: MOD_sizes.o
nblist.o: MOD_ani.o MOD_neigh.o MOD_pme.o MOD_utilgpu.o MOD_vdw.o
nblistgpu.o: MOD_neigh.o MOD_pme.o MOD_sizes.o MOD_utilgpu.o MOD_vdw.o
nblist_vdw.gpu.o: MOD_neigh.o MOD_utilgpu.o
newinduce_pmegpu.o: MOD_pme.o MOD_utilgpu.o
newinduce_pme2gpu.o: MOD_pme.o MOD_sizes.o MOD_utilgpu.o
newinduce_shortreal.o: MOD_neigh.o MOD_pme.o MOD_utilcomm.o
newinduce_shortrealgpu.o: MOD_pme.o MOD_sizes.o MOD_utilcomm.o MOD_utilgpu.o
orthogonalize.o: MOD_orthogonalize.o MOD_utilgpu.o
pimd.o: MOD_beads.o MOD_commstuffpi.o
pmestuff.o: MOD_pme.o MOD_sizes.o MOD_utilgpu.o
pmestuffgpu.o: MOD_neigh.o MOD_sizes.o MOD_pme.o MOD_utilgpu.o
prmkey.o: MOD_neigh.o MOD_vdw.o
prtdynbeads.o: MOD_beads.o
respa.o: MOD_ani.o MOD_utilgpu.o
respa1.o: MOD_ani.o MOD_utilgpu.o
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
verlet.o: MOD_ani.o MOD_utilgpu.o
version.o:
