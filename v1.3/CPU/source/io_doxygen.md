# I/O and Tinker-HP 

**I/O of trajectories** works as in Tinker 8: if no keyword is specified
then frames will printed individually at the selected frequency. By
adding the keyword **archive** then the frames will be directly appended
in an \*.arc file.

Moreover, Tinker-HP also supports **CHARMM dcd format**: adding the line
**dcdio** to the keyfile will print a trajectory with the dcd format.
When running the analyze or the bar postprocessing programs, if
**dcdio** is specified in the key file, then the programs will look for
the associated \*dcd file and process the frames it contains.


Other keywords control the writing of some quantities during a molecular dynamics trajectory:
  - SAVE-VELOCITY: print velocities in a separate file
  - SAVE-FORCES: print forces in a separate file
  - SAVE-INDUCED: print induced dipoles in a separate file
  - PBC-UNWRAP: print unwrappped trajectory, without this keywords atoms are wrapped by molecules


Tinker-HP deals with restart files for dynamic trajectories the same way
as Tinker-8 does by creating a \*.dyn file encompassing current
positions, velocities and accelerations of the system.
