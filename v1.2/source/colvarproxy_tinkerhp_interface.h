// Test linking Colvars through a C interface

#ifdef __cplusplus
extern "C" {
#endif

// C functions exposed by Colvars to Tinker

// This will be callable from C or Fortran
// Returns proxy pointer and stores it to static global unique_colvarproxy_object
void *allocate_colvars();

/// Called from Tinker: run the Colvars computation
void compute_colvars();

/// Called from Tinker: delete the Colvars proxy object, triggering post-run
void delete_colvars();


// FORTRAN routines exposed by Tinker to Colvars
// implemented in MOD_colvars.f

/// Get MPI communicator info from Tinker
void get_mpi_(int** commnuicator, int** rank, int** nrpoc);

void get_sim_temp_(cvm::real**);

void get_sim_boltzmann_(cvm::real**);

void get_sim_dt_(cvm::real**);

void rand_colvars_(cvm::real*);

/// Number of atoms
void get_system_size_(int* size);

/// Send IDs of Colvars atoms to Tinker
void set_cvatoms_ids_(int* ncvatoms, int * cvatoms_ids);

/// Get mass of given atom from Tinker
void get_mass_atom_(int* atom_number, double* mass);

/// Get position of given atom from Tinker
void get_pos_atom_(int* atom_number, double* xr, double* yr, double* zr);

/// Get total forces from Tinker
void get_forces_tot_(int* atom_number, double* fx, double* fy, double* fz);

/// Add Specified energy to Tinker
void add_energy_tinker_(double* energy);

/// Send given forces to Tinker on specified atoms
void add_forces_tinker_(int* index, cvm::rvector* force);

/// Apply minimum image convention using Tinker's PBC to a distance vector
void pbc_image_(cvm::rvector* r1, cvm::rvector* r2, cvm::rvector* dr);

/// Get the filename of the Tinker input
void get_input_filename_(char** input_name, int* len_name);

/// Check the existence of the colvars restart file through tinker
void get_restart_filename_(char** input_name, int* len_name);

/// Get PBC box dimensions, assuming orthorhombic cell
void get_pbc_(double* a, double *b, double *c);

/// Handle fatal error
void fatal_error_(void);

/// Handle colvars flag
void set_use_colvars_(bool* use_colvars);

// Lambda-dynamics interface

/// If lambda-dynamics enabled in back-end, special keyword in colvar config
/// Then the following is not necessary - Tinker doesn't need to know about the mass
/// Create special colvar tied to lambda with given mass and initial value
// void create_lambda_colvar_(double* mass, double* value);

/// Get value of alchemical lambda parameter from back-end
void get_alch_lambda_(double* lambda);

/// Set value of alchemical lambda parameter in back-end
void set_alch_lambda_(double* lambda);

/// Get (unbiased force field) energy derivative with respect to lambda
void get_delambda_(double* dE_dlambda);

/// Apply a scalar force on dE_dlambda (back-end distributes it onto atoms)
void apply_force_delambda_(double* force);

/// Get second derivative of the (unbiased force field) energy with respect to lambda
void get_d2edlambda2_(double* d2E_dlambda2);

#ifdef __cplusplus
}
#endif
