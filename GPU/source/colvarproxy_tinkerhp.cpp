#include "colvarproxy_tinkerhp.h"
#include "colvarproxy_tinkerhp_interface.h"
#include <colvarscript.h>
#include <mpi.h>

// This to avoid having to pass around proxy pointers
// This is unique to each MPI rank
static colvarproxy_tinkerhp* unique_colvarproxy_object;

// **********************************************
// The following three pure C functions are callable from Fortran
// technically they could live in a separate, tiny source file

void *allocate_colvars()
{
    colvarproxy_tinkerhp *proxy = new colvarproxy_tinkerhp();
    unique_colvarproxy_object = proxy;
    return proxy;
}

void compute_colvars ()
{
  unique_colvarproxy_object->compute();
}

void delete_colvars ()
{
  delete unique_colvarproxy_object;
}


// **********************************************
// Implementation of colvarproxy_tinkerhp follows

colvarproxy_tinkerhp::colvarproxy_tinkerhp()
{
  init();
}

colvarproxy_tinkerhp::~colvarproxy_tinkerhp()
{
  if ( rankcv != NULL && rankcv == 0) {
    // Finalize e.g. output from biases, and flush output streams
    post_run();
    // Delete the colvarmodule instance
    delete colvars;
  }
}

// Colvars Initialization
void colvarproxy_tinkerhp::init() {
  // Initialize colvars.
  get_system_size_(&n);
  get_sim_dt_(&sim_dt);
  cvm::real T;
  get_sim_temp_(&T);
  set_target_temperature(T);
  get_sim_boltzmann_(&sim_boltzmann);
  angstrom_value_ = 1.0;
  get_mpi_(&commcv,&rankcv,&nproccv);

  colvars_restart = false;
  // the input file must be "<input>.colvars"
  char* buff;
  int len_buff;
  get_input_filename_(&buff,&len_buff);
  if (len_buff>0) {
    bool use_colvars = true;
    set_use_colvars_(&use_colvars);
    config_file = std::string(buff,len_buff);
  // Output
    get_output_filename_(&buff,&len_buff);
    output_prefix_str = std::string(buff,len_buff);
  //
    get_restart_filename_(&buff,&len_buff);
    input_prefix_str = std::string(buff,len_buff);

//
  // initiate module: this object will be the communication proxy
    colvars = new colvarmodule (this);
    if (cvm::debug()) {
      cvm::log("Using TINKER-HP interface, version 1.0\n");
    }
    colvars->read_config_file(config_file.c_str());
    colvars->setup_input();
    colvars->setup_output();
  //Tinker-HP always starts at 1, so add 1 rightaway
    colvars->it++;

//    if (step != 0) {
//      cvm::log("Initializing step number to "+cvm::to_str(step)+".\n");
//      colvars->it = colvars->it_restart = step;
//    }
  
    if (cvm::debug()) {
      cvm::log ("colvars_atoms = "+cvm::to_str (atoms_ids)+"\n");
      cvm::log ("positions = "+cvm::to_str (atoms_positions)+"\n");
      cvm::log (cvm::line_marker);
    }

    int ncvatoms = atoms_ids.size();
    set_cvatoms_ids_(&ncvatoms, &atoms_ids[0]);
//  
    if (cvm::debug())
      log("done initializing the colvars proxy object.\n");
    {
    }
   // initialize multi-replica support, if available
  if (replica_enabled() == COLVARS_OK) {
    {
      int inter_comm_temp;
      get_root2root_(&inter_comm_temp);
      get_index_rep_(&inter_me);
      get_num_rep_(&inter_num);
      inter_comm = MPI_Comm_f2c( (MPI_Fint) inter_comm_temp );
      MPI_Comm_rank(inter_comm, &inter_me);
      MPI_Comm_size(inter_comm, &inter_num);
    }
  }
  }
  else {
    bool use_colvars = false;
    set_use_colvars_(&use_colvars);
  }
} // End colvars initialization.

// **** ATOMS ****
int colvarproxy_tinkerhp::init_atom(int atom_number)
{
  // save time by checking first whether this atom has been requested before
  // (this is more common than a non-valid atom number)
  int aid = (atom_number-1);

//  if ( (atom_number < 0) || (atom_number > n) ) {
//    cvm::fatal_error ("Error: invalid atom number specified, "+
//                      cvm::to_str (atom_number+1)+"\n");
//  }

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_ncopies[i] += 1;
      return i;
    }
  }

  if (aid < 0) {
    return COLVARS_INPUT_ERROR;
  }

  int const index = add_atom_slot(aid);
  // masses
  double mass;
  get_mass_atom_(&atom_number,&mass);
  atoms_masses.back() = mass;
  updated_masses_= true;
  return index;
}

void colvarproxy_tinkerhp::add_energy(cvm::real energy) {
  add_energy_tinker_(&energy);
  bias_energy += energy;
}

void colvarproxy_tinkerhp::request_total_force (bool yesno) {
  total_force_requested = yesno;
}

int colvarproxy_tinkerhp::check_atom_id(int atom_number)
{
  int const aid = atom_number-1;

  if (cvm::debug())
    log("Adding atom "+cvm::to_str(atom_number)+
        " for collective variables calculation.\n");

  // TODO add upper boundary check?
  if ( (aid < 0) ) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  return aid;
}

// **************** PERIODIC BOUNDARY CONDITIONS ****************
//  Get the PBC-aware distance vector between two positions
cvm::rvector colvarproxy_tinkerhp::position_distance (cvm::atom_pos const &pos1,
				cvm::atom_pos const &pos2) const {
  cvm::rvector r1, r2, dr;
  r1[0] = pos1.x;
  r1[1] = pos1.y;
  r1[2] = pos1.z;
  r2[0] = pos2.x;
  r2[1] = pos2.y;
  r2[2] = pos2.z;

  pbc_image_(&r1, &r2, &dr);
  return cvm::rvector( dr[0], dr[1], dr[2] );
}

int colvarproxy_tinkerhp::set_unit_system(std::string const &units_in, bool /*check_only*/)
{
  if (units_in != "real") {
    cvm::error("Error: Specified unit system \"" + units_in + "\" is unsupported in Tinker-HP. Supported units are \"real\" (A, kcal/mol).\n");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}

// trigger colvars computation
double colvarproxy_tinkerhp::compute()
{
// only the master does the calculation
   if (rankcv>0) {
     return 0.0;
   }

  double a,b,c;
  get_pbc_(&a,&b,&c);

  unit_cell_x.set(a, 0.0, 0.0);
  unit_cell_y.set(0.0, b, 0.0);
  unit_cell_z.set(0.0, 0.0, c);


  // only orthogonal unit cell in PBC with Tinker-HP
  boundaries_type = boundaries_pbc_ortho;
  colvarproxy_system::update_pbc_lattice();


  if (cvm::debug()) {
    cvm::log(std::string(cvm::line_marker)+
             "colvarproxy_tinker, step no. "+cvm::to_str(colvars->it)+"\n"+
             "Updating internal data.\n");
  }

  // zero the forces on the atoms, so that they can be accumulated by the colvars
  for (size_t i = 0; i < atoms_new_colvar_forces.size(); i++) {
    atoms_new_colvar_forces[i].reset();
  }

  // get the new positions
  for (size_t i = 0; i < atoms_ids.size(); i++) {
    cvm::rvector r1;
    int atom_number = i+1;
    get_pos_atom_(&atom_number,&r1[0],&r1[1],&r1[2]);
    atoms_positions[i] = r1;
  }

  bias_energy = 0.0;

  // fill the total forces array
  for (size_t i = 0;i <atoms_ids.size(); i++) {
    int atom_number = i+1;
    cvm::rvector forces;
    get_forces_tot_(&atom_number,&forces[0],&forces[1],&forces[2]);
    atoms_total_forces[i] = forces;
  }

  if (cvm::debug()) {
    cvm::log("atoms_ids = "+cvm::to_str(atoms_ids)+"\n");
    cvm::log("atoms_ncopies = "+cvm::to_str(atoms_ncopies)+"\n");
    cvm::log("atoms_positions = "+cvm::to_str(atoms_positions)+"\n");
    cvm::log("atoms_new_colvar_forces = "+cvm::to_str(atoms_new_colvar_forces)+"\n");
    cvm::log("atoms_total_forces = "+cvm::to_str(atoms_total_forces)+"\n");
    cvm::log("CVM timestep: " + cvm::to_str(cvm::step_absolute())+"\n");
    cvm::log("CVM temperature:  " + cvm::to_str(target_temperature()) + "\n");
  }

  // call the collective variable module
  colvars->calc();

  // Colvars maintains its own time step count (Tinker-HP always starts at 1)
  colvars->it++;

  if (cvm::debug()) {
    std::cout<<"dt "<<dt()<<std::endl;
    std::cout<<"bias "<<bias_energy<<std::endl;
  }

  // fill the forces array
  for (size_t i = 0;i <atoms_ids.size(); i++) {
    int atom_number = i+1;
    cvm::rvector forces;
    forces = atoms_new_colvar_forces[i];
    add_forces_tinker_(&atom_number,&forces);
  }


  if (cvm::debug()) {
    cvm::log("Value of first colvar: " + cvm::to_str((*colvars->variables())[0]->value()) + "\n");
    cvm::log("atoms_ids = "+cvm::to_str(atoms_ids)+"\n");
    cvm::log("atoms_ncopies = "+cvm::to_str(atoms_ncopies)+"\n");
    cvm::log("atoms_positions = "+cvm::to_str(atoms_positions)+"\n");
    cvm::log("atoms_new_colvar_forces = "+cvm::to_str(atoms_new_colvar_forces)+"\n");
    cvm::log("atoms_total_forces = "+cvm::to_str(atoms_total_forces)+"\n");
  }

  return bias_energy;
}

// multi-replica support

int colvarproxy_tinkerhp::replica_enabled()
{
  bool doreplica = false;
  get_use_reps_(&doreplica);
  return (doreplica) ? COLVARS_OK : COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_tinkerhp::replica_index()
{
  int index;
  get_index_rep_(&index);
  return index;
}

//
int colvarproxy_tinkerhp::num_replicas()
{
  int inter_num;
  get_num_rep_(&inter_num);
  return inter_num;
}
//
//
void colvarproxy_tinkerhp::replica_comm_barrier()
{
  MPI_Barrier(inter_comm);
}
//
//
int colvarproxy_tinkerhp::replica_comm_recv(char* msg_data,
                                          int buf_len, int src_rep)
{
  MPI_Status status;
  int retval;

  retval = MPI_Recv(msg_data,buf_len,MPI_CHAR,src_rep,0,inter_comm,&status);
  if (retval == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_CHAR, &retval);
  } else retval = 0;
  return retval;
}
//
//
int colvarproxy_tinkerhp::replica_comm_send(char* msg_data,
                                          int msg_len, int dest_rep)
{
  int retval;
  retval = MPI_Send(msg_data,msg_len,MPI_CHAR,dest_rep,0,inter_comm);
  if (retval == MPI_SUCCESS) {
    retval = msg_len;
  } else retval = 0;
  return retval;
}
/// Use embedded Tcl interpreter for callbacks if available
#ifdef COLVARS_TCL

int colvarproxy_tinkerhp::run_force_callback()
{
  int res = colvarproxy::tcl_run_force_callback();
  if (res == COLVARS_NOT_IMPLEMENTED)
    cvm::log("Tcl scripts are not enabled in this build of the Colvars library / Tinker-HP.\n");
  if (res != COLVARS_OK)
    cvm::error("Error running Tcl forces script.\n");
  return res;
}

int colvarproxy_tinkerhp::run_colvar_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         colvarvalue &value)
{
  int res = colvarproxy::tcl_run_colvar_callback(name, cvc_values, value);
  if (res == COLVARS_NOT_IMPLEMENTED)
    cvm::log("Tcl scripts are not enabled in this build of the Colvars library / Tinker-HP.\n");
  if (res != COLVARS_OK)
    cvm::error("Error running Tcl colvar script.\n");
  return res;
}

int colvarproxy_tinkerhp::run_colvar_gradient_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         std::vector<cvm::matrix2d<cvm::real> > &gradient)
{
  int res = colvarproxy::tcl_run_colvar_gradient_callback(name, cvc_values,
                                                       gradient);
  if (res == COLVARS_NOT_IMPLEMENTED)
    cvm::log("Tcl scripts are not enabled in this build of the Colvars library / Tinker-HP.\n");
  if (res != COLVARS_OK)
    cvm::error("Error running Tcl colvar gradient script.\n");
  return res;
}

#endif
