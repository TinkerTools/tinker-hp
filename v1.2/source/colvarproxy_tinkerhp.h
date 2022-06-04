/// Implementation of the colvarproxy class for Tinker-HP
/// Interfaces to Fortran via C functions declared in colvarproxy_tinkerhp_interface.h


#include <colvarproxy.h>
#include "colvarproxy_tinkerhp_interface.h"

/// \brief Communication between colvars and TINKER
/// (implementation of \link colvarproxy \endlink)
class colvarproxy_tinkerhp : public colvarproxy {
public:
  colvarproxy_tinkerhp();
  ~colvarproxy_tinkerhp();

  // perform colvars computation. returns biasing energy
  double compute();
  /// Print a message to the main log
  void log(std::string const &message) { std::cout << "colvars: " << message; }

  // methods for tinker to move data or trigger actions in the proxy
  bool total_forces_enabled() const { return total_force_requested; };

private:
  std::string config_file;

  bool total_force_requested;
  bool colvars_restart;
  double bias_energy;

  /// initialization
  void init();
  /// size of the system
  int n;
  /// \brief Value of the unit for atomic coordinates with respect to
  /// angstroms (used by some variables for hard-coded default values)
  cvm::real backend_angstrom_value()
  {
    return 1.0;
  }
  int set_unit_system(std::string const &units_in, bool check_only);

  /// \brief target Boltzmann constant
  cvm::real* sim_boltzmann;
  cvm::real boltzmann() { return *sim_boltzmann; }

  /// \brief Target temperature of the simulation (K units)
  cvm::real* sim_temperature;

  cvm::real temperature() {return *sim_temperature;}

  /// \brief Target Time step of the simulation (fs)
  cvm::real* sim_dt; 

  /// \brief Time step of the simulation (fs)
  cvm::real dt() { return *sim_dt;}

  /// \brief target mpi communicator, associated rank and number of procs
  int* commcv;
  int* rankcv;
  int* nproccv;

  /// \brief Pseudo-random number with Gaussian distribution
  cvm::real rand_gaussian(){
    double temp_rand;
    rand_colvars_(&temp_rand);
    return temp_rand;
    } 
  /// Pass restraint energy value for current timestep to MD engine
  void add_energy(cvm::real energy);

  void request_total_force (bool yesno);

  /// numeric index (1-based)
  int init_atom(int atom_number);

  /// Check that this atom number is valid, but do not initialize the
  /// corresponding atom yet
  int check_atom_id(int atom_number);

  // **************** PERIODIC BOUNDARY CONDITIONS ****************
  cvm::rvector position_distance (cvm::atom_pos const &pos1,
				                          cvm::atom_pos const &pos2) const;

  /// Print a message to the main log and let the rest of the program handle the error
  /// In Tinker-HP we consider all errors as fatal
  inline void error(std::string const &message) { std::cout << "colvars: " << message << std::endl; fatal_error_(); }
 
  /// Print a message to the main log and exit with error code
  inline void fatal_error(std::string const &message) { std::cout << "colvars: " <<  message << std::endl; fatal_error_(); }

  /// Get value of alchemical lambda parameter from back-end
  int get_alch_lambda(cvm::real* lambda) {
    // Call C/Fortran implementation
    get_alch_lambda_(lambda);
    return COLVARS_OK;
  }

  /// Set value of alchemical lambda parameter in back-end
  int send_alch_lambda() {
    // Call C/Fortran implementation
    set_alch_lambda_(&cached_alch_lambda);
    return COLVARS_OK;
  }

  /// Get energy derivative with respect to lambda
  int get_dE_dlambda(cvm::real* dE_dlambda) {
    // Call C/Fortran implementation
    get_delambda_(dE_dlambda);
    return COLVARS_OK;
  }

  /// Get energy derivative with respect to lambda
  int apply_force_dE_dlambda(cvm::real* force) {
    // Call C/Fortran implementation
    apply_force_delambda_(force);
    return COLVARS_OK;
  }

  /// Get energy derivative with respect to lambda
  int get_d2E_dlambda2(cvm::real* d2E_dlambda2) {
    // Call C/Fortran implementation
    get_d2edlambda2_(d2E_dlambda2);
    return COLVARS_OK;
  }
};
