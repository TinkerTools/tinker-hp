//          init_ml_ressources(   rank, devID,  nn_name, model_file, debug ):
extern int32_t init_ml_ressources(int32_t, int32_t, char*, char*, int32_t);

//          ml_models (   coord,   energy,  gradient,  cell,  atomic_species, neighl1, neighl2, dist,     dxyz,   istrict,    natm,  npairs,   nstrict, dograd)
extern int32_t ml_models (uint64_t, uint64_t, uint64_t, float *, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t,uint64_t,  int32_t, int32_t, int32_t, int32_t);

//   ml_disp ( array_ptr, size )
extern int32_t ml_debug(uint64_t, int32_t);

extern void nuke_context( int );
