global_load=0
use_custom_torchani=False
try:
  import traceback
  from typing import Optional
  from time import time
  import sys
  import os
  import numpy as np
  try:
    from mlplugin import ffi
  except:
    ffi = lambda : None
    ffi.def_extern = lambda : lambda a:a
    print("Warning: mlplugin not found, ffi was defined as dummy object")
    global_load=1

  try:
    print(sys.argv[0])
  except:
    sys.argv.append('mlinterface.py')
  
  _hartree2kcalmol = 627.5094738
  _ev2kcalmol      =  23.06054
  # _default_model_dir is modified by mlbuilder.py to point to the Tinker intallation directory
  _default_model_dir = "."
  try:
    _model_dir = os.environ["TINKER_ML_DIR"]
  except:
    #print("TINKER_ML_DIR environment variable not set."
    #  ," Using default model directory: "+_default_model_dir, flush=True)
    _model_dir = _default_model_dir
  _ctype2tensordtype = {}
  _ctype2dtype       = {}
  _ctype2tfdtype     = {}
  _tf_models = ["DEEPMD"]
  _torchani_models = ["ANI_GENERIC","ANI1X","ANI1CCX","ANI2X","ML_MBD"]

  def load_modules(ml_key:str,rank:int=0,debug:bool=False)->int:
    """ Load modules according to the model to initialize """
    key_upper=ml_key.upper()
    try:
      if key_upper in _tf_models:
        if debug:
          os.environ["TF_CPP_MIN_LOG_LEVEL"] = "0"
        else:
          os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
        #os.environ["KMP_AFFINITY"] = "granularity=fine,verbose,compact,1,0"
        #os.environ["OMP_NUM_THREADS"] = "1"
        #os.environ["KMP_BLOCKTIME"] = "1"
        #os.environ["TF_INTRA_OP_PARALLELISM_THREADS"] = "1"
        #os.environ["TF_INTER_OP_PARALLELISM_THREADS"] = "1"
        global tf
        if key_upper == "DEEPMD":
          if (rank==0): print("loading DEEPMD module",flush=True)
          global DeepPot,default_tf_session_config
          from deepmd.env import tf,default_tf_session_config
          tf.compat.v1.enable_eager_execution()
          from deepmd.infer import DeepPot
          from  deepmd.loggers import set_log_handles
          if debug:
            set_log_handles(10)
          else:
            set_log_handles(40)
        else:
          import tensorflow as tf
          tf.compat.v1.enable_eager_execution()
      elif key_upper in _torchani_models:
        global torch, use_custom_torchani
        import torch
        if use_custom_torchani:
          try:
            global models, MOD_neigh
            from torchanimulti_2x import models
            import MOD_neighborlist as MOD_neigh
            if (rank==0): print("loading torchanimulti_2x module",flush=True)
          except:
            if (rank==0): print("torchanimulti_2x module not found",flush=True)
            use_custom_torchani = False
        if not use_custom_torchani:
          if (rank==0): print("Loading standard torchani",flush=True)
          global models, torchani, NbList,MOD_neigh
          from torchani  import models
          from torchani.aev import NbList
          MOD_neigh = None
      try:
        global pcd, GPUArray
        import pycuda.driver   as pcd
        from pycuda.gpuarray import GPUArray
      except Exception as exp:
        print('Exception: Fail to load pycuda with exception:', exp, flush=True)
        return 2
    except Exception as exp:
      print('Exception: Fail to load modules with exception:', exp, flush=True)
      return 1
    return 0

  def build_ctypes_converters()->int:
    """ Build dictionary of type correspondance with ctypes 
       (for numpy, pytorch and tensorflow) """
    try:
      # Build numpy.dtype dict 
      for prefix in ('int', 'uint'):
          for log_bytes in range(4):
            ctype  = '%s%d_t' % (prefix, 8 * (2**log_bytes))
            dtype  = '%s%d' % (prefix[0], 2**log_bytes)
            _ctype2dtype[ctype] = np.dtype(dtype)
      _ctype2dtype['float'] = np.dtype('f4')
      _ctype2dtype['double'] = np.dtype('f8')

      # Build tensor.dtype dict
      if "torch" in globals():
        _ctype2tensordtype["int8_t"] = torch.int8
        _ctype2tensordtype["short"]  = torch.short
        _ctype2tensordtype["int16_t"]= torch.int16
        _ctype2tensordtype["int"]    = torch.int
        _ctype2tensordtype["int32_t"]= torch.int32
        _ctype2tensordtype["int64_t"]= torch.int64
        _ctype2tensordtype["long"]   = torch.long
        _ctype2tensordtype["float"]  = torch.float
        _ctype2tensordtype["float32"]= torch.float32
        _ctype2tensordtype["double"] = torch.double
        _ctype2tensordtype["float64"]= torch.float64

      # Build tf.dtype dict
      if "tf" in globals():
        _ctype2tfdtype["int8_t"] = tf.int8
        _ctype2tfdtype["int16_t"]= tf.int16
        _ctype2tfdtype["int32_t"]= tf.int32
        _ctype2tfdtype["int64_t"]= tf.int64
        _ctype2tfdtype["float16"]= tf.float16
        _ctype2tfdtype["float"]  = tf.float32
        _ctype2tfdtype["float32"]= tf.float32
        _ctype2tfdtype["double"] = tf.double
        _ctype2tfdtype["float64"]= tf.float64
      
      return 0
    except Exception as err:
      print('Exception: Fail to build ctypes with exception:', err, flush=True)
      #print(traceback.format_exc(),flush=True)
      return 1

  class AniRessources:
    nb_species   = 0
    mod_torch    = 0
    mod_tf       = 1

    def __init__(self,rank:int ,devID:int, mlpot_key:str
                    ,model_file:Optional[str]=None
                    ,verbose:bool=False
                    ,debug:bool=False)->None:

      self.verbose = debug and mlpot_key in _torchani_models
      if verbose:
        init_time = time()

      self.mlpot_key  = mlpot_key
      self.rank       = rank
      self.devID      = devID
      self.use_mod    = self.mod_torch if mlpot_key in _torchani_models else self.mod_tf

      try:
        # SELECT DEVICE
        if mlpot_key in _torchani_models:
            #pcd.init()
            torch.cuda.init()
            #self.cuDeviceId = pcd.Device(devID)
            #self.context    = self.cuDeviceId.retain_primary_context()
            device_name     = 'cuda:'+str(devID) if torch.cuda.is_available() else 'cpu'
            self.device     = torch.device(device_name)
            if self.verbose: print("I rank ",rank," select device "+device_name,flush=True)
            try:
               test = torch.tensor(np.empty(512),dtype=torch.float32,device=self.device)
            except Exception as exp:
               raise Exception('test torch '+str(exp))
        elif mlpot_key in _tf_models:
            default_tf_session_config.device_count['GPU']=1
            default_tf_session_config.gpu_options.visible_device_list=str(devID)
            physical_devices = tf.config.list_physical_devices('GPU')
            # Only use the GPU corresponding to devID
            tf.config.set_visible_devices(physical_devices[devID], 'GPU')
            logical_devices = tf.config.list_logical_devices('GPU')
            assert len(logical_devices) == 1
            self.device = tf.device('/GPU:0')
        
        # SELECT MODEL
        if mlpot_key in _torchani_models:
            if (self.verbose): print( f'Loading {mlpot_key} model...', flush=True)
            if use_custom_torchani: MOD_neigh.verbose = self.verbose
            self.pbc       = torch.tensor([True, True, True]).to(self.device)
            if self.mlpot_key == "ANI1CCX":
              self.model = models.ANI1ccx(periodic_table_index=True).to(self.device)
            elif self.mlpot_key == "ANI1X":
              self.model   = models.ANI1x(periodic_table_index=True).to(self.device)
            elif self.mlpot_key == "ANI2X":
              self.model   = models.ANI2x(periodic_table_index=True).to(self.device)
            elif self.mlpot_key == "ML_MBD":
              self.model = torch.jit.load(_model_dir+'/mlmbd_model.pt').to(self.device)
            elif self.mlpot_key=="ANI_GENERIC" and model_file is not None:
              self.model_file=model_file
              print(f'(from file "{model_file}")',flush=True)
              if model_file.endswith(".pt"):
                self.model = torch.jit.load(model_file).to(self.device)
              elif model_file.endswith(".json"):
                self.model = models.BuiltinModel.from_json(
                                model_file, periodic_table_index=True
                              ).to(self.device)
              elif model_file.endswith(".yaml"):
                self.model = models.BuiltinModel.from_yaml(
                                model_file, periodic_table_index=True
                              ).to(self.device)
              elif model_file.endswith(".pkl"):
                self.model = models.BuiltinModel.from_pickle(
                                model_file, periodic_table_index=True
                              ).to(self.device)
              else:
                raise ValueError("model file "+model_file+" not supported",flush=True)
            else:
              raise ValueError("ani model "+self.mlpot_key+" not supported",flush=True)
            self.model.train(False)
        elif mlpot_key in _tf_models:
            if (self.verbose): print(f'Loading {mlpot_key} model...', flush=True)
            if self.mlpot_key=="DEEPMD" and model_file is not None:
              self.model_file=model_file
              print(f'(from file "{model_file}")',flush=True)
              with self.device:
                self.model = DeepPot(model_file)
            else:
              print('Warning: Loading default DeePMD model...', flush=True)
              with self.device:
                self.model = DeepPot(_model_dir+'/deepmd_scanprl2021_watermodel.pb')
        else:
            raise ValueError("model "+mlpot_key+" not found", flush=True)
        
        

      except pcd.CompileError as cex:
        raise Exception('CompileError '+str(cex))
      except Exception as exp:
        raise Exception('arRessources init ',str(exp))
      except:
        raise Exception('Unknown Exception detected in AniRessources.__init__')
    
    def compute(self,coord_ptr, atm_ener_ptr, gradient_ptr, cell_ptr, atm_sp_ptr
        , neigh1Idx_ptr, neigh2Idx_ptr, dist_ptr, dxyz_ptr, istrict_ptr, nb_species, nb_pairs, nb_strict, dograd:bool)->int:

      ierr=0
      if self.verbose: init_time = time()

      if self.mlpot_key in _torchani_models:
        torchani_mbd_code = self.mlpot_key == 'ML_MBD'
        if not torchani_mbd_code:
          self.getmcache(nb_species)

        cell             = torch.from_numpy(asarray(ffi, cell_ptr, [3,3])).to(self.device)
        atomic_species   = asTensor(   atm_sp_ptr,'int64_t',[1,nb_species],self.device,order='F')
        coordinates      = asTensor(    coord_ptr,  'float',[1,nb_species,3],self.device,order='F').requires_grad_()
        atomic_energies  = asTensor( atm_ener_ptr,  'float',[1,nb_species],self.device)
        if dograd: gradient  = asTensor( gradient_ptr, 'float',[nb_species,3],self.device)

        if use_custom_torchani:
          MOD_neigh.atom_index1_temp = asTensor(neigh1Idx_ptr,'int32_t',[1,nb_pairs],self.device).type(torch.long)
          MOD_neigh.atom_index2_temp = asTensor(neigh2Idx_ptr,'int32_t',[1,nb_pairs],self.device).type(torch.long)
          MOD_neigh.distance_temp    = asTensor(dist_ptr ,  'float',[1,nb_pairs],self.device)
          with torch.no_grad():
            atom_index12 = torch.stack([MOD_neigh.atom_index1_temp.squeeze(),MOD_neigh.atom_index2_temp.squeeze()])
            selected_coordinates = coordinates.flatten(0, 1).index_select(0, atom_index12.view(-1)).view(2, -1, 3)
            shifts=asTensor(dxyz_ptr ,  'float',[nb_pairs,3],self.device,order='F')
            MOD_neigh.vec_temp  = selected_coordinates[0] - selected_coordinates[1] + shifts
        else:
          atom_index_1=asTensor(neigh1Idx_ptr,'int32_t',[1,nb_pairs],self.device).type(torch.long)
          atom_index_2 = asTensor(neigh2Idx_ptr,'int32_t',[1,nb_pairs],self.device).type(torch.long)
          nblist = NbList(
            atom_index12 = torch.stack([atom_index_1.squeeze(),atom_index_2.squeeze()])
            ,shifts = asTensor(dxyz_ptr ,  'float',[nb_pairs,3],self.device,order='F')
          ) 
        istrict = asTensor(istrict_ptr,'int32_t',[nb_species],self.device).type(torch.long)[:nb_strict]
        
        if self.verbose:
          torch.cuda.synchronize()
          print(' ML input GPU processes', flush=True)
          print(torch.cuda.list_gpu_processes(), flush=True)
          print(' calc ML input {:.6f} s'.format(time()-init_time), flush=True)

        if not torchani_mbd_code:
          if self.verbose:
            init_time_model = time()

          if use_custom_torchani:
            predictor = self.model((atomic_species, coordinates), cell=cell, pbc=self.pbc)
          else:
            predictor = self.model.atomic_energies((atomic_species, coordinates)
              , cell=cell, pbc=self.pbc ,nblist=nblist
              ,shift_energies=False)

          if self.verbose:
            torch.cuda.synchronize()
            print(' calc ML model {:.6f} s'.format(time()-init_time_model), flush=True)
            init_time_transf = time()

          if use_custom_torchani:
            atomic_energies[:] = _hartree2kcalmol*MOD_neigh.aeani_temp[:]
          else:
            atomic_energies[:] = _hartree2kcalmol*predictor.energies

          if self.verbose:
            torch.cuda.synchronize()
            init_time_grad = time()
            print(' calc ML transf {:.6f} s'.format(time()-init_time_transf), flush=True)

          if dograd:
            gradient[:] = _hartree2kcalmol*torch.autograd.grad(
                              predictor.energies[0][istrict].sum()
                              ,coordinates
                           )[0].squeeze().squeeze()

            if self.verbose:
              torch.cuda.synchronize()
              print(' calc ML grad {:.6f} s'.format(time()-init_time_grad), flush=True)

        else:
          if self.verbose:
            init_time_model = time()

          dict_species = {1: 0, 6: 1, 7: 2, 8: 3}
          atomic_species[0][:] = torch.tensor([dict_species[x[:]] for x in atomic_species[0]], device=self.device)

          if use_custom_torchani:
            predictor = self.model((atomic_species, coordinates))
          else:
            predictor = self.model.atomic_energies((atomic_species, coordinates) ,nblist=nblist)
            
          if self.verbose:
            torch.cuda.synchronize()
            print(' calc ML-MBD model {:.6f} s'.format(time()-init_time_model), flush=True)

          atomic_energies[:] = predictor.energies[:]

          if dograd:
            if self.verbose: init_time_grad = time()
            gradient[:] = torch.autograd.grad(
                            predictor.energies[0][istrict].sum()
                            ,coordinates
                          )[0].squeeze().squeeze()

            if self.verbose:
              torch.cuda.synchronize()
              print(' calc ML-MBD grad {:.6f} s'.format(time()-init_time_grad), flush=True)

      elif self.mlpot_key == 'DEEPMD':
        atomic_energies    = asGPUarray(atm_ener_ptr, 'float',[1,nb_species])
        if dograd: gradient= asGPUarray(gradient_ptr, 'float',[nb_species,3])

        with self.device:
          cell             = tf.convert_to_tensor(asarray(ffi, cell_ptr, [1,9]))
          atomic_species   = astfTensor(atm_sp_ptr,'int32_t',[nb_species], order='F').numpy()
          coordinates      = astfTensor(coord_ptr,  'float',[1,nb_species,3],order='F').numpy().reshape([1,-1])

          ########################################################
          ### can be useful if dict_species become complicated ###
          #dict_species     = {1: 1, 8: 0}   
          #for old,new in dict_species.items():
          #  if old!=new: atomic_species[atomic_species==old] = new
          ########################################################
          
          # much simpler here
          atomic_species[atomic_species==8] = 0

          e, f, v, ae, av = self.model.eval(coordinates, cell, atomic_species, True)
          atomic_energies[0,:] = _ev2kcalmol*ae.flatten().astype('float32')
          if dograd: gradient[:] = -_ev2kcalmol*f.reshape([-1,3]).astype('float32')

      else:
        print('''Codes other than Pytorch and TensorFlow are not yet available.
                 Please modify: mlbuilder.py.''', flush=True)
        ierr = 1

      if self.verbose:
        torch.cuda.synchronize()
        print(' calc ML total {:.6f} s'.format(time()-init_time), flush=True)
      
      return ierr

    def getmcache(self, ns ):
      """ create a big Tensor to trigger pytorch memory allocation"""
      if ns > self.nb_species:
        self.nb_species = int(ns * 1.10)
        allocSize  = np.int64(self.nb_species)*150*(8+4+32+20)*3
        allocSize += np.int64(self.nb_species)*384*3

        if self.verbose:
            torch.cuda.synchronize()
            print('getmcache',torch.cuda.list_gpu_processes(), flush=True)

        if (self.use_mod==self.mod_torch):
           Ten  = torch.empty([allocSize],dtype=torch.float32,device=self.device,requires_grad=False)

        if self.verbose:
            torch.cuda.synchronize()
            print('getmcache reserved', (allocSize*4)/1e6, 'MB on '
              , self.device, torch.cuda.list_gpu_processes(), flush=True)
      return 0

  def asarray(ffi, ptr, shape, **kwargs):
    length = np.prod(shape)
    T = ffi.getctype(ffi.typeof(ptr).item)
    if T not in _ctype2dtype:
      raise RuntimeError("Cannot create an array for element type: %s" % T)
    a = np.frombuffer(ffi.buffer(ptr, length * ffi.sizeof(T)), _ctype2dtype[T]).reshape(shape, **kwargs, order='F')
    return a

  def asGPUarray(ptr, ctype, shape, **kwargs):
    return GPUArray(gpudata=ptr,shape=shape,dtype=_ctype2dtype[ctype], **kwargs)

  def asTensor(ptr, ctype, shape, device, **kwargs):
    ga     = asGPUarray(ptr, ctype, shape, **kwargs)
    return torch.as_tensor(ga,dtype=_ctype2tensordtype[ctype],device=device)

  def astfTensor(ptr, ctype, shape, **kwargs):
    ga     = asGPUarray(ptr, ctype, shape, **kwargs)
    return tf.convert_to_tensor(ga.get())

except Exception as err:
  global_load = 1
  try: 
    import traceback
    print(traceback.format_exc(),flush=True)
  except:
    pass
  print('Failed global initialization with exception:', err, flush=True)

#-------------------------------------------------------------------------------
# define actual interface

@ffi.def_extern()
def init_ml_ressources(rank,devID,nn_name,model_file_,debug_int):
  try:
     if global_load != 0: return 1
     init_time  = time()
     ierr       = 0
     debug      = debug_int != 0
     ml_key     = ffi.string(nn_name).decode('UTF-8').strip().upper()
     model_file = ffi.string(model_file_).decode('UTF-8').strip()
     if (debug and rank==0): print('init ML ressources',rank,model_file,debug_int,ml_key,flush=True)

     load_err = load_modules(ml_key,rank,debug)
     if  load_err != 0: return 10+load_err
    
     ctype_err = build_ctypes_converters()
     if ctype_err != 0: return 20+ctype_err

     if debug: print(f'_model_dir="{_model_dir}"')

     global ar
     ar=AniRessources(rank,devID,ml_key,model_file=model_file,debug=debug)
     if(ar.use_mod==ar.mod_torch): torch.cuda.synchronize()

     if (ar.verbose): print(f'init ML ressources done {time()-init_time} s'.format(), flush=True)
     if rank==0:
        print('', flush=True)
        print(' *****    Using ML potential engine    ***** ', flush=True)
        print('', flush=True)
  except Exception as err:
     ierr = 1
     if (rank<4):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",flush=True)
        print(traceback.format_exc(),flush=True)
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<",flush=True)
     #raise Exception('init_ml_ressources:'+str(err))
  return ierr

@ffi.def_extern()
def ml_models(coord_ptr, atm_ener_ptr, gradient_ptr, cell_ptr, atm_sp_ptr
   ,neigh1Idx_ptr, neigh2Idx_ptr, dist_ptr, dxyz_ptr, istrict_ptr, nb_species, nb_pairs, nb_strict, dograd_int):
  try:
    dograd = dograd_int != 0
    ierr = ar.compute(coord_ptr, atm_ener_ptr, gradient_ptr, cell_ptr, atm_sp_ptr
         , neigh1Idx_ptr, neigh2Idx_ptr, dist_ptr, dxyz_ptr,istrict_ptr, nb_species, nb_pairs, nb_strict, dograd)
  except Exception as err:
    ierr = 1
    if (ar.rank<4):
       print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",flush=True)
       print(traceback.format_exc(),flush=True)
       print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<",flush=True)
    #raise Exception("ml_models "+str(err))
  return ierr

@ffi.def_extern()
def ml_debug(array_ptr, size):
  try:
    print('-----------mdi',ar.rank,size,array_ptr,type(array_ptr),flush=True)
    array = asTensor(array_ptr,'float',[1,size],ar.device,order='F')
    print('-----------mdi',array[0,0:9],flush=True)
    return 0
  except Exception as err:
    print(traceback.format_exc(),flush=True)
    return 1



@ffi.def_extern()
def nuke_context( i ):
  pass
  #try:
  #   if ar.mlpot_key in _torchani_models:
  #      if ar.verbose: print('ML potential pop CUDA Context', flush=True)
  #      #ar.context.pop()
  #except Exception as err:
  #   raise Exception('nuke_context'+str(err))
