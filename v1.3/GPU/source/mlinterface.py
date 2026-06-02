global_load = 0
try:
    import traceback
    from time import time, sleep
    import sys
    import os
    import numpy as np
    import gc
    import yaml

    try:
        from mlplugin import ffi
    except:
        ffi = lambda: None
        ffi.def_extern = lambda: lambda a: a
        print("Warning: mlplugin not found, ffi was defined as dummy object")
        global_load = 1

    try:
        print(sys.argv[0])
    except:
        sys.argv.append("mlinterface.py")

    # _default_model_dir is modified by mlbuilder.py to point to the Tinker intallation directory
    # _default_model_dir = "."
    # try:
    #     _model_dir = os.environ["TINKER_ML_DIR"]
    # except:
    #     # print("TINKER_ML_DIR environment variable not set."
    #     #  ," Using default model directory: "+_default_model_dir, flush=True)
    #     _model_dir = _default_model_dir
    _ctype2dtype = {}

    class GPUPointer:
        def __init__(self, ptr: int, size: int, dtype: np.dtype, device: int = 0):
            self.ptr = ptr
            self.size = size
            self.dtype = np.dtype(dtype)
            self.device = device
            self.__cuda_array_interface__ = {
                "shape": (self.size,),
                "typestr": self.dtype.str,
                "data": (int(self.ptr), False),
                "stream": None,
                "version": 3,
            }

    def load_modules(
        rank: int = 0,
        devId: int = 0,
        port: int = 45321,
    ) -> int:
        """Load modules according to the model to initialize"""

        #### load CuPy
        try:
            global cp, UnownedMemory,MemoryPointer
            import cupy as cp
            from cupy.cuda import UnownedMemory,MemoryPointer

        except Exception as exp:
            print("Exception: Failed to load cupy with exception:", exp, flush=True)
            return 2

        #### load jax
        try:
            global jax, jnp
            import jax
            import jax.numpy as jnp

            if port >= 0:
                if port == 0:
                    import socket
                    from contextlib import closing

                    sleep(np.abs(np.random.normal() + rank))
                    with closing(
                        socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                    ) as s:
                        s.bind(("localhost", 0))
                        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                        port = s.getsockname()[1]
                else:
                    port = port + rank

                jax.distributed.initialize(
                    coordinator_address=f"localhost:{port}",
                    num_processes=1,
                    process_id=0,
                    local_device_ids=devId,
                )

                device = jax.devices("gpu")[0]
                port_str = f", port {port}"

            else:
                device = jax.devices("gpu")[devId]
                port_str = ""

            print(f"rank {rank}, jax device", device, port_str, flush=True)
            jax.config.update("jax_default_device", device)
            jax.config.update("jax_default_matmul_precision", "highest")
        except Exception as exp:
            print("Exception: Failed to load jax with exception:", exp, flush=True)
            return 1

        ### load fennol
        try:
            global fennol, au, PERIODIC_TABLE
            import fennol
            from fennol.utils import AtomicUnits as au
            from fennol.utils.periodic_table import PERIODIC_TABLE
        except Exception as exp:
            print("Exception: Failed to load fennol with exception:", exp, flush=True)
            return 3

        return 0

    def build_ctypes_converters() -> int:
        """Build dictionary of numpy type correspondance with ctypes"""
        try:
            _ctype2dtype["int32_t"] = np.int32
            _ctype2dtype["int64_t"] = np.int64
            _ctype2dtype["uint32_t"] = np.uint32
            _ctype2dtype["uint64_t"] = np.uint64
            _ctype2dtype["float"] = np.float32
            _ctype2dtype["double"] = np.float64

            return 0
        except Exception as err:
            print("Exception: Failed to build ctypes with exception:", err, flush=True)
            # print(traceback.format_exc(),flush=True)
            return 1

    def asarray(ptr, shape, **kwargs):
        length = np.prod(shape)
        T = ffi.getctype(ffi.typeof(ptr).item)
        if T not in _ctype2dtype:
            raise RuntimeError("Cannot create an array for element type: %s" % T)
        a = np.frombuffer(
            ffi.buffer(ptr, length * ffi.sizeof(T)), _ctype2dtype[T]
        ).reshape(
            shape, **kwargs
        )  # , order="F")
        return a

    def writeGPUptr(ptr, array, wait=True):
        if wait:
            array.block_until_ready()
        shape = array.shape if array.ndim > 0 else (1,)
        mem = UnownedMemory(ptr, array.nbytes, owner=None)
        memptr = MemoryPointer(mem, offset=0)
        out = cp.ndarray(shape, dtype=array.dtype, memptr=memptr)

        mem = UnownedMemory(array.unsafe_buffer_pointer(), array.nbytes, owner=array)
        memptr = MemoryPointer(mem, offset=0)
        out[:] = cp.ndarray(shape, dtype=array.dtype, memptr=memptr)[:]
        # out = GPUArray(gpudata=ptr, shape=array.shape, dtype=array.dtype)
        # out.set(
        #     GPUArray(
        #         gpudata=array.unsafe_buffer_pointer(),
        #         shape=array.shape,
        #         dtype=array.dtype,
        #     )
        # )

    def asJaxArray(ptr, ctype, shape):
        length = np.prod(shape)
        T = GPUPointer(ptr, length, _ctype2dtype[ctype])
        return jnp.asarray(T, copy=False).reshape(shape)

except Exception as err:
    global_load = 1
    try:
        import traceback

        print(traceback.format_exc(), flush=True)
    except:
        pass
    print("Failed global initialization with exception:", err, flush=True)

# -------------------------------------------------------------------------------
# define actual interface


@ffi.def_extern()
def init_ml_ressources(rank_, devID, model_file_, debug_int, port, use_lambda_int, qtot_, qligand_, gc_stride_,
    nat, species_ptr, nligand, index_ligand_ptr):
    ierr = 0
    try:
        if global_load != 0:
            return 1
        init_time = time()
        debug = debug_int != 0
        model_file = ffi.string(model_file_).decode("UTF-8").strip()

        global verbose
        verbose = debug
        global rank
        rank = rank_
        if debug and rank == 0:
            print("init ML ressources", rank, model_file, debug_int, flush=True)

        load_err = load_modules(rank,devID, port)
        if load_err != 0:
            return 10 + load_err

        ctype_err = build_ctypes_converters()
        if ctype_err != 0:
            return 20 + ctype_err
        
        global ml_config
        ml_config = {}
        if os.path.exists("mlconfig.yaml"):
            with open("mlconfig.yaml", "r") as f:
                ml_config = yaml.safe_load(f)
            print("mlconfig.yaml", ml_config, flush=True)

        # if debug:
        #     print(f'TINKER_ML_DIR="{_model_dir}"')

        ## load model
        global model
        model = fennol.load(model_file)
        
        global energy_multiplier
        energy_multiplier = au.KCALPERMOL / model.Ha_to_model_energy

        global qtot, qligand
        qtot = float(qtot_)
        qligand = float(qligand_)

        global gc_stride
        gc_stride = gc_stride_

        global first_eval
        first_eval = True

        global ncalls
        ncalls = 0

        use_lambda = use_lambda_int != 0
        gradient_keys = ["coordinates", "cells"]
        if use_lambda:
            gradient_keys = gradient_keys + ["alch_elambda", "alch_vlambda"]
        energy_and_gradient = model.get_gradient_function(
            *gradient_keys, jit=False, variables_as_input=True
        )

        global energy_and_gradient_and_virial
        @jax.jit
        def energy_and_gradient_and_virial(variables,inputs):
            e, de, output = energy_and_gradient(variables, inputs)
            dedx = de["coordinates"]
            dedcells = de["cells"][0,:,:]
            cells = inputs["cells"][0,:,:]
            x = inputs["coordinates"]
            batch_index = inputs["batch_index"]
            vir = jax.ops.segment_sum(
                dedx[:, :, None] * x[:, None, :],
                batch_index,
                num_segments=len(inputs["natoms"]),
            ) + jnp.matmul(dedcells, cells.T)

            return e, de, vir, output

        nat = int(nat)
        species = asarray(species_ptr, [nat])

        global species_set, species_count
        species_set, species_count = np.unique(species, return_counts=True)

        print("species_set", species_set.tolist(), "species_count", species_count.tolist(), flush=True)

        if use_lambda:
            nligand = int(nligand)
            index_ligand = asarray(index_ligand_ptr, [nligand]) - 1
            print("ligand indices",index_ligand.tolist(), flush=True)
            species_ligand = species[index_ligand]
            global species_ligand_count
            species_ligand_set, species_ligand_count_ = np.unique(species_ligand, return_counts=True)
            species_ligand_count = np.zeros(len(species_set), dtype=int)
            for i, s in enumerate(species_set):
                for j,sl in enumerate(species_ligand_set):
                    if s == sl:
                        species_ligand_count[i] = species_ligand_count_[j]
                        break
            print("species_ligand_count", species_ligand_count.tolist(), flush=True)

        

        if verbose:
            print(f"init ML ressources done {time()-init_time} s".format(), flush=True)
        if rank == 0:
            print("", flush=True)
            print(" *****    Using ML potential engine    ***** ", flush=True)
            print("", flush=True)
    except Exception as err:
        ierr = 1
        if rank < 4:
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", flush=True)
            print(traceback.format_exc(), flush=True)
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<", flush=True)
        # raise Exception('init_ml_ressources:'+str(err))
    return ierr


@ffi.def_extern()
def ml_models(
    coord_ptr,
    atm_ener_ptr,
    gradient_ptr,
    vir_ptr,
    cell_ptr,
    atm_sp_ptr,
    neigh1Idx_ptr,
    neigh2Idx_ptr,
    dist_ptr,
    dxyz_ptr,
    trueat_ptr,
    nat,
    nb_pairs,
    nb_atoms_full,
    dograd_int,
    use_lambda_int,
    alch_elambda,
    alch_vlambda,
    alch_group_ptr,
    dedle_ptr,
    dedlv_ptr,
):
    ierr = 0
    try:
        dograd = dograd_int != 0
        use_lambda = use_lambda_int != 0

        cell = asarray(cell_ptr, [1, 3, 3])
        cell_inv = np.linalg.inv(cell)

        species = asJaxArray(atm_sp_ptr, "int32_t", [nb_atoms_full])
        coordinates = asJaxArray(coord_ptr, "float", [nb_atoms_full, 3])

        true_atoms = asJaxArray(trueat_ptr, "int32_t", [nb_atoms_full]).astype(bool)

        d12 = asJaxArray(dist_ptr, "float", [nb_pairs])
        edge_src = asJaxArray(neigh1Idx_ptr, "int32_t", [nb_pairs])
        edge_dst = asJaxArray(neigh2Idx_ptr, "int32_t", [nb_pairs])
        pbc_shifts = asJaxArray(dxyz_ptr, "float", [nb_pairs, 3])

        total_charge = jnp.asarray(qtot, dtype=jnp.float32)
        if use_lambda:
            alch_group = asJaxArray(alch_group_ptr, "int32_t", [nb_atoms_full])
            alch_elambda = jnp.asarray(float(alch_elambda), dtype=jnp.float32)
            alch_vlambda = jnp.asarray(float(alch_vlambda), dtype=jnp.float32)
            # print("alch_elambda", float(alch_elambda), "alch_vlambda", float(alch_vlambda))
            ligand_charge = jnp.asarray(qligand, dtype=jnp.float32)
        
        graph = {
            "edge_src": edge_src,
            "edge_dst": edge_dst,
            "d12": d12,
            "pbc_shifts": pbc_shifts,
            "overflow": False,
            "keep_graph": True,
        }

        assert nb_atoms_full >= nat, "nb_atoms_full < nat"
        natoms = jnp.asarray([nat], dtype=jnp.int32)
        if nb_atoms_full == nat:
            batch_index = jnp.zeros(nb_atoms_full, dtype=jnp.int32)
        else:
            # natoms = jnp.asarray([nat, nb_atoms_full - nat], dtype=jnp.int32)
            batch_index = jnp.concatenate(
                (
                    jnp.zeros(nat, dtype=jnp.int32),
                    jnp.ones(nb_atoms_full - nat, dtype=jnp.int32),
                )
            )
        
        global inputs, first_eval
        if first_eval:
            first_eval = False
            inputs = {
                **ml_config,
                "total_charge": total_charge,
                "species": species,
                "coordinates": coordinates,
                "natoms": natoms,
                "batch_index": batch_index,
                "cells": cell,
                "reciprocal_cells": cell_inv,
                "graph": graph,
                "recompute_species_index": True,
                "true_atoms": true_atoms,
                "species_set": species_set,
                "species_count": species_count,
                "flags": {"recompute_species_index": None},
            }
            if use_lambda:
                inputs["alch_vlambda"] = alch_vlambda
                inputs["alch_elambda"] = alch_elambda
                inputs["alch_group"] = alch_group
                inputs["alch_ligand_charge"] = ligand_charge
                inputs["species_ligand_count"] = species_ligand_count
            inputs = model.preprocess(**inputs)
            inputs_ = inputs

        else:
            inputs["batch_index"] = batch_index
            inputs["true_atoms"] = true_atoms
            inputs["species"] = species
            inputs["coordinates"] = coordinates
            inputs["cells"] = cell
            inputs["reciprocal_cells"] = cell_inv
            inputs["graph"] = graph
            inputs["total_charge"] = total_charge
            inputs["natoms"] = natoms
            if use_lambda:
                inputs["alch_vlambda"] = alch_vlambda
                inputs["alch_elambda"] = alch_elambda
                inputs["alch_group"] = alch_group
                inputs["alch_ligand_charge"] = ligand_charge
            inputs_ = model.preprocessing.process(model.preproc_state, inputs)
            model.preproc_state, state_up, inputs_, overflow = (
                model.preprocessing.check_reallocate(model.preproc_state, inputs_)
            )
            # print(rank,model.preproc_state,state_up,overflow)
            if overflow:
                print("rank", rank, ": nblist overflow => reallocating nblist")
                print("   size updates:", state_up)
        
        if dograd:
            _, de, vir, output = energy_and_gradient_and_virial(model.variables, inputs_)
            dedx = de["coordinates"]
            # _, f, output = model._energy_and_forces(model.variables, inputs)
            # jax.debug.print(
            #     f"rank {rank} nat {nat} nb_atoms_full {nb_atoms_full}" + " {a} \n {b}",
            #     a=f, b=true_atoms,
            # )
            # gradient = asGPUarray(gradient_ptr, "float", [3 * nat])
            writeGPUptr(gradient_ptr, energy_multiplier * dedx.flatten(), wait=True)
            writeGPUptr(vir_ptr, energy_multiplier * vir.flatten(), wait=True)
            if use_lambda:
                dedle = energy_multiplier * de["alch_elambda"]
                dedlv = energy_multiplier * de["alch_vlambda"]
                writeGPUptr(dedle_ptr, dedle, wait=True)
                writeGPUptr(dedlv_ptr, dedlv, wait=True)
                # jax.debug.print("dedle {deldle} dedlv {dedlv}", deldle=dedle, dedlv=dedlv)
            
            del de, vir, dedx

        else:
            _, output = model._total_energy(model.variables, inputs_)

        writeGPUptr(
            atm_ener_ptr, energy_multiplier * output["atomic_energies"], wait=True
        )

        del output, inputs_
        global ncalls
        ncalls += 1
        if ncalls % gc_stride == 0:
            gc.collect()
        
    except Exception as err:
        ierr = 1
        if rank ==0:
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", flush=True)
            print(traceback.format_exc(), flush=True)
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<", flush=True)
        # raise Exception("ml_models "+str(err))
    return ierr
