function(tinkerhp_validate_requested_configuration)
  if(NOT TINKERHP_BACKEND STREQUAL "GPU")
    message(FATAL_ERROR "TINKERHP_BACKEND=${TINKERHP_BACKEND} is reserved for a future CPU path. v1 only supports GPU.")
  endif()

  if(NOT TINKERHP_PRECISION STREQUAL "DOUBLE" AND NOT TINKERHP_PRECISION STREQUAL "MIXED")
    message(FATAL_ERROR "TINKERHP_PRECISION=${TINKERHP_PRECISION} is not supported. The maintained CMake path supports DOUBLE and MIXED.")
  endif()

  if(NOT TINKERHP_LINK_MODE STREQUAL "STATIC" AND NOT TINKERHP_LINK_MODE STREQUAL "SHARED")
    message(FATAL_ERROR "TINKERHP_LINK_MODE=${TINKERHP_LINK_MODE} is not supported. Use STATIC or SHARED.")
  endif()

  if(NOT TINKERHP_FFT_BACKEND STREQUAL "GENERIC")
    message(FATAL_ERROR "TINKERHP_FFT_BACKEND=${TINKERHP_FFT_BACKEND} is deferred. v1 only supports GENERIC.")
  endif()

  foreach(_unsupported
      TINKERHP_ENABLE_OPENMP
      TINKERHP_ENABLE_HIP
      TINKERHP_ENABLE_NVSHMEM)
    if(${_unsupported})
      message(FATAL_ERROR "${_unsupported}=ON is part of the deferred CMake roadmap and is intentionally blocked in v1.")
    endif()
  endforeach()
endfunction()

function(tinkerhp_preflight_toolchain out_cuda_version out_workaround_flags)
  if(TINKERHP_INTERNAL_VALIDATE_ONLY)
    if(DEFINED TINKERHP_INTERNAL_VALIDATE_CUDA_VERSION)
      set(_validate_cuda_version "${TINKERHP_INTERNAL_VALIDATE_CUDA_VERSION}")
    else()
      set(_validate_cuda_version "12.5")
    endif()
    set(${out_cuda_version} "${_validate_cuda_version}" PARENT_SCOPE)
    set(${out_workaround_flags} "" PARENT_SCOPE)
    message(STATUS "Tinker-HP internal validate-only mode: skipping NVHPC and CUDA preflight.")
    return()
  endif()

  if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    message(FATAL_ERROR "The maintained GPU CMake path requires an NVHPC Fortran compiler reached through mpif90.")
  endif()

  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    message(FATAL_ERROR "The maintained GPU CMake path requires an NVHPC C++ compiler reached through mpic++.")
  endif()

  if(NOT CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    message(FATAL_ERROR "The maintained GPU CMake path requires nvcc as CMAKE_CUDA_COMPILER.")
  endif()

  get_filename_component(_fc_name "${CMAKE_Fortran_COMPILER}" NAME)
  get_filename_component(_cxx_name "${CMAKE_CXX_COMPILER}" NAME)
  get_filename_component(_cuda_name "${CMAKE_CUDA_COMPILER}" NAME)

  if(NOT _fc_name MATCHES "^mpi")
    message(FATAL_ERROR "CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} is not an MPI wrapper. Use mpif90.")
  endif()
  if(NOT _cxx_name MATCHES "^mpi")
    message(FATAL_ERROR "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} is not an MPI wrapper. Use mpic++ or mpicxx.")
  endif()
  if(NOT _cuda_name STREQUAL "nvcc")
    message(FATAL_ERROR "CMAKE_CUDA_COMPILER=${CMAKE_CUDA_COMPILER} is not nvcc.")
  endif()

  set(_nvhpc_version "${CMAKE_Fortran_COMPILER_VERSION}")
  if(_nvhpc_version VERSION_GREATER_EQUAL "22.9" AND _nvhpc_version VERSION_LESS_EQUAL "24.5")
    message(FATAL_ERROR
      "NVHPC ${_nvhpc_version} is blocked for this CMake path. The documented unstable range is 22.9 through 24.5."
    )
  endif()

  set(_workaround_flags "")
  if(_nvhpc_version VERSION_GREATER_EQUAL "21.3" AND _nvhpc_version VERSION_LESS_EQUAL "21.9")
    if(TINKERHP_ENABLE_NVIDIA_21X_WORKAROUND)
      list(APPEND _workaround_flags "${TINKERHP_NVIDIA_21X_WORKAROUND_FLAG}")
      message(STATUS
        "Applying the NVHPC ${_nvhpc_version} workaround flag ${TINKERHP_NVIDIA_21X_WORKAROUND_FLAG}."
      )
    else()
      message(WARNING
        "NVHPC ${_nvhpc_version} is in the guarded 21.3-21.9 range and the automatic workaround is disabled."
      )
    endif()
  endif()

  set(_probe "${CMAKE_BINARY_DIR}/tinkerhp_cuda_probe.f90")
  file(WRITE "${_probe}" "program tinkerhp_cuda_probe\nend program tinkerhp_cuda_probe\n")
  execute_process(
    COMMAND "${CMAKE_Fortran_COMPILER}" -dryrun -cuda -c "${_probe}"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    RESULT_VARIABLE _probe_result
    OUTPUT_VARIABLE _probe_stdout
    ERROR_VARIABLE _probe_stderr
  )
  if(NOT _probe_result EQUAL 0)
    message(FATAL_ERROR "Unable to probe the NVHPC CUDA runtime version from ${CMAKE_Fortran_COMPILER}.")
  endif()

  set(_probe_output "${_probe_stdout}\n${_probe_stderr}")
  string(REGEX MATCH "cuda([0-9]+\\.[0-9]+)" _cuda_match "${_probe_output}")
  if(CMAKE_MATCH_1)
    set(_fortran_cuda_version "${CMAKE_MATCH_1}")
  else()
    string(REGEX MATCH "/cuda/([0-9]+\\.[0-9]+)" _cuda_path_match "${_probe_output}")
    if(CMAKE_MATCH_1)
      set(_fortran_cuda_version "${CMAKE_MATCH_1}")
    else()
      message(FATAL_ERROR "Unable to extract the NVHPC CUDA runtime version from the compiler dry-run output.")
    endif()
  endif()

  string(REGEX MATCH "^[0-9]+\\.[0-9]+" _nvcc_cuda_version "${CMAKE_CUDA_COMPILER_VERSION}")
  if(NOT _fortran_cuda_version STREQUAL _nvcc_cuda_version)
    message(FATAL_ERROR
      "CUDA version mismatch: NVHPC Fortran toolchain uses ${_fortran_cuda_version} while nvcc reports ${_nvcc_cuda_version}."
    )
  endif()

  set(${out_cuda_version} "${_fortran_cuda_version}" PARENT_SCOPE)
  set(${out_workaround_flags} "${_workaround_flags}" PARENT_SCOPE)
endfunction()
