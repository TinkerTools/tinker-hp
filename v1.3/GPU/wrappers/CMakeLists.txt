add_library(tinkerhp_wrapper STATIC
  thrust_cache.cu
  csort.cu
  thrust_module.f90
)

set_source_files_properties(thrust_module.f90 PROPERTIES
  Fortran_PREPROCESS OFF
)

set_target_properties(tinkerhp_wrapper PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  CUDA_ARCHITECTURES "${TINKERHP_CUDA_ARCHITECTURES_VALUE}"
  Fortran_MODULE_DIRECTORY "${TINKERHP_WRAPPER_MODULE_DIR}"
)

target_include_directories(tinkerhp_wrapper PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${TINKERHP_WRAPPER_MODULE_DIR}"
)

target_compile_options(tinkerhp_wrapper PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:${TINKERHP_CUDA_FORTRAN_FLAGS}>
  $<$<COMPILE_LANGUAGE:CUDA>:${TINKERHP_CUDA_FLAGS}>
)
