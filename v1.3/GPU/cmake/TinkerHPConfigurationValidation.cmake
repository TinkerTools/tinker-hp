cmake_minimum_required(VERSION 3.31)

if(NOT DEFINED TINKERHP_SOURCE_DIR)
  get_filename_component(TINKERHP_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
else()
  get_filename_component(TINKERHP_SOURCE_DIR "${TINKERHP_SOURCE_DIR}" ABSOLUTE)
endif()
if(NOT DEFINED TINKERHP_VALIDATION_DIR)
  set(TINKERHP_VALIDATION_DIR "${CMAKE_CURRENT_BINARY_DIR}/cmake-configuration-validation")
endif()

function(tinkerhp_validate_configuration name precision link_mode fpa ortho colvars nvtx no_mutation suffix required_pattern forbidden_pattern)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=${precision}
      -DTINKERHP_LINK_MODE=${link_mode}
      -DTINKERHP_ENABLE_FPA=${fpa}
      -DTINKERHP_ENABLE_ORTHO_BOX_ONLY=${ortho}
      -DTINKERHP_ENABLE_COLVARS=${colvars}
      -DTINKERHP_ENABLE_NVTX=${nvtx}
      -DTINKERHP_ENABLE_NO_MUTATION=${no_mutation}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} configuration failed")
  endif()

  file(READ "${build_dir}/tinkerhp-build-report.txt" report)
  string(FIND "${report}" "Link mode: ${link_mode}" link_mode_index)
  if(link_mode_index EQUAL -1)
    message(FATAL_ERROR "${name} did not report link mode ${link_mode}")
  endif()
  string(FIND "${report}" "Binary suffix: ${suffix}" suffix_index)
  if(suffix_index EQUAL -1)
    message(FATAL_ERROR "${name} did not report suffix ${suffix}")
  endif()

  file(READ "${build_dir}/compile_commands.json" commands)
  if(NOT required_pattern STREQUAL "")
    string(FIND "${commands}" "${required_pattern}" required_index)
    if(required_index EQUAL -1)
      message(FATAL_ERROR "${name} is missing ${required_pattern}")
    endif()
  endif()
  if(NOT forbidden_pattern STREQUAL "")
    string(FIND "${commands}" "${forbidden_pattern}" forbidden_index)
    if(NOT forbidden_index EQUAL -1)
      message(FATAL_ERROR "${name} unexpectedly contains ${forbidden_pattern}")
    endif()
  endif()

  if(link_mode STREQUAL "SHARED")
    file(READ "${build_dir}/build.ninja" ninja)
    string(FIND "${ninja}" "TARGET_FILE = lib/libtinkerhp_core${CMAKE_SHARED_LIBRARY_SUFFIX}" shared_core_index)
    if(shared_core_index EQUAL -1)
      message(FATAL_ERROR "${name} does not build tinkerhp_core as a shared library")
    endif()
    string(FIND "${ninja}" "Fortran_SHARED_LIBRARY_LINKER__tinkerhp_core" shared_linker_index)
    if(shared_linker_index EQUAL -1)
      message(FATAL_ERROR "${name} does not link the shared core with the Fortran MPI driver")
    endif()
    string(FIND "${ninja}" "-acc" shared_acc_index)
    string(FIND "${ninja}" "-cuda" shared_cuda_index)
    if(shared_acc_index EQUAL -1 OR shared_cuda_index EQUAL -1)
      message(FATAL_ERROR "${name} does not apply the GPU link flags to the shared core")
    endif()
    string(FIND "${commands}" "-fPIC" pic_index)
    if(pic_index EQUAL -1)
      message(FATAL_ERROR "${name} does not compile shared-core objects as PIC")
    endif()
    foreach(target
        tinkerhp_modules
        tinkerhp_fortran_common
        tinkerhp_fortran_gpu
        tinkerhp_cuda_fortran
        tinkerhp_cxx
        tinkerhp_2decomp_fft)
      string(FIND "${ninja}" "target ${target}\n" target_index)
      if(target_index EQUAL -1)
        message(FATAL_ERROR "${name} is missing the ${target} build rules")
      endif()
      string(SUBSTRING "${ninja}" ${target_index} 5000 target_rules)
      string(FIND "${target_rules}" "-fPIC" target_pic_index)
      if(target_pic_index EQUAL -1)
        message(FATAL_ERROR "${name} does not compile ${target} as PIC")
      endif()
    endforeach()
    file(READ "${build_dir}/cmake_install.cmake" install_script)
    string(FIND "${install_script}" "tinkerhp_core${CMAKE_SHARED_LIBRARY_SUFFIX}" shared_install_index)
    if(shared_install_index EQUAL -1)
      message(FATAL_ERROR "${name} does not install the shared core")
    endif()
    string(FIND "${ninja}" "ORIGIN/../${CMAKE_INSTALL_LIBDIR}" rpath_index)
    if(rpath_index EQUAL -1)
      message(FATAL_ERROR "${name} does not link executables with a relative shared-core RPATH")
    endif()
    if(colvars STREQUAL "ON")
      string(FIND "${ninja}" "CXXFLAGS=" colvars_cxxflags_index)
      if(colvars_cxxflags_index EQUAL -1)
        message(FATAL_ERROR "${name} does not retain the Colvars external archive build")
      endif()
      string(SUBSTRING "${ninja}" ${colvars_cxxflags_index} 1000 colvars_build_command)
      string(FIND "${colvars_build_command}" "-fPIC" colvars_pic_index)
      if(colvars_pic_index EQUAL -1)
        message(FATAL_ERROR "${name} does not compile the Colvars archive as PIC")
      endif()
    endif()
  endif()
endfunction()

function(tinkerhp_assert_shared_colvars_acceptance_registration)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/shared-colvars-acceptance")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=MIXED
      -DTINKERHP_LINK_MODE=SHARED
      -DTINKERHP_ENABLE_FPA=ON
      -DTINKERHP_ENABLE_COLVARS=ON
      -DTINKERHP_GPU_ARCHITECTURES=75
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "shared Colvars acceptance-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  foreach(expected
      "tinkerhp-shared-colvars-acceptance"
      "\"mpirun\" \"-n\" \"1\""
      "bin/dynamic_clv.fpa\" \"dhfr2dyn\" \"10\" \"2\" \"20\" \"2\" \"300\" \"-k\" \"dhfr2dyn_r\""
      "WORKING_DIRECTORY \"${TINKERHP_SOURCE_DIR}/examples\"")
    string(FIND "${tests}" "${expected}" expected_index)
    if(expected_index EQUAL -1)
      message(FATAL_ERROR "shared Colvars acceptance test is missing ${expected}")
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_parallel_scientific_registration name precision)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DBUILD_TESTING=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=${precision}
      -DTINKERHP_ENABLE_PARALLEL_SCIENTIFIC_TESTS=ON
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} parallel scientific-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  if(precision STREQUAL "DOUBLE")
    set(expected_tests
      "tinkerhp-parallel-scientific-dynamic-double-auto-4"
      "--precision\" \"DOUBLE\""
      "--mpi-ranks\" \"4\""
      "dhfr2dyn_r.double.dyn.ref"
      "dhfr2dyn_r.key")
    set(forbidden_test "tinkerhp-parallel-scientific-dynamic-mixed")
  else()
    set(expected_tests
      "tinkerhp-parallel-scientific-dynamic-mixed-auto-2"
      "tinkerhp-parallel-scientific-dynamic-mixed-1d-4"
      "--mpi-ranks\" \"2\""
      "--mpi-ranks\" \"4\""
      "dhfr2dyn_r_1d.key")
    set(forbidden_test "tinkerhp-parallel-scientific-dynamic-double")
  endif()
  list(APPEND expected_tests
    "\"mpirun\" \"-n\""
    "RUN_SERIAL \"TRUE\" TIMEOUT \"180\""
    "parallel-scientific-test-contract")
  foreach(expected IN LISTS expected_tests)
    string(FIND "${tests}" "${expected}" expected_index)
    if(expected_index EQUAL -1)
      message(FATAL_ERROR "${name} parallel scientific test is missing ${expected}")
    endif()
  endforeach()
  string(FIND "${tests}" "${forbidden_test}" forbidden_index)
  if(NOT forbidden_index EQUAL -1)
    message(FATAL_ERROR "${name} registered a parallel scientific test for another precision")
  endif()
endfunction()

function(tinkerhp_assert_bench6_scientific_registration name precision fpa expected reference)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DBUILD_TESTING=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=${precision}
      -DTINKERHP_ENABLE_FPA=${fpa}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} Bench6 scientific-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  string(REGEX MATCHALL
    "add_test\\(\\[=\\[tinkerhp-scientific-dynamic-bench6-[a-z]+\\]=\\]"
    registered "${tests}")
  list(LENGTH registered registration_count)
  if(expected STREQUAL "")
    if(NOT registration_count EQUAL 0)
      message(FATAL_ERROR "${name} unexpectedly registered a Bench6 scientific test")
    endif()
    return()
  endif()
  string(FIND "${registered}" "${expected}" expected_index)
  if(NOT registration_count EQUAL 1 OR expected_index EQUAL -1)
    message(FATAL_ERROR "${name} did not register only ${expected}")
  endif()
  foreach(required
      "--precision\" \"${precision}\""
      "--mpi-ranks\" \"1\""
      "bench6_r.key"
      "bench6.xyz"
      "bench6.dyn"
      "${reference}"
      "bench6\" \"50\" \"2\" \"20\" \"2\" \"300\" \"-k\" \"bench6_r")
    string(FIND "${tests}" "${required}" required_index)
    if(required_index EQUAL -1)
      message(FATAL_ERROR "${name} Bench6 scientific test is missing ${required}")
    endif()
  endforeach()
  if(fpa STREQUAL "ON")
    string(FIND "${tests}" "--feature\" \"FPA\"" feature_index)
    if(feature_index EQUAL -1)
      message(FATAL_ERROR "${name} Bench6 FPA test is missing FPA build evidence")
    endif()
    string(FIND "${tests}" "--dynamic-energy-threshold\" \"1e-4\"" threshold_index)
    if(threshold_index EQUAL -1)
      message(FATAL_ERROR "${name} Bench6 FPA test is missing its declared threshold")
    endif()
  endif()
endfunction()

function(tinkerhp_assert_performance_registration name precision fpa enabled expected)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DBUILD_TESTING=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=${precision}
      -DTINKERHP_ENABLE_FPA=${fpa}
      -DTINKERHP_ENABLE_PERFORMANCE_TESTS=${enabled}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} performance-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  string(REGEX MATCHALL
    "add_test\\(\\[=\\[tinkerhp-performance-dynamic-dhfr-b-fpa\\]=\\]"
    registered "${tests}")
  list(LENGTH registered registration_count)
  if(expected STREQUAL "")
    if(NOT registration_count EQUAL 0)
      message(FATAL_ERROR "${name} unexpectedly registered a performance test")
    endif()
    return()
  endif()
  if(NOT registration_count EQUAL 1)
    message(FATAL_ERROR "${name} did not register exactly one performance test")
  endif()
  set(test_marker "add_test([=[tinkerhp-performance-dynamic-dhfr-b-fpa]=]")
  set(properties_marker
    "set_tests_properties([=[tinkerhp-performance-dynamic-dhfr-b-fpa]=]")
  string(FIND "${tests}" "${test_marker}" test_start)
  string(FIND "${tests}" "${properties_marker}" properties_start)
  if(test_start EQUAL -1 OR properties_start EQUAL -1)
    message(FATAL_ERROR "${name} performance test registration is incomplete")
  endif()
  string(SUBSTRING "${tests}" ${test_start} -1 test_tail)
  string(FIND "${test_tail}" "\n" test_line_length)
  string(SUBSTRING "${test_tail}" 0 ${test_line_length} test_record)
  string(SUBSTRING "${tests}" ${properties_start} -1 properties_tail)
  string(FIND "${properties_tail}" "\n" properties_line_length)
  string(SUBSTRING "${properties_tail}" 0 ${properties_line_length} properties_record)
  set(registration "${test_record}\n${properties_record}")
  foreach(required
      "--profile" "dynamic-performance"
      "--precision" "MIXED"
      "--mpi-ranks" "1"
      "--feature" "FPA"
      "dhfr2dyn_B_perf.key"
      "dhfr2dyn_B_perf.fpa.dyn.ref"
      "5000" "-k" "dhfr2dyn_B_perf"
      "LABELS" "performance"
      "RUN_SERIAL" "TRUE" "TIMEOUT" "900")
    string(FIND "${registration}" "${required}" required_index)
    if(required_index EQUAL -1)
      message(FATAL_ERROR "${name} performance test is missing ${required}")
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_watersmall_scientific_registration name precision expected reference)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DBUILD_TESTING=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=${precision}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} Watersmall AMOEBA+ scientific-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  string(REGEX MATCHALL
    "add_test\\(\\[=\\[${expected}\\]=\\]"
    registered "${tests}")
  list(LENGTH registered registration_count)
  if(expected STREQUAL "")
    if(NOT registration_count EQUAL 0)
      message(FATAL_ERROR "${name} unexpectedly registered a Watersmall AMOEBA+ scientific test")
    endif()
    return()
  endif()
  if(NOT registration_count EQUAL 1)
    message(FATAL_ERROR "${name} did not register only ${expected}")
  endif()
  foreach(required
      "--precision\" \"${precision}\""
      "--mpi-launcher\" \"mpirun\""
      "--mpi-ranks\" \"1\""
      "watersmall_amoebap.key"
      "watersmall.xyz"
      "watersmall.dyn"
      "${reference}"
      "--energy-horizon\" \"50\""
      "--temperature-horizon\" \"20\""
      "watersmall\" \"50\" \"2\" \"200\" \"2\" \"300\" \"-k\" \"watersmall_amoebap.key")
    string(FIND "${tests}" "${required}" required_index)
    if(required_index EQUAL -1)
      message(FATAL_ERROR "${name} Watersmall AMOEBA+ scientific test is missing ${required}")
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_pimd_watersmall_scientific_registration
    name precision fpa parallel expected target_temperature)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DBUILD_TESTING=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=${precision}
      -DTINKERHP_ENABLE_FPA=${fpa}
      -DTINKERHP_ENABLE_PARALLEL_SCIENTIFIC_TESTS=${parallel}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} Watersmall PIMD scientific-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  string(REGEX MATCHALL
    "add_test\\(\\[=\\[tinkerhp(-parallel)?-scientific-pimd-watersmall-[a-z0-9-]+\\]=\\]"
    registered "${tests}")
  list(LENGTH registered registration_count)
  if(NOT registration_count EQUAL ${expected})
    message(FATAL_ERROR "${name} registered ${registration_count} Watersmall PIMD scientific tests; expected ${expected}")
  endif()
  if(expected EQUAL 0)
    return()
  endif()
  string(TOLOWER "${precision}" precision_lower)
  set(expected_tests "tinkerhp-scientific-pimd-watersmall-${precision_lower}")
  if(expected EQUAL 2)
    list(APPEND expected_tests "tinkerhp-parallel-scientific-pimd-watersmall-4")
  endif()
  foreach(test_name IN LISTS expected_tests)
    string(FIND "${tests}" "add_test([=[${test_name}]=]" test_start)
    string(FIND "${tests}" "set_tests_properties([=[${test_name}]=]" properties_start)
    if(test_start EQUAL -1 OR properties_start EQUAL -1)
      message(FATAL_ERROR "${name} is missing the ${test_name} registration or properties")
    endif()
    string(SUBSTRING "${tests}" ${test_start} -1 test_tail)
    string(FIND "${test_tail}" "\n" test_line_length)
    string(SUBSTRING "${test_tail}" 0 ${test_line_length} test_record)
    string(SUBSTRING "${tests}" ${properties_start} -1 properties_tail)
    string(FIND "${properties_tail}" "\n" properties_line_length)
    string(SUBSTRING "${properties_tail}" 0 ${properties_line_length} properties_record)
    if(test_name MATCHES "^tinkerhp-parallel-")
      set(expected_ranks 4)
    else()
      set(expected_ranks 1)
    endif()
    foreach(required
        "--profile\" \"pimd\""
        "--precision\" \"${precision}\""
        "--mpi-launcher\" \"mpirun\""
        "--mpi-ranks\" \"${expected_ranks}\""
        "watersmall.key"
        "watersmall.xyz"
        "watersmall.double.pimd.ref"
        "watersmall\" \"1000\" \"2\" \"20\" \"2\" \"300"
        "--centroid-temperature-target\" \"${target_temperature}")
      string(FIND "${test_record}" "${required}" required_index)
      if(required_index EQUAL -1)
        message(FATAL_ERROR "${name} ${test_name} is missing ${required}")
      endif()
    endforeach()
    foreach(required "PROCESSORS \"${expected_ranks}\"" "RUN_SERIAL \"TRUE\"" "TIMEOUT \"120\"")
      string(FIND "${properties_record}" "${required}" required_index)
      if(required_index EQUAL -1)
        message(FATAL_ERROR "${name} ${test_name} is missing ${required}")
      endif()
    endforeach()
    if(fpa STREQUAL "ON" AND precision STREQUAL "MIXED")
      string(FIND "${test_record}" "--feature\" \"FPA\"" feature_index)
      if(feature_index EQUAL -1)
        message(FATAL_ERROR "${name} ${test_name} is missing FPA build evidence")
      endif()
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_analyze_beads_watersmall_functional_registration name precision fpa expected)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DBUILD_TESTING=ON
      -DTINKERHP_PRECISION=${precision}
      -DTINKERHP_ENABLE_FPA=${fpa}
    RESULT_VARIABLE configure_result
    OUTPUT_QUIET
    ERROR_QUIET)
  if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR "${name} Watersmall analyze_beads functional-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  string(REGEX MATCHALL
    "add_test\\(\\[=\\[tinkerhp-functional-analyze-beads-watersmall-[a-z0-9-]+\\]=\\]"
    registered "${tests}")
  list(LENGTH registered registration_count)
  if(NOT registration_count EQUAL ${expected})
    message(FATAL_ERROR "${name} registered ${registration_count} Watersmall analyze_beads functional tests; expected ${expected}")
  endif()
  if(expected EQUAL 0)
    return()
  endif()
  foreach(required
      "--profile\" \"analyze-beads\""
      "--precision\" \"${precision}\""
      "--mpi-launcher\" \"mpirun\""
      "watersmall.key"
      "watersmall.xyz"
      "--bead-archive-source\" \"watersmall.xyz\""
      "--bead-archive-count\" \"32\""
      "watersmall\" \"E\" \"300")
    string(FIND "${tests}" "${required}" required_index)
    if(required_index EQUAL -1)
      message(FATAL_ERROR "${name} Watersmall analyze_beads functional test is missing ${required}")
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_pimd_piqtb_functional_registration)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/pimd-piqtb-parallel")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DBUILD_TESTING=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=MIXED
      -DTINKERHP_ENABLE_FPA=OFF
      -DTINKERHP_ENABLE_PARALLEL_SCIENTIFIC_TESTS=ON
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "PIQTB functional-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  foreach(required
      "tinkerhp-parallel-functional-pimd-piqtb-watersmall-4"
      "--profile\" \"pimd-piqtb\""
      "--mpi-ranks\" \"4\""
      "watersmall_pimd.key"
      "watersmall\" \"1\" \"2\" \"2\" \"2\" \"300\" \"-k\" \"watersmall_pimd")
    string(FIND "${tests}" "${required}" required_index)
    if(required_index EQUAL -1)
      message(FATAL_ERROR "PIQTB functional test is missing ${required}")
    endif()
  endforeach()
  string(FIND "${tests}" "--feature\" \"FPA\"" feature_index)
  if(NOT feature_index EQUAL -1)
    message(FATAL_ERROR "PIQTB functional test must not require FPA")
  endif()
endfunction()

function(tinkerhp_assert_iceii_minimize_scientific_registration name precision expected reference)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DBUILD_TESTING=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=${precision}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} IceII Minimize scientific-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  string(REGEX MATCHALL
    "add_test\\(\\[=\\[tinkerhp-scientific-minimize-iceii-double\\]=\\]"
    registered "${tests}")
  list(LENGTH registered registration_count)
  if(expected STREQUAL "")
    if(NOT registration_count EQUAL 0)
      message(FATAL_ERROR "${name} unexpectedly registered an IceII Minimize scientific test")
    endif()
    return()
  endif()
  if(NOT registration_count EQUAL 1)
    message(FATAL_ERROR "${name} did not register exactly one IceII Minimize scientific test")
  endif()
  foreach(required
      "--profile\" \"minimize\""
      "--precision\" \"DOUBLE\""
      "--mpi-launcher\" \"mpirun\""
      "--mpi-ranks\" \"1\""
      "iceII-large-amoeba_b.key"
      "iceII-large-amoeba.xyz"
      "${reference}"
      "iceII-large-amoeba.xyz\" \"1.0\" \"-k\" \"iceII-large-amoeba_b")
    string(FIND "${tests}" "${required}" required_index)
    if(required_index EQUAL -1)
      message(FATAL_ERROR "${name} IceII Minimize scientific test is missing ${required}")
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_parallel_bench6_fpa_registration name fpa parallel expected)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DBUILD_TESTING=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=MIXED
      -DTINKERHP_ENABLE_FPA=${fpa}
      -DTINKERHP_ENABLE_PARALLEL_SCIENTIFIC_TESTS=${parallel}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} parallel Bench6 FPA scientific-test configuration failed")
  endif()

  file(READ "${build_dir}/CTestTestfile.cmake" tests)
  string(REGEX MATCHALL
    "add_test\\(\\[=\\[tinkerhp-parallel-scientific-dynamic-bench6-fpa-2\\]=\\]"
    registered "${tests}")
  list(LENGTH registered registration_count)
  if(NOT expected)
    if(NOT registration_count EQUAL 0)
      message(FATAL_ERROR "${name} unexpectedly registered parallel Bench6 FPA")
    endif()
    return()
  endif()
  if(NOT registration_count EQUAL 1)
    message(FATAL_ERROR "${name} did not register parallel Bench6 FPA exactly once")
  endif()
  foreach(required
      "tinkerhp-parallel-scientific-dynamic-bench6-fpa-2"
      "--mpi-ranks\" \"2\""
      "--feature\" \"FPA\""
      "--dynamic-energy-threshold\" \"1e-4\""
      "--dynamic-temperature-target\" \"300\""
      "bench6_r.fpa.dyn.ref"
      "bench6\" \"50\" \"2\" \"20\" \"2\" \"300\" \"-k\" \"bench6_r"
      "parallel-scientific-test-contract")
    string(FIND "${tests}" "${required}" required_index)
    if(required_index EQUAL -1)
      message(FATAL_ERROR "parallel Bench6 FPA scientific test is missing ${required}")
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_fortran_ninja_preprocess_definitions name source)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  file(READ "${build_dir}/build.ninja" ninja)
  string(FIND "${ninja}" "${source}-pp.f90 |" source_index)
  if(source_index EQUAL -1)
    message(FATAL_ERROR "${name} is missing the preprocessing rule for ${source} in build.ninja")
  endif()
  string(SUBSTRING "${ninja}" ${source_index} 1000 preprocess_rule)
  foreach(definition IN LISTS ARGN)
    string(REGEX MATCH "DEFINES = [^\n]*${definition}" definition_match "${preprocess_rule}")
    if(NOT definition_match)
      message(FATAL_ERROR "${name} does not pass ${definition} to the Fortran preprocessing scan for ${source}")
    endif()
  endforeach()

  string(FIND "${ninja}" "${source}.o: Fortran_COMPILER" source_index)
  if(source_index EQUAL -1)
    message(FATAL_ERROR "${name} is missing the compile rule for ${source} in build.ninja")
  endif()
  string(SUBSTRING "${ninja}" ${source_index} 1000 compile_rule)
  string(REGEX MATCH "FLAGS = [^\n]*(-cpp|-D[A-Za-z0-9_])" leaked_preprocessor_flags "${compile_rule}")
  if(leaked_preprocessor_flags)
    message(FATAL_ERROR "${name} passes preprocessor flags to compilation of ${source}-pp.f90: ${leaked_preprocessor_flags}")
  endif()
endfunction()

function(tinkerhp_assert_fortran_ninja_preprocess_definitions_absent name source)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  file(READ "${build_dir}/build.ninja" ninja)
  string(FIND "${ninja}" "${source}-pp.f90 |" source_index)
  if(source_index EQUAL -1)
    message(FATAL_ERROR "${name} is missing the preprocessing rule for ${source} in build.ninja")
  endif()
  string(SUBSTRING "${ninja}" ${source_index} 1000 preprocess_rule)
  foreach(definition IN LISTS ARGN)
    string(FIND "${preprocess_rule}" "${definition}" definition_index)
    if(NOT definition_index EQUAL -1)
      message(FATAL_ERROR "${name} unexpectedly passes ${definition} to the Ninja preprocessing scan for ${source}")
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_fortran_makefile_preprocess_definitions)
  set(build_dir "${TINKERHP_VALIDATION_DIR}/fortran-makefile-preprocess")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G "Unix Makefiles"
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_PRECISION=MIXED
      -DTINKERHP_LINK_MODE=SHARED
      -DTINKERHP_ENABLE_FPA=ON
      -DTINKERHP_ENABLE_COLVARS=ON
      -DTINKERHP_GPU_ARCHITECTURES=75
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Unix Makefiles Fortran preprocessing configuration failed")
  endif()

  set(makefile_flag_files
    "2decomp_fft/CMakeFiles/tinkerhp_2decomp_fft.dir/flags.make"
    "CMakeFiles/tinkerhp_fortran_gpu.dir/flags.make"
    "CMakeFiles/tinkerhp_cuda_fortran.dir/flags.make"
  )
  set(makefile_definitions -D_OPENACC -D_CUDA -D_CUDA)
  foreach(flag_file definition IN ZIP_LISTS makefile_flag_files makefile_definitions)
    file(READ "${build_dir}/${flag_file}" flags)
    string(FIND "${flags}" "${definition}" definition_index)
    if(definition_index EQUAL -1)
      message(FATAL_ERROR
        "Unix Makefiles does not pass ${definition} to the Fortran preprocessing scan in ${flag_file}"
      )
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_nvhpc_makefile_fortran_launcher)
  set(fixture_dir "${TINKERHP_VALIDATION_DIR}/nvhpc-makefile-fortran-launcher")
  set(source_dir "${fixture_dir}/source")
  set(build_dir "${fixture_dir}/build")
  set(launcher "${fixture_dir}/launcher.sh")
  set(fake_compiler "${fixture_dir}/fake-fortran-compiler.sh")
  set(arguments_file "${fixture_dir}/compiler-arguments.txt")
  file(REMOVE_RECURSE "${fixture_dir}")
  file(MAKE_DIRECTORY "${source_dir}" "${build_dir}")

  set(TINKERHP_MAKEFILE_LAUNCHER_SOURCE_DIR "${source_dir}")
  configure_file(
    "${TINKERHP_SOURCE_DIR}/cmake/TinkerHPNVHPCMakeFortranLauncher.sh.in"
    "${launcher}"
    @ONLY
  )
  file(CHMOD "${launcher}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
  file(WRITE "${fake_compiler}" "#!/bin/sh\nprintf '%s\\n' \"$@\" > '${arguments_file}'\n")
  file(CHMOD "${fake_compiler}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )

  execute_process(
    COMMAND "${launcher}" "${fake_compiler}"
      -acc -D_OPENACC -c "${source_dir}/eliaison1.cu.f90" -o output.o
    WORKING_DIRECTORY "${build_dir}"
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "NVHPC Makefiles Fortran launcher fixture failed")
  endif()

  file(READ "${arguments_file}" arguments)
  string(FIND "${arguments}" "-D_OPENACC" openacc_definition_index)
  if(NOT openacc_definition_index EQUAL -1)
    message(FATAL_ERROR "NVHPC Makefiles Fortran launcher does not remove -D_OPENACC after -acc")
  endif()
  string(FIND "${arguments}" "../source/eliaison1.cu.f90" relative_source_index)
  if(relative_source_index EQUAL -1)
    message(FATAL_ERROR "NVHPC Makefiles Fortran launcher does not relativize project Fortran sources")
  endif()

  execute_process(
    COMMAND "${launcher}" "${fake_compiler}"
      -cuda -D_CUDA -c "${source_dir}/eliaison1.cu.f90" -o output.o
    WORKING_DIRECTORY "${build_dir}"
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "NVHPC Makefiles CUDA Fortran launcher fixture failed")
  endif()
  file(READ "${arguments_file}" arguments)
  string(FIND "${arguments}" "-D_CUDA" cuda_definition_index)
  if(cuda_definition_index EQUAL -1)
    message(FATAL_ERROR "NVHPC Makefiles Fortran launcher unexpectedly removes -D_CUDA")
  endif()
endfunction()

function(tinkerhp_assert_makefile_conditional_module_dependencies)
  set(fixture_dir "${TINKERHP_VALIDATION_DIR}/makefile-conditional-modules")
  set(source_dir "${fixture_dir}/source")
  set(build_dir "${fixture_dir}/build")
  file(REMOVE_RECURSE "${fixture_dir}")
  file(MAKE_DIRECTORY "${source_dir}")
  file(WRITE "${source_dir}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.31)
project(tinkerhp_makefile_conditional_modules LANGUAGES Fortran)
add_library(conditional_modules OBJECT
  cuda_provider.F90
  cuda_consumer.F90
  openacc_provider.F90
  openacc_consumer.F90
)
target_compile_definitions(conditional_modules PRIVATE _CUDA _OPENACC)
]=])
  file(WRITE "${source_dir}/cuda_provider.F90" [=[
#ifdef _CUDA
module utilcu
end module utilcu
#endif
]=])
  file(WRITE "${source_dir}/cuda_consumer.F90" [=[
module echargecu
#ifdef _CUDA
use utilcu
#endif
end module echargecu
]=])
  file(WRITE "${source_dir}/openacc_provider.F90" [=[
#ifdef _OPENACC
module utilacc
end module utilacc
#endif
]=])
  file(WRITE "${source_dir}/openacc_consumer.F90" [=[
module empolegpu
#ifdef _OPENACC
use utilacc
#endif
end module empolegpu
]=])

  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${source_dir}" -B "${build_dir}" -G "Unix Makefiles"
      -DCMAKE_Fortran_COMPILER=gfortran
    RESULT_VARIABLE configure_result
  )
  if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR "Unix Makefiles conditional-module fixture configuration failed")
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${build_dir}"
    RESULT_VARIABLE build_result
  )
  if(NOT build_result EQUAL 0)
    message(FATAL_ERROR "Unix Makefiles conditional-module fixture build failed")
  endif()

  file(READ "${build_dir}/CMakeFiles/conditional_modules.dir/depend.make" dependencies)
  foreach(expected
      "cuda_consumer.F90.o: CMakeFiles/conditional_modules.dir/utilcu.mod.stamp"
      "openacc_consumer.F90.o: CMakeFiles/conditional_modules.dir/utilacc.mod.stamp")
    string(FIND "${dependencies}" "${expected}" dependency_index)
    if(dependency_index EQUAL -1)
      message(FATAL_ERROR
        "Unix Makefiles does not retain the conditional module dependency ${expected}"
      )
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_compatibility_make_submodule_classification)
  set(fixture_dir "${TINKERHP_VALIDATION_DIR}/compatibility-make-submodules")
  file(REMOVE_RECURSE "${fixture_dir}")
  file(MAKE_DIRECTORY "${fixture_dir}")
  file(WRITE "${fixture_dir}/depend.mk" "# fixture dependencies\n")
  file(WRITE "${fixture_dir}/Makefile"
    "include ${TINKERHP_SOURCE_DIR}/source/Makefile.pgi\n\nprint-submodule-inventory:\n\t@:\n")
  foreach(source
      sm_ordinary.f90
      sm_openacc.gpu.f90
      sm_cuda.cu.f90)
    file(WRITE "${fixture_dir}/${source}" "! make-rule fixture\n")
  endforeach()

  execute_process(
    COMMAND make -n -f "${TINKERHP_SOURCE_DIR}/source/Makefile.pgi"
      "src_dir=${fixture_dir}"
      "MOD_SRCS=sm_ordinary.f90 sm_openacc.gpu.f90"
      "DEV_F_SRCS=sm_cuda.cu.f90"
      "RunF77=fortran-openacc"
      "RunCUFC=fortran-cuda"
      "GPUFLAGS=OPENACC_FLAGS"
      "CUFFLAGS=CUDA_FORTRAN_FLAGS"
      build_mod
    WORKING_DIRECTORY "${fixture_dir}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Compatibility Make canonical submodule fixture failed: ${error}")
  endif()

  foreach(expected
      "OPENACC_FLAGS -I${fixture_dir} -c ${fixture_dir}/sm_ordinary.f90 -o sm_ordinary.o"
      "OPENACC_FLAGS -I${fixture_dir} -c ${fixture_dir}/sm_openacc.gpu.f90 -o sm_openacc.gpu.o"
      "CUDA_FORTRAN_FLAGS -I${fixture_dir} -c ${fixture_dir}/sm_cuda.cu.f90 -o cuf.sm_cuda.o")
    string(FIND "${output}" "${expected}" expected_index)
    if(expected_index EQUAL -1)
      message(FATAL_ERROR "Compatibility Make does not select the expected canonical submodule rule: ${expected}")
    endif()
  endforeach()

  execute_process(
    COMMAND make -prRn -f "${fixture_dir}/Makefile"
      "src_dir=${fixture_dir}"
      print-submodule-inventory
    WORKING_DIRECTORY "${fixture_dir}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE make_database
    ERROR_VARIABLE error
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Compatibility Make submodule inventory fixture failed: ${error}")
  endif()
  string(REGEX MATCH
    "build_mod:[^\n]*sm_openacc\\.gpu\\.o[^\n]*sm_ordinary\\.o[^\n]*cuf\\.sm_cuda\\.o"
    module_phase_match "${make_database}")
  if(module_phase_match STREQUAL "")
    message(FATAL_ERROR "Compatibility Make does not place discovered canonical submodules in build_mod")
  endif()

  execute_process(
    COMMAND make -n -f "${TINKERHP_SOURCE_DIR}/source/Makefile.pgi"
      "src_dir=${TINKERHP_SOURCE_DIR}/source"
      "RunF77=fortran-openacc"
      "RunCUFC=fortran-cuda"
      "GPUFLAGS=OPENACC_FLAGS"
      "CUFFLAGS=CUDA_FORTRAN_FLAGS"
      mdstate.gpu.o cuf.MOD_util.o
    WORKING_DIRECTORY "${fixture_dir}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE legacy_output
    ERROR_VARIABLE error
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Compatibility Make legacy filename fixture failed: ${error}")
  endif()
  foreach(expected
      "OPENACC_FLAGS -I${TINKERHP_SOURCE_DIR}/source -c ${TINKERHP_SOURCE_DIR}/source/mdstate.gpu.f90"
      "CUDA_FORTRAN_FLAGS -I${TINKERHP_SOURCE_DIR}/source -c ${TINKERHP_SOURCE_DIR}/source/MOD_util.cu.f90 -o cuf.MOD_util.o")
    string(FIND "${legacy_output}" "${expected}" expected_index)
    if(expected_index EQUAL -1)
      message(FATAL_ERROR "Compatibility Make no longer supports the legacy filename rule: ${expected}")
    endif()
  endforeach()
endfunction()

function(tinkerhp_assert_compatibility_make_submodule_dependencies)
  set(fixture_dir "${TINKERHP_VALIDATION_DIR}/compatibility-make-submodule-dependencies")
  file(REMOVE_RECURSE "${fixture_dir}")
  file(MAKE_DIRECTORY "${fixture_dir}")
  file(COPY "${TINKERHP_SOURCE_DIR}/source/depend.py" DESTINATION "${fixture_dir}")
  file(WRITE "${fixture_dir}/MOD_memory.f90" "module memory_fixture\nend module memory_fixture\n")
  file(WRITE "${fixture_dir}/MOD_parent.f90" "module parent_module\nend module parent_module\n")
  file(WRITE "${fixture_dir}/MOD_shared.f90" "module shared_module\nend module shared_module\n")
  file(WRITE "${fixture_dir}/fixture.inc.f90" "! include dependency fixture\n")
  file(WRITE "${fixture_dir}/sm_ordinary.f90"
    "submodule (PARENT_MODULE) ordinary_implementation\nuse tinmemory\nuse shared_module\n#include \"fixture.inc.f90\"\ncontains\nend submodule ordinary_implementation\n")
  file(WRITE "${fixture_dir}/sm_openacc.gpu.f90"
    "submodule(parent_module) openacc_implementation\ncontains\nend submodule openacc_implementation\n")
  file(WRITE "${fixture_dir}/MOD_device_parent.cu.f90"
    "module device_parent\nend module device_parent\n")
  file(WRITE "${fixture_dir}/sm_cuda.cu.f90"
    "submodule (device_parent) cuda_implementation\ncontains\nend submodule cuda_implementation\n")

  find_program(python3_executable NAMES python3 REQUIRED)
  execute_process(
    COMMAND "${python3_executable}" "${fixture_dir}/depend.py" host
    WORKING_DIRECTORY "${fixture_dir}"
    RESULT_VARIABLE result
    ERROR_VARIABLE error
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Compatibility Make host dependency fixture failed: ${error}")
  endif()
  file(READ "${fixture_dir}/depend.mk" host_dependencies)
  foreach(expected
      "sm_ordinary.o: ${fixture_dir}/fixture.inc.f90 MOD_memory.o MOD_parent.o MOD_shared.o"
      "sm_openacc.gpu.o: MOD_parent.o")
    string(FIND "${host_dependencies}" "${expected}" expected_index)
    if(expected_index EQUAL -1)
      message(FATAL_ERROR "Compatibility Make host dependencies omit ${expected}")
    endif()
  endforeach()
  string(REGEX MATCHALL "sm_ordinary\\.o:[^\n]*" ordinary_rules "${host_dependencies}")
  if(NOT ordinary_rules STREQUAL
      "sm_ordinary.o: ${fixture_dir}/fixture.inc.f90 MOD_memory.o MOD_parent.o MOD_shared.o")
    message(FATAL_ERROR "Compatibility Make emits a duplicate or unexpected ordinary submodule prerequisite")
  endif()
  string(FIND "${host_dependencies}" "cuf." device_index)
  if(NOT device_index EQUAL -1)
    message(FATAL_ERROR "Compatibility Make host dependencies include CUDA Fortran sources")
  endif()

  execute_process(
    COMMAND "${python3_executable}" "${fixture_dir}/depend.py" device
    WORKING_DIRECTORY "${fixture_dir}"
    RESULT_VARIABLE result
    ERROR_VARIABLE error
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Compatibility Make device dependency fixture failed: ${error}")
  endif()
  file(READ "${fixture_dir}/depend.mk" device_dependencies)
  string(FIND "${device_dependencies}"
    "cuf.sm_cuda.o: cuf.MOD_device_parent.o" expected_index)
  if(expected_index EQUAL -1)
    message(FATAL_ERROR
      "Compatibility Make device dependencies omit the CUDA Fortran submodule parent prerequisite"
    )
  endif()
  string(REGEX MATCHALL "cuf\\.sm_cuda\\.o:[^\n]*" cuda_rules "${device_dependencies}")
  if(NOT cuda_rules STREQUAL "cuf.sm_cuda.o: cuf.MOD_device_parent.o")
    message(FATAL_ERROR "Compatibility Make emits a duplicate or unexpected CUDA Fortran submodule prerequisite")
  endif()
endfunction()

function(tinkerhp_assert_preset_entry name kind)
  string(JSON length LENGTH "${presets}" ${kind})
  set(found FALSE)
  math(EXPR last "${length} - 1")
  foreach(index RANGE ${last})
    string(JSON candidate GET "${presets}" ${kind} ${index} name)
    if(candidate STREQUAL name)
      set(found TRUE)
    endif()
  endforeach()
  if(NOT found)
    message(FATAL_ERROR "Missing ${kind} entry for ${name}")
  endif()
endfunction()

function(tinkerhp_assert_preset name)
  foreach(kind configurePresets buildPresets testPresets)
    tinkerhp_assert_preset_entry(${name} ${kind})
  endforeach()
endfunction()

function(tinkerhp_assert_no_preset_prefix prefix kind)
  string(JSON length LENGTH "${presets}" ${kind})
  math(EXPR last "${length} - 1")
  foreach(index RANGE ${last})
    string(JSON candidate GET "${presets}" ${kind} ${index} name)
    string(FIND "${candidate}" "${prefix}" prefix_index)
    if(prefix_index EQUAL 0)
      message(FATAL_ERROR "Unexpected ${kind} entry ${candidate} with private preset prefix ${prefix}")
    endif()
  endforeach()
endfunction()

function(tinkerhp_write_fake_python path version modules)
  file(WRITE "${path}" "#!/bin/sh\nif [ \"$1\" = \"-c\" ]; then\n  case \"$2\" in\n    *version_info*) echo \"${version}\"; exit 0 ;;\n    *sysconfig*) printf '%s\\n%s\\n' '/fake-deephp/lib' 'libpython3.11.so'; exit 0 ;;\n    *) exit ${modules} ;;\n  esac\nfi\nexit 0\n")
  file(CHMOD "${path}" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
endfunction()

function(tinkerhp_validate_nn_configuration name origin version modules cuda expect_success)
  set(fixture_dir "${TINKERHP_VALIDATION_DIR}/deephp-fixtures/${name}")
  set(project_dir "${fixture_dir}/deep-hp")
  set(fake_python "${fixture_dir}/python")
  set(fake_uv "${fixture_dir}/uv")
  set(build_dir "${TINKERHP_VALIDATION_DIR}/${name}")
  file(REMOVE_RECURSE "${fixture_dir}" "${build_dir}")
  file(MAKE_DIRECTORY "${project_dir}")
  file(WRITE "${project_dir}/pyproject.toml" "[project]\nname = 'fake-deep-hp'\n")
  tinkerhp_write_fake_python("${fake_python}" "${version}" "${modules}")

  if(origin STREQUAL "managed-local")
    file(MAKE_DIRECTORY "${project_dir}/.venv/bin")
    file(COPY_FILE "${fake_python}" "${project_dir}/.venv/bin/python" ONLY_IF_DIFFERENT)
    file(CHMOD "${project_dir}/.venv/bin/python" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  elseif(origin STREQUAL "managed-created")
    file(WRITE "${fake_uv}" "#!/bin/sh\nmkdir -p '${project_dir}/.venv/bin'\ncp '${fake_python}' '${project_dir}/.venv/bin/python'\n")
    file(CHMOD "${fake_uv}" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif()

  set(command
    "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${build_dir}" -G Ninja
    -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran
    -DTINKERHP_ENABLE_NN=ON
    -DTINKERHP_PRECISION=MIXED
    -DDEEPHP_PROJECT_DIR=${project_dir}
    -DTINKERHP_INTERNAL_VALIDATE_CUDA_VERSION=${cuda}
  )
  if(origin STREQUAL "explicit" OR origin STREQUAL "managed-local")
    list(APPEND command -DDEEPHP_PYTHON_INTERPRETER=${fake_python})
  endif()
  if(origin STREQUAL "managed-created")
    list(APPEND command -DDEEPHP_UV_EXECUTABLE=${fake_uv})
  endif()
  if(name STREQUAL "nn-colvars")
    list(APPEND command -DTINKERHP_ENABLE_COLVARS=ON)
  endif()
  if(name STREQUAL "nn-composition")
    list(APPEND command
      -DTINKERHP_ENABLE_ORTHO_BOX_ONLY=ON
      -DTINKERHP_ENABLE_NVTX=ON
      -DTINKERHP_ENABLE_NO_MUTATION=ON
    )
  endif()
  execute_process(COMMAND ${command} RESULT_VARIABLE result)
  if(expect_success AND NOT result EQUAL 0)
    message(FATAL_ERROR "${name} configuration failed")
  elseif(NOT expect_success AND result EQUAL 0)
    message(FATAL_ERROR "${name} unexpectedly configured")
  endif()
  if(NOT expect_success)
    return()
  endif()

  file(READ "${build_dir}/tinkerhp-build-report.txt" report)
  string(FIND "${report}" "Deep-HP environment origin: ${origin}" origin_index)
  if(origin_index EQUAL -1)
    message(FATAL_ERROR "${name} did not report ${origin}")
  endif()
  string(FIND "${report}" "Deep-HP Python version: ${version}" version_index)
  if(version_index EQUAL -1)
    message(FATAL_ERROR "${name} did not report Python ${version}")
  endif()
  if(origin STREQUAL "managed-created")
    string(REGEX MATCH "^([0-9]+)\\." _cuda_match "${cuda}")
    string(FIND "${report}" "Deep-HP CUDA profile: cuda${CMAKE_MATCH_1}" profile_index)
    if(profile_index EQUAL -1)
      message(FATAL_ERROR "${name} did not report the managed CUDA profile")
    endif()
  endif()
  if(name STREQUAL "nn-colvars")
    set(expected_suffix "_nn_clv.fpa")
  else()
    set(expected_suffix "_nn.fpa")
  endif()
  string(FIND "${report}" "Binary suffix: ${expected_suffix}" suffix_index)
  if(suffix_index EQUAL -1)
    message(FATAL_ERROR "${name} did not report suffix ${expected_suffix}")
  endif()
  file(READ "${build_dir}/compile_commands.json" commands)
  string(FIND "${commands}" "-DNN_SUPPORT" nn_index)
  if(nn_index EQUAL -1)
    message(FATAL_ERROR "${name} does not emit NN_SUPPORT")
  endif()
  if(name STREQUAL "nn-composition")
    foreach(definition -DORTHOGONAL_BOX_SHAPE_ONLY -DUSE_NVTX -DTINKER_NO_MUTATE)
      string(FIND "${commands}" "${definition}" definition_index)
      if(definition_index EQUAL -1)
        message(FATAL_ERROR "${name} does not compose ${definition} with NN_SUPPORT")
      endif()
    endforeach()
  endif()
  file(READ "${build_dir}/build.ninja" ninja)
  string(FIND "${ninja}" "mlinterface" plugin_index)
  if(plugin_index EQUAL -1)
    message(FATAL_ERROR "${name} does not generate the CFFI plugin in the build graph")
  endif()
  string(FIND "${ninja}" "/fake-deephp/lib/libpython3.11.so" python_library_index)
  if(python_library_index EQUAL -1)
    message(FATAL_ERROR "${name} does not link the Deep-HP Python library")
  endif()
  file(READ "${build_dir}/cmake_install.cmake" install_script)
  string(FIND "${install_script}" "mlinterface" install_index)
  if(install_index EQUAL -1)
    message(FATAL_ERROR "${name} does not install the CFFI plugin")
  endif()
endfunction()

function(tinkerhp_assert_deferred_feature_rejected)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${TINKERHP_VALIDATION_DIR}/deferred-openmp" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_ENABLE_OPENMP=ON
    RESULT_VARIABLE result
  )
  if(result EQUAL 0)
    message(FATAL_ERROR "Deferred OpenMP configuration unexpectedly succeeded")
  endif()
endfunction()

function(tinkerhp_assert_invalid_link_mode_rejected)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${TINKERHP_SOURCE_DIR}" -B "${TINKERHP_VALIDATION_DIR}/invalid-link-mode" -G Ninja
      -DTINKERHP_INTERNAL_VALIDATE_ONLY=ON
      -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran
      -DTINKERHP_LINK_MODE=MODULE
    RESULT_VARIABLE result
  )
  if(result EQUAL 0)
    message(FATAL_ERROR "Invalid link mode configuration unexpectedly succeeded")
  endif()
endfunction()

file(READ "${TINKERHP_SOURCE_DIR}/CMakePresets.json" presets)
foreach(name gpu-double gpu-fpa gpu-fpa-colvars gpu-fpa-nn gpu-fpa-nn-colvars)
  tinkerhp_assert_preset(${name})
endforeach()
foreach(name gpu-fpa-shared gpu-fpa-shared-colvars)
  tinkerhp_assert_preset_entry(${name} configurePresets)
  tinkerhp_assert_preset_entry(${name} buildPresets)
endforeach()
foreach(kind configurePresets buildPresets testPresets)
  tinkerhp_assert_no_preset_prefix(gpu-sm75- ${kind})
endforeach()

file(READ "${TINKERHP_SOURCE_DIR}/CMakeUserPresets.json" presets)
foreach(name
    gpu-sm75-fpa-shared
    gpu-sm75-fpa-shared-colvars
    gpu-sm75-fpa-nn
    gpu-sm75-double
    gpu-sm75-double-colvars
    gpu-sm75-double-plumed)
  tinkerhp_assert_preset_entry(${name} configurePresets)
  tinkerhp_assert_preset_entry(${name} buildPresets)
endforeach()
foreach(name
    gpu-sm75-fpa-shared-colvars
    gpu-sm75-double
    gpu-sm75-double-colvars
    gpu-sm75-double-plumed)
  tinkerhp_assert_preset_entry(${name} testPresets)
endforeach()
foreach(name gpu-sm75-fpa-shared-colvars-make gpu-sm75-double-shared-make)
  tinkerhp_assert_preset(${name})
endforeach()

tinkerhp_validate_configuration(double DOUBLE STATIC ON OFF OFF OFF OFF .gpu "" USE_DETERMINISTIC_REDUCTION)
tinkerhp_validate_configuration(mixed-fpa MIXED STATIC ON ON OFF OFF OFF .fpa ORTHOGONAL_BOX_SHAPE_ONLY "")
tinkerhp_validate_configuration(mixed-fpa-no-ortho MIXED STATIC ON OFF OFF OFF OFF .fpa USE_DETERMINISTIC_REDUCTION ORTHOGONAL_BOX_SHAPE_ONLY)
tinkerhp_validate_configuration(mixed MIXED STATIC OFF OFF OFF OFF OFF .mixed -DMIXED USE_DETERMINISTIC_REDUCTION)
tinkerhp_validate_configuration(mixed-fpa-colvars MIXED STATIC ON OFF ON OFF OFF _clv.fpa USE_DETERMINISTIC_REDUCTION "")
tinkerhp_validate_configuration(mixed-fpa-nvtx-no-mutation MIXED STATIC ON OFF OFF ON ON .fpa -DUSE_NVTX "")
tinkerhp_validate_configuration(double-shared DOUBLE SHARED ON OFF OFF OFF OFF .gpu "" USE_DETERMINISTIC_REDUCTION)
tinkerhp_validate_configuration(double-shared-colvars DOUBLE SHARED ON OFF ON OFF OFF _clv.gpu "" USE_DETERMINISTIC_REDUCTION)
tinkerhp_validate_configuration(mixed-fpa-shared MIXED SHARED ON OFF OFF OFF OFF .fpa USE_DETERMINISTIC_REDUCTION "")
tinkerhp_validate_configuration(mixed-fpa-shared-colvars MIXED SHARED ON OFF ON OFF OFF _clv.fpa USE_DETERMINISTIC_REDUCTION "")
tinkerhp_assert_shared_colvars_acceptance_registration()
tinkerhp_assert_parallel_scientific_registration(parallel-scientific-double DOUBLE)
tinkerhp_assert_parallel_scientific_registration(parallel-scientific-mixed MIXED)
tinkerhp_assert_parallel_bench6_fpa_registration(parallel-bench6-fpa-enabled ON ON TRUE)
tinkerhp_assert_parallel_bench6_fpa_registration(parallel-bench6-fpa-disabled OFF ON FALSE)
tinkerhp_assert_parallel_bench6_fpa_registration(parallel-bench6-option-disabled ON OFF FALSE)
tinkerhp_assert_bench6_scientific_registration(
  bench6-double DOUBLE OFF tinkerhp-scientific-dynamic-bench6-double
  bench6_r.double.dyn.ref)
tinkerhp_assert_bench6_scientific_registration(
  bench6-fpa MIXED ON tinkerhp-scientific-dynamic-bench6-fpa
  bench6_r.fpa.dyn.ref)
tinkerhp_assert_bench6_scientific_registration(bench6-mixed-no-fpa MIXED OFF "" "")
tinkerhp_assert_performance_registration(
  performance-fpa-enabled MIXED ON ON tinkerhp-performance-dynamic-dhfr-b-fpa)
tinkerhp_assert_performance_registration(performance-option-disabled MIXED ON OFF "")
tinkerhp_assert_performance_registration(performance-no-fpa MIXED OFF ON "")
tinkerhp_assert_performance_registration(performance-double DOUBLE ON ON "")
tinkerhp_assert_watersmall_scientific_registration(
  watersmall-double DOUBLE tinkerhp-scientific-dynamic-watersmall-amoebap-double
  watersmall_amoebap.double.dyn.ref)
tinkerhp_assert_watersmall_scientific_registration(watersmall-mixed MIXED "" "")
tinkerhp_assert_pimd_watersmall_scientific_registration(
  pimd-watersmall-double DOUBLE OFF OFF 1
  300)
tinkerhp_assert_pimd_watersmall_scientific_registration(
  pimd-watersmall-fpa MIXED ON OFF 1
  300)
tinkerhp_assert_pimd_watersmall_scientific_registration(
  pimd-watersmall-mixed-no-fpa MIXED OFF OFF 1
  300)
tinkerhp_assert_pimd_watersmall_scientific_registration(
  pimd-watersmall-parallel-mixed-no-fpa MIXED OFF ON 2
  300)
tinkerhp_assert_pimd_watersmall_scientific_registration(
  pimd-watersmall-parallel-fpa MIXED ON ON 2
  300)
tinkerhp_assert_pimd_piqtb_functional_registration()
tinkerhp_assert_analyze_beads_watersmall_functional_registration(
  analyze-beads-watersmall-double DOUBLE OFF 1)
tinkerhp_assert_analyze_beads_watersmall_functional_registration(
  analyze-beads-watersmall-fpa MIXED ON 1)
tinkerhp_assert_analyze_beads_watersmall_functional_registration(
  analyze-beads-watersmall-mixed-no-fpa MIXED OFF 0)
tinkerhp_assert_iceii_minimize_scientific_registration(
  iceii-minimize-double DOUBLE tinkerhp-scientific-minimize-iceii-double
  iceII-large-amoeba_b.double.min.ref)
tinkerhp_assert_iceii_minimize_scientific_registration(iceii-minimize-mixed MIXED "" "")
tinkerhp_assert_fortran_makefile_preprocess_definitions()
tinkerhp_assert_nvhpc_makefile_fortran_launcher()
tinkerhp_assert_makefile_conditional_module_dependencies()
tinkerhp_assert_compatibility_make_submodule_classification()
tinkerhp_assert_compatibility_make_submodule_dependencies()
tinkerhp_validate_nn_configuration(nn-managed-local managed-local 3.11.8 0 12.5 TRUE)
tinkerhp_validate_nn_configuration(nn-explicit explicit 3.10.14 0 12.5 TRUE)
tinkerhp_validate_nn_configuration(nn-created-cuda12 managed-created 3.11.8 0 12.5 TRUE)
tinkerhp_validate_nn_configuration(nn-created-cuda13 managed-created 3.11.8 0 13.1 TRUE)
tinkerhp_validate_nn_configuration(nn-colvars explicit 3.11.8 0 12.5 TRUE)
tinkerhp_validate_nn_configuration(nn-composition explicit 3.11.8 0 12.5 TRUE)
tinkerhp_validate_nn_configuration(nn-managed-wrong-version managed-local 3.10.14 0 12.5 FALSE)
tinkerhp_validate_nn_configuration(nn-explicit-missing-module explicit 3.11.8 1 12.5 FALSE)
tinkerhp_validate_nn_configuration(nn-created-unsupported-cuda managed-created 3.11.8 0 11.8 FALSE)
tinkerhp_assert_deferred_feature_rejected()
tinkerhp_assert_invalid_link_mode_rejected()

# CMake must receive definitions for the preprocessing scan, but compilation of
# the generated *-pp.f90 file must not invoke CPP again.
tinkerhp_assert_fortran_ninja_preprocess_definitions(mixed-fpa
  "source/newinduce_pmegpu.f90" -DMIXED -DUSE_DETERMINISTIC_REDUCTION -DORTHOGONAL_BOX_SHAPE_ONLY)
tinkerhp_assert_fortran_ninja_preprocess_definitions(mixed-fpa
  "src/fft_cufft.f90" -DOVERWRITE)
tinkerhp_assert_fortran_ninja_preprocess_definitions(mixed-fpa-colvars
  "source/newinduce_pmegpu.f90" -DMIXED -DUSE_DETERMINISTIC_REDUCTION -DCOLVARS -DCOLVARS_MPI -DCOLVARS_TCL)
tinkerhp_assert_fortran_ninja_preprocess_definitions_absent(mixed-fpa
  "src/fft_cufft.f90" -D_OPENACC)
tinkerhp_assert_fortran_ninja_preprocess_definitions_absent(mixed-fpa-shared-colvars
  "source/empole1gpu.f90" -D_OPENACC)
tinkerhp_assert_fortran_ninja_preprocess_definitions_absent(mixed-fpa-shared-colvars
  "source/echarge.cu.f90" -D_CUDA)
tinkerhp_assert_fortran_ninja_preprocess_definitions(mixed-fpa-nvtx-no-mutation
  "source/newinduce_pmegpu.f90" -DMIXED -DUSE_DETERMINISTIC_REDUCTION -DUSE_NVTX -DTINKER_NO_MUTATE)

file(READ "${TINKERHP_VALIDATION_DIR}/mixed-fpa-nvtx-no-mutation/build.ninja" nvtx_ninja)
string(FIND "${nvtx_ninja}" "-lnvToolsExt" nvtx_link_index)
if(nvtx_link_index EQUAL -1)
  message(FATAL_ERROR "NVTX configuration does not link nvToolsExt")
endif()

file(READ "${TINKERHP_VALIDATION_DIR}/mixed-fpa/compile_commands.json" ortho_commands)
file(READ "${TINKERHP_VALIDATION_DIR}/mixed-fpa-no-ortho/compile_commands.json" plain_commands)
string(REPLACE "${TINKERHP_VALIDATION_DIR}/mixed-fpa" "<build>" ortho_commands "${ortho_commands}")
string(REPLACE "${TINKERHP_VALIDATION_DIR}/mixed-fpa-no-ortho" "<build>" plain_commands "${plain_commands}")
string(REPLACE "-DORTHOGONAL_BOX_SHAPE_ONLY " "" ortho_commands "${ortho_commands}")
if(NOT ortho_commands STREQUAL plain_commands)
  message(FATAL_ERROR
    "Orthogonal-box-only mode changed compilation beyond ORTHOGONAL_BOX_SHAPE_ONLY."
  )
endif()
