if(NOT DEFINED TINKERHP_BUILD_DIR)
  message(FATAL_ERROR "TINKERHP_BUILD_DIR is required")
endif()

set(_required_files
  "${TINKERHP_BUILD_DIR}/compile_commands.json"
)
foreach(_required_file IN LISTS _required_files)
  if(NOT EXISTS "${_required_file}")
    message(FATAL_ERROR "Missing generated build file: ${_required_file}")
  endif()
endforeach()

set(_metadata_files)
foreach(_pattern
    "compile_commands.json"
    "build.ninja"
    "Makefile"
    "FortranModules.json"
    "Fortran.dd"
    "*.d"
    "*.ddi"
    "DependInfo.cmake"
    "depend.make"
    "depend.mk")
  file(GLOB_RECURSE _matches LIST_DIRECTORIES false
    "${TINKERHP_BUILD_DIR}/${_pattern}")
  list(APPEND _metadata_files ${_matches})
endforeach()
list(REMOVE_DUPLICATES _metadata_files)

set(_module_maps)
set(_dependency_files)
set(_dyndep_files)
set(_compiler_dependency_files)
set(_compatibility_dependency_files)
foreach(_metadata_file IN LISTS _metadata_files)
  if(_metadata_file MATCHES "/FortranModules\\.json$")
    list(APPEND _module_maps "${_metadata_file}")
  endif()
  if(_metadata_file MATCHES "(/Fortran\\.dd|\\.d$|\\.ddi$|/DependInfo\\.cmake$|/depend\\.(make|mk)$)")
    list(APPEND _dependency_files "${_metadata_file}")
  endif()
  if(_metadata_file MATCHES "/Fortran\\.dd$")
    list(APPEND _dyndep_files "${_metadata_file}")
  elseif(_metadata_file MATCHES "\\.(d|ddi)$")
    list(APPEND _compiler_dependency_files "${_metadata_file}")
  elseif(_metadata_file MATCHES "/(DependInfo\\.cmake|depend\\.(make|mk))$")
    list(APPEND _compatibility_dependency_files "${_metadata_file}")
  endif()
  file(READ "${_metadata_file}" _metadata)
  string(REGEX MATCH
    "(^|[^A-Za-z0-9])(mpistuff|utilcomm|utilgpu|utilcu)([^A-Za-z0-9]|$)"
    _retired_surface "${_metadata}")
  if(_retired_surface)
    message(FATAL_ERROR
      "Generated metadata restores a retired compatibility surface in "
      "${_metadata_file}: ${_retired_surface}")
  endif()
endforeach()

set(_validation_fixture_prefix
  "${TINKERHP_BUILD_DIR}/cmake-configuration-validation/")
set(_primary_module_maps)
set(_primary_dyndep_files)
set(_primary_compiler_dependency_files)
set(_primary_compatibility_dependency_files)
foreach(_metadata_file IN LISTS _module_maps)
  string(FIND "${_metadata_file}" "${_validation_fixture_prefix}"
    _fixture_index)
  if(_fixture_index EQUAL -1)
    list(APPEND _primary_module_maps "${_metadata_file}")
  endif()
endforeach()
foreach(_metadata_file IN LISTS _dyndep_files)
  string(FIND "${_metadata_file}" "${_validation_fixture_prefix}"
    _fixture_index)
  if(_fixture_index EQUAL -1)
    list(APPEND _primary_dyndep_files "${_metadata_file}")
  endif()
endforeach()
foreach(_metadata_file IN LISTS _compiler_dependency_files)
  string(FIND "${_metadata_file}" "${_validation_fixture_prefix}"
    _fixture_index)
  if(_fixture_index EQUAL -1)
    list(APPEND _primary_compiler_dependency_files "${_metadata_file}")
  endif()
endforeach()
foreach(_metadata_file IN LISTS _compatibility_dependency_files)
  string(FIND "${_metadata_file}" "${_validation_fixture_prefix}"
    _fixture_index)
  if(_fixture_index EQUAL -1)
    list(APPEND _primary_compatibility_dependency_files
      "${_metadata_file}")
  endif()
endforeach()

if(EXISTS "${TINKERHP_BUILD_DIR}/build.ninja")
  if(NOT _primary_module_maps)
    message(FATAL_ERROR
      "No generated FortranModules.json was found; build the configured targets "
      "before running this closure test")
  endif()
  if(NOT _primary_dyndep_files OR NOT _primary_compiler_dependency_files)
    message(FATAL_ERROR
      "Ninja build is missing Fortran dyndep or compiler dependency metadata")
  endif()
elseif(EXISTS "${TINKERHP_BUILD_DIR}/Makefile")
  if(NOT _primary_compatibility_dependency_files)
    message(FATAL_ERROR
      "Makefile build is missing CMake or Make dependency metadata")
  endif()
else()
  message(FATAL_ERROR "Unsupported or incomplete generated build directory")
endif()
message(STATUS
  "Checked ${_metadata_files} generated dependency metadata files; "
  "retired compatibility surfaces absent")
message(STATUS "Primary Fortran module metadata: ${_primary_module_maps}")
message(STATUS "Primary Fortran dyndep metadata: ${_primary_dyndep_files}")
message(STATUS
  "Primary compiler dependency metadata: ${_primary_compiler_dependency_files}")
message(STATUS
  "Primary compatibility dependency metadata: ${_primary_compatibility_dependency_files}")
