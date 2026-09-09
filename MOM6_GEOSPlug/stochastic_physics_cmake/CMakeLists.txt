esma_set_this (OVERRIDE stochastic_physics)

# This list is from the CMakeLists.txt in the stochastic_physics source
#
list(APPEND _stoch_phys_srcs
  kinddef.F90
  mpi_wrapper.F90
  halo_exchange.fv3.F90
  plumes.F90
  mersenne_twister.F90
  random_numbers.F90
  stochy_namelist_def.F90
  compns_stochy.F90
  stochy_internal_state_mod.F90
  spectral_transforms.F90
  stochastic_physics.F90
  stochy_patterngenerator.F90
  stochy_data_mod.F90
  get_stochy_pattern.F90
  cellular_automata_global.F90
  cellular_automata_sgs.F90
  update_ca.F90
  lndp_apply_perts.F90
)

esma_mepo_style(stochastic_physics stochastic_physics_rel_path REL_PATH ..)
set (stochastic_physics_path ${CMAKE_CURRENT_SOURCE_DIR}/${stochastic_physics_rel_path})

set (SRCS)
foreach (file ${_stoch_phys_srcs})
  list (APPEND SRCS ${stochastic_physics_path}/${file})
endforeach ()

find_package(LAPACK REQUIRED)

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES FMS::fms_r8 LAPACK::LAPACK
  TYPE SHARED
)

# The CESMCOUPLED is to prevent it trying to use esmf_dgemm
target_compile_definitions (${this} PRIVATE INTERNAL_FILE_NML CESMCOUPLED)

# This 'resets' the Intel DEBUG flags for MOM6. The stock debug flags use
# 'all,noarg_temp_created' which seem to be too aggressive for FMS/MOM6. This
# moves them back to the 'bounds,uninit' GEOS used to build with.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND CMAKE_BUILD_TYPE MATCHES Debug)
   string(REPLACE "all,noarg_temp_created" "bounds,uninit" _tmp "${GEOS_Fortran_FLAGS_DEBUG}")
   set (CMAKE_Fortran_FLAGS_DEBUG "${_tmp}")
endif ()

# We have to do something similar for ifx
if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND CMAKE_BUILD_TYPE MATCHES Debug)
   string(REPLACE "all,noarg_temp_created" "bounds" _tmp "${GEOS_Fortran_FLAGS_DEBUG}")
   set (CMAKE_Fortran_FLAGS_DEBUG "${_tmp}")
endif ()

# Specs for r8 version
string(REPLACE " " ";" tmp ${FREAL8})
foreach (flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach ()
