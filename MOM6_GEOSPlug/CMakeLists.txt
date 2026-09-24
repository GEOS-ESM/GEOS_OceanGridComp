esma_set_this ()

set (SRCS
  MOM6_GEOSPlug.F90
  )

esma_mepo_style(stochastic_physics stochastic_physics_rel_path FOUND stochastic_physics_found)
if (stochastic_physics_found)
  esma_add_subdirectories(stochastic_physics_cmake)
  set (MOM6_USE_STOCHASTIC_PHYSICS TRUE)
else ()
  set (MOM6_USE_STOCHASTIC_PHYSICS FALSE)
  message (WARNING "stochastic_physics package not found; using the MOM6 compile-only stubs. The external package is required for runs with stochastic physics enabled.")
endif ()

esma_add_subdirectories(mom6_cmake)

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL FMS::fms_r8 mom6
  TYPE SHARED
  )

# Specs for r8 version
string(REPLACE " " ";" tmp ${FREAL8})
foreach (flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach ()

install(
   DIRECTORY mom6_app
   DESTINATION etc/MOM6
   USE_SOURCE_PERMISSIONS
   )

mapl_acg (${this}   MOM6_GEOSPlug_StateSpecs.rc
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
          GET_POINTERS DECLARE_POINTERS)
