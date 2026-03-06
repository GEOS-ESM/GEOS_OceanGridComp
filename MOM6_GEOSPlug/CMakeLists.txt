esma_set_this ()

set (SRCS
  MOM6_GEOSPlug.F90
  )

esma_add_subdirectories(
  mom6_cmake
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL FMS::fms mom6
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
