esma_set_this ()

#target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

if (GEOS_OCEAN STREQUAL MOM5)
  set (ocean_plug MOM_GEOSPlugMod)
elseif(GEOS_OCEAN STREQUAL MOM6)
  set (ocean_plug MOM6_GEOSPlug)
elseif(GEOS_OCEAN STREQUAL MIT)
  set (ocean_plug MIT_GEOS5PlugMod)
elseif(GEOS_OCEAN STREQUAL datasea)
   # no real ocean model, just datasea
else()
  message (FATAL_ERROR "Invalid ocean model selection: ${GEOS_OCEAN}")
endif()
#option (BUILD_MIT_OCEAN "BUILD_MIT_OCEAN is turned off by default" OFF) #OFF by default
#if ( NOT GEOS_OCEAN EQUAL "MIT")
esma_add_library (${this}
SRCS GEOS_OceanGridComp.F90
SUBCOMPONENTS ${ocean_plug} GEOSdatasea_GridComp
DEPENDENCIES GEOSdatasea_GridComp MAPL esmf)

# Note to SA: If Fahad has problems then add back this logic.

mapl_acg (${this}  GEOS_Ocean_StateSpecs.rc
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
          GET_POINTERS DECLARE_POINTERS)
