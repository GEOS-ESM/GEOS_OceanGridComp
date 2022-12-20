esma_set_this ()

esma_add_library (${this} SRCS GEOS_DataSeaGridComp.F90 DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran)

mapl_acg (${this}   GEOS_DataSea_StateSpecs.rc
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
          GET_POINTERS DECLARE_POINTERS)
