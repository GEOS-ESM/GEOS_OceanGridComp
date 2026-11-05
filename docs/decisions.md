# GEOS_OceanGridComp: previous architectural decisions

**Audience:** coding agents and developers working in the `GEOS_OceanGridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOS_OceanGridComp` @ `8d91ab7` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Only what is easy to see in the code now. Entries grow over time. Context:
[`ARCHITECTURE.md`](ARCHITECTURE.md). `OGCM file` is
`GEOSogcm_GridComp/GEOS_OgcmGridComp.F90` in the parent repository `GEOSgcm_GridComp`.

### One thin wrapper component in front of every ocean
- Pattern: `OCEAN` holds exactly one active ocean child (`DATASEA`, `MOM`, `MOM6` or `MIT`), chosen in its `SetServices` from `USE_DATASEA:`
  and `OCEAN_NAME:`; every plug-in presents the same short names (`TW SW UW VW ...`) and `OCEAN` re-exports them as its own with
  `MAPL_AddExportSpec(..., CHILD_ID=OCN)` (`GEOS_OceanGridComp.F90`, `SetServices`).  [Verified]
- Likely reason: the parent `OGCM` and the sea-ice and radiation siblings are written once against `OCEAN`, whichever model runs.  [Inferred]
- Question for maintainers: the plug-ins differ in some names (`MOM_3D_MASK`, `MIT_3D_MASK`, `MOM_2D_MASK`; `UWC`/`VWC` only for `MOM6`; `T_Freeze`
  only for `MOM6`). Which of these differences are part of the `OCEAN` contract?

### Ocean models are plug-ins, not part of this wrapper
- Pattern: each ocean model is wrapped in its own directory and gridded component (`MOM_GEOS5PlugMod`, `MOM6_GEOSPlug`, `MIT_GEOS5PlugMod`); the
  model source is a separate nested repository cloned beneath it (`components.yaml`: `mom`, `mom6`, `stochastic_physics`, `mit` with `local:` under the plug-in
  directory; plug-in `CMakeLists.txt`: `SUBCOMPONENTS @mom`, `esma_mepo_style(mom6 ...)`, `SET (MITGCM_ROOTDIR ... /@mit)`).  [Verified]
- Likely reason: the models have their own upstream histories and build systems, and the plug-in is the only GEOS-specific code.  [Inferred]

### MOM and MOM6 are loaded as shared objects by name; MIT and data ocean are linked statically
- Pattern: `MOM` and `MOM6` are added with `MAPL_AddChild(name, 'setservices_', parentGC=GC, sharedObj=...)`, library names from `MOM_GEOS5PLUGMOD:` /
  `MOM6_GEOSPLUG:`; `MIT` and `DATASEA` are added with `SS=<SetServices>` (`GEOS_OceanGridComp.F90`, `SetServices`). The `MOM` and `MOM6` plug-ins are built
  `TYPE SHARED` (their `CMakeLists.txt`) and each ends in a free `SetServices` subroutine; `MIT_GEOS5PlugMod` is compiled in only with `BUILD_MIT_OCEAN`.  [Verified]
- Likely reason: both FMS-based models (MOM5, MOM6) may not be linkable into one executable, so each is a library selected at run time.  [Inferred; `gcm_setup`
  sets `OCEAN_PRELOAD` to `libmom` or `libmom6` for the chosen model (`src/Applications/@GEOSgcm_App/gcm_setup`, `OCNMODEL` branches), which points the same way]
- Question for maintainers: what led to a build-time switch (`BUILD_MIT_OCEAN`) for `MIT`, which also replaces the sea-ice composite?

### The wrapper weights fluxes by ocean fraction and fills orphan points
- Pattern: before running the child, `OCEAN` multiplies fluxes by `WGHT = FROCEAN/MASK` and copies them into the child's import state; after running, it writes
  child exports back only where `WGHT > 0` and advances the surface temperature of "orphan" points (ocean fraction but no child ocean point) with an explicit
  heat budget over `OGCM_TOP_LAYER:` (`GEOS_OceanGridComp.F90`, `Run`).  [Verified]
- Likely reason: the ocean models and the atmosphere use different land-sea masks; the wrapper makes the exchange conserve what the ocean sees and gives the
  atmosphere a surface temperature everywhere it has ocean.  [Inferred]
- Question for maintainers: which component is meant to own the orphan-point treatment?

### Imports are filled by hand; child imports are terminated
- Pattern: for a real ocean `OCEAN` calls `MAPL_TerminateImport` for the fluxes it will copy itself and fills the child's import state in `Run`
  (`GEOS_OceanGridComp.F90`, `SetServices`, `Run`). The parent `OGCM` likewise terminates most `OCEAN` imports and fills them with `MAPL_LocStreamTransform`
  (`OGCM file`, `SetServices`, `Run`).  [Verified]
- Likely reason: the fluxes must be scaled by `WGHT` and combined with sea-ice terms (for example `SHFLX - FHOCN`, `RAIN + FRESH`), which a plain connection cannot do.  [Inferred]

### Dual ocean: a second data ocean and a two-phase run
- Pattern: `DUAL_OCEAN: 1` adds a second `DATASEA` child, registers a second `Run` phase on `OCEAN`, on `OGCM` and on each of the `MOM`, `MOM6` and `MIT` plug-ins, and computes a
  temperature correction `DEL_TEMP` from the difference between the data-ocean and model surface temperature (`GEOS_OceanGridComp.F90`, `Run`;
  `MOM_GEOS5PlugMod.F90`, `MOM6_GEOSPlug.F90`, `MIT_GEOS5PlugMod.F90`, each `Run2`; `OGCM file`, `Run`).  [Verified]
- Likely reason: nudge the model ocean toward observed SST in the predictor/corrector steps of regular replay.  [Inferred]
- Question for maintainers: is `DUAL_OCEAN` supported for all three ocean models? The template mentions it only as a `MOM5` hint (`@MOM5#DUAL_OCEAN: 1` in
  `src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`).

### Plug-in state specs: ACG for MOM, MOM6 and data ocean; hand-written for MIT
- Pattern: `OCEAN`, `DATASEA`, `MOM` and `MOM6` declare state in `*_StateSpecs.rc` through `mapl_acg`; `MIT_GEOS5PlugMod` lists imports and exports in
  an array of `mstate` records and calls `MAPL_AddImportSpec` / `MAPL_AddExportSpec` in a loop, mirroring every import as an export (`CMakeLists.txt` files;
  `MIT_GEOS5PlugMod.F90`, `SetServices`).  [Verified]
- Likely reason: not apparent from the code.  [Inferred: the MIT plug-in may predate the use of ACG here]

### Private ocean-model state lives in an ESMF internal state, not in MAPL internal state
- Pattern: the plug-ins keep the model's own data structures (`ocean_public_type`, `ice_ocean_boundary_type`, the model state) in a pointer wrapper stored with
  `ESMF_UserCompSetInternalState` and restart through the model's own restart routine from the `Record` entry point
  (`MOM_GEOS5PlugMod.F90`, `Initialize`, `Record`; `MOM6_GEOSPlug.F90`, `Initialize`, `Record`).  [Verified]
- Likely reason: the model's restart format and I/O belong to the model; MAPL's internal-state restart handles only the few fields declared with
  `MAPL_AddInternalSpec` (for example `TS_FOUND`, `FRACICE`, `VOLICE`).  [Inferred]
