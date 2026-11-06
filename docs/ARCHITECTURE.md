# GEOS_OceanGridComp Architecture

**Audience:** coding agents and developers working in the `GEOS_OceanGridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOS_OceanGridComp` @ `8d91ab7` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Paths are relative to this repository's root unless they start with `src/` or `@`
(fixture-relative). Statements without a label are Verified and carry a
`(path, routine-or-key)` citation. *Inferred:* and *Unverified:* mark the weaker
cases. `OGCM file` means `src/Components/@GEOSgcm_GridComp/GEOSogcm_GridComp/GEOS_OgcmGridComp.F90`
(parent repository `GEOSgcm_GridComp`). `MAPL file` means
`src/Shared/@MAPL/generic/MAPL_Generic.F90`. Definitions at first use: MOM5 and MOM6 are the
Modular Ocean Model versions 5 and 6 (GFDL); MITgcm is the MIT General Circulation Model;
FMS is the GFDL Flexible Modeling System library; a *plug-in* here is a MAPL gridded component
that wraps one ocean model or a data ocean.

## 1. Placement

- The fixture lists this repository as `GEOS_OceanGridComp` with
  `local: ./src/Components/@GEOSgcm_GridComp/GEOSogcm_GridComp/@GEOS_OceanGridComp`,
  `tag: v3.13.0`, `develop: develop` (`components.yaml`, `GEOS_OceanGridComp`). This document was
  checked against the `develop` checkout; the tags below are the pinned versions in `components.yaml`.
- Parent repository: `GEOSgcm_GridComp`, whose ocean composite `OGCM`
  (`GEOSogcm_GridComp/GEOS_OgcmGridComp.F90`, module `GEOS_OgcmGridCompMod`) adds this
  repository's component as child `OCEAN` (`OGCM file`, `SetServices`:
  `OCEAN = MAPL_AddChild(GC, NAME='OCEAN', SS=OceanSetServices)`). Documents:
  `GEOSgcm_GridComp` architecture (`../../../docs/ARCHITECTURE.md` in the fixture) (the `OGCM` composite itself is
  marked "not yet documented" there), fixture architecture (`../../../../../../docs/ARCHITECTURE.md` in the fixture),
  application repository (`../../../../../Applications/@GEOSgcm_App/docs/ARCHITECTURE.md` in the fixture)
  (run-time templates and `gcm_setup`).
- Children (links down, all gitignored here, all listed in `components.yaml`):

| Fixture name | Local path in this repository | `tag:` / `develop:` | Used by |
|---|---|---|---|
| `mom` (remote `../MOM5.git`) | `MOM_GEOS5PlugMod/@mom` | `geos/5.1.0+2.1.0` / `geos5` | `MOM_GEOS5PlugMod` |
| `mom6` (remote `../MOM6.git`, `recurse_submodules`) | `MOM6_GEOSPlug/@mom6` | `geos/v3.9` / `geos/main` | `MOM6_GEOSPlug` |
| `stochastic_physics` | `MOM6_GEOSPlug/@stochastic_physics` | `geos/v1.0.0` / `geos/master` | `MOM6_GEOSPlug` (optional) |
| `mit` (remote `../MITgcm.git`) | `MIT_GEOS5PlugMod/@mit` | `checkpoint68o` / none | `MIT_GEOS5PlugMod` |


## 2. Directory map

| Path | What it is | Module / entry |
|---|---|---|
| `GEOS_OceanGridComp.F90` | The `OCEAN` wrapper: chooses one plug-in, runs it, converts and weights imports | `GEOS_OceanGridCompMod`, `SetServices` |
| `GEOS_Ocean_StateSpecs.rc` | Import, internal and export specs of `OCEAN` (ACG input) | `CMakeLists.txt`, `mapl_acg` |
| `GEOSdatasea_GridComp/` | Data ocean plug-in (`DATASEA`) | `GEOS_DataSeaGridCompMod` |
| `MOM_GEOS5PlugMod/` | Plug-in around MOM5 (`MOM`); `mom_app/mkMAPLTripolar` grid-file scripts | `MOM_GEOS5PlugMod` |
| `MOM6_GEOSPlug/` | Plug-in around MOM6 (`MOM6`); `mom6_app/` run-time configurations; `mom6_cmake/`, `stochastic_physics_cmake/` build the nested sources | `MOM6_GEOSPlugMod` |
| `MIT_GEOS5PlugMod/` | Plug-in around MITgcm (`MIT`); `mitgcm_setup/` driver code, `configs/` experiment setups, `instructions/` | `MIT_GEOS5PlugMod` |
| `.github/`, `.circleci/`, `README.md`, `CONTRIBUTING.md`, licence files | CI, code owners, project files | |

The build target is named `GEOS_OceanGridComp` (`esma_set_this`, `CMakeLists.txt`).

## 3. The `OCEAN` component

### 3.1 Purpose and place in the hierarchy
`OCEAN` is a thin wrapper (module header: "light-weight gridded component that serves an interface
to ocean/data_ocean components", `GEOS_OceanGridComp.F90`). It has exactly one active ocean child
(`OCN`), plus a second data-ocean child (`OCNd`) in dual-ocean mode. It is run by `OGCM` through
`MAPL_GenericRunChildren` in ordinary mode, or explicitly by phase in dual-ocean mode (`OGCM file`, `Run`).
It exchanges data with sibling components of `OGCM` (sea ice, ocean radiation, ocean biogeochemistry), not with
the atmosphere directly: the atmosphere-ocean exchange is made by `GEOS_GcmGridComp.F90` and `OGCM`
(see the `GEOSgcm_GridComp` architecture (`../../../docs/ARCHITECTURE.md` in the fixture), section 3.3).

### 3.2 `SetServices` (`GEOS_OceanGridComp.F90`)
1. Reads `USE_DATASEA:` (default 1), `OGCM_TOP_LAYER:` (default 10.0, called `OrphanDepth`) and
   `USE_DATA_ATM4OCN:` (default `.FALSE.`).
2. Adds the active child:
   - `USE_DATASEA` /= 0: child `DATASEA`, `MAPL_AddChild(GC, NAME=..., SS=DataSeaSetServices)`.
   - else reads `OCEAN_NAME:` (default `"MOM"`): `MOM` and `MOM6` are loaded from shared libraries,
     `MAPL_AddChild(OCEAN_NAME, 'setservices_', parentGC=GC, sharedObj=...)`, with the library name from
     `MOM_GEOS5PLUGMOD:` (default `libMOM_GEOS5PlugMod.so`) or `MOM6_GEOSPLUG:` (default
     `libMOM6_GEOSPlug.so`); `MIT` is added statically, `SS=MITSetServices`, only when compiled with
     `BUILD_MIT_OCEAN`; any other value writes a message and fails (`VERIFY_(999)`).
3. Reads `DUAL_OCEAN:` (default 0). If set, adds a second child named `DATASEA` (`OCNd`).
4. Declares state: `#include` of the ACG output for `GEOS_Ocean_StateSpecs.rc` (imports, one internal,
   exports), the hand-written export `MASKO`, and exports borrowed from the child
   (`MAPL_AddExportSpec(..., CHILD_ID=OCN)`): always `TW SW UW VW`; when not data ocean also
   `DH UWB VWB SLV T S`; `SSH PBO` only for `MOM`; `UWC VWC` only for `MOM6`.
5. For a non-data ocean, `MAPL_TerminateImport` removes from the child the imports `TAUX TAUY PENUVR PENPAR PENUVF
   PENPAF DRNIR DFNIR DISCHARGE CALVING LWFLX SHFLX QFLUX RAIN SNOW SFLX SWHEAT`: `OCEAN` fills them by hand
   in `Run`. `PEN_OCN` is deliberately not terminated "since it is not used in the plug" (source comment).
   With `DUAL_OCEAN`, the child's `DEL_TEMP` import is terminated too ("we will fill it here").
6. Sets Initialize and Run entry points (`Initialize`, `Run`); with `DUAL_OCEAN` a second `Run` entry point
   (phase 2, "only used in the predictor part of the replay for dual ocean"). Finalize and the record
   (checkpoint) entry point are not overridden, so the MAPL defaults apply (`MAPL file`,
   `register_generic_entry_points`: `MAPL_GenericFinalize`, `MAPL_GenericRecord`).
7. Calls `MAPL_GenericSetServices`, which wires the children to each other and to `OCEAN`'s own states
   (`MAPL file`, `MAPL_WireComponent`).

### 3.3 State
- **Imports** are the specs in `GEOS_Ocean_StateSpecs.rc` (`FROCEAN TAUX TAUY PENUVR PENPAR PENUVF PENPAF DRNIR
  DFNIR SWHEAT DISCHARGE LWFLX SHFLX QFLUX SNOW RAIN FRESH FSALT FHOCN PEN_OCN`; `CALVING` only if
  `DO_DATA_ATM4OCN`; bundles `TR` and `TRFLUX` only if `OCEAN_NAME == 'MOM'`; `FRACICEd` only if
  `dual_ocean`) **plus** any child import that is neither connected nor terminated: MAPL adds such child
  import specs to the parent's import spec (`MAPL file`, `MAPL_WireComponent`: "Imports that are not internally
  satisfied have their specs put in the GC's import spec to be externally satisfied"). *Inferred:* from that
  code, not checked in a run, this is how child-only imports such as `PS` (all real-ocean plug-ins), `PICE` (`MOM`, `MOM6`),
  `STROCNXB STROCNYB AICEU` (`MOM`), `TAUXBOT TAUYBOT AICE` (`MOM6`), `WGHT` and `ICESTATES` (`MIT`), and
  `FRACICE`, `DATA_SST`, `DATA_UW` (`DATASEA`) appear in `OCEAN`'s import state. `OGCM` then fills or
  connects them (3.5).
- **Internal:** `TS_FOUND` (K, default 280.0, `FRIENDLYTO = trim(COMP_NAME)`), so it is also an export
  (`GEOS_Ocean_StateSpecs.rc`; `MAPL file`, `MAPL_StateAddInternalSpec`).
- **Exports:** `SS_FOUND FRZMLT TAUX TAUY SWHEAT RFLUX DISCHARGE CALVING FROCEAN LWFLX SWFLX SHFLX QFLUX SFLX RAIN
  SNOW PEN_OCN` (diagnostic copies of what `Run` passed on), `T_Freeze` only if `OCEAN_NAME == 'MOM6'`
  (`GEOS_Ocean_StateSpecs.rc`), `MASKO`, and the borrowed child exports of 3.2.
- **Private state:** one clock, stored with `ESMF_UserCompSetInternalState` under the name
  `<OCEAN_NAME>_internal_state` (`GEOS_OceanGridComp.F90`, `Initialize`).
- **Restarts:** only `TS_FOUND` is declared here, with no `RESTART` option (MAPL default
  `MAPL_RestartOptional`, `MAPL file`, `MAPL_StateAddInternalSpec`). The template names its file with
  `OCEAN_INTERNAL_RESTART_FILE: ocean_internal_rst` on `@COUPLED` lines (`src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`).
  The ocean models restart themselves (section 5).

### 3.4 `Initialize` and `Run` (`GEOS_OceanGridComp.F90`)
`Initialize`:
- Creates the private clock named `<OCEAN_NAME>Clock`, step `OCEAN_DT:` (default `RUN_DT:`, the heartbeat),
  starting at the current time.
- Calls `MAPL_GenericInitialize` (which initializes the children).
- Non-data ocean only: allocates export `MASKO` and fills it from the child's mask: `<OCEAN_NAME>_3D_MASK` level 1 for
  `MOM`/`MIT`, `MOM_2D_MASK` for `MOM6`.

`Run` (phase 1), when the private clock is not past the supervisor's time (it rewinds the private clock if it is):
1. Gets the mask as above (data ocean: a field of 1.0).
2. Weight `WGHT = FROCEAN / MASK` where `MASK > 0` and `FROCEAN /= MAPL_UNDEF`, else 0.
3. Non-data ocean: copies imports into the child's import state multiplied by `WGHT`: `PENUVR PENPAR PENUVF PENPAF DRNIR
   DFNIR LWFLX QFLUX SNOW TAUX TAUY`; `SHFLX = (SHFLX - FHOCN)*WGHT`; `RAIN = (RAIN + FRESH)*WGHT`;
   `SFLX = FSALT*WGHT`; `DISCHARGE = DISCHARGE*WGHT` and `CALVING = 0`, or, with `USE_DATA_ATM4OCN`, `DISCHARGE`
   and `CALVING` taken unweighted and clipped at 0; `SWHEAT` weighted for `MOM`/`MIT` only (the
   `MOM6` branch does nothing: "No 3D Mask from MOM6. Do nothing for now!") and `RFLUX` accumulates the
   weighted heating in levels where the 3D mask is 0. With
   `BUILD_MIT_OCEAN` it also sets the child's `WGHT` import. Then fills the diagnostic exports.
4. Runs the child in a loop `do while (MyTime <= endTime)`, advancing the private clock by `OCEAN_DT`:
   ordinary mode `MAPL_GenericRunChildren(GC, IMPORT, EXPORT, PrivateState%CLOCK)`; dual mode in 3.6.
5. After the loop: `SS_FOUND` (34.0 psu where `WGHT = 0`, else child `SW`); `FRZMLT` from the child where `WGHT > 0`
   (data ocean: 0); `T_Freeze` likewise, only for `MOM6`; `TS_FOUND` from child `TW` where `WGHT > 0` (dual mode:
   from `TWd` where the ice fraction is below 0.05).
6. "Orphan points" (ocean fraction present, child mask 0): `TS_FOUND` is advanced with the surface heat budget over
   `OrphanDepth`, `FRZMLT` computed, and `TS_FOUND` limited below by the freezing temperature at salinity 34.0 (non-data
   ocean only).

### 3.5 Interface with the parent `OGCM` (`OGCM file`)
- Selection keys read in `OGCM`'s `SetServices`: `USE_DATASEA:` (default 1; comment "this get resource is
  repeated in Ocean - change both together!"), `OCEAN_NAME:` (default `MOM`), `USE_DATA_ATM4OCN:`, `DUAL_OCEAN:`,
  `OCEAN_EXT_DATA:`, `OCEAN_SSS_DATA:` and `SEAICE_THICKNESS_EXT_DATA:` (the last two only when `USE_DATASEA == 1`).
  `OCEAN`'s `SetServices` reads `USE_DATASEA:`, `OCEAN_NAME:`, `USE_DATA_ATM4OCN:`, `DUAL_OCEAN:` again, and each plug-in
  reads `DUAL_OCEAN:` once more. `OGCM` asserts: `USE_DATA_ATM4OCN` needs `USE_DATASEA = 0`; `USE_DATASEA /= 0` needs
  `USE_CICE_Thermo = 0`, `USE_DATASEAICE /= 0` and `USE_OCEANOBIOGEOCHEM = 0`.
- Connections `MAPL_AddConnectivity` in `SetServices`, to or from `OCEAN`:
  `ORAD -> OCEAN: SWHEAT` and `OCEAN -> ORAD: DH, MASKO` (non-data ocean); `SEAICE -> OCEAN: FRACICE` (always),
  `SEAICETHICKNESS` (if `SEAICE_THICKNESS_EXT_DATA`), `ICESTATES` (`BUILD_MIT_OCEAN`), `FRACICEd` (dual);
  `OCEAN -> SEAICE: SS_FOUND` (dual), `UWB VWB UW VW SLV` (non-data ocean), `UWC VWC` (non-data, `OCEAN_NAME` neither `MOM` nor `MIT`),
  `TS_FOUND SS_FOUND FRZMLT` renamed `SST SSS FRZMLT` (`USE_CICE_Thermo > 1`); `SEAICE -> OCEAN` for a non-data ocean:
  `TAUXBOT TAUYBOT HICE HSNO FRESH FSALT FHOCN` plus `STROCNXB STROCNYB AICEU` for `MOM` and `AICE` otherwise (the
  `else` branch is commented "MOM6" but also covers `MIT`); `OCEAN -> OBIO: DH T S MASKO` (with ocean biogeochemistry).
- Imports of `OCEAN` that `OGCM` fills by hand on the exchange grid from the surface skin (not by connection), with
  `MAPL_LocStreamTransform` into `GIM(OCEAN)` (`OGCM file`, `Run`): `TAUX TAUY OUSTAR3 PS` (pointers fetched with
  `notfoundOK`), `LWFLX SHFLX QFLUX SNOW RAIN PEN_OCN`, and for a real ocean `PENUVR PENUVF PENPAR PENPAF DRNIR DFNIR` and
  `DISCHARGE` (not with `USE_DATA_ATM4OCN`). `FROCEAN` is filled in `Initialize` for every child that has it
  (`MAPL_LocStreamFracArea`).
- `OGCM`'s `SetServices` ends with calls of `MAPL_TerminateImport(GC, names, child_ids)` with arrays. This is the "all but" form:
  every unconnected import is terminated *except* one whose name is in the list, for a child that is in the list (not positional pairs:
  `MAPL file`, `MAPL_TerminateImportAllBut`, `skip = any(names == short_name) .and. any(child_ids == I)`; for example `DATA_SIT`, listed for `seaice`, is also
  exempt for `ocean` when `ocean` is in the same call). The imports left unterminated for `OCEAN` are: `DATA_UW DATA_VW` when `OCEAN_EXT_DATA` is false (plus `DISCHARGE CALVING` with `USE_DATA_ATM4OCN`); `DATA_SST`
  (and `DATA_SSS` when `OCEAN_SSS_DATA`) when `OCEAN_EXT_DATA` is true and the ocean is a data ocean; none when it is true and the
  ocean is real (`OGCM file`, `SetServices`). *Inferred:* the kept imports are the ones left to ExtData; everything else that
  `OCEAN` imports and nobody connects (for example `PS`, `PICE`) is terminated at this level; `OGCM` fills some of those by hand (above).
  This agrees with `DataSea_ExtData.yaml`, which names `DATA_UW` and `DATA_VW` (5.1).
- Exports `OGCM` reads back from `GEOS_OceanGridComp`'s export state: `UW VW TS_FOUND SS_FOUND FRZMLT` and, non-data, `UWB VWB`
  (`OGCM file`, `Run`).
- Tracers: for `OCEAN_NAME == MOM`, `OGCM` collects ocean-biogeochemistry friendlies into the `TR` bundle of `OCEAN`'s import
  state (`OGCM file`, `Initialize`, `MAPL_GridCompGetFriendlies(GCS(OBIO),"OCEAN", BUNDLE)`).

### 3.6 Dual-ocean mode (`DUAL_OCEAN: 1`)
A second data ocean `DATASEA` (`OCNd`) runs beside the real ocean and the real ocean is nudged toward it, for the predictor part
of regular replay (`REPLAY_MODE: Regular`, asserted in `GEOS_GcmGridComp.F90`; see the
`GEOSgcm_GridComp` architecture (`../../../docs/ARCHITECTURE.md` in the fixture), section 3.4).
- `Run` phase 1 (corrector): run `OCNd` phase 1, its couplers, `OCN` phase 1, its couplers; then
  `DEL_TEMP = (TWd - TW)*DT/(DT + TAU_SST)` where mask > 0 and `FRACICE` of `OCNd` is below 0.05, and
  `DEL_TEMP = ((-0.054*SW + MAPL_TICE) - TW)*DT/(DT + TAU_SST_UNDER_ICE)` where mask > 0, that ice fraction is at least 0.05 and
  import `FRACICEd` exceeds it (`TAU_SST:` default 432000 s, `TAU_SST_UNDER_ICE:` default 86400 s); `DEL_TEMP` is 0 elsewhere.
  Then run `OCN` phase 2, which adds `DEL_TEMP` to the top temperature level.
- `Run` phase 2 (predictor): run `OCNd` phase 1 and its couplers only.

## 4. The plug-ins at a glance

Selection is made twice: `OGCM` reads the keys for connections, `OCEAN` for the child (3.2, 3.5).

| | `DATASEA` | `MOM` | `MOM6` | `MIT` |
|---|---|---|---|---|
| Directory, module | `GEOSdatasea_GridComp`, `GEOS_DataSeaGridCompMod` | `MOM_GEOS5PlugMod`, `MOM_GEOS5PlugMod` | `MOM6_GEOSPlug`, `MOM6_GEOSPlugMod` | `MIT_GEOS5PlugMod`, `MIT_GEOS5PlugMod` |
| Wraps | Nothing: reads SST (and optionally SSS) data | MOM5 (`@mom`) with FMS | MOM6 (`@mom6`) with FMS, optional `@stochastic_physics` | MITgcm (`@mit`) through generated driver code |
| Selected by | `USE_DATASEA:` /= 0 (default 1) | `USE_DATASEA: 0` and `OCEAN_NAME: MOM` (default) | `USE_DATASEA: 0` and `OCEAN_NAME: MOM6` | `USE_DATASEA: 0`, `OCEAN_NAME: MIT`, and build option `BUILD_MIT_OCEAN` |
| Added as child | static `SS=DataSeaSetServices` | shared object `libMOM_GEOS5PlugMod.so`, key `MOM_GEOS5PLUGMOD:`, symbol `setservices_` | shared object `libMOM6_GEOSPlug.so`, key `MOM6_GEOSPLUG:`, symbol `setservices_` | static `SS=MITSetServices` |
| State specs | `GEOS_DataSea_StateSpecs.rc` (ACG) | `MOM_GEOS5PlugMod_StateSpecs.rc` (ACG, import and export only) | `MOM6_GEOSPlug_StateSpecs.rc` (ACG) | hand-written `MAPL_Add*Spec` in `SetServices` |
| Restarts | None of its own | Model files via `ocean_model_restart` (Record), `ocean_model_end` (Finalize) | same pattern | Internal `FRACICE`, `VOLICE`; Record and Finalize do no model I/O |

Entry points, imports, exports and build details are in section 5 and 8; the `*_StateSpecs.rc` files (and, for `MIT`, `SetServices`) are the source of truth for field lists.

## 5. Plug-in subsections

### 5.1 `GEOSdatasea_GridComp` (`DATASEA`)
- Purpose: stands in for an ocean by passing sea surface temperature and, optionally, salinity and surface currents to the rest of the
  coupled system (`GEOS_DataSeaGridComp.F90`, `RUN` description: "Periodically refreshes the SST, sea ice and optionally SSS
  from dataset(s)"). It is the default ocean (`USE_DATASEA` defaults to 1), used for atmosphere-only (AMIP-style) runs and, in
  dual-ocean mode, as the second child.
- `SetServices`: sets only the `Run` entry point; reads `OCEAN_EXT_DATA:` and `OCEAN_SSS_DATA:` (both default `.FALSE.`) which are
  conditions of the ACG imports `DATA_SST` and `DATA_SSS` (`GEOS_DataSea_StateSpecs.rc`, `COND`); calls `MAPL_GenericSetServices`.
- `Run`: obtains `SST` from the ExtData import `DATA_SST` if `OCEAN_EXT_DATA`, else from the binary file named by `DATA_SST_FILE:`
  with `MAPL_ReadForcing` (`OGCM_IS_FCST:` selects `INIT_ONLY`); `SSS` from `DATA_SSS` if `OCEAN_SSS_DATA`, else `SW = 30.0`.
  `UW`, `VW` come from `DATA_UW`, `DATA_VW`, set to 0 where the speed is 10 m s-1 or more. `SST_ADJ_UND_ICE:` (default 0):
  1 limits SST to the freezing value and blends it with `FRACICE`; 2 relaxes toward it with `CTB:` (default 1.0e-4) and `DT:`.
  Exports `TW = SST`, `SW`.
- State: no internal state, nothing checkpointed (`GEOS_DataSea_StateSpecs.rc`, empty `category: INTERNAL`).
- ExtData: `DataSea_ExtData.yaml` points `DATA_UW` and `DATA_VW` at collection `/dev/null`; it is installed to `etc`
  (`CMakeLists.txt`) and listed in `src/Applications/@GEOSgcm_App/construct_extdata_yaml_list.py`. *Unverified:* where the
  ExtData entries for `DATA_SST` and `DATA_SSS` come from; none is in this repository or in `AGCM.rc.tmpl`/`gcm_setup`
  (grep of `OCEAN_EXT_DATA`, `DATA_SST`).
- Coupling: `FRACICE` import comes from `SEAICE` via `OGCM`'s connection. `OGCM` leaves `DATA_UW DATA_VW` unterminated when `OCEAN_EXT_DATA` is
  false, and `DATA_SST` (with `DATA_SSS` when `OCEAN_SSS_DATA`) when it is true; all other unconnected imports are terminated (3.5).

### 5.2 `MOM_GEOS5PlugMod` (`MOM`, MOM5)
- Purpose: wraps MOM5 (`@mom`). Header: "MAPL/ESMF Gridded Component that acts as a wrapper for MOM. ... heavy dependencies on FMS and
  MOM" (`MOM_GEOS5PlugMod.F90`). `gcm_setup` prints a warning that a MOM5 coupled set-up "is _known_ to have problems" and refers to
  an issue in the MOM5 repository (`src/Applications/@GEOSgcm_App/gcm_setup`, `OCNMODEL` branch `MOM5`); *Unverified:* the
  content of that issue.
- Selection and loading: section 4. A free `subroutine SetServices` at the end of the file calls the module's `SetServices` (*Inferred:* so the shared
  object has the symbol `setservices_` that `MAPL_AddChild` looks up).
- `SetServices`: reads `DUAL_OCEAN:`; ACG import and export specs; entry points `Initialize`, `Run`, `Finalize`, write-restart
  `Record`, and `Run2` as a second run phase when dual (`SetServices`).
- `Initialize`: takes the FMS layout from `MAPL_Get(NX, NY)`, the date from the clock, `fms_init` with the VM's MPI
  communicator, `ocean_model_init`, checks MOM's local size against the MAPL grid (assertions), allocates the flux "bulletin board"
  `Ice_ocean_boundary`, calls `MAPL_GenericInitialize`, then fills `TW SW DH MOM_3D_MASK AREA` and the other
  exports early because the parent reads them before the first `Run`.
- `Run`: converts to MOM units and sign: `salt_flux = -SFLX`, `lw_flux`, `t_flux`, `q_flux`, `lprec`, `fprec`, `runoff`, shortwave
  components (`PENUVR+PENPAR` as visible direct, ...), pressure `P = MOM_PICE_SCALING * PICE` (key default 1.0); builds the
  shortwave decay fraction `swheat_fr` from `SWHEAT` and surface shortwave, subtracts the surface flux from level 1 of `SWHEAT`; converts the
  A-grid stress to B-grid (`transformA2B`), rotates it with `cos_rot`/`sin_rot` and combines it with ice-base stress:
  `U_flux = (...)*(1 - AICEU) - STROCNXB`; copies the `TR` bundle into MOM tracers named `GENtracer_NNN`, calls
  `update_ocean_model` (skipped when `steady_state_ocean:` /= 0), copies tracers back, and fills exports (B-grid and A-grid
  surface currents, `TW`, `SW`, `SSH`, `SLV`, `FRZMLT`, `DH`, `T`, `S`, `RHO`, `HC`, `MLD`, ...). Header notes (*Unverified:*) that `MLD`
  is the KPP boundary layer depth and `PSI` is not filled.
- `Run2` (dual only): adds `DEL_TEMP` to the top level of MOM's temperature tracer (`mom4_put_prog_tracer`).
- `Finalize`: `ocean_model_end`, `diag_manager_end`, `field_manager_end`, `fms_io_exit` (`fms_end` is commented out), then `MAPL_GenericFinalize`.
- Restarts: `Record` calls `ocean_model_restart(Ocean_State, timeStamp)` when `MAPL_RecordAlarmIsRinging(MODE=MAPL_Write2Disk)` is
  true; `ocean_model_end` is also called in `Finalize`. No MAPL internal state. `gcm_run.j` copies `$EXPDIR/RESTART/*` to `INPUT`
  before the run and `RESTART/*` to `INPUT` after each segment (`@COUPLED` lines). *Inferred:* MOM5 reads its restart from `INPUT`
  during `ocean_model_init`.
- Run-time model configuration: `gcm_setup` copies `input.nml` and the `*table` files from `$ETCDIR/MOM5/geos5/<IM>x<JM>/INPUT`;
  `@mom/CMakeLists.txt` installs `data/geos5` to `etc/MOM5`. The grid file `MAPL_Tripolar.nc` (`OGCM.GRIDSPEC` in `AGCM.rc.tmpl`) is
  made with `mom_app/mkMAPLTripolar/mkMAPLTripolarMOM5.py`.
- Trap: `TR` is a bundle of references; `Run` copies its arrays into MOM tracers before `update_ocean_model` and back afterwards.

### 5.3 `MOM6_GEOSPlug` (`MOM6`, MOM6)
- Purpose: wraps MOM6 (`@mom6`) through its FMS driver (`ocean_model_MOM.F90`, `FMS_cap`); header: "acts as a coupler for MOM"
  (`MOM6_GEOSPlug.F90`). The CI coupled run uses `gcm_ocean_type: MOM6` (`.circleci/config.yml`).
- Selection and loading: section 4; the free `SetServices` at the end of the file (it `use`s `MOM6_GEOSPlugMod`).
- `SetServices`: reads `DUAL_OCEAN:`; ACG imports, exports (`MOM6_GEOSPlug_StateSpecs.rc`); entry points as for `MOM`.
- `Initialize`: `MOM_infra_init`, `field_manager_init`, `MOM_diag_manager_init`; checks `ocean_wind_stagger:` (default `AGRID`;
  `BGRID`/`CGRID` stop with a message, other values are errors); `ocean_model_init(..., iwind_stagger)`,
  `ocean_model_init_sfc`; saves `Ocean_State` in the private state type; checks the local domain against the MAPL grid;
  allocates `Boundary` over the data domain (`g_isd:g_ied, g_jsd:g_jed`); `MAPL_GenericInitialize`; fills `MOM_2D_MASK`, `TW`, `SW`, `AREA`, `SLV`,
  `DH`, `T`, `S`.
- `Run`: same conversion pattern as `MOM`, plus `calving`; stress on the water `U = TAUX*(1-AICE) - TAUXBOT*AICE` (A-grid), rotated
  with `cos_rot`, `sin_rot`; boundary pressure from `MOM_PICE_SCALING:` (default 1.0) times `PICE`, plus `PS` when `pres_loading:` is `YES`
  (default `NO`); `DTFREEZE_DS:` (default -0.054); `update_ocean_model` unless `steady_state_ocean:` /= 0 (the other branch calls it with
  two extra `.false.` arguments; source comment: "though this will cause MOM_FATAL error"). Exports `FRAZIL` and `MELT_POT` from
  heat fluxes, `FRZMLT` = `MELT_POT` or, where `FRAZIL > 0`, `FRAZIL`; `T_Freeze = dTf_dS*SW` (where mask > 0, else -1.8); `UWC`/`VWC`.
- `Run2` (dual only): adds `DEL_TEMP` to level 1 of the `temp` prognostic tracer.
- `Finalize`: `ocean_model_end` ("this also calls ocean_model_save_restart", source comment), `MOM_diag_manager_end`, `field_manager_end`,
  `io_infra_end`, deallocates `Boundary`. `Record`: `ocean_model_restart(Ocean_State, timeStamp)` on the same alarm as `MOM`.
- Restarts and configuration: MOM6 writes `RESTART/`; `gcm_run.j` copies `MOM_input` and `MOM_override` into the run directory
  and sets `DT` and `DT_THERM` in `MOM_override` to the heartbeat. `gcm_setup` copies `MOM_input`, `MOM_override`, `input.nml`,
  `*table` from `$ETCDIR/MOM6/mom6_app/<IM>x<JM>/`, plus `MOM_saltrestore` and `common/MOM_dataatm` if a data atmosphere is
  used (`CMakeLists.txt` installs `mom6_app` to `etc/MOM6`). Resolutions present: `72x36`, `540x458`, `720x576`, `1440x1080`,
  `2880x2240` (`mom6_app/`).
- Stochastic physics: `@stochastic_physics` is built into the MOM6 library when found, else MOM6's compile-only stubs are used with a warning
  (`CMakeLists.txt`, `MOM6_USE_STOCHASTIC_PHYSICS`); the plug-in source does not mention it.

### 5.4 `MIT_GEOS5PlugMod` (`MIT`, MITgcm)
- Purpose: wraps MITgcm (`@mit`) on an `llc` ocean grid (`OGCM_GRID_TYPE = llc` for `MIT` in `gcm_setup`; `OGCM.GRIDSPEC: mit.ascii` in `AGCM.rc.tmpl`). Header: "wrapper for MITgcm" (`MIT_GEOS5PlugMod.F90`). Compiled only with `BUILD_MIT_OCEAN`, and then
  `MOM` and `MOM6` are not built (`CMakeLists.txt`).
- Link to sea ice: with `BUILD_MIT_OCEAN` the sea-ice composite uses `MITSEAICEDYNA` (`GEOSseaice_GridComp/GEOS_SeaIceGridComp.F90`,
  `GEOSMITDyna_GridComp`), which exports the state container `ICESTATES` (connected `SEAICE -> OCEAN`, 3.5). `MIT_GEOS5PlugMod`'s `Initialize`
  reads the attribute `ICECOMPNAME` from it and gets `<name>_Imports` and `<name>_Exports`.
- `SetServices`: reads `DUAL_OCEAN:`; declares by hand 19 imports and 17 exports with a table of `mstate` entries; every import is also
  mirrored as an export "for diagnostic purposes"; adds `ICESTATES` (state item), `FRACICE`/`VOLICE` exports and internal specs
  (`MAPL_RestartOptional`, default 0.0); reads `CICE_N_ICE_CATEGORIES:` and `CICE_N_ICE_LAYERS:` unconditionally with
  `ESMF_ConfigGetAttribute`; entry points as for `MOM`.
- `Initialize`: reads `OCEAN_DIR:` (template: `mitocean_run`, `@MIT` line), changes into it (`mysetdir`), `DRIVER_INIT` with the MPI
  communicator, back (`popdir`), `MAPL_GenericInitialize`; sets `DH` from hard-wired values (50 m, 70 m, ... 690 m for level 15 and below) and
  fills `<name>_3D_MASK`, `TW`, `SW` from the driver.
- `Run`: reads `RUN_DT:`, `OCEAN_DT:`, `OCEAN_DIR:`, `ACTIVE_OCEAN:` (default 1; with 0 the driver is not run); passes imports and the ice state to the driver with
  `DRIVER_SET_IMPORT_STATE`, runs `DRIVER_RUN(istate, 1)`, reads exports with `DRIVER_GET_EXPORT_STATE` (`US VS TS SS MASK`, ice
  tendencies `DEL*`, `FRZMLT`). The weight `WGHT` is provided by `OCEAN` (`#ifdef BUILD_MIT_OCEAN`, 3.4).
- `Run2` (dual): adds `DEL_TEMP` to `TW` where the mask is positive and passes it with `DRIVER_SET_IMPORT_STATE(..., 'TS', T)` then
  `mitgcm_import_fill_ts`.
- `Record` evaluates the record alarm and does nothing else; `Finalize` only calls `MAPL_GenericFinalize` (`MIT_GEOS5PlugMod.F90`).
  *Unverified:* MITgcm's own pickup files are controlled by its namelists (`configs/<id>/input/data`, `pickupStrictlyMatch`).
- Build and setup: `MIT_CONFIG_ID` selects `configs/<id>/code` (c12_cs32_01, c90_llc90_02/03/05, c180_llc270_01); `genmake2` and custom
  targets `mitgcmuv`, `mitgcmuv_driver` build the model and the driver library in `mitgcm_setup/` (`CMakeLists.txt`);
  checkout recipes in `instructions/`.

## 6. Nested model repositories (link-only)

`@mom`, `@mom6`, `@stochastic_physics` and `@mit` are separate repositories (section 1); this repository holds only their wrappers
and configuration. Their own documents (present in a fixture checkout): `MOM6_GEOSPlug/@mom6/AGENTS.md`, `MIT_GEOS5PlugMod/@mit/README.md`,
`MOM_GEOS5PlugMod/@mom/README.md`, `MOM6_GEOSPlug/@stochastic_physics/README.md`. Change them only in their own repository.

## 7. Configuration that changes structure

| Key | File that sets it | Read in | Effect |
|---|---|---|---|
| `USE_DATASEA:` (default 1) | `AGCM.rc` (`@COUPLED USE_DATASEA:0` in `AGCM.rc.tmpl`) | `GEOS_OceanGridComp.F90`, `GEOS_OgcmGridComp.F90` | data ocean or real ocean; read twice |
| `OCEAN_NAME:` (default `MOM`) | `AGCM.rc` (`@OCEAN_NAME`, set in `gcm_setup`) | same two files | `MOM`, `MOM6` or `MIT` |
| `MOM_GEOS5PLUGMOD:`, `MOM6_GEOSPLUG:` | none found in the fixture (defaults used) | `GEOS_OceanGridComp.F90` | shared-library names |
| `DUAL_OCEAN:` (default 0) | commented hint `@MOM5#DUAL_OCEAN: 1` in `AGCM.rc.tmpl` | `OCEAN`, `OGCM`, each plug-in | second `DATASEA` child, phase 2 |
| `USE_DATA_ATM4OCN:` (default false) | `AGCM.rc` | `OCEAN`, `OGCM` | `CALVING` import, unweighted discharge |
| `OCEAN_EXT_DATA:`, `OCEAN_SSS_DATA:` | none found in `AGCM.rc.tmpl`, `gcm_setup` | `DATASEA`, `OGCM`, ocean radiation, data sea ice | data source (ExtData or file) and which imports are left unterminated |
| `OCEAN_DT:` (default `RUN_DT:`) | not in `AGCM.rc.tmpl` (it sets `OGCM_RUN_DT:`) | `OCEAN`, `MIT` | private clock step |
| `OGCM_TOP_LAYER:` (default 10.0), `TAU_SST:`, `TAU_SST_UNDER_ICE:` | none found | `GEOS_OceanGridComp.F90` | orphan-point depth, dual-ocean relaxation |
| `ocean_wind_stagger:` (default `AGRID`), `pres_loading:` (default `NO`), `DTFREEZE_DS:`, `MOM_PICE_SCALING:` (default 1.0), `steady_state_ocean:` (default 0) | only `steady_state_ocean: 0` and `OCEAN_PICE_SCALING: 0.0` are in `AGCM.rc.tmpl` (`@COUPLED` lines) | plug-ins (5.2, 5.3) | scheme options |

Resource lookup tries the component name as a prefix: `MAPL_GetResource` on a component passes its child name (`MOM`, `MOM6`, `MIT`,
`OCEAN`, `DATASEA`) and tries `<name>_<label>`, `<type>_<label>`, the bare label, then `<name>.<label>`
(`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_GetResource` with `component_name = state%compname`;
`src/Shared/@MAPL/base/MAPL_Resource.F90`, `get_labels_with_prefix`).

## 8. Build structure

- `CMakeLists.txt`: option `BUILD_MIT_OCEAN` (default OFF). OFF: library `GEOS_OceanGridComp` from `GEOS_OceanGridComp.F90` with
  `SUBCOMPONENTS MOM6_GEOSPlug MOM_GEOS5PlugMod GEOSdatasea_GridComp`. ON: defines `BUILD_MIT_OCEAN`, subcomponents
  `MIT_GEOS5PlugMod GEOSdatasea_GridComp` only. `mapl_acg` runs on `GEOS_Ocean_StateSpecs.rc`.
- `MOM_GEOS5PlugMod`: `TYPE SHARED`, `SUBCOMPONENTS @mom`, compiled with the real-8 flags (`FREAL8`); source list plus an optional
  `seawater.F90` if present; `mapl_acg` for imports and exports only.
- `MOM6_GEOSPlug`: `TYPE SHARED`, depends on `MAPL`, `FMS::fms_r8` and the target `mom6`, which `mom6_cmake/CMakeLists.txt` builds from an explicit
  list of `@mom6` sources (infrastructure `MOM6_INFRA` is `FMS1` or `FMS2`, default `FMS2`). Real-8 flags; installs `mom6_app` to `etc/MOM6`.
- `MIT_GEOS5PlugMod`: returns early with a message if `MIT_CONFIG_ID` is not defined; depends on custom targets `mitgcmuv` and `mitgcmuv_driver`;
  build options by compiler (`build_options.gnu`, `build_options.intel`).
- The shared-object plug-ins are also listed as `SUBCOMPONENTS`; whether the executable therefore links them is *Unverified:* in the
  fixture architecture (`../../../../../../docs/ARCHITECTURE.md` in the fixture), section 3 ("Plug-ins loaded by name").
- Tests and CI: no unit-test target in this repository's own `CMakeLists.txt` files. `.github/workflows/workflow.yml`: fixture builds
  (`ifort`, `gfortran-15`, `ifx`, Debug) and a spack build; `.circleci/config.yml`: fixture build, atmosphere-only run, coupled run with `gcm_ocean_type: MOM6`.

## 9. Known traps

- **Selection is read in several places.** `USE_DATASEA`, `OCEAN_NAME`, `USE_DATA_ATM4OCN` and `DUAL_OCEAN` are read by `OGCM`, by `OCEAN`
  and (for `DUAL_OCEAN`) by each plug-in. Defaults agree today (`GEOS_OgcmGridComp.F90`, `GEOS_OceanGridComp.F90`); the source comment
  asks that they change together.
- **`OCEAN`'s imports are not only its `.rc` file.** Unconnected child imports are promoted to it (3.3); adding or removing a plug-in
  import changes what `OGCM` must provide, and `OGCM`'s connections to names a child does not have are relaxed rather than failing
  (`src/Shared/@MAPL/generic/VarConn.F90`, `checkReq`).
- **Module state.** `OCEAN_NAME`, `DO_DATASEA`, `OCN`, `OCNd`, `DUAL_OCEAN`, `OrphanDepth` are module variables (also `DUAL_OCEAN` in each plug-in).
  *Inferred:* one `OCEAN` instance per executable.
- **Weighting.** Fluxes handed to the real ocean are multiplied by `FROCEAN/MASK` (3.4); exports are written back only where `WGHT > 0`, so
  other points keep their earlier or default values.
- **Mask names differ.** `MOM_3D_MASK`, `MIT_3D_MASK` (built from the component name) and `MOM_2D_MASK` (`MOM6`); `OCEAN` selects by `OCEAN_NAME`.
- **Precision.** The `MOM`/`MOM6` plug-ins are compiled with the real-8 flags (`FREAL8`, their `CMakeLists.txt`) and declare GEOS-side fields as
  `real(kind=4)` (`G5KIND`, `GeosKind`); conversions are explicit in `Run`.
- **Clocks.** Ordinary mode runs the child with `OCEAN`'s private clock (step `OCEAN_DT`); dual-ocean mode with the supervisor `CLOCK`.
- **Phase numbers.** `Run` reads `currentPhase` and subtracts 10 if it is 10 or more ("to be replaced by MAPL get_phase").

## 10. Decisions, observations and questions
Decisions: [`decisions.md`](decisions.md).
### Observations (noted, not fixed)
1. `OCEAN`'s `Run` fetches `CALVING` from the child's import state for every non-data ocean (`MAPL_GetPointer(GIM(OCN), CALVING, 'CALVING')`
   without `notFoundOK`); `MOM6_GEOSPlug_StateSpecs.rc` declares `CALVING`, but `MOM_GEOS5PlugMod_StateSpecs.rc` and `MIT_GEOS5PlugMod.F90` do not (grep).
   What happens for `OCEAN_NAME` `MOM` or `MIT`? (The CI coupled run uses `MOM6` only: `.circleci/config.yml`.)
2. The orphan-point block assigns `FRZMLTe` without an `associated` test, while the earlier `FRZMLT` block tests it
   (`GEOS_OceanGridComp.F90`, `Run`). Is the `FRZMLT` export always allocated at that point?
3. `MOM_GEOS5PlugMod.F90` declares `Ocean_State` as a local pointer in `Initialize`, `Run`, `Finalize` and `Record`, without `=> null()`. Only `Initialize` passes it, uninitialized, to `ocean_model_init` (which may set it); the other routines have no assignment
   to it. The private state type holds only `Ocean` and `Ice_ocean_boundary`. `MOM6_GEOSPlug.F90` stores `Ocean_State` in its private state. MOM5's `ocean_model_init`
   tests `associated(Ocean_state)` and otherwise allocates (`@mom/src/mom5/ocean_core/ocean_model.F90`). What is the intended lifetime of `Ocean_State` in the MOM5 plug-in?
4. `OCEAN_PICE_SCALING: 0.0` is set in `AGCM.rc.tmpl`; the plug-ins read `MOM_PICE_SCALING:` with default 1.0 and no source in the fixture reads
   `OCEAN_PICE_SCALING` (grep; `src/Applications/@UMD_Etc/UMD_rc/AGCM.rc` sets `MOM_PICE_SCALING`). No code in `GEOSogcm_GridComp` fills the `PICE` import (grep of
   `GEOS_OgcmGridComp.F90`, `GEOS_OceanGridComp.F90`, the sea-ice composite); the `MOM6` source comment says "PICE that is available here is all = 0".
   *Unverified:* what provides `PICE`. Which key is meant to scale ice pressure?
5. `MOM`/`MOM6`/`MIT` plug-in module headers say "It does not use the configuration" although each reads resources (`DUAL_OCEAN:`, ...).
6. `OCEAN` `Run`: the `select case (trim(OCEAN_NAME))` for the heating copy lists `"DATASEA"`, but sits inside the `DO_DATASEA == 0` branch where
   `OCEAN_NAME` is never `DATASEA`. For `MOM6` the `SWHEAT` copy is skipped, while `HEATe = HEAT` still runs when the `SWHEAT` export is requested; the
   `MOM6` plug-in's `SWHEAT` import is declared as `DUMMY_IMPORT_solar_heating_rate` (`MOM6_GEOSPlug_StateSpecs.rc`). What is the role of `SWHEAT` for `MOM6`?
7. `MIT_GEOS5PlugMod.F90`, `Run2`: `T` points to the export `TW` (`MAPL_GetPointer(EXPORT, T, 'TW')`) and `deallocate(T)` follows. What owns that memory?
8. `MIT_GEOS5PlugMod.F90`: `Initialize` overwrites `DH` with hard-wired thicknesses; `Record` and `Finalize` do no model I/O.
9. `MIT_GEOS5PlugMod/CMakeLists.txt` has `elseif()` with no condition after `if(CMAKE_BUILD_TYPE MATCHES "Debug")`, and sets `FC` to `ifort`/`icc` for any compiler id that
   matches `Intel`.
10. `OCEAN` `Run`: `T_Freeze_e` is looked up only when `OCEAN_NAME == "MOM6"` and nullified otherwise (with the data ocean `OCEAN_NAME` is `DATASEA`), so the
    branch `T_Freeze_e = -1.8` for `DO_DATASEA /= 0` appears unreachable. What is the intended `T_Freeze` export for the data ocean?
11. `OGCM file`, `SetServices`: the connection list in the `else` branch is commented `! MOM6` but also applies to `MIT`. In the `OCEAN_EXT_DATA` block (data
    ocean) the call labelled "both data_sss and data_sit" follows the two inner `if` blocks and is not inside an `else`, so it also runs when only one of the two is set;
    *Inferred:* each "all but" call terminates every import outside its own list, so the net set left unterminated is the intersection of the
    lists and the unconditional call changes nothing in the two single-flag cases.

### Questions for maintainers
- Which of `MOM` (a warning about known problems is printed by `gcm_setup`) and `MIT` are supported configurations today, and is either exercised by a test?
- Where are `OCEAN_EXT_DATA` and the ExtData entries for `DATA_SST` and `DATA_SSS` meant to be set?
