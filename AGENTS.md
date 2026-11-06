# AGENTS.md

**Audience:** coding agents and developers working in the `GEOS_OceanGridComp` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOS_OceanGridComp` @ `253669b` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

> **Build status on this line:** not built, not ported. On
> `release/MAPL-v3`, `src/Components/@GEOSgcm_GridComp/CMakeLists.txt` lists
> `GEOSogcm_GridComp` only as a commented-out entry ("not yet ported to MAPL3 -
> restore when ported"); the ocean composite `OGCM` (`GEOSogcm_GridComp`) that
> adds this repository's `OCEAN` component is therefore not compiled, and
> `GEOS_OceanGridComp` is not built. The code in this repository is the MAPL2
> code, kept for a later port. Relative to the MAPL2 line, the only source change is the
> `use MAPL_Constants` import in the `MOM` and `MOM6` plug-ins; this line is also
> behind `develop` (for example no optional stochastic physics in the MOM6 build
> and no `mom6_app/common/MOM_dataatm`), and the documents describe this line.
> The statements below about MAPL internals (`MAPL file`) and the parent `OGCM`
> describe the MAPL2 code the plug-ins were written for; they were not re-checked
> against MAPL3.

## What this is

The `OCEAN` gridded component of the GEOS coupled model and the four ocean
plug-ins it can run: a data ocean (`GEOSdatasea_GridComp`) and wrappers around
three ocean models, MOM5 (`MOM_GEOS5PlugMod`), MOM6 (`MOM6_GEOSPlug`) and
MITgcm (`MIT_GEOS5PlugMod`). `OCEAN` is a child of the ocean composite `OGCM`
(`GEOSogcm_GridComp/GEOS_OgcmGridComp.F90` in the parent repository
`GEOSgcm_GridComp`). In the GEOSgcm fixture the repository is at
`src/Components/@GEOSgcm_GridComp/GEOSogcm_GridComp/@GEOS_OceanGridComp`
(`components.yaml`, key `GEOS_OceanGridComp`).

The ocean models themselves are separate nested repositories (`@mom`, `@mom6`,
`@stochastic_physics`, `@mit`), gitignored here. `components.yaml` of the GEOSgcm fixture on this line lists `mom`, `mom6` and `mit`; it has no `stochastic_physics` entry.
They are link-only: this repository documents only the plug-ins around them.

## Documents

Read `docs/ARCHITECTURE.md` first; its section 4 is the plug-in table and
section 5 has one subsection per plug-in. There are no separate component files.

| Document | Content |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Placement, directory map, the `OCEAN` wrapper (selection, run loop, dual ocean), plug-in table and subsections, interface to the `OGCM` parent, nested model repositories, build, configuration keys, traps, observations |
| [`docs/decisions.md`](docs/decisions.md) | Previous architectural decisions (pattern, likely reason) |
| `../../docs/ARCHITECTURE.md` | Parent repository `GEOSgcm_GridComp` (root component tree; `GEOSogcm_GridComp` itself is not yet documented there) |
| `../../../../../docs/ARCHITECTURE.md` | The GEOSgcm fixture (repository tree, build, `components.yaml`) |

## Rules

- Note bugs you see; do not fix them unless asked. Record them under
  "Observations" in `docs/ARCHITECTURE.md`.
- Keep the documents current: if you change a child selection, a resource key,
  a plug-in import or export, an entry point or a build option described there,
  update the document in the same change. Cite files and routines, not line numbers.
- Do not edit generated files (`*_Import___.h`, `*_Export___.h`, ...); edit the
  `*_StateSpecs.rc` file. `MIT_GEOS5PlugMod` is the exception: it declares its
  state by hand in `SetServices`.
- A nested model repository is a separate Git repository; change it only there.
  `git status` here does not show it.
- Code owners (`.github/CODEOWNERS`): `@GEOS-ESM/ocean-team` owns everything;
  `CMakeLists.txt` files also notify `@GEOS-ESM/cmake-team`.

## Build and test

This repository does not build alone and is not built on this line (build status above); the GEOSgcm fixture (its `README.md`, `parallel_build.csh`) does not compile it. The option `BUILD_MIT_OCEAN`
(default OFF) selects the plug-ins (`docs/ARCHITECTURE.md` section 8). No unit-test target was found outside the nested model repositories.
CI: `.github/workflows/workflow.yml`, `.circleci/config.yml` (includes a coupled run with `gcm_ocean_type: MOM6`).
