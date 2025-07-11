esma_set_this (OVERRIDE mom6)

# This is for selecting the MOM6 infrastructure interface
## We default to FMS2 and look if something is passed in...
set (DEFAULT_MOM6_INFRA "FMS2")
if (NOT MOM6_INFRA)
  set (MOM6_INFRA ${DEFAULT_MOM6_INFRA})
endif ()

## ... and make sure it is allowed
set (ALLOWED_MOM6_INFRA "FMS1" "FMS2")
if (NOT MOM6_INFRA IN_LIST ALLOWED_MOM6_INFRA)
  message (FATAL_ERROR "The only allowed values for MOM6_INFRA are ${ALLOWED_MOM6_INFRA}")
endif ()

# got below list from following two commands:
# 1. find src/ -print | grep .F90 | sort -d
# 2. ls -1 pkg/CVMix-src/src/*0 ; ls -1 pkg/CVMix-src/src/drivers/*0 ; ls -1 pkg/CVMix-src/src/shared/*0 ; ls -1 pkg/GSW-Fortran/modules/*0 ; ls -1 pkg/GSW-Fortran/toolbox/*0
#
list( APPEND MOM6_SRCS
   src/ALE/coord_adapt.F90
   src/ALE/coord_hycom.F90
   src/ALE/coord_rho.F90
   src/ALE/coord_sigma.F90
   src/ALE/coord_zlike.F90
   src/ALE/MOM_ALE.F90
   src/ALE/MOM_hybgen_regrid.F90
   src/ALE/MOM_hybgen_remap.F90
   src/ALE/MOM_hybgen_unmix.F90
   src/ALE/MOM_regridding.F90
   src/ALE/MOM_remapping.F90
   src/ALE/P1M_functions.F90
   src/ALE/P3M_functions.F90
   src/ALE/PCM_functions.F90
   src/ALE/PLM_functions.F90
   src/ALE/polynomial_functions.F90
   src/ALE/PPM_functions.F90
   src/ALE/PQM_functions.F90
   src/ALE/Recon1d_EMPLM_CWK.F90
   src/ALE/Recon1d_EMPLM_WA.F90
   src/ALE/Recon1d_EMPLM_WA_poly.F90
   src/ALE/Recon1d_EPPM_CWK.F90
   src/ALE/Recon1d_MPLM_CWK.F90
   src/ALE/Recon1d_MPLM_WA.F90
   src/ALE/Recon1d_MPLM_WA_poly.F90
   src/ALE/Recon1d_PCM.F90
   src/ALE/Recon1d_PLM_CW.F90
   src/ALE/Recon1d_PLM_CWK.F90
   src/ALE/Recon1d_PLM_hybgen.F90
   src/ALE/Recon1d_PPM_CW.F90
   src/ALE/Recon1d_PPM_CWK.F90
   src/ALE/Recon1d_PPM_H4_2018.F90
   src/ALE/Recon1d_PPM_H4_2019.F90
   src/ALE/Recon1d_PPM_hybgen.F90
   src/ALE/Recon1d_type.F90
   src/ALE/regrid_consts.F90
   src/ALE/regrid_edge_values.F90
   src/ALE/regrid_interp.F90
   src/ALE/regrid_solvers.F90
   src/core/MOM_barotropic.F90
   src/core/MOM_boundary_update.F90
   src/core/MOM_check_scaling.F90
   src/core/MOM_checksum_packages.F90
   src/core/MOM_continuity.F90
   src/core/MOM_continuity_PPM.F90
   src/core/MOM_CoriolisAdv.F90
   src/core/MOM_density_integrals.F90
   src/core/MOM_dynamics_split_RK2.F90
   src/core/MOM_dynamics_split_RK2b.F90
   src/core/MOM_dynamics_unsplit.F90
   src/core/MOM_dynamics_unsplit_RK2.F90
   src/core/MOM.F90
   src/core/MOM_forcing_type.F90
   src/core/MOM_grid.F90
   src/core/MOM_interface_heights.F90
   src/core/MOM_isopycnal_slopes.F90
   src/core/MOM_open_boundary.F90
   src/core/MOM_PressureForce.F90
   src/core/MOM_PressureForce_FV.F90
   src/core/MOM_PressureForce_Montgomery.F90
   src/core/MOM_stoch_eos.F90
   src/core/MOM_transcribe_grid.F90
   src/core/MOM_unit_tests.F90
   src/core/MOM_variables.F90
   src/core/MOM_verticalGrid.F90
   src/core/MOM_porous_barriers.F90
   src/diagnostics/MOM_debugging.F90
   src/diagnostics/MOM_diagnose_KdWork.F90
   src/diagnostics/MOM_diagnose_MLD.F90
   src/diagnostics/MOM_diagnostics.F90
   src/diagnostics/MOM_harmonic_analysis.F90
   src/diagnostics/MOM_obsolete_diagnostics.F90
   src/diagnostics/MOM_obsolete_params.F90
   src/diagnostics/MOM_PointAccel.F90
   src/diagnostics/MOM_spatial_means.F90
   src/diagnostics/MOM_sum_output.F90
   src/diagnostics/MOM_wave_speed.F90
   src/equation_of_state/MOM_EOS.F90
   src/equation_of_state/MOM_EOS_base_type.F90
   src/equation_of_state/MOM_EOS_Jackett06.F90
   src/equation_of_state/MOM_EOS_Roquet_SpV.F90
   src/equation_of_state/MOM_EOS_Roquet_rho.F90
   src/equation_of_state/MOM_EOS_linear.F90
   src/equation_of_state/MOM_EOS_TEOS10.F90
   src/equation_of_state/MOM_EOS_UNESCO.F90
   src/equation_of_state/MOM_EOS_Wright.F90
   src/equation_of_state/MOM_EOS_Wright_full.F90
   src/equation_of_state/MOM_EOS_Wright_red.F90
   src/equation_of_state/MOM_TFreeze.F90
   src/equation_of_state/MOM_temperature_convert.F90
#  src/equation_of_state/TEOS10/gsw_chem_potential_water_t_exact.f90
#  src/equation_of_state/TEOS10/gsw_ct_freezing_exact.f90
#  src/equation_of_state/TEOS10/gsw_ct_freezing_poly.f90
#  src/equation_of_state/TEOS10/gsw_ct_from_pt.f90
#  src/equation_of_state/TEOS10/gsw_ct_from_t.f90
#  src/equation_of_state/TEOS10/gsw_entropy_part.f90
#  src/equation_of_state/TEOS10/gsw_entropy_part_zerop.f90
#  src/equation_of_state/TEOS10/gsw_gibbs.f90
#  src/equation_of_state/TEOS10/gsw_gibbs_ice.f90
#  src/equation_of_state/TEOS10/gsw_gibbs_pt0_pt0.f90
#  src/equation_of_state/TEOS10/gsw_mod_error_functions.f90
#  src/equation_of_state/TEOS10/gsw_mod_freezing_poly_coefficients.f90
#  src/equation_of_state/TEOS10/gsw_mod_gibbs_ice_coefficients.f90
#  src/equation_of_state/TEOS10/gsw_mod_kinds.f90
#  src/equation_of_state/TEOS10/gsw_mod_specvol_coefficients.f90
#  src/equation_of_state/TEOS10/gsw_mod_teos10_constants.f90
#  src/equation_of_state/TEOS10/gsw_mod_toolbox.f90
#  src/equation_of_state/TEOS10/gsw_pt0_from_t.f90
#  src/equation_of_state/TEOS10/gsw_pt_from_ct.f90
#  src/equation_of_state/TEOS10/gsw_pt_from_t.f90
#  src/equation_of_state/TEOS10/gsw_rho.f90
#  src/equation_of_state/TEOS10/gsw_rho_first_derivatives.f90
#  src/equation_of_state/TEOS10/gsw_rho_second_derivatives.f90
#  src/equation_of_state/TEOS10/gsw_specvol.f90
#  src/equation_of_state/TEOS10/gsw_specvol_first_derivatives.f90
#  src/equation_of_state/TEOS10/gsw_specvol_second_derivatives.f90
#  src/equation_of_state/TEOS10/gsw_sp_from_sr.f90
#  src/equation_of_state/TEOS10/gsw_sr_from_sp.f90
#  src/equation_of_state/TEOS10/gsw_t_deriv_chem_potential_water_t_exact.f90
#  src/equation_of_state/TEOS10/gsw_t_freezing_exact.f90
#  src/equation_of_state/TEOS10/gsw_t_freezing_poly.f90
#  src/equation_of_state/TEOS10/gsw_t_from_ct.f90
   src/framework/MOM_array_transform.F90
   src/framework/MOM_checksums.F90
   src/framework/MOM_coms.F90
   src/framework/MOM_coupler_types.F90
   src/framework/MOM_cpu_clock.F90
   src/framework/MOM_data_override.F90
   src/framework/MOM_diag_mediator.F90
   src/framework/MOM_diag_remap.F90
   src/framework/MOM_document.F90
   src/framework/MOM_domains.F90
   src/framework/MOM_dyn_horgrid.F90
   src/framework/MOM_ensemble_manager.F90
   src/framework/MOM_error_handler.F90
   src/framework/MOM_file_parser.F90
   src/framework/MOM_get_input.F90
   src/framework/MOM_hor_index.F90
   src/framework/MOM_horizontal_regridding.F90
   src/framework/MOM_interpolate.F90
   src/framework/MOM_intrinsic_functions.F90
   src/framework/MOM_io.F90
   src/framework/MOM_io_file.F90
   src/framework/MOM_murmur_hash.F90
   src/framework/MOM_netcdf.F90
   src/framework/MOM_random.F90
   src/framework/MOM_restart.F90
   src/framework/MOM_safe_alloc.F90
   src/framework/MOM_string_functions.F90
   src/framework/MOM_unique_scales.F90
   src/framework/MOM_unit_scaling.F90
   src/framework/MOM_write_cputime.F90
   src/framework/numerical_testing_type.F90
   src/framework/posix.F90
   src/framework/posix.h
   src/ice_shelf/MOM_ice_shelf_diag_mediator.F90
   src/ice_shelf/MOM_ice_shelf_dynamics.F90
   src/ice_shelf/MOM_ice_shelf.F90
   src/ice_shelf/MOM_ice_shelf_initialize.F90
   src/ice_shelf/MOM_ice_shelf_state.F90
   src/ice_shelf/MOM_marine_ice.F90
   src/ice_shelf/user_shelf_init.F90
   src/initialization/MOM_coord_initialization.F90
   src/initialization/MOM_fixed_initialization.F90
   src/initialization/MOM_grid_initialize.F90
   src/initialization/MOM_shared_initialization.F90
   src/initialization/MOM_state_initialization.F90
   src/initialization/MOM_tracer_initialization_from_Z.F90
   src/ocean_data_assim/MOM_oda_driver.F90
   src/ocean_data_assim/MOM_oda_incupd.F90
   src/parameterizations/lateral/MOM_hor_visc.F90
   src/parameterizations/lateral/MOM_internal_tides.F90
   src/parameterizations/lateral/MOM_interface_filter.F90
   src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90
   src/parameterizations/lateral/MOM_load_love_numbers.F90
   src/parameterizations/lateral/MOM_MEKE.F90
   src/parameterizations/lateral/MOM_MEKE_types.F90
   src/parameterizations/lateral/MOM_mixed_layer_restrat.F90
   src/parameterizations/lateral/MOM_self_attr_load.F90
   src/parameterizations/lateral/MOM_streaming_filter.F90
   src/parameterizations/lateral/MOM_thickness_diffuse.F90
   src/parameterizations/lateral/MOM_spherical_harmonics.F90
   src/parameterizations/lateral/MOM_tidal_forcing.F90
   src/parameterizations/lateral/MOM_wave_drag.F90
   src/parameterizations/lateral/MOM_Zanna_Bolton.F90
   src/parameterizations/vertical/MOM_ALE_sponge.F90
   src/parameterizations/vertical/MOM_bkgnd_mixing.F90
   src/parameterizations/vertical/MOM_bulk_mixed_layer.F90
   src/parameterizations/vertical/MOM_CVMix_conv.F90
   src/parameterizations/vertical/MOM_CVMix_ddiff.F90
   src/parameterizations/vertical/MOM_CVMix_KPP.F90
   src/parameterizations/vertical/MOM_CVMix_shear.F90
   src/parameterizations/vertical/MOM_diabatic_aux.F90
   src/parameterizations/vertical/MOM_diabatic_driver.F90
   src/parameterizations/vertical/MOM_diapyc_energy_req.F90
   src/parameterizations/vertical/MOM_energetic_PBL.F90
   src/parameterizations/vertical/MOM_entrain_diffusive.F90
   src/parameterizations/vertical/MOM_full_convection.F90
   src/parameterizations/vertical/MOM_geothermal.F90
   src/parameterizations/vertical/MOM_internal_tide_input.F90
   src/parameterizations/vertical/MOM_kappa_shear.F90
   src/parameterizations/vertical/MOM_opacity.F90
   src/parameterizations/vertical/MOM_regularize_layers.F90
   src/parameterizations/vertical/MOM_set_diffusivity.F90
   src/parameterizations/vertical/MOM_set_viscosity.F90
   src/parameterizations/vertical/MOM_sponge.F90
   src/parameterizations/vertical/MOM_tidal_mixing.F90
   src/parameterizations/vertical/MOM_vert_friction.F90
   src/parameterizations/stochastic/MOM_stochastics.F90
   src/tracer/advection_test_tracer.F90
   src/tracer/boundary_impulse_tracer.F90
   src/tracer/DOME_tracer.F90
   src/tracer/dyed_obc_tracer.F90
   src/tracer/dye_example.F90
   src/tracer/ideal_age_example.F90
   src/tracer/ISOMIP_tracer.F90
   src/tracer/MARBL_forcing_mod.F90
   src/tracer/MARBL_tracers.F90
   src/tracer/MOM_CFC_cap.F90
   src/tracer/MOM_hor_bnd_diffusion.F90
   src/tracer/MOM_neutral_diffusion.F90
   src/tracer/MOM_OCMIP2_CFC.F90
   src/tracer/MOM_offline_aux.F90
   src/tracer/MOM_offline_main.F90
   src/tracer/MOM_tracer_advect.F90
   src/tracer/MOM_tracer_advect_schemes.F90
   src/tracer/MOM_tracer_diabatic.F90
   src/tracer/MOM_tracer_flow_control.F90
   src/tracer/MOM_tracer_hor_diff.F90
   src/tracer/MOM_tracer_registry.F90
   src/tracer/MOM_tracer_types.F90
   src/tracer/MOM_tracer_Z_init.F90
   src/tracer/nw2_tracers.F90
   src/tracer/oil_tracer.F90
   src/tracer/pseudo_salt_tracer.F90
   src/tracer/RGC_tracer.F90
   src/tracer/tracer_example.F90
   src/user/adjustment_initialization.F90
   src/user/baroclinic_zone_initialization.F90
   src/user/basin_builder.F90
   src/user/benchmark_initialization.F90
   src/user/BFB_initialization.F90
   src/user/BFB_surface_forcing.F90
   src/user/circle_obcs_initialization.F90
   src/user/dense_water_initialization.F90
   src/user/DOME2d_initialization.F90
   src/user/DOME_initialization.F90
   src/user/dumbbell_initialization.F90
   src/user/dumbbell_surface_forcing.F90
   src/user/dyed_channel_initialization.F90
   src/user/dyed_obcs_initialization.F90
   src/user/external_gwave_initialization.F90
   src/user/Idealized_Hurricane.F90
   src/user/ISOMIP_initialization.F90
   src/user/Kelvin_initialization.F90
   src/user/lock_exchange_initialization.F90
   src/user/MOM_controlled_forcing.F90
   src/user/MOM_wave_interface.F90
   src/user/Neverworld_initialization.F90
   src/user/Phillips_initialization.F90
   src/user/RGC_initialization.F90
   src/user/Rossby_front_2d_initialization.F90
   src/user/SCM_CVMix_tests.F90
   src/user/seamount_initialization.F90
   src/user/shelfwave_initialization.F90
   src/user/sloshing_initialization.F90
   src/user/soliton_initialization.F90
   src/user/supercritical_initialization.F90
   src/user/tidal_bay_initialization.F90
   src/user/user_change_diffusivity.F90
   src/user/user_initialization.F90
   src/user/user_revise_forcing.F90
#  pkg/CVMix-src/src/cvmix_driver.F90
   pkg/CVMix-src/src/cvmix_io.F90
#  pkg/CVMix-src/src/dummy.F90
   pkg/CVMix-src/src/drivers/cvmix_bgrnd_BL.F90
   pkg/CVMix-src/src/drivers/cvmix_ddiff_drv.F90
   pkg/CVMix-src/src/drivers/cvmix_kpp_drv.F90
   pkg/CVMix-src/src/drivers/cvmix_shear_drv.F90
   pkg/CVMix-src/src/drivers/cvmix_tidal_Simmons.F90
   pkg/CVMix-src/src/shared/cvmix_background.F90
   pkg/CVMix-src/src/shared/cvmix_convection.F90
   pkg/CVMix-src/src/shared/cvmix_ddiff.F90
   pkg/CVMix-src/src/shared/cvmix_kinds_and_types.F90
   pkg/CVMix-src/src/shared/cvmix_kpp.F90
   pkg/CVMix-src/src/shared/cvmix_math.F90
   pkg/CVMix-src/src/shared/cvmix_put_get.F90
   pkg/CVMix-src/src/shared/cvmix_shear.F90
   pkg/CVMix-src/src/shared/cvmix_tidal.F90
   pkg/CVMix-src/src/shared/cvmix_utils.F90
   pkg/GSW-Fortran/modules/gsw_mod_baltic_data.f90
   pkg/GSW-Fortran/modules/gsw_mod_error_functions.f90
   pkg/GSW-Fortran/modules/gsw_mod_freezing_poly_coefficients.f90
   pkg/GSW-Fortran/modules/gsw_mod_gibbs_ice_coefficients.f90
   pkg/GSW-Fortran/modules/gsw_mod_kinds.f90
   pkg/GSW-Fortran/modules/gsw_mod_netcdf.f90
   pkg/GSW-Fortran/modules/gsw_mod_saar_data.f90
   pkg/GSW-Fortran/modules/gsw_mod_sp_coefficients.f90
   pkg/GSW-Fortran/modules/gsw_mod_specvol_coefficients.f90
   pkg/GSW-Fortran/modules/gsw_mod_teos10_constants.f90
   pkg/GSW-Fortran/modules/gsw_mod_toolbox.f90
#  pkg/GSW-Fortran/test/gsw_check_functions.f90
#  pkg/GSW-Fortran/test/gsw_poly_check.f90
   pkg/GSW-Fortran/toolbox/gsw_add_barrier.f90
   pkg/GSW-Fortran/toolbox/gsw_add_mean.f90
   pkg/GSW-Fortran/toolbox/gsw_adiabatic_lapse_rate_from_ct.f90
   pkg/GSW-Fortran/toolbox/gsw_adiabatic_lapse_rate_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_alpha.f90
   pkg/GSW-Fortran/toolbox/gsw_alpha_on_beta.f90
   pkg/GSW-Fortran/toolbox/gsw_alpha_wrt_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_alpha_wrt_t_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_beta_const_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_beta.f90
   pkg/GSW-Fortran/toolbox/gsw_cabbeling.f90
   pkg/GSW-Fortran/toolbox/gsw_c_from_sp.f90
   pkg/GSW-Fortran/toolbox/gsw_chem_potential_water_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_chem_potential_water_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_cp_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_first_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_first_derivatives_wrt_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_freezing_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_freezing.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_freezing_first_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_freezing_first_derivatives_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_freezing_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_from_enthalpy_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_from_enthalpy.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_from_entropy.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_from_pt.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_from_rho.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_from_t.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_maxdensity.f90
   pkg/GSW-Fortran/toolbox/gsw_ct_second_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_deltasa_atlas.f90
   pkg/GSW-Fortran/toolbox/gsw_deltasa_from_sp.f90
   pkg/GSW-Fortran/toolbox/gsw_dilution_coefficient_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_dynamic_enthalpy.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy_ct_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy_diff.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy_first_derivatives_ct_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy_first_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy_second_derivatives_ct_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy_second_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy_sso_0.f90
   pkg/GSW-Fortran/toolbox/gsw_enthalpy_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_entropy_first_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_entropy_from_pt.f90
   pkg/GSW-Fortran/toolbox/gsw_entropy_from_t.f90
   pkg/GSW-Fortran/toolbox/gsw_entropy_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_entropy_part.f90
   pkg/GSW-Fortran/toolbox/gsw_entropy_part_zerop.f90
   pkg/GSW-Fortran/toolbox/gsw_entropy_second_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_fdelta.f90
   pkg/GSW-Fortran/toolbox/gsw_frazil_properties.f90
   pkg/GSW-Fortran/toolbox/gsw_frazil_properties_potential.f90
   pkg/GSW-Fortran/toolbox/gsw_frazil_properties_potential_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_frazil_ratios_adiabatic.f90
   pkg/GSW-Fortran/toolbox/gsw_frazil_ratios_adiabatic_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_geo_strf_dyn_height.f90
   pkg/GSW-Fortran/toolbox/gsw_geo_strf_dyn_height_pc.f90
   pkg/GSW-Fortran/toolbox/gsw_gibbs.f90
   pkg/GSW-Fortran/toolbox/gsw_gibbs_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_gibbs_ice_part_t.f90
   pkg/GSW-Fortran/toolbox/gsw_gibbs_ice_pt0.f90
   pkg/GSW-Fortran/toolbox/gsw_gibbs_ice_pt0_pt0.f90
   pkg/GSW-Fortran/toolbox/gsw_gibbs_pt0_pt0.f90
   pkg/GSW-Fortran/toolbox/gsw_grav.f90
   pkg/GSW-Fortran/toolbox/gsw_helmholtz_energy_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_hill_ratio_at_sp2.f90
   pkg/GSW-Fortran/toolbox/gsw_ice_fraction_to_freeze_seawater.f90
   pkg/GSW-Fortran/toolbox/gsw_internal_energy.f90
   pkg/GSW-Fortran/toolbox/gsw_internal_energy_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_ipv_vs_fnsquared_ratio.f90
   pkg/GSW-Fortran/toolbox/gsw_kappa_const_t_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_kappa.f90
   pkg/GSW-Fortran/toolbox/gsw_kappa_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_kappa_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_latentheat_evap_ct.f90
   pkg/GSW-Fortran/toolbox/gsw_latentheat_evap_t.f90
   pkg/GSW-Fortran/toolbox/gsw_latentheat_melting.f90
   pkg/GSW-Fortran/toolbox/gsw_linear_interp_sa_ct.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_ice_equilibrium_sa_ct_ratio.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_ice_equilibrium_sa_ct_ratio_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_ice_into_seawater.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_ice_sa_ct_ratio.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_ice_sa_ct_ratio_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_seaice_equilibrium_sa_ct_ratio.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_seaice_equilibrium_sa_ct_ratio_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_seaice_into_seawater.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_seaice_sa_ct_ratio.f90
   pkg/GSW-Fortran/toolbox/gsw_melting_seaice_sa_ct_ratio_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_nsquared.f90
   pkg/GSW-Fortran/toolbox/gsw_pot_enthalpy_from_pt_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_pot_enthalpy_from_pt_ice_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_pot_enthalpy_ice_freezing.f90
   pkg/GSW-Fortran/toolbox/gsw_pot_enthalpy_ice_freezing_first_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_pot_enthalpy_ice_freezing_first_derivatives_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_pot_enthalpy_ice_freezing_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_pot_rho_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_pressure_coefficient_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_pressure_freezing_ct.f90
   pkg/GSW-Fortran/toolbox/gsw_pt0_cold_ice_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_pt0_from_t.f90
   pkg/GSW-Fortran/toolbox/gsw_pt0_from_t_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_pt_first_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_pt_from_ct.f90
   pkg/GSW-Fortran/toolbox/gsw_pt_from_entropy.f90
   pkg/GSW-Fortran/toolbox/gsw_pt_from_pot_enthalpy_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_pt_from_pot_enthalpy_ice_poly_dh.f90
   pkg/GSW-Fortran/toolbox/gsw_pt_from_pot_enthalpy_ice_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_pt_from_t.f90
   pkg/GSW-Fortran/toolbox/gsw_pt_from_t_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_pt_second_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_rho_alpha_beta.f90
   pkg/GSW-Fortran/toolbox/gsw_rho.f90
   pkg/GSW-Fortran/toolbox/gsw_rho_first_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_rho_first_derivatives_wrt_enthalpy.f90
   pkg/GSW-Fortran/toolbox/gsw_rho_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_rho_second_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_rho_second_derivatives_wrt_enthalpy.f90
   pkg/GSW-Fortran/toolbox/gsw_rho_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_rr68_interp_sa_ct.f90
   pkg/GSW-Fortran/toolbox/gsw_saar.f90
   pkg/GSW-Fortran/toolbox/gsw_saar_init.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_freezing_estimate.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_freezing_from_ct.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_freezing_from_ct_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_freezing_from_t.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_freezing_from_t_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_from_rho.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_from_sp_baltic.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_from_sp.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_from_sstar.f90
   pkg/GSW-Fortran/toolbox/gsw_sa_p_inrange.f90
   pkg/GSW-Fortran/toolbox/gsw_seaice_fraction_to_freeze_seawater.f90
   pkg/GSW-Fortran/toolbox/gsw_sigma0.f90
   pkg/GSW-Fortran/toolbox/gsw_sigma1.f90
   pkg/GSW-Fortran/toolbox/gsw_sigma2.f90
   pkg/GSW-Fortran/toolbox/gsw_sigma3.f90
   pkg/GSW-Fortran/toolbox/gsw_sigma4.f90
   pkg/GSW-Fortran/toolbox/gsw_sound_speed.f90
   pkg/GSW-Fortran/toolbox/gsw_sound_speed_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_sound_speed_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol_alpha_beta.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol_anom_standard.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol_first_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol_first_derivatives_wrt_enthalpy.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol_second_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol_second_derivatives_wrt_enthalpy.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol_sso_0.f90
   pkg/GSW-Fortran/toolbox/gsw_specvol_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_sp_from_c.f90
   pkg/GSW-Fortran/toolbox/gsw_sp_from_sa_baltic.f90
   pkg/GSW-Fortran/toolbox/gsw_sp_from_sa.f90
   pkg/GSW-Fortran/toolbox/gsw_sp_from_sk.f90
   pkg/GSW-Fortran/toolbox/gsw_sp_from_sr.f90
   pkg/GSW-Fortran/toolbox/gsw_sp_from_sstar.f90
   pkg/GSW-Fortran/toolbox/gsw_spiciness0.f90
   pkg/GSW-Fortran/toolbox/gsw_spiciness1.f90
   pkg/GSW-Fortran/toolbox/gsw_spiciness2.f90
   pkg/GSW-Fortran/toolbox/gsw_sr_from_sp.f90
   pkg/GSW-Fortran/toolbox/gsw_sstar_from_sa.f90
   pkg/GSW-Fortran/toolbox/gsw_sstar_from_sp.f90
   pkg/GSW-Fortran/toolbox/gsw_t_deriv_chem_potential_water_t_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_t_freezing_exact.f90
   pkg/GSW-Fortran/toolbox/gsw_t_freezing.f90
   pkg/GSW-Fortran/toolbox/gsw_t_freezing_first_derivatives.f90
   pkg/GSW-Fortran/toolbox/gsw_t_freezing_first_derivatives_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_t_freezing_poly.f90
   pkg/GSW-Fortran/toolbox/gsw_t_from_ct.f90
   pkg/GSW-Fortran/toolbox/gsw_t_from_pt0_ice.f90
   pkg/GSW-Fortran/toolbox/gsw_thermobaric.f90
   pkg/GSW-Fortran/toolbox/gsw_turner_rsubrho.f90
   pkg/GSW-Fortran/toolbox/gsw_util_indx.f90
   pkg/GSW-Fortran/toolbox/gsw_util_interp1q_int.f90
   pkg/GSW-Fortran/toolbox/gsw_util_sort_real.f90
   pkg/GSW-Fortran/toolbox/gsw_util_xinterp1.f90
   pkg/GSW-Fortran/toolbox/gsw_z_from_p.f90
#  choose above set MOM6_INFRA interface
   config_src/infra/${MOM6_INFRA}/MOM_coms_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_constants.F90
   config_src/infra/${MOM6_INFRA}/MOM_couplertype_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_cpu_clock_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_data_override_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_diag_manager_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_domain_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_ensemble_manager_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_error_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_interp_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_io_infra.F90
   config_src/infra/${MOM6_INFRA}/MOM_time_manager.F90
)

list( APPEND MOM6_SRCS
   config_src/drivers/FMS_cap/MOM_surface_forcing_gfdl.F90
   config_src/drivers/FMS_cap/ocean_model_MOM.F90
   # oda
   config_src/external/ODA_hooks/kdtree.f90
   config_src/external/ODA_hooks/ocean_da_core.F90
   config_src/external/ODA_hooks/ocean_da_types.F90
   config_src/external/ODA_hooks/write_ocean_obs.F90
   # tracer
   config_src/external/GFDL_ocean_BGC/FMS_coupler_util.F90
   config_src/external/GFDL_ocean_BGC/generic_tracer.F90
   config_src/external/GFDL_ocean_BGC/generic_tracer_utils.F90
   config_src/external/GFDL_ocean_BGC/MOM_generic_tracer.F90
   # drifters-particles
   config_src/external/drifters/MOM_particles.F90
   config_src/external/drifters/MOM_particles_types.F90
   # stochastic physics
   config_src/external/stochastic_physics/get_stochy_pattern.F90
   config_src/external/stochastic_physics/stochastic_physics.F90
   # database comms
   config_src/external/database_comms/database_client_interface.F90
   config_src/external/database_comms/MOM_database_comms.F90
   # MARBL
   config_src/external/MARBL/marbl_constants_mod.F90
   config_src/external/MARBL/marbl_interface.F90
   config_src/external/MARBL/marbl_interface_public_types.F90
   config_src/external/MARBL/marbl_logging.F90
)


esma_mepo_style(mom6 MOM6_rel_path REL_PATH ..)
set (MOM6_path ${CMAKE_CURRENT_SOURCE_DIR}/${MOM6_rel_path})

set (SRCS)
foreach (file ${MOM6_SRCS})
  list (APPEND SRCS ${MOM6_path}/${file})
endforeach ()

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES FMS::fms_r8
  INCLUDES
  $<BUILD_INTERFACE:${MOM6_path}/config_src/memory/dynamic_nonsymmetric>
#  choose above set MOM6_infra interface
  $<BUILD_INTERFACE:${MOM6_path}/config_src/infra/${MOM6_INFRA}>
  $<BUILD_INTERFACE:${MOM6_path}/src/framework>
  TYPE SHARED
)

target_compile_definitions (${this} PRIVATE use_netCDF)

# This 'resets' the Intel DEBUG flags for MOM6. The stock debug flags use
# 'all,noarg_temp_created' which seem to be too aggressive for FMS/MOM6. This
# moves them back to the 'bounds,uninit' GEOS used to build with.
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel" AND CMAKE_BUILD_TYPE MATCHES Debug)
   string(REPLACE "all,noarg_temp_created" "bounds,uninit" _tmp "${GEOS_Fortran_FLAGS_DEBUG}")
   set (CMAKE_Fortran_FLAGS_DEBUG "${_tmp}")
endif ()

# Specs for r8 version
string(REPLACE " " ";" tmp ${FREAL8})
foreach (flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach ()
