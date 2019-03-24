class VariableUnits:

    CANONICAL_UNITS = {};

    CANONICAL_UNITS['longitude'] = "degree_east"
    CANONICAL_UNITS['latitude'] = "degree_north"
    CANONICAL_UNITS['depth'] = "m"
    CANONICAL_UNITS['sea_binary_mask'] = "1"
    CANONICAL_UNITS['wet_binary_mask'] = "1"
    CANONICAL_UNITS['mesh_size'] = "m"

    # HYDRO
    # Sea Surface
    CANONICAL_UNITS['sea_surface_height_above_mean_sea_level'] = "m"
    CANONICAL_UNITS['sea_surface_height_above_geoid'] = "m"
    CANONICAL_UNITS['sea_surface_temperature'] = "C"
    CANONICAL_UNITS['sea_surface_salinity'] = "psu"
    CANONICAL_UNITS['sea_surface_pressure'] = "Pa"
    CANONICAL_UNITS['sea_surface_density'] = "kg m-3"
    CANONICAL_UNITS['eastward_sea_water_velocity_at_sea_water_surface'] = "m s-1"
    CANONICAL_UNITS['northward_sea_water_velocity_at_sea_water_surface'] = "m s-1"
    CANONICAL_UNITS['sea_water_speed_at_sea_water_surface'] = "m s-1"
    CANONICAL_UNITS['sea_water_from_direction_at_sea_water_surface'] = "degree"  # from North=0"
    CANONICAL_UNITS['sea_water_to_direction_at_sea_water_surface'] = "degree"  # from North=0"

    # 2D
    CANONICAL_UNITS['bathymetry'] = "m"
    CANONICAL_UNITS['barotropic_eastward_sea_water_velocity'] = "m s-1"
    CANONICAL_UNITS['barotropic_northward_sea_water_velocity'] = "m s-1"
    CANONICAL_UNITS['barotropic_sea_water_speed'] = "m s-1"
    CANONICAL_UNITS['barotropic_sea_water_from_direction'] = "degree"  # from North=0"
    CANONICAL_UNITS['barotropic_sea_water_to_direction'] = "degree"  # from North=0"

    # 3D
    CANONICAL_UNITS['sea_water_turbidity'] = "FTU"
    CANONICAL_UNITS['sea_water_pressure_at_sea_water_surface'] = "dbar"
    CANONICAL_UNITS['sea_water_electrical_conductivity'] = "S m-1"
    CANONICAL_UNITS['sea_water_temperature'] = "degree"
    CANONICAL_UNITS['sea_water_salinity'] = "psu"
    CANONICAL_UNITS['sea_water_density'] = "kg m-3"
    CANONICAL_UNITS['baroclinic_eastward_sea_water_velocity'] = "m s-1"
    CANONICAL_UNITS['baroclinic_northward_sea_water_velocity'] = "m s-1"

    CANONICAL_UNITS['water_volume_transport_into_sea_water_from_rivers'] = "m3 s-1"

    # WAVES
    CANONICAL_UNITS['sea_surface_wave_significant_height'] = "m"
    CANONICAL_UNITS['sea_surface_wave_breaking_height'] = "m"
    CANONICAL_UNITS['sea_surface_wave_mean_period'] = "s"
    CANONICAL_UNITS['sea_surface_wave_peak_period'] = "s"
    CANONICAL_UNITS['sea_surface_wave_from_direction'] = "degree" # from North=0 / East=90"
    CANONICAL_UNITS['sea_surface_wave_to_direction'] = "degree" # from North=0 / East=90"
    CANONICAL_UNITS['eastward_atmosphere_momentum_flux_to_waves'] = "m2 s-1"
    CANONICAL_UNITS['northward_atmosphere_momentum_flux_to_waves'] = "m2 s-1"
    CANONICAL_UNITS['eastward_waves_momentum_flux_to_ocean'] = "m2 s-1"
    CANONICAL_UNITS['northward_waves_momentum_flux_to_ocean'] = "m2 s-1"
    CANONICAL_UNITS['eastward_sea_surface_wave_stokes_drift_velocity'] = "m s-1"
    CANONICAL_UNITS['northward_sea_surface_wave_stokes_drift_velocity'] = "m s-1"
    CANONICAL_UNITS['radiation_pressure_bernouilli_head'] = "m2 s-2"
    CANONICAL_UNITS['sea_surface_wave_energy_flux_to_ocean'] = "W m-2"
    CANONICAL_UNITS['sea_surface_wave_energy_dissipation_at_ground_level'] = "W m-2"

    # METEO
    CANONICAL_UNITS['surface_air_pressure'] = "Pa"
    CANONICAL_UNITS['eastward_wind_10m'] = "m s-1"
    CANONICAL_UNITS['northward_wind_10m'] = "m s-1"
    CANONICAL_UNITS['wind_speed_10m'] = "m s-1"
    CANONICAL_UNITS['wind_to_direction_10m'] = "degree" # to North=0 / East=90"
    CANONICAL_UNITS['wind_from_direction_10m'] = "degree" # from North=0 / East=90"
    CANONICAL_UNITS['rainfall_amount'] = "kg m-1"
    CANONICAL_UNITS['eastward_wind_stress'] = "W m2"
    CANONICAL_UNITS['northward_wind_stress'] = "W m2"
    CANONICAL_UNITS['topography'] = "m"



