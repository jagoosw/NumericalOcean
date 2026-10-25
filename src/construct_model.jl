using Oceananigans.Units
using WorldOceanAtlasTools
using OMIPSimulations: henyey_diffusivity_field, νhb, woa_to_teos10!, salinity_surface_restoring, NormalizeTotalWater,
                       omip_progress_callback
using Oceananigans.TurbulenceClosures: TriadIsopycnalSkewSymmetricDiffusivity, FluxTapering
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.OutputReaders: Cyclical
using OceanBioME.Models.GasExchangeModel.ScaledGasTransferVelocity: JRA55
using ClimaSeaIce: IncrementalRemapping
using Oceananigans.Architectures: architecture
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength
using NumericalEarth.DataWrangling.SeaWiFS: SeaWiFSMonthly
using NumericalEarth.DataWrangling.WOA: WOAClimatology
using NumericalEarth.DataWrangling: ConservativeSurfaceFluxRestoringCallback
using NumericalEarth.EarthSystemModels.InterfaceComputations: COARELogarithmicSimilarityProfile, ConvectiveGustiness,
                                                              WindDependentWaveFormulation, TemperatureDependentAirViscosity,
                                                              atmosphere_sea_ice_stability_functions,
                                                              ImpureSaturationSpecificHumidity, AtmosphericThermodynamics,
                                                              ThreeEquationHeatFlux, MomentumBasedFrictionVelocity,
                                                              InterfaceTemperatureMeltwater, GillSaturationEnhancement

"""
$(TYPEDSIGNATURES)

Build a JRA55-forced ocean--sea-ice `Simulation` on `grid`,
optionally carrying `biogeochemistry` with OAEMIP's initial conditions, rivers and iron dust.
Extra passive `tracers` are added to `(:T, :S)`.

Model time zero is `start_date`. The JRA55 forcing cycles every `forcing_cycle` seconds (see [`forcing_period`](@ref):
21915 days for 1958--2018, 365 for repeat-year forcing) and the climatologies (WOA salinity, SeaWiFS chlorophyll,
iron dust) every `climatology_cycle` (see [`climatology_period`](@ref): 365.25 days, or 365), so all stay in phase
with each other and the calendar from cycle to cycle, the climatologies to within 0.75 days. Cycle any
`atmosphere_tracers` the same way, e.g. `MaunaLoaCO₂(arch; dates, start_date, period = forcing_cycle)`.

Every `progress_frequency` iterations (`nothing` for never) the time, extrema of the state and wall time per step are
logged, followed by a `PROBE` line of memory and GC use (unless `ENV["OMIP_PROBE"] = "0"`). Its reductions synchronise
the GPU, so keep it infrequent.
"""
function forced_ocean_simulation(grid; 
                                 arch = architecture(grid),
                                 backend_size = 50,
                                 start_date = DateTime(1958, 1, 1),
                                 end_date = DateTime(2018, 1, 1),
                                 atmosphere_tracers = NamedTuple(),
                                 staging = false,
                                 jra55_dataset = MultiYearJRA55(),
                                 jra55_dir = staging ? stage_jra55!(jra55_dataset) : forcing_dir[],
                                 forcing_cycle = forcing_period(jra55_dataset, start_date, end_date),
                                 climatology_cycle = climatology_period(forcing_cycle),
                                 land = JRA55PrescribedLand(grid; 
                                                            dir = jra55_dir, 
                                                            dataset = jra55_dataset,
                                                            start_date, 
                                                            end_date, 
                                                            time_indices_in_memory = backend_size, 
                                                            time_indexing = Cyclical(forcing_cycle),
                                                            prefetch = true,
                                                            maximum_search_radius = default_maximum_search_radius(grid),
                                                            spread_radius = default_river_spread_radius(grid),
                                                            maximum_spread_cells = default_river_spread_cells(grid)),

                                 eddy_closure = grid isa ORCA1GRID ? TriadIsopycnalSkewSymmetricDiffusivity(VerticallyImplicitTimeDiscretization();
                                                                                                            κ_skew = 800, κ_symmetric = 800, 
                                                                                                            slope_limiter = FluxTapering(1e-2))
                                                                   : nothing,
                                 background_diffusivity = VerticalScalarDiffusivity(κ = henyey_diffusivity_field(grid), ν = 3e-5),
                                 vertical_closure = CATKEVerticalDiffusivity(VerticallyImplicitTimeDiscretization(); 
                                                                             mixing_length = CATKEMixingLength(Cᵇ = 0.01), 
                                                                             maximum_viscosity = 3, 
                                                                             maximum_tracer_diffusivity = 3, 
                                                                             maximum_tke_diffusivity = 3, 
                                                                             negative_tke_damping_time_scale = 10),
                                 horizontal_viscosity = grid isa ORCA1GRID ? 
                                                        HorizontalScalarBiharmonicDiffusivity(ν = νhb, κ = nothing, discrete_form = true,
                                                                                              parameters = 50days) :
                                                        nothing,
                                 closures = (vertical_closure, eddy_closure, horizontal_viscosity, background_diffusivity),
                                 salt_restoring = salinity_surface_restoring(grid, WOAMonthly(); 
                                                                             restoring_dir = restoring_dir[], 
                                                                             piston_velocity = 0.254,
                                                                             time_indexing = Cyclical(climatology_cycle)),
                                 chlorophyll = seawifs_chlorophyll_climatology(grid; climatology_cycle),
                                 radiative_forcing = TwoColorRadiation(grid; chlorophyll),
                                 momentum_advection = WENOVectorInvariant(; vorticity_order = grid isa ORCA1GRID ? 5 : 9,
                                                                            time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)),
                                 biogeochemistry = nothing,
                                 tracers = nothing,
                                 forcing = NamedTuple(),
                                 river_tracers = filter(tracer -> haskey(DEFAULT_RIVER_INPUTS, replicated_tracer(tracer)), required_biogeochemical_tracers(biogeochemistry)),
                                 river_fluxes = isempty(river_tracers) ? NamedTuple() :
                                                RivR2OSurfaceFlux(grid; forced_tracers = river_tracers, climatology = true, start_year = 1985, end_year = 1995, land),
                                 dust_deposition = :Fe ∈ required_biogeochemical_tracers(biogeochemistry) ?
                                                   (; Fe = InterpolatedIronDust(BGCInitDustIronDeposition(grid; time_indexing = Cyclical(climatology_cycle));
                                                                                hard_fraction = 0.987)) : NamedTuple(),
                                 time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5),
                                 tracer_advection = tracer_advection_schemes(tracers, biogeochemistry, time_discretization),
                                 ocean = ocean_simulation(grid;
                                                          Δt = 1minutes,
                                                          tracers = (:T, :S, something(tracers, ())...),
                                                          radiative_forcing,
                                                          momentum_advection,
                                                          tracer_advection,
                                                          materialize_buoyancy_gradients = true,
                                                          free_surface = SplitExplicitFreeSurface(grid; substeps = default_barotropic_substeps(grid)),
                                                          additional_surface_fluxes = merge((; S = salt_restoring), river_fluxes),
                                                          closure = filter(!isnothing, closures),
                                                          river_routing = land.river_routing,
                                                          biogeochemistry,
                                                          forcing = merge(forcing, dust_deposition)),

                                 sea_ice_immersed_latitudes = (-50, 35),
                                 sea_ice_grid = build_sea_ice_grid(grid, sea_ice_immersed_latitudes, GridFittedBottom),
                                 sea_ice = sea_ice_simulation(sea_ice_grid, ocean;
                                                              advection = IncrementalRemapping(),
                                                              thickness_categories = 4),
                                 atmosphere = JRA55PrescribedAtmosphere(arch; 
                                                                        dataset = jra55_dataset,
                                                                        tracers = atmosphere_tracers,
                                                                        dir = jra55_dir,
                                                                        start_date,
                                                                        end_date,
                                                                        time_indexing = Cyclical(forcing_cycle),
                                                                        time_indices_in_memory = backend_size),
                                 sea_ice_albedo = SeaIceAlbedo(sea_ice.model.ice_thickness,
                                                               sea_ice.model.snow_thickness,
                                                               sea_ice.model.snow_thermodynamics.top_surface_temperature),
                                 radiation = JRA55PrescribedRadiation(arch;
                                                                      dataset = jra55_dataset,
                                                                      dir = jra55_dir,
                                                                      start_date,
                                                                      end_date,
                                                                      time_indices_in_memory = backend_size,
                                                                      time_indexing = Cyclical(forcing_cycle),
                                                                      prefetch = true,
                                                                      ocean_surface   = SurfaceRadiationProperties(0.06, 1.00),
                                                                      sea_ice_surface = SurfaceRadiationProperties(sea_ice_albedo, 1.0)),
                                 atmosphere_correction = nothing,
                                 radiation_correction = nothing,
                                 biogeochemistry_interface_kwargs = (; base_transfer_velocity = JRA55()),
                                 air_kinematic_viscosity = TemperatureDependentAirViscosity(),
                                 similarity_form = COARELogarithmicSimilarityProfile(),
                                 atmosphere_ocean_fluxes = SimilarityTheoryFluxes(; similarity_form,
                                                                                    subgrid_velocities           = ConvectiveGustiness(minimum_gustiness = 0.5),
                                                                                    momentum_roughness_length    = MomentumRoughnessLength(; wave_formulation = WindDependentWaveFormulation(), air_kinematic_viscosity),
                                                                                    temperature_roughness_length = ScalarRoughnessLength(; air_kinematic_viscosity),
                                                                                    water_vapor_roughness_length = ScalarRoughnessLength(; air_kinematic_viscosity)),
                                 atmosphere_sea_ice_fluxes = SimilarityTheoryFluxes(; stability_functions = atmosphere_sea_ice_stability_functions(),
                                                                                      similarity_form,
                                                                                      subgrid_velocities = ConvectiveGustiness(minimum_gustiness = 0.2),
                                                                                      momentum_roughness_length = 5e-4,
                                                                                      temperature_roughness_length = 5e-5,
                                                                                      water_vapor_roughness_length = 5e-5),
                                 sea_ice_ocean_heat_flux = ThreeEquationHeatFlux(; heat_transfer_coefficient = 0.0057,
                                                                                   friction_velocity = MomentumBasedFrictionVelocity()),
                                 ice_meltwater_enthalpy = InterfaceTemperatureMeltwater(),
                                 saturation_enhancement = GillSaturationEnhancement(),
                                 humidity = (ocean = ImpureSaturationSpecificHumidity(AtmosphericThermodynamics.Liquid(), 0.98; saturation_enhancement),
                                             sea_ice = ImpureSaturationSpecificHumidity(AtmosphericThermodynamics.Ice(); saturation_enhancement)),
                                 interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice; radiation, land,
                                                                  atmosphere_ocean_fluxes,
                                                                  atmosphere_sea_ice_fluxes,
                                                                  sea_ice_ocean_heat_flux,
                                                                  ice_meltwater_enthalpy,
                                                                  atmosphere_ocean_interface_specific_humidity = humidity.ocean,
                                                                  atmosphere_sea_ice_interface_specific_humidity = humidity.sea_ice,
                                                                  exchanger_correction = atmosphere_correction, radiation_correction,
                                                                  biogeochemistry_interface_kwargs),
                                 Δt = default_Δt(grid),
                                 stop_time = 5 * 21915days,
                                 progress_frequency = 100,
                                 T_and_S_initialisation_dataset = WOAAnnual(),
                                 initialise_bgc = true)

    if T_and_S_initialisation_dataset isa WOAClimatology
        T_init = Field(Metadatum(:temperature; dir = restoring_dir[], dataset = T_and_S_initialisation_dataset), grid)
        S_init = Field(Metadatum(:salinity;    dir = restoring_dir[], dataset = T_and_S_initialisation_dataset), grid)
        woa_to_teos10!(T_init, S_init)
        set!(ocean.model, T=T_init, S=S_init)
    elseif !isnothing(T_and_S_initialisation_dataset)
        @warn "Initial temperature and salinity from $(nameof(typeof(T_and_S_initialisation_dataset))) are set without conversion, " *
              "but the model uses TEOS-10 Conservative Temperature and Absolute Salinity (only WOA datasets are converted automatically)"
        set!(ocean.model, T = Metadatum(:temperature; dir = restoring_dir[], dataset = T_and_S_initialisation_dataset),
                          S = Metadatum(:salinity; dir = restoring_dir[], dataset = T_and_S_initialisation_dataset))
    end

    if !isnothing(biogeochemistry) && initialise_bgc
        for tracer in required_biogeochemical_tracers(ocean.model.biogeochemistry)
            tracer ∈ (:T, :S) && continue
            if haskey(BGC_INIT_TRACERS, replicated_tracer(tracer))
                set!(ocean.model.tracers[tracer], 
                     Metadatum(BGC_INIT_TRACERS[replicated_tracer(tracer)]; 
                               dataset = BGCInitClimatology(), 
                               date = start_date))
            else
                @warn "BGCInitClimatology has no initial condition for $tracer"
            end
        end
    end

    set!(sea_ice.model, h = Metadatum(:sea_ice_thickness;     dir = restoring_dir[], dataset = ECCO4Monthly(), date = DateTime(1993, 1, 1)),
                        ℵ = Metadatum(:sea_ice_concentration; dir = restoring_dir[], dataset = ECCO4Monthly(), date = DateTime(1993, 1, 1)))

    model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, land, interfaces)

    simulation = Simulation(model; Δt, stop_time)

    if staging && !(jra55_dataset isa RepeatYearJRA55)
        staging_callback = JRA55DataStagingCallback(; source_dir = forcing_dir[],
                                                      staging_dir = staging_dir[],
                                                      start_date,
                                                      cycle_length = forcing_cycle)
        add_callback!(simulation, staging_callback, TimeInterval(30days))
    end

    add_callback!(simulation, ConservativeSurfaceFluxRestoringCallback(salt_restoring, ocean.model), IterationInterval(1))
    add_callback!(simulation, NormalizeTotalWater(model, Δt, nothing), IterationInterval(1))

    if !isnothing(progress_frequency)
        add_callback!(simulation, omip_progress_callback(Ref(time_ns())), IterationInterval(progress_frequency); name = :progress)
    end

    return simulation
end

    
