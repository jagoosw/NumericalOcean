using JLD2: ZstdFilter
using Oceananigans: fields
using Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields
using Oceananigans.Grids: znodes, znode
using Oceananigans.Utils: on_architecture
using NumericalEarth.Diagnostics: MixedLayerDepthField
using OMIPSimulations: AsyncJLD2Writer, FlushAsyncOutput, flush_output!
using OMIPSimulations: zonal_volume_flux, meridional_volume_flux, zonal_tracer_transport, meridional_tracer_transport,
                       uu_at_ccc, vv_at_ccc, one_at_ccc
using OceanBioME.Models.NutrientsPlanktonDetritusModels.PlanktonModels: community_productivity
using OceanBioME.Models.GasExchangeModel: surface_value
using OceanBioME.Models.CarbonChemistryModel: silicate_concentration, phosphate_concentration
using OceanBioME.Models.NutrientsPlanktonDetritusModels.InorganicCarbonModels: AbstractInorganicCarbon

"""
$(TYPEDSIGNATURES)

Add `writer` to `simulation.output_writers` as `name`. With `async = true` it is wrapped in an `AsyncJLD2Writer`,
which compresses and writes on a background task (so needs `julia --threads ≥ 2`), and a callback is added to
write anything still queued at the end of `run!`. Call `flush_output!(simulation)` in a `finally` around `run!`
too, so that an error does not lose queued writes. A `Checkpointer` waits for the queued writes before it
saves, so add it after the writers (see [`checkpointer!`](@ref)).
"""
function add_output_writer!(simulation, name, writer; async = false)
    if async
        writer = AsyncJLD2Writer(writer)
        simulation.callbacks[:flush_async_output] = Callback(FlushAsyncOutput(), IterationInterval(typemax(Int)))
    end

    simulation.output_writers[name] = writer

    return nothing
end

#####
##### OMIP physical diagnostics (OMIPSimulations' `add_omip_diagnostics!` without the checkpointer)
#####

"""
$(TYPEDSIGNATURES)

Add OMIPSimulations' surface, 3-D field and global-average output writers to `simulation`, written in the
background with `async = true` (see [`add_output_writer!`](@ref)).
"""
function omip_diagnostics!(simulation;
                           field_mean_interval = 5days,
                           surface_averaging_interval = 5days,
                           field_averaging_interval = 15days,
                           averaging_stride = 1,
                           output_dir = NumericalOcean.output_dir[],
                           filename_prefix = "omip",
                           file_splitting_interval = 360days,
                           async = false)

    model   = simulation.model
    ocean   = model.ocean
    sea_ice = model.sea_ice
    grid    = ocean.model.grid
    Nz      = size(grid, 3)

    T, S = ocean.model.tracers.T, ocean.model.tracers.S
    u, v, w = ocean.model.velocities
    η = ocean.model.free_surface.displacement

    hi = sea_ice.model.ice_thickness
    ℵi = sea_ice.model.ice_concentration
    hs = sea_ice.model.snow_thickness
    ui, vi = sea_ice.model.velocities

    tos = view(T, :, :, Nz)
    sos = view(S, :, :, Nz)

    surface_outputs = (tos       = tos,
                       sos       = sos,
                       zos       = η,
                       uos       = view(u, :, :, Nz),
                       vos       = view(v, :, :, Nz),
                       tossq     = tos * tos,
                       sossq     = sos * sos,
                       zossq     = Field(η * η),
                       mlotst    = MixedLayerDepthField(ocean.model.buoyancy, grid, ocean.model.tracers),
                       tauuo     = model.interfaces.net_fluxes.ocean.u,
                       tauvo     = model.interfaces.net_fluxes.ocean.v,
                       hfds      = model.interfaces.net_fluxes.ocean.T,
                       wfo       = model.interfaces.net_fluxes.ocean.S,
                       hfss      = model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat,
                       hfls      = model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat,
                       siconc    = ℵi,
                       sithick   = hi,
                       siu       = ui,
                       siv       = vi,
                       sitemptop = sea_ice.model.ice_thermodynamics.top_surface_temperature,
                       sisnthick = hs,
                       JTf       = NumericalEarth.Diagnostics.frazil_temperature_flux(model),
                       JTn       = NumericalEarth.Diagnostics.net_ocean_temperature_flux(model),
                       JTio      = NumericalEarth.Diagnostics.sea_ice_ocean_temperature_flux(model),
                       JTao      = NumericalEarth.Diagnostics.atmosphere_ocean_temperature_flux(model),
                       JSn       = NumericalEarth.Diagnostics.net_ocean_salinity_flux(model),
                       JSio      = NumericalEarth.Diagnostics.sea_ice_ocean_salinity_flux(model))

    # Each writer is handed `ocean.model`, but Oceananigans serializes `including` against the coupled model
    writer = JLD2Writer(ocean.model, surface_outputs;
                        including = Symbol[],
                        schedule = AveragedTimeInterval(surface_averaging_interval; stride = averaging_stride),
                        dir = output_dir,
                        filename = filename_prefix * "_surface",
                        file_splitting = TimeInterval(file_splitting_interval),
                        overwrite_files = true,
                        jld2_kw = Dict(:compress => ZstdFilter()))

    add_output_writer!(simulation, :surface, writer; async)

    bop = Oceananigans.Models.buoyancy_operation(ocean.model)

    field_outputs = (to     = T,
                     so     = S,
                     uo     = u,
                     vo     = v,
                     wo     = w,
                     bo     = bop,
                     uosq   = KernelFunctionOperation{Center, Center, Center}(uu_at_ccc, grid, u),
                     vosq   = KernelFunctionOperation{Center, Center, Center}(vv_at_ccc, grid, v),
                     uvol   = KernelFunctionOperation{Face,   Center, Center}(zonal_volume_flux,      grid, u),
                     vvol   = KernelFunctionOperation{Center, Face,   Center}(meridional_volume_flux, grid, v),
                     uvolto = KernelFunctionOperation{Face,   Center, Center}(zonal_tracer_transport,      grid, u, T),
                     vvolto = KernelFunctionOperation{Center, Face,   Center}(meridional_tracer_transport, grid, v, T),
                     uvolso = KernelFunctionOperation{Face,   Center, Center}(zonal_tracer_transport,      grid, u, S),
                     vvolso = KernelFunctionOperation{Center, Face,   Center}(meridional_tracer_transport, grid, v, S))

    writer = JLD2Writer(ocean.model, field_outputs;
                        including = Symbol[],
                        schedule = AveragedTimeInterval(field_averaging_interval; stride = averaging_stride),
                        dir = output_dir,
                        filename = filename_prefix * "_fields",
                        file_splitting = TimeInterval(file_splitting_interval),
                        overwrite_files = true,
                        jld2_kw = Dict(:compress => ZstdFilter()))

    add_output_writer!(simulation, :fields, writer; async)

    # Global volume means are stored as `Integral`s to be divided by `voco` offline, because `Average`
    # freezes its volume denominator at construction, which is wrong on a z-star grid. They are saved in
    # double precision: Float32 resolves only ~1e-7 of a global integral, coarser than the drifts they track
    average_outputs = (zosga = Average(η),
                       to_h  = Average(T,   dims = (1, 2)),
                       so_h  = Average(S,   dims = (1, 2)),
                       bo_h  = Average(bop, dims = (1, 2)),
                       voco  = Integral(KernelFunctionOperation{Center, Center, Center}(one_at_ccc, grid)),
                       soco  = Integral(S),
                       hoco  = Integral(T),
                       sivol = Integral(hi * ℵi),
                       snvol = Integral(hs * ℵi))

    writer = JLD2Writer(ocean.model, average_outputs;
                        including = Symbol[],
                        array_type = Array{Float64},
                        schedule = AveragedTimeInterval(field_mean_interval; stride = averaging_stride),
                        dir = output_dir,
                        filename = filename_prefix * "_averages",
                        file_splitting = TimeInterval(file_splitting_interval),
                        overwrite_files = true)

    add_output_writer!(simulation, :averages, writer; async)

    return nothing
end

#####
##### OAEMIP biogeochemical diagnostics
#####

@inline carbon_chemistry_diagnostic(i, j, k, grid, cc, DIC, Alk, fields, output) =
    @inbounds cc(; DIC = DIC[i, j, k],
                   Alk = Alk[i, j, k],
                   T   = fields.T[i, j, k],
                   S   = fields.S[i, j, k],
                   silicate  =  silicate_concentration(grid, i, j, k, fields),
                   phosphate = phosphate_concentration(grid, i, j, k, fields),
                   output)

@inline stored_pH_diagnostic(i, j, k, grid, pH) = @inbounds pH[i, j, 1]

@inline stored_pH_fCO₂_diagnostic(i, j, k, grid, cc, DIC, fields, pH) =
    @inbounds cc(; DIC = DIC[i, j, k],
                   T   = fields.T[i, j, k],
                   S   = fields.S[i, j, k],
                   pH  = pH[i, j, 1],
                   output = Val(:fCO₂))

surface_carbon_chemistry_diagnostics(::Nothing, grid, cc, DIC, Alk, fields) =
    (pH   = KernelFunctionOperation{Center, Center, Center}(carbon_chemistry_diagnostic, grid, cc, DIC, Alk, fields, Val(:pHᶠ)),
     fCO₂ = KernelFunctionOperation{Center, Center, Center}(carbon_chemistry_diagnostic, grid, cc, DIC, Alk, fields, Val(:fCO₂)))

surface_carbon_chemistry_diagnostics(stored_pH, grid, cc, DIC, Alk, fields) =
    (pH   = KernelFunctionOperation{Center, Center, Center}(stored_pH_diagnostic, grid, stored_pH),
     fCO₂ = KernelFunctionOperation{Center, Center, Center}(stored_pH_fCO₂_diagnostic, grid, cc, DIC, fields, stored_pH))

# the (DIC, Alk) tracer names of each carbonate system replicate, and the suffix its diagnostics are saved with
# (`CarbonateSystem(N)` names them `DIC1`/`Alk1`, ... when `N > 1`)
carbonate_replicates(::AbstractInorganicCarbon{1}) = ((DIC = :DIC, Alk = :Alk, suffix = ""), )
carbonate_replicates(::AbstractInorganicCarbon{N}) where N =
    ntuple(n -> (DIC = Symbol(:DIC, n), Alk = Symbol(:Alk, n), suffix = string(n)), N)

suffixed(names, suffix) = map(name -> Symbol(name, suffix), names)

function carbonate_replicate_diagnostics(replicate, ocean, interface, chemistry_fields)
    (; DIC, Alk, suffix) = replicate

    water_concentration = interface[DIC].water_concentration
    cc = water_concentration.carbon_chemistry
    stored_pH = hasproperty(water_concentration, :pH) ? water_concentration.pH : nothing

    DIC_field = ocean.tracers[DIC]
    Alk_field = ocean.tracers[Alk]

    pH, fCO₂ = surface_carbon_chemistry_diagnostics(stored_pH, ocean.grid, cc, DIC_field, Alk_field, chemistry_fields)

    surface_carbon_flux = DIC_field.boundary_conditions.top.condition.func.flux_field

    surface = NamedTuple{suffixed((:pH, :fCO₂, :surface_carbon_flux), suffix)}((pH, fCO₂, surface_carbon_flux))

    integrals = NamedTuple{suffixed((:surface_carbon_flux, :DIC, :Alk), suffix)}((Integral(surface_carbon_flux),
                                                                                   Integral(DIC_field),
                                                                                   Integral(Alk_field)))

    return surface, integrals
end

# community productivity in cells whose centre lies above `z_export` (the depth the POP export is saved at)
@inline function upper_productivity(i, j, k, grid, z_export, args...)
    above = znode(i, j, k, grid, Center(), Center(), Center()) > z_export
    return ifelse(above, community_productivity(i, j, k, grid, args...), zero(grid))
end

@inline oxygen_saturation_diagnostic(i, j, k, grid, oc, clock, fields) =
    surface_value(oc, i, j, grid, clock, fields)

"""
$(TYPEDSIGNATURES)

Add OAEMIP's biogeochemical output writers (carbonate system, O₂ saturation, productivity, POP export and
global carbon integrals) for a `NutrientsPlanktonDetritus` model with `ImplicitProductivity`, written in the
background with `async = true` (see [`add_output_writer!`](@ref)).
"""
function oaemip_diagnostics!(simulation;
                             surface_averaging_interval = 5days,
                             field_averaging_interval = 15days,
                             averaging_stride = 1,
                             output_dir = NumericalOcean.output_dir[],
                             filename_prefix = "omip",
                             file_splitting_interval = 360days,
                             async = false)

    ocean = simulation.model.ocean.model
    biogeochemistry = ocean.biogeochemistry
    grid = ocean.grid

    interface = simulation.model.interfaces.properties
    inorganic_carbon = biogeochemistry.underlying_biogeochemistry.inorganic_carbon

    # the potential temperature and practical salinity the gas exchange is computed with
    chemistry_fields = merge(fields(ocean), (T = simulation.model.interfaces.exchanger.ocean.state.T,
                                             S = interface.gas_exchange_state.practical_salinity))

    carbonate_diagnostics = map(r -> carbonate_replicate_diagnostics(r, ocean, interface, chemistry_fields), carbonate_replicates(inorganic_carbon))

    carbonate_surface   = mapreduce(first, merge, carbonate_diagnostics)
    carbonate_integrals = mapreduce(last,  merge, carbonate_diagnostics)

    oc = interface.O₂.air_concentration
    O₂_sat = KernelFunctionOperation{Center, Center, Center}(oxygen_saturation_diagnostic, grid, oc, ocean.clock, chemistry_fields)

    productivity_args = (biogeochemistry.underlying_biogeochemistry.plankton,
                         biogeochemistry.underlying_biogeochemistry,
                         fields(ocean),
                         biogeochemical_auxiliary_fields(biogeochemistry))

    PP = Field(Integral(KernelFunctionOperation{Center, Center, Center}(community_productivity, grid, productivity_args...),
                        dims = 3))

    POP_export = biogeochemistry.underlying_biogeochemistry.detritus.sinking.flux.POP

    z_export = on_architecture(CPU(), znodes(POP_export))
    k_100 = findmin(abs, z_export .+ 100)[2]

    # productivity above the export depth, so the export across it can be budgeted against it
    PP_100 = Field(Integral(KernelFunctionOperation{Center, Center, Center}(upper_productivity, grid,
                                                                            convert(eltype(grid), z_export[k_100]),
                                                                            productivity_args...),
                            dims = 3))

    output_fields      = ocean.tracers[(required_biogeochemical_tracers(inorganic_carbon)..., :PO₄, :Fe, :O₂, :DOP)]
    output_surface     = merge(carbonate_surface, (; O₂_sat))
    vertical_integrals = (; PP, PP_100)
    export_flux        = (; POP_export)
    spatial_integrals  = carbonate_integrals

    writer = JLD2Writer(ocean, output_fields;
                        schedule = AveragedTimeInterval(field_averaging_interval; stride = averaging_stride),
                        dir = output_dir,
                        filename = filename_prefix * "_bgc_fields",
                        file_splitting = TimeInterval(file_splitting_interval),
                        overwrite_files = true,
                        jld2_kw = Dict(:compress => ZstdFilter()))

    add_output_writer!(simulation, :bgc_fields, writer; async)

    writer = JLD2Writer(ocean, output_surface;
                        indices = (:, :, grid.Nz),
                        schedule = AveragedTimeInterval(surface_averaging_interval; stride = averaging_stride),
                        dir = output_dir,
                        filename = filename_prefix * "_bgc_surfaces",
                        file_splitting = TimeInterval(file_splitting_interval),
                        overwrite_files = true,
                        jld2_kw = Dict(:compress => ZstdFilter()))

    add_output_writer!(simulation, :bgc_surfaces, writer; async)

    writer = JLD2Writer(ocean, vertical_integrals;
                        indices = (:, :, 1),
                        schedule = AveragedTimeInterval(field_averaging_interval; stride = averaging_stride),
                        dir = output_dir,
                        filename = filename_prefix * "_bgc_vertical_integrals",
                        file_splitting = TimeInterval(file_splitting_interval),
                        overwrite_files = true,
                        jld2_kw = Dict(:compress => ZstdFilter()))

    add_output_writer!(simulation, :bgc_vertical_integrals, writer; async)

    writer = JLD2Writer(ocean, export_flux;
                        indices = (:, :, k_100),
                        schedule = AveragedTimeInterval(surface_averaging_interval; stride = averaging_stride),
                        dir = output_dir,
                        filename = filename_prefix * "_bgc_exports",
                        file_splitting = TimeInterval(file_splitting_interval),
                        overwrite_files = true,
                        jld2_kw = Dict(:compress => ZstdFilter()))

    add_output_writer!(simulation, :bgc_exports, writer; async)

    # global inventories are ~3e21 mmol, where Float32 only resolves ~3e14 mmol (0.3 Tmol): too coarse to
    # see a perturbation in the difference between replicates, so these are saved in double precision
    writer = JLD2Writer(ocean, spatial_integrals;
                        indices = (1, 1, 1),
                        array_type = Array{Float64},
                        schedule = AveragedTimeInterval(field_averaging_interval; stride = averaging_stride),
                        dir = output_dir,
                        filename = filename_prefix * "_bgc_integrals",
                        file_splitting = TimeInterval(file_splitting_interval),
                        overwrite_files = true,
                        jld2_kw = Dict(:compress => ZstdFilter()))

    add_output_writer!(simulation, :bgc_spatial_integrals, writer; async)

    return nothing
end

#####
##### Checkpointing
#####

"""
$(TYPEDSIGNATURES)

Add a checkpointer of the coupled model, for `run!(simulation; pickup = true)`. Call it after the diagnostics:
writers run in the order they were added, so at an iteration where an output and a checkpoint coincide the
output is written first.
"""
function checkpointer!(simulation;
                       checkpoint_interval = 360days,
                       output_dir = NumericalOcean.output_dir[],
                       filename_prefix = "omip")

    simulation.output_writers[:checkpointer] = Checkpointer(simulation.model;
                                                            schedule = TimeInterval(checkpoint_interval),
                                                            dir      = output_dir,
                                                            prefix   = filename_prefix * "_checkpoint",
                                                            cleanup  = false,
                                                            verbose  = true)

    return nothing
end
