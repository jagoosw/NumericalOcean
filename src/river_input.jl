# RivR2O river loads, copied from OAEMIP/src/river_input.jl

using NCDatasets
using JLD2
using Statistics
using Oceananigans: Forcing, CPU, Field, KernelFunctionOperation
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Grids: Face
using Oceananigans.Grids: λnode, φnode
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.OutputReaders: interpolating_time_indices, Cyclical
using Oceananigans.Operators: volume, Az
using Oceananigans.Utils: on_architecture
using NumericalEarth.DataWrangling.JRA55: JRA55PrescribedLand
using Oceananigans.Units: days
using Adapt: Adapt, adapt


# we only have todo this once so its going to be ugly
function map_river_load(grid, lon, lat; land = nothing)
    grid = on_architecture(CPU(), grid)
    isnothing(land) && (land = JRA55PrescribedLand(grid)) 

    river_target_i = on_architecture(CPU(), land.river_routing.rivers.target_i)
    river_target_j = on_architecture(CPU(), land.river_routing.rivers.target_j)

    river_target_lon = map(n->λnode(river_target_i[n], river_target_j[n], 1, grid, C, C, C), 1:length(river_target_i))
    river_target_lat = map(n->φnode(river_target_i[n], river_target_j[n], 1, grid, C, C, C), 1:length(river_target_i))

    distance = @. sqrt(abs(lon - river_target_lon')^2 + abs(lat - river_target_lat')^2)

    nearest_target = argmin.(eachrow(distance))
    target_indices = sortperm(nearest_target)

    unique_targets = sort(unique(nearest_target))
    source_indices = ones(Int, length(unique_targets))

    for n in 2:length(unique_targets)
        source_indices[n] = source_indices[n-1] + count(m -> m == unique_targets[n-1], nearest_target)
    end

    # one past the last mouth closes the final target's range, and again for the empty range of cells without a river
    N = length(nearest_target)
    source_indices = [source_indices..., N + 1, N + 1]

    target_map = Field{Center, Center, Nothing}(grid, Int)

    for (n, idx) in enumerate(unique_targets)
        i = river_target_i[idx]
        j = river_target_j[idx]
        @inbounds target_map[i, j, 1] = n
    end

    interior(target_map)[interior(target_map) .== 0] .= length(unique_targets) + 1

    return target_map, target_indices, source_indices
end

function read_file(file_path)
    ds = NCDataset(file_path)
    lon = Array(ds["lon"])
    lat = Array(ds["lat"])
    variables = tuple([v for v in keys(ds) if ndims(ds[v]) == 2]...)
    FT = eltype(lon)

    # mouths are the cells where any load is nonzero: RivR2O puts DIP on cells without DIC (37 % of it, e.g. a
    # second Amazon cell holding more DIP than the Amazon's DIC mouth), which a DIC-only mask drops
    loads = map(name -> coalesce.(ds[name][:, :], zero(FT)), variables)
    idx = findall(.|(map(load -> .!iszero.(load), loads)...))
    mouth_lon = lon[getindex.(idx, 1)]
    mouth_lat = lat[getindex.(idx, 2)]

    data = NamedTuple{Symbol.(variables)}(map(load -> convert.(FT, load[idx]), loads))

    return data, mouth_lon, mouth_lat, idx
end

function construct_timeseries(dir = "data/r2o_river_inputs_1901_2024")
    files = filter(endswith(".nc"), readdir(dir))
    years = sort(parse.(Int, [match(r"\d{4}", f).match for f in files]))

    data1, mouth_lon, mouth_lat, mouth_map = read_file(joinpath(dir, files[1]))

    variable_names = keys(data1)
    river_count = length(mouth_lon)
    FT = eltype(data1.DIC)

    data = NamedTuple{variable_names}(map(_->zeros(FT, river_count, length(years)), variable_names))

    for (n, yr) in enumerate(years)
        ds = NCDataset(joinpath(dir, "rivr2o_riverinputs_$yr.nc"))

        for var in variable_names
            vals = ds[var][mouth_map]
            vals[ismissing.(vals)] .= zero(FT) # not sure this is the correct interpretation of this data yet
            vals[vals .< 0] .= zero(FT)
            getproperty(data, var)[:, n] .= vals
        end
    end

    return data, years, mouth_lon, mouth_lat
end

const yr = 365days

grid_name(grid) = "$(size(grid, 1))_$(size(grid, 2))" 
grid_name(grid::ImmersedBoundaryGrid) = "immersed_"*grid_name(grid.underlying_grid)

const DEFAULT_RIVER_TRACERS = (:DIC, :Alk, :PO₄, :DOP)

# RivR2O loads are in 10⁶ g yr⁻¹, converted to mmol s⁻¹ with 10⁶ g × 10³ mmol/mol / molar mass. Following the
# RivR2O guidance (zenodo.org/records/14889524):
# PO₄ from DIP, which already includes the labile DOP (terrestrial DOC is C:P ≈ 2583:1, so carries no P of its own)
# DIC from DIC + labile DOC + semi-labile DOC (the latter should decay over ~1.5 yr; we have no tracer for it)
# DOP from the recycled POC, as marine organic matter at the model's C:P of 117
# Alk from DIC (weathering bicarbonate, 1:1)
const DEFAULT_RIVER_INPUTS = (DIC = (DIC = 1e9/12/yr, DOC_l = 1e9/12/yr, DOC_sl = 1e9/12/yr),
                              Alk = (DIC = 1e9/12/yr, ),
                              PO₄ = (DIP = 1e9/30.97/yr, ),
                              DOP = (POC = 1e9/12/117/yr, ))

"""
    RivR2OForcing(grid; forced_tracers = (:DIC, :Alk, :PO₄, :DOP), kwargs...)

The RivR2O river loads as a 3D `Forcing` per tracer, nonzero only in the surface layer but evaluated
at every level. Prefer [`RivR2OSurfaceFlux`](@ref), which applies the same input as a top boundary
flux; this is kept for comparison. Carbonate replicates (`DIC1`, `Alk2`, ...) in `forced_tracers` get
the `DIC`/`Alk` input. See `river_load_data` for the keyword arguments.
"""
function RivR2OForcing(grid; forced_tracers = DEFAULT_RIVER_TRACERS, kwargs...)
    input_tracers = replicated_tracers(forced_tracers)
    data, times, source_map, offsets = river_load_data(grid; forced_tracers = input_tracers, kwargs...)
    forcing = build_river_forcing(architecture(grid), input_tracers, data, times, source_map, offsets)
    return replicate_river_inputs(forcing, forced_tracers)
end

"""
    RivR2OSurfaceFlux(grid; forced_tracers = (:DIC, :Alk, :PO₄, :DOP), kwargs...)

The RivR2O river loads as a top `FluxBoundaryCondition` per tracer, `-load / area` at the river mouth
cells (negative is into the ocean), for the `additional_surface_fluxes` of the ocean. With
`climatology = true` the flux is constant and is computed once into a 2D `Field`; otherwise the loads
are time interpolated at the surface each evaluation. Identical in the top cell to [`RivR2OForcing`](@ref).
Carbonate replicates (`DIC1`, `Alk2`, ...) in `forced_tracers` get the `DIC`/`Alk` input.
See `river_load_data` for the keyword arguments.
"""
function RivR2OSurfaceFlux(grid; forced_tracers = DEFAULT_RIVER_TRACERS, kwargs...)
    input_tracers = replicated_tracers(forced_tracers)
    data, times, source_map, offsets = river_load_data(grid; forced_tracers = input_tracers, kwargs...)
    loads = build_river_loads(architecture(grid), input_tracers, data, times, source_map, offsets)
    fluxes = map(load -> river_surface_flux(grid, load), loads)
    return replicate_river_inputs(fluxes, forced_tracers)
end

# `CarbonateSystem(N)` names its tracers `DIC1`, `Alk1`, ... when `N > 1`, and each replicate takes the river input
# (and initial condition) of the tracer it replicates. The inputs are only read, so the replicates share them.
function replicated_tracer(name)
    base = Symbol(rstrip(isdigit, string(name)))
    return base ∈ (:DIC, :Alk) ? base : name
end

replicated_tracers(forced_tracers) = Tuple(unique(map(replicated_tracer, forced_tracers)))

replicate_river_inputs(inputs, forced_tracers) =
    NamedTuple{Tuple(forced_tracers)}(map(name -> inputs[replicated_tracer(name)], Tuple(forced_tracers)))

function river_load_data(grid;
                         dir = bgc_dir[],
                         full_dir = joinpath(dir, "r2o_river_inputs_1901_2024"),
                         forced_tracers = DEFAULT_RIVER_TRACERS,
                         input_tracers = DEFAULT_RIVER_INPUTS,
                         land = nothing,
                         climatology = false,
                         start_year = 1958,
                         end_year = 2017,
                         load_cache = true)

    climatology_suffix = climatology ? "_climatology" : ""

    # only plain CPU arrays are cached - JLD2 can't reliably round trip `RiverLoad`/`Forcing`/`Field`/GPU arrays
    # (it silently returns reconstructed types which then fail to compile in kernels)
    save_name = joinpath(full_dir, "mapped_$(grid_name(grid))_$(start_year)_$(end_year)$(climatology_suffix)_v6.jld2")

    if isfile(save_name) & load_cache
        return jldopen(save_name) do f
            f["data"], f["times"], f["source_map"], f["offsets"]
        end
    end

    if isnothing(land)
        land = JRA55PrescribedLand(grid)
    end

    data, years, lon, lat = construct_timeseries(full_dir)
    times = @. (years - start_year + 0.5) * yr
    time_inds = (times .> 0) .& (times .< (end_year + 1 - start_year) * yr)

    source_map, order, offsets = map_river_load(grid, lon, lat; land)

    data = NamedTuple{forced_tracers}(
        map(dn -> total_load(data, keys(input_tracers[dn]), values(input_tracers[dn]))[order, time_inds],
            forced_tracers)
    )

    if climatology
        data = NamedTuple{forced_tracers}(map(d -> mean(d, dims = 2)[:, 1], values(data)))
        times = nothing
    else
        times = times[time_inds]
    end

    source_map = source_map.data # CPU OffsetArray, indexed as map[i, j, 1] in the kernel

    jldsave(save_name; data, times, source_map, offsets)

    return data, times, source_map, offsets
end

build_river_forcing(arch, forced_tracers, args...) =
    map(load -> Forcing(load, discrete_form = true), build_river_loads(arch, forced_tracers, args...))

function build_river_loads(arch, forced_tracers, data, times, source_map, offsets)
    if isnothing(times)
        time_interpolation = nothing
    else
        time_interpolation = Cyclical(times[end] - times[1] + (times[end] - times[end-1]))
    end

    source_map = on_architecture(arch, source_map)
    times = on_architecture(arch, times)
    offsets = on_architecture(arch, offsets)

    return NamedTuple{forced_tracers}(map(tn -> RiverLoad(on_architecture(arch, data[tn]), time_interpolation, times, source_map, offsets), forced_tracers))
end

total_load(data, names, scales) =
    sum(getproperty(data, name) .* scale for (name, scale) in zip(names, scales))

const C = Center()
const CCC = (C, C, C)

struct RiverLoad{L, TI, T, M, O} <: Function # DiscreteForcing only calls `func` if it is a `Function`
         source_load :: L
  time_interpolation :: TI
               times :: T
                 map :: M
       target_offset :: O
end

Adapt.adapt_structure(to, rl::RiverLoad) =
    RiverLoad(adapt(to, rl.source_load),
              adapt(to, rl.time_interpolation),
              adapt(to, rl.times),
              adapt(to, rl.map),
              adapt(to, rl.target_offset))

# the total load of the river mouths mapped to column (i, j)
@inline function accumulated_river_load(ri::RiverLoad, i, j, grid, clock)
    accumulated_load = zero(eltype(grid))

    @inbounds begin
        c = ri.map[i, j, 1]
        load = ri.source_load
        f, n₁, n₂ = interpolating_time_indices(ri.time_interpolation, ri.times, clock.time)

        for m in ri.target_offset[c]:(ri.target_offset[c+1]-1)
            accumulated_load += load[m, n₁] * (1 - f) + load[m, n₂] * f
        end
    end

    return accumulated_load
end

@inline function accumulated_river_load(ri::RiverLoad{<:Any, <:Any, Nothing}, i, j, grid, clock)
    accumulated_load = zero(eltype(grid))

    @inbounds begin
        c = ri.map[i, j, 1]
        load = ri.source_load

        for m in ri.target_offset[c]:(ri.target_offset[c+1]-1)
            accumulated_load += load[m]
        end
    end

    return accumulated_load
end

# as a `Forcing`: the load into the top cell
@inline (ri::RiverLoad)(i, j, k, grid, clock, fields) =
    accumulated_river_load(ri, i, j, grid, clock) / volume(i, j, k, grid, CCC...) * (k == grid.Nz)

# as a top boundary flux (positive upwards, so negative into the ocean): with the boundary condition's
# `Az / volume` this is the same tendency as the `Forcing`
@inline (ri::RiverLoad)(i, j, grid, clock, fields) =
    - accumulated_river_load(ri, i, j, grid, clock) / Az(i, j, grid.Nz+1, grid, Center(), Center(), Face())

# time constant loads only need evaluating once
river_surface_flux(grid, load::RiverLoad{<:Any, <:Any, Nothing}) =
    FluxBoundaryCondition(Field(KernelFunctionOperation{Center, Center, Nothing}(river_surface_flux_kernel, grid, load)))

river_surface_flux(grid, load) = FluxBoundaryCondition(load; discrete_form = true)

@inline river_surface_flux_kernel(i, j, k, grid, load) = load(i, j, grid, nothing, nothing)

