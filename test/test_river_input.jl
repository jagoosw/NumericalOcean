# RivR2O river loads: the RiverLoad kernels on hand-built arrays, and the loading/mapping pipeline on a
# tiny fake RivR2O dataset with a fake `land` (only `land.river_routing.rivers.target_i/j` is used).

using Oceananigans.OutputReaders: Cyclical
using Oceananigans.Forcings: DiscreteForcing
using Oceananigans.BoundaryConditions: BoundaryCondition, DiscreteBoundaryFunction

const yr = 365days

@testset "DEFAULT_RIVER_INPUTS / total_load" begin
    @test keys(GO.DEFAULT_RIVER_INPUTS) == GO.DEFAULT_RIVER_TRACERS == (:DIC, :Alk, :PO₄, :DOP)
    @test GO.DEFAULT_RIVER_INPUTS.DIC.DIC ≈ 1e9 / 12 / yr       # 10⁶ g C/yr -> mmol C/s
    @test GO.DEFAULT_RIVER_INPUTS.PO₄.DIP ≈ 1e9 / 30.97 / yr
    @test keys(GO.DEFAULT_RIVER_INPUTS.PO₄) == (:DIP,) # labile DOP is already in DIP
    @test GO.DEFAULT_RIVER_INPUTS.DOP.POC ≈ GO.DEFAULT_RIVER_INPUTS.DIC.DIC / 117 # marine C:P
    @test keys(GO.DEFAULT_RIVER_INPUTS.DIC) == (:DIC, :DOC_l, :DOC_sl)

    data = (a = [1.0 2.0; 3.0 4.0], b = [10.0 20.0; 30.0 40.0])
    @test GO.total_load(data, (:a,), (2.0,)) == 2 .* data.a
    @test GO.total_load(data, (:a, :b), (2.0, 0.5)) == 2 .* data.a .+ 0.5 .* data.b
    @test GO.total_load(merge(data, (; c = data.a)), (:a, :b, :c), (2.0, 0.5, 1.0)) == 3 .* data.a .+ 0.5 .* data.b
end

# Columns: (1, 1) -> target 1 (mouths 1:2), (2, 2) -> target 2 (mouth 3), everything else -> 3 (empty range)
const RIVER_GRID = RectilinearGrid(size = (3, 3, 2), extent = (3, 3, 4)) # cells 1 × 1 × 2 m
function test_river_map()
    map = fill(3, 3, 3, 1)
    map[1, 1, 1] = 1
    map[2, 2, 1] = 2
    return map
end
const RIVER_OFFSETS = [1, 3, 4, 4]
const RIVER_TIMES = [0.0, 10.0, 20.0]
const RIVER_LOAD = [1.0 2.0 3.0; 10.0 20.0 30.0; 100.0 200.0 300.0] # mouths × times

clock_at(t) = (; time = t)

@testset "RiverLoad: climatology (time constant)" begin
    load = GO.build_river_loads(CPU(), (:X,), (; X = [1.0, 10.0, 100.0]), nothing, test_river_map(), RIVER_OFFSETS).X
    @test load isa GO.RiverLoad
    @test isnothing(load.time_interpolation)
    @test isnothing(load.times)

    @test GO.accumulated_river_load(load, 1, 1, RIVER_GRID, nothing) == 11  # mouths 1 and 2
    @test GO.accumulated_river_load(load, 2, 2, RIVER_GRID, nothing) == 100 # mouth 3
    @test GO.accumulated_river_load(load, 3, 3, RIVER_GRID, nothing) == 0   # no river
    @test GO.accumulated_river_load(load, 1, 2, RIVER_GRID, nothing) == 0

    # as a Forcing: load / volume, only in the top cell
    @test load(1, 1, 2, RIVER_GRID, nothing, nothing) ≈ 11 / 2
    @test load(1, 1, 1, RIVER_GRID, nothing, nothing) == 0
    # as a flux: -load / area (negative into the ocean)
    @test load(1, 1, RIVER_GRID, nothing, nothing) ≈ -11
    @test load(2, 2, RIVER_GRID, nothing, nothing) ≈ -100

    # the time constant flux is precomputed into a 2D Field
    bc = GO.river_surface_flux(RIVER_GRID, load)
    @test bc isa BoundaryCondition
    flux = interior(bc.condition, :, :, 1)
    @test flux ≈ [-11 0 0; 0 -100 0; 0 0 0]

    # flux × area == forcing × volume: the same tracer tendency in the top cell
    Az, V = 1.0, 2.0
    @test -load(1, 1, RIVER_GRID, nothing, nothing) * Az ≈ load(1, 1, 2, RIVER_GRID, nothing, nothing) * V
end

@testset "RiverLoad: time varying" begin
    loads = GO.build_river_loads(CPU(), (:X, :Y), (; X = RIVER_LOAD, Y = 2 .* RIVER_LOAD),
                                 RIVER_TIMES, test_river_map(), RIVER_OFFSETS)
    load = loads.X
    @test keys(loads) == (:X, :Y)

    # cyclic over the record plus one spacing
    @test load.time_interpolation isa Cyclical
    @test load.time_interpolation.period == 30

    acc(i, j, t) = GO.accumulated_river_load(load, i, j, RIVER_GRID, clock_at(t))

    @test acc(1, 1, 0.0) ≈ 11
    @test acc(2, 2, 0.0) ≈ 100
    @test acc(3, 3, 7.0) == 0
    @test acc(1, 1, 5.0) ≈ 16.5   # midpoint of 11 and 22
    @test acc(1, 1, 25.0) ≈ 22    # midway between t = 20 (33) and the wrap to t = 30 ≡ 0 (11)
    @test GO.accumulated_river_load(loads.Y, 1, 1, RIVER_GRID, clock_at(5.0)) ≈ 33

    # Off-midpoint times expose the interpolation weights: they are swapped in
    # `accumulated_river_load` (src/river_input.jl), so these are currently wrong.
    @test acc(1, 1, 10.0) ≈ 22  # exactly on the second record
    @test acc(1, 1, 2.5) ≈ 0.75 * 11 + 0.25 * 22
    @test acc(2, 2, 20.0) ≈ 300

    @test load(1, 1, 2, RIVER_GRID, clock_at(5.0), nothing) ≈ 16.5 / 2
    @test load(1, 1, 1, RIVER_GRID, clock_at(5.0), nothing) == 0
    @test load(1, 1, RIVER_GRID, clock_at(5.0), nothing) ≈ -16.5

    bc = GO.river_surface_flux(RIVER_GRID, load)
    @test bc.condition isa DiscreteBoundaryFunction

    forcing = GO.build_river_forcing(CPU(), (:X,), (; X = RIVER_LOAD), RIVER_TIMES, test_river_map(), RIVER_OFFSETS)
    @test forcing.X isa DiscreteForcing
    @test forcing.X.func isa GO.RiverLoad
end

@testset "RiverLoad in a model" begin
    # the flux form runs through a real model's tracer tendency
    load = GO.build_river_loads(CPU(), (:c,), (; c = [1.0, 10.0, 100.0]), nothing, test_river_map(), RIVER_OFFSETS).c
    c_bcs = FieldBoundaryConditions(top = GO.river_surface_flux(RIVER_GRID, load))
    model = HydrostaticFreeSurfaceModel(RIVER_GRID; tracers = :c, buoyancy = nothing,
                                        boundary_conditions = (; c = c_bcs))
    time_step!(model, 1.0)
    c = interior(model.tracers.c)
    @test c[1, 1, 2] ≈ 11 / 2
    @test c[2, 2, 2] ≈ 100 / 2
    @test c[1, 1, 1] == 0
    @test sum(c) * 2 ≈ 111 # mass conserved: Σ c V = Σ load × Δt
end

# A tiny fake RivR2O dataset: mouths at (lon, lat) grid points with load 1, 10, 100, 1000 × year index
function write_fake_rivr2o(dir, years)
    lon = [5.0, 6.0, 26.0, 34.0, 20.0]
    lat = [-15.0, -14.0, -5.0, 16.0]
    for (n, year) in enumerate(years)
        NCDataset(joinpath(dir, "rivr2o_riverinputs_$year.nc"), "c") do ds
            defDim(ds, "lon", length(lon)); defDim(ds, "lat", length(lat))
            defVar(ds, "lon", lon, ("lon",)); defVar(ds, "lat", lat, ("lat",))
            for var in ("DIC", "DOC_l", "DIP", "DOC_sl", "POC")
                A = zeros(length(lon), length(lat))
                A[1, 1] = n; A[2, 2] = 10n; A[3, 3] = 100n; A[4, 4] = 1000n
                defVar(ds, var, A, ("lon", "lat"))
            end
        end
    end
end

# targets at the centres of cells (1, 1) = (5°, -15°), (3, 2) = (25°, -5°) and (4, 4) = (35°, 15°)
const RIVER_LL_GRID = LatitudeLongitudeGrid(size = (4, 4, 2), longitude = (0, 40), latitude = (-20, 20), z = (-100, 0))
const FAKE_LAND = (; river_routing = (; rivers = (; target_i = [1, 3, 4], target_j = [1, 2, 4])))

@testset "map_river_load" begin
    lon = [5.0, 6.0, 26.0, 34.0]
    lat = [-15.0, -14.0, -5.0, 16.0]
    target_map, order, offsets = GO.map_river_load(RIVER_LL_GRID, lon, lat; land = FAKE_LAND)

    @test order == [1, 2, 3, 4]
    @test target_map[1, 1, 1] == 1
    @test target_map[3, 2, 1] == 2
    @test offsets[1:3] == [1, 3, 4]

    columns(c) = offsets[c]:(offsets[c + 1] - 1)
    @test columns(target_map[1, 1, 1]) == 1:2 # mouths 1 and 2 -> target 1
    @test columns(target_map[3, 2, 1]) == 3:3 # mouth 3 -> target 2
    @test isempty(columns(target_map[2, 3, 1])) # cells without a river get nothing

    # Mouth 4 -> target 3 at (4, 4), the last target. `map_river_load` closes the offsets with a
    # repeat of the last start instead of one-past-the-end, so the last target gets an empty range
    # (and shares its index with the cells that have no river): its rivers are dropped.
    @test columns(target_map[4, 4, 1]) == 4:4
    @test target_map[4, 4, 1] != target_map[2, 3, 1]

    # unsorted mouths are reordered by target
    target_map, order, offsets = GO.map_river_load(RIVER_LL_GRID, [26.0, 5.0, 34.0, 6.0], [-5.0, -15.0, 16.0, -14.0];
                                                   land = FAKE_LAND)
    @test order == [2, 4, 1, 3]
    @test offsets[1:3] == [1, 3, 4]
end

@testset "river_load_data / RivR2OSurfaceFlux (fake RivR2O data)" begin
    dir = mktempdir()
    full_dir = joinpath(dir, "r2o_river_inputs_1901_2024")
    mkpath(full_dir)
    write_fake_rivr2o(full_dir, 1985:1987)

    kw = (; full_dir, land = FAKE_LAND, start_year = 1985, end_year = 1987, load_cache = false)
    Alk_scale = GO.DEFAULT_RIVER_INPUTS.Alk.DIC

    data, times, source_map, offsets = GO.river_load_data(RIVER_LL_GRID; forced_tracers = (:Alk, :DIC), kw...)
    @test keys(data) == (:Alk, :DIC)
    @test times ≈ [0.5, 1.5, 2.5] .* yr
    @test size(data.Alk) == (4, 3) # mouths × years
    @test data.Alk[:, 1] ≈ [1, 10, 100, 1000] .* Alk_scale
    @test data.DIC[:, 1] ≈ 3 .* data.Alk[:, 1] # DIC + DOC_l + DOC_sl, all at 1e9/12/yr
    @test source_map[1, 1, 1] == 1

    # the loads grow with the year index
    @test data.Alk[:, 2] ≈ 2 .* data.Alk[:, 1]
    @test data.Alk[:, 3] ≈ 3 .* data.Alk[:, 1]

    # a cache file is written
    @test any(startswith("mapped_4_4_1985_1987"), readdir(full_dir))
    @test isfile(joinpath(full_dir, "mapped_4_4_1985_1987_v6.jld2"))

    clim, clim_times, _, _ = GO.river_load_data(RIVER_LL_GRID; forced_tracers = (:Alk,), climatology = true, kw...)
    @test isnothing(clim_times)
    @test clim.Alk isa AbstractVector
    @test length(clim.Alk) == 4
    @test isfile(joinpath(full_dir, "mapped_4_4_1985_1987_climatology_v6.jld2"))
    @test clim.Alk ≈ 2 .* [1, 10, 100, 1000] .* Alk_scale # mean over 1985-1987 (see above)

    # load_cache = true round-trips through the JLD2 cache
    cached = GO.river_load_data(RIVER_LL_GRID; forced_tracers = (:Alk,), climatology = true,
                                kw..., load_cache = true)
    @test cached[1].Alk == clim.Alk
    @test isnothing(cached[2])
    @test cached[4] == offsets

    # end-to-end flux: -load / Az at the river mouth cells
    fluxes = GO.RivR2OSurfaceFlux(RIVER_LL_GRID; forced_tracers = (:Alk,), climatology = true, kw...)
    @test keys(fluxes) == (:Alk,)
    flux = fluxes.Alk.condition
    Az(i, j) = Oceananigans.Operators.Az(i, j, 2, RIVER_LL_GRID, Center(), Center(), Center())
    @test flux[1, 1, 1] ≈ -(clim.Alk[1] + clim.Alk[2]) / Az(1, 1)
    @test flux[3, 2, 1] ≈ -clim.Alk[3] / Az(3, 2)
    @test flux[2, 3, 1] == 0
    @test all(interior(flux) .<= 0)

    time_varying = GO.RivR2OSurfaceFlux(RIVER_LL_GRID; forced_tracers = (:Alk,), kw...)
    @test time_varying.Alk.condition isa DiscreteBoundaryFunction

    forcing = GO.RivR2OForcing(RIVER_LL_GRID; forced_tracers = (:Alk,), climatology = true, kw...)
    @test forcing.Alk isa DiscreteForcing

    # carbonate replicates share the input of the tracer they replicate
    replicated = GO.RivR2OSurfaceFlux(RIVER_LL_GRID; forced_tracers = (:Alk1, :Alk2), climatology = true, kw...)
    @test keys(replicated) == (:Alk1, :Alk2)
    @test replicated.Alk1 === replicated.Alk2
    @test interior(replicated.Alk2.condition) == interior(flux)

    replicated_forcing = GO.RivR2OForcing(RIVER_LL_GRID; forced_tracers = (:Alk1, :Alk2), climatology = true, kw...)
    @test keys(replicated_forcing) == (:Alk1, :Alk2)
    @test replicated_forcing.Alk2 isa DiscreteForcing
end

@testset "carbonate replicates" begin
    @test GO.replicated_tracer(:DIC2) == :DIC
    @test GO.replicated_tracer(:Alk1) == :Alk
    @test GO.replicated_tracer(:Alk12) == :Alk
    @test GO.replicated_tracer(:DIC) == :DIC
    @test GO.replicated_tracer(:PO₄) == :PO₄
    @test GO.replicated_tracer(:O₂) == :O₂
    @test GO.replicated_tracer(:T2) == :T2 # only DIC and Alk are replicated

    @test GO.replicated_tracers((:DIC1, :Alk1, :DIC2, :Alk2, :PO₄)) == (:DIC, :Alk, :PO₄)

    inputs = (DIC = 1, Alk = 2, PO₄ = 3)
    @test GO.replicate_river_inputs(inputs, (:DIC1, :Alk1, :DIC2, :Alk2, :PO₄)) ==
          (DIC1 = 1, Alk1 = 2, DIC2 = 1, Alk2 = 2, PO₄ = 3)
    @test GO.replicate_river_inputs(inputs, (:DIC, :Alk)) == (DIC = 1, Alk = 2)

    # the tracers `forced_ocean_simulation` gives river inputs and initial conditions to
    bgc = MITgcmDIC(RectilinearGrid(size = (1, 1, 2), extent = (1, 1, 2)); inorganic_carbon = CarbonateSystem(2))
    tracers = Oceananigans.Biogeochemistry.required_biogeochemical_tracers(bgc)
    @test filter(t -> haskey(GO.DEFAULT_RIVER_INPUTS, GO.replicated_tracer(t)), tracers) ⊇ (:DIC1, :Alk1, :DIC2, :Alk2, :PO₄, :DOP)
    @test all(t -> haskey(GO.BGC_INIT_TRACERS, GO.replicated_tracer(t)), (:DIC1, :Alk1, :DIC2, :Alk2))
end
