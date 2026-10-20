using Oceananigans.Units: minutes, days

function default_maximum_search_radius(grid)
    Nx, Ny, _ = size(grid)
    return max(5, ceil(Int, 3 / ((360 / Nx + 180 / Ny) / 2)))
end

default_river_spread_radius(grid)        = 1.2
default_river_spread_radius(::ORCA1GRID) = nothing

default_river_spread_cells(grid)        = nothing
default_river_spread_cells(::ORCA1GRID) = 8

default_κ_skew(grid)        = nothing
default_κ_skew(::ORCA1GRID) = 800

default_κ_symmetric(grid)        = nothing
default_κ_symmetric(::ORCA1GRID) = 800

default_momentum_advection_order(grid)        = nothing
default_momentum_advection_order(::ORCA1GRID) = 5

default_biharmonic_timescale(grid)             = 50days
default_biharmonic_timescale(::ORCATPivotGRID) = nothing

default_barotropic_substeps(grid)             = 100
default_barotropic_substeps(::ORCA1GRID)      = 300
default_barotropic_substeps(::ORCATPivotGRID) = 200

default_Δt(grid)                 = 30minutes
default_Δt(::ORCA1GRID)          = 90minutes
default_Δt(grid::ORCATPivotGRID) = is_orca_twelfth(grid) ? 6minutes : 20minutes

# T and S at 7th order and the extra and biogeochemical tracers at 5th, as in OMIPSimulations. Every tracer is named,
# because `ocean_simulation` gives any it is not told about an explicit `WENO(order = 5)`.
# Next to land the reconstruction ends in first-order upwind (OMIPSimulations' `boundary_scheme = :upwind`) rather than
# the default `Centered(order = 2)`, which grows spurious extrema in coastal corner and dead-end cells (e.g. S ≈ 43 and
# PO₄ ≈ -2.3 mmol m⁻³ next to the Bosphorus).
function tracer_advection_schemes(tracers, biogeochemistry, time_discretization)
    temperature_salinity = WENO(; order = 7, time_discretization, boundary_scheme = UpwindBiased(order = 1))
    other = WENO(; order = 5, time_discretization, boundary_scheme = UpwindBiased(order = 1))
    other_names = filter(name -> name ∉ (:T, :S), (something(tracers, ())..., required_biogeochemical_tracers(biogeochemistry)...))
    return merge((T = temperature_salinity, S = temperature_salinity),
                 NamedTuple(name => other for name in other_names))
end
