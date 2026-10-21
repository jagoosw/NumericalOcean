# preset grids, but any other grid can be dropped in
using KernelAbstractions: @kernel, @index
using OMIPSimulations: omip_vertical_discretization
using Oceananigans.Architectures: architecture
using Oceananigans.Grids: OrthogonalSphericalShellGrid, RightFaceFolded, RightCenterFolded
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.OrthogonalSphericalShellGrids: Tripolar
using Oceananigans.BoundaryConditions: FPivot, TPivot, fill_halo_regions!
using Oceananigans.DistributedComputations: Distributed, global_size

# eORCA1 is the only ORCA mesh folded about an F-point pivot (eORCA025/eORCA12 pivot on T points and
# `TripolarGrid` defaults to a U pivot), so the fold parameters of its `Tripolar` mapping identify it.
# A `TripolarGrid` built with `fold_topology = RightFaceFolded` would also match.
const ORCA1Mapping        = Tripolar{<:Any, <:Any, <:Any, RightFaceFolded, FPivot}
const ORCA1UnderlyingGrid = OrthogonalSphericalShellGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:ORCA1Mapping}
const ORCA1GRID           = Union{ORCA1UnderlyingGrid, ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:ORCA1UnderlyingGrid}}

# eORCA025 and eORCA12 both fold about a T-point pivot, so they share a type and only their size
# (1440 and 4320 distinct columns) tells them apart.
const ORCATPivotMapping        = Tripolar{<:Any, <:Any, <:Any, RightCenterFolded, TPivot}
const ORCATPivotUnderlyingGrid = OrthogonalSphericalShellGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:ORCATPivotMapping}
const ORCATPivotGRID           = Union{ORCATPivotUnderlyingGrid, ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:ORCATPivotUnderlyingGrid}}

# On a distributed grid `size` is the local size, so recover the global one
global_Nx(grid) = global_Nx(architecture(grid), grid)
global_Nx(arch, grid) = size(grid)[1]
global_Nx(arch::Distributed, grid) = global_size(arch, size(grid))[1]

is_orca_quarter(grid) = false
is_orca_quarter(grid::ORCATPivotGRID) = global_Nx(grid) == 1440

is_orca_twelfth(grid) = false
is_orca_twelfth(grid::ORCATPivotGRID) = global_Nx(grid) == 4320

# eORCA1's narrowed straits, retuned so their transports approach observations. Each entry is matched to
# the channel cell nearest (longitude, latitude); `width` (m) sets both of its faces (`Inf` removes the
# narrowing) and `depth` (m) its bottom, and anything left out keeps the mesh/bathymetry files' value.
# First-pass linear estimates from a 10-year mean of the spin-up with the files' geometry:
#
#   * Lombok: 0.50 Sv against 2.6 Sv (INSTANT). The files' 13 km by 166 m channel is shallower than the
#     ~240 m sill in ETOPO 2022; at that depth the transport scales to 2.6 Sv at ~49 km.
#   * Ombai: 0.54 Sv against 4.9 Sv (INSTANT), which would need more than a whole cell, so unnarrowed.
#   * Gibraltar: two-layer exchange of 1.24 Sv against 0.80 Sv (Soto-Navarro et al. 2010), scaled from
#     20 km; the exchange is hydraulically controlled, so it scales with width.
#   * Dardanelles: exchange of 0.25 Sv against 0.033 Sv (Ünlüata et al. 1990) through the files' 10 km by 73 m
#     channel. Narrowing it to the 1.3 km a linear scaling asks for NaNs within ~200 iterations (the explicit
#     biharmonic viscosity at the velocity cell between the channel and the full cell beside it), so it keeps
#     the 10 km and is made shallower instead: a two-layer exchange scales as width × depth^(3/2), giving
#     73 m / 7.6^(2/3) ≈ 20 m.
#
# Halmahera (126.5ᵒE, 129.5ᵒE) carries a closed ±1.1 Sv loop round the island, which a width cannot fix,
# and the 119.5ᵒE strait's channel cell is land, so both keep the files' geometry.
const ORCA1_STRAITS = (lombok      = (longitude = 115.5, latitude = -8.5, width = 49e3, depth = 240),
                       ombai       = (longitude = 124.5, latitude = -8.5, width = Inf),
                       gibraltar   = (longitude = -5.4,  latitude = 36.2, width = 13e3),
                       dardanelles = (longitude = 26.9,  latitude = 40.3, depth = 20))

function ORCA1(arch; Nz = 70, Δz_top = 1.5, depth = 5500, Δzmax = nothing,
                     immersed_bottom = GridFittedBottom, minimum_depth = 0,
                     straits = ORCA1_STRAITS)

    z_faces = omip_vertical_discretization(Nz, depth; surface_grid_size = Δz_top,
                                                       maximum_grid_size = Δzmax)

    grid = ORCAGrid(arch;
                    dataset = ORCAOne(),
                    Nz,
                    z = z_faces,
                    halo = (8, 8, 8),
                    with_bathymetry = true,
                    immersed_bottom,
                    major_basins = 1,
                    minimum_depth,
                    straits,
                    active_cells_map = true)

    return grid
end

function ORCA025(arch; Nz = 100, Δz_top = 1.5, depth = 5500, Δzmax = nothing,
                     immersed_bottom = GridFittedBottom, minimum_depth = 20)

    z_faces = omip_vertical_discretization(Nz, depth; surface_grid_size = Δz_top,
                                                       maximum_grid_size = Δzmax)

    grid = ORCAGrid(arch;
                    dataset = ORCAQuarter(),
                    Nz,
                    z = z_faces,
                    halo = (8, 8, 8),
                    with_bathymetry = true,
                    immersed_bottom,
                    major_basins = 1,
                    minimum_depth,
                    active_cells_map = true)

    return grid
end

function ORCA12(arch; Nz = 100, Δz_top = 1.5, depth = 5500, Δzmax = nothing,
                     immersed_bottom = GridFittedBottom, minimum_depth = 20)

    z_faces = omip_vertical_discretization(Nz, depth; surface_grid_size = Δz_top,
                                                       maximum_grid_size = Δzmax)

    grid = ORCAGrid(arch;
                    dataset = ORCATwelfth(),
                    Nz,
                    z = z_faces,
                    halo = (8, 8, 8),
                    with_bathymetry = true,
                    immersed_bottom,
                    major_basins = 1,
                    minimum_depth,
                    active_cells_map = true)

    return grid
end

build_sea_ice_grid(grid, ::Nothing, immersed_bottom) = grid

@kernel function _immerse_latitude_band!(bottom_height, grid, south, north)
    i, j = @index(Global, NTuple)
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    @inbounds z = bottom_height[i, j, 1]
    @inbounds bottom_height[i, j, 1] = ifelse((φ > south) & (φ < north), oftype(z, 100), z)
end

function build_sea_ice_grid(grid, latitudes::Tuple, immersed_bottom)
    south, north = latitudes

    arch       = architecture(grid)
    underlying = grid.underlying_grid

    bottom = Field{Center, Center, Nothing}(underlying)
    parent(bottom) .= parent(bottom_height_field(grid))

    FT = eltype(grid)
    launch!(arch, underlying, :xy, _immerse_latitude_band!, bottom, underlying, convert(FT, south), convert(FT, north))
    fill_halo_regions!(bottom)

    return ImmersedBoundaryGrid(underlying, GridFittedBottom(bottom); active_cells_map = true)
end