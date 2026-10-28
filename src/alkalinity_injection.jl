# Alkalinity releases (from OAEMIP's src/alkalinity_injection.jl). Both forms give identical results: the first
# stores the release pattern in a field (more memory, less computation), the second recomputes it every evaluation.
using GPUArraysCore: @allowscalar
using KernelAbstractions: @kernel, @index
using Oceananigans.AbstractOperations: Integral
using Oceananigans.Grids: λnode, φnode
using Oceananigans.ImmersedBoundaries: immersed_cell
using Oceananigans.OrthogonalSphericalShellGrids: spherical_distance

import Adapt

"""
$(TYPEDSIGNATURES)

A release into the surface cells with a Gaussian footprint, as a `forcing` (in tracer units per second) on the tracer
it is passed for, e.g. `forcing = (; Alk2 = GaussianInjection(grid, (-49.5, 0), 5 * 365days))`.

`total_injection` (tracer units × m³, e.g. 1000e12 mmol = 1 Tmol of alkalinity) is spread over the surface cells
within `cutoff` (m) of `center = (λ, φ)` (degrees), weighted by `exp(-r² / 2σ²)` with `σ` in m, at a constant rate
from `start_time` to `start_time + duration` (s of model time).
"""
struct GaussianInjection{M, R, S, D}
          mask :: M
          rate :: R
    start_time :: S
      duration :: D
end

# the mask is only nonzero in the surface layer, so this is the injection into the surface cell only
@inline (gi::GaussianInjection)(i, j, k, grid, clock, fields) =
    @inbounds gi.mask[i, j, k] * gi.rate * (gi.start_time <= clock.time < gi.start_time + gi.duration)

Adapt.adapt_structure(to, gi::GaussianInjection) =
    GaussianInjection(Adapt.adapt(to, gi.mask),
                      Adapt.adapt(to, gi.rate),
                      Adapt.adapt(to, gi.start_time),
                      Adapt.adapt(to, gi.duration))

function GaussianInjection(grid, center, start_time;
                           σ = 250000,
                           cutoff = 3σ,
                           total_injection = 1000e12, # mmol, i.e. 1 Tmol
                           duration = 365days)

    mask = CenterField(grid)
    set_injection_site!(mask, grid, center, σ, cutoff)

    return GaussianInjection(mask, injection_rate(mask, total_injection, duration), start_time, duration)
end

"""
$(TYPEDSIGNATURES)

[`GaussianInjection`](@ref) with the footprint recomputed every evaluation rather than stored in a field.
"""
struct GaussianInjectionLowMemory{C, W, O, P, R, S, D}
        center :: C
         width :: W
        cutoff :: O
        radius :: P
          rate :: R
    start_time :: S
      duration :: D
end

@inline function (gi::GaussianInjectionLowMemory)(i, j, k, grid, clock, fields)
    mask = injection_mask(i, j, grid, gi.center, gi.width, gi.cutoff, gi.radius)
    return mask * gi.rate * (gi.start_time <= clock.time < gi.start_time + gi.duration) * (k == grid.Nz)
end

function GaussianInjectionLowMemory(grid, center, start_time;
                                    σ = 250000,
                                    cutoff = 3σ,
                                    total_injection = 1000e12, # mmol, i.e. 1 Tmol
                                    duration = 365days)

    mask = CenterField(grid)
    set_injection_site!(mask, grid, center, σ, cutoff)

    radius = planet_radius(grid)

    return GaussianInjectionLowMemory(center, σ, cutoff, radius, injection_rate(mask, total_injection, duration), start_time, duration)
end

# the rate (per m³ of footprint per second) at which `total_injection` is released over `duration`
function injection_rate(mask, total_injection, duration)
    volume = @allowscalar first(Field(Integral(mask)))
    return total_injection / volume / duration
end

# `Oceananigans.defaults` is a mutable global, which a GPU kernel can't read, so the radius is passed in
planet_radius(grid) = convert(eltype(grid), Oceananigans.defaults.planet_radius)

@inline function injection_mask(i, j, grid, center, σ, cutoff, radius)
    λ₀, φ₀ = center
    λ = λnode(i, j, 1, grid, Center(), Center(), Center())
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())

    r = spherical_distance(λ₀, φ₀, λ, φ, radius)

    raw_mask = exp(-r^2 / (2σ^2))
    return ifelse(r > cutoff, zero(raw_mask), raw_mask) * !immersed_cell(i, j, grid.Nz, grid)
end

set_injection_site!(mask, grid, center, σ, cutoff) =
    launch!(architecture(grid), grid, :xy, _set_injection_site!, mask, grid, center, σ, cutoff, planet_radius(grid))

@kernel function _set_injection_site!(mask, grid, center, σ, cutoff, radius)
    i, j = @index(Global, NTuple)
    @inbounds mask[i, j, grid.Nz] = injection_mask(i, j, grid, center, σ, cutoff, radius)
end
