# # Spinning up the ``1\degree`` configuration
#
# This example builds the `NumericalOcean1-JRA55_v0.2.0` spin-up (see [Spun-up configurations](@ref)): the
# default ocean--sea ice configuration on eORCA1, run from rest through 5 cycles of 1958--2018 JRA55 forcing
# (5 × 21915 days, about 300 years). At around 140 simulated years per day on an H100 this takes a little over
# two days, so it is not run when the documentation is built.
#
# ## Packages

using NumericalOcean
using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Dates
using CUDA

auto_config_directories!()

# ## Grid
#
# The nominal ``1\degree`` eORCA1 tripolar grid on the GPU.

grid = ORCA1(GPU())

# ## Simulation
#
# Everything is left at its default: multi-year JRA55 forcing from 1 January 1958, temperature and salinity
# initialised from the WOA annual climatology (converted to TEOS-10), sea ice from ECCO4 (January 1993), the
# ocean at rest, and `stop_time = 5 * 21915days`. There is no [`pickup_restart!`](@ref): this run is what
# the restarts come from.

simulation = forced_ocean_simulation(grid)

# ## Output and checkpoints
#
# The OMIP diagnostics, and a checkpoint every 360 days. The `"orca"` prefix gives the
# `orca_checkpoint_iteration*.jld2` files that [`pickup_restart!`](@ref) looks for: to restart from your own
# spin-up set `NumericalOcean.OMIP_RESTART_DIR[]` to the output directory.

omip_diagnostics!(simulation; async = true)
checkpointer!(simulation; filename_prefix = "orca")

# ## Running
#
# A spin-up this long will usually span several jobs: after the first, resume from the latest checkpoint
# with `run!(simulation; pickup = true)`.

try
    run!(simulation)
finally
    flush_output!(simulation) # write anything still queued
end

#md # ## Results
#md #
#md # ![Surface fields in the final 10 years](../assets/figures/spinup/physics_surface.mp4)
