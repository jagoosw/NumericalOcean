using JLD2: JLD2, jldopen
using Oceananigans: restore_prognostic_state!
using Oceananigans.OutputWriters: load_checkpoint_state
using Oceananigans.TimeSteppers: Clock, update_state!
using Oceananigans.Units: days
using NumericalEarth.EarthSystemModels: EarthSystemModel, components
using NumericalEarth.DataWrangling: Metadata, native_times, all_dates, compute_native_date_range

const OMIP_RESTART_DIR = Ref(on_bouchet ? "/nfs/roberts/pi/pi_ey239/js5256/NumericalOceanSpinups/NumericalOcean1-JRA55_v0.2.0" : "data/spinups")

#####
##### Finding the checkpoint for a date
#####

"""
$(TYPEDSIGNATURES)

Return the `(path, time)` of every `prefix_iteration*.jld2` checkpoint in `dir`, sorted by `time`
(model seconds, read from the checkpointed clock).
"""
function checkpoint_times(dir, prefix)
    regex = Regex("^" * prefix * "_iteration[0-9]+\\.jld2\$")
    paths = [joinpath(dir, file) for file in readdir(dir) if occursin(regex, file)]
    isempty(paths) && throw(ArgumentError("no $(prefix)_iteration*.jld2 checkpoints in $dir"))
    times = [jldopen(path, "r") do file; checkpointed_time(file["simulation/model/clock"]); end for path in paths]
    order = sortperm(times)
    return paths[order], times[order]
end

# The `Checkpointer` writes the whole `Clock` as one dataset, but a group holding `time` is also accepted
checkpointed_time(clock::JLD2.Group) = Float64(clock["time"])
checkpointed_time(clock) = Float64(clock.time)

"""
$(TYPEDSIGNATURES)

The period (seconds) over which the JRA55 `dataset` forcing between `start_date` and `end_date` cycles, as its
`FieldTimeSeries` infers it: the last node is included, so 1958--2018 cycles every 21915 days + 3 hours, not 21915 days.
This is for the 3-hourly atmospheric state; the daily river fluxes cycle every 21916 days.

This is the cycle of runs built before `forced_ocean_simulation` set the period explicitly; later runs cycle
every [`forcing_period`](@ref) (21915 days for 1958--2018).
"""
function forcing_cycle_length(dataset, start_date, end_date; name = :temperature)
    dates = compute_native_date_range(all_dates(dataset, name), start_date, end_date)
    times = native_times(Metadata(name; dataset, dates))
    return times[end] - times[1] + (times[end] - times[end-1])
end

"""
$(TYPEDSIGNATURES)

The model time (seconds since `start_date`) to give a checkpoint of the state on `checkpoint_date`.
With `repeat_year`, `checkpoint_date` is first moved to `start_date`'s year (by calendar date, so the 29th of
February becomes the 28th), and the time wrapped into `[0, 365days)`.
"""
function restart_clock_time(checkpoint_date, start_date; repeat_year = false)
    repeat_year || return period_to_seconds(checkpoint_date - start_date)

    same_year_date = checkpoint_date + Year(year(start_date) - year(checkpoint_date))
    return mod(period_to_seconds(same_year_date - start_date), 365days)
end

"""
$(TYPEDSIGNATURES)

Choose the checkpoint (from `paths` at model `times`) whose date in the forcing cycle is closest to `date`.
The checkpointed run started on `restart_start_date` and cycled the forcing every `cycle_length` seconds;
only checkpoints from `cycle` (0-based, default: the last one with any checkpoints) whose clock time for
a run starting on `start_date` would be non-negative are considered.

Return `(path, checkpoint_date, cycle)`.
"""
function select_restart(paths, times, date;
                        start_date,
                        repeat_year = false,
                        restart_start_date = DateTime(1958, 1, 1),
                        cycle_length,
                        cycle = nothing)

    cycles = [fld(t, cycle_length) for t in times]
    dates  = [restart_start_date + seconds_to_period(t - c * cycle_length) for (t, c) in zip(times, cycles)]
    cycle  = something(cycle, maximum(cycles))

    candidates = [n for n in eachindex(paths)
                  if cycles[n] == cycle && restart_clock_time(dates[n], start_date; repeat_year) >= 0]

    isempty(candidates) && throw(ArgumentError("no checkpoint in cycle $cycle from on or after $start_date"))

    n = candidates[argmin([abs(dates[n] - date) for n in candidates])]

    return paths[n], dates[n], cycle
end

#####
##### Restoring it
#####

set_restart_clock!(::Nothing, time) = nothing
set_restart_clock!(component::Simulation, time) = set_restart_clock!(component.model, time)
set_restart_clock!(component, time) = hasproperty(component, :clock) ? set_restart_clock!(component.clock, time) : nothing

function set_restart_clock!(clock::Clock, time)
    clock.time = time
    clock.iteration = 0
    return nothing
end

function set_restart_clock!(model::EarthSystemModel, time)
    set_restart_clock!(model.clock, time)

    for component in components(model)
        set_restart_clock!(component, time)
    end

    # the restore loaded the prescribed forcing's in-memory windows at the checkpoint's own clock, which is
    # elsewhere in the forcing record unless the new run cycles it identically; the coupled `update_state!`
    # doesn't refresh them, so it would read outside the windows
    update_state!(model.atmosphere)
    update_state!(model.land)
    update_state!(model.radiation)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Set `simulation`'s model to the checkpoint of the OMIP spin-up closest to `date`, and its clock to that
checkpoint's date as seconds since `start_date` (the `start_date` the simulation was built with).

With `repeat_year = true` (for `RepeatYearJRA55` runs, with `start_date = DateTime(1990, 1, 1)`), the clock is
instead set to the checkpoint's time of year since `start_date`, so e.g. `pickup_restart!(simulation, DateTime(1990, 1, 1);
start_date = DateTime(1990, 1, 1), repeat_year = true)` starts from the spin-up's state closest to 1 Jan 1990 at the
forcing's matching time of year.

The checkpoints are `dir/prefix_iteration*.jld2`, written by a run forced by `MultiYearJRA55` from `restart_start_date`
to `restart_end_date`, which it cycled every `cycle_length` seconds. The default, 21915 days + 3 hours for 1958--2018,
is the period that run's forcing inferred for itself (see [`forcing_cycle_length`](@ref)); pass
`cycle_length = forcing_period(MultiYearJRA55(), restart_start_date, restart_end_date)` (21915 days) for checkpoints
from runs whose forcing is phase-locked by `forced_ocean_simulation`. A checkpoint's date is that of the atmospheric
forcing it was last forced with. Only checkpoints from `cycle` (0-based, default: the last cycle present) are considered.

The new run's forcing is all phase-locked to the calendar, but in a checkpointed run with inferred periods only the
3-hourly JRA55 atmosphere and radiation sit at the dates assigned here: the rest drifted against them, so a state from its `n`th
cycle was spun up under rivers `n × 0.875` days behind where the new run puts them on the same date, WOA salinity
restoring `n × 15.125` days ahead and SeaWiFS chlorophyll `n × 44.875` days behind, the climatologies a further day
ahead for each 29 February (WOA) or behind for each 365-day year (SeaWiFS) since `restart_start_date`.

Only the model state is restored, and the iteration is set to zero, so the simulation initialises as a fresh run
(output writer and callback schedules start from the new time, and the first step is forward Euler). Tracers and
other state missing from the checkpoint (e.g. biogeochemical tracers when restarting from a physics-only spin-up)
keep the values they were given when the simulation was built. `simulation.stop_time` is not changed: it is model
time, so with `repeat_year = false` it must exceed the restored clock time.

Return the path of the checkpoint used.
"""
function pickup_restart!(simulation, date;
                         start_date = DateTime(1958, 1, 1),
                         repeat_year = false,
                         dir = OMIP_RESTART_DIR[],
                         prefix = "orca_checkpoint",
                         restart_start_date = DateTime(1958, 1, 1),
                         restart_end_date = DateTime(2018, 1, 1),
                         cycle_length = forcing_cycle_length(MultiYearJRA55(), restart_start_date, restart_end_date),
                         cycle = nothing)

    paths, times = checkpoint_times(dir, prefix)
    path, checkpoint_date, cycle = select_restart(paths, times, date;
                                                  start_date, repeat_year, restart_start_date, cycle_length, cycle)

    time = restart_clock_time(checkpoint_date, start_date; repeat_year)

    @info "Picking up from $path (cycle $cycle, $checkpoint_date), model time $(prettytime(time)) since $start_date"

    state = load_checkpoint_state(path)
    restore_prognostic_state!(simulation.model, state.model)
    set_restart_clock!(simulation.model, time)
    simulation.initialized = false

    return path
end
