module NumericalOcean

export auto_config_directories!,
       ORCA1, ORCA025, ORCA12,
       forced_ocean_simulation,
       MaunaLoaCO₂,
       PrescribedChlorophyllAttenuation,
       seawifs_chlorophyll_climatology,
       GaussianInjection, GaussianInjectionLowMemory,
       omip_diagnostics!, oaemip_diagnostics!, checkpointer!,
       add_output_writer!, flush_output!,
       pickup_restart!

using Dates
using DocStringExtensions
using NumericalEarth
using Oceananigans

const on_bouchet = isdir("/nfs/roberts/pi/pi_ey239")
const forcing_dir = Ref(on_bouchet ? "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/forcing_data" : "data/forcing_data")
const restoring_dir = Ref(on_bouchet ? "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/climatology" : "data/climatology")
const staging_dir = Ref(on_bouchet ? joinpath("/tmp", "jra55_" * get(ENV, "SLURM_JOB_ID", string(getpid()))) : "data/staging")
const bgc_dir = Ref(on_bouchet ? "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data" : "data")
const output_dir = Ref(on_bouchet ? (joinpath("/nfs/roberts/pi/pi_ey239/", get(ENV, "USER", "no_user_this_should_error_as_long_as_no_one_makes_this_folder"), "output")) :
                                    "data/output")

include("forcing_cycles.jl")
include("data_management.jl")
include("grids.jl")
include("defaults.jl")
include("unified_bgc.jl")
include("river_input.jl")
include("chlorophyll_climatology.jl")
include("interpolated_iron_dust.jl")
include("biogeochemistry.jl")
include("light_attenuation.jl")
include("mauna_loa.jl")
include("alkalinity_injection.jl")
include("construct_model.jl")
include("diagnostics.jl")
include("restart.jl")

function __init__()
    auto_config_directories!()
    @info "Setup directories. Forcing: $(forcing_dir[]), restoring: $(restoring_dir[]), output: $(output_dir[])."
end

end # module NumericalOcean
