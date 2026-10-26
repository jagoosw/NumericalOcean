# NumericalOcean: global configurations in Oceananigans 
[![DOI](https://zenodo.org/badge/1407843229.svg)](https://doi.org/10.5281/zenodo.23375796) 

This configuration package assumes familiarity with [Oceananigans](https://github.com/clima/oceananigans.jl/) and to some extent [NumericalEarth](https://github.com/numericalearth/numericalearth.jl/).
## Setup
This repo depends on un-released branches of Oceananigans, ClimaSeaIce, OceanBioME, SeawaterPolynomials, and NumericalEarth. To use install create a julia project:
```console
~$ mkdir AnOceanProject
~$ cd AnOceanProject
~$ touch Project.toml
~$ julia --project -t auto
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.12.7 (2026-08-15)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org release
|__/                   |

julia> using Pkg

julia> Pkg.add([PackageSpec(url="https://github.com/NumericalEarth/NumericalEarth.jl.git", rev="jsw/omip-bgc"),
                PackageSpec(url="https://github.com/NumericalEarth/NumericalEarth.jl.git", rev="jsw/omip-bgc", subdir="experiments/OMIPSimulations"),
                PackageSpec(url="https://github.com/Clima/Oceananigans.jl.git", rev="jsw/for-oaemip"),
                PackageSpec(url="https://github.com/jagoosw/ClimaSeaIce.jl.git", rev="jsw/for-oaemip"),
                PackageSpec(url="https://github.com/OceanBioME/OceanBioME.jl.git", rev="jsw/for-oaemip"),
                PackageSpec(url="https://github.com/clima/SeawaterPolynomials.jl.git", rev="ss/conversion-functions"),
                PackageSpec(url="https://github.com/jagoosw/GlobalOcean", rev="v0.1.0")])
```
The branches have to be added first, in the same call or before `NumericalOcean`, because Pkg ignores the `[sources]` in NumericalOcean's `Project.toml` when it is installed as a package, and `OMIPSimulations` and some dependencies are not registered.
This puts Oceananigans, OceanBioME and NumericalEarth in your project, so they can be used in scripts; you will also need to add CUDA to run on the GPU:
```console
julia> Pkg.add("CUDA")
```
Alternatively, clone this repo and use it as the project (`julia --project=. -e 'using Pkg; Pkg.instantiate()'`), where its `[sources]` are used directly.
## Physics only setup
As shown in `examples/1degree.jl`, the default configuration is an ocean-sea ice model model forced with JRA55 atmosphere and land, configured to setup (semi)-validated closures.
The model is initialised with WOA climatology temperature and salinity, and ECCO sea ice concentration and thickness.
```julia
using NumericalOcean
using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Dates
using CUDA

auto_config_directories!()

grid = ORCA1(GPU())

simulation = forced_ocean_simulation(grid)

pickup_restart!(simulation, DateTime(1958, 1, 1))

omip_diagnostics!(simulation)
checkpointer!(simulation)

run!(simulation)
```
The `grid` can be any global grid, but automatic setup is only checked for `ORCA1`, `ORCA025`, and `ORCA12` which are nominally $1\degree$, $1/4\degree$, and $1/12\degree$.
The whole model is setup through `forced_ocean_simulation` so everything can be modified by passing arguments. 
For example, to release add a passive tracer `c` and release it in a circle:
```julia
c_release = Forcing((λ, φ, z, t) -> (sqrt(λ^2 + φ^2) < 1)&(t < 86400))
simulation = forced_ocean_simulation(grid; 
                                     tracers = (:c, ),
                                     forcing = (; c = c_release))
```
Please see the `forced_ocean_simulation` docstring for details, documentation to follow.

## Biogeochemistry
`forced_ocean_simulation` tried to automatically setup the necessary additional features for biogeochemistry: rivers, dust, and gas exchange.
If you `using OceanBioME` a `NumericalEarth` extension is loaded that automatically sets up $CO_2$ and $O_2$ gas exchange.
River loads are added from RivR2O co-located with the land model river mouths.
Dust is supplied by [UnifiedBGC](https://roms-tools.readthedocs.io/en/latest/datasets_overview.html#unified-bgc-dataset) from ROMS-tools.
BGC tracers are also initialised from UnifiedBGC.
For example setting up an equivalent to the [MITgcm DIC package](https://mitgcm.readthedocs.io/en/latest/phys_pkgs/dic.html) which tracks $PO_4$, $Fe$, $DOP$, $DIC$, $Alk$ and $O_2$ with Mona Loa atmospheric $pCO_2$:
```julia
using OceanBioME
surface_PAR = PARFromShortwave(grid)
light_attenuation = PrescribedAttenuationPAR(grid, surface_PAR;
                                             attenuation = PrescribedChlorophyllAttenuation(first_color_fraction = 0.0),
                                             attenuation_discrete_form = true)

biogeochemistry = MITgcmDIC(grid;
                            light_attenuation,
                            open_bottom = false,
                            implicit_sinking = true,
                            store_flux = true)

pCO₂ = MaunaLoaCO₂(arch; dates = DateTime(1990, 1, 15):Month(1):DateTime(1990, 12, 15),
                         start_date = DateTime(1990, 1, 1),
                         period = 365days)

simulation = forced_ocean_simulation(grid;
                                     biogeochemistry,
                                     atmosphere_tracers = (; pCO₂),
                                     jra55_dataset = RepeatYearJRA55(),
                                     start_date = DateTime(1990, 1, 1))
```
`MaunaLoaCO₂(arch)` with no arguments gives the monthly record for 1958--2017, cycled with the multi-year JRA55 forcing (see below); here we instead repeat 1990 to go with repeat-year forcing.
We have to use `PARFromShortwave` to get the downwelling radiation from the atmosphere component of the model and, since this model doesn't prognose phytoplankton, we have to use the chlorophyll climatology from the radiation model to attenuation the light through `PrescribedAttenuationPAR` with `PrescribedChlorophyllAttenuation`.

### Alkalinity release
You can prescribe a Guassian alkalinity release following the OAEMIP protocol by adding a forcing:
```julia
alkalinity_release = GaussianInjection(grid, (-49.5, 0), 0.0; duration = 20days) # grid, center, start time; duration
simulation = forced_ocean_simulation(grid;
                                     biogeochemistry,
                                     atmosphere_tracers = (; pCO₂),
                                     jra55_dataset = RepeatYearJRA55(),
                                     start_date = DateTime(1990, 1, 1),
                                     forcing = (; Alk2 = release))
```
Default is 1Tmol over 1 year in a 250km standard deviation patch cutoff at 3 standard deviations over the surface cells only.

## Forcing data
To change the automatically configured data and output directories you can change the value of `NumericalOcean.forcing_dir`, `NumericalOcean.restoring_dir`, `NumericalOcean.staging_dir`, `NumericalOcean.bgc_dir`, and `NumericalOcean.output_dir`.
They are [`Ref`s ](https://discourse.julialang.org/t/what-is-ref/47610) so you set it like `NumericalOcean.bgc_dir[] = "~/bgc_data/"`.
`staging_dir` is only needed when you have a slow `forcing_dir` so that data can be asynchronously loaded to a faster drive before it is needed.

### Normally...
When you first load this package it will create directories `data/forcing_data`, `data/climatology`, `data/output`, and `data/staging`.
When you run the model for the first time it is going to try and download the whole JRA55 surface dataset which is around 700GiB, it is probably best to have arranged for this before running `scripts/download_forcing.jl`. 
Repeat year (1990) is much smaller, to use pass `jra55_dataset = RepeatYearJRA55()`.
It is also sufficiently small to be fully loaded into GPU memory on a H100 if you pass `backend_size = 2920` which confers meaningful speed up.
To use the biogeochemistry setup you will also need to download [UnifiedBGC](https://roms-tools.readthedocs.io/en/latest/datasets_overview.html#unified-bgc-dataset) and the [RivR2O](https://zenodo.org/records/14889524) dataset.

### On Bouchet
If you are running on Bouchet the package will point to a predownloaded version of all the data, if you can not access it please email [me](mailto:jago.strong-wright@yale.edu).
On Bouchet `staging_dir` will automatically configure to use the fast NVMe storage on the nodes and will be used if you pass `staging=true` to `forced_ocean_simulation` (although I have found this results in minimal performance gain).

## Restart files
`pickup_restart!(simulation, date)` starts a simulation from the checkpoint of the OMIP spin-up (`orca_checkpoint_iteration*.jld2` in `NumericalOcean.OMIP_RESTART_DIR[]`, which points at the spin-up on Bouchet) whose forcing date is closest to `date`, and sets the model clock so the forcing is picked up where the checkpointed run had it.
The spin-up was forced with 1958-2018 JRA55 and checkpointed every 360 days, so checkpoints fall on different dates in each cycle; by default only the last cycle is searched (pass `cycle = n`, 0-based, for another).
```julia
simulation = forced_ocean_simulation(grid)

# the last-cycle state nearest 1 January 1958, with the clock at its time since 1958-01-01
pickup_restart!(simulation, DateTime(1958, 1, 1))
```
For a run with a different `start_date` pass the same `start_date` to `pickup_restart!`: only checkpoints dated on or after it can be used.
For repeat-year runs pass `repeat_year = true` to pick the checkpoint nearest `date` (in any year) and set the clock to its time of year since `start_date`:
```julia
simulation = forced_ocean_simulation(grid; jra55_dataset = RepeatYearJRA55(), start_date = DateTime(1990, 1, 1))

pickup_restart!(simulation, DateTime(1990, 1, 1); start_date = DateTime(1990, 1, 1), repeat_year = true)
```
Only the model state is restored: the iteration is reset to zero and output and callback schedules start afresh from the new time.
State not in the checkpoint (e.g. biogeochemical tracers when restarting from the physics-only spin-up) keeps the values it was initialised with.

## Performance
Currently I get around 90 simulated years per day with the MITgcm BGC, and 140 simulated years per day with physics only on a H100 on Bouchet.

## Miscellaneous
I have been requesting 8 CPUs with 16G of memory each (probably excessive) and running julia with `--threads=4 --gcthreads=4,1` which helped performance a bit.

# Spun-up configurations
Spinup checkpoints are available on [Bouchet](https://docs.ycrc.yale.edu/clusters/bouchet/), if you would like access elsewhere we can share them (each one is ~4.5GiB).
## NumericalOcean1-JRA55_v0.1.0
As built by `forced_ocean_simulation(ORCA1(GPU()))` with default arguments.

| | |
|---|---|
| **Grid** | eORCA1 tripolar ($1\degree$), 70 z-star levels to 5500 m (1.5 m at the surface) |
| **Time step** | 90 minutes, split Runge--Kutta 3 ([Silvestri et al., 2026](https://doi.org/10.22541/essoar.15002225/v1)), split-explicit free surface |
| **Equation of state** | TEOS-10 (Conservative Temperature and Absolute Salinity) |
| **Momentum advection** | WENO vector invariant, 5th order ([Silvestri et al., 2024](https://doi.org/10.1029/2023MS004130)) |
| **Tracer advection** | WENO, 7th order |
| **Vertical mixing** | CATKE ([Wagner et al., 2025](https://doi.org/10.1029/2024MS004522)) |
| **Background diffusivity** | κ = max(2×10⁻⁶, 10⁻⁵ \|sin φ\|) m² s⁻¹ ([Henyey et al., 1986](https://doi.org/10.1029/JC091iC07p08487)), ν = 3×10⁻⁵ m² s⁻¹ |
| **Mesoscale eddies** | GM ([Gent & McWilliams, 1990](https://doi.org/10.1175/1520-0485(1990)020%3C0150:IMIOCM%3E2.0.CO;2)) and Redi ([Redi, 1982](https://doi.org/10.1175/1520-0485(1982)012%3C1154:OIMBCR%3E2.0.CO;2)) on triads ([Griffies et al., 1998](https://doi.org/10.1175/1520-0485(1998)028%3C0805:IDIAZC%3E2.0.CO;2)), κ_GM = 800 m² s⁻¹, κ_Redi = 1200 m² s⁻¹ |
| **Horizontal viscosity** | biharmonic, ν₄ = (cell area)² / 50 days |
| **Sea ice** | ClimaSeaIce, 4 thickness categories, EVP dynamics ([Kimmritz et al., 2017](https://doi.org/10.1016/j.ocemod.2017.05.006)), landfast basal stress ([Lemieux et al., 2015](https://doi.org/10.1002/2014JC010678)), incremental remapping ([Lipscomb & Hunke, 2004](https://doi.org/10.1175/1520-0493(2004)132%3C1341:MSITUI%3E2.0.CO;2)) |
| **Atmosphere and runoff** | JRA55-do 1958--2018 ([Tsujino et al., 2018](https://doi.org/10.1016/j.ocemod.2018.07.002)) |
| **Air-sea fluxes** | COARE similarity theory ([Edson et al., 2013](https://doi.org/10.1175/JPO-D-12-0173.1)) |
| **Ice-ocean heat flux** | three-equation ([Holland & Jenkins, 1999](https://doi.org/10.1175/1520-0485(1999)029%3C1787:MTIOIA%3E2.0.CO;2)) |
| **Shortwave penetration** | two-band, SeaWiFS chlorophyll monthly climatology (1998–2010) |
| **Salinity restoring** | to monthly WOA, piston velocity 0.254 m day⁻¹ |
| **Initial conditions** | WOA annual temperature and salinity, ECCO4 sea ice (January 1993), ocean at rest |
| **Spin-up** | 5 cycles of 1958--2018 JRA55 |
