using JLD2
using Oceananigans: Field, FieldTimeSeries, Center, interior, set!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.OutputReaders: Cyclical
using Oceananigans.Units: days
using NumericalEarth.DataWrangling: Metadata, Metadatum, native_times
using NumericalEarth.DataWrangling.JRA55: MultiYearJRA55
using NumericalEarth.DataWrangling.SeaWiFS: SeaWiFSMonthly

"""
$(TYPEDSIGNATURES)

SeaWiFS surface chlorophyll (mg m⁻³) on `grid` as a monthly climatology: each month the mean over `years` of that
month's composite, each regridded and inpainted as `Field(Metadatum(...), grid)` does. `SeaWiFSMonthly` holds monthly
composites of single months (1997-09 to 2010-12), not a climatology, so a single year carries that year's anomalies and
gaps; the default `years` are the whole years of the record. The climatology sits on the times of the year 2000, as
`forced_ocean_simulation` used the single year 2000 before, so that restarts (`restart.jl`) keep its phase, and cycles
every `climatology_cycle` (the default is `forced_ocean_simulation`'s for 1958--2018 JRA55).

The 12 monthly means are cached in `dir` as a plain `Nx × Ny × 12` array and read back unless `load_cache = false`;
computing them reads `12 × length(years)` monthly files, which must already be in `dir` (or be downloadable).
The cache is keyed by the grid's size only, so delete it after changing the land mask.
"""
function seawifs_chlorophyll_climatology(grid;
                                         years = 1998:2010,
                                         climatology_cycle = climatology_period(forcing_period(MultiYearJRA55(), DateTime(1958, 1, 1), DateTime(2018, 1, 1))),
                                         dir = restoring_dir[],
                                         load_cache = true)

    save_name = joinpath(dir, "seawifs_chlorophyll_climatology_$(grid_name(grid))_$(first(years))_$(last(years))_v1.jld2")

    if isfile(save_name) & load_cache
        means = jldopen(save_name) do file
            file["means"]
        end
    else
        monthly_means = [Field{Center, Center, Nothing}(grid) for month in 1:12]
        for month in 1:12, year in years
            snapshot = Field(Metadatum(:chlorophyll; dataset = SeaWiFSMonthly(), date = DateTime(year, month, 1), dir), grid)
            interior(monthly_means[month]) .+= interior(snapshot) ./ length(years)
        end
        means = cat((Array(interior(monthly_means[month])) for month in 1:12)...; dims = 3)
        jldsave(save_name; means, years = collect(years))
    end

    monthly_means = [Field{Center, Center, Nothing}(grid) for month in 1:12]
    for month in 1:12
        set!(monthly_means[month], means[:, :, month:month])
    end
    fill_halo_regions!.(monthly_means)

    times = native_times(Metadata(:chlorophyll; dataset = SeaWiFSMonthly(), dates = (DateTime(2000, 1, 1), DateTime(2000, 12, 1)), dir))
    chlorophyll = FieldTimeSeries{Center, Center, Nothing}(grid, times; time_indexing = Cyclical(climatology_cycle))
    set!(chlorophyll, monthly_means)

    return chlorophyll
end
