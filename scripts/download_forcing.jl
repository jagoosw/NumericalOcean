using NumericalOcean
using Dates
using Downloads: Downloads
using Oceananigans: CPU

using NumericalEarth.DataWrangling: Metadata, Metadatum, WOAAnnual, WOAMonthly, ECCO4Monthly
using NumericalEarth.DataWrangling.JRA55: RepeatYearJRA55, MultiYearJRA55
using NumericalEarth.DataWrangling.SeaWiFS: SeaWiFSMonthly
using NumericalEarth.DataWrangling.ORCA: ORCAOne
using WorldOceanAtlasTools  

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

forcing_dir       = NumericalOcean.forcing_dir[]
climatology       = NumericalOcean.restoring_dir[]
sea_ice_date      = DateTime(1993, 1, 1)
chlorophyll_dates = (DateTime(1998, 1, 1), DateTime(2010, 12, 1)) # the whole years of SeaWiFS, for its climatology

jra55 = MultiYearJRA55()
jra55_variables = (:eastward_velocity, :northward_velocity, :temperature, :specific_humidity,
                   :sea_level_pressure, :rain_freshwater_flux, :snow_freshwater_flux,
                   :downwelling_shortwave_radiation, :downwelling_longwave_radiation,
                   :river_freshwater_flux, :iceberg_freshwater_flux)

for name in jra55_variables
    @info "JRA55 $name"
    Downloads.download(Metadata(name; dataset = jra55, dir = forcing_dir))
end

for name in (:temperature, :salinity)
    @info "WOA Annual $name"
    Downloads.download(Metadata(name; dataset = WOAAnnual(), dir = climatology))
    @info "WOA Monthly $name"
    Downloads.download(Metadata(name; dataset = WOAMonthly(), dir = climatology))
end

for name in (:sea_ice_thickness, :sea_ice_concentration)
    @info "ECCO4Monthly $name"
    Downloads.download(Metadata(name; dataset = ECCO4Monthly(), dates = (sea_ice_date, sea_ice_date), dir = climatology))
end

@info "SeaWiFS chlorophyll"
Downloads.download(Metadata(:chlorophyll; dataset = SeaWiFSMonthly(), dates = chlorophyll_dates, dir = climatology))

@info "ORCA eORCA1 mesh + bathymetry"
Downloads.download(Metadatum(:mesh_mask;     dataset = ORCAOne()))
Downloads.download(Metadatum(:bottom_height; dataset = ORCAOne()))

@info "SeaWiFS chlorophyll climatology on ORCA1 (cached in the climatology directory)"
seawifs_chlorophyll_climatology(ORCA1(CPU()))

@info "done"
