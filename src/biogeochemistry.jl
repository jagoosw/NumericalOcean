using OceanBioME: DiscreteBiogeochemistry, NutrientsPlanktonDetritus
using NumericalEarth.Oceans: RiverConcentration

import NumericalEarth.Oceans: default_freshwater_tracer_content

# Riverine iron is only carried by river runoff (not rain, snow or melt), at 0.01 mmol Fe / m³
default_freshwater_tracer_content(::Val{:Fe}, ::DiscreteBiogeochemistry) =
    RiverConcentration(convert(Oceananigans.defaults.FloatType, 0.01))

# Compute every NutrientsPlanktonDetritus transition in its own kernel rather than inline in the tracer tendency
Oceananigans.Biogeochemistry.separate_tracer_transitions(bgc::DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}) =
    Oceananigans.Biogeochemistry.required_biogeochemical_tracers(bgc)
