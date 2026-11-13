# River iron content and the separate-transition override for NutrientsPlanktonDetritus

using NumericalEarth.Oceans: RiverConcentration, default_freshwater_tracer_content
using Oceananigans.Biogeochemistry: separate_tracer_transitions, required_biogeochemical_tracers
using Oceananigans.Fields: ZeroField

@testset "default_freshwater_tracer_content" begin
    grid = RectilinearGrid(size = (2, 2, 2), extent = (1, 1, 1))
    bgc = MITgcmDIC(grid)

    content = default_freshwater_tracer_content(Val(:Fe), bgc)
    @test content isa RiverConcentration
    @test content.concentration == convert(Oceananigans.defaults.FloatType, 0.01)
    @test content.concentration isa Oceananigans.defaults.FloatType

    # other tracers keep NumericalEarth's default (no content), and so does Fe without DiscreteBiogeochemistry
    @test default_freshwater_tracer_content(Val(:DIC), bgc) isa ZeroField
    @test default_freshwater_tracer_content(Val(:PO₄), bgc) isa ZeroField
    @test default_freshwater_tracer_content(Val(:Fe), nothing) isa ZeroField
end

@testset "separate_tracer_transitions" begin
    grid = RectilinearGrid(size = (2, 2, 2), extent = (1, 1, 1))
    for bgc in (MITgcmDIC(grid), LOBSTER(grid), NPZD(grid), ImplicitBiology(grid))
        @test bgc isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}
        @test separate_tracer_transitions(bgc) == required_biogeochemical_tracers(bgc)
    end
    @test separate_tracer_transitions(nothing) == ()
end
