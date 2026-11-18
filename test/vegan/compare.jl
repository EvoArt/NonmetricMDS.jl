# Compare NonmetricMDS.jl with vegan::metaMDS (monoMDS) on the dune data.
# Usage: julia --project=. test/vegan/compare.jl OUTDIR
using NonmetricMDS, DelimitedFiles, Test, Random
dir = ARGS[1]
rd(f) = readdlm(joinpath(dir, f), ','; header=true)[1]
D = Float64.(rd("D.csv"))
Xv = Float64.(rd("points.csv"))
sv = rd("stress.csv")[1, 1]
r = nmds(D; k=2, tries=50, rng=Xoshiro(1))
println("NonmetricMDS stress ", r.stress, "  vegan ", sv)
@testset "vegan metaMDS" begin
    @test r.stress <= sv + 1e-3                    # at least as good a solution
    @test procrustes(Xv, r.points).rmse < 0.05 * sqrt(sum(abs2, Xv) / size(Xv, 1))
end
