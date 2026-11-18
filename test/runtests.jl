using NonmetricMDS, Test, LinearAlgebra, Random, Statistics, Distances, StableRNGs

eucl(X) = pairwise(Euclidean(), X; dims=1)

@testset "NonmetricMDS" begin
    @testset "isotonic regression" begin
        rng = StableRNG(1)
        for _ in 1:50
            y = randn(rng, 30)
            f = NonmetricMDS.isotonic(y)
            @test issorted(f)
            @test sum(f) ≈ sum(y)                       # PAVA preserves the mean
            # optimal: no small monotone perturbation lowers the squared error
            for _ in 1:20
                g = sort(f .+ 0.01 .* randn(rng, 30))
                @test sum(abs2, f .- y) <= sum(abs2, g .- y) + 1e-12
            end
        end
        @test NonmetricMDS.isotonic([3.0, 1.0, 2.0]) == [2.0, 2.0, 2.0]
        @test NonmetricMDS.tie_blocks([1, 1, 2, 3, 3, 3]) == [1:2, 3:3, 4:6]
    end

    @testset "PCoA recovers Euclidean configurations" begin
        rng = StableRNG(2)
        X = randn(rng, 25, 3)
        p = pcoa(eucl(X); k=3)
        @test eucl(p.points) ≈ eucl(X)
        @test sum(p.explained) ≈ 1
        @test issorted(p.eigenvalues; rev=true)
        @test size(pcoa(eucl(X); k=2).points) == (25, 2)
    end

    @testset "NMDS recovers configurations" begin
        rng = StableRNG(3)
        X = randn(rng, 30, 2)
        D = eucl(X)
        r = nmds(D; k=2, tries=5, rng=StableRNG(4))
        @test r.stress < 1e-3
        @test procrustes(X, r.points).rmse < 0.02
        # non-metric: any monotone transform gives the same solution
        r3 = nmds(D .^ 3; k=2, tries=5, rng=StableRNG(4))
        @test r3.stress < 1e-3
        @test procrustes(X, r3.points).rmse < 0.02
        @test r.repeats >= 1
        @test size(r.points) == (30, 2)
        # centred and on principal axes
        @test maximum(abs, mean(r.points; dims=1)) < 1e-10
        C = cov(r.points)
        @test abs(C[1, 2]) < 1e-8 && C[1, 1] >= C[2, 2]
        s = shepard(r)
        @test issorted(s.dissimilarity)
        @test length(s.fitted) == 30 * 29 ÷ 2
    end

    @testset "stress on noisy data, ties, positional API" begin
        rng = StableRNG(5)
        Y = rand(rng, 0:5, 20, 8)                  # counts → Bray-Curtis with ties
        D = pairwise(BrayCurtis(), Y; dims=1)
        r = nmds(D; tries=10, rng=StableRNG(6))
        @test 0 < r.stress < 0.3
        X, S = nmds(D, 2, 3)
        @test size(X) == (20, 2) && S isa Float64
        @test size(nmds(D, 3).points) == (20, 3)
        @test occursin("stress", sprint(show, MIME("text/plain"), r))
        @test occursin("explained", sprint(show, MIME("text/plain"), pcoa(D)))
    end

    @testset "input checks" begin
        @test_throws ArgumentError nmds(rand(3, 4))
        @test_throws ArgumentError nmds([0 1 2; 1 0 1; 2 5 0])
        @test_throws ArgumentError nmds([0 -1 2; -1 0 1; 2 1 0])
    end
end
