```@meta
CurrentModule = NonmetricMDS
```

# NonmetricMDS.jl

NonmetricMDS.jl finds a low-dimensional representation of observations from their pairwise dissimilarities. `nmds` uses monotone regression with SMACOF to preserve the order of dissimilarities; `pcoa` provides a metric principal coordinates ordination. The input must be a square, symmetric matrix with non-negative entries. Row and column `i` refer to the same observation.

## From dissimilarities to an ordination

The matrix below contains distances between five points. Use `nmds` to place them in two dimensions:

```julia
using NonmetricMDS, Random

D = [0.0 1.0 1.0 sqrt(2) 2.0;
     1.0 0.0 sqrt(2) 1.0 1.0;
     1.0 sqrt(2) 0.0 1.0 sqrt(5);
     sqrt(2) 1.0 1.0 0.0 sqrt(2);
     2.0 1.0 sqrt(5) sqrt(2) 0.0]

result = nmds(D; k=2, tries=10, rng=MersenneTwister(42))
result.points   # five rows, two coordinate columns
result.stress   # Kruskal's stress-1; lower means a closer fit
```

`tries` controls the number of starts. The first uses a PCoA solution and the others use random configurations. The returned `NMDSResult` also records whether the best start converged, its iteration count, and how many starts reached a similar solution. For metric principal coordinates instead, call `pcoa(D; k=2)` and inspect its `points`, `eigenvalues`, and `explained` fields.

## Shepard plot

[`shepard`](@ref) returns the original dissimilarities, distances in the ordination, and the fitted monotone values, sorted by dissimilarity. For example, with the optional [Plots.jl](https://docs.juliaplots.org/) package installed:

```julia
using Plots

s = shepard(result)
scatter(s.dissimilarity, s.distance;
        xlabel="Dissimilarity", ylabel="Ordination distance",
        label="Pairs", legend=:topleft)
plot!(s.dissimilarity, s.fitted; label="Monotone fit", linewidth=2)
```

The points show how well the ordination preserves each pairwise dissimilarity; the line shows the fitted monotone relationship. `Plots.jl` is needed only to draw the figure.

## API reference

```@docs
NonmetricMDS
nmds
NMDSResult
shepard
pcoa
PCoAResult
procrustes
```
