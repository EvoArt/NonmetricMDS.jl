# NonmetricMDS.jl

[![Dev docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://EvoArt.github.io/NonmetricMDS.jl/dev/)
[![CI](https://github.com/EvoArt/NonmetricMDS.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/EvoArt/NonmetricMDS.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/EvoArt/NonmetricMDS.jl/graph/badge.svg)](https://codecov.io/gh/EvoArt/NonmetricMDS.jl)

NonmetricMDS.jl turns a square, symmetric dissimilarity matrix into a low-dimensional ordination. It provides non-metric multidimensional scaling (`nmds`) using monotone regression and SMACOF, principal coordinates analysis (`pcoa`), Procrustes alignment, and data for Shepard plots. Each row of the returned coordinates represents one observation.

Install it with Julia's package manager (`] add NonmetricMDS`) once the package is registered. To use the current repository before registration, run `] add https://github.com/EvoArt/NonmetricMDS.jl`.

```julia
using NonmetricMDS, Random

# Pairwise dissimilarities between five observations.
D = [0.0 1.0 1.0 sqrt(2) 2.0;
     1.0 0.0 sqrt(2) 1.0 1.0;
     1.0 sqrt(2) 0.0 1.0 sqrt(5);
     sqrt(2) 1.0 1.0 0.0 sqrt(2);
     2.0 1.0 sqrt(5) sqrt(2) 0.0]

result = nmds(D; k=2, tries=10, rng=MersenneTwister(42))
result.points   # 5 × 2 coordinates, one row per observation
result.stress   # Kruskal's stress-1; smaller is better

s = shepard(result)  # dissimilarity, distance, and fitted vectors
```

For a Shepard plot, plot `s.distance` against `s.dissimilarity` and overlay `s.fitted`. A complete plotting example and the [API reference](https://EvoArt.github.io/NonmetricMDS.jl/dev/) are in the documentation. Plotting is optional; the package does not require a plotting library.

[![Schematic of non-metric MDS](docs/schematic.PNG)](https://link.springer.com/content/pdf/10.1007/BF02289694.pdf)
