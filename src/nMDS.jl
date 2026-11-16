"""
    nMDS

Ordination of a dissimilarity matrix: non-metric multidimensional scaling
([`nmds`](@ref)) and principal coordinates analysis ([`pcoa`](@ref)).

`nmds` minimises Kruskal's stress-1 by SMACOF majorisation with monotone
(isotonic) regression, from several starts, as `vegan::metaMDS` does with
`monoMDS`. The solution is centred and rotated to its principal axes.
"""
module nMDS

using LinearAlgebra, Random, Statistics
using Printf: @sprintf

export nmds, pcoa, NMDSResult, PCoAResult, procrustes, shepard

include("pava.jl")
include("pcoa.jl")
include("procrustes.jl")
include("smacof.jl")
include("show.jl")

end
