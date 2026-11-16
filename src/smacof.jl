"""
    NMDSResult

A non-metric MDS solution.

- `points`: one row per observation, centred and rotated to principal axes
- `stress`: Kruskal's stress-1 of the best solution
- `tries`, `repeats`: number of starts, and how many reached the best
  solution again (same stress within 1e-4 and Procrustes RMSE < 0.01)
- `converged`: whether the best start met the stress tolerance
- `dissimilarities`, `distances`, `fitted`: per pair (lower triangle,
  column-major), for a Shepard plot (see [`shepard`](@ref))
- `labels`: observation names (may be empty)

`X, S = nmds(D, 2)` still works: iterating gives `points` then `stress`.
"""
struct NMDSResult
    points::Matrix{Float64}
    stress::Float64
    tries::Int
    repeats::Int
    converged::Bool
    iterations::Int
    dissimilarities::Vector{Float64}
    distances::Vector{Float64}
    fitted::Vector{Float64}
    labels::Vector{String}
end

Base.iterate(r::NMDSResult, s=1) = s == 1 ? (r.points, 2) : s == 2 ? (r.stress, 3) : nothing
Base.length(::NMDSResult) = 2

"""
    nmds(D; k=2, tries=20, maxit=500, tol=1e-7, rng=Random.default_rng(),
         labels=String[], threaded=true)

Non-metric multidimensional scaling of the symmetric dissimilarity matrix `D`
in `k` dimensions.

Each start minimises Kruskal's stress-1 by SMACOF majorisation, refitting a
monotone regression of distances on dissimilarities at every step (ties in
`D` are treated by the primary approach, as `monoMDS` does by default).
The first start is the PCoA solution; the others are random. The best
solution is centred and rotated to its principal components.

Stress below 0.05 is an excellent representation, below 0.1 good, below 0.2
usable, above 0.2 poor (Kruskal's rule of thumb; Clarke 1993).

Positional forms of earlier versions still work: `nmds(D, k)` and
`nmds(D, k, tries)`.
"""
function nmds(D; k::Integer=2, tries::Integer=20, maxit::Integer=500, tol::Real=1e-7,
              rng::AbstractRNG=Random.default_rng(), labels=String[], threaded::Bool=true)
    n = _check_dissimilarity(D)
    k >= 1 || throw(ArgumentError("k must be at least 1"))
    k < n || throw(ArgumentError("k = $k dimensions need more than $k observations"))
    tries >= 1 || throw(ArgumentError("tries must be at least 1"))

    prob = _problem(D)
    starts = Vector{Matrix{Float64}}(undef, tries)
    starts[1] = _pcoa_start(D, k, n, rng)
    for t in 2:tries
        starts[t] = rand(rng, n, k) .- 0.5
    end
    runs = Vector{Any}(undef, tries)
    if threaded && Threads.nthreads() > 1
        Threads.@threads for t in 1:tries
            runs[t] = _smacof(prob, starts[t], maxit, tol)
        end
    else
        for t in 1:tries
            runs[t] = _smacof(prob, starts[t], maxit, tol)
        end
    end

    best = argmin([r.stress for r in runs])
    b = runs[best]
    X = _principal_axes(b.X)
    repeats = count(runs) do r
        abs(r.stress - b.stress) < 1e-4 && procrustes(X, r.X).rmse < 0.01 * _scale(X)
    end
    d = _pair_distances(X, prob)
    dhat = similar(d)
    monotone_fit!(dhat, d, copy(prob.order), prob.blocks, similar(d), similar(d))
    return NMDSResult(X, b.stress, tries, repeats, b.converged, b.iterations,
                      prob.delta, d, dhat, string.(collect(labels)))
end

nmds(D, k::Integer; kw...) = nmds(D; k=k, kw...)
nmds(D, k::Integer, tries::Integer, args...; kw...) = nmds(D; k=k, tries=tries, kw...)

_scale(X) = sqrt(sum(abs2, X .- mean(X; dims=1)) / size(X, 1))

# pairs (i > j) in column-major order of the lower triangle, their
# dissimilarities, and the dissimilarity order with its tied runs
function _problem(D)
    n = size(D, 1)
    m = n * (n - 1) ÷ 2
    I = Vector{Int}(undef, m)
    J = Vector{Int}(undef, m)
    delta = Vector{Float64}(undef, m)
    p = 0
    for j in 1:n-1, i in j+1:n
        p += 1
        I[p], J[p], delta[p] = i, j, D[i, j]
    end
    order = sortperm(delta)
    blocks = tie_blocks(delta[order])
    return (n=n, I=I, J=J, delta=delta, order=order, blocks=blocks)
end

function _pcoa_start(D, k, n, rng)
    p = pcoa(D; k=k)
    X = p.points
    size(X, 2) < k && (X = hcat(X, 1e-3 .* randn(rng, n, k - size(X, 2))))
    return X
end

function _pair_distances(X, prob)
    d = Vector{Float64}(undef, length(prob.I))
    _distances!(d, X, prob.I, prob.J)
    return d
end

function _distances!(d, X, I, J)
    k = size(X, 2)
    @inbounds for p in eachindex(d)
        i, j = I[p], J[p]
        s = 0.0
        for c in 1:k
            s += (X[i, c] - X[j, c])^2
        end
        d[p] = sqrt(s)
    end
    return d
end

stress1(d, dhat) = sqrt(sum(abs2, d .- dhat) / sum(abs2, d))

# One start of non-metric SMACOF (de Leeuw 1977; Borg & Groenen 2005, ch. 9).
# Disparities are rescaled so their sum of squares is the number of pairs;
# the Guttman transform then lowers raw stress for fixed disparities.
function _smacof(prob, X0, maxit, tol)
    n, I, J = prob.n, prob.I, prob.J
    m = length(I)
    X = copy(X0)
    Xn = similar(X)
    k = size(X, 2)
    d = Vector{Float64}(undef, m)
    dhat = similar(d)
    ybuf, fbuf = similar(d), similar(d)
    order = copy(prob.order)
    _distances!(d, X, I, J)
    X .*= sqrt(m / sum(abs2, d))      # unit mean squared distance
    s_old = Inf
    s = Inf
    converged = false
    it = 0
    while it < maxit
        it += 1
        _distances!(d, X, I, J)
        monotone_fit!(dhat, d, order, prob.blocks, ybuf, fbuf)
        s = stress1(d, dhat)
        if s < 1e-6 || s_old - s < tol * s_old
            converged = true
            break
        end
        s_old = s
        dhat .*= sqrt(m / sum(abs2, dhat))
        # Guttman transform: X ← B(X) X / n
        fill!(Xn, 0.0)
        @inbounds for p in 1:m
            d[p] > 0 || continue
            b = dhat[p] / d[p]
            i, j = I[p], J[p]
            for c in 1:k
                t = b * (X[i, c] - X[j, c])
                Xn[i, c] += t
                Xn[j, c] -= t
            end
        end
        X .= Xn ./ n
    end
    return (X=X, stress=s, converged=converged, iterations=it)
end

# centre and rotate to principal components (vegan::postMDS, pc = TRUE)
function _principal_axes(X)
    Xc = X .- mean(X; dims=1)
    F = svd(Xc)
    Y = Xc * F.V
    return _fix_signs!(Y)
end

"""
    shepard(r::NMDSResult)

Data for a Shepard diagram: a named tuple of the original dissimilarities,
the ordination distances and the monotone-regression fit, one entry per
pair of observations, sorted by dissimilarity.
"""
function shepard(r::NMDSResult)
    o = sortperm(r.dissimilarities)
    return (dissimilarity=r.dissimilarities[o], distance=r.distances[o], fitted=r.fitted[o])
end
