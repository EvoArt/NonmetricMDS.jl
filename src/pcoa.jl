"""
    PCoAResult

Principal coordinates of a dissimilarity matrix.

- `points`: one row per observation, one column per axis kept
- `eigenvalues`: all eigenvalues of the double-centred matrix, largest first
  (negative ones mean the dissimilarity is not Euclidean)
- `explained`: share of the total positive eigenvalue sum on each kept axis
- `labels`: observation names (may be empty)
"""
struct PCoAResult
    points::Matrix{Float64}
    eigenvalues::Vector{Float64}
    explained::Vector{Float64}
    labels::Vector{String}
end

"""
    pcoa(D; k=2, labels=String[])

Principal coordinates analysis (classical, metric MDS) of the symmetric
dissimilarity matrix `D`, keeping `k` axes (at most the number of positive
eigenvalues). The same as `cmdscale(D, k, eig=TRUE)` in R, up to the sign of
each axis.
"""
function pcoa(D; k::Integer=2, labels=String[])
    n = _check_dissimilarity(D)
    G = Symmetric(gower_center(D))
    F = eigen(G)
    ord = sortperm(F.values; rev=true)
    vals = F.values[ord]
    vecs = F.vectors[:, ord]
    tol = sqrt(eps()) * max(abs(vals[1]), 1.0)
    npos = count(>(tol), vals)
    k = min(k, npos)
    pts = vecs[:, 1:k] .* sqrt.(vals[1:k])'
    _fix_signs!(pts)
    expl = vals[1:k] ./ sum(v for v in vals if v > tol)
    return PCoAResult(pts, vals, expl, string.(collect(labels)))
end

# -1/2 J D² J
function gower_center(D)
    A = -0.5 .* float.(D) .^ 2
    A .-= mean(A; dims=1)
    A .-= mean(A; dims=2)
    return A
end

# make the largest-magnitude coordinate of each axis positive, so results
# do not flip between runs or platforms
function _fix_signs!(X)
    for j in axes(X, 2)
        i = argmax(abs.(view(X, :, j)))
        X[i, j] < 0 && (X[:, j] .*= -1)
    end
    return X
end

function _check_dissimilarity(D)
    n = size(D, 1)
    size(D, 2) == n || throw(ArgumentError("dissimilarity matrix must be square, got $(size(D))"))
    n >= 3 || throw(ArgumentError("need at least 3 observations, got $n"))
    any(x -> ismissing(x) || isnan(x), D) && throw(ArgumentError("dissimilarity matrix has missing or NaN values"))
    any(<(0), D) && throw(ArgumentError("dissimilarities must be non-negative"))
    for j in 1:n, i in j+1:n
        isapprox(D[i, j], D[j, i]; rtol=1e-8, atol=1e-12) ||
            throw(ArgumentError("dissimilarity matrix is not symmetric (D[$i,$j] ≠ D[$j,$i])"))
    end
    return n
end
