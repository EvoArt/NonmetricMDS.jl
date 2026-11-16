"""
    procrustes(X, Y; scale=true)

Rotate (and, with `scale=true`, scale) configuration `Y` to best match `X`
after centring both, as `vegan::procrustes(X, Y, symmetric=FALSE)`.
Returns a named tuple `(Yrot, rmse, ss)`: the rotated `Y`, the root mean
squared difference per point, and the residual sum of squares.
"""
function procrustes(X, Y; scale::Bool=true)
    size(X) == size(Y) || throw(DimensionMismatch("configurations differ in size: $(size(X)) vs $(size(Y))"))
    Xc = X .- mean(X; dims=1)
    Yc = Y .- mean(Y; dims=1)
    F = svd(Xc' * Yc)
    A = F.V * F.U'                     # rotation of Y onto X
    c = scale ? sum(F.S) / sum(abs2, Yc) : 1.0
    Yrot = c .* (Yc * A) .+ mean(X; dims=1)
    ss = sum(abs2, (Yrot .- mean(X; dims=1)) .- Xc)
    return (Yrot=Yrot, rmse=sqrt(ss / size(X, 1)), ss=ss)
end
