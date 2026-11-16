# Monotone (isotonic) regression by pool-adjacent-violators.

"""
    isotonic!(fit, y)

Least-squares non-decreasing fit to `y` (unit weights), written to `fit`.
"""
function isotonic!(fit, y)
    n = length(y)
    n == 0 && return fit
    # blocks as a stack of (sum, count); merged while they violate the order
    sums = Vector{Float64}(undef, n)
    cnts = Vector{Int}(undef, n)
    top = 0
    for i in 1:n
        top += 1
        sums[top] = y[i]
        cnts[top] = 1
        while top > 1 && sums[top-1] / cnts[top-1] >= sums[top] / cnts[top]
            sums[top-1] += sums[top]
            cnts[top-1] += cnts[top]
            top -= 1
        end
    end
    i = 1
    for b in 1:top
        v = sums[b] / cnts[b]
        for _ in 1:cnts[b]
            fit[i] = v
            i += 1
        end
    end
    return fit
end

isotonic(y) = isotonic!(similar(y, Float64), y)

# Kruskal's primary approach to ties: within each run of equal
# dissimilarities, order the pairs by their current distance before the
# monotone regression, so tied dissimilarities need not get equal fits.
# `order` sorts pairs by dissimilarity; `blocks` are the tied runs in it.
function monotone_fit!(dhat, d, order, blocks, ybuf, fbuf)
    for r in blocks
        length(r) > 1 && sort!(view(order, r); by=i -> d[i])
    end
    @inbounds for (k, i) in enumerate(order)
        ybuf[k] = d[i]
    end
    isotonic!(fbuf, ybuf)
    @inbounds for (k, i) in enumerate(order)
        dhat[i] = fbuf[k]
    end
    return dhat
end

# runs of equal values in sorted `x`
function tie_blocks(x)
    blocks = UnitRange{Int}[]
    i = 1
    n = length(x)
    while i <= n
        j = i
        while j < n && x[j+1] == x[i]
            j += 1
        end
        push!(blocks, i:j)
        i = j + 1
    end
    return blocks
end
