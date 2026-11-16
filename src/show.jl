# REPL displays

_sty(io, s, color; bold=false) = get(io, :color, false) ?
    string(bold ? "\e[1m" : "", "\e[38;5;$(color)m", s, "\e[0m") : string(s)

function _stress_word(s)
    s < 0.05 ? ("excellent", 114) : s < 0.1 ? ("good", 114) : s < 0.2 ? ("usable", 221) : ("poor", 203)
end

function Base.show(io::IO, r::NMDSResult)
    n, k = size(r.points)
    print(io, "NMDSResult($n points, $k dims, stress = ", @sprintf("%.4f", r.stress), ")")
end

function Base.show(io::IO, ::MIME"text/plain", r::NMDSResult)
    n, k = size(r.points)
    println(io, _sty(io, "NMDSResult", 117; bold=true), _sty(io, " · ", 240), "$n points in $k dimensions")
    word, col = _stress_word(r.stress)
    println(io, _sty(io, "  stress      ", 245), @sprintf("%.4f", r.stress), "  ", _sty(io, "($word)", col))
    rep = r.repeats >= 2 ? _sty(io, "✓", 114) : _sty(io, "only once: try more starts (tries=…)", 221)
    println(io, _sty(io, "  solution    ", 245), "best of $(r.tries) starts, reached $(r.repeats)× ", rep)
    r.converged || println(io, _sty(io, "  warning     ", 245),
                           _sty(io, "best start hit maxit before converging", 221))
    print(io, _sty(io, "  fields      ", 245), "points, stress, shepard(r) for a Shepard plot")
end

function Base.show(io::IO, r::PCoAResult)
    print(io, "PCoAResult($(size(r.points, 1)) points, $(size(r.points, 2)) axes)")
end

function Base.show(io::IO, ::MIME"text/plain", r::PCoAResult)
    n, k = size(r.points)
    println(io, _sty(io, "PCoAResult", 117; bold=true), _sty(io, " · ", 240), "$n points, $k axes")
    ex = join((@sprintf("PCo%d %.1f%%", i, 100e) for (i, e) in enumerate(r.explained)), "  ")
    println(io, _sty(io, "  explained   ", 245), ex)
    nneg = count(<(-sqrt(eps()) * abs(r.eigenvalues[1])), r.eigenvalues)
    nneg > 0 && println(io, _sty(io, "  note        ", 245),
                        "$nneg negative eigenvalues (the dissimilarity is not Euclidean)")
    print(io, _sty(io, "  fields      ", 245), "points, eigenvalues, explained")
end
