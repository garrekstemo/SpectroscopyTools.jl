"""
Baseline correction algorithms for spectroscopy.

Provides:
- `arpls_baseline` — Asymmetrically Reweighted PLS (Baek et al., 2015)
- `snip_baseline` — Statistics-sensitive Non-linear Iterative Peak-clipping (Ryan et al., 1988)
- `rubberband_baseline` — Lower convex hull (rubber band) baseline
"""

# =============================================================================
# Asymmetrically Reweighted Penalized Least Squares (arPLS)
# =============================================================================

"""
    arpls_baseline(y; λ=1e5, maxiter=100, tol=1e-6) -> Vector

Asymmetrically Reweighted Penalized Least Squares baseline correction.
"""
function arpls_baseline(y::AbstractVector{<:Real};
                        λ::Real=1e5,
                        maxiter::Int=100,
                        tol::Real=1e-6)
    n = length(y)

    λ > 0 || throw(ArgumentError("λ must be positive, got $λ"))
    n ≥ 3 || throw(ArgumentError("Need at least 3 points, got $n"))

    D = _diff_matrix(n, 2)
    DtD = D' * D

    w = ones(n)
    z = similar(y, Float64)

    for iter in 1:maxiter
        w_prev = copy(w)

        W = spdiagm(0 => w)
        z .= (W + λ * DtD) \ (w .* y)

        d = y - z

        d_neg = d[d .< 0]
        if isempty(d_neg)
            break
        end

        m = mean(d_neg)
        σ = std(d_neg; corrected=false)

        if σ < eps()
            break
        end

        threshold = -m + 2σ
        for i in eachindex(w)
            w[i] = 1 / (1 + exp(2 * (d[i] - threshold) / σ))
        end

        δ = norm(w - w_prev) / (norm(w) + eps())
        δ < tol && break
    end

    return z
end

# =============================================================================
# SNIP
# =============================================================================

"""
    snip_baseline(y; iterations=40, decreasing=true) -> Vector

SNIP baseline correction using iterative peak clipping.
"""
function snip_baseline(y::AbstractVector{<:Real};
                       iterations::Int=40,
                       decreasing::Bool=true)
    n = length(y)

    iterations > 0 || throw(ArgumentError("iterations must be positive"))
    n ≥ 3 || throw(ArgumentError("Need at least 3 points, got $n"))

    z = convert(Vector{Float64}, copy(y))

    _lls_transform!(z)

    window_sizes = decreasing ? (iterations:-1:1) : (1:iterations)

    for k in window_sizes
        for i in (k+1):(n-k)
            neighbor_avg = (z[i-k] + z[i+k]) / 2
            z[i] = min(z[i], neighbor_avg)
        end
    end

    _lls_inverse!(z)

    return z
end

function _lls_transform!(z::Vector{Float64})
    for i in eachindex(z)
        v = max(z[i], 0.0)
        z[i] = log(log(sqrt(v + 1) + 1) + 1)
    end
end

function _lls_inverse!(z::Vector{Float64})
    for i in eachindex(z)
        v = z[i]
        z[i] = (exp(exp(v) - 1) - 1)^2 - 1
    end
end

# =============================================================================
# Rubber band baseline
# =============================================================================

"""
    rubberband_baseline(x, y)

Rubber band baseline correction using the lower convex hull.

Equivalent to stretching a rubber band under the spectrum -- the baseline
follows the lower envelope. Used in OPUS for FTIR baseline correction.

Returns the baseline y-values (same length as input).
"""
function rubberband_baseline(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have the same length"))
    n >= 2 || throw(ArgumentError("Need at least 2 points"))

    xv = Float64.(x)
    yv = Float64.(y)

    # Compute lower convex hull indices using Andrew's monotone chain
    # We want the lower hull of points (x[i], y[i])
    # Sort by x, then build lower hull
    order = sortperm(xv)
    hull = Int[]  # indices into the original arrays (via order)

    for idx in order
        while length(hull) >= 2
            h1 = hull[end-1]
            h2 = hull[end]
            # Cross product to check left turn (lower hull wants right turns)
            cross = (xv[h2] - xv[h1]) * (yv[idx] - yv[h1]) -
                    (yv[h2] - yv[h1]) * (xv[idx] - xv[h1])
            if cross <= 0
                pop!(hull)
            else
                break
            end
        end
        push!(hull, idx)
    end

    # Interpolate hull back to original x-grid
    hull_x = xv[hull]
    hull_y = yv[hull]
    itp = Interpolations.linear_interpolation(hull_x, hull_y, extrapolation_bc=Interpolations.Flat())
    return itp.(xv)
end

# =============================================================================
# Unified API
# =============================================================================

"""
    correct_baseline(y; method=:arpls, kwargs...) -> NamedTuple

Correct baseline and return both corrected spectrum and baseline.

Returns `(y=corrected, baseline=baseline)`.
"""
function correct_baseline(y::AbstractVector{<:Real};
                          method::Symbol=:arpls,
                          kwargs...)
    baseline = if method == :arpls
        arpls_baseline(y; kwargs...)
    elseif method == :snip
        snip_baseline(y; kwargs...)
    elseif method == :rubberband
        # rubberband_baseline needs x-values; create a dummy index grid
        rubberband_baseline(collect(1.0:length(y)), y)
    else
        available = (:arpls, :snip, :rubberband)
        throw(ArgumentError("Unknown method :$method. Available: $available"))
    end

    corrected = y - baseline
    return (y=corrected, baseline=baseline)
end

"""
    correct_baseline(x, y; kwargs...) -> NamedTuple

Version that also returns x values for convenience.
Methods that use x-values (rubberband) receive real x, not dummy indices.
"""
function correct_baseline(x::AbstractVector, y::AbstractVector{<:Real};
                          method::Symbol=:arpls, kwargs...)
    if method == :rubberband
        baseline = rubberband_baseline(x, y)
        corrected = y - baseline
        return (x=collect(x), y=corrected, baseline=baseline)
    end
    result = correct_baseline(y; method=method, kwargs...)
    return (x=collect(x), y=result.y, baseline=result.baseline)
end

# =============================================================================
# Internal utilities
# =============================================================================

function _diff_matrix(n::Int, d::Int)
    D = sparse(1.0I, n, n)

    for _ in 1:d
        m = size(D, 1)
        D = D[2:m, :] - D[1:m-1, :]
    end

    return D
end
