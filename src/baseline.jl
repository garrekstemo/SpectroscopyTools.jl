"""
Baseline correction algorithms for spectroscopy.

Provides:
- `als_baseline` — Asymmetric Least Squares (Eilers & Boelens, 2005)
- `arpls_baseline` — Asymmetrically Reweighted PLS (Baek et al., 2015)
- `snip_baseline` — Statistics-sensitive Non-linear Iterative Peak-clipping (Ryan et al., 1988)
"""

using SparseArrays
using LinearAlgebra

# =============================================================================
# Asymmetric Least Squares (ALS)
# =============================================================================

"""
    als_baseline(y; λ=1e5, p=0.01, maxiter=10, tol=1e-6) -> Vector

Asymmetric Least Squares baseline correction.
"""
function als_baseline(y::AbstractVector{<:Real};
                      λ::Real=1e5,
                      p::Real=0.01,
                      maxiter::Int=10,
                      tol::Real=1e-6)
    n = length(y)

    0 < p < 1 || throw(ArgumentError("p must be in (0, 1), got $p"))
    λ > 0 || throw(ArgumentError("λ must be positive, got $λ"))
    n ≥ 3 || throw(ArgumentError("Need at least 3 points, got $n"))

    D = _diff_matrix(n, 2)
    DtD = D' * D

    w = ones(n)
    z = similar(y, Float64)
    z_prev = similar(z)

    for iter in 1:maxiter
        copyto!(z_prev, z)

        W = spdiagm(0 => w)
        z .= (W + λ * DtD) \ (w .* y)

        for i in eachindex(w)
            w[i] = y[i] > z[i] ? p : (1 - p)
        end

        if iter > 1
            δ = norm(z - z_prev) / (norm(z) + eps())
            δ < tol && break
        end
    end

    return z
end

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
    baseline = if method == :als
        als_baseline(y; kwargs...)
    elseif method == :arpls
        arpls_baseline(y; kwargs...)
    elseif method == :snip
        snip_baseline(y; kwargs...)
    else
        available = (:als, :arpls, :snip)
        throw(ArgumentError("Unknown method :$method. Available: $available"))
    end

    corrected = y - baseline
    return (y=corrected, baseline=baseline)
end

"""
    correct_baseline(x, y; kwargs...) -> NamedTuple

Version that also returns x values for convenience.
"""
function correct_baseline(x::AbstractVector, y::AbstractVector{<:Real}; kwargs...)
    result = correct_baseline(y; kwargs...)
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
