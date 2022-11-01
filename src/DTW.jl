import ProgressMeter
using Statistics
using StatsBase
using Distances
using FillArrays
using ProgressMeter
using DataStructures
using SlidingDistancesBase
import SlidingDistancesBase: floattype, lastlength, setup_normalizer
using LoopVectorization
using Plots, Plots.PlotMeasures # Plots required both for @layout and for the margins
using Requires
using UnPack


#### utils.jl ####

@inline function indmin3(a,b,c,i,j)
    if a <= b
        if a <= c
            return 1,i-1,j-1
        else
            return 3,i,j-1
        end
    else
        if b <= c
            return 2,i-1,j
        else
            return 3,i,j-1
        end
    end
end


"""
    imin, imax = radiuslimits(r, n::Int, m::Int)
    imin, imax = radiuslimits(r, seq1, seq2)
"""
function radiuslimits(r, n::Int, m::Int)
    d = abs(m-n)
    if m >= n
        imin, imax = max.((1:n) .- r,1), min.((1:n) .+ (r+d), m)
    else
        imin, imax = max.((1:n) .- (r+d),1), min.((1:n) .+ r, m)
    end

    imin, imax
end

radiuslimits(r,seq1, seq2) = radiuslimits(r, lastlength(seq1), lastlength(seq2))

"""
    inds = align_signals(s::AbstractVector{<:AbstractArray}, master = argmax(length.(s)); method = :dtw, output=:indices)

Compute a set of indices such that `s[i][inds[i]]` is optimally aligned to `s[master]`. All elements of `inds` will have the same length.

# Arguments:
- `s`: A vector of signals to align
- `master`: Index of the signal used as reference. All the other signals will be aligned to this one.
- `method`: `:dtw` uses the warping paths from dtw between `s[master]` and `s[i]`. `:xcorr` uses `DSP.finddelay` which internally computes the cross correlation between signals, which often results in a slight misalignment.  
- `output`: The default is to output the aligning indices, alternatively, the aligned `:signals` themselves can be outputted.
"""
function align_signals(s::AbstractVector{<:AbstractArray}, master::Integer=argmax(length.(s)); method=:dtw, output=:indices, postprocess=nothing)
    inds = [1:lastlength(s) for s in s]
    # find delays to align with master
    d = map(s) do si
        si === s[master] && (return 0)
        if method ∈ (:xcorr, :crosscorr, :dsp)
            SlidingDistancesBase.DSP.finddelay(s[master], si) # suboptimal because xcorr does not do exactly what we want
        elseif method ∈ (:dtw, :DTW)
            d,i1,i2 = dtw(si, s[master], postprocess=postprocess)
            round(Int, median(i2-i1))
        else
            throw(ArgumentError("Unknown method"))
        end
    end

    # find left and right (virtual) zero padding
    lp = maximum(d)
    rp = maximum(length(s[master]) .- (length.(s) .+ d))
    
    # New window length
    wl = length(inds[master]) - lp - rp
    
    # trim individual index sets to fit into new master window
    for i in eachindex(inds)
        start = max(1, 1+lp-d[i])
        stop = min(length(s[i]),start+wl-1)
        inds[i] = start : stop
    end
    @assert all(length.(inds) .== length(inds[1]))
    if output === :indices
        return inds 
    else
        return [s[!,i] for (s,i) in zip(s,inds)]
    end
end



#### distance_interface.jl ####

# methods for estimating dtw #

abstract type DTWDistance{D <: Union{Function, Distances.PreMetric}} <: Distances.SemiMetric end


"""
    struct DTW{D,N} <: DTWDistance{D}

# Keyword arguments:
- `radius`: The maximum allowed deviation of the matching path from the diagonal
- `dist`: Inner distance
- `transportcost` If >1, an additional penalty factor for non-diagonal moves is added.
- `normalizer`: defaults to `Nothing`

If the two time-series are of equal length, [`dtw_cost`](@ref) is called, if not, [`dtwnn`](@ref) is called.
"""
struct DTW{D,N,FK} <: DTWDistance{D}
    "The maximum allowed deviation of the matching path from the diagonal"
    radius::Int
    dist::D
    "If >1, an additional penalty factor for non-diagonal moves is added."
    transportcost::Float64
    postprocess::FK
end
DTW(r,dist::D=SqEuclidean(); transportcost=1, postprocess::FK=nothing, normalizer::Type{N}=Nothing) where {D,N,FK} = DTW{D, N, FK}(r,dist,transportcost,postprocess)
DTW(;radius,dist=SqEuclidean(), transportcost=1, postprocess::FK=nothing, normalizer::Type{N}=Nothing) where {N,FK} = DTW{typeof(dist), N, FK}(radius,dist,transportcost,postprocess)

"""
    struct SoftDTW{D, T} <: DTWDistance{D}

# Arguments:
- `γ`: smoothing parameter default to 1. Smaller value makes the distance closer to standard DTW.
- `dist`: Inner distance, defaults to `SqEuclidean()`.
- `transportcost`
"""
Base.@kwdef struct SoftDTW{D,T,R} <: DTWDistance{D}
    γ::T
    "The maximum allowed deviation of the matching path from the diagonal"
    dist::D = SqEuclidean()
    "If >1, an additional penalty factor for non-diagonal moves is added."
    transportcost::Float64 = 1.0
    radius::R = nothing
    SoftDTW(γ, dist=SqEuclidean(),transportcost=1,radius=nothing) = new{typeof(dist), typeof(γ), typeof(radius)}(γ,dist,transportcost,radius)
end


"""
    struct FastDTW{D} <: DTWDistance{D}

- `radius`
- `dist` inner distance
"""
Base.@kwdef struct FastDTW{D} <: DTWDistance{D}
    radius::Int
    dist::D = SqEuclidean()
    FastDTW(r, dist=SqEuclidean()) = new{typeof(dist)}(r, dist)
end



function Distances.evaluate(d::DTW{<:Any,N}, x, y; kwargs...) where N
    if lastlength(x) == lastlength(y)
        x, y = normalize(N, x), normalize(N, y)
        return dtw_cost(x, y, d.dist, d.radius; transportcost = d.transportcost, kwargs...)
    end
    if lastlength(x) > lastlength(y)
        x, y = y, x
    end
    dtwnn(
        x,
        y,
        d.dist,
        d.radius,
        N;
        prune_endpoints = false,
        transportcost = d.transportcost,
        kwargs...,
    ).cost
end

Distances.evaluate(d::SoftDTW, x, y) = soft_dtw_cost(x, y, d.dist, γ=d.γ, radius=d.radius)
Distances.evaluate(d::FastDTW, x, y) =
    fastdtw(x, y, d.dist, d.radius)[1]

distpath(d::DTW, x, y) = dtw(x, y, d.dist; transportcost=d.transportcost, postprocess = d.postprocess)
distpath(d::DTW, x, y, i2min::AbstractVector, i2max::AbstractVector; transportcost=d.transportcost) =
    dtw(x, y, i2min, i2max, d.dist)
distpath(d::FastDTW, x, y) = fastdtw(x, y, d.dist, d.radius)

(d::DTWDistance)(x,y;kwargs...) = Distances.evaluate(d,x,y;kwargs...)

"""
    distance_profile(d::DTWDistance, Q::AbstractArray{S}, T::AbstractArray{S}; kwargs...) where S

Optimized method for computing the distance profile using DTW distances. kwargs are sent to [`dtwnn`](@ref).
"""
function SlidingDistancesBase.distance_profile(d::DTW{<:Any,N}, Q::AbstractArray{S}, T::AbstractArray{S}; kwargs...) where {N,S}
    m = lastlength(Q)
    n = lastlength(T)
    n >= m || throw(ArgumentError("Q cannot be longer than T"))
    l = n-m+1
    res = dtwnn(Q, T, d.dist, d.radius, N; saveall=true, kwargs...)
    res.dists
end

#### filters.jl ####

function imfilter(A, kern)
    size(kern,1) == size(kern,2) || throw(ArgumentError("Only square kernels supported, you provided $(size(kern))."))
    isodd(size(kern,1)) || throw(ArgumentError("Only odd-sized kernels supported, you provided $(size(kern))."))
    A2 = LocalFilters.convolve(A, kern)
    pad = size(kern, 1) 
    @views A2[1:pad, :] .= A[1:pad, :]
    @views A2[end-pad+1:end, :] .= A[end-pad+1:end, :]
    @views A2[:, 1:pad] .= A[:, 1:pad]
    @views A2[:, end-pad+1:end] .= A[:, end-pad+1:end]
    A2
end


function gaussian2(n)
    t = range(-2,stop=2, length=n)
    k = exp.(-t.^2)
    K = k*k'
    K ./= sum(K)
end

function gaussian(n)
    t = range(-2,stop=2, length=n)
    K = [exp(-(x-y)^2) for x in t, y in t]
    K ./= sum(K)
end

#### dtw.jl ####

#####################################
#     Basic interface functions     #
#####################################

"""
    cost,i1,i2 = dtw(seq1, seq2, [dist=SqEuclidean, postprocess=nothing])
    cost,i1,i2 = dtw(seq1, seq2, dist, i2min, i2max)

Perform dynamic-time warping to measure the distance between two sequences.

Find a set of indices (`i1`,`i2`) that align two series (`seq1`,`seq2`) by
dynamic axis warping. Also returns the distance (after warping) according to
the SemiMetric `dist`, which defaults to squared Euclidean distance (see
Distances.jl). If `seq1` and `seq2` are matrices, each column is considered
an observation.

If `i2min/max` are provided, do DTW to align `seq1` and `seq2` confined to a window. Vectors `i2min` and
`i2max` specify (inclusive) lower and upper bounds for `seq2` for each index in
`seq1`. Thus, `i2min` and `i2max` are required to be the same length as `seq1`.

If `filternernel::AbstractMatrix` is provided, it's used to filter the cost matrix. Create a suitable kerlen using, e.g., `ImageFiltering.Kernel.gaussian(3)`. The filtering of the cost matrix makes the warping smoother, effectively penalizing small-scale warping.

See also [`dtw_cost`](@ref) and [`dtwnn`](@ref).
"""
function dtw(args...; kwargs...)
    D = dtw_cost_matrix(args...; kwargs...)
    return trackback(D)
end

##############################
#  Cost matrix computations  #
##############################

Distances.pairwise(d::PreMetric, s1::AbstractVector, s2::AbstractVector; dims=2) = evaluate.(Ref(d), s1, transpose(s2))

function Distances.pairwise(d::PreMetric, s1::AbstractArray, s2::AbstractArray; dims=2)
    [evaluate(d, s1[!,i], s2[!,j]) for i in 1:lastlength(s1), j in lastlength(s2)]
end

@inbounds function dtw_cost_matrix(seq1::AbstractArray{T}, seq2::AbstractArray{T}, dist::SemiMetric = SqEuclidean();
    transportcost=1,
    postprocess = nothing) where T
    # Build the cost matrix
    m = lastlength(seq2)
    n = lastlength(seq1)

    # Initialize first column and first row
    D = Distances.pairwise(dist, seq2, seq1, dims=2)
    @assert size(D) == (m,n)

    for r=2:m
        D[r,1] += D[r-1,1]
    end
    for c=2:n
        D[1,c] += D[1,c-1]
    end

    # Complete the cost matrix
    for c = 2:n
        for r = 2:m
            best_neighbor_cost = min(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1])
            D[r, c] += best_neighbor_cost
        end
    end

    if postprocess !== nothing
        D = postprocess(D)
    end

    return D
end


Base.@propagate_inbounds function dtw_cost_matrix(
    seq1::AbstractArray{T},
    seq2::AbstractArray{T},
    dist::SemiMetric,
    i2min::AbstractVector{U},
    i2max::AbstractVector{U};
    transportcost = 1
) where {T,U<:Integer}
    n = lastlength(seq1) # of columns in cost matrix
    m = lastlength(seq2) # of rows in cost matrix
    Base.@boundscheck begin
        n == length(i2min) || throw(ArgumentError("i2min does not match length of seq1."))
        n == length(i2max) || throw(ArgumentError("i2max does not match length of seq1."))
        1 == i2min[1]      || throw(ArgumentError("i2min must start at 1."))
        m == i2max[end]    || throw(ArgumentError("i2max must end at length(seq2), was $(i2max[end]) ≂̸ $(m)"))
    end

    # Build the (n x m) cost matrix into a WindowedMatrix, because it's ragged.
    # That type gives efficient storage with convenient [r,c] indexing and returns
    # Inf when accessed outside the window.
    D = WindowedMatrix(i2min, i2max, Inf)

    # First column first
    D[1, 1] = evaluate(dist, seq1[!, 1], seq2[!, 1])
    for r = 2:i2max[1]
        D[r, 1] = D[r-1, 1] + evaluate(dist, seq1[!, 1], seq2[!, r])
    end

    # Complete the cost matrix from columns 2 to m.
    for c = 2:n
        for r = i2min[c]:i2max[c]
            best_neighbor_cost = min(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1])
            D[r, c] = best_neighbor_cost + evaluate(dist, seq1[!, c], seq2[!, r])
        end
    end

    return D
end

########################################
#  Find Best Path through Cost Matrix  #
########################################

"""
    cost,cols,rows = trackback(D::Matrix)

Given the cost matrix `D`, computes the optimal track from end to beginning.
Returns `cols` and `rows` which are vectors respectively holding the track.
"""
function trackback(D::AbstractMatrix{T}) where {T<:Number}

    # initialize trackback throught rows/columns
    r, c       = size(D)
    rows, cols = Int[r], Int[c]

    # estimate that we'll need N⋅logN elements
    N  = max(r, c)
    sz = 2 * N
    sizehint!(rows, sz)
    sizehint!(cols, sz)

    # do trackback
    @inbounds while r > 1 && c > 1
        tb, r, c = indmin3(D[r-1, c-1], D[r-1, c], D[r, c-1], r, c)
        push!(rows, r)
        push!(cols, c)
    end
    # Possibly either r>1 or c>1 at this point (but not both).
    # Add the unfinished part of the track to reach [1,1]
    for r = r-1:-1:1
        push!(rows, r)
        push!(cols, 1)
    end
    for c = c-1:-1:1
        push!(rows, 1)
        push!(cols, c)
    end
    return D[end, end], reverse(cols), reverse(rows)
end





"""
    dtw_cost(a::AbstractArray, b::AbstractArray, dist::Distances.SemiMetric, r::Int; best_so_far = Inf, cumulative_bound = Zeros(length(a)))

Perform dynamic time warping to measure the distance between two sequences.

Calculate the DTW cost between `a` and `b` with maximum warping radius `r`. You may provide values of `best_so_far` and `cumulative_bound` in order to enable early stopping.

# Keyword arguments:
- `best_so_far`: The best cost value obtained so far (optional)
- `cumulative_bound`: A vector the same length as a and b (optional)
- `s1`: Optional storage vector of length 2r+1, can be used to save allocations.
- `s2`: Optional storage vector of length 2r+1, can be used to save allocations.

Providing the two vectors `s1, s2` does not save very much time, but it makes the function completely allocation free. Can be useful in a threaded context.

See also [`dtw`](@ref) and [`dtwnn`](@ref).
This code was inspired by https://www.cs.ucr.edu/~eamonn/UCRsuite.html
"""
function dtw_cost(
    a::AbstractArray{QT},
    b::AbstractArray,
    dist::F,
    r::Int;
    transportcost = 1,
    best_so_far = typemax(floattype(QT)),
    cumulative_bound = Zeros(lastlength(a)),
    s1 = fill(typemax(floattype(QT)), 2r + 1),
    s2 = fill(typemax(floattype(QT)), 2r + 1),
    kwargs...
) where {QT, F <: Union{Distances.SemiMetric, Function}}

    T = floattype(QT)
    cost      = s1 # just change the name of these variables
    cost_prev = s2

    # Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    m, mb = lastlength(a), lastlength(b)
    mb == m || throw(ArgumentError("a and b must have the same length, got $m and $mb. To compare two series of different lengths, use function dtw"))
    length(cumulative_bound) == m || throw(ArgumentError("cumulative_bound and a must have the same length."))
    length(s1) == 2r+1 || throw(ArgumentError("s1 must be length 2r+1."))
    length(s2) == 2r+1 || throw(ArgumentError("s2 must be length 2r+1."))

    local k

    for i = 0:m-1
        k = max(0, r - i)
        min_cost = typemax(T)

        for j = max(0, i - r):min(m - 1, i + r)
            if j == 0 && i == 0
                cost[k+1] = dist(a[!,1], b[!,1]; kwargs...)
                min_cost = cost[k+1]
                k += 1
                continue
            end
            y = (j - 1 < 0) || (k - 1 < 0)     ? typemax(T) : cost[k]
            x = (i - 1 < 0) || (k + 1 > 2 * r) ? typemax(T) : transportcost*cost_prev[k+2]
            z = (i - 1 < 0) || (j - 1 < 0)     ? typemax(T) : transportcost*cost_prev[k+1]

            cost[k+1] = min(x, y, z) + dist(a[!,i+1], b[!,j+1]; kwargs...)

            # Find minimum cost in row for early stopping
            if cost[k+1] < min_cost
                min_cost = cost[k+1]
            end
            k += 1
        end

        # We can abandon early if the current cumulative distace with lower bound together are larger than best_so_far
        if ((i + r) < (m - 1)) && (min_cost + cumulative_bound[i+r+1] >= best_so_far)
            return min_cost + cumulative_bound[i+r+1]
        end

        cost_prev, cost = cost, cost_prev
    end

    # the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    final_dtw = cost_prev[k]
    return T(final_dtw)
end




@inbounds function soft_dtw_cost_matrix(seq1::AbstractArray, seq2::AbstractArray, dist::SemiMetric = SqEuclidean(); γ = 1,
    transportcost=1, radius=nothing)
    # Build the cost matrix
    m = lastlength(seq2)
    n = lastlength(seq1)

    # Initialize first column and first row
    D = Distances.pairwise(dist, seq2, seq1, dims=2)
    @assert size(D) == (m,n)

    for r=2:m
        D[r,1] += D[r-1,1]
    end
    for c=2:n
        D[1,c] += D[1,c-1]
    end

    # Complete the cost matrix
    if radius === nothing
        for c = 2:n
            for r = 2:m
                D[r, c] += softmin(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1], γ)
            end
        end
    else
        for c = 2:n
            for r = 2:m
                if abs(c-r) > radius
                    D[r, c] += 1/γ
                    # continue
                end    
                D[r, c] += softmin(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1], γ)
            end
        end
    end

    return D
end



"""
    soft_dtw_cost(args...; γ = 1, kwargs...)

Perform Soft DTW. This is a differentiable version of DTW. The "distance" returned by this function is quite far from a true distance and can be negative. A smaller value of `γ` makes the distance closer to the standard DTW distance.


To differentiate w.r.t. the first argument, try
```julia
using ReverseDiff
da = ReverseDiff.gradient(a->soft_dtw_cost(a,b), a)
```

Ref: "Soft-DTW: a Differentiable Loss Function for Time-Series" https://arxiv.org/pdf/1703.01541.pdf

#Arguments:
- `args`: same as for [`dtw`](@ref)
- `γ`: The smoothing factor. A small value means less smoothing and a result closer to [`dtw_cost`](@ref)
- `kwargs`: same as for [`dtw`](@ref)
"""
function soft_dtw_cost(args...; γ = 1, kwargs...)
    D = soft_dtw_cost_matrix(args...; γ = γ, kwargs...)
    D[end,end]
end


@fastmath @inline function softmin(a, b, c, γ)
    γ = -γ
    a,b,c = a/γ, b/γ, c/γ
    maxv = max(a,b,c)
    ae,be,ce = exp(a - maxv), exp(b - maxv), exp(c - maxv)
    γ*(log(ae+be+ce) + maxv)
end

const LVB = LoopVectorization.VectorizationBase
@inline function softmin(a::T, b::T, c::T, γ) where T <: Union{Float64, Float32}
    γ = -γ
    ninvγ = one(T) / γ
    v = LVB.Vec{4,T}(a, b, c, typemax(T))
    v = v * ninvγ
    @fastmath maxv = min(a,b,c) * ninvγ
    ve = exp(v - maxv) * LVB.Vec{4,T}(one(T), one(T), one(T), zero(T))
    γ*(log(LVB.vsum(ve)) + maxv)
end


#### gdtw.jl ####

# A cache to preallocate everything for the GDTW distance
struct GDTWWorkspace{T1, T2, T3}
    τ::T2
    l::T1
    u::T1
    l_prev::T1
    u_prev::T1
    l₀::T1
    u₀::T1
    min_costs::T2
    costs::T3
end

"""
    GDTWWorkspace{T}(M, N)

Creates a cache of numeric type `T` for use in [`gdtw`](@ref).
"""
function GDTWWorkspace(::Type{T}, M, N) where {T}
    GDTWWorkspace{Vector{T}, Matrix{T}, Array{T, 3}}(
        zeros(T, M, N), zeros(T, N), zeros(T, N),
        zeros(T, N), zeros(T, N), zeros(T, N),
        zeros(T, N), zeros(T, M, N), zeros(T, M, M, N)
    )
end

GDTWWorkspace(M, N) = GDTWWorkspace(Float64, M, N)

# refine the bounds, as described in Section 4.2 of DB19
function refine!(l_current, u_current, l_prev, u_prev, l₀, u₀, warp; η)
    @avx for i in eachindex(l_current, u_current, l_prev, u_prev, l₀, u₀, warp)
        δ = η * (u_prev[i] - l_prev[i]) / 2
        l_current[i] = max(warp[i] - δ, l₀[i])
        u_current[i] = min(warp[i] + δ, u₀[i])
    end
    return nothing
end

# inital choices of `l` and `u`, modified from Eq (7) of DB19
function inital_bounds!(l, u, t, smin, smax, symmetric)
    # We need to loosen the bounds to account for floating point error
    # Otherwise these bounds can be too tight and disallow valid moves
    # which can lead to very wrong results.
    smin = .99*smin
    smax = 1.01*smax

    @inbounds for i in eachindex(t, l, u)
        # You must be able to get to `warp[i]` in time `t[i]`, so
        #   `smax >= warp[i] / t[i] >= smin`
        # This gives a lower and upper bound on `warp[i]`, i.e., on `τ[:, i]`.
        # You must be able to get to `1` from `warp[i]` in time `1-t[i]`, so
        #   `smin <= (1-warp[i])/(1-t[i]) <= smax`
        # this gives another lower and upper bound.
        # The resulting bounds:
        lower = max(smin * t[i], 1 - smax * (1 - t[i]))
        upper = min(smax * t[i], 1 - smin * (1 - t[i]))

        if symmetric
            # We need to apply the above bounds to `ψ(s) = 2s - ϕ(s) = 2t[i] - warp[i]`
            # as well Which leads to the following.
            l[i] = max(lower, 2*t[i] - upper)
            u[i] = min(upper, 2*t[i] - lower)
        else
            l[i] = lower
            u[i] = upper
        end
    end
    return nothing
end

# an inplace update to `τ` with the bounds `l` and `u`.
# produces `τ` as described in Section 4.1 of DB19.
function update_τ!(τ, t, M, l, u)
    N = length(t)
    @assert size(τ) == (M, N)
    @inbounds for t = 1:N, j = 1:M
        τ[j, t] = l[t] + ((j - 1) / (M - 1)) * (u[t] - l[t])
    end
    return nothing
end

"""
    gdtw(
        x,
        y,
        ::Type{T}  = Float64;
        symmetric::Bool = true,
        M::Int     = 100,
        N::Int     = 100,
        t          = range(T(0), stop = T(1), length = N),
        cache::GDTWWorkspace = GDTWWorkspace(T, M, length(t)),
        λcum       = T(0.01),
        λinst      = T(0.01),
        η          = T(1 / 8),
        max_iters  = 3,
        metric     = (x, y) -> norm(x - y),
        Rcum       = abs2,
        smin::Real = T(0.001),
        smax::Real = T(5.0),
        Rinst      = symmetric  ?
                        ϕ′ -> ( (smin <= ϕ′ <= smax)
                            && (smin <= 2 - ϕ′ <= smax) ) ? (ϕ′-1)^2 : typemax(T)
                                :
                        ϕ′ -> (smin <= ϕ′ <= smax) ? (ϕ′-1)^2 : typemax(T),
        verbose    = false,
        warp       = zeros(T, length(t)),
        callback   = nothing,
    ) where T -> cost, ϕ, ψ

Computes a general DTW distance following [DB19](https://arxiv.org/abs/1905.12893).

Aims to find `ϕ(s)` to minimize

    ∫ metric(x(ϕ(s)), y(ψ(s))) + λinst*Rinst(ϕ'(s) - 1) + λcum*Rcum(ϕ(s) - s) ds

over the interval `s ∈ [0,1]`, where `ψ(s) = 2s - ϕ(s)` (if `symmetric=true`) or `ψ(s) = s`
(if `symmetric = false`). The integral is discretized in time into `N` points (or according
to the times `t`, if `t` is specified). Additionally, the possible values obtained by `ϕ`
(and hence `ψ`) at each possible time `s` are discretized into `M` points.

If `max_iters > 1`, then after solving the doubly-discretized problem to obtain the optimal `ϕ`,
the problem is solved again by choosing a new discretization of `M` possible values
of `ϕ(s)` in an interval (whose width is governed by the parameter `η`) around the
previous optimal value. This is repeated until the problem has been solved `max_iters`
times in total. Setting `verbose=true` prints the cost at each iteration; a "high enough"
value of `max_iters` can be chosen by inspecting when the cost stabilizes sufficiently.

The parameters are:

* `x`: the continuous time signal to warp (see [`LinearInterpolation`](@ref) for generating such a signal from discrete data)
* `y`: the continuous-time signal to warp to
* `T`: the numeric type to be used in the problem
* `symmetric`: if true, `ψ(s) = 2s - ϕ(s)`, otherwise `ψ(s) = s`.
* `t`: the discretization of time on `[0,1]`; either `t` or `N` should be specified
* `M`: the discretization of the values obtained by the warping path
* `metric`:  a function `metric(u,v) -> ℝ` to compute differences between the signals at a time point (such as a Distances.jl distance)
* `Rcum`: penalty function on the cumulative warp
* `Rinst`: penalty function on the instantaenous warping. Should be infinite outside of `[smin, smax]`.
* `smin`, `smax`: minimum and maximum allowed instantaenous warping. Should have `smin > 0` and `smin < smax`.
* `λcum`, `λinst`: the regularization constants for `Rcum` and `Rinst`, respectively

The following may be pre-allocated and reused between distance computations with the same `M` and `N` (or `length(t)`).

* `cache`: a cache of matrices and vectors, generated by `GDTW.GDTWWorkspace{T}(N,M)`

"""
function gdtw(args...; kwargs...)
    data = prepare_gdtw(args...; kwargs...)
    cost = iterative_gdtw!(data)
    return cost, gdtw_warpings(data)...
end

"""
    prepare_gdtw(x, y; kwargs...)

Creates a NamedTuple of parameters, using the same keyword argments as `dist`.
A preprocessing step before calling `iterative_gdtw!`.
"""
function prepare_gdtw(
    x,
    y,
    ::Type{T}  = Float64;
    symmetric::Bool = true,
    M::Int     = 100,
    N::Int     = 100,
    t::AbstractVector{T} = range(T(0), stop = T(1), length = N),
    cache::GDTWWorkspace = GDTWWorkspace(T, M, length(t)),
    λcum::T    = T(0.01),
    λinst::T   = T(0.01),
    η::T       = T(1 / 8),
    max_iters  = 3,
    metric     = (x, y) -> norm(x - y),
    Rcum       = abs2,
    smin::T    = T(0.001),
    smax::T    = T(5.0),
    Rinst      = symmetric  ?
                    ϕ′ -> ( (smin <= ϕ′ <= smax)
                        && (smin <= 2 - ϕ′ <= smax) ) ? (ϕ′-1)^2 : typemax(T)
                            :
                    ϕ′ -> (smin <= ϕ′ <= smax) ? (ϕ′-1)^2 : typemax(T),
    verbose::Bool = false,
    warp       = zeros(T, length(t)),
    callback   = nothing,
) where T
    N = length(t)

    (M > N / smax) || @warn "`M <= N / smax`; problem may be infeasible" M N smax


    @unpack l₀, l_prev, l,  u₀, u_prev, u, τ = cache
    inital_bounds!(l₀, u₀, t, smin, smax, symmetric)
    l_prev .= l₀
    u_prev .= u₀
    u .= u₀
    l .= l₀
    update_τ!(τ, t, M, l, u)

    function node_weight(j, s)
        s == length(t) && return zero(T)
        Rval = Rcum(τ[j, s] - t[s])
        yval = symmetric ? 2*t[s] - τ[j, s] : t[s]
        (t[s+1] - t[s])*(metric(x(τ[j, s]), y(yval)) + λcum * Rval)
    end

    @inline function edge_weight((j, s), (k, s2))
        s + 1 ≠ s2 && return typemax(T)
        ϕ′ = (τ[k, s+1] - τ[j, s]) / (t[s+1] - t[s])
        (t[s+1] - t[s]) * (λinst * Rinst(ϕ′))
    end

    return (
        iter        = Ref(1),
        N           = N,
        M           = M,
        τ           = τ,
        node_weight = node_weight,
        edge_weight = edge_weight,
        η           = η,
        max_iters   = max_iters,
        t           = t,
        smin        = smin,
        smax        = smax,
        callback    = callback,
        verbose     = verbose,
        metric      = metric,
        cache       = cache,
        warp        = warp,
        symmetric   = symmetric,
    )
end

"""
    iterative_gdtw!(data; max_iters = data.max_iters, verbose = data.verbose) -> cost

Runs the GDTW algorithm with iterative refinement for `max_iters` iterations,
returning the resulting cost. Uses the [`GDTWWorkspace`](@ref) in `data.cache` and
updates the `data.warp` vector. Here, `data` is usually obtained by [`prepare_gdtw`](@ref).

This can be called multiple times on the same `data` with a higher value of `max_iters`
to refine a calculation without starting over, which can be useful for checking convergence.

## Example

```julia
data = prepare_gdtw(x,y; max_iters = 3)
cost3 = iterative_gdtw!(data) # performs 3 iterations, returns the cost
cost10 = iterative_gdtw!(data; max_iters = 10) # performs 7 more iterations, returns the cost
ϕ, ψ = gdtw_warpings(data)

# Equivalently,
cost10, ϕ, ψ = gdtw(x,y; max_iters = 10)
```
"""
function iterative_gdtw!(data; max_iters = data.max_iters, verbose = data.verbose)
    @unpack N, M, τ, η, iter, t, callback = data
    @unpack cache, warp, symmetric = data
    if iter[] > max_iters
        @warn "`iter[] > max_iters`; no iterations performed." iter[] max_iters
        return zero(eltype(warp))
    end
    local cost

    # First iteration is special cased because
    # we can't send the cost and warp to the callback yet,
    # and we can quit early if we don't need to do refinement.
    if iter[] == 1
        if callback !== nothing
            callback((iter=1, t=t, τ=τ))
        end

        cost = single_gdtw!(data)
        verbose && @info "Iteration" iter[] cost
        iter[] += 1
        max_iters == 1 && return cost
    end

    @unpack l_prev, u_prev, l, u, l₀, u₀ = cache

    while iter[] <= max_iters
        l_prev .= l
        u_prev .= u
        refine!(l, u, l_prev, u_prev, l₀, u₀, warp; η=η)
        update_τ!(τ, t, M, l, u)
        cost = single_gdtw!(data)
        if callback !== nothing
            callback((iter = iter, t = t, τ = τ, warp = warp, cost = cost))
        end
        verbose && @info "Iteration" iter[] cost

        iter[] += 1
    end

    return cost
end

"""
    gdtw_warpings(data) -> ϕ, ψ

Computes the interpolations from a `data` `NamedTuple`
with entries for the time points `t`, warping points `warp`,
and a boolean `symmetric`.
"""
function gdtw_warpings(data)
    @unpack t, warp, symmetric = data
    if symmetric
        ψ = LinearInterpolation(2*t - warp, t)
    else
        ψ = LinearInterpolation(t, t)
    end
    ϕ = LinearInterpolation(warp, t)
    return ϕ, ψ
end

## Dynamic programming to compute the distance

function single_gdtw!(data::T) where {T}
    @unpack N, M, node_weight, cache, edge_weight, τ, warp = data
    @unpack min_costs, costs = cache
    calc_costs!(min_costs, costs, N, M, node_weight, edge_weight)
    cost = min_costs[end, end]
    trackback!(warp, costs, τ)
    return cost
end

function calc_costs!(min_costs, costs, N, M, node_weight::F1, edge_weight::F2) where {F1,F2}
    @boundscheck checkbounds(min_costs, 1:M, 1:N)
    @boundscheck checkbounds(costs, 1:M, 1:M, 1:N)

    @inbounds begin
        min_costs .= node_weight.(1:M, permutedims(1:N))
        # t = 2 case
        for j = 1:M
            costs[1, j, 2] = min_costs[1, 1] + edge_weight((1, 1), (j, 2))
            min_costs[j, 2] += costs[1, j, 2]
        end
        for t = 3:N
            for j = 1:M
                mi = typemax(eltype(costs))
                for k = 1:M
                    c = min_costs[k, t-1] + edge_weight((k, t - 1), (j, t))
                    costs[k, j, t] = c
                    mi = ifelse(c < mi, c, mi)
                end
                min_costs[j, t] += mi
            end
        end
    end
    return nothing
end

function trackback!(warp, costs, τ)
    (M, _, N) = size(costs)
    @boundscheck checkbounds(costs, 1:M, 1:M, 1:N)
    c = M
    @inbounds for t = N:-1:3
        warp[t] = τ[c, t]
        c = argmin(@views costs[:, c, t])
    end
    warp[2] = τ[c, 2]
    warp[1] = τ[1, 1]

    return nothing
end

## Interpolations


"""
LinearInterpolation(x::AbstractVector) -> Function

Provides a linear interpolation of `x` on the interval `[0,1]`.
"""
struct LinearInterpolation{Tx,Tt} <: Function
    x::Tx
    t::Tt
    function LinearInterpolation(x::Tx, ts::Ts) where {Tx,Ts}
        issorted(ts) || throw(ArgumentError("Time parameter `ts` must be sorted in increasing order."))
        T = eltype(Tx)
        t = (ts .- T(first(ts))) ./ T( last(ts) - first(ts))
        Tt = typeof(t)
        new{Tx,Tt}(x, t)
    end
end

LinearInterpolation(x) = LinearInterpolation(x, axes(x, ndims(x)))

function (xt::LinearInterpolation)(s)
    x = xt.x
    0 <= s <= 1 || return zero(x[!, 1])
    t = xt.t
    i = searchsortedlast(t, s)
    (i == 0) && return x[!, 1]
    (i == lastlength(x)) && return x[!, lastlength(x)]
    (s == t[i]) && return x[!, i]
    weight = (s - t[i]) / (t[i+1] - t[i])
    omw = 1 - weight
    x[!, i] .* omw .+ x[!, i+1] .* weight
end

#### dtwnn.jl ####

struct DTWWorkspace{T,N,AT<:AbstractArray,D}
    q::AT
    dist::D
    r::Int
    l::Vector{T}
    u::Vector{T}
    l_buff::Vector{T}
    u_buff::Vector{T}
    cb::Vector{T}
    c1::Vector{T}
    c2::Vector{T}
    function DTWWorkspace(q::AbstractArray{QT}, dist, r::Int, ::Type{N}=Nothing) where {QT, N}
        T      = floattype(QT)
        m      = lastlength(q)
        n      = 2r + 1
        l      = zeros(T, m)
        u      = zeros(T, m)
        l_buff = zeros(T, m)
        u_buff = zeros(T, m)
        cb     = zeros(T, m)
        c1     = zeros(T, n)
        c2     = zeros(T, n)
        new{T, N, typeof(q), typeof(dist)}(q, dist, r, l, u, l_buff, u_buff, cb, c1, c2)
    end
end


struct DTWSearchResult{QT,C,D} <: AbstractSearchResult{QT}
    q::QT
    cost::C
    loc::Int
    prunestats
    dists::D
end

SlidingDistancesBase.value(r::DTWSearchResult) = r.cost
SlidingDistancesBase.location(r::DTWSearchResult) = r.loc
SlidingDistancesBase.payload(r::DTWSearchResult) = r.dists
SlidingDistancesBase.target(r::DTWSearchResult) = r.q

Base.findmin(results::Vector{<:DTWSearchResult}) = (i=argmin(results); (results[i].cost,i))
Base.findmax(results::Vector{<:DTWSearchResult}) = _findres(results, >)
Base.minimum(results::Vector{<:DTWSearchResult}) = findmin(results)[1]
Base.maximum(results::Vector{<:DTWSearchResult}) = findmax(results)[1]

function _findres(results::Vector{<:DTWSearchResult}, comp)
    mapreduce((a,b)->comp(a[1], b[1]) ? a : b, enumerate(results)) do (i,r)
        minimum(r.dists), i
    end
end

function lower_upper_envs!(w::DTWWorkspace{T}, q, bsf, query = false) where {T}
    du, dl = Deque{Int}(), Deque{Int}()
    push!(du, 0)
    push!(dl, 0)
    r = w.r
    if query
        u,l = w.u, w.l
    else
        u,l = w.u_buff, w.l_buff
    end
    m = lastlength(w.q)
    for i = 1:m-1
        if i > r
            u[i-r] = q[first(du)+1]
            l[i-r] = q[first(dl)+1]
        end
        if q[i+1] > q[i]
            pop!(du)
            while (!isempty(du) && q[i+1] > q[last(du)+1])
                pop!(du)
            end
        else
            pop!(dl)
            while (!isempty(dl) && q[i+1] < q[last(dl)+1])
                pop!(dl)
            end
        end
        push!(du, i)
        push!(dl, i)
        if i == 2r + 1 + first(du)
            popfirst!(du)
        elseif (i == 2r + 1 + first(dl))
            popfirst!(dl)
        end
    end
    for i = m:m+r
        u[i-r] = q[first(du)+1]
        l[i-r] = q[first(dl)+1]
        if i - first(du) >= 2r + 1
            popfirst!(du)
        end
        if i - first(dl) >= 2r + 1
            popfirst!(dl)
        end
    end
end

function lb_endpoints(dist, q, buffer, best_so_far; kwargs...)
    m = lastlength(q)

    x1 = buffer[!,1]
    y1 = buffer[!,m]
    lb = dist(q[!,1], x1; kwargs...) + dist(q[!,m], y1; kwargs...)
    lb >= best_so_far && return lb

    x2 = buffer[!,2]
    d = min(dist(x2, q[!,1]; kwargs...), dist(x1,q[!,2]; kwargs...), dist(x2,q[!,2]); kwargs...)
    lb += d
    lb >= best_so_far && return lb

    y2 = buffer[!,m-1]
    d = min(dist(y2, q[!,m]; kwargs...), dist(y1,q[!,m-1]; kwargs...), dist(y2,q[!,m-1]); kwargs...)
    lb += d
    lb >= best_so_far && return lb

    return lb
    # TODO: can add more comparisons here
end

function lb_env!(w::DTWWorkspace{T}, buffer, best_so_far; kwargs...) where T
    lb = zero(T)
    q, dist, u, l = w.q, w.dist, w.u, w.l
    for i in 1:lastlength(q)
        x = buffer[!,i] # This function only supports data with natural ordering
        d = zero(T)
        if x > u[i]
            d = dist(x, u[i]; kwargs...)
        elseif x < l[i]
            d = dist(x, l[i]; kwargs...)
        end
        lb += d
        w.cb[i] = d
        lb > best_so_far && return lb
    end
    return lb
end

function rev_cumsum!(cb)
    @inbounds for k = length(cb)-1:-1:1
        cb[k] = cb[k+1] + cb[k]
    end
end

"""
    search_result = dtwnn(q, y, dist, rad, [normalizer::Type{Nothing}]; kwargs...)

Compute the nearest neighbor to `q` in `y`. An optinal normalizer type can be supplied, see, `ZNormalizer, DiagonalZNormalizer, NormNormalizer`.

# Arguments:
- `q`: query (the short time series)
- `y`: data ( the long time series)
- `dist`: distance
- `rad`: radius
- `showprogress`: Defaults to true
- `prune_endpoints = true`: use endpoint heuristic
- `prune_envelope  = true`: use envelope heuristic
- `bsf_multiplier  = 1`: If > 1, require lower bound to exceed `bsf_multiplier*best_so_far`.
- `saveall = false`: compute a dense result (takes longer, no early stopping methods used). If false, then a vector of lower bounds on the distance is stored in `search_result.dists`, if true, all distances are computed and stored.
- `avoid`: If an integer index (or set of indices) is provided, this index will be avoided in the search. This is useful in case `q` is a part of `y`.
"""
function dtwnn(q::AbstractArray{QT}, y, dist, rad, ::Type{N}=Nothing; kwargs...) where {QT,N}
    q, y = setup_normalizer(N, q, y)
    w = DTWWorkspace(q, dist, rad, N)::DTWWorkspace{floattype(QT),N,typeof(q),typeof(dist)}
    dtwnn(w, y; kwargs...)
end

function dtwnn(w::DTWWorkspace{T,normalizer}, y::AbstractArray;
    prune_endpoints = true,
    prune_envelope  = true,
    saveall         = false,
    bsf_multiplier  = 1,
    transportcost   = 1,
    showprogress    = true,
    avoid           = nothing,
    kwargs...) where {T,normalizer}


    normalizer !== Nothing && !isa(y, AbstractNormalizer) && @warn("Normalizer in use but `y` is not wrapped in a normalizer object. This will result in highly suboptimal performance", maxlog=10)
    bsf_multiplier >= 1 || throw(DomainError("It does not make sense to have the bsf_multiplier < 1"))
    best_so_far = typemax(T)
    best_loc    = 1
    q           = w.q
    m           = lastlength(q)
    my          = actuallastlength(y)
    my >= m || throw(ArgumentError("q must be shorter than y, swap inputs."))
    onedim      = ndims(q) == 1 && eltype(q) <: Real
    onedim && prune_envelope && lower_upper_envs!(w, q, best_so_far, true) # Result stored in w

    # Counters to keep track of how many times lb helps
    prune_end   = 0
    prune_env   = 0
    dists       = fill(typemax(T), my-m+1)

    prog = Progress((my-m)÷max(my÷100,1), dt=1, desc="DTW NN")
    # @inbounds @showprogress 1.5 "DTW NN" for it = 1:my-m
    @inbounds for it = 1:my-m+1
        showprogress && it % max(my÷100,1) == 0 && next!(prog)
        advance!(y)
        avoid !== nothing && it ∈ avoid && continue
        bsf = bsf_multiplier*best_so_far
        ym = getwindow(y, m, it)
        if prune_endpoints && !saveall
            lb_end = lb_endpoints(w.dist, w.q, ym, bsf; kwargs...)
            if lb_end > bsf
                prune_end += 1
                continue
            end
        end
        if onedim && prune_envelope && !saveall # This bound only works when there is a natural ordering
            # lower_upper_envs!(w, ym, bsf) # This step is only required for reverse bound
            lb_env = lb_env!(w, ym, bsf; kwargs...) # updates w.cb
            rev_cumsum!(w.cb)
            if lb_env > bsf
                prune_env += 1
                continue
            end
        end
        # If we get here, we must normalize the entire y
        buffern = normalize(normalizer, ym) # This only normalizes what's not already normalized

        newdist = dtw_cost(q, buffern, w.dist, w.r;
            cumulative_bound = w.cb,
            best_so_far      = saveall ? typemax(T) : bsf,
            s1               = w.c1,
            s2               = w.c2,
            transportcost    = transportcost,
            kwargs...
        )
        dists[it] = newdist
        if newdist < best_so_far
            best_so_far = newdist
            best_loc = it
        end
    end
    prunestats = (prune_end=prune_end, prune_env=prune_env)
    DTWSearchResult(q, best_so_far, best_loc, prunestats, dists)
end

struct Neighbor{T}
    i::Int
    d::T
end

Base.isless(n1::Neighbor,n2::Neighbor) = isless(n1.d, n2.d)
Base.isless(n1, n2::Neighbor) = isless(n1, n2.d)
Base.isless(n1::Neighbor,n2) = isless(n1.d, n2)

"""
    dists, inds = sparse_distmat(y::Vector{<:AbstractVector{T}}, k, dist, radius; kwargs...) where T

Compute the `k` nearest neighbors between signals in `y`, corresponding to the `k` smallest entries in each row of the pairwise distance matrix. The return values are vectors of length-k vectors with the calculated distances and neighbor indices.

#Arguments:
- `y`: Vector of vectors containing the signals
- `k`: number of neighbors
- `dist`: the inner metric, e.g., `SqEuclidean()`
- `kwargs`: these are sent to `dtw_cost`.
- `showprogress = true`
"""
function sparse_distmat(
    y::AbstractVector{<:AbstractArray{S}},
    k,
    dist,
    rad;
    showprogress::Bool = true,
    kwargs...,
) where {S}
    T     = floattype(S)
    N     = length(y)
    INDS  = [zeros(Int, k) for _ = 1:N]
    DISTS = [zeros(T, k) for _ = 1:N]
    showprogress && (p = Progress(N^2, 1, "sparse_distmat"))
    for i = 1:N
        bsf   = typemax(T)
        dists = BinaryMaxHeap{Neighbor{T}}()
        for j = 1:N
            j == i && (showprogress && next!(p); continue)
            d = lb_endpoints(dist, y[i], y[j], bsf; kwargs...)
            if d < bsf
                d = dtw_cost(y[i], y[j], dist, rad; best_so_far = bsf, kwargs...)
            end
            push!(dists, Neighbor(j, T(d)))
            if length(dists) > k
                bsf = pop!(dists)
            end
            showprogress && next!(p)
        end

        for j = k:-1:1
            n = pop!(dists)
            INDS[i][j] = n.i
            DISTS[i][j] = n.d
        end
    end
    DISTS, INDS
end


#### dba.jl ####

using LinearAlgebra

"""
    DBAResult(cost,converged,iterations,cost_trace)

Holds results of a DTW Barycenter Averaging (DBA) fit.
"""
mutable struct DBAResult
    cost::Float64
    converged::Bool
    iterations::Int
    cost_trace::Vector{Float64}
end

"""
    avgseq, results = dba(sequences, dist::DTWDistance; kwargs...)

Perfoms DTW Barycenter Averaging (DBA) given a collection of `sequences`
and the current estimate of the average sequence.

Example usage:

    x = [1., 2., 2., 3., 3., 4.]
    y = [1., 3., 4.]
    z = [1., 2., 2., 4.]
    avg,result = dba([x,y,z], DTW(3))
"""
function dba(
    sequences::AbstractVector,
    dtwdist::DTWDistance;
    init_center = rand(sequences),
    iterations::Int = 1000,
    rtol::Float64 = 1e-5,
    store_trace::Bool = false,
    show_progress::Bool = true,
    i2min::AbstractVector = [],
    i2max::AbstractVector = [],
)

    # initialize dbavg as a random sample from the dataset
    nseq = length(sequences)
    dbavg = deepcopy(init_center)

    # storage for each iteration
    newavg = zeros(size(dbavg))
    counts = zeros(Int, lastlength(dbavg))

    # variables storing optimization progress
    converged = false
    iter = 0
    cost, newcost = Inf, Inf
    cost_trace = Float64[]

    # display optimization progress
    if show_progress
        p = ProgressMeter.ProgressThresh(rtol)
    end

    ## main loop ##
    while !converged && iter < iterations

        # do an iteration of dba
        newcost = dba_iteration!(
            newavg,
            dbavg,
            counts,
            sequences,
            dtwdist;
            i2min = i2min,
            i2max = i2max,
        )
        iter += 1

        # store history of cost while optimizing (optional)
        store_trace && push!(cost_trace, newcost)

        # check convergence
        Δ = (cost - newcost) / newcost
        if Δ < rtol
            converged = true
        else
            # update estimate
            cost = newcost
            dbavg = deepcopy(newavg)
        end

        # update progress bar
        if show_progress
            ProgressMeter.update!(
                p,
                Δ;
                showvalues = [
                    (:iteration, iter),
                    (Symbol("max iteration"), iterations),
                    (:cost, cost),
                ],
            )
        end
    end

    return newavg, DBAResult(newcost, converged, iter, cost_trace)
end


"""

Performs one iteration of DTW Barycenter Averaging (DBA) given a collection of
`sequences` and the current estimate of the average sequence, `dbavg`. Returns
an updated estimate, and the cost/loss of the previous estimate
"""
function dba_iteration!(
    newavg::T,
    oldavg::T,
    counts::Array{Int,1},
    sequences::AbstractVector{T},
    d::DTWDistance;
    i2min::AbstractVector = [],
    i2max::AbstractVector = [],
) where {T}

    # sum of dtw dist of all sequences to center
    total_cost = 0.0

    # store stats for barycenter averages
    counts .= 0
    newavg .= 0

    # main ploop
    for seq in sequences
        # time warp signal versus average
        # if one of the two is empty, use unconstrained window. If both are nonempty, but not the same length, distpath will throw error
        if isempty(i2min) && isempty(i2max)
            cost, i1, i2 = distpath(d, oldavg, seq)
        else
            cost, i1, i2 = distpath(d, oldavg, seq, i2min, i2max)
        end
        total_cost += cost

        # store stats for barycentric average
        for j = 1:length(i2)
            counts[i1[j]] += 1
            newavg[!,i1[j]] += seq[!,i2[j]]
        end
    end

    # compute average and return total cost
    for i in eachindex(counts)
        newavg[!,i] = newavg[!,i] / counts[i]
    end

    return total_cost
end


# weirdly enought, this works for the dtw_dba_miniexample,  but does not work for dtw_dbaclust
#@generated function _sequentize{T,N}(s::AbstractArray{T,N})
#    :( Sequence[ Sequence(@ncall($N, view, s, n-> n==$N ? i : Colon())) for i = 1:size(s,2) ] )
#end


#### dbaclust.jl ####

"""
    DBAclustResult(centers,clustids,result)

Holds results of a DBAclust run.

"""
mutable struct DBAclustResult{T}
    centers::T
    clustids::Array{Int}
    converged::Bool
    iterations::Int
    dbaresult::DBAResult
end

"""

    optional_threaded(cond, ex)

Execute ex with multi-threading if cond is true,
otherwise execute with one core
"""
macro optional_threaded(cond, ex)
    quote
        if $(esc(cond))
            $(esc(:(Threads.@threads $ex)))
        else
            $(esc(ex))
        end
    end
end


"""
    dbaclust(
        sequences,
        nclust::Int,
        dtwdist::DTWDistance;
        n_init::Int           = 1,
        iterations::Int       = 100,
        inner_iterations::Int = 10,
        rtol::Float64         = 1e-4,
        rtol_inner::Float64   = rtol,
        threaded::Bool        = false,
        show_progress::Bool   = true,
        store_trace::Bool     = true,
        i2min::AbstractVector = [],
        i2max::AbstractVector = [],
    )


# Arguments:
- `nclust`: Number of clsuters
- `n_init`: Number of initialization tries
- `inner_iterations`: Number of iterations in the inner alg.
- `i2min`: Bounds on the warping path
- `threaded`: Use multi-threading 
"""
function dbaclust(
    sequences,
    nclust::Int,
    dtwdist::DTWDistance;
    n_init::Int           = 1,
    iterations::Int       = 100,
    inner_iterations::Int = 10,
    rtol::Float64         = 1e-4,
    rtol_inner::Float64   = rtol,
    threaded::Bool        = false,
    show_progress::Bool   = true,
    store_trace::Bool     = true,
    i2min::AbstractVector = [],
    i2max::AbstractVector = [],
)

    n_init < 1 && throw(ArgumentError("n_init must be greater than zero"))

    T = typeof(sequences)

    results = Array{DBAclustResult{T}}(undef, n_init)
    show_progress && (p = Progress(n_init))

    @optional_threaded threaded for i = 1:n_init
        results[i] = dbaclust_single(
            sequences,
            nclust,
            dtwdist;
            iterations = iterations,
            inner_iterations = inner_iterations,
            rtol = rtol,
            rtol_inner = rtol_inner,
            threaded = threaded,
            show_progress = false,
            store_trace = store_trace,
            i2min = i2min,
            i2max = i2max,
        )
        show_progress && next!(p)
    end
    best = results[1]

    for i = 2:n_init
        if results[i].dbaresult.cost < best.dbaresult.cost
            best = results[i]
        end
    end

    return best
end


"""
    avgseq, results = dbaclust_single(sequences, dist; kwargs...)

Perfoms a single DTW Barycenter Averaging (DBA) given a collection of `sequences`
and the current estimate of the average sequence.

Example usage:

    x = [1,2,2,3,3,4]
    y = [1,3,4]
    z = [1,2,2,4]
    avg,result = dba([x,y,z], DTW(3))
"""
function dbaclust_single(
    sequences::AbstractVector,
    nclust::Int,
    dtwdist::DTWDistance;
    threaded::Bool        = false,
    init_centers::AbstractVector = dbaclust_initial_centers(
        sequences,
        nclust,
        dtwdist;
        threaded
    ),
    iterations::Int       = 100,
    inner_iterations::Int = 10,
    rtol::Float64         = 1e-4,
    rtol_inner::Float64   = rtol,
    show_progress::Bool   = true,
    store_trace::Bool     = true,
    i2min::AbstractVector = [],
    i2max::AbstractVector = [],
)

    T = floattype(eltype(sequences))
    # rename for convienence
    avgs   = init_centers
    N = length(avgs[1])

    # check initial centers have the same length
    if !all(length(a) == N for a in avgs)
        throw(ArgumentError("all initial centers should be the same length"))
    end

    # dimensions
    nseq      = length(sequences)
    maxseqlen = maximum([length(s) for s in sequences])


    # TODO switch to ntuples?
    counts    = [zeros(Int, N) for _ = 1:nclust]
    sums      = [Array{T}(undef,N) for _ = 1:nclust]

    # cluster assignments for each sequence
    clus_asgn = Array{Int}(undef,nseq)
    c         = 0

    # arrays storing path through dtw cost matrix
    i1, i2    = Int[], Int[]

    # variables storing optimization progress
    converged       = false
    iter            = 0
    inner_iter      = 0
    converged_inner = false
    last_cost       = Inf
    total_cost      = 0.0
    cost_trace      = Float64[]
    costs           = Array{Float64}(undef,nseq)

    # main loop ##
    if show_progress
        prog = ProgressMeter.ProgressThresh(rtol, 2)
        ProgressMeter.update!(
            prog,
            Inf;
            showvalues = [
                (:iteration, iter),
                (Symbol("max iteration"), iterations),
                (:cost, total_cost),
            ],
        )
    end#showprogress

    while !converged && iter < iterations

        # first, update cluster assignments based on nearest
        # centroid (measured by dtw distance). Keep track of
        # total cost (sum of all distances to centers).
        total_cost = 0.0
        for s = 1:nseq
            # process sequence s
            seq = sequences[s]

            # find cluster assignment for s
            costs[s] = Inf
            if threaded
                # using multi-core
                cluster_dists_      = Array{Float64}(undef, nclust)
                cluster_i1s_        = Array{Vector{Int64}}(undef, nclust)
                cluster_i2s_        = Array{Vector{Int64}}(undef, nclust)

                Threads.@threads for c_ = 1:nclust
                    # if one of the two is empty, use unconstrained window. If both are nonempty, but not the same lenght, distpath will throw error
                    if isempty(i2min) && isempty(i2max)
                        cost, i1_, i2_ = distpath(dtwdist, avgs[c_], seq)
                    else
                        cost, i1_, i2_ = distpath(dtwdist, avgs[c_], seq, i2min, i2max)
                    end
                    cluster_dists_[c_]  = cost
                    cluster_i1s_[c_]    = i1_
                    cluster_i2s_[c_]    = i2_
                end
                cost, c = findmin(cluster_dists_)
                i1      = cluster_i1s_[c]
                i2      = cluster_i2s_[c]
                costs[s] = cost
            else
                # using single-core
                for c_ = 1:nclust
                    # if one of the two is empty, use unconstrained window. If both are nonempty, but not the same lenght, distpath will throw error
                    if isempty(i2min) && isempty(i2max)
                        cost, i1_, i2_ = distpath(dtwdist, avgs[c_], seq)
                    else
                        cost, i1_, i2_ = distpath(dtwdist, avgs[c_], seq, i2min, i2max)
                    end
                    if cost < costs[s]
                        # store cluster, and match indices
                        c           = c_
                        i1          = i1_
                        i2          = i2_
                        costs[s]    = cost
                    end
                end
            end

            # s was assigned to cluster c
            clus_asgn[s] = c
            cnt          = counts[c]
            sm           = sums[c]
            avg          = avgs[c]

            # update stats for barycentric average for
            # the assigned cluster
            for t = 1:length(i2)
                cnt[i1[t]] += 1
                sm[i1[t]]  += seq[i2[t]]
            end
        end

        # if any centers are unused, and reassign them to the sequences
        # with the highest cost
        unused = setdiff(1:nclust, unique(clus_asgn))
        if !isempty(unused)
            # reinitialize centers
            @optional_threaded threaded for c in unused
                avgs[c] = deepcopy(sequences[argmax(costs)])
                @optional_threaded threaded for s = 1:nseq
                    seq = sequences[s]
                    if isempty(i2min) && isempty(i2max)
                        cost, = distpath(dtwdist, avgs[c], seq)
                    else
                        cost, = distpath(dtwdist, avgs[c], seq, i2min, i2max)
                    end
                    if costs[s] > cost
                        costs[s] = cost
                    end
                end
            end
            # we need to reassign clusters, start iteration over
            continue
        end

        # store history of cost while optimizing (optional)
        total_cost = sum(costs)
        store_trace && push!(cost_trace, total_cost)

        # check convergence
        Δ = (last_cost - total_cost) / total_cost
        if Δ < rtol
            converged = true
        else
            last_cost = total_cost
        end

        # update barycenter estimates
        for (a, s, c) in zip(avgs, sums, counts)
            for t = 1:N
                c[t] == 0 && continue
                a[t] = s[t] / c[t]
            end
            # zero out sums and counts for next iteration
            s .= 0
            c .= 0
        end

        # add additional inner dba iterations

        for i = 1:nclust
            seqs            = view(sequences, clus_asgn .== i)
            inner_iter      = 0
            converged_inner = false
            oldcost         = 1.0e100
            while !converged_inner && inner_iter < inner_iterations
                newcost = dba_iteration!(
                    sums[i],
                    avgs[i],
                    counts[i],
                    seqs,
                    dtwdist;
                    i2min = i2min,
                    i2max = i2max,
                )
                copy!(avgs[i], sums[i])
                inner_iter += 1
                δ = (oldcost - newcost) / oldcost
                if δ < rtol_inner
                    converged_inner = true
                else
                    oldcost = newcost
                end
            end
        end

        # update progress bar
        iter += 1
        show_progress && ProgressMeter.update!(
            prog,
            Δ;
            showvalues = [
                (:iteration, iter),
                (Symbol("max iteration"), iterations),
                (:cost, total_cost),
            ],
        )
    end

    return DBAclustResult(
        avgs,
        clus_asgn,
        converged,
        iter,
        DBAResult(total_cost, converged_inner, inner_iter, cost_trace),
    )
end


"""
   dbaclust_initial_centers(sequences, nclust, dtwdist::DTWDistance; threaded::Bool = false)

Uses kmeans++ (but with dtw distance) to initialize the centers
for dba clustering.
"""
function dbaclust_initial_centers(
    sequences::AbstractVector,
    nclust::Int,
    dtwdist::DTWDistance;
    threaded::Bool        = false
)
    # number of sequences in dataset
    nseq          = length(sequences)
    # distance of each datapoint to each center
    dists         = zeros(nclust, nseq)
    # distances to closest center
    min_dists     = zeros(1, nseq)
    # choose a center uniformly at random
    center_ids    = zeros(Int, nclust)
    center_ids[1] = rand(1:nseq)

    # assign the rest of the centers
    for c = 1:(nclust-1)

        # first, compute distances for the previous center
        cent = sequences[center_ids[c]]
        @optional_threaded threaded for i = 1:nseq
            # this distance will be zero
            i == center_ids[c] && continue
            # else, compute dtw distance
            seq = sequences[i]
            dists[c, i], = distpath(dtwdist, seq, cent)
        end

        # for each sequence, find distance to closest center
        minimum!(min_dists, dists[1:c, :])

        min_dists .= abs2.(min_dists)

        # sample the next center
        center_ids[c+1] = sample(1:nseq, Weights(view(min_dists, :)))
    end

    # return list of cluster centers
    return [deepcopy(sequences[c]) for c in center_ids]
end


#### windowed_matrix.jl ####

#
# Define a "windowed matrix" type, where a maximum and minimum row is specified for
# each column.  Trying to read the matrix outside that range yields a default value
# (this defaults to Inf, though it can be specified when constructing the object).
# Trying to write outside that range throws a BoundsError().
#
# This is provided for efficient implementation of restricted Dynamic Time Warping
# algorithms.
#
# Joe Fowler
# NIST Boulder Laboratories
# December 2014
#

struct WindowedMatrix{T<:Real} <: AbstractArray{T,2}
    nrow::Int
    ncol::Int
    ncells::Int

    cost::Vector{T}
    rowmin::Vector{Int}
    rowmax::Vector{Int}
    rowspercol::Vector{Int}
    idxcol::Vector{Int}  # Index (in cost) of 1st element in each column
    defaultval::T

    function WindowedMatrix{T}(
        rmin::Vector{Int},
        rmax::Vector{Int},
        default::T,
    ) where {T<:Real}
        rowmin = copy(rmin)
        rowmax = copy(rmax)
        rowspercol = rowmax .- rowmin .+ 1
        ncells = sum(rowspercol)
        nrow = maximum(rowmax)
        ncol = length(rowmax)
        cost = zeros(T, ncells)
        idxcol = 1 .+ vcat([0], cumsum(rowspercol[1:end-1])) # index of 1st element per column

        new(nrow, ncol, ncells, cost, rowmin, rowmax, rowspercol, idxcol, default)
    end
end


# If the default matrix value is given, then the matrix takes on its type
function WindowedMatrix(rmin::Vector{Int}, rmax::Vector{Int}, default::T) where {T<:Real}
    WindowedMatrix{T}(rmin, rmax, default)
end

# If no default matrix value is given, it will be Inf, and matrix will hold
# Float64 values.
WindowedMatrix(rmin::Vector{Int}, rmax::Vector{Int}) =
    WindowedMatrix{Float64}(rmin, rmax, Inf)

Base.size(W::WindowedMatrix) = W.nrow, W.ncol
Base.size(W::WindowedMatrix, i) = size(W)[i]

@inline function Base.getindex(W::WindowedMatrix, r::Integer, c::Integer)
    if c < 1 || c > W.ncol || r < W.rowmin[c] || r > W.rowmax[c]
        return W.defaultval
    end

    offset = r - W.rowmin[c]
    W.cost[W.idxcol[c]+offset]
end


@inline function Base.setindex!(W::WindowedMatrix, val, r::Integer, c::Integer)
    Base.@boundscheck if c < 1 || c > W.ncol || r < W.rowmin[c] || r > W.rowmax[c]
        throw(BoundsError())
    end

    offset = r - W.rowmin[c]
    W.cost[W.idxcol[c]+offset] = val
    return W
end



#### fastdtw.jl ####


"""
    cost,i1,i2 = fastdtw(seq1,seq2,dist,radius)

Perform dynamic-time warping using the FastDTW algorithm to measure the distance between two sequences. Note that regular DTW often performs better than FastDTW https://arxiv.org/abs/2003.11246

Computes FastDTW approximation to the DTW, described in Salvador & Chan,
Intelligent Data Analysis (2007).

See also [`dtw`](@ref), [`dtw_cost`](@ref), [`dtwnn`](@ref).
"""
function fastdtw(
        seq1::AbstractArray,
        seq2::AbstractArray,
        dist::Union{SemiMetric, Function},
        radius::Int,
    )

    MinSize = max(radius + 2, 10)
    N1 = lastlength(seq1)
    N2 = lastlength(seq2)
    if N1 <= MinSize || N2 <= MinSize
        return (dtw(seq1, seq2, dist))
    end

    # Call recursively on a pair of sequences half this length
    compressed1 = compress2(seq1)
    compressed2 = compress2(seq2)
    _cost, lowrescol, lowresrow = fastdtw(compressed1, compressed2, dist, radius)

    # Now resample that path to the finer resolution, find the correct
    # window around it, and get the DTW given that window.
    hirescol, hiresrow = expandpath(lowrescol, lowresrow, N1, N2)
    idx2min, idx2max = computewindow(hirescol, hiresrow, radius)
    cost1, newcol, newrow = dtw(seq1, seq2, dist, idx2min, idx2max)
end


# Given a path through low-res space, generate an approximate path
# through high-res space. It should have dimension Ncol x Nrow

@inbounds function expandpath(lowrescol, lowresrow, Ncol, Nrow)
    @assert div(Ncol+1,2) == lowrescol[end]
    @assert div(Nrow+1,2) == lowresrow[end]
    Np = length(lowrescol)
    @assert Np == length(lowresrow)

    hirescol = zeros(eltype(lowrescol), 2*Np)
    hiresrow = zeros(eltype(lowresrow), 2*Np)
    hirescol[1] = hiresrow[1] = c = r = 1
    for i=1:Np-1
        # Select plan according to the next move in lowres path.
        if lowrescol[i+1] == lowrescol[i]  # Next move is up
            r += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            r += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r

        elseif lowresrow[i+1] == lowresrow[i] # Next move is sideways
            c += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            c += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r

        else  # Next move is diagonal.
            c += 1; r += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            c += 1; r += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r
        end
    end
    hirescol[end] = Ncol
    hiresrow[end] = Nrow
    # When expanding to an odd numbered size, it's possible to repeat
    # the last step.  Fix that:
    if hirescol[end]==hirescol[end-1] && hiresrow[end]==hiresrow[end-1]
        hirescol = hirescol[1:end-1]
        hiresrow = hiresrow[1:end-1]
    end
    hirescol, hiresrow
end

# yshort = compress2(y)
#   Returns a shortened time series that is half the length of the input sequence.
#   The length of the compressed sequence is always even.
function compress2(seq::AbstractArray)
    # Navg = div(length(seq), 2)
    n = lastlength(seq)
    evenseq = 0.5*(seq[!,1:2:n-1]+seq[!,2:2:n])
    if lastlength(seq)%2 == 1
        return cat(evenseq, seq[!, n], dims=ndims(seq))
    end
    evenseq
end



# Given the lists of (col,row) indices for the optimal path, compute a "window"
# around that path of the given radius.
# Returns (rowmin, rowmax), each a vector of length pathcols[end], representing
# for each column, the minimum and maximum row numbers used in that column.

@inbounds function computewindow(pathcols, pathrows, radius)
    Np = length(pathcols)
    @assert Np == length(pathrows)
    Ncol = pathcols[end]
    Nrow = pathrows[end]

    # Find the min/max row at each column in the path.
    pathmin = zeros(Int, Ncol)
    pathmax = zeros(Int, Ncol)
    for i=1:Np
        c,r = pathcols[i], pathrows[i]
        pathmax[c] = r
        if pathmin[c] == 0
            pathmin[c] = r
        end
    end

    # The window in each column for "radius" r starts at the pathmin
    # of the rth-previous column and ends at the pathmax of the
    # rth-next column, plus (in each case) the radius.
    if radius < Ncol-1 && radius < Nrow-1
        rowmin = vcat(fill(1,radius), pathmin[1:end-radius] .- radius)
        rowmax = vcat(pathmax[radius+1:end] .+ radius, fill(Nrow,radius))

        # Window values must be in the range [1:Nrow].
        for c=1:Ncol
            if rowmin[c]<1; rowmin[c]=1; end
            if rowmax[c]>Nrow; rowmax[c]=Nrow; end
        end
    else
        rowmin = fill(1,Ncol)
        rowmax = fill(Nrow,Ncol)
    end
    rowmin, rowmax
end

#### plots.jl ####

export dtwplot


"""
    dtwplot(seq1, seq2, [dist=SqEuclidean()]; transportcost=1, diagonal=false)
    dtwplot(seq1, seq2, D, i1, i2; transportcost=1, diagonal=false)

Given two sequences, perform dynamic time warping and plot
the results. If alignment has already been computed, pass
the indices `i1` and `i2` to make the plot.

`diagonal = true` plots a diagonal marker as visual aid.
"""
dtwplot

function handleargs(seq1, seq2, dist::SemiMetric = SqEuclidean(); kwargs...)
    D = dtw_cost_matrix(seq1, seq2, dist; kwargs...)
    cost, i1, i2 = trackback(D)
    seq1, seq2, D, i1, i2
end

function handleargs(seq1, seq2, dist, i2min, i2max; kwargs...)
    D = dtw_cost_matrix(seq1, seq2, dist, i2min, i2max; kwargs...)
    cost, i1, i2 = trackback(D)
    seq1, seq2, D, i1, i2
end

handleargs(h; kwargs...) = handleargs(h.args...; kwargs...)

@userplot DTWPlot

@recipe function f(h::DTWPlot; transportcost=1, diagonal=false, postprocess=nothing)
    seq1, seq2, D, i1, i2 = handleargs(h; transportcost=transportcost, postprocess=postprocess)

    n1, n2 = lastlength(seq1), lastlength(seq2)

    all = ndims(seq1) ∈ (1,2)
    # set up the subplots
    legend --> false
    link := :both
    grid --> false
    if all
        layout --> @layout [
            a b{0.8w,0.8h}
            _ c
        ]
    else
        layout --> 1
    end

    left_margin --> 0mm
    bottom_margin --> 0mm
    top_margin --> 0mm
    right_margin --> 0mm

    # heatmap
    @series begin
        clims --> (0, 3 * D[end, end])
        seriestype := :heatmap
        formatter --> (z) -> ""
        subplot := (all ? 2 : 1)
        D
    end

    # the rest of the plots are paths

    # main plot
    s1 = @series begin
        seriestype := :path
        linecolor --> :auto
        linewidth --> 3
        subplot := (all ? 2 : 1)
        formatter --> (z) -> ""
        i1, i2
    end

    if diagonal
        m2 = max(n1, n2)
        m1 = min(n1, n2)
        d = m2-m1
        seriestype := :path
        subplot := (all ? 2 : 1)
        imi, ima = radiuslimits(d,seq1, seq2)
        if d == 0
            @series 1:n1
        else
            @series begin
                [imi ima]
            end
        end
    end

    if all
        if ndims(seq1) == 1
            # left line plot
            @series begin
                subplot := 1
                seq2, 1:n2
            end

            # bottom line plot
            @series begin
                subplot := 3
                1:n1, seq1
            end
        else
            # left line plot
            @series begin
                seriestype := :heatmap
                subplot := 1
                seq2'
            end

            # bottom line plot
            @series begin
                seriestype := :heatmap
                subplot := 3
                seq1
            end

        end
    end
end



@userplot MatchPlot

znorm(x) = (x = x.- mean(x); x ./= std(x))
using Statistics
@recipe function f(h::MatchPlot; transportcost=1, separation=2, ds=1, postprocess=nothing)
    x, y, D, i1, i2 = handleargs(h; transportcost=transportcost, postprocess=postprocess)
    x,y = znorm.((x,y))
    s1 = x .- separation
    s2 = y .+ separation

    @series begin
        s1
    end
    @series begin
        s2
    end
    @series begin
        primary := false
        linecolor --> :black
        seriesalpha --> 0.2
        i = fill(Inf, 1, length(i1))
        vec([i1'; i2'; i][:,1:ds:end]), vec([s1[i1]'; s2[i2]'; i][:,1:ds:end])
    end
end

@recipe function plot(r::DTWSearchResult)
    title --> "DTW-NN Search result"
    yguide --> "Distance"
    label --> round(r.cost, sigdigits=4)
    if length(r.dists) == 1
        @series begin
            seriestype := :scatter
            makersize --> 15
            markershape --> :x
            group := 1
            r.dists
        end
        @series begin
            seriestype := :hline
            linestyle := :dash
            primary := false
            group := 1
            r.dists
        end
    else
        @series r.dists
    end

    @series begin
        seriestype := :vline
        linestyle := :dash
        linecolor := :black
        primary := false
        [r.loc]
    end

end
