using JuMP
import MathOptInterface

include("types.jl")
include("preprocessor.jl")
include("positions_and_covering.jl")

"""
Construct a delta-map of the space, such that:
    runningsum(deltamap(positions)) == covering(positions)

By factoring the problem this way, we can establish the same objective with
fewer conditions.
"""
function delta_transform(model::Model, parts::Vector{Rect}, houghmap)
    (np, bins, ht, wd) = size(houghmap)
    # Now we compute a delta-map of shape [np, bins, ht, wd], where
    # iff object k is at (i, j) in bin q, the value:
    #    at (k, q, i, j) is 1,
    #    at (k, q, i+part[k].ht, j) is -1,
    #    at (k, q, i, j+part[k].wd) is -1,
    #    at (k, q, i+part[k].ht, j+part[k].wd) is 1
    # This is the _inverse_ transformation to that of the running sum.
    # (You can think of this is one dimension, and the solution here is the convolution of
    # that along two dimensions.)
    @variable(model, -1 <= deltamap[1:np, 1:bins, 1:ht, 1:wd] <= 1, integer=true, base_name="deltamap")
    # While technically the size of deltamap should be [np, bins, ht+1, wd+1], here we
    # can ignore the last column because the runsum of the output must be zero there.
    for k in 1:np
        part = parts[k]
        for q in 1:bins
            for i in 1:ht
                # The start coordinates an object must be for the value in this cell to be -1
                for j in 1:wd
                    start_i = i - part.h
                    start_j = j - part.w

                    # The inverse of the running sum transform.
                    # We are computing the position at which the object must begin to 
                    # impart +1 or -1 to this cell.
                    cond_bottom = (start_i > 0) ? houghmap[k,q,start_i,j] : 0
                    cond_right = (start_j > 0) ? houghmap[k,q,i,start_j] : 0
                    cond_bottom_right = ((start_i > 0) && (start_j > 0)) ? houghmap[k,q,start_i,start_j] : 0

                    @constraint(model, deltamap[k,q,i,j] == houghmap[k,q,i,j] - cond_bottom - cond_right + cond_bottom_right)
                end
            end
        end
    end

    return deltamap
end


abstract type RunningSumMode end
struct Naïve <: RunningSumMode end
struct Incremental <: RunningSumMode end

"""
At each position (i, j), compute the running sum of all elements above and left of it. (That is, the sum of [1:i,1:j]).

The number inside each cell in runsum is the number of objects that span over
that cell. We wish to find an assignment of objects that are non-overlapping,
and we achieve that by constraining all cells to be at most one.
"""
function runningsum(model::Model, deltamap; mode::Naïve)
    @assert false # Don't use this method!

    # Compute the running sum using the naïve method.
    _, bins, ht, wd = size(deltamap)
    # The <= 1 constraint is what prevents objects overlapping:
    @variable(model, 0 <= runsum[1:bins, 1:ht, 1:wd] <= 1, integer=true, base_name="runningsummap")
    for q in 1:bins
        for i in 1:ht
            for j in 1:wd
                @constraint(model, runsum[q,i,j] == sum(deltamap[:,q,1:i,1:j]))
            end
        end
    end

    return runsum
end

function runningsum(model::Model, deltamap; mode::Incremental)
    np, bins, ht, wd = size(deltamap)
    # Compute this incrementally, using:
    #      sum[i, j] := vals[i,j] + sum[i-1, j] + sum[i, j-1] - sum[i-1, j-1]
    # This is a well-known technique, but I first read it in the 1999 Viola-Jones
    # Face-recognition paper. There it was used to efficiently compute cascades
    # of Haarlike filters.
    @variable(model, 0 <= runsum[1:bins, 1:ht, 1:wd] <= 1, integer=true, base_name="runningsummap")
    for q in 1:bins
        for i in 1:ht
            for j in 1:wd
                this_cell = sum(deltamap[:,q,i,j])

                above = (i == 1) ? 0 : ( runsum[q, i - 1, j] )
                left  = (j == 1) ? 0 : ( runsum[q, i, j - 1] )
                aboveleft = ((i == 1) || (j == 1)) ? 0 : ( runsum[q, i - 1, j - 1] )

                @constraint(model, runsum[q,i,j] == this_cell + above + left - aboveleft)
            end
        end
    end

    return runsum
end

"""
This implements Hough and Cover (H&C), which is a derivative of Positions and Covering (P&C).

We use the Hough transform (equivalent to Positions step of P&C) and use a cumulative-sum map to ensure no collisions exist. The use of a running-sum map (and appropriate pre-transformation) allows us to greatly reduce the number of conditions.
"""
function hough_and_cover(model::Model, problem::Problem, bins::Integer;
        runsummode::RunningSumMode=Incremental(), timeout::Float64=Inf,
        prep::Bool=true)

    parts = problem.parts
    if problem.rotations
        # Expand parts to include rotations:
        parts = [parts; [Rect(p.h, p.w) for p in parts]]
    end

    houghmap, supply = positions(model, parts, bins, problem.bin_h, problem.bin_w)
    if problem.rotations
        np = length(problem.parts)
        # We want exactly one of the oriented object and its rotation:
        @constraint(model, supply[1:np] .+ supply[np+1:2*np] .== 1)
    else
        @constraint(model, supply .== 1)
    end

    deltamap = delta_transform(model, parts, houghmap)
    runsum = runningsum(model, deltamap; mode=runsummode)

    if isfinite(timeout)
        println("Setting timeout $timeout")
        MOI.set(model, MOI.TimeLimitSec(), timeout)
    end
    # Now that we have created the problem, we solve it:
    optimize!(model)
    if termination_status(model) == OPTIMAL && primal_status(model) == FEASIBLE_POINT
        valmap = positions⁻¹(value.(houghmap))
        if problem.rotations
            return gather_rotations(problem, valmap)
        else
            return valmap, nothing
        end
    end
    return nothing
end

function gather_rotations(problem::Problem, valmap)
    np = length(problem.parts)
    vm = Vector{CartesianIndex{3}}(undef, np)
    rot= zeros(Bool, (np,))

    for k in 1:np
        if !isnothing(valmap[k])
            vm[k] = valmap[k]
            rot[k] = false
        else
            vm[k] = valmap[k + np]
            rot[k] = true
        end
    end

    return vm, rot
end
