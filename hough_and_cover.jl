using JuMP
import Cbc
include("types.jl")

# This implements our Hough and Cover technique, which is a derivative of
# the Positions and Covering (P&C) technique. In this technique, we use the Hough
# transform (equivalent to Positions step of P&C) and use a cumulative-sum map to
# ensure no collisions exist.
# This latter technique is taken from the 1999 Viola-Jones face detection paper that
# uses Haar Cascades to detect faces.
function hac_transform(model::Model, problem::Problem, bins::Int;
        rotations::Bool=false, verbose::Bool=true,
        incremental_runsum=true)
    @assert !rotations # Don't support rotations for now.

    np = length(problem.parts)
    wd = problem.bin_w
    ht = problem.bin_h

    # This is the first step in the transformation:
    @variable(model, houghmap[1:np, 1:bins, 1:ht, 1:wd], binary=true, base_name="houghmap")

    for k in 1:np
        maxy = (ht - problem.parts[k].h + 1)
        maxx = (wd - problem.parts[k].w + 1)
        # The object exists at exactly one position:
        @constraint(model, sum(houghmap[k,:,1:maxy,1:maxx]) == 1)
        # The object cannot exist at these positions
        @constraint(model, houghmap[k,:,(maxy + 1):ht,:] .== 0)
        @constraint(model, houghmap[k,:,:,(maxx + 1):wd] .== 0)
    end

    # Hough transform complete!
    # houghmap[k,i,j] is 1 iff object k is at position (i, j).
    # Now we compress this space to get a delta-map of shape [k, wd + 1, ht + 1], where
    # iff object k is at (i, j), the value at (k, i, j) is 1 and the position at (k, i+ht, j+wd) is -1.
    @variable(model, -1 <= deltamap[1:np, 1:bins, 1:ht, 1:wd] <= 1, integer=true, base_name="deltamap")
    for k in 1:np
        part = problem.parts[k]
        for q in 1:bins
            for i in 1:ht
                # The start coordinates an object must be for the value in this cell to be -1
                for j in 1:wd
                    start_i = i - part.h
                    start_j = j - part.w

                    # The inverse of the running sum transform:
                    cond_bottom = (start_i > 0) ? houghmap[k,q,start_i,j] : 0
                    cond_right = (start_j > 0) ? houghmap[k,q,i,start_j] : 0
                    cond_bottom_right = ((start_i > 0) && (start_j > 0)) ? houghmap[k,q,start_i,start_j] : 0
                    @constraint(model, deltamap[k,q,i,j] == houghmap[k,q,i,j] - cond_bottom - cond_right + cond_bottom_right)
                end
            end
        end
    end

    # Final step in the transform is to construct a running-sum map, similar to that used to
    # efficiently calculate Haarlike features in the 1999 Viola-Jones Face Detection paper.
    # The idea is that, at each position (i, j), to compute the running sum of all elements above- and to the left of it. (That is, the sum of [1:i,1:j]).
    # There are two strategies we can use to compute this:
    #  (1) Naive (which may be faster for this application?)
    #  (2) Incrementally, using:
    #      sum[i, j] := vals[i,j] + sum[i-1, j] + sum[i, j-1] - sum[i-1, j-1]

    # The number inside each cell in runsum is the number of objects that span over
    # that cell. We wish to find an assignment of objects that are non-overlapping,
    # so all cells must be constrained to be at most one.
    @variable(model, 0 <= runsum[1:bins, 1:ht, 1:wd] <= 1, integer=true, base_name="runningsummap")
    for q in 1:bins
        for i in 1:ht
            for j in 1:wd
                if incremental_runsum
                    this_cell = sum(deltamap[:,q,i,j])

                    above = (i == 1) ? 0 : ( runsum[q, i - 1, j] )
                    left  = (j == 1) ? 0 : ( runsum[q, i, j - 1] )
                    aboveleft = ((i == 1) || (j == 1)) ? 0 : ( runsum[q, i - 1,j - 1] )

                    @constraint(model, runsum[q,i,j] == this_cell + above + left - aboveleft)
                else
                    @constraint(model, runsum[q,i,j] == sum(deltamap[:,q,1:i,1:j]))
                end
            end
        end
    end

    # Now that we have created the problem, we solve it:
    optimize!(model)
    if verbose
        println("Solved $(termination_status(model)), primal $(primal_status(model))")
        println("Objective: ", objective_value(model))
    end
    Test.@test termination_status(model) == OPTIMAL
    Test.@test primal_status(model) == FEASIBLE_POINT

    # Extract the result:
    valmap = value.(houghmap)
    coords = Vector{CartesianIndex{3}}()

    for k in 1:np
        idxes = findall(x -> x == 1, valmap[k,:,:,:])
        @assert length(idxes) == 1
        idx = first(idxes)
        println("element $(k) at $(idx)")

        push!(coords, idx)
    end
    return coords
end

function solver_hac(prob::Problem; prior_num::Int = 1;
                    bin_search::string = "increment",
                    optimizer = Cbc.Optimizer)

    model = Model(optimizer)
    set_optimizer_attribute(model, "logLevel", 0)
    retval = hac_transform(model, prob, 1)
    return retval
end
