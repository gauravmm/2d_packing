using JuMP
import Cbc
include("types.jl")

# This implements our Hough and Cover technique, which is a derivative of
# the Positions and Covering (P&C) technique. In this technique, we use the Hough
# transform (equivalent to Positions step of P&C) and use a cumulative-sum map to
# ensure no collisions exist.
# This latter technique is taken from the 1999 Viola-Jones face detection paper that
# uses Haar Cascades to detect faces.
function hac_transform(model::Model, problem::Problem ; rotations::Bool=false, verbose::Bool=true)
    # This function performs the hough-and-cover transform, producing

    # This function constructs a positions map where the cell at (k, y, x) == 1
    # when positions[k].y == y and positions[k].x == x.
    @assert !rotations # Don't support rotations for now.

    np = length(problem.parts)
    wd = problem.bin_w
    ht = problem.bin_h

    # This is the first step in the transformation:
    @variable(model, houghmap[1:np, 1:ht, 1:wd], binary = true, base_name = "houghmap")

    for k in 1:np
        maxy = (ht - problem.parts[k].h + 1)
        maxx = (wd - problem.parts[k].w + 1)
        # The object exists at exactly one position:
        @constraint(model, sum(houghmap[k,1:maxy,1:maxx]) == 1)
        # The object cannot exist at these positions
        @constraint(model, houghmap[k,(maxy + 1):ht,:] .== 0)
        @constraint(model, houghmap[k,:,(maxx + 1):wd] .== 0)
    end

    # Hough transform complete!
    # houghmap[k,i,j] is 1 iff object k is at position (i, j).
    # Now we compress this space to get a delta-map of shape [k, wd + 1, ht + 1], where
    # iff object k is at (i, j), the value at (k, i, j) is 1 and the position at (k, i+ht, j+wd) is -1.
    @variable(model, -1 <= deltamap[1:np, 1:(ht + 1), 1:(wd + 1)] <= 1, integer=true, base_name = "deltamap")
    for k in 1:np
        part = problem.parts[k]
        for i in 1:(ht + 1)
            # The start coordinates an object must be for the value in this cell to be -1
            for j in 1:(wd + 1)
                start_i = i - part.h
                start_j = j - part.w
                # This is part of the constraint we are adding:
                cond_minus = ((start_i > 0) && (start_j > 0)) ? houghmap[k,start_i,start_j] : 0
                cond_plus = ((i <= ht) && (j <= wd)) ? houghmap[k,i,j] : 0
                @constraint(model, deltamap[k,i,j] == cond_plus - cond_minus)
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
    @variable(model, 0 <= runsum[1:(ht + 1), 1:(wd + 1)] <= np, integer = true, base_name = "runningsummap")
    for i in 1:(ht + 1)
        for j in 1:(wd + 1)
            """
            this_cell = sum(deltamap[:,i,j])

            above = (i == 1) ? 0 : ( runsum[i - 1, j] )
            left  = (j == 1) ? 0 : ( runsum[i, j - 1] )
            aboveleft = ((i == 1) || (j == 1)) ? 0 : ( runsum[i - 1,j - 1] )

            @constraint(model, runsum[i,j] == this_cell + above + left - aboveleft)
            """
            @constraint(model, runsum[i,j] == sum(deltamap[:,1:i,1:j]))
        end
    end

    #@objective(model, Min, sum(runsum));

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
    coords = Vector{CartesianIndex{2}}()

    for k in 1:np
        idxes = findall(x -> x == 1, valmap[k,:,:])
        @assert length(idxes) == 1
        idx = first(idxes)
        println("element $(k) at $(idx)")

        push!(coords, idx)
    end
    return coords
end

function hac_solve(prob::Problem)
    model = Model(GLPK.Optimizer)
    # model = Model(Cbc.Optimizer)
    retval = hac_transform(model, prob)

    return retval
end
