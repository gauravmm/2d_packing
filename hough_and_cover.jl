using JuMP
using types

# This implements our Hough and Cover technique, which is a derivative of
# the Positions and Covering (P&C) technique. In this technique, we use the Hough
# transform (equivalent to Positions step of P&C) and use a cumulative-sum map to
# ensure no collisions exist.
# This latter technique is taken from the 1999 Viola-Jones face detection paper that
# uses Haar Cascades to detect faces.

function parameter_space_transform(model :: Model, problem :: Problem ; rotations :: Bool=false)
    # This function constructs a positions map where the cell at (k, y, x) == 1
    # when positions[k].y == y and positions[k].x == x.
    assert(!rotations) # Don't support rotations for now.

    np = length(positions)
    wd = problem.wbin
    ht = problem.hbin

    # This is the first step in the transformation:
    @variable(model, houghmap[1:np, 1:ht, 1:wd], binary=true, base_name="houghmap")

    for k in 1:np
        maxy = (ht - problem.parts[k].h + 1)
        maxx = (wd - problem.parts[k].w + 1)
        # The object exists at exactly one position:
        @constraint(model, sum(houghmap[k,1:maxy,1:maxx]) == 1)
        # The object cannot exist at these positions
        @constraint(model, houghmap[k,(maxy+1):ht,:] == 0)
        @constraint(model, houghmap[k,:,(maxx+1):wd] == 0)
    end

    # Hough transform complete!
    # houghmap[k,i,j] is 1 iff object k is at position (i, j).
    # Now we compress this space to get a delta-map of shape [k, wd + 1, ht + 1], where
    # iff object k is at (i, j), the value at (k, i, j) is 1 and the position at (k, i+ht, j+wd) is -1.
    @variable(model, deltamap[1:np, 1:(ht+1), 1:(wd+1)], binary=true, base_name="deltamap")
    for k in 1:np
        part = problem.parts[k]
        for i in 1:(ht+1)
            for j in 1:(wd+1)
                # The start coordinates an object must be for the value in this cell to be -1
                start_i = i - ht
                start_j = j - wd
                # This is part of the constraint we are adding:
                cond_minus = (start_i > 0 && start_j > 0) ? houghmap[k,start_i,start_j] : 0
                cond_plus = (i <= ht && j <= wd) ? houghmap[k,i,j] : 0
                @constraint(model, deltamap[k,i,j] == cond_plus-cond_minus)
            end
        end
    end

    # Final step in the transform is to construct a subset-sum map, similar to that used to
    # efficiently calculate Haarlike features in the 1999 Viola-Jones Face Detection paper.
    # The idea is that, at each position (i, j), to compute the running sum of all elements above- and to the left of it. (That is, the sum of [1:i,1:j]).
    # There are two strategies we can use to compute this:
    #  (1) Naive, which may be faster for this application.
    #  (2) If we do this in a rowwise manner, we can get
    #      sum[i, j] as vals[i,j] + sum[i-1, j] + sum[i, j-1] - sum[i-1, j-1]

    return 
end


function example_knapsack(; verbose = true)
    profit = [5, 3, 2, 7, 4]
    weight = [2, 8, 4, 2, 5]
    capacity = 10
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:5], Bin)
    # Objective: maximize profit
    @objective(model, Max, profit' * x)
    # Constraint: can carry all
    @constraint(model, weight' * x <= capacity)
    # Solve problem using MIP solver
    optimize!(model)
    if verbose
        println("Objective is: ", objective_value(model))
        println("Solution is:")
        for i in 1:5
            print("x[$i] = ", value(x[i]))
            println(", p[$i]/w[$i] = ", profit[i] / weight[i])
        end
    end
    Test.@test termination_status(model) == OPTIMAL
    Test.@test primal_status(model) == FEASIBLE_POINT
    Test.@test objective_value(model) == 16.0
    return
end