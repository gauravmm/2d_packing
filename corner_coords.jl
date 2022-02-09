using JuMP
import MathOptInterface

include("types.jl")


function corner_coords_construct(model::Model, parts::Vector{Rect}, ht::Int32, wd::Int32)
    np = length(parts)

    # Allocate variables
    @variable(model, 1 <= xmin[1:np] <= wd, integer = true, base_name = "xmin")
    @variable(model, 1 <= xmax[1:np] <= wd, integer = true, base_name = "xmax")
    @variable(model, 1 <= ymin[1:np] <= ht, integer = true, base_name = "ymin")
    @variable(model, 1 <= ymax[1:np] <= ht, integer = true, base_name = "ymax")

    println(parts)

    # 1(a) variables encode sizes
    for q in 1:np
        @assert parts[q].w + 1 == parts[q].h
        n = parts[q].w # Objects are of size (n, n+1)
        @constraint(model, n <= xmax[q] - xmin[q] + 1 <= n + 1)
        @constraint(model, n <= ymax[q] - ymin[q] + 1 <= n + 1)
        @constraint(model, (xmax[q] - xmin[q] + 1) + (ymax[q] - ymin[q] + 1) == 2n + 1)
    end

    M = 2 * (ht + wd) # We need a large constant for the disjunction.

    # 1(b) no-overlap constraint
    for q in 1:(np-1)
        for p in (q+1):np
            # anonymous variable syntax:
            z = @variable(model, [1:4], binary = true)
            # At least one of the following needs to hold:
            # xmin[q] > xmax[p] # p is left of q
            # xmin[p] > xmax[q] # p is right of q
            # ymin[q] > ymax[p] # p is above q
            # ymin[p] > ymax[q] # p is below q
            #
            # We implement this by creating an auxiliary variable z[i]
            # which is 1 if the corresponding constraint holds.
            # Since we can't use <, we instead add one to the LHS and use <=
            @constraint(model, -xmin[q] + xmax[p] + 1 <= M * (1 - z[1]))
            @constraint(model, -xmin[p] + xmax[q] + 1 <= M * (1 - z[2]))
            @constraint(model, -ymin[q] + ymax[p] + 1 <= M * (1 - z[3]))
            @constraint(model, -ymin[p] + ymax[q] + 1 <= M * (1 - z[4]))

            # We then constrain there to be at least one z that holds.
            @constraint(model, sum(z) >= 1)
        end
    end

    return (xmin, ymin, xmax, ymax)
end

function corner_coords_output(xmin, ymin, xmax, ymax)
    np = length(xmin)
    # Extract the result:
    coords = Vector{Union{CartesianIndex{3},Nothing}}(undef, np)
    rot = zeros(Bool, (np,))

    println()
    println("n\txmin)\txmax)\tymin)\tymax)")

    for k in 1:np
        coords[k] = CartesianIndex{3}((1, round(Int, ymin[k]), round(Int, xmin[k])))
        rot[k] = (xmax[k] - xmin[k]) > (ymax[k] - ymin[k])
        println("$k\t$(xmin[k])\t$(xmax[k])\t$(ymin[k])\t$(ymax[k])")
    end

    return coords, rot
end


"""
This implements the corner-coordinate algorithm.

It is a super simple algorithm that uses corner coordinates to prevent overlaps between pairs of rectangles. This is based on https://www.cs.cmu.edu/~mheule/15217-f21/assignments/assignment11.pdf
"""
function corner_coords(model::Model, problem::Problem, bins::Integer; timeout::Float64 = Inf, dump_model::Bool = false)
    construc_ns = time_ns()

    parts = problem.parts
    @assert problem.rotations
    @assert bins == 1


    xmin, ymin, xmax, ymax = corner_coords_construct(model, parts, problem.bin_h, problem.bin_w)
    if dump_model
        println("|>     Problem constructed in $((time_ns() - construc_ns)*10^-9) s")
        write_to_file(model, "model-cc-$(problem.problem_class)-$(problem.abs_inst).lp.gz")
        throw(OutputToFileDone())
    end

    if isfinite(timeout)
        println("Setting timeout $timeout")
        MOI.set(model, MOI.TimeLimitSec(), timeout)
    end

    # Now that we have created the problem, we solve it:
    optimize!(model)
    if termination_status(model) == OPTIMAL && primal_status(model) == FEASIBLE_POINT
        return corner_coords_output(value.(xmin), value.(ymin), value.(xmax), value.(ymax))
    end
    return nothing
end
