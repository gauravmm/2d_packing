"""Implements the Positions and Covering (P&C) technique as an evaluation baseline.

We observe that the _Positions_ in P&C is equivalent to the Hough transform, and so we
re-use that primitive transform (and its inverse) in our Hough and Cover (H&C) approach.
"""

using JuMP
import Gurobi

include("types.jl")


"""
Constructs and constraints a starting position map.
"""
function positions_constraints(model, houghmap, k, ht, wd, part_h, part_w)
    # Maximum positions at which the top-left of the object can be placed:
    maxy = (ht - part_h + 1)
    maxx = (wd - part_w + 1)
    # The object cannot exist at these positions
    @constraint(model, houghmap[k,:,(maxy + 1):ht,:] .== 0)
    @constraint(model, houghmap[k,:,:,(maxx + 1):wd] .== 0)

    # Instead of using this constraint:
    # @constraint(model, sum(houghmap[k,:,1:maxy,1:maxx]) == 1)
    # We compute the number of positions the object exists at:
    return sum(houghmap[k,:,1:maxy,1:maxx])
end

function positions(model::Model, parts::Vector{Rect}, bins::Integer, ht::Integer, wd::Integer; rotations::Bool=false, strengthen::Bool=true)
    np = length(parts)

    # Generate the positions map
    @variable(model, houghmap[1:np, 1:bins, 1:ht, 1:wd], binary=true, base_name="houghmap")

    for k in 1:np
        num_occur = positions_constraints(model, houghmap, k, ht, wd, parts[k].h, parts[k].w)
        @constraint(model, num_occur == 1)
    end

    # If strengthen is true, then we implement the "strengthening the convex polytope"
    # tricks from P&C. Of the three presented tricks, the only one not already covered in
    # our implementation is the constraint limiting the total area of all items in a bin
    # to be no more than the total area of a bin.
    if strengthen
        for b in 1:bins
            area_of_objects_in_bin = [sum(houghmap[k,b,:,:])*(parts[k].h*parts[k].w) for k in 1:np]
            @constraint(model, sum(area_of_objects_in_bin) <= ht*wd)
        end
    end

    return houghmap
end

function positions(model::Model, parts::Vector{Rect}, bins::Integer, ht::Integer, wd::Integer; rotations::Bool=true, strengthen::Bool=true)
    np = length(parts)

    # Generate the positions map
    @variable(model, houghmap[1:(rotations ? 2*np : np), 1:bins, 1:ht, 1:wd], binary=true, base_name="houghmap")

    for k in 1:np
        num_occur_orient = positions_constraints(model, houghmap, 2*k, ht, wd, parts[k].h, parts[k].w)
        num_occur_rotate = positions_constraints(model, houghmap, 2*k+1, ht, wd, parts[k].w, parts[k].h)
        @constraint(model, (num_occur_orient + num_occur_rotate) == 1)
    end

    # If strengthen is true, then we implement the "strengthening the convex polytope"
    # tricks from P&C.
    if strengthen
        for b in 1:bins
            area_of_objects_in_bin = [sum(houghmap[((2*k):(2*k+1)),b,:,:])*(parts[k].h*parts[k].w) for k in 1:np]
            @constraint(model, sum(area_of_objects_in_bin) <= ht*wd)
        end
    end

    return houghmap
end

"""
Given a solved position map, recover the indices of each object.
"""
function positions⁻¹(valmap::Array{Float64,4}; rotations::Bool=false)
    np, _, _, _ = size(valmap)
    # Extract the result:
    @assert !rotations
    coords = Vector{CartesianIndex{3}}(undef, np)

    # TODO: Write this using a single call to findall and sorting by first coordinate.
    for k in 1:np
        idxes = findall(x -> x >= .999, valmap[k,:,:,:])
        @assert length(idxes) == 1
        coords[k] = first(idxes)
    end

    return coords
end

function positions⁻¹(valmap::Array{Float64,4}; rotations::Bool=true)
    np, _, _, _ = size(valmap)
    # Extract the result:
    @assert !rotations
    coords = Vector{CartesianIndex{4}}(undef, np)

    for k in 1:np
        idxes = findall(x -> x >= .999, valmap[((2*k):(2*k+1)),:,:,:])
        @assert length(idxes) == 1
        coords[k] = first(idxes)
    end

    return coords
end

"""
Given a problem, finds initial bounds on the number of bins:
"""
function bin_bounds(problem::Problem)
    area_per_bin = problem.bin_h * problem.bin_w
    area_objects = sum([r.w * r.h for r in problem.parts])

    # This is the lower bound imposed by the total covered area:
    # cld is ceiling(division(nanny, donkey)):
    lower_bound = Int(cld(area_objects, area_per_bin))

    # Upper bound:
    # A trivial upper bound is four times the lower bound, because it is
    # guaranteed to be higher than the case where (under optimal packing)
    # each bin is slightly more than half full along each axis.
    # In higher-dimensional problems, this gets exponentially worse.
    upper_bound = 4*lower_bound

    return (lower_bound, upper_bound)
end

"""Implements the original covering part of the P&C algorithm, as in the paper.

This method solves very slowly due to a completely avoidable combinatorial explosion
caused by computing the covering set for each object separately.
"""
function covering(model::Model, parts::Vector{Rect}, houghmap; rotations::Bool=false)
    @assert rotations == false

    (np, bins, ht, wd) = size(houghmap)
    # Now we compute a covering of shape [np, bins, ht, wd], where
    # iff object k is at (i, j) in bin q, the value at (k, q, i:i+h-1, j:j+w-1) is 1
    @variable(model, covering[1:np, 1:bins, 1:ht, 1:wd] <= 1, integer=true, base_name="covermap")

    for k in 1:np
        part = parts[k]
        maxy = (ht - part.h + 1)
        maxx = (wd - part.w + 1)

        for q in 1:bins
            for i in 1:maxy
                for j in 1:maxx
                    # We use >= so that the cell in the covering is 1 if the original position
                    # is anywhere such that the object overlaps.
                    @constraint(model, covering[k,q,i:(i+part.h-1),j:(j+part.w-1)] .>= houghmap[k,q,i,j])
                end
            end
        end
    end

    # Now we add the no-overlapping covering:
    for q in 1:bins
        for i in 1:ht
            for j in 1:wd
                @constraint(model, sum(covering[:,q,i,j]) <= 1)
            end
        end
    end
    return covering
end

"""
This implements Positions and Covering (P&C)
"""
function positions_and_covering(model::Model, problem::Problem, bins::Int;
        rotations::Bool=false, timeout::Float64=Inf)
    @assert !rotations # Don't support rotations for now.

    # houghmap contains x^i_{jk} from the P&C paper
    houghmap = positions(model, problem.parts, bins, problem.bin_h, problem.bin_w)
    # covermap contains C, the correspondence matrix
    covermap = covering(model, problem.parts, houghmap)

    # Now that we have created the problem, we solve it:
    if isfinite(timeout)
        println("Setting timeout $timeout")
        MOI.set(model, MOI.TimeLimitSec(), timeout)
    end
    optimize!(model)
    if termination_status(model) == OPTIMAL && primal_status(model) == FEASIBLE_POINT
        return positions⁻¹(value.(houghmap))
    end
    return nothing
end

"""
Run P&C (or H&C) across a range of bins, increasing the bin size on each failure. 
"""
function solver_incremental(prob::Problem; known_bins::Int = 0,
                    solver_func=positions_and_covering,
                    optimizer=Gurobi.Optimizer,
                    timeout=Inf)
    lb, ub = bin_bounds(prob)
    if known_bins > 0
        lb, ub = known_bins
    end

    start_time = time_ns()
    bins = lb
    while bins <= ub
        last_time = time_ns()

        time_spent = (last_time - start_time)*10^-9
        if isfinite(timeout)
            if time_spent >= timeout
                println("|>     Timed out after $time_spent sec! At bin: $bins - 1")
                return nothing
            end
        end

        model = Model(optimizer)
        retval = solver_func(model, prob, bins; timeout=timeout<=0 ? -1 : timeout - time_spent)

        if !isnothing(retval)
            # We have a solution!
            end_time = time_ns()
            return Solution(true, bins, retval, end_time - start_time, end_time - last_time)
        end
        bins += 1
    end
    println("Oops! We failed to find a solution with bins in ($lb, $ub)")
    return nothing
end
