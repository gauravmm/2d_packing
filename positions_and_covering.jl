"""Implements the Positions and Covering (P&C) technique as an evaluation baseline.

We observe that the _Positions_ in P&C is equivalent to the Hough transform, and so we
re-use that primitive transform (and its inverse) in our Hough and Cover (H&C) approach.
"""

using JuMP

include("types.jl")


"""
Constructs and constraints a starting position map.
"""
function positions(model::Model, parts::Vector{Rect}, bins::Integer, ht::Integer, wd::Integer; rotations::Bool=false)
    @assert !rotations # We don't support rotations yet
    np = length(parts)

    # Generate the positions map
    @variable(model, houghmap[1:np, 1:bins, 1:ht, 1:wd], binary=true, base_name="houghmap")

    for k in 1:np
        # Maximum positions at which the top-left of the object can be placed:
        maxy = (ht - parts[k].h + 1)
        maxx = (wd - parts[k].w + 1)
        # The object exists at exactly one position:
        @constraint(model, sum(houghmap[k,:,1:maxy,1:maxx]) == 1)
        # The object cannot exist at these positions
        @constraint(model, houghmap[k,:,(maxy + 1):ht,:] .== 0)
        @constraint(model, houghmap[k,:,:,(maxx + 1):wd] .== 0)
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
