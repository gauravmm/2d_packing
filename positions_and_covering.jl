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
end