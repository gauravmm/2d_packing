# Preprocess the problem based on Section 4 of Combinatorial Benders Decomposition (BD) paper.

using Test
using Graphs

include("types.jl")

#
# Bin and Object scaling clique
#

function preprocess(prob::Problem)
    # Repeat until no further changes.
    while true
        is_changed, prob = shrink_and_expand(prob)
        if !is_changed
            return prob
        end
    end
end

function preprocess⁻¹(prob::Problem, soln::Solution) :: Solution
    # The shrink_and_expand step does not require any changes to the solution.
    return soln
end

function shrink_and_expand(prob::Problem)::Tuple{Bool, Problem}
    is_changed = false

    # First we compute the possible widths or heights achievable using any combination
    # of objects
    pos_wd = compute_possible_lengths([p.w for p in prob.parts], prob.bin_w)
    pos_ht = compute_possible_lengths([p.h for p in prob.parts], prob.bin_h)
    # pos_wd is a [len(parts) + 1, bin_w] array where [i, j] is true if the value j can be reached using
    # a sum of the widths of any parts except i. The final row does not exclude any element.

    # Alvarez-Valdes et al. (2009), via BD paper.
    # We shrink the bin size (along each axis) to the smallest possible sum of object sizes:
    bin_w = findlast(pos_wd[end,:])
    bin_h = findlast(pos_ht[end,:])
    is_changed |= (prob.bin_w != bin_w) | (prob.bin_h != bin_h)

    # Carlier et al. (2007), via BD paper.
    # We grow object j along each dimension if there are no sums of other object lengths such that
    # it can exactly occupy the entire width of the (shrunk) bin.
    parts::Vector{Rect} = copy(prob.parts)
    for (i, p) in enumerate(parts)
        w = compute_growth(pos_wd[i,:], bin_w, p.w)
        h = compute_growth(pos_ht[i,:], bin_h, p.h)
        parts[i] = Rect(w, h)
        is_changed |= (w != p.w && h != p.h)
    end

    # If anything is changed, then return an updated object.
    return (is_changed, Problem(prob.seq, prob.problem_class, prob.num,
            prob.rel_inst, prob.abs_inst,
            bin_w, bin_h, parts,
            prob.rotations))
end

function compute_growth(sizes::BitVector, bin_len::Integer, obj_len::Integer)
    if obj_len == bin_len
        return obj_len
    end

    # Z is the difference between the current object size and the maximum size it can be
    # without changing what objects can be fit along a single line in the bin.
    z = findlast(sizes[1:(bin_len - obj_len)])
    if z === nothing
        # If there are no other objects that can come between this object and the
        # edge of the bin, then set the object to fill the entire bin.
        z = 0
    end
    return bin_len - z
end

function compute_possible_lengths(values::Vector{T}, max_val::Integer) where T <: Integer
    # Compute the possible lengths that can be achieved up to (and including) max_val
    # Return a list of length(values) + 1 elements where the element at position i is a mask
    # excluding values[i]

    rv = falses((length(values) + 1, max_val))

    for (i, v) in enumerate(values)
        # mask controls which rows to write to:
        mask = trues(length(values) + 1)
        mask[i] = false

        # Any existing possible sums can be extended by adding v to it
        if v < max_val
            rv[mask,(v+1):end] .|= rv[mask,1:(end-v)]
        end
        # And we can also directly reach v.
        rv[mask, v] .= true
    end

    return rv
end

# @testset "possible lengths" begin
#     @test all(compute_possible_lengths([1, 2, 4], 8) .== [0 1 0 1 0 1 0 0; 1 0 0 1 1 0 0 0; 1 1 1 0 0 0 0 0; 1 1 1 1 1 1 1 0])
#     @test all(compute_possible_lengths([2, 3, 5], 5) .== [0 0 1 0 1; 0 1 0 0 1; 0 1 1 0 1; 0 1 1 0 1])
# end;

#
# Incompatible Items and Maximal Incompatible Clique
#

function conflicts(parts::Vector{Rect}, bin_h::Integer, bin_w::Integer)
    g = SimpleGraph(length(parts))
    np = length(parts)

    if np <= 1
        return (g, [])
    end

    for i in 1:(np - 1)
        for j in (i+1):np
            if (parts[i].w + parts[j].w) > bin_w && (parts[i].h + parts[j].h) > bin_h
                # Undirected graph:
                add_edge!(g, i, j)
            end
        end
    end

    cliques = sort!(maximal_cliques(g), by=length, rev=true)
    return (g, cliques[1])
end
