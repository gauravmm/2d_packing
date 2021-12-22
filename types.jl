# Types for this.
using JuMP

struct Rect
    w :: Int32
    h :: Int32
end

struct Problem
    seq::Int32
    problem_class::Int32
    num::Int32
    rel_inst::Int32
    abs_inst::Int32
    bin_w::Int32
    bin_h::Int32
    parts::Vector{Rect}
    rotations::Bool
end

struct Solution
    solved::Bool
    bins::Int32
    positions::Vector{CartesianIndex{3}}
    rotations::Union{Vector{Int}, Nothing}
    total_time::Float64
    last_time::Float64
end
