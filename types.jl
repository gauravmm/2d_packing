# Types for this.
using JuMP

struct Rect
    w :: Int32
    h :: Int32
end

struct Problem
    problem_class::Int32
    num::Int32
    rel_inst::Int32
    abs_inst::Int32
    bin_w::Int32
    bin_h::Int32
    parts::Vector{Rect}
end
