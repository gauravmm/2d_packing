using JuMP
import GLPK
import Test

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

function popline(f::IOStream)
    line = readline(f)
    if isempty(strip(line))
        return
    end

    left = strip(line[1:5])
    right = strip(line[6:10])

    return [parse(Int32, left) (isempty(right) ? -1 : parse(Int32, right))]
end

function build_problems_unibo(fn; basepath = "data/unibo/")
    # Build a problem from the UNIBO files
    pat = normpath(joinpath(basepath, fn))

    problems = Vector{Problem}()
    # next = Problem()

    open(pat) do f
        while !eof(f)
            problem_class, _ = popline(f)
            num, _ = popline(f)
            rel_ins, abs_ins = popline(f)
            hbin, wbin = popline(f)
            objs = Vector{Rect}()

            for _ in 1:num
                h, w = popline(f)
                push!(objs, Rect(w, h))
            end
            push!(problems, Problem(problem_class, num, rel_ins, abs_ins, wbin, hbin, objs))

            if !eof(f)
                popline(f)
            end
        end
    end

    return problems
end

function main(; from_files=["Class_01.2bp"], do_first=1)
    for filename in from_files
        problems = build_problems_unibo(filename)
        print("Loaded $(length(problems)) problems from $filename")

    end
end

main()
