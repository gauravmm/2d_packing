import Pkg
Pkg.activate(".")

using JuMP
import GLPK
import Test

include("hough_and_cover.jl")

function popline(f::IOStream)
    line = readline(f)
    if isempty(strip(line))
        return
    end

    left = strip(line[1:5])
    right = strip(line[6:10])

    return [parse(Int32, left) (isempty(right) ? -1 : parse(Int32, right))]
end

function build_problems_basic()
    return [
        Problem(1, 1, 0, 0, 4, 4, [Rect(2, 2)])
        Problem(1, 1, 0, 0, 4, 4, [Rect(3, 3)])
        Problem(1, 1, 0, 0, 4, 4, [Rect(4, 4)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(2, 2), Rect(2, 2)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(3, 1), Rect(2, 2)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(2, 3), Rect(2, 2)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(2, 3), Rect(2, 3)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(3, 3), Rect(2, 3)])
    ]
end

function build_problems_unibo(fn; basepath = "data/unibo/")
    # Build a problem from the UNIBO files
    pat = normpath(joinpath(basepath, fn))

    problems = Vector{Problem}()
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

function main(solver, problems)
    println("Burn-in test")
    soln = solver_incremental(Problem(1, 1, 0, 0, 4, 4, [Rect(2, 2)]), solver_func=solver)
    @assert !isnothing(soln)
    println("Compilation time ~$(round(soln.total_time / 10^6)) ms")

    for (i, prob) in enumerate(problems)
        soln = solver_incremental(problems[i], solver_func=solver)
        if isnothing(soln)
            println("|> $(i)\tNO SOLUTION")
        else
            verify = check_solution(prob, soln)
            println("|> $(i)\t$(verify)\t$(soln.bins)\t$(round(soln.total_time / 10^6)) ms")
        end
    end

    println("Done!")
end

function problems_from_unibo(;filenames::Vector{String}=["Class_02.2bp"], do_first=0)
    problems = Vector{}()
    for filename in filenames
        np = build_problems_unibo(filename)
        if isnothing(np)
            println("Error loading from $filename")
        else
            append!(problems, np)
            println("Loaded $(length(np)) problems from $filename")
        end
    end
    if do_first > 0 && do_first < length(problems)
        return problems[1:do_first]
    end
    return problems
end

function check_solution(prob, soln)
    return true
end

problems = problems_from_unibo(do_first=25)
println("|> POSITIONS AND COVERING")
main(positions_and_covering, problems)
println("|> HOUGH AND COVER")
main(hough_and_cover, problems)
