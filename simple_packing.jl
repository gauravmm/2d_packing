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

function main(; from_unibo=["Class_02.2bp"], do_first=1)
    for filename in from_unibo
        problems = build_problems_unibo(filename)
        println("Loaded $(length(problems)) problems from $filename")

        for i in 1:( do_first > 0 ? do_first : length(problems))
            println("Solving $(i)...")
            solver_hac(problems[i])
            println("Done!")
        end
    end
end

#main()

function check_solution(prob, soln)
    println(soln)
end

function do_basic_tests(solver)
    # Create and test some simple problems:
    problems = [
        Problem(1, 1, 0, 0, 4, 4, [Rect(2, 2)])
        Problem(1, 1, 0, 0, 4, 4, [Rect(3, 3)])
        Problem(1, 1, 0, 0, 4, 4, [Rect(4, 4)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(2, 2), Rect(2, 2)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(3, 1), Rect(2, 2)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(2, 3), Rect(2, 2)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(2, 3), Rect(2, 3)])
        Problem(1, 2, 0, 0, 4, 4, [Rect(3, 3), Rect(2, 3)])
    ]

    for (i, prob) in enumerate(problems)
        soln = solver(prob)
        if isnothing(soln)
            println("|> $(i)\tNO SOLUTION")
        else
            verify = check_solution(prob, soln)
            println("|> $(i)\t$(verify)\t$(soln.bins)\t$(round(soln.total_time / 10^6)) ms")
        end
    end

    println("Done!")
end

do_basic_tests(solver_hac)
