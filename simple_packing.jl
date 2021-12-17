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

function main(solvers, problems; timeout_factor=5, initial_timeout=1000)
    println("Burn-in test")
    for solver in solvers
        soln = solver_incremental(Problem(1, 1, 0, 0, 4, 4, [Rect(2, 2)]), solver_func=solver)
        @assert !isnothing(soln)
        println("Compilation time $(String(Symbol(solver))): ~$(round(soln.total_time / 10^6)) ms")
    end

    for (i, prob) in enumerate(problems)
        best_time = Inf
        for solver in solvers
            soln = solver_incremental(problems[i], solver_func=solver; timeout= isfinite(best_time) ? best_time*timeout_factor : initial_timeout)
            if isnothing(soln)
                println("\n|> $(String(Symbol(solver)))\t$(problems[i].num)\tNO SOLUTION")
            else
                verify = check_solution(prob, soln)
                println("\n|> $(String(Symbol(solver)))\t$(problems[i].num)\t$(verify)\t$(soln.bins)\t$(round(soln.total_time / 10^6)/1000) s")
                best_time = min(best_time, soln.total_time*10^-9)
            end
        end
    end

    println("Done!")
end

function problems_from_unibo(;filenames::Vector{String}=["Class_02.2bp"], do_first=0, skip_first=0)
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
    if skip_first > 0 && skip_first < length(problems)
        return problems[skip_first+1:length(problems)]
    end

    return problems
end

function check_solution(prob::Problem, soln::Solution)
    if isnothing(soln)
        return false
    end
    bins = soln.bins
    testarr = zeros(Int32, bins, prob.bin_h, prob.bin_w)
    ht = prob.bin_h
    wd = prob.bin_w

    for (part::Rect, pos::CartesianIndex{3}) in zip(prob.parts, soln.positions)
        bin, i, j = Tuple(pos)
        if bin <= 0 || bin > bins || i <= 0 || i + part.h - 1 > ht || j <= 0 || j + part.w - 1 > wd
            println("|>     Coordinate out of bounds $bin, $i, $j")
        end

        testarr[bin,i:(i+part.h-1),j:(j+part.w-1)] .+= 1
    end

    if length(soln.positions) != prob.num
        println("|>     Missing solutions!")
        return false
    elseif maximum(testarr) > 1
        println("|>     Overlapping bins.")
        return false
    else
        return true
    end
end

if true
    timeout_factor=10
    println("|> UNIBO")
    println("|>     Set timeout_factor=$timeout_factor")
    problems = problems_from_unibo(; skip_first=11)
    main([hough_and_cover, positions_and_covering], problems)
else
    println("|> TEST PROBLEMS")
    main([hough_and_cover, positions_and_covering], build_problems_basic())
end