import Pkg
Pkg.activate(".")

using JuMP
using ArgParse
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
        Problem(1, 1, 1, 0, 0, 4, 4, [Rect(2, 2)], false)
        Problem(2, 1, 1, 0, 0, 4, 4, [Rect(3, 3)], false)
        Problem(3, 1, 1, 0, 0, 4, 4, [Rect(4, 4)], false)
        Problem(4, 1, 2, 0, 0, 4, 4, [Rect(2, 2), Rect(2, 2)], false)
        Problem(5, 1, 2, 0, 0, 4, 4, [Rect(3, 1), Rect(2, 2)], false)
        Problem(6, 1, 2, 0, 0, 4, 4, [Rect(2, 3), Rect(2, 2)], false)
        Problem(7, 1, 2, 0, 0, 4, 4, [Rect(2, 3), Rect(2, 3)], false)
        Problem(8, 1, 2, 0, 0, 4, 4, [Rect(3, 3), Rect(2, 3)], false)
    ]
end

function build_problems_unibo(fn; basepath="./")
    # Build a problem from the UNIBO files
    pat = normpath(joinpath(basepath, fn))
    seq = 1

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
            push!(problems, Problem(seq, problem_class, num, rel_ins, abs_ins, wbin, hbin, objs, false))
            seq+=1

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
        soln = solver_incremental(Problem(-1, 1, 1, 0, 0, 4, 4, [Rect(2, 2)], false), solver_func=solver)
        @assert !isnothing(soln)
        println("Compilation time $(String(Symbol(solver))): ~$(round(soln.total_time / 10^6)) ms")
    end

    println("|>     Set timeout_factor=$timeout_factor")
    for (i, prob) in enumerate(problems)
        best_time = Inf
        for solver in solvers
            soln = solver_incremental(problems[i], solver_func=solver; timeout= isfinite(best_time) ? best_time*timeout_factor : initial_timeout)
            if isnothing(soln)
                println("\n|> $(String(Symbol(solver)))\t$(problems[i].seq)\t$(problems[i].problem_class)\tNO SOLUTION")
            else
                verify = check_solution(prob, soln)
                println("\n|> $(String(Symbol(solver)))\t$(problems[i].seq)\t$(problems[i].problem_class)\t$(verify)\t$(soln.bins)\t$(round(soln.total_time / 10^6)/1000) s")
                best_time = min(best_time, soln.total_time*10^-9)
            end
        end
    end

    println("Done!")
end

function problems_from_unibo(;filenames::Vector{String}=[], do_first=0, skip_first=0)
    problems = Vector{}()
    for filename in filenames
        np = build_problems_unibo(filename)

        if isnothing(np)
            println("Error loading from $filename")
        else
            if do_first > 0 && do_first < length(np)
                np = np[1:do_first]
            end
            if skip_first > 0 && skip_first < length(np)
                np = np[skip_first+1:length(np)]
            end

            append!(problems, np)
            println("Loaded $(length(np)) problems from $filename")
        end
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

settings = ArgParseSettings()
@add_arg_table settings begin
    "filename"
        help="2bp file to load"
        required=true
        action=:store_arg
    "number"
        help="which question(s) to evaluate in the file"
        nargs='*'
        action = :store_arg
end

if true
    parsed_args = parse_args(ARGS, settings)
    problems = problems_from_unibo(;filenames=[parsed_args["filename"]])
    idxes = parsed_args["number"]
    if length(idxes) == 0
        println("No number specified, running all")
        idxes = collect(1:length(problems))
    end

    println("Requested problems $idxes")
    main([hough_and_cover, positions_and_covering], problems[idxes])

else
    println("|> TEST PROBLEMS")
    main([hough_and_cover, positions_and_covering], build_problems_basic())
end
