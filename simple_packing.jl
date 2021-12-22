import Pkg
Pkg.activate(".")

using JuMP
using ArgParse
import Test

include("problems.jl")
include("hough_and_cover.jl")

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

function check_solution(prob::Problem, soln::Solution)
    if isnothing(soln)
        return false
    end
    bins = soln.bins
    testarr = zeros(Int32, (bins, prob.bin_h, prob.bin_w))
    ht = prob.bin_h
    wd = prob.bin_w

    for (part::Rect, pos::CartesianIndex{3}) in zip(prob.parts, soln.positions)
        bin, i, j = Tuple(pos)
        if bin <= 0 || bin > bins || i <= 0 || i + part.h - 1 > ht || j <= 0 || j + part.w - 1 > wd
            println("|>     Coordinate out of bounds $bin, $i, $j")
        end

        testarr[bin,i:(i+part.h-1),j:(j+part.w-1)] .+= 1
    end

    if length(soln.positions) != length(prob.parts)
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
    "basic"
        help="basic tests"
        action=:command
    "basicrot"
        help="basic rotation tests"
        action=:command
    "unibo"
        help="run examples from the UNIBO dataset"
        action=:command
    "nqsq"
        help="not quite squares (with rotations!)"
        action=:command
end

@add_arg_table settings["unibo"] begin
    "filename"
        help="2bp file to load"
        required=true
        action=:store_arg
    "number"
        help="which question(s) to evaluate in the file"
        nargs='*'
        arg_type=Int
        action = :store_arg
end

@add_arg_table settings["nqsq"] begin
    "binsize"
        help="the size of the bin, as (binsize-1, binsize)"
        arg_type=Int
        required=true
        action=:store_arg
    "maxn"
        help="the largest not quite square, as (maxn-1, maxn)"
        arg_type=Int
        required=true
        action = :store_arg
end


parsed_args = parse_args(ARGS, settings)
if parsed_args["%COMMAND%"] == "basic"
    println("|> TEST PROBLEMS")
    main([hough_and_cover, positions_and_covering], build_problems_basic())

elseif parsed_args["%COMMAND%"] == "basicrot"
    println("|> TEST ROTATION PROBLEMS")
    main([hough_and_cover], build_problems_basicrot())

elseif parsed_args["%COMMAND%"] == "unibo"
    problems = build_problems_unibo(;filenames=[parsed_args["filename"]])
    idxes = parsed_args["number"]
    if length(idxes) == 0
        println("No number specified, running all")
        idxes = collect(1:length(problems))
    end

    println("Requested problems $idxes")
    main([hough_and_cover, positions_and_covering], problems[idxes])

elseif parsed_args["%COMMAND%"] == "nqsq"
    binsize = parsed_args["binsize"]
    binsize = parsed_args["maxn"]
    problems = build_problems_nqsq(parsed_args["binsize"], parsed_args["maxn"])

    main([hough_and_cover, positions_and_covering], problems[idxes])

end
