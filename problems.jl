include("types.jl")

#
# Basic problem
#
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
        Problem(8, 1, 2, 0, 0, 3, 3, [Rect(1, 1) for _ in 1:9], false)
    ]
end

function build_problems_basicrot()
    return [
        Problem(8, 1, 2, 0, 0, 3, 3, [Rect(3, 1), Rect(1, 3)], true)
        Problem(8, 1, 2, 0, 0, 3, 3, [Rect(3, 2), Rect(1, 3)], true)
        Problem(8, 1, 2, 0, 0, 3, 3, [Rect(1, 3), Rect(1, 3), Rect(3, 1)], true)
    ]
end


#
# UNIBO reference problems
#
function popline(f::IOStream)
    line = readline(f)
    if isempty(strip(line))
        return
    end

    left = strip(line[1:5])
    right = strip(line[6:10])

    return [parse(Int32, left) (isempty(right) ? -1 : parse(Int32, right))]
end

function unibo_parse(fn; basepath="./")
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

function build_problems_unibo(;filenames::Vector{String}=[], do_first=0, skip_first=0)
    problems = Vector{}()
    for filename in filenames
        np = unibo_parse(filename)

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


#
# Not-quite Square problems
#
function build_problems_nqsq(b::Int, n::Int)
    @assert n > 1
    @assert b > 1

    return [
        Problem(0, b, n, 0, 0, b-1, b, [Rect(i-1, i) for i in 2:n], true)
    ]
end
