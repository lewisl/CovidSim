using BenchmarkTools

function build_dict_of_arrays(dims, n)
    return Dict(i => rand(dims...) for i in 1:n)
end

function build_dict_of_indices(dims,n)
    dinds = Dict(zip(collect(1:n),randperm(n)))
    arr = rand(dims...,n)
    return dinds, arr
end


function create_and_time(dims, n)
    numdims = length(dims)
    @assert 2 <= numdims <= 3 "ndims must be 2 or 3"

    darrs = build_dict_of_arrays(dims, n)
    (dinds, arr) = build_dict_of_indices(dims, n)

    getone = rand(1:n, 1)[1]
    lst = dinds[getone]
    println(lst)

    t1 = @benchmark ans1 = $darrs[$getone]  # 9.5 nanoseconds for n = 3000
    if numdims == 2
        t2 = @benchmark ans2 = $arr[:,:,$lst]  # 
    elseif numdims == 3
        t2 = @benchmark ans2 = $arr[:,:,:, $lst]  # 2.7 microseconds for n = 3000
    end

    println(mean(t1), " ", mean(t2))
    # return darrs, dinds, arr
end