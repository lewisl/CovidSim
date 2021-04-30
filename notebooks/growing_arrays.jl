# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# ## Experiments Growing Arrays

# %% [markdown]
# #### Conclusion
# *for growing:* We want to pre-allocate a bucket big enough to hold the largest expected number of items and hope we make it big enough. Then, insert new items in place.
#
# *for deleting:* We want to have a boolean that matches the maximum size of the bucket. Values need to be initialized to *false*. Then, we set each item to true for deleting from the bucket. Then we will have to expand the bucket again and expand the filt of the bucket back to maximum length. But, expanding by resize! is really slow.  So, deleting really doesn't work.  Instead, we want to zero out the elements we don't want any more. But, then we need to sort to put the zeros at the end (or the beginning--but the end is easier).
#
# *for zeroing and sorting:* Both operations are pretty slow.
#
# *conclusion:* We have to conclude that the bucket strategy doesn't work if the bucket is big and/or the number of items we have in the bucket changes each time. The fastest thing is to change values of the vector in place and then use findall to generate indexes to each value.

# %%
using BenchmarkTools
using StatsBase
using Random

# %%
# naive push!

function grow_by_push(arr, item, n_add=100)
#     arr = copy(arr)
    for i in 1:n_add
        push!(arr, item)
    end
    
    return arr
end
        
    

# %%
# append and return, append!

function append_and_return(arr1, arr2)
    arr_in = copy(arr1)
    append!(arr_in, arr2)
    return arr_in
end

# %%
function append(arr1, arr2)
    append!(arr1, arr2)
    return
end

# %%
function insertat_loop!(arr1, arr2, position)
    insat = position
    for i in 1:length(arr2)
        arr1[insat] = arr2[i]
        insat += 1
    end
end

# %%
function insertat!(arr1, arr2, position)
    pos_end = position + length(arr2) -1
    arr1[position:pos_end] = arr2
end

# %%
arr, item = setup(10_000, Float64)
arr = copy(arr)
@benchmark grow_by_push(copy(arr), item, 5_000) setup=(arr)

# %%

@benchmark append_and_return(arr1, arr2) setup=(arr1=zeros(Float64, 10_000); arr2=ones(Float64, 5_000))

# %%

@benchmark append(arr1, arr2) setup=(arr1=zeros(Float64, 10_000); arr2=ones(Float64, 5_000))

# %%

@benchmark append!(arr1, arr2) setup=(arr1=zeros(Float64, 10_000); arr2=ones(Float64, 5_000))

# %%
@benchmark zeros(Float64, 15_000)

# %%
@benchmark Vector{Float64}(undef, 15_000) # a lot faster to define the memory w/o filling it

# %%
@benchmark insertat!(arr1, arr2, 10_001) setup=(arr1 = zeros(Float64,15_000); arr2=ones(Float64, 5_000))


# %%
@benchmark insertat_loop!(arr1, arr2, 10_001) setup=(arr1 = zeros(Float64,15_000); arr2=ones(Float64, 5_000))

# %% [markdown]
# ## Deleting strategies

# %%
d = @benchmarkable deleteat!(a, delbool) setup=(a = collect(1:15_000); 
                                 delbool = rand([false, false, false, false, false, true],15_000))

run(d)

# %%
# delete by index without sorting (set indices in sorted order...   
# the application will need to sort, but this is for one benchmark

e = @benchmarkable deleteat!(a, delidx) setup=(a = collect(1:15_000); 
                                 delidx = collect(2500:2:7500))
run(e)

# %%
# delete by index with sorting--try different ways to sort
f = @benchmarkable deleteat!(a, sort!(delidx)) setup=(a = collect(1:15_000); 
                                 delidx = shuffle(collect(2500:2:7500)))
run(f)

# %%
# evaluate different ways to sort: sort and assign
g = @benchmarkable x = sort(delidx) setup=(delidx = shuffle(collect(2500:2:7500)))
run(g)

# %%
# evaluate different ways to sort: sort in place--> this is better than allocating a new vector
g = @benchmarkable sort!(delidx) setup=(delidx = shuffle(collect(2500:2:7500)))
run(g)

# %%

h = @benchmarkable delidx[sortperm(delidx)] setup=(delidx = shuffle(collect(2500:2:7500)))

run(h)

# %%
ix = Vector{Int}(undef, 2501)
i = @benchmarkable delidx[sortperm!($ix,delidx)] setup=(delidx = shuffle(collect(2500:2:7500)))

run(i)

# %% [markdown]
# ### Resizing to grow a bucket

# %%
j = @benchmarkable resize!(bucket, 15_000) setup = (bucket = rand([1,2,3],12500))
run(j)

# %%
k = @benchmarkable append!(bucket, zeros(Int,2500)) setup = (bucket = rand([1,2,3],12500))
run(k)

# %%
l = @benchmarkable append!(bucket, Vector{Int}(undef,2500)) setup = (bucket = rand([1,2,3],12500))
run(l)

# %% [markdown]
# ### zeroing and sorting the bucket

# %%
a = collect(1:15_000)
zz = collect(2500:2:7500)
a[zz] .= 0
sort!(a, rev=true)

# %%
m = @benchmarkable sort!(a, rev=true)  setup=(a = collect(1:15_000); zz = collect(2500:2:7500); a[zz].= 0)
run(m)

# %% [markdown]
# ### Instead of a bucket, change the values in place and then filter

# %%
n = @benchmarkable (findall(a .== 2)) setup=(a = fill(1, 15_000);idx = collect(2500:2:7500);a[idx] .= 2; )


run(n)

# %%
