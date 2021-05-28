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

# %%
using BenchmarkTools
using Distributions

# %%
rs = rand(50)

# %%
rand()

# %%
catdist = Dict(:out=>[1,5,7], :pr=>[.1, .7, .2])

# %%
function catsim(catdist, rands=[])
    if isempty(rands)
        return rand()
            
            

# %%
prs = catdist[:pr]

# %%
@btime append!([$prs[1]],[$prs[i] + $prs[i-1] for i = 2:length($prs)])

# %%
function categorical_sim(prs::Vector{Float64}, do_assert=true)
    do_assert && @assert isapprox(sum(prs), 1.0)
    x = rand()
    cumpr = 0.0
    i = 0
    for pr in prs
        cumpr += pr
        i += 1
        if x <= cumpr 
            break
        end
    end
    i
end

function categorical_sim(prs::Vector{Float64}, n::Int, do_assert=true)
    do_assert && @assert isapprox(sum(prs), 1.0)
    ret = Vector{Int}(undef, n)
    
    @inbounds for i in 1:n
        ret[i] = categorical_sim(prs, false)
    end
    ret
end
        
        
    

# %%
@btime categorical_sim($catdist[:pr], false);

# %%
categorical_sim(catdist[:pr])

# %%
@btime rand(Categorical(catdist[:pr]));

# %%
24/207

# %%
typeof(Int[1,2,3])

# %%
for l in eachindex([10,20])
    println(l)
end

# %%
@btime rand(Categorical($catdist[:pr]),2000);

# %%
@btime categorical_sim($catdist[:pr],2000);

# %%
@btime zeros(Int, 2000);

# %%
@btime Vector{Int}(undef, 2000);

# %%
