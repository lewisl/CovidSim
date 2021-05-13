# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# %%
using CovidSim_ilm

# %%
using StatsBase
using TypedTables
using BenchmarkTools
using Distributions
using YAML

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid/ilm-src"))


# %%
struct Spreadcase
    day::Int
    cfdelta::Tuple{Float64,Float64}  # (4,5)  # contact factors
    tfdelta::Tuple{Float64,Float64}  # (6,5)  # touch factors
    comply::Float64             # compliance percentage
    cfcase::Dict{Int64, Dict{String, Float64}}
    tfcase::Dict{Int64, Dict{String, Float64}}
end
# %%
testdict = Dict(1 => Dict("foo" => 1e6))
testcase = Spreadcase(1, (.2, .4), (.18,.6), .9, testdict, testdict)

# %%
cases = Dict(Symbol("case_", i) => testcase for i in 1:10)

# %%
Base.summarysize(cases)

# %%
@btime $cases[:case_1];

# %%
@btime $cases[:case_1].tfcase;

# %%
case_1 = cases[:case_1]

# %%
@btime $case_1.tfcase;

# %%
cfstr = "contact_factors:                                        
  1:                                                   
    {nil: 1.1, mild: 1.1, sick: 0.7, severe: 0.5}       
  2:
    {nil: 2.1, mild: 2.0, sick: 1.0, severe: 0.6}
  3:
    {nil: 2.1, mild: 2.0, sick: 1.0, severe: 0.6}
  4:
    {nil: 1.7, mild: 1.6, sick: 0.7, severe: 0.5}
  5:
    {nil: 1.0, mild: 0.9, sick: 0.6, severe: 0.5}"

# %%
cf = YAML.load(cfstr)

# %%
cf = cf["contact_factors"]

# %%
function shifter(x::Array, oldmin, oldmax, newmin, newmax)
    newmin .+ (newmax - newmin) / (oldmax - oldmin) .* (x .- oldmin)
end


function shifter(x::Array, newmin=0.9, newmax=1.5)
    oldmin = minimum(x)
    oldmax = maximum(x)
    shifter(x, oldmin, oldmax, newmin, newmax)
end

function shifter!(d::Dict, newmin=0.9, newmax=1.5)
    oldmin = limdict(d, <)
    oldmax = limdict(d, >)
    
    for k1 in keys(d)
        for k2 in keys(d[k1])
            x = d[k1][k2]
            d[k1][k2] = newmin + (newmax - newmin) / (oldmax - oldmin) * (x - oldmin)
        end
    end
end
            

# %%
shifter([ 1 2; 3 4])

# %%
function limdict(dct, op)
    minop = <
    cv = op == minop ? Inf : -Inf
    for v1 in values(dct)
        for v2 in values(v1)
            cv = op(v2, cv) ? v2 : cv
        end
    end
    return cv
end

# %%
limdict(cf,<)

# %%
cfcopy = deepcopy(cf)

# %%
shifter!(cfcopy, 0.2, 2.0)

# %%
cfcopy

# %%
function add1(dct)
    for k in keys(dct)
        dct[k] += 1
    end
end

# %%
dct = Dict("one"=>1, "two"=>2)

# %%
add1(dct)

# %%
dct

# %%
typeof(>)

# %%
supertype(typeof(add1))

# %%
?Base.@kwdef

# %%
