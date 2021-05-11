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

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid/ilm-src"))

# %% [markdown]
# # Test setup and population matrix

# %%
# set locale
locale = 38015

# %%
alldict = setup(150, [locale])

# %%
alldict["dat"]

# %%
ilmat = alldict["dat"]["popdat"][locale]

# %%
columnnames(ilmat)

# %%
countmap(ilmat.agegrp)

# %%
sum(ilmat.status)  # everyone begins as unexposed

# %%
geodf = alldict["geo"]   # the date for all locales has been read into a dataframe

# %%
density_factor = geodf[geodf[!, :fips] .== locale, :density_factor][]

# %%
alldict["sp"]  # the spread parameters are loaded as a dict of float arrays

# %%
send_risk = alldict["sp"][:send_risk]

# %%
recv_risk = alldict["sp"][:recv_risk]

# %%
contact_factors = alldict["sp"][:contact_factors]

# %%
contact_factors[5]

# %%
touch_factors =  alldict["sp"][:touch_factors]

# %%
touch_factors[1]

# %%
alldict["sp"][:touch_factors]

# %%
dectree = alldict["dt_dict"]["dt"] # the decision trees for all age groups are loaded

# %%
typeof(dectree)

# %% [markdown]
# Dict{Int64, OrderedCollections.OrderedDict{Int, Dict{String, Vector{T} where T}

# %%
dectree[5][25]

# %%
typeof(dectree[5][25][7]["outcomes"])

# %%
function get_node(dectree, agegrp, sickday, fromcond)
    dectree[agegrp][sickday][fromcond]
end

# %%
@time node = get_node(dectree, 5, 25, 7)

# %%
@btime get_node(dectree, 5, 25, 7)["probs"]

# %% [markdown]
# # Create a seed case

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# %% [markdown]
# # Run a simulation

# %%
result_dict, series = run_a_sim(180, locale, showr0=false, silent=true, spreadcases=[], runcases=[seed_1_6]);

# %%
result_dict

# %%
popdat = result_dict["dat"]["popdat"][locale]

# %%
countmap(popdat.cond)

# %%
countmap(popdat.status)

# %%
virus_outcome(series, locale)

# %%
Int(5.0)

# %% [markdown]
# # Plotted results

# %%
cumplot(series, locale)

# %% [markdown]
# Note that the orangle line labeled Infectious that shows the number of infected people is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: Some people got sick today. Some people get better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. Newspaper tracking shows the new active infections of each day--who got sick today? The next day, if no one new got sick the line would be at zero--even though the people who got sick aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. This is the least common way to show the data.

# %%
818/95626

# %% [markdown]
# ### Ways to reduce allocations for indexing

# %%
a = rand(1000); b = rand(["this", "is", "it"], 1000); c = rand(1:7,1000); i = 1:1000;
t = Table(a=a, b=b, c=c, i=i)

# %%
function optfindall(p, X, maxlen=0)
    if maxlen==0
        out = Vector{Int}(undef, length(X))
    elseif isa(maxlen, Int)
        out = Vector{Int}(undef, maxlen)
    else
        out = Vector{Int}(undef, floor(Int, maxlen * length(X)))
    end
    ind = 0
    @inbounds for (i, x) in pairs(X)
        if p(x)
            out[ind+=1] = i
        end
    end
    resize!(out, ind)
    return out
end

# %%
floor(Int,.23  * length(t))

# %% tags=[]
@time optfindall(==(6),t.c,.2);

# %% jupyter={"outputs_hidden": true} tags=[]
@time findall(t.c .== 6)

# %% jupyter={"outputs_hidden": true} tags=[]
for k in eachindex(t.c .== 3)
    println(k, " ", t.i[k])
end
