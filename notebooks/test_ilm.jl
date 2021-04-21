# ---
# jupyter:
#   jupytext:
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
using CovidSim_ilm

# %%
using StatsBase
using TypedTables

# %%
cd("/Users/lewislevin/Dropbox/Online Coursework/Covid/ilm-src")

# %% [markdown]
# # Test setup and population matrix

# %%
alldict = setup(150, [38015])

# %%
alldict["dat"]

# %%
alldict["dat"]["popdat"]

# %%
Base.summarysize(alldict["dat"]["popdat"][38015])

# %%
ilmat = alldict["dat"]["popdat"][38015]

# %%
columnnames(ilmat)

# %%
countmap(ilmat.agegrp)

# %%
sum(ilmat.status)  # everyone begins as unexposed

# %%
geodf = alldict["geo"]   # the date for all locales has been read into a dataframe

# %%
density_factor = geodf[geodf[!, :fips] .== 38015, :density_factor][]

# %%
alldict["sp"]  # the spread parameters are loaded as a dict of float arrays

# %%
alldict["dt_dict"]["dt"] # the decision trees for all age groups are loaded

# %%
alldict["dt_dict"]["decpoints"]  # the decpoints for all agegrps are loaded as array of day values

# %% [markdown]
# # Create a seed case

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# %% [markdown]
# # Run a simulation

# %%
result_dict, env, series = run_a_sim(180, 38015, showr0=false, silent=true, spreadcases=[], runcases=[seed_1_6]);

# %%
result_dict

# %%
virus_outcome(series, 38015)  # has errors

# %%
cumplot(series, 38015)

# %%
67327/95626

# %%
815/95626

# %%
