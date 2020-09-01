# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.5.1
#     language: julia
#     name: julia-1.5
# ---

# %%
using CovidSim_ilm

# %%
using StatsBase

# %%
cd("ilm-src")

# %% [markdown]
# # Test setup

# %%
alldict = setup(150, [38015])

# %%
ilmat = alldict["dat"]["openmx"][38015]

# %%
countmap(ilmat[:,2])

# %%
sum(ilmat[:,cpop_status])  # everyone begins as unexposed

# %%
geodf = alldict["geo"]   # the date for all locales has been read into a dataframe

# %%
density_factor = geodf[geodf[:fips] .== 38015, :density_factor][]

# %%
alldict["sp"]  # the spread parameters are loaded as a dict of float arrays

# %%
alldict["dt"] # the decision trees for all age groups are loaded

# %%
alldict["decpoints"]  # the decpoints for all agegrps are loaded as array of day values

# %% [markdown]
# # Create a seed case

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# %% [markdown]
# # Run a sim without any summary tracking

# %%
alldict = run_a_sim(180, 38015, showr0=false, silent=true, spreadcases=[], runcases=[seed_1_6]);

# %%
