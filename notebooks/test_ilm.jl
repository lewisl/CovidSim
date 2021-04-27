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
virus_outcome(series, 38015)

# %% [markdown]
# # Plotted results

# %%
cumplot(series, 38015)

# %% [markdown]
# Note that the orangle line labeled Infectious that shows the number of infected people is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: Some people got sick today. Some people get better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. Newspaper tracking shows the new active infections of each day--who got sick today? The next day, if no one new got sick the line would be at zero--even though the people who got sick aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. This is the least common way to show the data.