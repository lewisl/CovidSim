# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.4.0
#     language: julia
#     name: julia-1.4
# ---

# %%
using CovidSim

# %% jupyter={"outputs_hidden": true}
using Plots
pyplot()

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# %%
alldict, env, series = run_a_sim(180,10, showr0=false, 
       silent=true,spreadcases=[],
       runcases=[seed_1_6]);
geo = alldict["geo"];

# %%
cumplot(series,10,geo=geo)

# %%
761/62238

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
deadvals = series[10][:cum][end,[map2series.dead]...]
pctvals = round.([deadvals[i] / deadvals[6] for i in 1:length(deadvals)], digits=3)
deadtbl = hcat(agelabels, deadvals, pctvals)

# %% [markdown]
# ### Worldometers Death Demographics for New York City
#
# <img src=attachment:image.png width="500" height="500">

# %% [markdown]
# #### CDC Age Demographics for Covid-19 Deaths
# <img src=attachment:image.png width="200" height="1000">

# %%
deadseries = series[10][:cum][:,[map2series.dead]...]
n = size(deadseries,1)

# %%
ageserieslabels = [agelabels[1] agelabels[2] agelabels[3] agelabels[4] agelabels[5]]
areaplot(1:n, deadseries[:,1:5],labels=ageserieslabels, title="Deaths by Age Group")

# %% [markdown]
# ##### Some random documentation for me

# %%
plotattr()

# %%
plotattr(:Subplot)

# %% jupyter={"source_hidden": true}
plotattr(:Series)
