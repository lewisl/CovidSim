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

# %%
using DataFrames
using Plots
pyplot()

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# %%
alldict, env, series = run_a_sim(180,10, showr0=false, 
       dtfilename="../parameters/dec_tree_all_25.csv",
       silent=true,spreadcases=[],
       runcases=[seed_1_6]);
geo = alldict["geo"];

# %%
cumplot(series,10,geo=geo)

# %%
1312/62400

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
deadvals = series[10][:cum][end,[map2series.dead]...]
pctvals = round.([deadvals[i] / deadvals[6] for i in 1:length(deadvals)], digits=3)
deadtbl = hcat(agelabels, deadvals, pctvals)

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
recovals = series[10][:cum][end,[map2series.recovered]...]
pctvals = round.([recovals[i] / recovals[6] for i in 1:length(recovals)], digits=3)
deadtbl = hcat(agelabels, recovals, pctvals)

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
unexvals = series[10][:cum][end,[map2series.unexposed]...]
pctvals = round.([unexvals[i] / unexvals[6] for i in 1:length(unexvals)], digits=3)
deadtbl = hcat(agelabels, unexvals, pctvals)

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

# %%
cumhistmx = alldict["dat"]["cumhistmx"]
newhistmx = alldict["dat"]["newhistmx"]
openmx = alldict["dat"]["openmx"];

# %%
summary = (total_infected = series[10][:cum][1, 6] - series[10][:cum][180,6],
total_pop = series[10][:cum][180,6] + series[10][:cum][180,54],
whos_left = series[10][:cum][180,map2series.dead[6]] + series[10][:cum][180,map2series.recovered[6]]
     + series[10][:cum][180,map2series.infectious[6]] + series[10][:cum][180,map2series.unexposed[6]],
end_unexposed = series[10][:cum][180,map2series.unexposed[6]],
end_infected = series[10][:cum][180,map2series.infectious[6]],
end_recovered = series[10][:cum][180,map2series.recovered[6]],
end_dead = series[10][:cum][180,map2series.dead[6]])

# %%
transeries = DataFrame(transq)
trans = (dead = sum(transeries[:,:dead]), recovered = sum(transeries[:,:recovered]))

# %%
err = summary.total_infected - (trans.recovered + trans.dead + summary.end_infected)

# %%
spreadseries = day2df(spreadq)
check_infected = sum(spreadseries[:,:infected])

# %% [markdown]
# end_exposed is ok (off by 2 from age rounding and 6 seeds)
# total infected is ok; matches check_infected

# %% [markdown]
# ##### Some random documentation for me

# %%
plotattr()

# %%
plotattr(:Subplot)

# %% jupyter={"source_hidden": true}
plotattr(:Series)

# %%
