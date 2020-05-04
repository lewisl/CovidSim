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
str_50 = sd_gen(start=50, comply=.8, cf=(.2,1.3), tf=(.18,.45))

# %%
# working with specific locale
locale = 2

# %%
alldict, env, series = run_a_sim(180,locale, showr0=false, 
       dtfilename="../parameters/dec_tree_all_25.csv",
       silent=true,spreadcases=[str_50],
       runcases=[seed_1_6]);
geo = alldict["geo"];

# %%
cumplot(series,locale,geo=geo)

# %% [markdown]
# #### Death Percentage by Age Group

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
deadvals = series[locale][:cum][end,[map2series.dead]...]
pctvals = round.([deadvals[i] / deadvals[6] for i in 1:length(deadvals)], digits=3)
deadtbl = hcat(agelabels, deadvals, pctvals)

# %% [markdown]
# ### Recovered Percentage by Age Group

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
recovals = series[locale][:cum][end,[map2series.recovered]...]
pctvals = round.([recovals[i] / recovals[6] for i in 1:length(recovals)], digits=3)
deadtbl = hcat(agelabels, recovals, pctvals)

# %% [markdown]
# ### Unexposed Percentage by Age Group

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
unexvals = series[locale][:cum][end,[map2series.unexposed]...]
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
deadseries = series[locale][:cum][:,[map2series.dead]...]
n = size(deadseries,1)

# %%
ageserieslabels = [agelabels[1] agelabels[2] agelabels[3] agelabels[4] agelabels[5]]
areaplot(1:n, deadseries[:,1:5],labels=ageserieslabels, title="Deaths by Age Group")

# %% [markdown]
# ## Plots by Disease Condition

# %%
condseries = series[locale][:cum][:,[map2series.nil[6], map2series.mild[6], map2series.sick[6], 
            map2series.severe[6]]]
n = size(condseries,1);

# %%
condlabels = ["nil", "mild", "sick", "severe"]
day = 180
condday = series[locale][:cum][day,[map2series.nil[6], map2series.mild[6], map2series.sick[6], 
            map2series.severe[6]]]
condend = series[locale][:cum][end,[map2series.nil[6], map2series.mild[6], map2series.sick[6], 
            map2series.severe[6]]]
condpct = round.(condday ./ sum(condday), digits=2)
println("Approximate Percentage Disease Condition\n(across all ages)")
condtbl = hcat(condlabels, condday, condpct)

# %%
condserieslabels = [condlabels[1] condlabels[2] condlabels[3] condlabels[4]]
areaplot(1:n, condseries[:,:],labels=condserieslabels, 
    title="Disease Conditions Over Time\n(across all ages)")

# %%
condserieslabels = [condlabels[4]]
areaplot(1:n, condseries[:,4],labels="Severe", title="Potential Hospital Burden")
maxsevere = maximum(condseries[:, 4])
half_yscale = floor(Int, maxsevere * 0.7)
annotate!((6,half_yscale,Plots.text("Burden: $maxsevere", 10, :left)))

# %%

# %%

# %%

# %% [markdown]
# ### Check the Basic Identities

# %%
cumhistmx = alldict["dat"]["cumhistmx"]
newhistmx = alldict["dat"]["newhistmx"]
openmx = alldict["dat"]["openmx"];

# %% jupyter={"source_hidden": true}

summary = (
           total_infected = series[locale][:cum][1, 6] - series[locale][:cum][180,6],
           total_pop = series[locale][:cum][180,6] + series[locale][:cum][180,54],
           whos_left = series[locale][:cum][180,map2series.dead[6]] + series[locale][:cum][180,map2series.recovered[6]]
              + series[locale][:cum][180,map2series.infectious[6]] + series[locale][:cum][180,map2series.unexposed[6]],
           end_unexposed = series[locale][:cum][180,map2series.unexposed[6]],
           end_infected = series[locale][:cum][180,map2series.infectious[6]],
           end_recovered = series[locale][:cum][180,map2series.recovered[6]],
           end_dead = series[locale][:cum][180,map2series.dead[6]]
       )

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

# %%

# %%

# %%

# %% [markdown]
# ##### Some random documentation for me

# %%
plotattr()

# %%
plotattr(:Subplot)

# %%
plotattr(:Series)

# %%
plotattr(:Plot)

# %%
plotattr("size")

# %%
plotattr(:Axis)

# %%
