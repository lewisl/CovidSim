# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.5.0
#     language: julia
#     name: julia-1.5.0-1.5
# ---

# %%
using CovidSim

# %%
using DelimitedFiles
using DataFrames
using Plots
using Dates
pyplot()

# %% [markdown]
# ## Comparing Simulation to Johns Hopkins Reported Data
# ### Seattle (Really, King County)

# %%
jhcases,jhfirst, jhlast = get_real_data()
jhdead, _, _ = get_real_data(series="dead")

# %%
n = jhcases.last.col - jhcases.first.col + 1; println(jhcases.first, " ", jhcases.last); println(n)

# %% [markdown]
# Get our locales.

# %%
seattle = (;fips=53033); newyork=(;fips=36061); bismarck=(;fips=38015); println(newyork)

# %%
sea = loc2df(confdat=jhcases.dat, deaddat=jhdead.dat, loc=seattle.fips)
rename!(sea, [:sea_infected, :sea_dead])
nyc = loc2df(confdat=jhcases.dat, deaddat=jhdead.dat, loc=newyork.fips)
rename!(nyc, [:nyc_infected, :nyc_dead])
bis = loc2df(confdat=jhcases.dat, deaddat=jhdead.dat, loc=bismarck.fips)
rename!(bis, [:bis_infected, :bis_dead])
tricities = hcat(sea, nyc, bis)

# %%
dt = CovidSim.setup_dt("../parameters/dec_tree_all.csv");

# %% [markdown]
# #### Run a simulation for the tricities

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# %%
alldict, env, series = run_a_sim(n, seattle.fips, showr0=true, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[],
       runcases=[seed_1_6]);

# %%
geo = alldict["geo"]

# %%
cumplot(series,seattle.fips,geo=geo)

# %%
cumplot(series,seattle.fips,[infectious, recovered, dead],geo=geo)

# %%
infection_outcome(series, seattle.fips)

# %% [markdown]
# #### Let's put some moderate social distancing in place around March 23
# This reduces R0 from around 2.0 to around 1.2. Still, some growth but much slower...

# %%
str_50 = sd_gen(start=50, comply=.75, cf=(.2,1.3), tf=(.18,.44))

# %%
alldict, env, series = run_a_sim(n, seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50],
       runcases=[seed_1_6]);

# %%
sim_r0(env, alldict["dt"], alldict["decpoints"])

# %%
cumplot(series,seattle.fips,[infectious, recovered, dead],geo=geo)

# %%
infection_outcome(series, seattle.fips)

# %% [markdown]
# #### What is the difference between active cases and reported cases?

# %%
cumplot(series,seattle.fips,[infectious, totinfected], geo=geo)

# %% [markdown]
# #### What's going on here?
# Infectious shows the *active* cases: the people on each day who are, on that day, infected. This means that people who have recovered or, sadly, have died are not included--they're not active any more.
#
# Totinfected shows the cumulative total of all the people who have ever become infected. This is what popular web dashboards and news sites show to depict growth and compare different places. Why don't they show *active* infected cases?  ...because they can't... They could subtract out the reported deaths, but that wouldn't make much difference, except in Italy or New York City. They don't know how many people have recovered because no public health jurisdiction accurately tracks this data.  When "anti-body" tests are reliably accurate and broadly available, then a few jurisdictions may begin to sample the population with anti-body tests. They will obtain an estimate of the percentage of the population that contracted and recovered from Coronavirus.
#
# Total cumulative infected is easier to track, if not accurately, and easier to compare across locations. In the early stages of the spread of infection, it focuses on what everyone wants to know: is the virus growing and how fast is it growing.  It tends to exaggerate the state of things because the line goes up steeply and never goes down--it can only become flat. 
#
# Active infected cases *will* go down and it shows the effect of so-called herd immunity: When there are many recovered and temporarily and partially immune people in the population, the virus spread slows down because there are fewer hosts available. Eventually, this will protect many other people as the spread becomes very slow and active virus in hosts dies out.

# %% [markdown]
# #### Align simulation and reported series on the day when 50 deaths were reached

# %%
sim20dead = Date("2020-01-22") + Day(findfirst(series[seattle.fips][:cum][:,map2series.dead[6]] .>= 50))
rpt20dead = Date("2020-01-22") + Day(findfirst(tricities[:, :sea_dead] .>= 50))
adjdays = Dates.value(sim20dead - rpt20dead)

# %%
println(tricities[49,:sea_dead], " ", series[seattle.fips][:cum][103,map2series.dead[6]])

# %%
alldict, env, series = run_a_sim(n+adjdays, seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50],
       runcases=[seed_1_6]);

# %% [markdown]
# #### Simulated vs. Reported Cases

# %%
plot(1:n, series[seattle.fips][:cum][1+adjdays:n+adjdays,map2series.totinfected[6]], label="Simulation",dpi=150, size=(400,260), tickfontsize=6)
plot!(1:n, tricities[:,:sea_infected], label="Reported Cases",dpi=150, size=(400,260), tickfontsize=6)
title!("King County Simulation vs. Reported Cases", titlefontsize=8)
xlabel!("Days: Feb. 1 to April 30", guidefontsize=8)


# %% [markdown]
# #### Which do you believe?
# Is this plausible?  Let's think about it.  The simulation with social distancing results in about 91,000 cumulative total infected and 772 deaths for a death rate of 8 tenths of one percent. The reported results show only 6207 cumulative total infected and 447 deaths for a death rate of over 7%. 
#
# Could we have unreported deaths in Seattle?  Even in the place with possibly the best and most honest approach to public health in the US, sure... Is it reasonable that Seattle's death rate is comparable to Italy's? Probably not. Let's split the difference on deaths and call it 600. Let's assume a moderately high death rate of 2%. This would mean Seattle really has 30,000 cases--so Seattle is underreporting by a factor of 5 and the simulation would be showing 3X the number cases rather than around 12X.
#
# Even so, doesn't the simulation still seem too high?  The simulated infection rate is only at 4% of King County's population. That doesn't seem outrageous--in fact, it is probably still too low. The good news is that the real death rate is probably significantly lower than Seattle's reported death rate and King County is building towards "herd immunity" more rapidly than the reported cases would suggest.
#
# This phenomenon is happening around the world with a respected article based estimating cases based on more accurate reported deaths. The article suggests that under-reporting in countries varies from 3X to more than 12X. 

# %%
plot(1:n, series[seattle.fips][:cum][1+adjdays:(n+adjdays),map2series.dead[6]], label="Simulation",dpi=150, size=(400,260), tickfontsize=6)
plot!(1:n, vcat(zeros(0),tricities[1:n-0, :sea_dead]), label="Reported Deaths",dpi=150, size=(400,260), tickfontsize=6)
title!("King County Simulation vs Reported Deaths", titlefontsize=10)
xlabel!("Days: Jan. 22 to May 23", guidefontsize=10)

# %% [markdown]
# ### Virtual New York City

# %%
alldict, env, series = run_a_sim(122, newyork.fips, showr0=true, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[],
       runcases=[seed_1_6]);

# %%
cumplot(series,newyork.fips,geo=geo)

# %% [markdown]
# Why does New York grow much faster than Seattle? Because it has a much higher population density.  The simulation uses population density as an input that influences the spread of the virus.

# %%
infection_outcome(series, newyork.fips)

# %% [markdown]
# #### Let's do some social distancing

# %%
alldict, env, series = run_a_sim(n+adjdays, newyork.fips, showr0=true, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50],
       runcases=[seed_1_6]);

# %% [markdown]
# #### Align simulation and reported series on the day when 50 deaths were reached

# %%
sim20dead = Date("2020-01-22") + Day(findfirst(series[newyork.fips][:cum][:,map2series.dead[6]] .>= 50))
rpt20dead = Date("2020-01-22") + Day(findfirst(tricities[:, :nyc_dead] .>= 50))
println("sim ", sim20dead," rept ", rpt20dead)
adjdays = Dates.value(sim20dead - rpt20dead)

# %% [markdown]
# ### With Social Distancing

# %%
cumplot(series,newyork.fips,[infectious, recovered, dead], geo=geo)

# %% [markdown]
# ### Active Cases vs. Cumulative Total Cases

# %%
cumplot(series,newyork.fips,[infectious, totinfected], geo=geo)

# %% [markdown]
# ### Simulation vs. Reported Total Cases

# %%
# plot(1:122, vcat(zeros(abs(adjdays)),series[newyork.fips][:cum][1:n+adjdays,map2series.totinfected[6]]), label="Simulation",dpi=150, size=(400,260), tickfontsize=6)
plot(1:n, series[newyork.fips][:cum][1+adjdays:n+adjdays,map2series.totinfected[6]], label="Simulation",dpi=150, size=(400,260), tickfontsize=6)
plot!(1:n, tricities[:,:nyc_infected], label="Reported Cases",dpi=150, size=(400,260), tickfontsize=7)
title!("New York City Simulation vs. Reported Cases", titlefontsize=7)
xlabel!("Days: Feb. 1 to April 30", guidefontsize=8)


# %%
plot(1:n, vcat(zeros(abs(adjdays)),series[newyork.fips][:cum][1:n+adjdays,map2series.dead[6]]), label="Simulation",dpi=150, size=(400,260), tickfontsize=6)
plot!(1:n, tricities[:,:nyc_dead], label="Reported Deaths",dpi=150, size=(400,260), tickfontsize=6)
title!("NYC Simulation vs Reported Deaths", titlefontsize=10)
xlabel!("Days: Feb. 1 to April 30", guidefontsize=10)

# %%
