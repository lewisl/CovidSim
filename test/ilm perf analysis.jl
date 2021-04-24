

##
# using Traceur
using Revise
using DataStructures

## setupl
cd("ilm-src")
using CovidSim_ilm
loc = 38015
alldict = setup(180, [loc])

## data structures
dt_dict = alldict["dt_dict"]  # decision trees for transition
popdat = alldict["dat"]["popdat"]
agegrp_idx = alldict["dat"]["agegrp_idx"]
cumhistmx = alldict["dat"]["cumhistmx"]
newhistmx = alldict["dat"]["newhistmx"]
geodf = alldict["geo"]
spread_params = alldict["sp"]

env = initialize_sim_env(geodf; spread_params...)
density_factor = geodf[geodf[!, :fips] .== loc, :density_factor][]


## prep for first run
reset!(ctr, :day)  # return and reset key to 0 :day leftover from prior runs
inc!(ctr, :day)  # increment the simulation day counter
loc = 38015

## seed population with sick people
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, CovidSim_ilm.nil, CovidSim_ilm.agegrps)

##
inc!(ctr, :day)  # increment the simulation day counter
spread!(popdat, loc, [], env, density_factor)  

##
transition!(popdat, loc, dt_dict) 



##

result_dict, env, series = run_a_sim(180, 38015, showr0=false, silent=true, spreadcases=[], runcases=[seed_1_6]);


##
cumplot(series, 38015)
