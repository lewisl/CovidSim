# ---
# jupyter:
#   jupytext:
#     formats: jl:percent,ipynb
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
alldict = setup(180, [locale])

# %%
alldict["dat"]

# %%
ilmat = alldict["dat"]["popdat"][locale]

# %%
ages = alldict["dat"]["agegrp_idx"][38015]

# %%
columnnames(ilmat)

# %%
countmap(ilmat.agegrp)

# %%
countmap(ilmat.status)  # everyone begins as unexposed

# %%
geodf = alldict["geo"]   # the date for all locales has been read into a dataframe

# %%
density_factor = geodf[geodf[!, :fips] .== locale, :density_factor][]

# %%
spreadparams = alldict["sp"]  # the spread parameters are loaded as a dict of float arrays

# %%
keys(spreadparams)

# %%
typeof(spreadparams.shape)

# %%
typeof(spreadparams.riskmx)

# %%
contact_factors = spreadparams.contact_factors

# %%
typeof(contact_factors)

# %%
contact_factors[5]

# %%
touch_factors =  spreadparams.touch_factors

# %%
touch_factors[3]

# %%
limdict = CovidSim_ilm.limdict
limdict(touch_factors, <)

# %%
# is shifter working?
shifter(touch_factors, (.18, .3)...)[3]

# %%
dectree = alldict["dectree"] # the decision trees for all age groups are loaded

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
result_dict, series = run_a_sim(180, locale, showr0=false, silent=true, runcases=[seed_1_6]);

# %%
result_dict

# %%
popdat = result_dict["dat"]["popdat"][locale]

# %%
countmap(popdat.cond)

# %%
countmap(popdat.status)

# %%
virus_outcome(series, locale, base=:pop)

# %%
series[locale][:cum]

# %% [markdown]
# # Plotted results

# %%
cumplot(series, locale)

# %% [markdown]
# Note that the orangle line labeled Infectious that shows the number of infected people is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: Some people got sick today. Some people get better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. Newspaper tracking shows the new active infections of each day--who got sick today? The next day, if no one new got sick the line would be at zero--even though the people who got sick aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. This is the least common way to show the data.

# %% [markdown]
# ## Test a social distancing case

# %%
sd1 = sd_gen(startday = 45, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, include_ages=[])    

# %%
sd1_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, include_ages=[])

# %%
result_dict, series = run_a_sim(180, locale, showr0=false, silent=true, runcases=[seed_1_6, sd1, sd1_end]);

# %%
virus_outcome(series, locale, base=:pop)

# %%
cumplot(series, locale)

# %%
cumplot(series, locale,[:infectious, :dead])

# %%
outdat = result_dict["dat"]["popdat"][locale]

# %% [markdown]
# ## Social distancing only among those age40_59, age60_79, age80_plus

# %%
sdolder = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, 
    include_ages=[age40_59, age60_79, age80_up])    

# %%
sdolder_end = sd_gen(startday = 80, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age40_59, age60_79, age80_up])    

# %%
result_dict, series = run_a_sim(180, locale, showr0=false, silent=true, 
    runcases=[seed_1_6, sdolder, sdolder_end]);


# %%
olderdat = result_dict["dat"]["popdat"][locale]

sd = findall(olderdat.s_d_comply .!= :none)

@Select(agegrp, s_d_comply)(olderdat[sd])

# %%
cumplot(series, locale)

# %%
cumplot(series, locale, [:infectious, :dead])

# %% [markdown]
# ## Social Distancing starts with everyone and then the younger folks party

# %%
sdyoung_end = sd_gen(startday = 100, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age0_19, age20_39])    

# %%
result_dict, series = run_a_sim(180, locale, showr0=false, silent=true, 
    runcases=[seed_1_6, sd1, sdyoung_end]);

# %%
mixdat = result_dict["dat"]["popdat"][locale]

sd = findall(mixdat.s_d_comply .!= :none)
young = findall((mixdat.agegrp .== age0_19) .| (mixdat.agegrp .== age20_39))
sd_young_idx = intersect(sd, young)
old = findall((mixdat.agegrp .== age40_59) .| (mixdat.agegrp .== age60_79) .| (mixdat.agegrp .== age80_up));

# %%
youngtab = Table(mixdat[young])
count(youngtab.s_d_comply .== :none)

# %%
oldtab = Table(mixdat[old])
count(oldtab.s_d_comply .!= :none)

# %%
typeof(mixdat.agegrp)

# %%
cumplot(series, locale)

# %%
cumplot(series, locale, [:infectious, :dead])

# %%
@Select(status, agegrp, cond, s_d_comply)(ilmat)

# %% [markdown]
# alldict

# %%
alldict

# %%
ages = alldict["dat"]["agegrp_idx"][38015]

# %%
include_ages = [age0_19, age20_39]

# %%
union((ages[i] for i in include_ages)...)

# %%
ilmat.s_d_comply[collect(1:5:95000)] .= :test

# %%
incase_idx = findall(ilmat.s_d_comply .== :test)

# %%
byage_idx = intersect(incase_idx, union((ages[i] for i in include_ages)...))

# %%
