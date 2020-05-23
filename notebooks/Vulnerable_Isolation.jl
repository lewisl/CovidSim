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
seattle = (fips=53033, id=2); nyc=(fips=36061, id=3); bismarck=(fips=38015,id=10)

# %%
geo = CovidSim.readgeodata("../data/geo2data.csv")
geo[:,1:7]

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# %%
alldict, env, series = run_a_sim(180,seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[],
       runcases=[seed_1_6]);

# %%
cumplot(series,seattle.fips,geo=geo)

# %% [markdown]
# ### Strong Social Distancing

# %% [markdown]
# Reset the model to defaults.

# %%
alldict, env, series = run_a_sim(180,seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[],
       runcases=[seed_1_6]);

# %%
str_50 = sd_gen(start=50, comply=.8, cf=(.2,1.3), tf=(.18,.45))

# %%
alldict, env, series = run_a_sim(180,seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50],
       runcases=[seed_1_6]);

# %%
cumplot(series,seattle.fips,[recovered, infectious, dead],geo=geo)

# %% [markdown]
# ### Open Totally (which won't happen...)
# This uses opening back to essentially no social distancing and an R0 between 1.9 and 2.0. People will voluntarily be more prudent and government recommendations and policies will provide for more limited opening. So, this shows why complete opening is't possible:  the full force of the virus does return with only a slight delay.

# %%
# Reset to defaults
alldict, env, series = run_a_sim(180,seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[],
       runcases=[seed_1_6]);

# %%
open = sd_gen(start=105, comply=0.0, cf=(.3,1.8), tf=(.18,.62))

# %%
alldict, env, series = run_a_sim(180,seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50, open],
       runcases=[seed_1_6]);

# %%
cumplot(series,seattle.fips,[infectious, dead],geo=geo)

# %%
r0_sim(;sa_pct=[1.0,0.0,0.0], density_factor=1.0, dt=[], cf=[], tf=[],
                compliance=[1.0], shift_contact=(.2,1.8), shift_touch=(.18,.62), pri=false, env=env)

# %% [markdown]
# The preceding estimates R0 based on equal representation in all demographic groups of the simulation.  But, the groups are not equally respresented so this is a slight underestimate of the socially determined R0.

# %% [markdown]
# ### Limited Opening

# %%
# reset the model to defaults
alldict, env, series = run_a_sim(180,seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[],
       runcases=[seed_1_6]);

# %%
mod_105 = sd_gen(start=105,cf=(.2,1.45), tf=(.18,.5),comply=.7)

# %%
r0_sim(;sa_pct=[1.0,0.0,0.0], density_factor=1.0, dt=[], cf=[], tf=[],
                compliance=[.75], shift_contact=(.2,1.45), shift_touch=(.18,.5), pri=false)

# %%
alldict, env, series = run_a_sim(180,seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50, mod_105],
       runcases=[seed_1_6]);

# %%
cumplot(series,seattle.fips,[recovered,infectious, dead],geo=geo)

# %% [markdown]
# Even with meaningful easing or restrictions, this moderate opening is very effective because it comes on the base of a lower curve that resulted from earlier strict social distancing.

# %% [markdown]
# This alternative would *not* be recommended. Widespread testing of asymptomatic people with contact tracing is very hard to achieve, but is a preferred alternative.  (Modeling effort is more significant: it is coming.)

# %% [markdown]
# ### An Alternative: Fewer Restrictions with Isolation of the Vulnerable
# Note that this cases models vulnerable people as the age groups from 60 to 80 and over 80. Other people are vulnerable: people with diabetes, hypertension, immuno-compromise, cancer patients, smokers, and others across age groups. It's beyond this model to attempt to represent these vulnerabilities across age groups at this poing.

# %%
mod_less_105 = sd_gen(start=105,cf=(.2,1.55), tf=(.18,.52),comply=.55)

# %%
r0_sim(;sa_pct=[1.0,0.0,0.0], density_factor=1.0, dt=[], cf=[], tf=[],
                compliance=[.65], shift_contact=(.2,1.55), shift_touch=(.18,.55), pri=false)

# %%
function isolate_vulnerable(locale; opendat=openmx, isodat=isolatedmx,testdat=openmx, env=env)
    if ctr[:day] == 105
        isolate!(.70,[unexposed, nil,mild,sick, severe],[5],1:laglim, locale; opendat=opendat, isodat=isodat)
        isolate!(.50,[unexposed,nil,mild,sick, severe],[4],1:laglim, locale; opendat=opendat, isodat=isodat)
    end
end

# %%
# reset the model to defaults
alldict, env, series = run_a_sim(180,seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[],
       runcases=[]);
laglim=25

# %%
alldict, env, series = run_a_sim(180,seattle.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50, mod_less_105],
       runcases=[seed_1_6, isolate_vulnerable]);

# %%
cumplot(series,seattle.fips,[recovered,infectious, dead],geo=geo)

# %% [markdown]
# **Assessment:**
# This was put together quickly and requires more careful review. The restrictions are considerably looser than the limited opening above, which is clear from the high infection rate. But, the number of deaths is lower. Since over 75% of deaths occur in agegroups over 60 combined with over 80, the moderation of deaths will primarily benefit those age groups. This warrants more work.

# %%
