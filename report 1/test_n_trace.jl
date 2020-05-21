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
cd("/Users/lewis/Dropbox/Online Coursework/Covid/src"); 

# %%
using CovidSim

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)
alldict, env, series = run_a_sim(180,38015, showr0=false, silent=true,spreadcases=[],
       runcases=[seed_1_6]);
geo = alldict["geo"];

# %% [markdown]
# # No social distancing; no isolation; all mingling...

# %%
cumplot(series,38015,geo=geo)

# %%
str_50 = sd_gen(start=40, comply=.9, cf=(.2,1.1), tf=(.18,.4))
alldict, env, series = run_a_sim(180,38015, showr0=false, silent=true,
    spreadcases=[str_50],
    runcases=[seed_1_6]);

# %% [markdown]
# # Moderately Strong Social Distancing starts on Day 50

# %%
cumplot(series,38015,geo=geo)

# %%
open = sd_gen(start=70, comply=0.9, cf=(.3,1.5), tf=(.25,.50))

# %%
alldict, env, series = run_a_sim(180,38015, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6]);

# %% [markdown]
# # It's OVER--Let's Open Up.   (Uh-Oh, the Dread "Double Bump")

# %%
cumplot(series,38015,geo=geo)

# %%
t_n_t100_160=CovidSim.t_n_t_case_gen(100,160,tc_perday=1000,breakout_pct=0.20, q_comply=0.9,past_contacts=false)

# %%
alldict, env, series = run_a_sim(180,38015, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6, t_n_t100_160]);

# %% [markdown]
# # Wait!  Test and Trace and Isolate Is Supposed to Work!

# %%
cumplot(series,38015,geo=geo)

# %%
t_n_t100_160=CovidSim.t_n_t_case_gen(100,160,tc_perday=10000,breakout_pct=0.20, q_comply=0.9,past_contacts=false)

# %%
alldict, env, series = run_a_sim(180,38015, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6, t_n_t100_160]);

# %%
cumplot(series,38015,geo=geo)

# %% [markdown]
# # If only we had more test kits!

# %%
