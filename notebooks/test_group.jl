# %% [markdown]
# ## Quick Test of Group Model ###

# %%
using Markdown
using InteractiveUtils

# %%
cd(joinpath(homedir(),"Dropbox/Online Coursework/Covid/ilm-src"))

# %%
using CovidSim_group

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# %%
result_dict, env, series = run_a_sim(180, 38015, showr0=false, silent=true, spreadcases=[], runcases=[seed_1_6]);

# %%
cumplot(series, 38015)

# %%
keys(result_dict)

# %%
result_dict["dat"]["openmx"][38015]

# %%
