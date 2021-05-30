# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
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

# %% [markdown]
# # New approach to decision trees for transition

# %%
using CovidSim_ilm

# %%
using StatsBase
using DelimitedFiles
using Distributions
using PrettyPrint
using JSON2
using YAML

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))

# %% [markdown]
# ### Current YAML Approach as of 4/28/2021

# %%
dectreefilename="../parameters/dec_tree_all_25.yml"

# %%
dectree_dict = setup_dt(dectreefilename)

# %%
dectree = dectree_dict["dt"]

# %%
pprint(sort(dectree[5]))

# %% tags=[]
display_tree(dectree)

# %% [markdown]
# ## Experiments with YAML 

# %%
newdectree_fname="../parameters/new.yml"

# %%
newtree = YAML.load_file(newdectree_fname)

# %%
newtree[5]

# %%
newtree[5][9]

# %%
newtree[5][9][5]

# %% [markdown]
# ## Test and Re-write sanity check for transition phases

# %%
dt = sort(dectree[5])

# %% tags=[]
CovidSim_ilm.walksequence

# %%
for i = 1:5
    done = CovidSim_ilm.walksequence(dectree[i])
    probs,allpr = CovidSim_ilm.verifyprobs(done)
    println("for agegroup ", i)
    println(probs)
    println(allpr)
end

# %%
