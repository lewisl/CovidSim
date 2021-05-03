# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
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
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid/ilm-src"))

# %%
ilmat = readdlm("../data/ilmtestdata.csv", ',',Int, header=true)[1]
ilmat = repeat(ilmat, 10_000)
refresh = copy(ilmat)

# %% [markdown]
# ### Current YAML Approach as of 4/28/2021

# %%
dectreefilename="../parameters/dec_tree_all_25.yml"

# %%
dectree_dict = setup_dt(dectreefilename)

# %%
dectree = dectree_dict["dt"]

# %% jupyter={"outputs_hidden": true} tags=[]
display_tree(dectree)

# %%
YAML.write(dectree)

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

# %%
retree = YAML.load_file(newdectree_fname)

# %%
