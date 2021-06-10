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
using JSON
using YAML

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))

# %% [markdown]
# ## Use an input format closer to the intended output format

# %%
newdtfname = "../parameters/new.yml"

# %%
newdt = YAML.load_file(newdtfname)

# %%
newdt[1]

# %%
newdt[1][5]

# %%
trees = setup_dt(newdtfname)

# %% tags=[]
display_tree(trees)

# %%
trees[age80_up]

# %%
trees[age80_up][5]

# %%
trees[age80_up][5][nil]

# %% [markdown]
# #### does it check out?

# %% tags=[]
seqs = CovidSim_ilm.getseqs(trees[age80_up])

# %%
CovidSim_ilm.verifyprobs(seqs)

# %% [markdown]
# ## Test and Re-write sanity check for transition phases

# %%
CovidSim_ilm.sanitycheck(trees)

# %%
