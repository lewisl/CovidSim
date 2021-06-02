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

# %%
using CovidSim_ilm

# %%
using StatsBase
using TypedTables
using BenchmarkTools
using Distributions
using YAML

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))

# %% [markdown]
# ## Load data and params

# %% tags=[]
# set locale
locale = 38015

# %% tags=[]
alldict = setup(180, [locale])

# %%
spreadparams = alldict["sp"]  # the spread parameters are loaded as a dict of float arrays

# %%
keys(spreadparams)

# %%
spfilename="../parameters/spread_params.yml"

# %%
spread_inputs = YAML.load_file(spfilename)

# %%
spr_day = 6; recv_age = 4

# %%
recvrisk = spread_inputs["recv_risk"]

# %%
sendrisk = spread_inputs["send_risk"]

# %%
@btime $sendrisk[$spr_day] * $recvrisk[$recv_age]

# %%
Statcond = Union{status, condition}

# %%
foo = Array{Statcond, 1}

# %%
foo = [nil, dead]

# %%
