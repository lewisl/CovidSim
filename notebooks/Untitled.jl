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
using StatsBase
using DelimitedFiles

# %%
ilmat = readdlm("../data/ilmtestdata.csv", ',',Int, header=true)[1]
ilmat = repeat(ilmat, 10_000)
r
