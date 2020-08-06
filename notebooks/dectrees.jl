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
dectrees, all_decpoints = setup_dt("../parameters/dec_tree_all_25.csv");

# %%
res = CovidSim.sanity_test_all(dectrees)

# %%
age_dist

# %%
wgt_deaths = age_dist .* res[:,4]

# %%
deathpct_by_age = wgt_deaths ./ sum(wgt_deaths)

# %% [markdown]
#  ```
#  "0-20"    27  0.039
#  "20-40"  101  0.147
#  "40-60"  222  0.322
#  "60-80"  257  0.373
#  "80+"     82  0.119
#  "Total"  689  1.0
#  ```

# %%
sum(deathpct_by_age)

# %%
