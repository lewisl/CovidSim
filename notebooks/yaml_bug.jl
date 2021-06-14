# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
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
using YAML
using PrettyPrint

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid/ilm-src"))

# %%
newdectree_fname="../parameters/transition.yml"

# %%
newtree = YAML.load_file(newdectree_fname)

# %%
mini2_fname = "../parameters/mini2.yml"

# %%
mini2 = YAML.load_file(mini2_fname)

# %% [raw]
#

# %%
mini1_fname = "../parameters/mini1.yml"

# %%
mini1 = YAML.load_file(mini1_fname)

# %%
pprint(mini1)

# %%
new_block_fname = "../parameters/new_yml_block.yml"

# %%
new_block = YAML.load_file(new_block_fname)

# %%
new_block["\ufeff1"]

# %%
