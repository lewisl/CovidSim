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
using YAML

# %%
cd(joinpath(homedir(), "Dropbox/Covid Modeling/Covid"))

# %%
sp = YAML.load_file("parameters/spread_params.yml")

# %%
sp["touch_factors"]

# %%
sp[:send_risk]

# %%
sp_sym = YAML.load_file("parameters/sp_sym.yml", dicttype=Dict{Symbol,Any})

# %%
sp_sym[:contact_factors]

# %%
convert(Dict{Symbol, Dict{Symbol, Float64}}, sp_sym[:contact_factors])

# %%
sp_sym[:contact_factors][:age0_19]

# %%
sp_sym[:shape]

# %%
Symbol('1')

# %%
Int(Symbol("1"))

# %%
typeof(ans)

# %%
