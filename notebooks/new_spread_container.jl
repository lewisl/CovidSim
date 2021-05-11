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
using BenchmarkTools

# %%
struct SpreadParams
    send_risk::Vector{Float64}
    recv_risk::Vector{Float64}
    shape::Float64
    touch_factors::Dict{Int64, NamedTuple{(:sick, :unexposed, :nil, :severe, :recovered, :mild), NTuple{6, Float64}}}
    contact_factors::Dict{Int64, NamedTuple{(:sick, :nil, :severe, :mild), NTuple{4, Float64}}}
    riskmx::Array{Float64, 2}
end

# %%
cd(joinpath(homedir(), "Dropbox/Covid Modeling/Covid/parameters"))

# %%
sps = YAML.load_file("spread_params.yml")

# %%
function send_risk_by_recv_risk(send_risk, recv_risk)
    recv_risk' .* send_risk  # (sickdaylim, agegrps)
end

# %% [markdown]
# ## function to build spread params as a struct
# The dicts for touch_factors and contact_factors are maps from agegrp to a named tuple for each.
# This eliminates a function barrier to breakup the dictionary, establish concrete types (no {Any}) 
# for Julia, which significantly improves performance. If we get rid of the dictionary, we can get rid 
# of the function barrier (which is a hideous thing). But, there are other complexities to this
# struct and the necessary named tuples.
# %%
function build_spread_params(spfilename)

    spread_params = YAML.load_file(spfilename)

    required_params = ["send_risk", "recv_risk", "contact_factors", "touch_factors", "shape"]
    has_all = true
    missing = []
    for p in required_params
        if !haskey(spread_params, p)
            push!(missing, p)
            has_all = false
        end
    end
    @assert has_all "required keys: $missing not in $(spfilename)"

    send_risk = send_risk_by_recv_risk(spread_params["send_risk"], spread_params["recv_risk"])

   spread_struct = SpreadParams(
       spread_params["send_risk"],
       spread_params["recv_risk"],
       spread_params["shape"],
       Dict(i => (; (Symbol(k) => v for (k,v) in spread_params["touch_factors"][i])...) 
            for i in keys(spread_params["touch_factors"])),
        Dict(i => (; (Symbol(k) => v for (k,v) in spread_params["contact_factors"][i])...) 
            for i in keys(spread_params["contact_factors"])),
       send_risk
   )
    
    return spread_struct
end

# %%
sps = build_spread_params("spread_params.yml")

# %%
sps.contact_factors

# %%
