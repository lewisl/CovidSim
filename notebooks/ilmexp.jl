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

# %% [markdown]
# # Prototype Spread and Transition with ILM Approach

# %%
using CovidSim

# %%
using StatsBase
using DelimitedFiles
using Distributions

# %%
ilmat = readdlm("../data/ilmtestdata.csv", ',',Int, header=true)[1]
ilmat = repeat(ilmat, 10_000)
refresh = copy(ilmat)

# %% [markdown]
# ## Spread

# %%
spfilename="../parameters/spread_params.toml"
spread_params = CovidSim.read_spread_params(spfilename)
contact_factors = spread_params[:contact_factors]
touch_factors = spread_params[:touch_factors]
send_risk = spread_params[:send_risk]
recv_risk = spread_params[:recv_risk]
riskmx = CovidSim.send_risk_by_recv_risk(send_risk, recv_risk) # (lags, agegrp);

# %%
spreadidx = ilmat[:,1] .== 2
count(spreadidx)

# %%
# for this test case lets have fewer people who are infected
resetidx = sample(findall(spreadidx),72000, replace=false)
ilmat[resetidx,1] .= 1
# index of all spreaders
spreadidx = findall(ilmat[:, cpop_status] .== 2)
n_spreaders = size(spreadidx, 1);

# %%
# no of contacts for each spreader
@show n_spreaders

contacts = zeros(Int, size(spreadidx,1),2) # second column for lag of the spreader
for i in 1:size(contacts, 1)  # cond, agegrp
    scale = contact_factors[ilmat[spreadidx[i], cpop_cond]-4, ilmat[spreadidx[i], cpop_agegrp]]
    contacts[i, 1] = round(Int,rand(Gamma(1.0,scale))) # assume density_factor = 1.0
    contacts[i, 2] = ilmat[spreadidx[i], cpop_lag]   # lag of the spreader who made this contact
end
n_contacts = sum(contacts[:,1])
@show n_contacts

# assign the contacts 
contactidx = findall(ilmat[:, cpop_status] .!= dead)
n_contactable = size(contactidx, 1)
@show n_contactable
choose_contacts = sample(contactidx, min(n_contacts, n_contactable), replace=false)

# determine contacts are consequential touches and if newly infected
n_touched = 0
n_newly_infected = 0

for i in 1:size(contacts,1)
    person = choose_contacts[i]
    status = ilmat[person, cpop_status]
    characteristic =  status in [1,3] ? [1,0,2][status] : max(0,ilmat[person, cpop_cond]-2)
    agegrp = ilmat[person, cpop_agegrp]
    touched = rand(Binomial(1, touch_factors[characteristic, agegrp]))
    # println(status, " ", characteristic, " ", agegrp, " ", touched)
    n_touched += touched
    if touched == 1 && characteristic == 1
        prob = riskmx[contacts[i, 2], agegrp]
        newly_infected = rand(Binomial(1, prob))
        # println(prob, " ", newly_infected)
        if newly_infected == 1
            ilmat[person, cpop_cond] = nil # nil === asymptomatic or pre-symptomatic
            ilmat[person, cpop_status] = infectious
        end
        n_newly_infected += newly_infected
    end
end


@show n_touched
@show n_newly_infected

# %%
spreadidx = findall(ilmat[:, cpop_status] .== 2)
size(spreadidx,1)

# %% [markdown]
# ## Transition

# %%
