######################################################################################
# setup and initialization functions
######################################################################################


function setup(n_days, locales; 
    geofilename="../data/geo2data.csv", 
    dectreefilename="../parameters/dec_tree_all_25.yml",
    spfilename="../parameters/spread_params.toml")

    # geodata
        geodata = buildgeodata(geofilename)

    # simulation data matrix
        datadict = build_data(locales, geodata, n_days)
        openmx = datadict["openmx"]   # alias to make it easier to do initialization

    # spread parameters
        spread_params = read_spread_params(spfilename)

    # transition decision trees 
        dt, all_decpoints = setup_dt(dectreefilename)

    # isolation probabilities: not sure we need this
        # iso_pr = build_iso_probs()

    return Dict("dat"=>datadict, "dt"=>dt, "decpoints"=>all_decpoints,
                "geo"=>geodata, "sp"=>spread_params)  # , "iso_pr"=>iso_pr
end


"""
Convert a vector of dates from a csv file in format "mm/dd/yyyy"
to a vector of Julia numeric Date values in format yyyy-mm-dd
"""
function quickdate(strdates)  # 20x faster than the built-in date parsing, e.g.--runs in 5% the time
    ret = [parse.(Int,i) for i in split.(strdates, '/')]
    ret = [Date.(i[3], i[1], i[2]) for i in ret]
end


function build_data(locales, geodata, n_days)

    pop = [geodata[geodata[:, "fips"] .== loc, "pop"][1] for loc in locales]

    openmx = Dict(loc => pop_data(geodata[geodata[:, "fips"] .== loc, "pop"][1]) for loc in locales)
    # isolatedmx = data_dict(locales, lags=size(lags,1), conds=size(conditions,1), agegrps=size(agegrps,1))
    # testmx = data_dict(locales, lags=size(lags,1), conds=size(conditions,1), agegrps=size(agegrps,1))

    cumhistmx = hist_dict(locales, n_days)
    newhistmx = hist_dict(locales, n_days)
    # return Dict("openmx"=>openmx, "isolatedmx"=>isolatedmx, "testmx"=>testmx, "cumhistmx"=>cumhistmx, "newhistmx"=>newhistmx)
    return Dict("openmx"=>openmx, "cumhistmx"=>cumhistmx, "newhistmx"=>newhistmx)
end


"""
Pre-allocate and initialize population data for one locale in the simulation.
"""
function pop_data(pop; age_dist=age_dist, intype=Int16)

    status = fill(intype(unexposed), pop) # Array{Int,1}(undef, popsize)
    agegrp = convert.(intype,rand(Categorical(age_dist), pop))  # Array{Int,1}(undef, popsize)
    cond = zeros(intype, pop)  # Array{Int,1}(undef, popsize) 
    lag = zeros(intype, pop)  # Array{Int,1}(undef, popsize) 
    recov_day = zeros(intype, pop)  # Array{Int,1}(undef, popsize) 
    dead_day = zeros(intype, pop)  # Array{Int,1}(undef, popsize) 
    cluster = zeros(intype, pop)  # Array{Int,1}(undef, popsize) 
    vax = zeros(intype, pop)  # Array{Int,1}(undef, popsize) 
    vax_day = zeros(intype, pop)  # Array{Int,1}(undef, popsize) 
    test = zeros(intype, pop)  # Array{Int,1}(undef, popsize) 
    test_day = zeros(intype, pop)  # Array{Int,1}(undef, popsize) 
    quar = zeros(intype, pop)
    quar_day = zeros(intype, pop)

    dat = hcat(status, agegrp, cond, lag, cluster, recov_day, dead_day, vax, 
        vax_day, test, test_day, quar, quar_day)

    return dat       
end


function hist_dict(locales, n_days; conds=length(conditions), agegrps=n_agegrps)
    dat = Dict{Int64, Array{T_int[]}}()
    for loc in locales
        dat[loc] = zeros(T_int[], conds, agegrps+1, n_days) # (conds, agegrps + 1, n_days) => (8, 6, 150)
    end
    return dat       
end


function buildgeodata(filename)
    geo = DataFrame(CSV.File(filename))
    insertcols!(geo, "density_factor" => shifter(geo[:, "density"],0.9,1.25))

    # fix dates   
    insertcols!(geo, "anchor2" => quickdate(geo[:, "anchor"]))
    insertcols!(geo, "limit2" => quickdate(geo[:, "limit"]))
    select!(geo, Not(["anchor", "limit"]))
    rename!(geo, "anchor2" => "anchor")
    rename!(geo, "limit2" => "limit")

    return geo
end


function read_spread_params(spfilename)
    spread_params = YAML.load_file(spfilename)

    required_params = ["send_risk", "recv_risk", "contact_factors", "touch_factors"]
    has_all = true
    missing = []
    for p in required_params
        if !haskey(spread_params, p)
            push!(missing, p)
            has_all = false
        end
    end
    @assert has_all "required keys: $missing not in $(spfilename)"

    # reshape and flip contact_factors
        cf = copy(spread_params["contact_factors"])
        spread_params["contact_factors"] = permutedims(reshape(cf,5,4), (2,1))
    # reshape and flip touch_factors
        tf = copy(spread_params["touch_factors"])
        spread_params["touch_factors"] = permutedims(reshape(tf,5,6), (2,1))
    # change keys to symbols--so we can use this as keyword arguments
    return Dict(Symbol(k)=>v for (k,v) in spread_params)
end


# calculate density_factor in setup, per locale

function shifter(x::Array, oldmin, oldmax, newmin, newmax)
    newmin .+ (newmax - newmin) / (oldmax - oldmin) .* (x .- oldmin)
end


function shifter(x, newmin=0.9, newmax=1.5)
    oldmin = minimum(x)
    oldmax = maximum(x)
    shifter(x, oldmin, oldmax, newmin, newmax)
end


######################################################################################
# SimEnv: simulation environment
######################################################################################


"""
Struct for variables used by many functions = the simulation environment
    
- pre-allocate large arrays, accessed and modified frequently
- hold complex parameter sets
"""
struct SimEnv{T<:Integer}      # the members are all mutable so we can change their values
    geodata::Array{Any, 2}
    spreaders::Array{T, 3} # laglim,4,5
    all_accessible::Array{T, 3} # laglim,6,5
    contacts::Array{T, 3} # laglim,4,5
    simple_accessible::Array{T, 2} # 6,5
    peeps::Array{T, 2} # 6,5
    touched::Array{T, 3} # laglim,6,5
    lag_contacts::Array{T, 1} # laglim,
    riskmx::Array{Float64, 2} # laglim,5
    contact_factors::Array{Float64, 2}  # 4,5 parameters for spread!
    touch_factors::Array{Float64, 2}  #  6,5  parameters for spread!
    send_risk_by_lag::Array{Float64, 1}  # laglim,  parameters for spread!
    recv_risk_by_age::Array{Float64,1}  # 5,  parameters for spread!
    sd_compliance::Array{Float64, 2} # (6,5) social_distancing compliance unexp,recov,nil:severe by age

    # constructor with keyword arguments and type compatible fillins--not suitable as defaults, see initialize_sim_env
    # T_int[] should be one of Int64, Int32 when calling the constructor
    function SimEnv{T}(; 
                                geodata=[T(0) "" ], # geodata
                                spreaders=zeros(T, 0,0,0),   # semicolon for all keyword (named) arguments)
                                all_accessible=zeros(T, 0,0,0),
                                contacts=zeros(T, 0,0,0),
                                simple_accessible=zeros(T, 0,0),
                                peeps=zeros(T, 0,0),
                                touched=zeros(T, 0,0,0),
                                lag_contacts=zeros(T, laglim),
                                riskmx=zeros(Float64, 0,0),
                                contact_factors=zeros(Float64, 0,0),
                                touch_factors=zeros(Float64, 0,0),
                                send_risk_by_lag=zeros(Float64,laglim),
                                recv_risk_by_age=zeros(Float64, 5),
                                sd_compliance=ones(Float64, 6,5)    
                            ) where T<:Integer
        return new(geodata, spreaders, all_accessible, contacts, simple_accessible, peeps,
                   touched, lag_contacts, riskmx, contact_factors,
                   touch_factors, send_risk_by_lag, recv_risk_by_age, sd_compliance)
    end
end



function initialize_sim_env(geodata; contact_factors, touch_factors, send_risk, recv_risk)

    ret = SimEnv{T_int[]}(
                geodata=geodata,
                spreaders=zeros(T_int[], laglim, 4, agegrps),
                all_accessible=zeros(T_int[], laglim, 6, agegrps),
                contacts=zeros(T_int[], laglim, 4, agegrps),
                simple_accessible=zeros(T_int[], 6, agegrps),
                peeps=zeros(T_int[], 6, agegrps),
                touched=zeros(T_int[], laglim, 6, agegrps),
                lag_contacts=zeros(T_int[], laglim),
                riskmx = send_risk_by_recv_risk(send_risk, recv_risk), # zeros(Float64,laglim,5),
                contact_factors = contact_factors,
                touch_factors = touch_factors,
                send_risk_by_lag = send_risk,
                recv_risk_by_age = recv_risk,
                sd_compliance = zeros(6, agegrps))

    return ret
end

    # contact_factors and touch_factors look like:
    #=
        contact_factors = 
                [ 1    1.8    1.8     1.5     1.0;    # nil
                  1    1.7    1.7     1.4     0.9;    # mild
                0.7    1.0    1.0     0.7     0.5;   # sick
                0.5    0.8    0.8     0.5     0.3]  # severe

      # agegrp    1     2      3       4       5

        touch_factors = 
                [.55    .62     .58     .4    .35;    # unexposed
                 .55    .62     .58     .4    .35;    # recovered
                 .55    .62     .58     .4    .35;    # nil
                 .55    .6      .5      .35   .28;    # mild
                 .28   .35      .28     .18   .18;    # sick
                 .18   .18      .18     .18   .18]   # severe
    =#
