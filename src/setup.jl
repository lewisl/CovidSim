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


