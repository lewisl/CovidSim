######################################################################################
# setup and initialization functions: ILM Model
######################################################################################


function setup(n_days, locales;  # must provide following inputs
    geofilename="../data/geo2data.csv", 
    dectreefilename="../parameters/dec_tree_all_25.yml",
    spfilename="../parameters/spread_params.yml")



    # geodata
        geodata = buildgeodata(geofilename)

    # simulation data matrix
        datadict = build_data(locales, geodata, n_days)

    # spread parameters
        spread_params = read_spread_params(spfilename)

    # transition decision trees     
        dt_dict = setup_dt(dectreefilename)

    # isolation probabilities: not sure we need this
        # iso_pr = build_iso_probs()

    return Dict("dat"=>datadict, "dt_dict"=>dt_dict, "geo"=>geodata, "sp"=>spread_params)  
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

    popdat = Dict(loc => pop_data(geodata[geodata[:, "fips"] .== loc, "pop"][1]) for loc in locales)

    # precalculate agegrp indices
    agegrp_idx = Dict(loc => precalc_agegrp_filt(popdat[loc])[2] for loc in locales)

    cumhistmx = hist_dict(locales, n_days)
    newhistmx = hist_dict(locales, n_days)
    # return Dict("popdat"=>popdat, "isolatedmx"=>isolatedmx, "testmx"=>testmx, "cumhistmx"=>cumhistmx, "newhistmx"=>newhistmx)
    return Dict("popdat"=>popdat, "agegrp_idx"=>agegrp_idx, "cumhistmx"=>cumhistmx, "newhistmx"=>newhistmx)
end


"""
Pre-allocate and initialize population data for one locale in the simulation.
"""
function pop_data(pop; age_dist=age_dist, intype=T_int[], cols="all")

    if cols == "all"
        parts = apportion(pop, age_dist)
        dat = Table(
            status = fill(intype(unexposed), pop),    
            agegrp = reduce(vcat,[fill(i, parts[i]) for i in agegrps]),
            cond = zeros(intype, pop),
            lag = zeros(intype, pop),   
            recov_day = zeros(intype, pop),  
            dead_day = zeros(intype, pop),   
            cluster = zeros(intype, pop),   
            vax = falses(pop),   
            vax_day = zeros(intype, pop),  
            test = falses(pop),  
            test_day = zeros(intype, pop),  
            quar = falses(pop),
            quar_day = zeros(intype, pop))

    elseif cols == "track"
        parts = apportion(pop, age_dist)
        dat = Table(
            status = fill(intype(unexposed), pop),        
            agegrp = reduce(vcat,[fill(i, parts[i]) for i in agegrps]),
            cond = zeros(intype, pop),  
            lag = zeros(intype, pop))  

    else
        @error "Wrong choice of cols in pop_data: $cols"
    end    

    return dat       
end


function hist_dict(locales, n_days; conds=length(conditions), agegrps=n_agegrps)
    dat = Dict{Int64, Array{T_int[]}}()
    for loc in locales
        dat[loc] = zeros(T_int[], n_days, last(last(map2series))) # (conds, agegrps + 1, n_days) => (8, 6, 150)
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

    # # reshape and flip contact_factors
    #     cf = copy(spread_params["contact_factors"])
    #     spread_params["contact_factors"] = permutedims(reshape(cf,5,4), (2,1))
    # # reshape and flip touch_factors
    #     tf = copy(spread_params["touch_factors"])
    #     spread_params["touch_factors"] = permutedims(reshape(tf,5,6), (2,1))
    # # change keys to symbols--so we can use this as keyword arguments

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


function apportion(x::Int, splits::Array)
    @assert isapprox(sum(splits), 1.0)
    maxidx = argmax(splits)
    parts = round.(Int, splits .* x)
    diff = sum(parts) - x
    parts[maxidx] -= diff
    return parts
end


######################################################################################
# precalculate agegrp indices--these do not change during the simulation
######################################################################################


function precalc_agegrp_filt(dat)  # dat for a single locale
    agegrp_filt_bit = Dict(agegrp => dat.agegrp .== agegrp for agegrp in agegrps)
    agegrp_filt_idx = Dict(agegrp => findall(agegrp_filt_bit[agegrp]) for agegrp in agegrps)
    return agegrp_filt_bit, agegrp_filt_idx
end
# agegrp_filt_bit, agegrp_filt_idx = precalc_agegrp_filt(ilmat);




######################################################################################
# SimEnv: simulation environment
######################################################################################


"""
Struct for variables used by many functions = the simulation environment
    
- pre-allocate large arrays, accessed and modified frequently
- hold complex parameter sets
"""
struct SimEnv{T<:Integer}      # the members are all mutable so we can change their values
    geodata::DataFrames.DataFrame
    riskmx::Array{Float64, 2}            # laglim,5
    contact_factors::Dict{Int64, Dict{Any, Any}}   # 4,5 parameters for spread!
    touch_factors::Dict{Int64,Dict{Any, Any}}     #  6,5  parameters for spread!
    send_risk::Array{Float64, 1}  # laglim,  parameters for spread!
    recv_risk::Array{Float64,1}   # 5,  parameters for spread!

    shape::Float64                       # parameter for spread!

    # constructor with keyword arguments and type compatible fillins--not suitable as defaults, see initialize_sim_env
    # T_int[] should be one of Int64, Int32 when calling the constructor
    function SimEnv{T}(; 
            geodata=DataFrame, # geodata
            riskmx=zeros(Float64, 0,0),
            contact_factors=Dict(),
            touch_factors=Dict(),
            send_risk=zeros(Float64,laglim),
            recv_risk=zeros(Float64, 5),
            shape=1.0
        ) where T<:Integer
        return new(geodata, riskmx, contact_factors, touch_factors, send_risk, recv_risk, shape)

    end
end



function initialize_sim_env(geodata; contact_factors, touch_factors, send_risk, recv_risk, shape)

    ret = SimEnv{T_int[]}(
        geodata=geodata,
        riskmx = send_risk_by_recv_risk(send_risk, recv_risk), # zeros(Float64,laglim,5),
        contact_factors = contact_factors,
        touch_factors = touch_factors,
        send_risk = send_risk,
        recv_risk = recv_risk,
        shape = shape)

    return ret
end

    # contact_factors and touch_factors look like:
    #=
    contact_factors = [
        1    1.8    1.8     1.5     1.0;     # nil
        1    1.7    1.7     1.4     0.9;     # mild
        0.7  1.0    1.0     0.7     0.5;   # sick
        0.5  0.8    0.8     0.5     0.3]   # severe

    # agegrp    1     2      3       4       5

    touch_factors = [
        .55    .62     .58     .4    .35;    # unexposed
        .55    .62     .58     .4    .35;    # recovered
        .55    .62     .58     .4    .35;    # nil
        .55    .6      .5      .35   .28;    # mild
        .28   .35      .28     .18   .18;    # sick
        .18   .18      .18     .18   .18]    # severe
    =#
