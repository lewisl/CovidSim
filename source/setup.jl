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
        spreadparams = build_spread_params(spfilename)

    # transition decision trees     
        dt_dict = setup_dt(dectreefilename)
        dectree = dt_dict["dt"]

    # isolation probabilities: not sure we need this
        # iso_pr = build_iso_probs()

    return Dict("dat"=>datadict, "dectree"=>dectree, "geo"=>geodata, "sp"=>spreadparams)  
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


function convert_to_enumkeys(dectree)

    for (k_age, v_age) in dectree
        for (k_sickday, v_sickday) in v_age
        end
    end
end

"""
Pre-allocate and initialize population data for one locale in the simulation.
"""
function pop_data(pop; age_dist=age_dist, cols="all")

    if cols == "all"
        parts = apportion(pop, age_dist)
        dat = Table(
            status = fill(unexposed, pop),    
            agegrp = reduce(vcat,[fill(age, parts[Int(age)]) for age in agegrps]), 
            cond = fill(notsick, pop),
            sickday = zeros(Int, pop),   
            recov_day = zeros(Int, pop),  
            dead_day = zeros(Int, pop),   
            cluster = zeros(Int, pop), 
            sdcomply = fill(:none, pop),  
            vax = zeros(Int, pop),   
            vax_day = zeros(Int, pop),  
            test = falses(pop),  
            test_day = zeros(Int, pop),  
            quar = falses(pop),
            quar_day = zeros(Int, pop))

    elseif cols == "track"
        parts = apportion(pop, age_dist)
        dat = Table(
            status = fill(unexposed, pop),        
            agegrp=reduce(vcat,[fill(age, parts[Int(age)]) for age in instances(agegrps)]), 
            cond = fill(notsick, pop),  
            sickday = zeros(Int, pop))  

    else
        @error "Wrong choice of cols in pop_data: $cols"
    end    

    return dat       
end


function hist_dict(locales, n_days; conds=allconds, agegrps=n_agegrps)
    dat = Dict{Int64, Array{Int}}()
    for loc in locales
        dat[loc] = zeros(Int, n_days, last(last(map2series))) # (conds, agegrps + 1, n_days) => (8, 6, 150)
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


function send_risk_by_recv_risk(send_risk, recv_risk)
    recv_risk' .* send_risk  # (sickdaylim, agegrps)
end


function build_spread_params(spfilename)

    spread_inputs = YAML.load_file(spfilename)

    required_params = ["send_risk", "recv_risk", "contact_factors", "touch_factors", "shape"]
    has_all = true
    missing = []
    for p in required_params
        if !haskey(spread_inputs, p)
            push!(missing, p)
            has_all = false
        end
    end
    @assert has_all "required keys: $missing not in $(spfilename)"

    send_risk = send_risk_by_recv_risk(spread_inputs["send_risk"], spread_inputs["recv_risk"])

    # named tuple doesn't result in type instability of Dict that requires "function barrier" to fix
    spreadparams = (
        send_risk          = spread_inputs["send_risk"]::Vector{Float64},
        recv_risk          = spread_inputs["recv_risk"]::Vector{Float64},
        contact_factors    = Dict(agegrp(Int(k1)) => 
                                Dict(symcond[Symbol(k2)] => Float64(v2) for (k2, v2) in v1)  for (k1, v1) in spread_inputs["contact_factors"]),
        touch_factors      = Dict(agegrp(Int(k1)) => 
                                Dict(symallconds[Symbol(k2)] => Float64(v2) for (k2, v2) in v1)  for (k1, v1) in spread_inputs["touch_factors"]),
        shape              = spread_inputs["shape"],
        # riskmx             = send_risk
        )
    
    return spreadparams
end


#####################################################################################
# dodgy math helper functions
#####################################################################################

@inline @fastmath function shifter(x::Array, oldmin, oldmax, newmin, newmax)
    newmin .+ (newmax - newmin) / (oldmax - oldmin) .* (x .- oldmin)
end

@inline @fastmath function shifter(x::Float64, oldmin, oldmax, newmin, newmax)
    newmin + (newmax - newmin) / (oldmax - oldmin) * (x - oldmin)
end

@inline function shifter(x::Array, newmin, newmax)
    oldmin = minimum(x)
    oldmax = maximum(x)
    shifter(x, oldmin, oldmax, newmin, newmax)
end


"""
    limdict(dct::Dict, op::Function)

Finds minimum or maxium value of the leaves of a dict.
Warning: not general! works on dict with 2 levels and 
numerical values at the lower level.
"""
function limdict(dct::AbstractDict, op::Function)
    minop = <
    cv = op == minop ? Inf : -Inf
    for v1 in values(dct)
        for v2 in values(v1)
            cv = op(v2, cv) ? v2 : cv
        end
    end
    return cv
end


"""
Warning: not general! works on dict with 2 levels and 
numerical values at the lower level.
"""
@inline function shifter(d::AbstractDict, newmin, newmax)
    ret = deepcopy(d)
    oldmin = limdict(d, <)
    oldmax = limdict(d, >)
    
    for k1 in keys(ret)
        for k2 in keys(ret[k1])
            x = ret[k1][k2]
            # ret[k1][k2] = newmin + (newmax - newmin) / (oldmax - oldmin) * (x - oldmin)
            ret[k1][k2] = shifter(x, oldmin, oldmax, newmin, newmax)
        end
    end

    return ret
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
    agegrp_filt_bit = Dict(age => dat.agegrp .== age for age in agegrps)
    agegrp_filt_idx = Dict(age => findall(agegrp_filt_bit[age]) for age in agegrps)
    return agegrp_filt_bit, agegrp_filt_idx
end
# agegrp_filt_bit, agegrp_filt_idx = precalc_agegrp_filt(ilmat);
