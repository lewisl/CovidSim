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

    # isolation probabilities: not sure we need this
        # iso_pr = build_iso_probs()

    return Dict("dat"=>datadict, "dt_dict"=>dt_dict, "geo"=>geodata, "sp"=>spreadparams)  
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
            sickday = zeros(intype, pop),   
            recov_day = zeros(intype, pop),  
            dead_day = zeros(intype, pop),   
            cluster = zeros(intype, pop), 
            s_d_comply = fill(:none, pop),  
            vax = zeros(intype, pop),   
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
            sickday = zeros(intype, pop))  

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


function send_risk_by_recv_risk(send_risk, recv_risk)
    recv_risk' .* send_risk  # (sickdaylim, agegrps)
end


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

    # named tuple doesn't result in type instability of Dict that requires "function barrier" to fix
    spreadparams = (
        # send_risk          = spread_params["send_risk"]::Vector{Float64},
        # recv_risk          = spread_params["recv_risk"]::Vector{Float64},
        contact_factors    = Dict(Int(k1) => Dict(string(k2) => Float64(v2) for (k2,v2) in v1) 
                                    for (k1, v1) in spread_params["contact_factors"]),
        touch_factors      = Dict(Int(k1) => Dict(string(k2) => Float64(v2) for (k2,v2) in v1) 
                                    for (k1, v1) in spread_params["touch_factors"]),
        shape              = spread_params["shape"]::Float64,
        riskmx             = send_risk::Array{Float64, 2}
        )
    
    return spreadparams
end


#####################################################################################
# dodgy math helper functions
#####################################################################################

function shifter(x::Array, oldmin, oldmax, newmin, newmax)
    newmin .+ (newmax - newmin) / (oldmax - oldmin) .* (x .- oldmin)
end


function shifter(x::Array, newmin, newmax)
    oldmin = minimum(x)
    oldmax = maximum(x)
    shifter(x, oldmin, oldmax, newmin, newmax)
end

"""
Warning: not general! works on dict with 2 levels and 
numerical values at the lower level.
"""
function limdict(dct::Dict, op::Function)
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
@inline function shifter(d::Dict, newmin, newmax)
    ret = deepcopy(d)
    oldmin = limdict(d, <)
    oldmax = limdict(d, >)
    
    @fastmath for k1 in keys(ret)
        for k2 in keys(ret[k1])
            x = ret[k1][k2]
            ret[k1][k2] = newmin + (newmax - newmin) / (oldmax - oldmin) * (x - oldmin)
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
    agegrp_filt_bit = Dict(agegrp => dat.agegrp .== agegrp for agegrp in agegrps)
    agegrp_filt_idx = Dict(agegrp => findall(agegrp_filt_bit[agegrp]) for agegrp in agegrps)
    return agegrp_filt_bit, agegrp_filt_idx
end
# agegrp_filt_bit, agegrp_filt_idx = precalc_agegrp_filt(ilmat);
