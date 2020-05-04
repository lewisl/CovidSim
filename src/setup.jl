######################################################################################
# setup and initialization functions
######################################################################################


function setup(n_days; geolim=15,  
    geofilename="../data/geo2data.csv", 
    dectreefilename="../parameters/dec_tree_all.csv",
    spfilename="../parameters/spread_params.toml")

    # geodata
        geodata = readgeodata(geofilename)
        numgeo = size(geodata,1)
        if geolim < numgeo
            numgeo = geolim
        end
        density_factor = shifter(geodata[:, density],0.9,1.25)
        geodata = [geodata density_factor]
        # fix dates   
        geodata[:, anchor] .= quickdate(geodata[:, anchor])
        geodata[:, restrict] .= quickdate(geodata[:, restrict])

    # simulation data matrix
        datadict = build_data(numgeo, n_days)
        openmx = datadict["openmx"]   # alias to make it easier to do initialization
        init_unexposed!(openmx, geodata, numgeo)

    # small parameters
        spread_params = read_spread_params(spfilename)

    # transition decision trees 
        dt = setup_dt(dectreefilename)

    # isolation probabilities: not sure we need this
        # iso_pr = build_iso_probs()

    return Dict("dat"=>datadict, "dt"=>dt, "geo"=>geodata, "sp"=>spread_params)  # , "iso_pr"=>iso_pr
end


"""
Convert a vector of dates from a csv file in format "mm/dd/yyyy"
to a vector of Julia numeric Date values in format yyyy-mm-dd
"""
function quickdate(strdates)  # 20x faster than the built-in date parsing, e.g.--runs in 5% the time
    ret = [parse.(Int,i) for i in split.(strdates, '/')]
    ret = [Date.(i[3], i[1], i[2]) for i in ret]
end


function init_unexposed!(dat, geodata, numgeo)
    for locale in 1:numgeo
        for agegrp in agegrps
            dat[locale][1, unexposed, agegrp] = floor(Int,age_dist[agegrp] * geodata[locale, popsize])
        end
    end
end


function build_data(numgeo, n_days)
    # openmx = zeros(Int, size(lags,1), size(conditions,1),  size(agegrps,1), numgeo)
    openmx = data_dict(numgeo, lags=size(lags,1), conds=size(conditions,1), agegrps=size(agegrps,1))
    isolatedmx = data_dict(numgeo, lags=size(lags,1), conds=size(conditions,1), agegrps=size(agegrps,1))

    cumhistmx = hist_dict(numgeo, n_days)
    newhistmx = hist_dict(numgeo, n_days)
    # isolatedhistmx = zeros(Int, size(conditions,1),  size(agegrps,1), numgeo, 1) # initialize for 1 day
    # return Dict("openmx"=>openmx, "isolatedmx"=>isolatedmx, "openhistmx"=>openhistmx, "isolatedhistmx"=>isolatedhistmx)
    return Dict("openmx"=>openmx, "isolatedmx"=>isolatedmx, "cumhistmx"=>cumhistmx, "newhistmx"=>newhistmx)
end

# one environment at a time
function data_dict(numgeo; lags=laglim, conds=8, agegrps=5)
    dat = Dict()
    for i = 1:numgeo
        dat[i] = zeros(Int, lags, conds, agegrps)
    end
    return dat       
end


function hist_dict(numgeo, n_days; conds=8, agegrps=5)
    dat = Dict()
    for i = 1:numgeo
        dat[i] = zeros(Int, conds, agegrps+1, n_days) # (conds, agegrps + 1, n_days) => (8, 6, 150)
    end
    return dat       
end


function readgeodata(filename)
    geodata = readdlm(filename, ','; header=true)[1]
end


function read_spread_params(spfilename)
    spread_params = TOML.parsefile(spfilename)

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

shifter(x::Array,a,b,c,d) = c .+ (d-c)/(b-a) .* (x .- a)

function shifter(x, newmin=0.9, newmax=1.5)
    oldmin = minimum(x)
    oldmax = maximum(x)
    shifter(x, oldmin, oldmax, newmin, newmax)
end

# not used
function build_iso_probs()
    # these don't need to sum to one in either direction!  independent events
    default_iso_pr = zeros(6,5)
    #                    enexp recover nil mild sick severe   
    default_iso_pr[:, a1] = [0.3, 0.3, 0.3, 0.4, 0.8, 0.9]  # this is for column a1
    default_iso_pr[:, a2] = [0.3, 0.3, 0.2, 0.3, 0.8, 0.9]
    default_iso_pr[:, a3] = [0.3, 0.3, 0.3, 0.4, 0.8, 0.9]
    default_iso_pr[:, a4] = [0.5, 0.5, 0.4, 0.6, 0.9, 0.9]
    default_iso_pr[:, a5] = [0.5, 0.5, 0.4, 0.6, 0.9, 0.9]

    iso_pr = Dict("default"=>default_iso_pr)

end
