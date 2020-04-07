######################################################################################
# setup and initialization functions
######################################################################################


function setup(geofilename; dectreefilename="dec_tree_all.csv", node_starts_filename="dec_tree_starts.csv", geolim=10)

    geodata = readgeodata(geofilename)
    numgeo = size(geodata,1)
    if geolim <= numgeo
        numgeo = geolim
    end
    datsdict = build_data(numgeo)
    openmx = datsdict["openmx"]
    init_unexposed!(openmx, geodata, numgeo)

    # transition decision trees 
    dt = setup_dt(dectreefilename; node_starts_filename=node_starts_filename)

    iso_pr = build_iso_probs()

    return Dict("dat"=>datsdict, "dt"=>dt, "iso_pr"=>iso_pr)
end


function init_unexposed!(dat, geodata, numgeo)
    for locale in 1:numgeo
        for agegrp in agegrps
            dat[locale][1, unexposed, agegrp] = floor(Int,age_dist[agegrp] * geodata[locale, popsize])
        end
    end
end


function build_data(numgeo)
    # openmx = zeros(Int, size(lags,1), size(conditions,1),  size(agegrps,1), numgeo)
    openmx = data_dict(numgeo, lags=size(lags,1), conds=size(conditions,1), agegrps=size(agegrps,1))
    isolatedmx = data_dict(numgeo, lags=size(lags,1), conds=size(conditions,1), agegrps=size(agegrps,1))

    # isolatedmx = zeros(Int, size(lags,1), size(conditions,1),  size(agegrps,1), numgeo)
    # openhistmx = zeros(Int, size(conditions,1), size(agegrps,1), numgeo, 1) # initialize for 1 day
    # isolatedhistmx = zeros(Int, size(conditions,1),  size(agegrps,1), numgeo, 1) # initialize for 1 day
    # return Dict("openmx"=>openmx, "isolatedmx"=>isolatedmx, "openhistmx"=>openhistmx, "isolatedhistmx"=>isolatedhistmx)
    return Dict("openmx"=>openmx, "isolatedmx"=>isolatedmx)
end

# one environment at a time
function data_dict(numgeo; lags=19, conds=10, agegrps=5)
    dat = Dict()
    for i = 1:numgeo
        dat[i] = zeros(Int, lags, conds, agegrps)
    end
    return dat       
end

"""
    Build container for data series of simulation outcomes by day.
    Top index is an integer that is the ID of a locale or 0 for total across locales.
    2nd index is either :new or :cum.
    Values at the leaf of :new is a DataFrame of new additions to a series.
    Values at the leaf of :cum is a DataFrame of the cumulative values of a series.
    Current series columns are:
        Travelers
        Infected
        Nil
        Mild
        Sick
        Severe
        Dead
        Recovered
        Isolated
    Rows are days of the simulation.
"""
function build_series(locales)
#=
columnnames = ["x", "y", "z"]
columns = [Symbol(col) => Float64[] for col in columnnames]
df1 = DataFrame(columns...)
df2 = DataFrame(columns...)
series_colnames = Dict( 1=>:Unexposed,  2=>:Infectious, 3=>:Recovered, 4=>:Dead, 5=>:Nil, 6=>:Mild, 7=>:Sick,
        8=>:Severe,  9=>:Travelers, 10=>:Isolated)

=#


    dseries = Dict{Int,Dict}()
    # columns = [series_colnames[i] => Int[] for i in 1:length(series_colnames)]

    template = DataFrame([series_colnames[i] => Int[] for i in 1:length(series_colnames)]...)

    # new = DataFrame(Travelers=Int[], Unexposed=Int[], Infected=Int[], Nil=Int[], Mild=Int[], Sick=Int[],
    #     Severe=Int[], Dead=Int[], Recovered=Int[], Isolated=Int[])
    # cum = DataFrame(Travelers=Int[], Unexposed=Int[],Infected=Int[], Nil=Int[], Mild=Int[], Sick=Int[],
    #     Severe=Int[], Dead=Int[], Recovered=Int[], Isolated=Int[]) # do by agegrp, total, by gender?
    for i in [0, locales...]
        dseries[i]=Dict{Symbol,DataFrame}()
        dseries[i][:new] = deepcopy(template)
        dseries[i][:cum] = deepcopy(template)
    end
    return dseries
end


function readgeodata(filename)
    geodata = readdlm(filename, ','; header=true)[1]
end


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
