######################################################################################
# setup and initialization functions
######################################################################################

const mean_density = 3316  # see Excel for backup


function setup(geofilename; dectreefilename="dec_tree_all.csv", geolim=15)

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
    datadict = build_data(numgeo)
    openmx = datadict["openmx"]
    isolatedmx = datadict["isolatedmx"]
    init_unexposed!(openmx, geodata, numgeo)

    # transition decision trees 
    dt = setup_dt(dectreefilename)

    # # isolation probabilities: not sure we need this
    # iso_pr = build_iso_probs()

    return Dict("dat"=>datadict, "dt"=>dt, "geo"=>geodata)  # , "iso_pr"=>iso_pr
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


function build_data(numgeo)
    # openmx = zeros(Int, size(lags,1), size(conditions,1),  size(agegrps,1), numgeo)
    openmx = data_dict(numgeo, lags=size(lags,1), conds=size(conditions,1), agegrps=size(agegrps,1))
    isolatedmx = data_dict(numgeo, lags=size(lags,1), conds=size(conditions,1), agegrps=size(agegrps,1))

    # openhistmx = zeros(Int, size(conditions,1), size(agegrps,1), numgeo, 1) # initialize for 1 day
    # isolatedhistmx = zeros(Int, size(conditions,1),  size(agegrps,1), numgeo, 1) # initialize for 1 day
    # return Dict("openmx"=>openmx, "isolatedmx"=>isolatedmx, "openhistmx"=>openhistmx, "isolatedhistmx"=>isolatedhistmx)
    return Dict("openmx"=>openmx, "isolatedmx"=>isolatedmx)
end

# one environment at a time
function data_dict(numgeo; lags=laglim, conds=8, agegrps=5)
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
        Unexposed
        Infectious
        Recovered
        Dead
        Nil
        Mild
        Sick
        Severe
        Travelers
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

    template = DataFrame([series_colnames[i] => Int[] for i in keys(series_colnames)]...)

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


# calculate density_factor in setup, per locale
test_density = rand((5000:3_000_000),20)  # for real, use US Census data

function minmax_norm(x)
    x_max = maximum(x, dims=1)
    x_min = minimum(x, dims=1)
    minmax = (x .- x_min) ./ (x_max .- x_min .+ 1e-08)
end

function scale_minmax(x, newmin, newmax) 
    round.(x .* (newmax - newmin) .+ newmin, digits=2)
end

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
