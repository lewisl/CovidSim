using DelimitedFiles
using DataStructures
using Random
using Distributions
using StatsBase


mutable struct People # indexes to population data
    a1::Int64  # age group 1: 0 to 19
    a2::Int64  # age group 1: 20 to 39
    a3::Int64  # age group 1: 40 to 59
    a4::Int64  # age group 1: 60 to 79
    a5::Int64  # age group 1: 80+
    lag::Bool
    name::String
end
    
    # alternate creator method to build array of indices
    People(name, lag::Bool) = People(
            0, 0, 0, 0, 0, lag, name)

# control constants
const age_dist = [0.251, 0.272,   0.255,   0.184,   0.039]
const agegrps = [1,2,3,4,5]
const travprobs = [2.0, 3.0, 5.0, 5.0, 0.6] # by age group
const locales = [1,2,3,4,5]
const lagrange = 1:19

# row indices of data matrix won't change
const id = 1
const size_cat = 2
const popsize =  3

# population conditions
const unexposed = People("unexposed", false) # active or cumulative
const exposed = People("exposed", true)  # aging over 18 days; means not_infectious
const infectious = People("infectious", true) # aging over 18 days
const recovered = People("recovered", false) # active or cumulative

# population outcomes
const mild = People("mild", true)
const sick = People("sick", true)
const severe = People("severe", true)
const dead = People("dead", false) # active or cumulative

# population actions
const isolated_exp = People("isolated_exp", true) # needs to be associated with condition
const isolated_un  = People("isolated_un", true)  # must be recovered or unexposed

# inbound travelers: exogenous, domestic, outbound (assumed domestic) travelers
const intravel = People("intravel", false) # exogenous travelers in
const crosstravel = People("crosstravel", false)
const outtravel = People("outtravel", false)
const travelq = Queue{NamedTuple{(:cnt, :from, :to, :agegrp, :cond),
    Tuple{Int64,Int64,Int64,Int64,String}}}()

function travitem(cnt, from, to, agegrp, cond)
    return (cnt=cnt, from=from, to=to, agegrp=agegrp, cond=cond)
end



######################################################################################
# setup and initialization functions
######################################################################################

function push_to_unexposed!(dat)
    for locale in locales
        for agegrp in agegrps
            dat[getfield(unexposed,agegrp), locale] = floor(Int,age_dist[agegrp] * dat[popsize, locale])
        end
    end
end


function build_matrix(symdata)
    list = [unexposed, exposed, infectious, isolated, recovered, mild, sick,
            severe, dead, intravel, crosstravel, outtravel]
    m,n = size(symdata)
    for item in list
        symdata = cat(symdata,zeros(Int64,5,n), dims=1)
        len = size(symdata,1)
        item.a1 = len-4; item.a2 = len-3; item.a3 = len - 2; item.a4 = len -1; item.a5 = len;
        if item.lag == true
            beg = item.a1
            fin = beg + 4
            symdata[beg:fin,:] .= [zeros(Int64,19)]
        end
    end
    return symdata
end


function setup(symdata_filename)
    symdata = readdlm(symdata_filename, ',')
    symdata = convert(Array{Any}, symdata)
    symdata = build_matrix(symdata)
    push_to_unexposed!(symdata)
    return symdata
end


#function daystep()  # one day of simulation

    # travel: some people travel in; distribution of them are infectious
        # pull item from queue:  add to destination: totalpop, agegrp/condition, agegrp/outcome
        #                        remove from source: totalpop, agegrp/condition, agegrp/outcome
        # add to log


    # evolve residents across all lags
            # start from last lag, by age group:
                    # severe distribute to sick, severe or die remain at lag 18
                    # sick distribute to severe, mild remain at lag 18
                    # mild distribute to sick, recovered (no nil at this point) at lag 18
                    # exposed distribute to exposed, mild, sick, severe at days 1-6 
            # for earlier lags:
                    # severe distribute to sick, severe or die => advance one lag
                    # sick distribute to sick, severe, mild => advance one lag
                    # mild distribute to sick, recovered => advance one lag
                    # nil distribute to nil, mild, recovered => advance one lag

            # go dead => remove from previous condition; remove from resident, remove from totalpop:  make sure we can count # dead today (or more than yesterday)
            # go recovered => remove from previous condition; remove from exposed? remove from infectious        


    # contact:  distribution of people(not isolated) contact distribution of residents across age&condition&lags


    # of those contacted, a distribution become exposed & nil at lag 0


    # DONE = travelout! some residents travel out
        # choose distribution of people traveling by age and condition
        # add to queue    

    # gather stats

#end



####################################################################################
#   functions for simulation events
####################################################################################


function evolve(locale)

    for agegrp = agegrps


end


function evolve_one(source, to_conds, locale, agegrp)   # severe to severe, sick, die, mild

    for cond in to_conds  #

end

# 16 microseconds for 5 locales
"""
For a locale, randomly choose the number of people from each agegroup with
condition of {unexposed, exposed, infectious, recovered} who travel to each
other locale. Add to the travelq.
"""
function travelout!(locale, rules=[])
    # choose distribution of people traveling by age and condition:
        # unexposed, exposed, infectious, recovered -> ignore lag for now
    # ret = []
    travdests = locales
    deleteat!(travdests,findfirst(isequal(locale), travdests))
    for peeps in [exposed, unexposed, infectious, recovered]
        for agegrp = agegrps
            if peeps.lag == true
                numfolks = sum(grab(peeps, agegrp, lagrange, locale)) # this locale, all lags
            else
                numfolks = grab(peeps, agegrp, 1, locale) # this locale, all lags
            end                
            travcnt = floor(Int,prob(travprobs[agegrp]) * numfolks)
            x = rand(travdests, travcnt)
            bydest = bucket(x, lim=length(travdests), bins=length(travdests))
            for dest in 1:length(bydest)
                isempty(bydest) && continue
                cnt = bydest[dest]
                iszero(cnt) && continue
                enqueue!(travelq, travitem(cnt, locale, dest, agegrp, peeps.name))
            end
        end
    end

end

# this could still be faster by memoizing x--don't look again at items already counted.
function bucket(x::Array; lim=5, bins = 5)
    ret = zeros(Int, bins)
    comp = collect(1:lim)
    for j in comp
        n = count(isequal(j),x)
        ret[j] += n
    end
    return ret
end

function prob(target)
    @assert 0.0 <= target <= 99.0 "target must be between 0.0 and 99.0"
    dgamma = Gamma(1.2,target)
    pr = rand(dgamma, 1)[1] / 100.0
end

####################################################################################
#   convenience functions for reading and inputting population statistics
####################################################################################

# single age, single lag, one locale
# example: grab(exposed, 1, 1, 1:3)
function grab(item::People, agegrp::Int, lag::Int, locale::Int; dat=symdata)
    return dat[getfield(item, agegrp), locale][lag]
end

# single age, multiple lags, one locale
# example: grab("exposed", 1, 1:10, 1)  for lag 0 to 9 days
function grab(item::People, age::Int, lag::UnitRange{Int64}, locale::Int; dat=symdata)
    return dat[getfield(item, age), locale][lag]
end 

# 1 values to single age, single lag, one locale
function input!(val, item::People, age::Int, lag::Int, locale::Int; dat=symdata)
    dat[getfield(item, age), locale][lag] = val
end

# TODO--doesn't work right.  assigns to all locales
# 1 or more values to single age, multiple lags, one locale
function input!(val, item::People, age::Int, lag::UnitRange{Int64}, locale; dat=symdata)
    dat[getfield(item, age), locale][lag] .= val
end