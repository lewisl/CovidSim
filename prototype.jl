using DelimitedFiles
using DataStructures
using Random
using Distributions
using StatsBase


# control constants
const age_dist = [0.251, 0.271,   0.255,   0.184,   0.039]
const travprobs = [2.0, 3.0, 5.0, 5.0, 0.6] # by age group
const lags = 1:19

# geo data columns
const id = 1
const state = 2
const city = 3
const size_cat = 4
const popsize =  5
const locales = [1,2,3,4,5]  # rows

# condition_outcome rows
const unexposed = 1
const infectious = 2
const recovered = 3
const dead = 4
const nil = 5
const mild = 6
const sick = 7
const severe = 8
const conditions = [unexposed, infectious, recovered, dead, nil, mild, sick, severe]
const condnames = Dict(1=>"unexposed",2=>"infectious",3=>"recovered", 4=>"dead",
                       5=>"nil",6=>"mild",7=>"sick",8=>"severe")

# agegrp channels at dimension 3
const a1 = 1 # 0-19
const a2 = 2 # 20-39
const a3 = 3 # 40-59
const a4 = 4 # 60-79
const a5 = 5 # 80+
const agegrps = [1,2,3,4,5]
const ages = length(agegrps)
const contact_risk_by_age = [.05, .10, .15, .30, .45]

# inbound travelers: exogenous, domestic, outbound (assumed domestic) travelers
const travelq = Queue{NamedTuple{(:cnt, :from, :to, :agegrp, :cond),
    Tuple{Int64,Int64,Int64,Int64,String}}}()

function travitem(cnt, from, to, agegrp, cond)
    return (cnt=cnt, from=from, to=to, agegrp=agegrp, cond=cond)
end



######################################################################################
# setup and initialization functions
######################################################################################

function init_unexposed!(dat, geodata)
    for locale in locales
        for agegrp in agegrps
            dat[1, unexposed, agegrp, locale] = floor(Int,age_dist[agegrp] * geodata[locale, popsize])
        end
    end
end


function build_data(numgeo)
    openmx = zeros(Int, size(lags,1), size(conditions,1),  size(agegrps,1), numgeo)
    isolatedmx = zeros(Int, size(lags,1), size(conditions,1),  size(agegrps,1), numgeo)
    openhistmx = zeros(Int, size(conditions,1), size(agegrps,1), numgeo, 1) # initialize for 1 day
    isolatedhistmx = zeros(Int, size(conditions,1),  size(agegrps,1), numgeo, 1) # initialize for 1 day
    return (openmx, isolatedmx, openhistmx, isolatedhistmx)
end

function readgeodata(filename)
    geodata = readdlm(filename, ','; header=true)[1]
end

function setup(geofilename)

    geodata = readgeodata(geofilename)
    numgeo = size(geodata,1)
    openmx, isolatedmx, openhistmx, isolatedhistmx = build_data(numgeo)
    init_unexposed!(openmx, geodata)

    return (openmx, isolatedmx, openhistmx, isolatedhistmx)
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


    # DONE:  distribution of people(not isolated) contact distribution of residents across age&condition&lags
            # DONE = spread!  touchers -> num_touched -> split by age/condition
            # DONE in spread!  of those contacted, a distribution become exposed & nil at lag 0


    # DONE = travelout! some residents travel out
        # choose distribution of people traveling by age and condition
        # add to queue    

    # gather stats

#end



####################################################################################
#   functions for simulation events
####################################################################################


function evolve(locale)

    for agegrp in agegrps
    end

end

# TODO:  use predetermined distribute probs instead of equally likely
function distribute(from_cond, to_conds, cnt)   # severe to severe, sick, die, mild
    numcells = length(to_conds)
    probvec = fill(1.0/float(numcells), numcells)
    x = catprob(probvec, cnt)
end


"""
How far do the infectious people spread the virus into previously unexposed people, by agegrp?
"""
function spread!(locale)
    # 110 microseconds
    # start with the number of infectious people      
    # now we ignore what their condition and age are is: TODO fix
    # should spread less for conditions in order: nil, mild, sick, severe
    # need to loop by condition
    spreaders = Int(sum(grab(infectious, 1:5,1:11, locale)))
    if spreaders == 0
        return nothing
    end
    # how many people are touched 
    touched = how_many_touched(spreaders)

    # by age group
    byage = split_by_age(touched)
    for i in 1:length(byage)
        byage[i] = ceil.(Int, rand(Binomial(byage[i], contact_risk_by_age[i])))
    end

    # move the people from unexposed:agegrp to infectious:agegrp
    plus!.(byage, infectious, agegrps, 1, locale)  
    
    return nothing
end


function split_by_age(cnt)::Array{Int64,1}
    numcells = ages
    probvec = age_dist
    x = catprob(probvec, cnt)
    # nums, _ = histo(x)
    nums = bucket(x,lim=5, bins=5)
    return nums
end

# 5.7 microseconds for 100 touchers
function how_many_touched(touchers, scale=6)
    dgamma = Gamma(1.2, scale)  #shape, scale
    x = rand(dgamma,touchers);
    nums, bounds = histo(x)
    return sum(nums .* bounds)
end

# 10.5 microseconds for 5 locales
"""
For a locale, randomly choose the number of people from each agegroup with
condition of {unexposed, exposed, infectious, recovered} who travel to each
other locale. Add to the travelq.
"""
function travelout!(locale, rules=[])
    # choose distribution of people traveling by age and condition:
        # unexposed, infectious, recovered -> ignore lag for now
    travdests = copy(locales)
    deleteat!(travdests,findfirst(isequal(locale), travdests))
    bins = lim = length(travdests) + 1
    for cond in [unexposed, infectious, recovered]
        name = condnames[cond]
        for agegrp in agegrps
            numfolks = sum(grab(cond, agegrp, lags, locale)) # this locale, all lags
            travcnt = floor(Int,prob(travprobs[agegrp]) * numfolks)  # TODO use probability
            x = rand(travdests, travcnt)  # randomize across destinations
            bydest = bucket(x, lim=lim, bins=bins)
            for dest in 1:length(bydest)
                isempty(bydest) && continue
                cnt = bydest[dest]
                iszero(cnt) && continue
                enqueue!(travelq, travitem(cnt, locale, dest, agegrp, name))
            end
        end
    end

end

# this could still be faster by memoizing x--don't look again at items already counted.
# discrete integer histogram
function bucket(x::Array; lim=5, bins = 5)
    ret = zeros(Int, bins)
    comp = collect(1:lim)
    for j in comp
        n = count(isequal(j),x)
        ret[j] += n
    end
    return ret
end

function histo(x)
    big = ceil(maximum(x))
    bins = Int(big)
    sm = floor(minimum(x))
    ret = zeros(Int, bins)
    binbounds = collect(1:bins)
    for i = 1:bins
        n = count(x -> i-1 < x <= i,x)
        ret[i] = n
    end
    return ret, binbounds
end


function prob(target)
    @assert 0.0 <= target <= 99.0 "target must be between 0.0 and 99.0"
    dgamma = Gamma(1.2,target)
    pr = rand(dgamma, 1)[1] / 100.0
end

function catprob(probvec, trials)
    @assert sum(probvec) == 1.0 "target vector must sum to 1.0"
    x = rand(Categorical(probvec), trials)
end

####################################################################################
#   convenience functions for reading and inputting population statistics
####################################################################################

# single age, single lag, one locale
# example: grab(exposed, 1, 1, 1:3)
function grab(condition, agegrp, lag, locale; dat=openmx)
    return dat[lag, condition, agegrp, locale]
end


function input!(val, condition, agegrp, lag, locale; dat=openmx)
    dat[lag, condition, agegrp, locale] = val
end


function plus!(val, condition, agegrp, lag, locale; dat=openmx)
    dat[lag, condition, agegrp,locale] += val
end


function minus!(val, condition, agegrp, lag, locale; dat=openmx)
    dat[lag, condition, agegrp, locale] -= val
end