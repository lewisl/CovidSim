using DelimitedFiles
using DataStructures
using Random
using Distributions
using StatsBase
using Printf


# control constants
const age_dist = [0.251, 0.271,   0.255,   0.184,   0.039]
const travprobs = [2.0, 3.0, 5.0, 5.0, 0.6] # by age group
const lags = 1:19   # rows

# geo data cubes (4th dim)
const id = 1
const state = 2
const city = 3
const size_cat = 4
const popsize =  5
# const locales = [1,2,3,4,5]  # rows

# condition_outcome columns
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
const infectious_cases = [nil, mild, sick, severe]

# evolve_prob_rows
const to_recovered = 1
const to_nil = 2
const to_mild = 3
const to_sick = 4
const to_severe = 5
const to_dead = 6

# agegrp channels at dimension 3
const a1 = 1 # 0-19
const a2 = 2 # 20-39
const a3 = 3 # 40-59
const a4 = 4 # 60-79
const a5 = 5 # 80+
const agegrps = [1,2,3,4,5]  
const ages = length(agegrps)
const contact_risk_by_age = [.05, .10, .15, .30, .45]

# this is too funky and hard to understand
# const age_amplify = Dict(1=>(1.3,.7), 2=>(1.1,.9), 3=>(1.0, 1.0), 4=>(0.9,1.1), 5=>(.8,1.2))

# inbound travelers: exogenous, domestic, outbound (assumed domestic) travelers
const travelq = Queue{NamedTuple{(:cnt, :from, :to, :agegrp, :cond),
    Tuple{Int64,Int64,Int64,Int64,String}}}()

function travitem(cnt, from, to, agegrp, cond)
    return (cnt=cnt, from=from, to=to, agegrp=agegrp, cond=cond)
end


######################################################################################
# setup and initialization functions
######################################################################################


function setup(geofilename, lim=10)

    geodata = readgeodata(geofilename)
    numgeo = size(geodata,1)
    if lim <= numgeo
        numgeo = lim
    end
    dats = build_data(numgeo)
    init_unexposed!(openmx, geodata, numgeo)

    ev_pr = build_evolve_probs()

    return (dats, ev_pr)
end


function init_unexposed!(dat, geodata, numgeo)
    for locale in 1:numgeo
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


function seed!(cnt, lag, status, cond, agegrp, locale; dat=openmx)
    if status == infectious
        input!.(cnt, status, agegrp,lag, locale)
        minus!.(cnt, unexposed, agegrp, lag, locale)
        input!.(cnt, cond, agegrp, lag, locale)
    elseif status == recovered
        dat[1, status, agegrp, locale] .= cnt
        # need to figure out what condition and lag get decremented!!
        dat[1, infectious, agegrp, locale] .-= cnt
    elseif status == dead 
        dat[lag, status, agegrp, locale] .= cnt
        # need to figure out what condition and lag get decremented!!
        dat[lag, infectious, agegrp, locale] .-= cnt        
    end
end

function build_evolve_probs()
    # evolve_probs will be of dims (4,6,5) => source conditions, destination conditions, agegrp planes
    # array labels

        # rows
        # to_recovered = 1
        # to_nil = 2
        # to_mild = 3
        # to_sick = 4
        # to_severe = 5
        # to_dead = 6
        # agegroups are columns: use a1,a2,a3,a4,a5 above
    # NOTE!!!: the acual arrays are what you see below transposed!

    nil_ev_pr = zeros(6,5)
    #               recover nil mild  sick severe dead
    nil_ev_pr[:, a1] = [0.7, 0.0, 0.2, 0.1, 0.0, 0.0]
    nil_ev_pr[:, a2] = [0.6, 0.0, 0.2, 0.2, 0.0, 0.0]
    nil_ev_pr[:, a3] = [0.5, 0.0, 0.3, 0.2, 0.0, 0.0]
    nil_ev_pr[:, a4] = [0.4, 0.0, 0.3, 0.3, 0.0, 0.0]
    nil_ev_pr[:, a5] = [0.3, 0.0, 0.4, 0.3, 0.0, 0.0]

    mild_ev_pr = zeros(6,5)
    mild_ev_pr[:, a1] = [0.3, 0.4, 0.0, 0.25, 0.05, 0.0]
    mild_ev_pr[:, a2] = [0.2, 0.4, 0.0, 0.3, 0.1, 0.0]
    mild_ev_pr[:, a3] = [0.2, 0.3, 0.0, 0.4, 0.1, 0.0]
    mild_ev_pr[:, a4] = [0.1, 0.2, 0.0, 0.5, 0.2, 0.0]
    mild_ev_pr[:, a5] = [0.05, 0.1, 0.0, 0.65, 0.2, 0.0]

    sick_ev_pr = zeros(6,5)
    sick_ev_pr[:, a1] = [0.2, 0.0, 0.7, 0.0, 0.1, 0.0]
    sick_ev_pr[:, a2] = [0.15, 0.0, 0.65, 0.0, 0.2, 0.0]
    sick_ev_pr[:, a3] = [0.1, 0.0, 0.55, 0.0, 0.3, 0.05]
    sick_ev_pr[:, a4] = [0.1, 0.0, 0.45, 0.0, 0.3, 0.15]
    sick_ev_pr[:, a5] = [0.0, 0.0, 0.3, 0.0, 0.4, 0.3]

    severe_ev_pr = zeros(6,5)
    severe_ev_pr[:, a1] = [0.0, 0.0, 0.0, 0.9, 0.0, 0.1]  # very few a1's ever get here
    severe_ev_pr[:, a2] = [0.0, 0.0, 0.0, 0.85, 0.0, 0.15]
    severe_ev_pr[:, a3] = [0.0, 0.0, 0.0, 0.85, 0.0, 0.15]
    severe_ev_pr[:, a4] = [0.0, 0.0, 0.0, 0.8, 0.0, 0.2]
    severe_ev_pr[:, a5] = [0.0, 0.0, 0.0, 0.0, 0.75, 0.25]

    ev_pr = Dict(nil=>nil_ev_pr, mild=>mild_ev_pr, sick=>sick_ev_pr, severe=>severe_ev_pr)

    return ev_pr
end




#function daystep()  # one day of simulation

    # loop across days of the simulation
        # function hooks for:...
        # seeding events--especially inbound international travel
        # outbreaks--when we don't have pockets implemented
        # rules--that restrict movement and affect touches and infection contacts

    # travel: some people travel in; distribution of them are infectious
        # pull item from queue:  add to destination: totalpop, agegrp/condition, agegrp/outcome
        #                        remove from source: totalpop, agegrp/condition, agegrp/outcome
        # add to log

    # DONE = evolve! distribute residents across all lags and conditions
            # start from last lag, by age group:
                    # severe distribute to sick, severe or die remain at lag 18
                    # sick distribute to severe, mild remain at lag 18
                    # mild distribute to sick, recovered (no nil at this point) at lag 18
                    # nil distribute to nil, mild, sick, severe at lags 2-6 
            # for earlier lags:
                    # severe distribute to sick, severe or die => advance one lag
                    # sick distribute to sick, severe, mild => advance one lag
                    # mild distribute to sick, recovered => advance one lag
                    # nil distribute to nil, mild, recovered => advance one lag

            # go dead => remove from previous condition; remove from resident, remove from totalpop:  make sure we can count # dead today (or more than yesterday)
            # go recovered => remove from previous condition; remove from exposed? remove from infectious    

    # loop across all locales to evolve all    


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


function evolve_all(locale)

    for agegrp in agegrps
    end

end

# TODO:  use predetermined distribute probs instead of equally likely
# function distribute(from_cond, to_conds, cnt)   # severe to severe, sick, die, mild
#     numcells = length(to_conds)
#     probvec = fill(1.0/float(numcells), numcells)
#     x = catprob(probvec, cnt)
# end

# TODO    This depends on what lag day it is!
function evolve!(ev_pr, conds, locale, lag19error; dat=openmx)  # also need to run for isolatedmx
    for cond in conds
        for agegrp in agegrps
                lastlag_distribute!(ev_pr, cond, agegrp, 19, locale, lag19error; dat=openmx)
        end
    end

    for lag = 18:-1:1
        for cond in conds
            for agegrp in agegrps

                # get the number of folks to distribute for the 
                folks = grab(cond,agegrp,lag,locale) # scalar
                # @printf("%12s %5d %5d %5d ",condnames[cond], agegrp,  lag, folks)
                # get the dist vector of folks to each outcome (6 outcomes)  
                    # distvec:   1: recovered 2: nil 3: mild 4: sick 5: severe 6: dead
                distvec = bucket(categorical_sample(ev_pr[cond][:,agegrp],folks), 6, 6)
                # @printf("%1s  %5d %5d %5d %5d %5d %5d \n", "|", distvec...)

                # for infectious states 
                plus!(distvec[2:5], infectious_cases, agegrp, lag+1,locale) # add to infectious cases for next lag
                minus!(sum(distvec[2:5]), cond, agegrp, lag, locale)  # everyone came from one cond from this lag

                # for recovered--add to recovered, substract from infectious and the specific infectious case
                out_of_infectious(distvec[to_recovered], cond, recovered, agegrp, lag, locale)

                # for dead--add to dead, substract from infectious and the specific infectious case
                out_of_infectious(distvec[to_dead], cond, dead, agegrp, lag, locale)
            end
        end
    end

end


function out_of_infectious(x, from, to, agegrp, lag, locale; dat=openmx)
    plus!(x, to, agegrp, 1, locale) # add to the target condition (either recovered or dead)
    minus!(x, from, agegrp, lag, locale)  # subtract from the infectious case where it came from
    minus!(x, infectious, agegrp, 1, locale)  # subtract from infectious
end


function lastlag_distribute!(ev_pr, cond, agegrp, lag, locale, lag19error; dat=openmx)
    condfolks = grab(cond, agegrp,19,locale)

    # not a distribution:  just 2 outcomes: either recovered or dead--not probabilistic
    recoverpr = ev_pr[cond][to_recovered, agegrp]
    recoverfolks = ceil(Int, recoverpr * condfolks)
    out_of_infectious(recoverfolks, cond, recovered, agegrp, lag, locale)

    deadpr = ev_pr[cond][to_dead, agegrp]
    deadfolks = ceil(Int, deadpr * condfolks)
    out_of_infectious(deadfolks, cond, dead, agegrp, lag, locale)

    residual = condfolks - recoverfolks - deadfolks
    @warn residual >= 0 "Uh-oh, some people infectious for more than 18 days" 
    lag19error(residual)  # accumulator

end


"""
    How far do the infectious people spread the virus to 
    previously unexposed people, by agegrp?
"""
function spread!(locale)
    # 180 microseconds
    # start with the number of infectious people      
    # now we ignore what their condition and age are is: TODO fix
    # should spread less for conditions in order: nil, mild, sick, severe
    # need to loop by condition

    # how many spreaders  TODO grab their condition.  Separate probs by condition
    spreaders = Int(sum(grab(infectious, agegrps, 1:19, locale)))
    if spreaders == 0
        return nothing
    end

    # how many people are touched   
    touched = how_many_touched(spreaders)

    # contact that causes infection
    # amplify probabilities with population density for locale
    byage = split_by_age(touched)
    for i in 1:length(byage) # this probabilisticly determines if contact resulted in contracting the virus
        byage[i] = ceil.(Int, rand(Binomial(byage[i], contact_risk_by_age[i])))
    end

    # move the people from unexposed:agegrp to infectious:agegrp and nil
    plus!.(byage, infectious, agegrps, 1, locale)  
    plus!.(byage, nil, agegrps, 1, locale)  
    
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
# TODO use population density as probability amplifier
function how_many_touched(touchers, scale=6)
    dgamma = Gamma(1.2, scale)  #shape, scale
    x = rand(dgamma,touchers);
    nums, bounds = histo(x)
    return sum(nums .* bounds)
end


"""
    For a locale, randomly choose the number of people from each agegroup with
    condition of {unexposed, infectious, recovered} who travel to each
    other locale. Add to the travelq.
"""
function travelout!(locale, numgeo, rules=[])
    # 10.5 microseconds for 5 locales
    # choose distribution of people traveling by age and condition:
        # unexposed, infectious, recovered -> ignore lag for now
    # TODO: more frequent travel to and from Major and Large cities
    locales =1:numgeo
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
function bucket(x, lim=5, bins = 5)
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


function categorical_sample(probvec, trials)
    @assert isapprox(sum(probvec), 1.0) "target vector must sum to 1.0"
    x = rand(Categorical(probvec), trials)
end

function accumulator(i)  # call the first time with zero?
               f(n) = i += n
               return f
       end

lag19error = accumulator(0)  # do judgment parameters result in people still sick after 19 days?

####################################################################################
#   convenience functions for reading and inputting population statistics
####################################################################################


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