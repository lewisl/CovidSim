# TODO
# more info
    # get fatality rate by age and co-morbidity CDC, Italian NIH
    # by agegroup, hospitalization %, ICU admission %, fatality %
    # UW virology, expansion of deaths by state on log chart
# do we really need hist matrices that copy everything?


using DelimitedFiles
using DataStructures
using DataFrames
using Random
using Distributions
using StatsBase
using Printf

mutable struct SimEnvironment  # maybe, maybe not...?
    simdday::Int
    opendat::Array
    isodat::Array
    newstats::DataFrame
    cumstats::DataFrame
end


# control constants
const age_dist = [0.251, 0.271,   0.255,   0.184,   0.039]
const lags = 1:19   # rows

# geo data cubes (4th dim)
const id = 1
const state = 2
const city = 3
const size_cat = 4
const popsize =  5
# const locales = [1,2,3,4,5]  # rows

# population centers sizecats
const major = 1
const large = 2
const medium = 3
const small = 4
const smaller = 5
const rural = 6

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
const evolve_cases = [recovered, nil, mild, sick, severe, dead]
const series_colnames = Dict(1=>:Unexposed,  2=>:Infectious, 3=>:Recovered, 4=>:Dead, 5=>:Nil, 6=>:Mild, 7=>:Sick,
        8=>:Severe,  9=>:Travelers, 10=>:Isolated)

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

# traveling constants
    # inbound travelers: exogenous, domestic, outbound (assumed domestic) travelers
    const travelq = Queue{NamedTuple{(:cnt, :from, :to, :agegrp, :lag, :cond),
        Tuple{Int64,Int64,Int64,Int64,Int64,String}}}()

    function travitem(cnt, from, to, agegrp, lag, cond)
        return (cnt=cnt, from=from, to=to, agegrp=agegrp, lag=lag, cond=cond)
    end

    const travprobs = [1.0, 2.0, 3.0, 3.0, 0.4] # by age group


# isolation
    const isolatedq = Queue{NamedTuple{(:cnt, :cond, :agegrp, :lag, :locale), 
         Tuple{Int64,Int64,Int64,Int64,Int64,}}}()

    function iso_item(cnt, cond, agegrp, lag, locale)
        return (cnt=cnt, cond=cond, agegrp=agegrp, lag=lag, locale=locale)
    end

# tracking statistics

    const newstatq = DataFrame(simday=Int[], cnt=Int[], locale=Int[], tocond=Int[])

    TravelStat = typeof((;cnt=0, locale=0, cond=0, to=0))
    function travelstat(;cnt=0, locale=0, cond=0, to=0)
        (cnt=cnt, locale=locale, cond=cond, to=to)
    end

    EvolveStat = typeof((;cnt=0, locale=0, tocond=5))
    function evolvestat(; cnt=0, locale=0, tocond=5)
        (cnt=cnt, locale=locale, tocond=tocond)
    end

    IsolateStat = typeof((;cnt=0, locale=0))
    function isolatestat(; cnt=0, locale=0)
        (cnt=cnt, locale=locale)
    end

    SpreadStat = typeof((;cnt=0, locale=0, tocond=5))
    function spreadstat(;cnt=0,locale=0,tocond=5)
       (cnt=cnt, locale=locale, tocond=tocond)
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
    datsdict = build_data(numgeo)
    openmx = datsdict["openmx"]
    init_unexposed!(openmx, geodata, numgeo)

    ev_pr = build_evolve_probs()

    iso_pr = build_iso_probs()

    return Dict("dat"=>datsdict, "ev_pr"=>ev_pr, "iso_pr"=>iso_pr)
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
    return Dict("openmx"=>openmx, "isolatedmx"=>isolatedmx, "openhistmx"=>openhistmx, "isolatedhistmx"=>isolatedhistmx)
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
function build_series(numgeo)
#=
columnnames = ["x", "y", "z"]
columns = [Symbol(col) => Float64[] for col in columnnames]
df1 = DataFrame(columns...)
df2 = DataFrame(columns...)
series_colnames = Dict(1=>"Travelers", 2=>"Unexposed", 3=>"Infected", 4=>"Nil", 5=>"Mild", 6=>"Sick",
        7=>"Severe", 8=>"Dead", 9=>"Recovered", 10=>"Isolated")

=#

    dseries = Dict{Int,Dict}()
    columns = [series_colnames[i] => Int[] for i in 1:length(series_colnames)]
    new = DataFrame(columns...)
    cum = DataFrame(columns...)
    # new = DataFrame(Travelers=Int[], Unexposed=Int[], Infected=Int[], Nil=Int[], Mild=Int[], Sick=Int[],
    #     Severe=Int[], Dead=Int[], Recovered=Int[], Isolated=Int[])
    # cum = DataFrame(Travelers=Int[], Unexposed=Int[],Infected=Int[], Nil=Int[], Mild=Int[], Sick=Int[],
    #     Severe=Int[], Dead=Int[], Recovered=Int[], Isolated=Int[]) # do by agegrp, total, by gender?
    for i in 0:numgeo
        dseries[i]=Dict{Symbol,DataFrame}()
        dseries[i][:new] = new
        dseries[i][:cum] = cum
    end
    return dseries
end

function readgeodata(filename)
    geodata = readdlm(filename, ','; header=true)[1]
end


function seed!(cnt, lag, conds, agegrps, locales; dat=openmx)
    @assert length(lag) == 1 "input only one lag value"    @warn "Seeding is for testing and may results in case counts out of balance"
    for loc in locales
        for cond in conds
            if cond in [nil, mild, sick, severe]
                input!.(cnt, cond, agegrps, lag, loc, dat=dat)
                minus!.(cnt, unexposed, agegrps, lag, loc, dat=dat)
            elseif cond in [recovered, dead]
                input!.(cnt, cond, agegrp, lag, loc, dat=dat)
                # need to figure out what condition and lag get decremented!!
            end
        end
        update_infectious!(loc, dat = dat)
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

        # columns  must sum to 1: for a given condition and age, all who evolve go somewhere
                # people don't all evolve--they stay in place
    # NOTE!!!: the acual arrays are what you see below transposed like the one for nil_ev_pr

 
 
    #                     a1   a2   a3   a4   a5 
        nil_ev_pr =     [0.2  0.2  0.1  0.1  0.1; # recover
                         0.5  0.5  0.4  0.3  0.3; # nil
                         0.2  0.2  0.3  0.3  0.3; # mild 
                         0.1  0.1  0.2  0.3  0.3; # sick 
                         0.0  0.0  0.0  0.0  0.0; # severe 
                         0.0  0.0  0.0  0.0  0.0] # dead


    mild_ev_pr = zeros(6,5)
    #                 recover nil mild  sick severe dead
    mild_ev_pr[:, a1] = [0.1, 0.3, 0.45, 0.15, 0.0, 0.0]
    mild_ev_pr[:, a2] = [0.1, 0.25, 0.45, 0.15, 0.05, 0.0]
    mild_ev_pr[:, a3] = [0.05, 0.25, 0.4, 0.2, 0.1, 0.0]
    mild_ev_pr[:, a4] = [0.05, 0.2, 0.4, 0.25, 0.1, 0.0]
    mild_ev_pr[:, a5] = [0.05, 0.2, 0.3, 0.35, 0.1, 0.0]

    sick_ev_pr = zeros(6,5)
    sick_ev_pr[:, a1] = [0.0, 0.2, 0.3, 0.4, 0.1, 0.0]
    sick_ev_pr[:, a2] = [0.0, 0.17, 0.32, 0.4, 0.11, 0.0]
    sick_ev_pr[:, a3] = [0.0, 0.1, 0.35, 0.35, 0.2, 0.0]
    sick_ev_pr[:, a4] = [0.0, 0.1, 0.34, 0.35, 0.2, 0.01]
    sick_ev_pr[:, a5] = [0.0, 0.07, 0.26, 0.35, 0.3, 0.02]

    severe_ev_pr = zeros(6,5)
    severe_ev_pr[:, a1] = [0.0, 0.0, 0.25, 0.425, 0.32, 0.005]  # very few a1's ever get here
    severe_ev_pr[:, a2] = [0.0, 0.0, 0.25, 0.38, 0.36, 0.01]
    severe_ev_pr[:, a3] = [0.0, 0.0, 0.22, 0.35, 0.42, 0.01]
    severe_ev_pr[:, a4] = [0.0, 0.0, 0.2, 0.285, 0.5, 0.015]
    severe_ev_pr[:, a5] = [0.0, 0.0, 0.15, 0.265, 0.57, 0.015]

    ev_pr = Dict(nil=>nil_ev_pr, mild=>mild_ev_pr, sick=>sick_ev_pr, severe=>severe_ev_pr)

    return ev_pr
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


#function daystep()  # one day of simulation

    # loop across days of the simulation
        # function hooks for:...
        # seeding events--especially inbound international travel
        # outbreaks--when we don't have pockets implemented
        # rules--that restrict movement and affect touches and infection contacts


    # need rules for how many people travel across locales; then rules for travel restrictions
         # DONE = travelin!: some people travel in; distribution of them are infectious
        # pull item from queue:  add to destination: agegrp/condition, agegrp/outcome
        #                        remove from source: agegrp/condition, agegrp/outcome
        # TODO add to log

    # some people will isolate: need to come up with rules for isolation, leaks from isolation 
        #   DONE = isolate_in put people into and out of isolation by agegroup, locale, condition, lag

    # DONE = evolve! distribute residents across all lags and conditions
        # openmx
        # isolatedmx
        # loop across all locales to evolve all    
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

            # TODO log new dead, new recovered  to change state 

    

    # DONE:  spread from infectious people to unexposed|recovered people(not isolated) 
            #contact distribution of residents across age&condition&lags
            # TODO: handle partial immunity of recovered
            # DONE = spread!  touchers -> num_touched -> split by age/condition
            # DONE in spread!  of those contacted, a distribution become exposed & nil at lag 0

            # TODO log new infectious


    # DONE = travelout! some residents travel out
        # choose distribution of people traveling by age and condition
        # add to queue    

    # TODO summarize current state values for unexposed, infectious, dead, recovered
        # add to series using change_state, summarize function

#end

# TODO Build locale data in a dataframe
#   cols: id, state, locale, pop, density, travelprobs, gender split, age dist (maybe)
# 
# 
# 
# 
# 

####################################################################################
#   functions for simulation events
####################################################################################
# things that cause condition changes: travel, spread, evolve, isolate

function evolve_all(locales=locales)

    # open environment
    for locale in locales
    end

    # isolated environment
    for locale in locales
    end

end


"""
    People who have become infectious evolve through cases from
    nil (asymptomatic) to mild to sick to severe, depending on their
    agegroup, days of being exposed, and statistics. The final outcomes
    are recovery or dying.
"""
function evolve!(ev_pr, locale; conds=infectious_cases, case="open", dat=openmx)  # TODO also need to run for isolatedmx

    # special handling because this is the last day back we track outcomes
    # lag = 19 # [recovered, nil, mild, sick, severe, dead]
    # for cond in conds
    #     for agegrp in agegrps
    #             lastlag_distribute!(ev_pr, cond, agegrp, lag, locale; dat=dat)
    #     end
    # end

    # for day 19
    #     implement decision point

    for lag = 18:-1:14
        for cond in conds  # [recovered, nil, mild, sick, severe, dead]
            for agegrp in agegrps

                dist_to_new_conditions!(cond, to_recovered:to_dead, agegrp, lag, locale, case=case, dat=dat)

            end
        end
    end    

    # for day 14
        # implement decision point

    # next 10 days to outcomes
    for lag = 13:-1:9
        for cond in conds  # [recovered, nil, mild, sick, severe, dead]
            for agegrp in agegrps

                dist_to_new_conditions!(cond, to_recovered:to_dead, agegrp, lag, locale, case=case, dat=dat)

            end
        end
    end
    # for day 9
        # implement decision point

    # next 4 days to increasing symptoms
    for lag = 8:-1:5
        for cond in conds  # [nil, mild, sick]
            for agegrp in agegrps

                dist_to_new_conditions!(cond, to_nil:to_sick, agegrp, lag, locale, case=case, dat=dat)

            end
        end
    end
    # for day 5
        # implement decision point
    # initial day + 4 days
    for lag = 4:-1:1
        for cond in conds  # [nil, mild, sick]
            for agegrp in agegrps

                dist_to_new_conditions!(cond, to_nil:to_mild, agegrp, lag, locale, case=case, dat=dat)
                # input!(grab(nil:severe,1:5,lag,1),nil:severe,1:5,lag+1,1, dat=dat)
                # minus!(grab(nil:severe,1:5,lag,1),nil:severe,1:5,lag,1, dat=dat)
            end
        end
    end

    update_infectious!(locale, dat = dat)
end


function dist_to_new_conditions!(fromcond, tocases, agegrp, lag, locale; case = "open", dat=openmx)

    # get the number of folks to be distributed 
    folks = grab(fromcond,agegrp,lag,locale, dat=dat) # scalar
    # @debug @sprintf("%12s %5d %5d %5d ",condnames[fromcond], agegrp,  lag, folks)


    # select probs for tocases from the raw evolve probs
    toprobs = zeros(6)
    toprobs[tocases] = snorm(ev_pr[fromcond][tocases, agegrp])

    # get the dist vector of folks to each outcome (6 outcomes)  
        # distvec:   1: recovered 2: nil 3: mild 4: sick 5: severe 6: dead
    distvec = bucket(categorical_sample(toprobs,folks), 6, 6)
    @debug begin
        cond = condnames[fromcond];recovered=distvec[1]; nil=distvec[2]; mild=distvec[3]; sick=distvec[4]; severe=distvec[5]; dead=distvec[6];
        "distribute $cond age $agegrp lag $lag CNT $folks to $recovered $nil $mild $sick $severe $dead" 
        end

    # distribute to infectious cases,recovered, dead
    moveout = sum(distvec)
    plus!(distvec[to_nil:to_severe], infectious_cases, agegrp, lag+1,locale, dat=dat) # add to infectious cases for next lag
    plus!(distvec[to_recovered], recovered, agegrp, 1, locale, dat=dat)  # recovered to lag 1
    plus!(distvec[to_dead], dead, agegrp, 1, locale, dat=dat)  # dead to lag 1
    minus!(moveout, fromcond, agegrp, lag, locale, dat=dat)  # subtract from cond in current lag

    queuestats(distvec, [recovered, nil,mild,sick,severe, dead], locale, evolvestat)
    queuestats(-moveout, fromcond, locale, evolvestat)
end


function queuestats(vals, conds, locale, func::typeof(evolvestat); case="open")
    @assert length(vals) == length(conds) "lengths of vals and indices don't match"
    worldday = simday(0)
    for i in eachindex(vals)
        vals[i] == 0 && continue
        additem = func(cnt = vals[i], locale=locale, tocond=conds[i])
        push!(newstatq, (worldday, additem.cnt, additem.locale, additem.tocond))
        # enqueue!(newstatq,(simday=simday(),additem...))
    end
end


function update_infectious!(locale; dat=openmx) # by single locale
    for agegrp in agegrps
        tot = total!([nil, mild, sick, severe],agegrp,:,locale,dat=dat) # sum across cases and lags per locale and agegroup
        input!(tot, infectious, agegrp, 1, locale, dat=dat) # update the infectious total for the locale and agegroup
    end
end


function lastlag_distribute!(ev_pr, cond, agegrp, lag, locale; dat=openmx)
    condfolks = grab(cond, agegrp,19,locale, dat=dat)

    # not a distribution:   2 final outcomes: either recovered or dead--not probabilistic


    recoverpr = ev_pr[cond][to_recovered, agegrp]
    recoverfolks = ceil(Int, recoverpr * condfolks)
    plus!(recoverfolks, recovered, agegrp, 1, locale, dat=dat)  # recovered to lag 1

    deadpr = ev_pr[cond][to_dead, agegrp]
    deadfolks = ceil(Int, deadpr * condfolks)
    plus!(deadfolks, dead, agegrp, 1, locale, dat=dat)  # recovered to lag 1

    minus!(recoverfolks+deadfolks, cond, agegrp, lag, locale, dat=dat)  # subtract from cond in current lag

    residual = condfolks - recoverfolks - deadfolks
    if residual > 0
        @warn "Uh-oh, some people infectious for more than 18 days " cnt=residual cond=condnames[cond] agegrp
    end
    lag19error(residual)  # accumulator
end


"""
    How far do the infectious people spread the virus to 
    previously unexposed people, by agegrp?
"""
function spread!(locale, dat=openmx)
    # 180 microseconds
    # start with the number of infectious people      
    # now we ignore what their condition and age are is: TODO fix
    # should spread less for conditions in order: nil, mild, sick, severe
    # need to loop by condition

    # how many spreaders  TODO grab their condition.  Separate probs by condition
    spreaders = Int(sum(grab(infectious, agegrps, 1:19, locale, dat=dat)))
    if spreaders == 0
        return nothing
    end

    # how many people are touched   
    touched = how_many_touched(spreaders)

    # contact that causes infection
    # amplify probabilities with population density for locale
    byage = split_by_age(touched)
    for i in 1:length(byage) # this probabilisticly determines if contact resulted in contracting the virus
        byage[i] = binomial_one_sample(byage[i], contact_risk_by_age[i])
    end

    # move the people from unexposed:agegrp to infectious:agegrp and nil
    plus!.(byage, infectious, agegrps, 1, locale, dat=dat)  
    plus!.(byage, nil, agegrps, 1, locale, dat=dat)  
    
    return nothing
end


function split_by_age(cnt)::Array{Int64,1}
    numcells = ages
    probvec = age_dist
    x = categorical_sample(probvec, cnt)
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
    # TODO: should the caller do the loop across locales?   YES
    locales =1:numgeo
    travdests = collect(locales)
    deleteat!(travdests,findfirst(isequal(locale), travdests))
    bins = lim = length(travdests) + 1
    for agegrp in agegrps
        for cond in [unexposed, infectious, recovered]
            name = condnames[cond]
            for lag in lags
                numfolks = sum(grab(cond, agegrp, lag, locale)) # this locale, all lags
                travcnt = floor(Int, gamma_prob(travprobs[agegrp]) * numfolks)  # interpret as fraction of people who will travel
                x = rand(travdests, travcnt)  # randomize across destinations
                bydest = bucket(x, lim, bins)
                for dest in 1:length(bydest)
                    isempty(bydest) && continue
                    cnt = bydest[dest]
                    iszero(cnt) && continue
                    enqueue!(travelq, travitem(cnt, locale, dest, agegrp, lag, name))
                end
            end
        end
    end
end


"""
    Assuming a daily cycle, at the beginning of the day 
    process the queue of travelers from the end of the previous day. 
    Remove groups of travelers by agegrp, lag, and condition
    from where they departed.  Add them to their destination.
"""
function travelin!(dat=openmx)
    while !isempty(travelq)
        g = dequeue!(travelq)
        cond = eval(Symbol(g.cond))
        minus!(g.cnt, cond, g.agegrp, g.lag, g.from, dat=dat)
        plus!(g.cnt, cond, g.agegrp, g.lag, g.to, dat=dat)
    end
end

"""
    Move people into isolation and out of the open movement environment.
"""
function isolate(iso_pr, locales = locales,  rules=[])
    for locale in locales
        iso_pr_locale = get(iso_pr, locale, iso_pr["default"])
        isolate_out!(iso_pr_locale, locale)
    end
    isolate_in!()
end


"""
    Process the isolated queue.
    Remove people from open environment in array openmx.
    Add people to the insolated environment in array isolatedmx.
"""
function isolate_in!()
    while !isempty(isolatedq)
        g = dequeue!(isolatedq)
        minus!(g.cnt,g.cond, g.agegrp, g.lag, g.locale, dat=openmx)
        plus!(g.cnt, g.cond, g.agegrp, g.lag, g.locale, dat=isolatedmx)
    end
end


"""
    Place people who are isolating into the isolated queue: isolatedq
"""
function isolate_out!(iso_pr_locale, locale)
    for (idx, cond) in enumerate([unexposed, recovered, nil, mild, sick, severe])
        for lag in lags
            for agegrp in agegrps
                cnt = binomial_one_sample(grab(cond, agegrp, lag, locale), 
                                        iso_pr_locale[idx, agegrp])
                cnt == 0 && continue
                enqueue!(isolatedq, iso_item(cnt, cond, agegrp, lag, locale))
            end
        end
    end
end

function unisolate(locale, rules=[])  # rules want to be functions
end


#######################################################################################
#  probability
#######################################################################################

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

# range counts to discretize PDF of continuous outcomes
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

"""
    Returns continuous value that represents gamma outcome for a given
    approximate center point (scale value of gamma).  We can interpret this 
    as a funny sort of probability or as a number outcome from a gamma
    distributed sample. 
    1.2 provides a good shape with long tail right and big clump left
"""
function gamma_prob(target; shape=1.2)
    @assert 0.0 <= target <= 99.0 "target must be between 0.0 and 99.0"
    dgamma = Gamma(1.2,target)
    pr = rand(dgamma, 1)[1] / 100.0
end


"""
    Returns a single number of successes or number of members or a group for whom
    the sampled outcome is a success given the probability of success.
"""
function binomial_one_sample(cnt, pr)::Int
    return rand(Binomial(cnt, pr))
end


function snorm(arr)
    return arr./(sum(arr))
end

function categorical_sample(probvec, trials)
    @assert isapprox(sum(probvec), 1.0) "target vector must sum to 1.0"
    x = rand(Categorical(probvec), trials)
end


#########################################################################################
# tracking
#########################################################################################

function accumulator(i)  # input for the first time initializes the counter
               f(n) = i += n
               return f
       end

lag19error = accumulator(0)  # judgment parameters result in how many people still sick after 19 days?


# the day in the world
simday = accumulator(1)  # simday(1) to advance by a day; 

function queue_to_stats(qdf, seriesdf, numgeo)
    size(qdf,1) == 0 && return   # empty queue

    rowinit = Dict(:Unexposed => 0,  :Infectious => 0, :Recovered=> 0, :Dead=> 0, :Nil=> 0, 
        :Mild=> 0, :Sick=> 0, :Severe=> 0,  :Travelers=> 0, :Isolated=> 0)

    thisday = qdf[1, :simday]
    @assert all(qdf[!, :simday] .== thisday)  "queue doesn't contain all the same days"
    agg = by(qdf, [:locale, :tocond], cnt = :cnt => sum)
    println(agg)
    for l in 1:numgeo
        filt = agg[agg.locale .== l, :]
        rowfill = copy(rowinit)
        for r in eachrow(filt)
            rowfill[Symbol(series_colnames[r.tocond])] = r.cnt 
        end
        push!(seriesdf[l][:new], rowfill)
    end
    # purge the queue
    deleterows!(qdf,1:size(qdf,1))
end

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
    dat[lag, condition, agegrp, locale] += val
end


function minus!(val, condition, agegrp, lag, locale; dat=openmx)
    dat[lag, condition, agegrp, locale] -= val
end

# avoid overloading base.sum 
function total!(condition, agegrp, lag, locale; dat=openmx)
    sum(dat[lag, condition, agegrp, locale])
end

function showq(qname)
    for item in qname
        println(item)
    end
end

function printsp(xs...)
    for x in xs
       print(x," ")
    end
   println()
end
