# TODO
# more info
    # get fatality rate by age and co-morbidity CDC, Italian NIH
    # by agegroup, hospitalization %, ICU admission %, fatality %
    # UW virology, expansion of deaths by state on log chart
# do we really need hist matrices that copy everything?
# probably need to raise effective death rate for people in 60-80 and 80+  per Vox article


using DelimitedFiles
using DataStructures
using DataFrames
using Random
using Distributions
using StatsBase
using Printf
using PyPlot
import Seaborn

include("dec_tree.jl")
include("setup.jl")



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
const major = 1  # 20
const large = 2  # 50
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
const transition_cases = [recovered, nil, mild, sick, severe, dead]
const series_colnames = Dict( 1=>:Unexposed,  2=>:Infectious, 3=>:Recovered, 4=>:Dead, 5=>:Nil, 6=>:Mild, 7=>:Sick,
        8=>:Severe,  9=>:Travelers, 10=>:Isolated)

# transition_prob_rows
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
const agegrps = 1:5  
const ages = length(agegrps)
const contact_risk_by_age = [.05, .1, .15, .20, .25]

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

    const newstatq = DataFrame(day=Int[], cnt=Int[], locale=Int[], tocond=Int[])

    TravelStat = typeof((;cnt=0, locale=0, cond=0, to=0))
    function travelstat(;cnt=0, locale=0, cond=0, to=0)
        (cnt=cnt, locale=locale, cond=cond, to=to)
    end

    TransitionStat = typeof((;cnt=0, locale=0, tocond=5))
    function transitionstat(; cnt=0, locale=0, tocond=5)
        (cnt=cnt, locale=locale, tocond=tocond)
    end

    IsolateStat = typeof((;cnt=0, locale=0))
    function isolatestat(; cnt=0, locale=0)
        (cnt=cnt, locale=locale)
    end

    SpreadStat = typeof((;cnt=0, locale=0, tocond=nil))
    function spreadstat(;cnt=0,locale=0,tocond=nil)
       (cnt=cnt, locale=locale, tocond=tocond)
    end

    # for some reason empty!(df) is deprecated
    empty!(df::DataFrame) = select!(df, Int[])  # select only a null column in place effectively empties


function seed!(cnt, lag, conds, agegrps, locales; dat=openmx)
    @assert length(lag) == 1 "input only one lag value"    
    # @warn "Seeding is for testing and may result in case counts out of balance"
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


"""
Map a condition index from rows in the data matrix to 
indices for the transition probabilities:

```
               unexposed  infectious  recovered  dead   nil  mild  sick  severe
matrix            1          2            3        4     5     6     7     8
transition pr    -1         -1            1        6     2     3     4     5
```

Tranition pr indices that return -1 are not used and will raise an error.
"""
function cond2tran_idx(cond)
    idx =   if cond == recovered  # value 3
                1
            elseif cond == nil    # value 5
                2
            elseif cond == mild   # value 6
                3
            elseif cond == sick   # value 7
                4
            elseif cond == severe # value 8
                5
            elseif cond == dead   # value 4
                6
            else
                -1                # for values 1, 2
            end
end
# better way to represent mapping for this example
    # cmap = (recovered=1, nil=2, mild=3, sick=4, severe=5, dead=6)
    # use as cmap.recovered
    # same speed:  named tuple easier to create as needed
    # can name them for specific purposes but that could be more confusing



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
        #   DONE = isolate_in put people into isolation by agegroup, locale, condition, lag
                    # and isolate_out=>remove from open   
        #  transition people in isolation across conditions
        # add isolation people to queue
        # MUST have unisolate:  remove people from isolatedmx and put them back into openmx

    # DONE = transition! distribute residents across all lags and conditions
        # openmx
        # isolatedmx
        # loop across all locales to transition all    
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
            # DONE: handle partial immunity of recovered
            # DONE = spread!  touchers -> num_touched -> split by age/condition
            # DONE in spread!  of those contacted, a distribution become exposed & nil at lag 0

            # DONE log new infectious


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
#   simulation runner
####################################################################################


function run_a_sim(geofilename, n_days, locales = [1]; silent=true)

    locales = locales   # force local scope to the loop
    simdict = setup(geofilename)
    dt = simdict["dt"]  # decision trees for transition
    # get iso_pr here
    openmx = simdict["dat"]["openmx"]
    # get isolatedmx here  
    dseries = build_series(locales)   # should use numgeo here

    # start the counter at zero
    reset!(ctr, :day)  # remove key :day leftover from prior runs


    starting_unexposed = reduce(hcat, [grab(unexposed, agegrps, 1, i, dat=openmx) for i in locales])
    starting_unexposed = vcat(locales', starting_unexposed)

    for i = 1:n_days
        inc!(ctr, :day)  # update the simulation day counter
        silent || println("simulation day: ", ctr[:day])
        for locale in locales
            # seed some people arriving as carriers
            if ctr[:day] == 1
                println("first seed locale $locale....")
                seed!([0,3,3,0,0],1,nil, agegrps, locale, dat=openmx)
                queuestats([0,3,3,0,0], locale, spreadstat; case="open")
            elseif ctr[:day] == 10
                println("second seed locale $locale...")
                seed!([0,3,3,0,0],1,nil, agegrps, locale, dat=openmx)
                queuestats([0,3,3,0,0], locale, spreadstat; case="open")
            end
            spread!(locale, dat=openmx)
            transition!(dt, locale, dat=openmx)
        end
        queue_to_newseries!(newstatq, dseries, locales)
    end
    println("Simulation completed for $(ctr[:day]) days.")

    for locale in locales
        new_to_cum!(dseries, locale, starting_unexposed)
    end
    return openmx, dseries, starting_unexposed
end


####################################################################################
#   functions for simulation events
####################################################################################
# things that cause condition changes: travel, spread, transition, isolate



"""
    transition!(dt, locale; case="open", dat=openmx)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability. The final 
outcomes are recovered or dead.
"""
function transition!(dt, locale; case="open", dat=openmx)  # TODO also need to run for isolatedmx

    decision_points = [5,9,14,19]               # in lag days
    nodes = Dict(19 => [(4,4)],                 # nodes in decision trees for each agegrp 
                 14 => [(3,2), (3,3), (3,4)],
                  9 => [(2,2), (2,3)],
                  5 => [(1,1)]
                 )

    for lag = 19:-1:1
        if lag in decision_points 
            # decision points determine how people's conditions change 
            for agegrp in agegrps # get the params from the dec_tree
                for node in nodes[lag]
                    toprobs = zeros(6)
                    for branch in dt[agegrp][node]  # agegroup index in array, node key in agegroup dict
                        toprobs[cond2tran_idx(branch.tocond)] = branch.pr
                    end
                    fromcond = dt[agegrp][node][1].fromcond
                    @debug @sprintf("%12s %3f %3f %3f %3f %3f %3f",condnames[fromcond], toprobs...)
                    dist_to_new_conditions!(fromcond, toprobs, agegrp, lag, locale, case=case, dat=dat)
                end
            end
        else
            # bump people up a day without changing their conditions
            input!(grab(nil:severe, agegrps,lag,locale, dat=dat),nil:severe, agegrps,lag+1, locale, dat=dat)
            minus!(grab(nil:severe, agegrps,lag,locale, dat=dat),nil:severe, agegrps,lag, locale, dat=dat)
        end
    end

    # total of all people who are nil, mild, sick, or severe across all lag days
    update_infectious!(locale, dat = dat)  
end


function dist_to_new_conditions!(fromcond, toprobs, agegrp, lag, locale; case = "open", dat=openmx, lastlag=19)
    @assert length(locale) == 1  "Assertion failed: length locale was not 1"
    transition_cases = [recovered, nil,mild,sick,severe, dead]

    # get the number of folks to be distributed 
    folks = grab(fromcond,agegrp,lag,locale, dat=dat) # scalar
    # @debug "folks $folks lag $lag age $agegrp cond $fromcond"
        # println("folks $folks lag $lag age $agegrp cond $fromcond")


    # get the dist vector of folks to each outcome (6 outcomes)  
        # distvec:   1: recovered 2: nil 3: mild 4: sick 5: severe 6: dead
    @assert isapprox(sum(toprobs), 1.0, atol=1e-3) "target vector must sum to 1.0; submitted $toprobs"
    distvec = bucket(categorical_sample(toprobs,folks), lim=6, bins=6)
    @debug begin
        cond = condnames[fromcond];rec=distvec[1]; ni=distvec[2]; mi=distvec[3]; si=distvec[4]; se=distvec[5]; de=distvec[6];
        "distribute $cond age $agegrp lag $lag CNT $folks to $rec $nil $mi $si $se $dead" 
        end


    # distribute to infectious cases,recovered, dead
    lag != lastlag && (plus!(distvec[to_nil:to_severe], infectious_cases, agegrp, lag+1,locale, dat=dat)) # add to infectious cases for next lag
    plus!(distvec[to_recovered], recovered, agegrp, 1, locale, dat=dat)  # recovered to lag 1
    plus!(distvec[to_dead], dead, agegrp, 1, locale, dat=dat)  # dead to lag 1

    leaveout = cond2tran_idx(fromcond)
    deleteat!(distvec, leaveout) # exclude the source condition 
    deleteat!(transition_cases, leaveout)

    moveout = sum(distvec)
    minus!(moveout, fromcond, agegrp, lag, locale, dat=dat)  # subtract from cond in current lag

    queuestats(distvec, transition_cases, locale, transitionstat)
    queuestats(-moveout, fromcond, locale, transitionstat)
end



# TODO this has to work for case = "isolated"
# for transition
function queuestats(vals, conds, locale, func::typeof(transitionstat); case="open")
    @assert length(vals) == length(conds) "lengths of vals and conds don't match"
    thisday = ctr[:day]
    for i in eachindex(vals)
        vals[i] == 0 && continue

        additem = func(cnt = vals[i], locale=locale, tocond=conds[i])
        push!(newstatq, (day=thisday, additem...))
    end
end


# for spread
function queuestats(vals, locale, func::typeof(spreadstat); case="open")
    @assert length(vals) >= 0 "lengths of vals must be greater than 0"
    thisday = ctr[:day]
    cnt = sum(vals)
    if cnt == 0
    else
        # net addition to nil
        additem = func(cnt = cnt, locale=locale)  # defaults to nil
        push!(newstatq, (day=thisday, additem...))
        # net subtraction from unexposed
        additem = func(cnt = -cnt, tocond=unexposed, locale = locale)
        push!(newstatq, (day=thisday, additem...))
    end
end


function update_infectious!(locale; dat=openmx) # by single locale
    for agegrp in agegrps
        tot = total!([nil, mild, sick, severe],agegrp,:,locale,dat=dat) # sum across cases and lags per locale and agegroup
        input!(tot, infectious, agegrp, 1, locale, dat=dat) # update the infectious total for the locale and agegroup
    end
end


"""
How far do the infectious people spread the virus to 
previously unexposed people, by agegrp?  For a single locale...
"""
function spread!(locale; dat=openmx)

    # how many spreaders  TODO grab their condition.  Separate probs by condition
    spreaders = grab(infectious_cases, agegrps, lags, locale, dat=dat) # 19 x 4 x 5 lag x cond x agegrp

    println(" all the spreaders ", sum(spreaders))

    # spreaders = permutedims(spreaders,[1,3,2]) # 19 x 5 x 4: lag x agegrp x cond  Not sure we need this
    if sum(spreaders) == 0
        return
    end

    all_accessible = grab([unexposed,recovered, nil, mild, sick, severe],agegrps,lags, locale, dat=dat)  #   19 x 6 x 5  lag x cond by agegrp
    all_unexposed = grab(unexposed, agegrps, 1, locale, dat=dat)  # 5,  agegrp for lag 1

    # how many people are contacted based on characteristics of spreader   
    numcontacted = how_many_contacted(spreaders, contact_factors)

    println("  all the contacted ", sum(numcontacted))

    # how many people are touched based on characteristics of recipient and potential contacts?
    numtouched = how_many_touched(numcontacted, all_accessible, all_unexposed)   # cond x agegrp

    println("  all the touched ", sum(numtouched))


    newinfected = how_many_infected(numtouched)    # x agegrp (only condition is nil, assumed lag = 1)

    println("  all the newly infected ", sum(newinfected))



    # a not ok way to estimate R0:  need each subsequent day of infected for the originating group; this is only 1 day
    # starters = sum(spreaders)
    # newlyinfected = sum(byage)
    # onestep_r0 = newlyinfected / starters
    # @debug "Spreaders $starters to newly infected: $newlyinfected for r0: $onestep_r0"

    # move the people from unexposed:agegrp to infectious:agegrp and nil
    lag = 1
    plus!.(newinfected, infectious, agegrps, lag, locale, dat=dat)  
    plus!.(newinfected, nil, agegrps, lag, locale, dat=dat)  
    minus!.(newinfected, unexposed, agegrps, lag, locale, dat=dat)

    # add to stats queue for today
    queuestats(sum(newinfected), locale, spreadstat) # sum(5 agegroups), nil is the default, single locale
    
    return 
end

# contact factors for the spreaders

#                         # nil mild sick severe
# const contact_factors = [  2   2    1     .5;      # agegrp 1
#                            4   4    2     1.0;     # agegrp 2
#                            4   4    2     1.0;     # agegrp 3
#                            3   3    1     0.5;     # agegrp 4
#                            2   2    1     0.2;     # agegrp 5
#                         ]

            # agegrp     1   2  3   4    5
const contact_factors = [2   4  4   3    2;  # nil
                         2   4  4   3    2;  # mild
                         1   2  2   1    1;  # sick
                         0.5 1  1  0.5  0.2  # severe
                        ]
# TODO calculate density_factor in setup, per locale
# TODO fix logistic shift and scale
test_density = rand((5000:3_000_000),20)  # use US Census data
function minmax(x) 
    x_max = maximum(x, dims=1)
    x_min = minimum(x, dims=1)
    minmax_density = (x .- x_min) ./ (x_max .- x_min .+ 1e-08)
end
scale_minmax(x, newmin, newmax) = x .* (newmax - newmin) .+ newmin 

"""
Based only on the characteristics of the spreaders, how many of the accessible
are potentially contacted?
"""
function how_many_contacted(spreaders, contact_factors, density_factor=1.1; scale=6)
    #=  This originally ignores the conditions of the touched--assumes they are all equally likely to be touched
        how_many_touched corrects this.
        We assume spreaders is small compared to all_accessible. At some point this might not be true:
        how_many_touched also handles this.
    =#
    sp_lags, sp_conds, sp_ages = size(spreaders)
    numcontacted = zeros(Int, sp_lags, sp_conds, sp_ages)  # 19 x 4 x 5 lag x cond x agegrp

    # SPREADERS: how many people are contacted by each spreader?
        # numcontacted is the potential number of people contacted by a spreader in each 
        # cell by lag (19), infectious cond (4), and agegrp(5)
    for agegrp in sp_ages
        for cond in sp_conds
            for lag in sp_lags
                scale = contact_factors[cond, agegrp]
                dgamma = Gamma(1.2, density_factor * scale)  #shape, scale
                println("dgamma ", dgamma)
                x = rand(dgamma,spreaders[lag, cond, agegrp]);
                println("x ", x)
                if isempty(x)
                else
                    nums, bounds = histo(x)
                    numcontacted[lag, cond, agegrp] = sum(nums .* bounds)  # Float64
                end
            end
        end
    end

    # numcontacted = clamp(floor(Int,numtouched), 0, floor(Int, .8 * tot_unexposed))

    # @assert numcontacted < tot_unexposed "Assertion Error: number touched $numtouched exceeds total unexposed $tot_unexposed"
    return numcontacted
end


"""
For all the spreaders and given how many accessible are potentially contacted, 
how many of the accessible and susceptible actually are "touched" by a spreader? 
"""
function how_many_touched(numcontacted, all_accessible, all_unexposed)  
#= 
    - who is accessible: start with all to capture effect of "herd immunity", when it arises
    - all_accessible 19 x 6 x 5  lag x cond x agegrp, includes unexposed, recovered, nil, mild, sick, severe
    - numcontacted 19 x 4 x 5  lag x cond x agegrp, includes nil, mild, sick, severe = the infectious

    There is a lot of setup before we get to business here.
=#

    # parameters for accessibility of the accessible--not who gets sick, just who is touched!
    mapx = (unexposed= 1, infectious=-1, recovered= 2, dead=-1, nil= 3, mild=  4, sick= 5, severe= 6)
    low, low2, low3 = .6, .5, .3
    lowest, lowest2, lowest3, vlowest = .28, .24, .18, .15
    highest, high, high2 = .9, .8, .7

    access_table = zeros(count(x->x>0, mapx),length(agegrps))  #     6 x 5    cond x agegrp
    access_table[mapx.unexposed, :] .= [low, highest, high, low, lowest]
    access_table[mapx.recovered, :] .= [low, highest, high, low, lowest]


    # not varying access pct for lag of infectious states because we are ignoring viral load by repeat exposures for now
    access_table[mapx.nil, :] .= [low, highest, high, low, lowest]
    access_table[mapx.mild, :] .= [low, highest, high2, low2, lowest2]
    access_table[mapx.sick, :] .= [lowest2, lowest, lowest2, lowest3, lowest3]
    access_table[mapx.severe,:] .= vlowest

    compress_numcontacted = sum(numcontacted,dims=1)[1,:,:] # sum by lag result: 4 x 5  => conds are nil, mild, sick, severe
    totcontacted = sum(compress_numcontacted)
    totaccessible = sum(all_accessible)

    # simplify accessible to unexposed, infectious, recovered by agegrps
    simple_accessible = sum(all_accessible, dims=1)[1,:,:]
    simple_accessible = vcat(simple_accessible[1:2, :], sum(simple_accessible[3:6,:],dims=1)) # 3 x 5
    s_a_pop_pct = reshape(simple_accessible ./ totaccessible, 15)


    # who gets touched in unexposed by agegrp?   
        # we'll normalize down and contact 80% this day
    #=
        - folks in  each cell of numcontacted touch a sample of the accessible reduced by the accessibility factor
        - we can draw a categorical sample for each cell to distribute them across the contact categories
        - we consider the folks who get contacted in infectious and recovered, because that happens--reduce 
           the number of unexposed who are touched
        - we throw out the folks in infectious and recovered and keep those in unexposed
        - use the accessibility factor in a binomial sample:  was the contact successful?
        - we should care about not touching more than the number of accessible but we really only care about unexposed
        - we'll test to be sure we don't exceed the unexposed and reduce touches to 80% of unexposed by agegrp
    =#

    mapi = (unexposed= -1, infectious=-1, recovered= -1, dead=-1, nil= 1, mild=  2, sick= 3, severe= 4)
                    # map to infectious conditions as 1,2,3,4
    dcat = Categorical(s_a_pop_pct)
    touched = zeros(length(agegrps))
    # loop over compress_numcontacted
    for agegrp in agegrps
        for cond in [nil, mild, sick, severe]
            # draw sample across 3 conds and 5 agegrps and only keep the unexposed for 5 age groups
            peeps = reshape(bucket(rand(dcat, compress_numcontacted[mapi[cond], agegrp]), lim=15, bins=15), 3,5)[1,:]
            for a in agegrps
                cnt = binomial_one_sample(peeps[a], access_table[mapx.unexposed, a])
                cnt = cnt < all_unexposed[a] ? cnt : .8 * all_unexposed[a]
                touched[a] += cnt
            end
        end
    end
    
    return touched
end
    

function how_many_infected(numtouched)
    # transmissibility by agegrp of recipient
    newinfected = zeros(size(numtouched,1))  # agegrps
    for i in eachindex(numtouched) # this probabilisticly determines if contact resulted in contracting the virus
        newinfected[i] = binomial_one_sample(numtouched[i], contact_risk_by_age[i])
    end
    return newinfected
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
Add people to the isolated environment in array isolatedmx.
"""
function isolate_move!()
    while !isempty(isolatedq)
        g = dequeue!(isolatedq)
        minus!(g.cnt,g.cond, g.agegrp, g.lag, g.locale, dat=openmx)
        plus!(g.cnt, g.cond, g.agegrp, g.lag, g.locale, dat=isolatedmx)
    end
end


"""
Place people who are isolating into the isolated queue: isolatedq
"""
function isolate_queue!(iso_pr_locale, locale)
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
function bucket(x; lim=5, bins = 5)
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

# unsafe
function snorm(arr)
    return arr ./ (sum(arr))
end

function categorical_sample(probvec, trials)
    x = rand(Categorical(probvec), trials)
end


#########################################################################################
# tracking
#########################################################################################

"""
- use incr!(ctr, :day) for day of the simulation:  creates and adds 1
- use reset!(ctr, :day) to remove :day and return its current value, set it to 0
- use ctr[:day] to return current value of day
"""
const ctr = counter(Symbol) # from package DataStructures

# all locales
function queue_to_newseries!(newstatq, dseries, locales)

    size(newstatq,1) == 0 && return   # empty queue

    thisday = newstatq[1, :day]
    @assert all(newstatq[!, :day] .== thisday)  "Assertion failure: queue doesn't contain all the same days"

    @debug "Size of queue: $(size(newstatq,1))"

    # rowinit = Dict(:Unexposed => 0,  :Infectious => 0, :Recovered=> 0, :Dead=> 0, :Nil=> 0, 
        # :Mild=> 0, :Sick=> 0, :Severe=> 0,  :Travelers=> 0, :Isolated=> 0)

    agg = by(newstatq, [:locale, :tocond], cnt = :cnt => sum)
    for l in locales
        filt = agg[agg.locale .== l, :]
        # rowfill = copy(rowinit)
        rowinit = zeros(10)
        for r in eachrow(filt)
            # rowfill[series_colnames[r.tocond]] = r.cnt 
            rowinit[r.tocond] = r.cnt
        end
        # rowfill[:Infectious] = rowfill[:Nil] + rowfill[:Mild] + rowfill[:Sick] + rowfill[:Severe]
        rowinit[infectious] = sum(rowinit[nil:severe])
        # push!(dseries[l][:new], rowfill)

        push!(dseries[l][:new], rowinit)
    end
    # purge the queue
    deleterows!(newstatq,1:size(newstatq,1))
end

# one locale at a time
function new_to_cum!(dseries, locale, starting_unexposed)
     # add starting unexposed

     println("Updating cumulative statistics for locale $locale.")

     locale_idx = findall(isequal(locale),starting_unexposed[1,:])
     startunexp = sum(starting_unexposed[2:6, locale_idx])
     newseries = dseries[locale][:new]
     cumseries = dseries[locale][:cum]

    # cum first row
    r1_new = collect(newseries[1,:])
    r1_cum = zeros(10) .+ r1_new
    r1_cum[1] += startunexp
    # r1_cum[2] += sum(r1_new[nil:severe])

    push!(cumseries, r1_cum)

    for r_idx in 2:size(newseries,1)
        r1 = collect(newseries[r_idx,:])
        r0 = collect(cumseries[r_idx-1,:])
        push!(cumseries, r0 .+ r1)
    end
end


function cumplot(dseries, locale; sb=false)
    sb && Seaborn.set()

        cumseries = dseries[locale][:cum]
        pldat = Matrix(cumseries[!,[:Unexposed,:Infectious,:Recovered, :Dead]]);
        labels = ["Unexposed", "Infectious","Recovered", "Dead"];
        n = size(cumseries,1)
        people = cumseries[1,:Unexposed] + cumseries[1,:Infectious]

        figure()
        plot(1:n,pldat)
        legend(labels)
        title("Covid for $people people over $n days")
  
end

function newplot(dseries, locale, item; sb=false)
    sb && Seaborn.set()

    if !(item in([:Unexposed,:Infectious,:Recovered, :Dead, :Nil, :Mild, :Sick, :Severe]))
        error("item must be one of :Unexposed,:Infectious,:Recovered, :Dead, :Nil, :Mild, :Sick, :Severe")
    end
    newseries = dseries[locale][:new]
    figure()
    pldat = Vector(newseries[!,item])
    labels = [string(item)]
    n = size(newseries,1)
    # people = cumseries[1,:Unexposed] + cumseries[1,:Infectious]

    bar(1:n,pldat)
    legend(labels)
    title("Covid new cases of $item over $n days")

end


####################################################################################
#   convenience functions for reading and inputting population statistics
#       these all work with one locale at a time EXCEPT grab
####################################################################################


function grab(condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "Assertion Error: locale must be a scalar"
    return dat[locale][lag, condition, agegrp]
end


function input!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "Assertion Error: locale must be a scalar"
    dat[locale][lag, condition, agegrp] = val
end


function plus!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "Assertion Error: locale must be a scalar"
    dat[locale][lag, condition, agegrp] += val
end


function minus!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "Assertion Error: locale must be a scalar"
    dat[locale][lag, condition, agegrp] -= val
end

# avoid overloading base.sum 
function total!(condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "Assertion Error: locale must be a scalar"
    sum(dat[locale][lag, condition, agegrp])
end


#############################################################
#  other convenience functions
#############################################################


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

sparsify!(x, eps=1e-8) = x[abs.(x) .< eps] .= 0.0;
