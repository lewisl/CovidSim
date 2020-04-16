
# pre-allocate large arrays, accessed and modified frequently
# hold complex parameter sets
mutable struct Env
    spreaders::Array{Int64,3} # 19,4,5
    all_accessible::Array{Int64,3} # 19,6,5
    numcontacts::Array{Int64,3} # 19,4,5
    simple_accessible::Array{Int64,2} # 6,5
    numtouched::Array{Int64,2} # 19,5
    lag_contacts::Array{Int,1} # 19,
    riskmx::Array{Float64,2} # 19,5
    contact_factors::Array{Float64,2}  # 4,5 parameters for spread!
    touch_factors::Array{Float64,2}  #  6,5  parameters for spread!
    send_risk_by_lag::Array{Float64,1}  # 19,  parameters for spread!
    recv_risk_by_age::Array{Float64,1}  # 5,  parameters for spread!
    sd_compliance::Array{Float64,2} # (6,5) social_distancing compliance unexp,recov,nil:severe by age

    # constructor with keyword arguments and very basic type compatible defaults
    function Env(;          spreaders=zeros(Int,0,0,0),   # semicolon for all keyword (named) arguments)
                            all_accessible=zeros(Int,0,0,0),
                            numcontacts=zeros(Int,0,0,0),
                            simple_accessible=zeros(Int,0,0),
                            numtouched=zeros(Int,0,0),
                            lag_contacts=zeros(Int,19),
                            riskmx=zeros(Float64,0,0),
                            contact_factors=zeros(Float64,0,0),
                            touch_factors=zeros(Float64,0,0),
                            send_risk_by_lag=zeros(Float64,19),
                            recv_risk_by_age=zeros(Float64,5),
                            sd_compliance=ones(Float64,6,5)        )
        return new(spreaders, all_accessible, numcontacts, simple_accessible, 
                   numtouched, lag_contacts, riskmx, contact_factors, 
                   touch_factors, send_risk_by_lag, recv_risk_by_age, sd_compliance)
    end
        
end

#  stash for temporary values changed during simulation cases
#      to change just once and then get the originals back
#      it is the users responsibility to get rid of stuff
const stash = Dict{Symbol, Array}()

# control constants
const age_dist = [0.251, 0.271,   0.255,   0.184,   0.039]
const lags = 1:19   # rows

# geo data: id,fips,county,city,state,sizecat,pop,density
const id = 1
const fips = 2
const county = 3
const city = 4
const state = 5
const sizecat = 6
const popsize =  7
const density = 8
const density_fac = 9


# population centers sizecats  
const major = 1  # 20
const large = 2  # 50
const medium = 3
const small = 4
const smaller = 5
const rural = 6

# condition_outcome columns and stat series columns
const unexposed = 1
const infectious = 2
const recovered = 3
const dead = 4
const nil = 5
const mild = 6
const sick = 7
const severe = 8
const travelers = 9
const isolated = 10

const conditions = [unexposed, infectious, recovered, dead, nil, mild, sick, severe]
const condnames = Dict(1=>"unexposed",2=>"infectious",3=>"recovered", 4=>"dead",
                       5=>"nil",6=>"mild",7=>"sick",8=>"severe")
const infectious_cases = [nil, mild, sick, severe]
const transition_cases = [recovered, nil, mild, sick, severe, dead]

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


# traveling constants
const travprobs = [1.0, 2.0, 3.0, 3.0, 0.4] # by age group



"""
Map a condition index from rows in the data matrix to
indices for the transition probabilities:

```
               unexposed  infectious  recovered  dead   nil  mild  sick  severe
data rows         1          2            3        4     5     6     7     8
transition pr    -1         -1            1        6     2     3     4     5
```

Tranition pr indices that return -1 are not used and will raise an error.

- Use with text literal in code as mapcond2tran.nil => 2
- Use with variables that stand for the data rows as mapcond2tran[nil] => 2
"""
const mapcond2tran = (unexposed=-1, infectious=-1, recovered=1, dead=6, nil=2, mild=3, sick=4, severe=5)


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
        #   DONE = isolate_queue! put people into isolation by agegroup, locale, condition, lag
                    # and isolate_move!=>remove from open, unisolate_queue!, unisolate_move!
                    # add isolation people to queue
                    # unisolate:  remove people from isolatedmx and put them back into openmx
        #  transition people in isolation across conditions
        

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
            # DONE: handle partial (treating as full) immunity of recovered
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


function run_a_sim(geofilename, n_days, locales; runcases=[], spreadcases=[], dtfilename = "dec_tree_all.csv",
                   nsfilename="nodstarts.csv", silent=true)
#=
    see cases.jl for runcases and spreadcases
=#

    !isempty(dayq) && (deleteat!(dayq, 1:length(dayq)))   # empty it

    locales = locales   # force local scope to the loop
    alldict = setup(geofilename; dectreefilename="dec_tree_all.csv",
                    node_starts_filename="dec_tree_starts.csv", geolim=10)
    dt_set = alldict["dt"]  # decision trees for transition
    # get iso_pr here
    openmx = alldict["dat"]["openmx"]
    isolatedmx = alldict["dat"]["isolatedmx"]
    dseries = build_series(locales)   # should use numgeo here
    geodata = alldict["geo"]

    # pre-allocate arrays and initialize parameters for spread!
    env = initialize_sim_env()

    # start the counter at zero
    reset!(ctr, :day)  # remove key :day leftover from prior runs

    starting_unexposed = reduce(hcat, [grab(unexposed, agegrps, 1, i, dat=openmx) for i in locales])
    starting_unexposed = vcat(locales', starting_unexposed)

    # spreadtime = 0
    # trantime = 0

    for i = 1:n_days
        inc!(ctr, :day)  # update the simulation day counter
        silent || println("simulation day: ", ctr[:day])
        for locale in locales
            # seed some people arriving as carriers
            if ctr[:day] == 1
                println("first seed locale $locale....")
                seed!([0,3,3,0,0],5,nil, agegrps, locale, dat=openmx)
                queuestats([0,3,3,0,0], locale, spreadstat; case="open")
            elseif ctr[:day] == 10
                println("second seed locale $locale...")
                seed!([0,3,3,0,0],5,nil, agegrps, locale, dat=openmx)
                queuestats([0,3,3,0,0], locale, spreadstat; case="open")
            end
            for case in runcases
                case(locale, opendat=openmx, isodat=isolatedmx, env=env)
            end
            density_factor = geodata[locale,density_fac]
            spread!(locale, density_factor, dat=openmx, env=env, spreadcases=spreadcases)
            transition!(dt_set, locale, dat=openmx)   # transition all infectious cases "in the open"
            transition!(dt_set, locale, dat=isolatedmx)  # transition all infectious cases in isolation
        end
        queue_to_newseries!(newstatq, dseries, locales) # update tracking stats for the day
    end
    println("Simulation completed for $(ctr[:day]) days.")

    for locale in locales
        new_to_cum!(dseries, locale, starting_unexposed)
    end

    return alldict, dseries, starting_unexposed
end



####################################################################################
#   functions for simulation events
####################################################################################
# things that cause condition changes: travel, spread, transition, isolate


function initialize_sim_env()
    Env(spreaders=zeros(Int64,19,4,5),  
        all_accessible=zeros(Int64,19,6,5),  
        numcontacts=zeros(Int64,19,4,5),  
        simple_accessible=zeros(Int64,6,5),  
        numtouched=zeros(Int64,19,5),   
        lag_contacts=zeros(Int64,19),      
        riskmx = zeros(Float64,19,5),    
        contact_factors =       [ 1    2.5    2.5     1.5    1;    # nil
                                  1    2.5    2.5     1.5    1;    # mild
                                0.7    1.0    1.0     0.7   0.5;   # sick
                                0.5    0.8    0.8     0.5   0.2],  # severe
                      # agegrp    1     2      3       4     5
        touch_factors =         [.6    .9      .8     .6   .35;    # unexposed
                                 .6    .9      .8     .6   .35;    # recovered
                                 .6    .9      .8     .6   .35;    # nil
                                 .6    .9      .7     .5   .28;    # mild
                                 .28   .35     .28    .18  .18;    # sick
                                 .18   .18     .18    .18  .18],   # severe                               
        send_risk_by_lag = [.1, .3, .7, .9, .9, .9, .8, .7, .5, .4, .2, .1, .1, 0.05, 0.05, 0.5, 0, 0, 0],
        recv_risk_by_age = [.1, .3, .3, .4, .5],
        sd_compliance = zeros(6,5))
end


function seed!(cnt, lag, conds, agegrps, locales; dat=openmx)
    @assert length(lag) == 1 "input only one lag value"
    # @warn "Seeding is for testing and may result in case counts out of balance"
    for loc in locales
        for cond in conds
            if cond in [nil, mild, sick, severe]
                input!.(cnt, cond, agegrps, lag, loc, dat=dat)
                minus!.(cnt, unexposed, agegrps, 1, loc, dat=dat)
            elseif cond in [recovered, dead]
                input!.(cnt, cond, agegrps, lag, loc, dat=dat)
                # need to figure out what condition and lag get decremented!!
            end
        end
        update_infectious!(loc, dat = dat)
    end
end


"""
    transition!(dt, locale; case="open", dat=openmx)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability. The final
outcomes are recovered or dead.
"""
function transition!(dt_set, locale; case="open", dat=openmx)  # TODO also need to run for isolatedmx

    nodestarts = dt_set.starts # day that a node takes effect
    dt = dt_set.dt             # decision trees

    for lag = 19:-1:1
        if lag in keys(nodestarts)
            # decision points determine how people's conditions change
            for agegrp in agegrps # get the params from the dec_tree
                for node in nodestarts[lag]
                    toprobs = zeros(6)
                    for branch in dt[agegrp][node]  # agegroup index in array, node key in agegroup dict
                        toprobs[mapcond2tran[branch.tocond]] = branch.pr
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
    if folks == 0
        return   # save some time if there is nothing to do
    end
    @debug "folks $folks lag $lag age $agegrp cond $fromcond"


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

    leaveout = mapcond2tran[fromcond]
    deleteat!(distvec, leaveout) # exclude the source condition
    deleteat!(transition_cases, leaveout)

    moveout = sum(distvec)
    minus!(moveout, fromcond, agegrp, lag, locale, dat=dat)  # subtract from cond in current lag

    queuestats(distvec, transition_cases, locale, transitionstat)
    queuestats(-moveout, fromcond, locale, transitionstat)
end


function update_infectious!(locale; dat=openmx) # by single locale
    for agegrp in agegrps
        tot = total!([nil, mild, sick, severe],agegrp,:,locale,dat=dat) # sum across cases and lags per locale and agegroup
        input!(tot, infectious, agegrp, 1, locale, dat=dat) # update the infectious total for the locale and agegroup
    end
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
function gamma_prob(target; shape=1.0)
    @assert 0.0 <= target <= 99.0 "target must be between 0.0 and 99.0"
    dgamma = Gamma(shape,target)
    pr = rand(dgamma, 1)[1] / 100.0
end


"""
Returns a single number of successes for a
sampled outcome of cnt tries with the input pr of success.
"""
function binomial_one_sample(cnt, pr)::Int
    return rand.(Binomial.(cnt, pr))
end


function categorical_sample(probvec, trials)
    x = rand(Categorical(probvec), trials)
end


function cumplot(dseries, locale, plseries=[:Unexposed,:Infectious,:Recovered, :Dead];
    geo=[])

    pyplot()
    theme(:ggplot2, foreground_color_border =:black, reuse = false)

    !(typeof(plseries) <: Array) && (plseries = [plseries])

    # the data
    n = size(dseries[locale][:cum],1)
    cumseries = Matrix([DataFrame(Day = 1:n) dseries[locale][:cum][!,plseries]])
    labels = string.(plseries)
    labels = reshape([labels...], 1, length(labels))
    people = dseries[locale][:cum][1,:Unexposed] + dseries[locale][:cum][1,:Infectious]
    cityname = !isempty(geo) ? geo[locale, city] : ""
    died = dseries[locale][:cum][end,:Dead]
    infected = dseries[locale][:cum][1,:Unexposed] - dseries[locale][:cum][end,:Unexposed]
    firstseries = plseries[1]
    half_yscale = floor(Int, maximum(dseries[locale][:cum][!,firstseries]) * 0.5)

    # the plot
    plot(   cumseries[:,1], cumseries[:,2:end], 
            size = (700,500),
            label = labels, 
            lw=2.3,
            title = "Covid for $people people in $cityname over $n days\nActive Cases for Each Day",
            xlabel = "Simulation Days",
            yaxis = ("People"),
            legendfontsize = 10,
            reuse = false
        )
    annotate!((6,half_yscale,Plots.text("Died: $died\nInfected: $infected", 10, :left)))

end

function newplot(dseries, locale, plseries=[:Infectious])

    pyplot()
    theme(:ggplot2, foreground_color_border =:black)

    !(typeof(plseries) <: Array) && (plseries = [plseries])

    # the data
    n = size(dseries[locale][:new],1)
    newseries = Matrix([DataFrame(Day = 1:n) dseries[locale][:new][!,plseries]])
    labels = string.(plseries)
    labels = reshape([labels...], 1, length(labels))
    people = dseries[locale][:cum][1,:Unexposed] + dseries[locale][:cum][1,:Infectious]

    # the plot
    bar(    newseries[:,1], newseries[:,2:end], 
            size = (700,500),
            label = labels, 
            lw=0,
            title = "Covid Daily Change for $people people over $n days",
            xlabel = "Simulation Days",
            yaxis = ("People"),
            reuse =false
        )

end


function day2df(dayq::Array)
    dayseries = DataFrame(dayq)

    dayseries[!, :cuminfected] .= zeros(Int, size(dayseries,1))
    dayseries[1, :cuminfected] = copy(dayseries[1,:infected])
    for i = 2:size(dayseries,1)
       dayseries[i,:cuminfected] = dayseries[i-1,:cuminfected] + dayseries[i,:infected]
    end

    return dayseries
end

function dayplot(dayseries::DataFrame)
    plot(dayseries[!,:day], dayseries[!,:spreaders],label="Spreaders", dpi=200,lw=2,
         xlabel="Simulation Days", ylabel="People", title="Daily Spread of Covid",
         bg_legend=:white)
    
    plot!(dayseries[!,:day], dayseries[!,:contacts],label="Contacts", dpi=200,lw=2)
    plot!(dayseries[!,:day], dayseries[!,:touched],label="Touched", dpi=200,lw=2)
    plot!(dayseries[!,:day], dayseries[!,:infected],label="Infected", dpi=200,lw=2)
    plot!(dayseries[!,:day], dayseries[!,:cuminfected],label="Cum Infected", dpi=200,lw=2)

end



function day_animate2(dayseries)
    n = size(dayseries,1)
    # daymat = Matrix(dayseries)

    xd = dayseries[1:5,:]

    topy = max(maximum(dayseries[!,:spreaders]),maximum(dayseries[!,:contacts]),
                maximum(dayseries[!,:touched]),maximum(dayseries[!,:infected]) )

    @df xd plot(:day, [:spreaders :contacts :touched :infected], color=^([:red :blue :green :orange]),
                labels=^(["Spreaders" "Contacts" "Touched" "Infected"]),dpi=200, lw=2,ylim=(0,topy))

    for i = 5:2:n
        xd = dayseries[i-2:i,:]

        @df xd plot!(:day, [:spreaders :contacts :touched :infected], color=^([:red :blue :green :orange]),
                 labels=false, dpi=200, lw=2, ylim=(0,3e4))
        gui()

        if i < round(Int, n/4)
            sleep(0.3)
        elseif i < round(Int,n/2)
            sleep(0.1)
        else
            sleep(.001)
        end
        # print("\nPress enter to continue, q enter to quit.> ");
        # ans = chomp(readline()) 
        # if ans == "q"
        #     break
        # end    
    end
end

# Plots.AnimatedGif("/var/folders/mf/73qj_8c91dzg4sw459_7mchm0000gn/T/jl_Js4px6.gif")

####################################################################################
#   convenience functions for reading and inputting population statistics
#                in the simulation data matrices
####################################################################################


function grab(condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    return dat[locale][lag, condition, agegrp]
end


function input!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    dat[locale][lag, condition, agegrp] = val
end


function plus!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    dat[locale][lag, condition, agegrp] += val
end


function minus!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    dat[locale][lag, condition, agegrp] -= val
end


function total!(condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    sum(dat[locale][lag, condition, agegrp])
end


#############################################################
#  other convenience functions
#############################################################


function printsp(xs...)
    for x in xs
       print(x," ")
    end
   println()
end

sparsify!(x, eps=1e-8) = x[abs.(x) .< eps] .= 0.0;
