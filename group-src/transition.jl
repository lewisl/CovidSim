####################################################
# transition.jl
#     change status of folks in simulation:
#           transition
#           travel
#           seed
#           isolate
####################################################


"""
Map a condition index from rows in the data matrix to
indices for the transition probabilities:

```
               unexposed  infectious  recovered  dead   nil  mild  sick  severe
data rows         1          2            3        4     5     6     7     8
transition pr    -1         -1            1        6     2     3     4     5
```

Transition probability indices that return -1 are not used and will raise an error.

- Use with text literal in code as map2pr.nil => 2
- Use with variables that stand for the data rows as map2pr[nil] => 2
"""
const map2pr = (unexposed=-1, infectious=-1, recovered=1, dead=6, nil=2, mild=3, sick=4, severe=5)

function seed!(day, cnt, lag, conds, agegrps, locale, dat)
    @assert length(lag) == 1 "input only one lag value"
    # @warn "Seeding is for testing and may result in case counts out of balance"
    if day == ctr[:day]
        println("*** seed day $(ctr[:day]) locale $locale....")
        for loc in locale
            for cond in conds
                @assert (cond in [nil, mild, sick, severe]) "Seed cases must have conditions of nil, mild, sick, or severe" 
                plus!(cnt, cond, agegrps, lag, loc, dat)
                minus!(cnt, unexposed, agegrps, 1, loc, dat)
                update_infectious!(loc, dat)
            end
        end
    end
end


# method to run through all existing locales in isolation
function transition!(dt, all_decpoints, dat)
    for locale in keys(dat)
        transition!(dt, all_decpoints, locale, dat)
    end
end


"""
    transition!(dt, all_decpoints, locale, dat)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability; then to 
recovered or dead.

Works for a single locale.
"""
function transition!(dt, all_decpoints, locale, dat)  

    # @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single integer or NamedTuple"
    iszero(dat[locale]) && (return)

    #pre-allocate variables updated in loop
    toprobs = @MVector zeros(Float64, 6)
    distvec = @MVector zeros(T_int[], 6)  #  
    # tree = Dict{Tuple{Int64, Int64}, Array{CovidSim.Branch, 1}}() 
    # age_decpoints = Dict{Int64, Array{Tuple{Int64, Int64}, 1}}()

    @inbounds for agegrp in agegrps
        for lag = laglim:-1:1
            if lag in all_decpoints[agegrp] # check if a decision tree applies to this lag           
                tree = dt[agegrp]    # .tree
                # age_decpoints = all_decpoints[agegrp]
                age_bump = copy(infectious_cases)
                for node in [node for node in keys(tree) if first(node) == lag] # get(age_decpoints, lag, []) # skip the loop is this agegrp doesn't have this decpoint
                    # println(ctr[:day],": node: ", node)
                    age_bump = copy(infectious_cases)
                    toprobs = @MVector zeros(Float64, 6)
                    for branch in tree[node]["branches"]  # agegroup index in array, node key in agegroup dict
                        toprobs[map2pr[branch["tocond"]]] = branch["pr"]
                    end
                    @assert isapprox(sum(toprobs), 1.0, atol=1e-6) """\ntoprobs not equal 1.0, 
                                                                        got $(sum(toprobs)) \nfor $node
                                                                        \nbranches: $(tree[node]["branches"])"""
                    fromcond = node[2]   # .fromcond  # all branches of a node MUST have the same fromcond
                    
                    # age_bump = filter(x->x!=fromcond,age_bump)   # remove fromcond distributed to new condition
                    removeit!(age_bump, fromcond)
                    folks = grab(fromcond,agegrp,lag,locale, dat) 

                    if folks > T_int[](0)
                        distribute_to_new_conditions!(folks, fromcond, toprobs, agegrp, lag, locale, node, distvec, dat)
                    end
                end  
                if !isempty(age_bump)  # bump people in conds that didn't get distributed above
                    bump_up!(age_bump, agegrp, lag, locale, dat) 
                end
            else 
                # when no decision tree, bump up every infected person one day within the same condition
                bump_up!(infectious_cases, agegrps, lag, locale, dat)
            end
        end
    end

    update_infectious!(locale, dat) # total all people who are nil, mild, sick, severe across all lags
    return
end


function removeit!(x,y)  # performs ok for tiny vectors
   kill = 1
   for i in eachindex(x)
       if x[i] == y; kill = i; break; end
   end
   deleteat!(x, kill)
end


"""
function bump-up!

Bump people from one lag to lag + 1 in the same disease condition.
"""
function bump_up!(to_bump, agegrp, lag, locale, dat)
    bump = grab(to_bump, agegrp, lag, locale, dat)

    if sum(bump) > T_int[](0)
        plus!(bump, to_bump, agegrp, lag+1, locale, dat)
        minus!(bump, to_bump, agegrp, lag,   locale, dat)
    end
end


"""
function distribute_to_new_conditions!

Based on decision trees for each age group, at specific decision points (in days), change
people's disease condition, or move them to recovered or dead.
"""
function distribute_to_new_conditions!(folks, fromcond, toprobs, agegrp, lag, locale, node, distvec, dat, lastlag=laglim)
    
    @inbounds begin
        @debug "day $(ctr[:day])  folks $folks lag $lag age $agegrp cond $fromcond"

        # set vector of folks to each outcome (6 outcomes): 1: recovered 2: nil 3: mild 4: sick 5: severe 6: dead
        @assert isapprox(sum(toprobs), 1.0, atol=1e-4) "target vector must sum to 1.0; submitted $toprobs"
        x = categorical_sample(toprobs, folks)  # integer results

        # distvec = static_bucket(x, vals=1:length(toprobs))   # toprobs ALWAYS = 6     # [count(x .== i) for i in 1:size(toprobs,1)]

        distvec[:] = [count(x .== i) for i in 1:6]  # length(toprobs)

        # @assert sum(distvec) == folks "someone got lost $res != $folks"

        if lag != lastlag  # infectious cases to next lag
            @views plus!(distvec[map2pr.nil:map2pr.severe], infectious_cases, agegrp, lag+1, locale, dat) # @views 
        end
        @views plus!(distvec[map2pr.recovered], recovered, agegrp, 1, locale, dat)  # recovered to lag 1
        @views plus!(distvec[map2pr.dead], dead, agegrp, 1, locale, dat)  # dead to lag 1
        @views minus!(folks, fromcond, agegrp, lag, locale, dat)  # subtract what we moved from the current lag

        push!(transq, (day=ctr[:day], lag=lag, agegrp=agegrp, fromcond=fromcond,  # @views primarily for debugging; can do some cool plots
                   newcond=distvec[map2pr.nil:map2pr.severe], recovered=distvec[map2pr.recovered],
                   dead=distvec[map2pr.dead], node=node, locale=locale))
    end
    return
end


function update_infectious!(locale, dat) # by single locale
    for agegrp in agegrps
        tot = sum(grab([nil, mild, sick, severe],agegrp,:,locale,dat)) # sum across cases and lags per locale and agegroup
        input!(tot, infectious, agegrp, 1, locale, dat) # update the infectious total for the locale and agegroup
    end
end


"""
For a locale, randomly choose the number of people from each agegroup with
condition of {unexposed, infectious, recovered} who travel to each
other locale. Add to the travelq.
"""
function travelout!(fromloc, locales, rules=[])    # TODO THIS WON'T WORK ANY MORE!
    # 10.5 microseconds for 5 locales
    # choose distribution of people traveling by age and condition:
        # unexposed, infectious, recovered -> ignore lag for now
    # TODO: more frequent travel to and from Major and Large cities
    # TODO: should the caller do the loop across locales?   YES
    travdests = collect(locales)
    deleteat!(travdests,findfirst(isequal(fromloc), travdests))
    bins = lim = length(travdests) + 1
    for agegrp in agegrps
        for cond in [unexposed, infectious, recovered]
            name = condnames[cond]
            for lag in lags
                numfolks = sum(grab(cond, agegrp, lag, fromloc)) # the from locale, all lags
                travcnt = floor(Int, gamma_prob(travprobs[agegrp]) * numfolks)  # interpret as fraction of people who will travel
                x = rand(travdests, travcnt)  # randomize across destinations
                bydest = bucket(x, vals=1:length(travdests))
                for dest in 1:length(bydest)
                    isempty(bydest) && continue
                    cnt = bydest[dest]
                    iszero(cnt) && continue
                    enqueue!(travelq, travitem(cnt, fromloc, dest, agegrp, lag, name))
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
Place people into isolation.

You can enter a percentage (as a fraction in [0.0, 1.0]), an array
of percentages, a number, or an array of numbers.

Use a dot after the function name to apply the same pct or number
input to several conditions, agegrps, lags, or locales.

Use a dot after the function name to apply an array: one or more
of agegrp, cond, or locale must have the same number of elements as the input.
"""
function isolate!(pct::Float64,cond,agegrp,lag,locale, opendat, isodat)
    for c in cond
        for age in agegrp
            for l in lag
                isolate_by!(pct::Float64,c,age,l,locale, opendat, isodat)
            end
        end
    end
end


function isolate_by!(pct::Float64,cond,agegrp,lag,locale, opendat, isodat)
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, lag, locale, opendat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, T_int[](0), T_int[](available))  # limit to max
    cnt < scnt && (@warn "Attempt to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale, opendat, isodat)
end


function isolate_by!(num, cond, agegrp, lag, locale, opendat, isodat)
    @assert sum(num) >= 0 "num must be greater than zero"
    if typeof(locale) <: Quar_Loc
        available = grab(cond, agegrp, lag, locale.locale, opendat)  # max
    else
        available = grab(cond, agegrp, lag, locale, opendat)  # max
    end
    cnt = clamp.(num, T_int[](0), T_int[](available))  # limit to max
    sum(cnt) < sum(num) && (@warn "Attempt to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale, opendat, isodat)
    return nothing
end  # this one works


function _isolate!(cnt, cond, agegrp, lag, locale::Integer, opendat, isodat)
    minus!(cnt, cond, agegrp, lag, locale, opendat)  # move out 
    update_infectious!(locale, opendat)
    plus!(cnt, cond, agegrp, lag, locale, isodat)  # move in
    update_infectious!(locale, isodat)
    return nothing  # this one works!
end

# for test and trace or any isolate that records the day of isolation
function _isolate!(cnt, cond, agegrp, lag, qloc::Quar_Loc, opendat, isodat)
    minus!(cnt, cond, agegrp, lag, qloc.locale, opendat)  # move out 
    update_infectious!(qloc.locale, opendat)
    plus!(cnt, cond, agegrp, lag, qloc, isodat)  # move in
    update_infectious!(qloc, isodat)
    return nothing  # this one works!
end

function unisolate!(pct::Float64,cond,agegrp,lag,locale, opendat, isodat)
    for c in cond
        for age in agegrp
            for l in lag
                unisolate_by!(pct::Float64,c,age,l,locale, opendat, isodat)
            end
        end
    end
end


function unisolate_by!(pct::Float64,cond,agegrp,lag,locale, opendat, isodat)
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, lag, locale, isodat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, T_int[](0), T_int[](available))  # limit to max
    cnt < scnt && (@warn "Attempt to unisolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale, opendat, isodat)
    return nothing  # this one works!
end


function unisolate_by!(num, cond, agegrp, lag, locale, opendat, isodat, mode=:both)
    @assert sum(num) >= 0 "sum(num) must be greater than zero"

    available = grab(cond, agegrp, lag, locale, isodat)  # max

    # println("day $(ctr[:day]) request to unisolate   ", sum(num))
    # println("day $(ctr[:day]) available to unisolate ", sum(available))

    cnt = clamp.(num, T_int[](0), T_int[](available))  # limit to max
    sum(cnt) < sum(num) && (@warn "Attempt to unisolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale,  opendat, isodat, :both)
    return nothing
end  # this one works


function _unisolate!(cnt, cond, agegrp, lag, locale,  opendat, isodat, mode=:both)
    if mode != :plus  # this is when mode = :minus or :both
        minus!(cnt, cond, agegrp, lag, locale, isodat)
        update_infectious!(locale, isodat)
    end
    if mode != :minus  # this is when mode = :plus or :both
        if typeof(locale) <: Quar_Loc
            locale = locale.locale
        end

        # println("day $(ctr[:day])  unquarantine is unisolating this many ", sum(cnt))

        plus!(cnt, cond, agegrp, lag, locale, opendat)
        update_infectious!(locale, opendat)
    end
    return 
end