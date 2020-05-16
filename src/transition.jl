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



function seed!(day, cnt, lag, conds, agegrps, locale; dat=openmx)
    @assert length(lag) == 1 "input only one lag value"
    # @warn "Seeding is for testing and may result in case counts out of balance"
    if day == ctr[:day]
        println("*** seed day $(ctr[:day]) locale $locale....")
        for loc in locale
            for cond in conds
                @assert (cond in [nil, mild, sick, severe]) "Seed cases must have conditions of nil, mild, sick, or severe" 
                input!.(cnt, cond, agegrps, lag, loc, dat=dat)
                minus!.(cnt, unexposed, agegrps, 1, loc, dat=dat)
                update_infectious!(loc, dat = dat)
            end
        end
    end
end


# method to run through all existing locales in isolation
function transition!(dt; dat)
    for locale in keys(dat)
        transition!(dt, locale, dat=dat)
    end
end


"""
    transition!(dt, locale; case="open", dat=openmx)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability; then to 
recovered or dead.

Works for a single locale.
"""
function transition!(dt, locale; dat)  

    # @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single integer or NamedTuple"
    iszero(dat[locale]) && (return)

    all_decpoints = reduce(merge,dt[agegrp].dec_points for agegrp in agegrps)
    toprobs = zeros(6)
    for lag = laglim:-1:1
        if lag in keys(all_decpoints) # check if a decision tree applies to this lag
            for agegrp in agegrps
                tree = dt[agegrp].tree
                age_decpoints = dt[agegrp].dec_points
                age_bump = copy(infectious_cases)
                for node in get(age_decpoints, lag, []) # skip the loop is this agegrp doesn't have this decpoint
                    toprobs[:] = zeros(6)
                    for branch in tree[node]  # agegroup index in array, node key in agegroup dict
                        toprobs[map2pr[branch.tocond]] = branch.pr
                    end
                    @assert isapprox(sum(toprobs), 1.0, atol=1e-6) "toprobs not equal 1.0, got $(sum(toprobs))"
                    fromcond = tree[node][1].fromcond  # all branches of a node MUST have the same fromcond
                    age_bump = filter(x->x!=fromcond,age_bump)   # remove fromcond distributed to new condition
                    folks = grab(fromcond,agegrp,lag,locale, dat=dat) # integer
                    if folks > 0
                        distribute_to_new_conditions!(folks, fromcond, toprobs, agegrp, lag, locale, node, dat=dat)
                    end
                end  
                if !isempty(age_bump)  # bump people in conds that didn't get distributed above
                    bump_up!(age_bump, agegrp, lag, locale, dat=dat) 
                end
            end
        else 
            # when no decision tree, bump up every infected person one day within the same condition
            bump_up!(infectious_cases, agegrps, lag, locale, dat=dat)
        end
    end

    update_infectious!(locale, dat = dat) # total all people who are nil, mild, sick, severe across all lags
    return
end


"""
function bump-up!

Bump people from one lag to lag + 1 in the same disease condition.
"""
function bump_up!(to_bump, agegrp, lag, locale; dat=openmx)
    bump = grab(to_bump, agegrp, lag, locale, dat=dat)
    if sum(bump) > 0
        plus!(bump, to_bump, agegrp, lag+1, locale, dat=dat)
        minus!(bump, to_bump, agegrp, lag,   locale, dat=dat)
    end
end


"""
function distribute_to_new_conditions!

Based on decision trees for each age group, at specific decision points (in days), change
people's disease condition, or move them to recovered or dead.
"""
function distribute_to_new_conditions!(folks, fromcond, toprobs, agegrp, lag, locale, node; dat=openmx, lastlag=laglim)

    @debug "day $(ctr[:day])  folks $folks lag $lag age $agegrp cond $fromcond"

    # set vector of folks to each outcome (6 outcomes): 1: recovered 2: nil 3: mild 4: sick 5: severe 6: dead
    @assert isapprox(sum(toprobs), 1.0, atol=1e-4) "target vector must sum to 1.0; submitted $toprobs"
    x = categorical_sample(toprobs, folks)  # integer results
    distvec = bucket(x, vals=1:length(toprobs))   # toprobs ALWAYS = 6     # [count(x .== i) for i in 1:size(toprobs,1)]
    @assert sum(distvec) == folks "someone got lost $res != $folks"

    if lag != lastlag  # infectious cases to next lag
        @views plus!(distvec[map2pr.nil:map2pr.severe], infectious_cases, agegrp, lag+1, locale, dat=dat) # @views 
    end
    @views plus!(distvec[map2pr.recovered], recovered, agegrp, 1, locale, dat=dat)  # recovered to lag 1
    @views plus!(distvec[map2pr.dead], dead, agegrp, 1, locale, dat=dat)  # dead to lag 1
    @views minus!(folks, fromcond, agegrp, lag, locale, dat=dat)  # subtract what we moved from the current lag

    @views push!(transq, (day=ctr[:day], lag=lag, agegrp=agegrp,   # @views primarily for debugging; can do some cool plots
                   newcond=distvec[map2pr.nil:map2pr.severe], recovered=distvec[map2pr.recovered],
                   dead=distvec[map2pr.dead], node=node, locale=locale))
    return
end


function update_infectious!(locale; dat=openmx) # by single locale
    for agegrp in agegrps
        tot = sum([nil, mild, sick, severe],agegrp,:,locale,dat=dat) # sum across cases and lags per locale and agegroup
        input!(tot, infectious, agegrp, 1, locale, dat=dat) # update the infectious total for the locale and agegroup
    end
end


"""
For a locale, randomly choose the number of people from each agegroup with
condition of {unexposed, infectious, recovered} who travel to each
other locale. Add to the travelq.
"""
function travelout!(fromloc, locales, rules=[])
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
function isolate!(pct::Float64,cond,agegrp,lag,locale; opendat=openmx, isodat=isolatedmx)
    for c in cond
        for age in agegrp
            for l in lag
                isolate_by!(pct::Float64,c,age,l,locale; opendat=opendat, isodat=isodat)
            end
        end
    end
end


function isolate_by!(pct::Float64,cond,agegrp,lag,locale; opendat=openmx, isodat=isolatedmx)
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, lag, locale, dat=opendat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, 0, available)  # limit to max
    cnt < scnt && (@warn "Attempt to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
end


function isolate_by!(num, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    @assert sum(num) >= 0 "num must be greater than zero"
    if typeof(locale) <: Quar_Loc
        available = grab(cond, agegrp, lag, locale.locale, dat=opendat)  # max
    else
        available = grab(cond, agegrp, lag, locale, dat=opendat)  # max
    end
    cnt = clamp.(num, 0, available)  # limit to max
    sum(cnt) < sum(num) && (@warn "Attempt to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing
end  # this one works


function _isolate!(cnt, cond, agegrp, lag, locale::Int; opendat=openmx, isodat=isolatedmx)
    minus!.(cnt, cond, agegrp, lag, locale, dat=opendat)  # move out 
    update_infectious!(locale, dat=opendat)
    plus!.(cnt, cond, agegrp, lag, locale, dat=isodat)  # move in
    update_infectious!(locale, dat=isodat)
    return nothing  # this one works!
end

# for test and trace or any isolate that records the day of isolation
function _isolate!(cnt, cond, agegrp, lag, qloc::Quar_Loc; opendat=openmx, isodat=isolatedmx)
    minus!(cnt, cond, agegrp, lag, qloc.locale, dat=opendat)  # move out 
    update_infectious!(qloc.locale, dat=opendat)
    plus!(cnt, cond, agegrp, lag, qloc, dat=isodat)  # move in
    update_infectious!(qloc, dat=isodat)
    return nothing  # this one works!
end

function unisolate!(pct::Float64,cond,agegrp,lag,locale; opendat=openmx, isodat=isolatedmx)
    for c in cond
        for age in agegrp
            for l in lag
                unisolate_by!(pct::Float64,c,age,l,locale; opendat=opendat, isodat=isodat)
            end
        end
    end
end


function unisolate_by!(pct::Float64,cond,agegrp,lag,locale; opendat = openmx, isodat=isolatedmx)
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, lag, locale, dat=isodat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, 0, available)  # limit to max
    cnt < scnt && (@warn "Attempt to unisolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing  # this one works!
end


function unisolate_by!(num, cond, agegrp, lag, locale; mode=:both, opendat=openmx, isodat=isolatedmx)
    @assert sum(num) >= 0 "sum(num) must be greater than zero"

    available = grab(cond, agegrp, lag, locale, dat=isodat)  # max

    println("day $(ctr[:day]) request to unisolate   ", sum(num))
    println("day $(ctr[:day]) available to unisolate ", sum(available))

    cnt = clamp.(num, 0, available)  # limit to max
    sum(cnt) < sum(num) && (@warn "Attempt to unisolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale; mode=:both, opendat=opendat, isodat=isodat)
    return nothing
end  # this one works


function _unisolate!(cnt, cond, agegrp, lag, locale; mode=:both, opendat=openmx, isodat=isolatedmx)
    if mode != :plus  # this is when mode = :minus or :both
        minus!(cnt, cond, agegrp, lag, locale, dat=isodat)
        update_infectious!(locale, dat=isodat)
    end
    if mode != :minus  # this is when mode = :plus or :both
        if typeof(locale) <: Quar_Loc
            locale = locale.locale
        end

        println("day $(ctr[:day])  unquarantine is unisolating this many ", sum(cnt))

        plus!(cnt, cond, agegrp, lag, locale, dat=opendat)
        update_infectious!(locale, dat=opendat)
    end
    return 
end