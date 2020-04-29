####################################################
# transition.jl
#     change status of folks in simulation:
#           transition
#           travel
#           seed
#           isolate
####################################################


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


"""
    transition!(dt, locale; case="open", dat=openmx)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability. The final
outcomes are recovered or dead.
"""
function transition!(dt, locale; dat=openmx)  # TODO also need to run for isolatedmx

    for lag = laglim:-1:1
        dec_tree_applied = false
        for agegrp in agegrps
            tree = dt[agegrp].tree
            node_decpoints = dt[agegrp].dec_points
            if lag in keys(node_decpoints) # check if a decision tree applies today
                for node in node_decpoints[lag]  
                    toprobs = zeros(6)
                    for branch in tree[node]  # agegroup index in array, node key in agegroup dict
                        toprobs[mapcond2tran[branch.tocond]] = branch.pr
                    end
                    fromcond = tree[node][1].fromcond  # all branches of a node MUST have the same fromcond
                    @debug @sprintf("%12s %3f %3f %3f %3f %3f %3f",condnames[fromcond], toprobs...)
                    distribute_to_new_conditions!(fromcond, toprobs, agegrp, lag, locale, dat=dat)
                end
                dec_tree_applied = true
            end
        end

        if lag == laglim  # there must be a node that starts on laglim and clears everyone left on the last lag!
            if sum(grab(infectious_cases, agegrps, laglim, locale, dat=dat)) !== 0
                @warn  "infectious conditions not zero at lag $laglim, day $(ctr[:day])--found $(sum(grab(infectious_cases, agegrps, laglim, locale, dat=dat)))"
            end
            if !dec_tree_applied # e.g., dec_tree was NOT applied for final lag
                @warn "decision node not applied at last lag--probably wrong outcomes"
            end
        end

        if dec_tree_applied
            continue  # EITHER move people with a decision node OR by just bumping them up
        end

        # on a day with no decision tree, bump every infected person up one day within the same condition
        input!(grab(nil:severe, agegrps, lag, locale, dat=dat), nil:severe, agegrps, lag+1, locale, dat=dat)
        minus!(grab(nil:severe, agegrps, lag, locale, dat=dat), nil:severe, agegrps, lag,   locale, dat=dat)
    end

    # total of all people who are nil, mild, sick, or severe across all lag days
    update_infectious!(locale, dat = dat)
end


function distribute_to_new_conditions!(fromcond, toprobs, agegrp, lag, locale; dat=openmx, lastlag=laglim)
    @assert length(locale) == 1  "Assertion failed: length locale was not 1"
    transition_cases = [recovered, nil,mild,sick,severe, dead]

    # get the number of folks to be distributed
    folks = grab(fromcond,agegrp,lag,locale, dat=dat) # scalar
    if folks == 0
        return   # save some time if there is nothing to do
    end
    @debug "day $(ctr[:day])  folks $folks lag $lag age $agegrp cond $fromcond"
    # println("on day $(ctr[:day]) folks $folks lag $lag from $fromcond to $toprobs")


    map2pr = (unexposed = -1, infectious = -1, recovered = 1, dead = 6, nil = 2, mild = 3, sick = 4, severe = 5)

    # get the distvect of folks to each outcome (6 outcomes): 1: recovered 2: nil 3: mild 4: sick 5: severe 6: dead
    @assert isapprox(sum(toprobs), 1.0, atol=1e-3) "target vector must sum to 1.0; submitted $toprobs"
    x = categorical_sample(toprobs, folks)
    # println(x)
    distvec = bucket(x, vals=1:length(toprobs))        # [count(x .== i) for i in 1:size(toprobs,1)]
    @debug begin
        cond = condnames[fromcond];rec=distvec[1]; ni=distvec[2]; mi=distvec[3]; si=distvec[4]; se=distvec[5]; de=distvec[6];
        "distribute $cond age $agegrp lag $lag CNT $folks to $rec $nil $mi $si $se $dead"
        end

    lag != lastlag && (plus!(distvec[map2pr.nil:map2pr.severe], infectious_cases, agegrp, lag+1, locale, dat=dat)) # add to infectious cases for next lag
    plus!(distvec[map2pr.recovered], recovered, agegrp, 1, locale, dat=dat)  # recovered to lag 1
    plus!(distvec[map2pr.dead], dead, agegrp, 1, locale, dat=dat)  # dead to lag 1
    minus!(sum(distvec), fromcond, agegrp, lag, locale, dat=dat)  # subtract what we moved from the current lag

    return
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
                bydest = bucket(x, vals=1:length(travdests))
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
Place people who are isolating into the isolated queue: isolatedq

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
    cnt < scnt && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
end


function isolate_by!(num::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    @assert num > 0 "num must be greater than zero"
    available = grab(cond, agegrp, lag, locale, dat=opendat)  # max
    cnt = clamp(num, 0, available)  # limit to max
    cnt < num && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing
end  # this one works


function _isolate!(cnt::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    minus!(cnt, cond, agegrp, lag, locale, dat=opendat)  # move out 
    update_infectious!(locale, dat=opendat)
    plus!(cnt, cond, agegrp, lag, locale, dat=isodat)  # move in
    update_infectious!(locale, dat=isodat)
    cnt !== 0 && queuestats(cnt=cnt, locale=locale, conds=cond, agegrp=agegrp, event=:isolate)
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
    cnt < scnt && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing  # this one works!
end


function unisolate_by!(num::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    @assert num > 0 "num must be greater than zero"
    available = grab(cond, agegrp, lag, locale, dat=isodat)  # max
    cnt = clamp(num, 0, available)  # limit to max
    cnt < num && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing
end  # this one works


function _unisolate!(cnt::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    minus!(cnt, cond, agegrp, lag, locale, dat=isodat)
    update_infectious!(locale, dat=isodat)
    plus!(cnt, cond, agegrp, lag, locale, dat=opendat)
    update_infectious!(locale, dat=opendat)
    cnt !== 0 && queuestats(cnt=-cnt, locale=locale, conds=cond, agegrp=agegrp, event=:isolate)
    return nothing
end