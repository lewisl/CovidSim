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

function seed!(day, cnt, sickday, conds, agegrps, locale, dat)
    @assert length(sickday) == 1 "input only one sickday value"
    # @warn "Seeding is for testing and may result in case counts out of balance"
    if day == day_ctr[:day]
        println("*** seed day $(day_ctr[:day]) locale $locale....")
        for loc in locale
            for cond in conds
                @assert (cond in [nil, mild, sick, severe]) "Seed cases must have conditions of nil, mild, sick, or severe" 
                plus!(cnt, cond, agegrps, sickday, loc, dat)
                minus!(cnt, unexposed, agegrps, 1, loc, dat)
                update_infectious!(loc, dat)
            end
        end
    end
end


# method to run through all existing locales in isolation
function transition!(dat, dt_dict)
    for locale in keys(dat)
        transition!(dat, dt_dict, locale)
    end
end



"""
    transition!(dat, dt_dict, locale)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability; then to 
recovered or dead.

Works for a single locale.
"""
function transition!(dat, dt_dict, locale)  

    # @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single integer or NamedTuple"
    # iszero(dat[locale]) && (return)

    for agegrp in agegrps
        agetree = dt_dict["dt"][agegrp]  # tree for a single agegrp
        for sickday in sort(collect(keys(agetree)), rev=true)  # tree for a sickday value
            sickdaytree = agetree[sickday]
            for fromcond in sort(collect(keys(sickdaytree)))
                condtree = sickdaytree[fromcond]
                folks = grab(fromcond, agegrp, sickday, locale, dat)

                if folks > 0
                    pr = condtree["probs"] # pr for all branches at the node
                    outcomes = condtree["outcomes"] # outcomes for all branches at the node

                    distrib = countmap(rand(Categorical(pr), folks))

                    for i in keys(distrib)
                        if outcomes[i] in [recovered, dead]
                            plus!(distrib[i], outcomes[i], agegrp, 1, locale, dat) # sickday is 1
                            minus!(distrib[i], fromcond, agegrp, sickday, locale, dat)
                        else  # in infectious conditions nil:severe
                            plus!(distrib[i], outcomes[i], agegrp, sickday, locale, dat)
                            minus!(distrib[i], fromcond, agegrp, sickday, locale, dat)
                        end
                    end
                    push!(transq, 
                            (day=day_ctr[:day], sickday=sickday, agegrp=agegrp, fromcond=fromcond, locale=locale,   # @views primarily for debugging; can do some cool plots
                             recovered = get(distrib, indexin(recovered,outcomes)[], 0),
                             dead      = get(distrib, indexin(dead,outcomes)[], 0),
                             nil       = get(distrib, indexin(nil,outcomes)[], 0),
                             mild      = get(distrib, indexin(mild,outcomes)[], 0),
                             sick      = get(distrib, indexin(sick,outcomes)[], 0),
                             severe    = get(distrib, indexin(severe,outcomes)[], 0)))
                end  # if folks
            end  # for fromcond
        end  # for sickday
    end  # for agegrp

    for sickday in sickdaylim-1:-1:1
        bump_up!(infectious_cases, agegrps, sickday, locale, dat)
    end

    update_infectious!(locale, dat)
end


"""
function bump-up!

Bump people from one sickday to sickday + 1 in the same disease condition.
"""
function bump_up!(to_cond, agegrp, sickday, locale, dat)
    bump = grab(to_cond, agegrp, sickday, locale, dat)

    if sum(bump) > Int(0)
        plus!(bump, to_cond, agegrp, sickday+1, locale, dat)
        minus!(bump, to_cond, agegrp, sickday,   locale, dat)
    end
end


function update_infectious!(locale, dat) # by single locale
    for agegrp in agegrps
        tot = sum(grab([nil, mild, sick, severe],agegrp,:,locale,dat)) # sum across cases and sickdays per locale and agegroup
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
        # unexposed, infectious, recovered -> ignore sickday for now
    # TODO: more frequent travel to and from Major and Large cities
    # TODO: should the caller do the loop across locales?   YES
    travdests = collect(locales)
    deleteat!(travdests,findfirst(isequal(fromloc), travdests))
    bins = lim = length(travdests) + 1
    for agegrp in agegrps
        for cond in [unexposed, infectious, recovered]
            name = condnames[cond]
            for sickday in sickdays
                numfolks = sum(grab(cond, agegrp, sickday, fromloc)) # the from locale, all sickdays
                travcnt = floor(Int, gamma_prob(travprobs[agegrp]) * numfolks)  # interpret as fraction of people who will travel
                x = rand(travdests, travcnt)  # randomize across destinations
                bydest = bucket(x, vals=1:length(travdests))
                for dest in 1:length(bydest)
                    isempty(bydest) && continue
                    cnt = bydest[dest]
                    iszero(cnt) && continue
                    enqueue!(travelq, travitem(cnt, fromloc, dest, agegrp, sickday, name))
                end
            end
        end
    end
end


"""
Assuming a daily cycle, at the beginning of the day
process the queue of travelers from the end of the previous day.
Remove groups of travelers by agegrp, sickday, and condition
from where they departed.  Add them to their destination.
"""
function travelin!(dat=openmx)
    while !isempty(travelq)
        g = dequeue!(travelq)
        cond = eval(Symbol(g.cond))
        minus!(g.cnt, cond, g.agegrp, g.sickday, g.from, dat=dat)
        plus!(g.cnt, cond, g.agegrp, g.sickday, g.to, dat=dat)
    end
end


"""
Place people into isolation.

You can enter a percentage (as a fraction in [0.0, 1.0]), an array
of percentages, a number, or an array of numbers.

Use a dot after the function name to apply the same pct or number
input to several conditions, agegrps, sickdays, or locales.

Use a dot after the function name to apply an array: one or more
of agegrp, cond, or locale must have the same number of elements as the input.
"""
function isolate!(pct::Float64,cond,agegrp,sickday,locale, opendat, isodat)
    for c in cond
        for age in agegrp
            for l in sickday
                isolate_by!(pct::Float64,c,age,l,locale, opendat, isodat)
            end
        end
    end
end


function isolate_by!(pct::Float64,cond,agegrp,sickday,locale, opendat, isodat)
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, sickday, locale, opendat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, Int(0), Int(available))  # limit to max
    cnt < scnt && (@warn "Attempt to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, sickday, locale, opendat, isodat)
end


function isolate_by!(num, cond, agegrp, sickday, locale, opendat, isodat)
    @assert sum(num) >= 0 "num must be greater than zero"
    if typeof(locale) <: Quar_Loc
        available = grab(cond, agegrp, sickday, locale.locale, opendat)  # max
    else
        available = grab(cond, agegrp, sickday, locale, opendat)  # max
    end
    cnt = clamp.(num, Int(0), Int(available))  # limit to max
    sum(cnt) < sum(num) && (@warn "Attempt to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, sickday, locale, opendat, isodat)
    return nothing
end  # this one works


function _isolate!(cnt, cond, agegrp, sickday, locale::Integer, opendat, isodat)
    minus!(cnt, cond, agegrp, sickday, locale, opendat)  # move out 
    update_infectious!(locale, opendat)
    plus!(cnt, cond, agegrp, sickday, locale, isodat)  # move in
    update_infectious!(locale, isodat)
    return nothing  # this one works!
end

# for test and trace or any isolate that records the day of isolation
function _isolate!(cnt, cond, agegrp, sickday, qloc::Quar_Loc, opendat, isodat)
    minus!(cnt, cond, agegrp, sickday, qloc.locale, opendat)  # move out 
    update_infectious!(qloc.locale, opendat)
    plus!(cnt, cond, agegrp, sickday, qloc, isodat)  # move in
    update_infectious!(qloc, isodat)
    return nothing  # this one works!
end

function unisolate!(pct::Float64,cond,agegrp,sickday,locale, opendat, isodat)
    for c in cond
        for age in agegrp
            for l in sickday
                unisolate_by!(pct::Float64,c,age,l,locale, opendat, isodat)
            end
        end
    end
end


function unisolate_by!(pct::Float64,cond,agegrp,sickday,locale, opendat, isodat)
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, sickday, locale, isodat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, Int(0), Int(available))  # limit to max
    cnt < scnt && (@warn "Attempt to unisolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, sickday, locale, opendat, isodat)
    return nothing  # this one works!
end


function unisolate_by!(num, cond, agegrp, sickday, locale, opendat, isodat, mode=:both)
    @assert sum(num) >= 0 "sum(num) must be greater than zero"

    available = grab(cond, agegrp, sickday, locale, isodat)  # max

    # println("day $(day_ctr[:day]) request to unisolate   ", sum(num))
    # println("day $(day_ctr[:day]) available to unisolate ", sum(available))

    cnt = clamp.(num, Int(0), Int(available))  # limit to max
    sum(cnt) < sum(num) && (@warn "Attempt to unisolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, sickday, locale,  opendat, isodat, :both)
    return nothing
end  # this one works


function _unisolate!(cnt, cond, agegrp, sickday, locale,  opendat, isodat, mode=:both)
    if mode != :plus  # this is when mode = :minus or :both
        minus!(cnt, cond, agegrp, sickday, locale, isodat)
        update_infectious!(locale, isodat)
    end
    if mode != :minus  # this is when mode = :plus or :both
        if typeof(locale) <: Quar_Loc
            locale = locale.locale
        end

        # println("day $(day_ctr[:day])  unquarantine is unisolating this many ", sum(cnt))

        plus!(cnt, cond, agegrp, sickday, locale, opendat)
        update_infectious!(locale, opendat)
    end
    return 
end