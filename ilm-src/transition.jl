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
                make_sick!(dat[loc]; cnt=cnt, fromage=agegrps, tocond=nil, tolag=lag)
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
function transition!(dt, all_decpoints, locale, dat, agegrp_idx)  

    # @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single integer or NamedTuple"
    locdat = dat[locale]
    locagegrp_idx = agegrp_idx[locale]
    filt_agegrp = [falses(size(locagegrp_idx[age],1)) for age in agegrps]


    for agegrp in agegrps
        tree = dt[agegrp]
        for node in sort(collect(keys(tree)), rev=true)
            nodelag, fromcond = node
            filt = ( (locdat[locagegrp_idx[agegrp],cpop_cond] .== fromcond) .& 
                     (locdat[locagegrp_idx[agegrp],cpop_status] .== infectious) .&
                     (locdat[locagegrp_idx[agegrp], cpop_lag] .== nodelag)  )

            # indices for the entire population by agegrp and the anded logical filters
            filt = locagegrp_idx[agegrp][filt]
            cnt_transition = size(filt,1)      

            if cnt_transition > 0
                probs = tree[node]["probs"]
                outcomes = tree[node]["outcomes"]
                choices = rand(Categorical(probs), cnt_transition) # returns ordinal Int index to outcomes     
                distrib = countmap(choices)
                for (i, person) in enumerate(filt)       # findall(filt)                 
                    tocond = outcomes[choices[i]]
                    if tocond in (dead, recovered)  # change status
                        locdat[person, cpop_status] = tocond
                        # leave disease condition and lag alone so we can see last state before death or recovery
                    else   # change disease condition
                        locdat[person, cpop_cond] = tocond
                    end  # if/else
                end  # for (idx, person)
                push!(transq, (day=ctr[:day], lag=nodelag, agegrp=agegrp, node=node, locale=locale,   # @views primarily for debugging; can do some cool plots
                                recovered=get(distrib, indexin(recovered,outcomes)[], 0),
                                dead=   get(distrib, indexin(dead,outcomes)[], 0),
                                nil=    get(distrib, indexin(nil,outcomes)[], 0),
                                mild=   get(distrib, indexin(mild,outcomes)[], 0),
                                sick=   get(distrib, indexin(sick,outcomes)[], 0),
                                severe= get(distrib, indexin(severe,outcomes)[], 0))
                       )

            end  # if cnt_transition
        end  # for node
    end  #for agegrp

    # bump everyone who is still infectious all at once in one go
    filt = findall(locdat[:,cpop_status] .== infectious)
    locdat[filt, cpop_lag] .+= 1   
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