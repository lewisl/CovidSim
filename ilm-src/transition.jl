####################################################
# transition.jl for ilm model
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
function transition!(dat, dt_dict)
    for locale in keys(dat)
        transition!(dat, locale, dt_dict)
    end
end


    
"""
    transition!(dat, locale, dt_dict)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability. Finally,  
they move to recovered or dead.

Works for a single locale.
"""
function transition!(dat, locale::Int, dt_dict)

    locdat = locale == 0 ? dat : dat[locale]

    dectree = dt_dict["dt"]  # decision tree for illness changes over time to recovery or death

    @inbounds for p in findall(locdat.status .== infectious)  # p for person          

                dtkey = (locdat.lag[p], locdat.cond[p])
                p_agegrp = locdat.agegrp[p]  # agegroup of person p = agegrp column of locale data, row p 
                node = get(dectree[p_agegrp], dtkey, ())
                @inbounds if isempty(node)  # no transition for this lag/day of the disease and current condition
                    # @assert p_tup.lag < laglim "Person made it to last day and was not removed:\n     $p_tup\n"
                    locdat.lag[p] += 1
                else  # change of the person p's state--a transition
                    choice = rand(Categorical(node["probs"]), 1) # which branch...?
                    tocond = node["outcomes"][choice][]
                    if tocond == dead  # change status, leave cond and lag as last state before death or recovery                        
                        locdat.status[p] = dead  # change the status
                        locdat.dead_day[p] = ctr[:day]
                        locdat.cond[p] = notsick
                    elseif tocond == recovered
                        locdat.recov_day[p] = ctr[:day]
                        locdat.status[p] = recovered
                        locdat.cond[p] = notsick
                    else   # change disease condition
                        locdat.cond[p] = tocond   # change the condition
                        locdat.lag[p] += 1  
                    end    
                end
    end  # for p
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
function travelin!(dat=popdat)
    while !isempty(travelq)
        g = dequeue!(travelq)
        cond = eval(Symbol(g.cond))
        minus!(g.cnt, cond, g.agegrp, g.lag, g.from, dat=dat)
        plus!(g.cnt, cond, g.agegrp, g.lag, g.to, dat=dat)
    end
end


"""
    Struct for supplying terms that make up a query of the population matrix.
    Fields of the struct are:
    - col: Symbol that is the column name in the population matrix
    - op: logical test operation, one of: ==, !=, <, >, <=, >=, etc.
    - val: the value to be compared, either Int or bit type for Bool: true, false
"""
struct Popquery
    col::Symbol
    op::Function
    val::Union{Int64, Bool}
    # row::Int
end


function doquery(dat, pq::Popquery, row)
    return (pq.op)(pq.val, getproperty(dat, pq.col)[row])
end


# TODO THIS NO LONGER WORKS!    

"""
    function isolate!(locdat, qty, filters::Array, val::Bool=true)

Place people into isolation or quarantine when val is the default: true.
Remove people from quarantine by setting val to false.

"""
function isolate!(locdat, qty, filters::Array, val=true)  
    # val = true adds to quarantine; val = false removes from quarantine

    thisday = ctr[:day]
    filts = copy(filters)

    push!(filts, Popquery(:quar, ==, !val)) # don't isolate people already in quarantine

    # break apart filters into tuples of values, operations, columns
    qvals = Tuple(q.val for q in filts)
    qops = Tuple(q.op for q in filts)
    qcols = Tuple(q.col for q in filts)

    isoq = 0
    rand_idx = randperm(size(locdat.status, 1))
    for samp in Iterators.partition(rand_idx, qty), p in samp
        conditions = true
        for (j,col) in enumerate(qcols)
            # conditions &= (qops[j])(qvals[j],  getproperty(locdat, col)[p]) 
            conditions &= doquery(locdat, Popquery(col, qops[j], qvals[j]), p)
        end
        
        if conditions
            isoq += 1
            locdat.quar[p] = val
            locdat.quar_day[p] = val ? thisday : 0
        end

        if isoq >= qty
            break
        end
    end

    sim_stash[:quar_cnt] = (val ? get(sim_stash, :quar_cnt, 0) + isoq 
                                :  get(sim_stash, :quar_cnt, 0) - isoq) 

    return isoq
end


"""
Return boolean filter of people currently in quarantine less
daily leakage, if applicable.
"""
function current_quar(locdat, leakage = .05)
    @assert 0.0 <= leakage <= 1.0 "leakage value must be between 0.0 and 1.0 inclusive"

    iq_filt = copy(locdat.quar)
    cnt = round(Int, sum(locdat.quar) * leakage)
    
    iq_filt[sample(findall(locdat.quar), cnt, replace=false)] .= false
    
    return iq_filt
end
