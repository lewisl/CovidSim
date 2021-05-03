####################################################
# transition.jl for ilm model
#     change status of folks in simulation:
#           transition
#           travel
#           seed
#           isolate
####################################################


function seed!(day, cnt, sickday, conds, agegrps, locale, dat)
    @assert length(sickday) == 1 "input only one sickday value"
    # @warn "Seeding is for testing and may result in case counts out of balance"
    if day == day_ctr[:day]
        println("*** seed day $(day_ctr[:day]) locale $locale....")
        for loc in locale
            for cond in conds
                @assert (cond in [nil, mild, sick, severe]) "Seed cases must have conditions of nil, mild, sick, or severe" 
                make_sick!(dat[loc]; cnt=cnt, fromage=agegrps, tocond=nil, tosickday=sickday)
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
@inline function transition!(locdat, infect_idx, dectree)

    # locdat = locale == 0 ? dat : dat[locale] #allocates 112 bytes

    # dectree = dt_dict["dt"]  # decision tree for illness changes over time to recovery or death

    # infect_idx = optfindall(==(infectious), locdat.status, 0.5)  # allocates thousands of bytes

    for p in infect_idx  # p for person    
        p_sickday = locdat.sickday[p] 
        p_cond = locdat.cond[p]
        p_agegrp = locdat.agegrp[p]  # agegroup of person p = agegrp column of locale data, row p 
        if !haskey(dectree[p_agegrp], p_sickday) || !haskey(dectree[p_agegrp][p_sickday], p_cond)
            # @assert p_tup.sickday < sickdaylim "Person made it to last day and was not removed:\n     $p_tup\n"
            locdat.sickday[p] += 1
        else  # change of the person p's state--a transition
            node = dectree[p_agegrp][p_sickday][p_cond]
            choice = categorical_sim(node["probs"]) # rand(Categorical(node["probs"])) # which branch...?
            tocond = node["outcomes"][choice]
            if tocond == dead  # change status, leave cond and sickday as last state before death or recovery                        
                locdat.status[p] = dead  # change the status
                locdat.dead_day[p] = day_ctr[:day]
                locdat.cond[p] = notsick
            elseif tocond == recovered
                locdat.recov_day[p] = day_ctr[:day]
                locdat.status[p] = recovered
                locdat.cond[p] = notsick
            else   # change disease condition
                locdat.cond[p] = tocond   # change the condition
                locdat.sickday[p] += 1  
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
function travelin!(dat=popdat)
    while !isempty(travelq)
        g = dequeue!(travelq)
        cond = eval(Symbol(g.cond))
        minus!(g.cnt, cond, g.agegrp, g.sickday, g.from, dat=dat)
        plus!(g.cnt, cond, g.agegrp, g.sickday, g.to, dat=dat)
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

    thisday = day_ctr[:day]
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


function in_quarantine(locdat, p, leakage = 0.05)::Bool
    @assert 0.0 <= leakage <= 1.0 "leakage value must be between 0.0 and 1.0 inclusive"
    if leakage == 1.0
        return false  # everyone leaks quarantine is always false
    elseif leakage == 0.0
        return locdat.quar[p] # no one leaks--depends on the person's status
    else
        rand() < leakage ? false : locdat.quar[p]
    end
end