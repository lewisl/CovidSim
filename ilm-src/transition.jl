####################################################
# transition.jl for ilm model
#     change status of folks in simulation:
#           transition
#           travel
####################################################

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

    for p in infect_idx  # p for person    
        p_sickday = Int(locdat.sickday[p]) 
        p_cond = Int(locdat.cond[p])
        p_agegrp = Int(locdat.agegrp[p])  # agegroup of person p = agegrp column of locale data, row p 
        if haskey(dectree[p_agegrp], p_sickday) && haskey(dectree[p_agegrp][p_sickday], p_cond)
            # change the person p's state--a transition
            node = dectree[p_agegrp][p_sickday][p_cond]
            choice = categorical_sim(node["probs"]) # rand(Categorical(node["probs"])) # which branch...?
            tocond = node["outcomes"][choice]
            if tocond == Int(dead)  # change status, leave cond and sickday as last state before death or recovery                        
                locdat.status[p] = dead  # change the status
                locdat.dead_day[p] = day_ctr[:day]
                locdat.cond[p] = notsick
            elseif tocond == Int(recovered)
                locdat.recov_day[p] = day_ctr[:day]
                locdat.status[p] = recovered
                locdat.cond[p] = notsick
            else   # change disease condition
                locdat.cond[p] = condition(tocond)   # change the condition = degree of sickness
                locdat.sickday[p] += 1  
            end    
        else # just increment sickday: one more day feeling the same kind of sick
            # @assert p_tup.sickday < sickdaylim "Person made it to last day and was not removed:\n     $p_tup\n"
            locdat.sickday[p] += 1
        end  # if haskey ... else
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