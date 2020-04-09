#########################################################################################
# tracking.jl
#########################################################################################

const series_colnames = Dict( 1=>:Unexposed,  2=>:Infectious, 3=>:Recovered, 4=>:Dead, 5=>:Nil, 6=>:Mild, 7=>:Sick,
        8=>:Severe,  9=>:Travelers, 10=>:Isolated)


# traveling 
    # inbound travelers: exogenous, domestic, outbound (assumed domestic) travelers
    const travelq = Queue{NamedTuple{(:cnt, :from, :to, :agegrp, :lag, :cond),
        Tuple{Int64,Int64,Int64,Int64,Int64,String}}}()

    function travitem(cnt, from, to, agegrp, lag, cond)
        return (cnt=cnt, from=from, to=to, agegrp=agegrp, lag=lag, cond=cond)
    end


# isolation   TODO  maybe not using this
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
    function isolatestat(; cnt=0, locale=0, tocond=isolated)
        (cnt=cnt, locale=locale, tocond=tocond)
    end

    SpreadStat = typeof((;cnt=0, locale=0, tocond=nil))
    function spreadstat(;cnt=0,locale=0,tocond=nil)
       (cnt=cnt, locale=locale, tocond=tocond)
    end


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


# for isolation
function queuestats(cnt, locale, func::typeof(isolatestat); case="open")
    thisday = ctr[:day]
    if cnt == 0
    else
        # net addition to isolated; no change to any condition
        additem = func(cnt = cnt, locale=locale)  # defaults to isolated
        push!(newstatq, (day=thisday, additem...))
    end
end



"""
- use incr!(ctr, :day) for day of the simulation:  creates and adds 1
- use reset!(ctr, :day) to remove :day and return its current value, set it to 0
- use ctr[:day] to return current value of day
"""
const ctr = counter(Symbol) # from package DataStructures

# all locales
function queue_to_newseries!(newstatq, dseries, locales)

    thisday = ctr[:day]
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


function showq(qname)
    for item in qname
        println(item)
    end
end
