#########################################################################################
# tracking.jl
#########################################################################################

const series_colnames = Dict( 1=>:Unexposed,  2=>:Infectious, 3=>:Recovered, 4=>:Dead, 5=>:Nil, 6=>:Mild, 7=>:Sick,
        8=>:Severe,  9=>:Travelers, 10=>:Isolated)


"""
- use incr!(ctr, :day) for day of the simulation:  creates and adds 1
- use reset!(ctr, :day) to remove :day and return its current value, set it to 0
- use ctr[:day] to return current value of day
"""
const ctr = counter(Symbol) # from package DataStructures

# for debugging whole simulations
const spreadq = []


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



# for Johns Hopkins US actual data
struct Col_ref
    date::String
    col::Int64
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

    # add :Total_infected to match Johns Hopkins statistics
    cumseries[!, :Total_infected] = cumseries[!, :Infectious] + cumseries[!, :Recovered] + cumseries[!, :Dead]
end


function showq(qname)
    for item in qname
        println(item)
    end
end


function reviewdays(q=spreadq)
    for it in q
        println(it)
        print("\nPress enter to continue, q enter to quit.> ");
        ans = chomp(readline())
        if ans == "q"
            break
        end
    end
end

function reviewdays(df::DataFrame)
    for it in eachrow(df)
        display(it)
        print("\nPress enter to continue, q enter to quit.> ");
        ans = chomp(readline())
        if ans == "q"
            break
        end
    end
end


###########################################################################################
#  Plotting
###########################################################################################


function cumplot(dseries, locale, plseries=[:Unexposed,:Infectious,:Recovered, :Dead];
    days="all",geo=[])

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
    half_yscale = floor(Int, maximum(dseries[locale][:cum][!,firstseries]) * 0.7)
    days = days == "all" ? (1:n) : days

    # the plot
    plot(   cumseries[days,1], cumseries[days,2:end], 
            size = (700,500),
            label = labels, 
            lw=2.3,
            title = "Covid for $people people in $cityname over $n days\nActive Cases for Each Day",
            xlabel = "Simulation Days",
            ylabel = "People",
            legendfontsize = 10,
            reuse = false
        )
    annotate!((6,half_yscale,Plots.text("Died: $died\nInfected: $infected", 10, :left)))
    gui()
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
    groupedbar(    newseries[:,1], newseries[:,2:end], 
            size = (700,500),
            label = labels, 
            lw=0,
            title = "Covid Daily Change for $people people over $n days",
            xlabel = "Simulation Days",
            yaxis = ("People"),
            reuse =false
        )
    gui()
end


function day2df(spreadq::Array)
    spreadseries = DataFrame(spreadq)

    spreadseries[!, :cuminfected] .= zeros(Int, size(spreadseries,1))
    spreadseries[1, :cuminfected] = copy(spreadseries[1,:infected])
    for i = 2:size(spreadseries,1)
       spreadseries[i,:cuminfected] = spreadseries[i-1,:cuminfected] + spreadseries[i,:infected]
    end

    return spreadseries
end


function dayplot(spreadq, plseries=[])
    dayplot(DataFrame(spreadq), plseries)
end

function dayplot(spreadseries::DataFrame, plseries=[])
    pyplot()
    plot(spreadseries[!,:day], spreadseries[!,:infected],label="Infected", dpi=200,lw=2,
         xlabel="Simulation Days", ylabel="People", title="Daily Spread of Covid",
         bg_legend=:white)
    
    for addlseries in plseries
        lbl = titlecase(string(addlseries))
        plot!(spreadseries[!,:day], spreadseries[!,addlseries],label=lbl, dpi=200,lw=2)
        # plot!(spreadseries[!,:day], spreadseries[!,:touched],label="Touched", dpi=200,lw=2)
        # plot!(spreadseries[!,:day], spreadseries[!,:spreaders],label="Spreaders", dpi=200,lw=2)
    end
    gui()  # force instant plot window
    # plot!(spreadseries[!,:day], spreadseries[!,:cuminfected],label="Cum Infected", dpi=200,lw=2)
end


function day_animate2(spreadseries)
    n = size(spreadseries,1)
    # daymat = Matrix(spreadseries)

    xd = spreadseries[1:5,:]

    topy = max(maximum(spreadseries[!,:spreaders]),maximum(spreadseries[!,:contacts]),
                maximum(spreadseries[!,:touched]),maximum(spreadseries[!,:infected]) )

    @df xd plot(:day, [:spreaders :contacts :touched :infected], color=^([:red :blue :green :orange]),
                labels=^(["Spreaders" "Contacts" "Touched" "Infected"]),dpi=200, lw=2,ylim=(0,topy))

    for i = 5:2:n
        xd = spreadseries[i-2:i,:]

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