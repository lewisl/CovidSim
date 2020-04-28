#########################################################################################
# tracking.jl
#########################################################################################


const map2series = (unexposed=1:6, infectious=7:12, recovered=13:18, dead=19:24, 
                    nil=25:30, mild=31:36, sick=37:42, severe=43:48)
const total = 6

"""
- use incr!(ctr, :day) for day of the simulation:  creates and adds 1
- use reset!(ctr, :day) to remove :day and return its current value, set it to 0
- use ctr[:day] to return current value of day
"""
const ctr = counter(Symbol) # from package DataStructures

# for debugging whole simulations
const spreadq = []



# for Johns Hopkins US actual data
struct Col_ref
    date::String
    col::Int64
end



# tracking statistics

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


function cumplot(series, locale, plcols=[unexposed, infectious, recovered, dead]; days="all",geo=[])

    pyplot()
    theme(:ggplot2, foreground_color_border =:black, reuse = false)

    !(typeof(plcols) <: Array) && (plcols = [plcols])

    # the data
    n = size(series[locale][:cum],1)
    days = days == "all" ? (1:n) : days
    cumseries = series[locale][:cum][days, [map2series[i][total] for i in plcols]]

    # labels and annotations
    labels = [titlecase(condnames[i]) for i in plcols]
    labels = reshape([labels...], 1, length(labels))
    people = series[locale][:cum][1, map2series[unexposed][total]] + series[locale][:cum][1,map2series[infectious][total]]
    cityname = !isempty(geo) ? geo[locale, city] : ""
    died = series[locale][:cum][end, map2series[dead][total]]
    infected = series[locale][:cum][1,map2series[unexposed][total]] - series[locale][:cum][end,map2series[unexposed][total]]
    firstseries = plcols[1]
    half_yscale = floor(Int, maximum(series[locale][:cum][:,map2series[firstseries][total]]) * 0.7)
 

    # the plot
    plot(   days, cumseries[days,1:end], 
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
    theme(:ggplot2, foreground_color_border =:black)
    
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