
####################################################################################
#   simulation runner: ILM Model
####################################################################################


function run_a_sim(n_days, locales; runcases=[], spreadcases=[], showr0 = true, silent=true, set_int_type=Int64,
            geofilename="../data/geo2data.csv", 
            dectreefilename="../parameters/dec_tree_all_25.yml",
            spfilename="../parameters/spread_params.yml")

    empty_all_caches!() # from previous runs

    T_int[] = set_int_type # update the global type of ints with the input value

    # access input data and pre-allocate storage
    alldict = setup(n_days, locales; geofilename=geofilename, 
                    dectreefilename=dectreefilename, spfilename=spfilename)

        dectree = alldict["dt_dict"]["dt"]  # decision trees for transition
        popdat = alldict["dat"]["popdat"]   # first key locale
        agegrp_idx = alldict["dat"]["agegrp_idx"]   # first key locale
        cumhistmx = alldict["dat"]["cumhistmx"]   # first key locale
        newhistmx = alldict["dat"]["newhistmx"]   # first key locale
        geodf = alldict["geo"]
        spreaddict = alldict["sp"]

    # start the day counter at zero
    reset!(day_ctr, :day)  # return and reset key to 0 :day leftover from prior runs

    locales = locales   # force local scope to be visible in the loop

    ######################
    # simulation loop
    ######################
    sprtime = 0
    trtime = 0
    idxtime = 0
    histtime = 0

    for i = 1:n_days
        inc!(day_ctr, :day)  # increment the simulation day counter
        silent || println("simulation day: ", day_ctr[:day])

        for loc in locales     # @inbounds
            
            locdat = popdat[loc]
            
            density_factor = geodf[geodf[!, :fips] .== loc, :density_factor][]
            for case in runcases
                case(loc, popdat, spreaddict)   
            end
            idxtime += @elapsed begin
                infect_idx = findall(locdat.status .== infectious)
                contactable_idx = findall(locdat.status .!= dead)
            end
            sprtime += @elapsed  spread!(locdat, infect_idx, contactable_idx, spreadcases, 
                spreaddict, density_factor)  # sptime += @elapsed 
            trtime += @elapsed  transition!(locdat, infect_idx, dectree)                       # trtime += @elapsed 

            # r0 displayed every 10 days
            if showr0 && (mod(day_ctr[:day],10) == 0)   # do we ever want to do this by locale -- maybe
                current_r0 = r0_sim(age_dist, popdat, loc, dectree, spreaddict, density_factor)
                println("day $(day_ctr[:day]), locale $loc: rt = $current_r0")
            end

        end

        histtime += @elapsed do_history!(locales, popdat, cumhistmx, newhistmx, agegrp_idx)
        silent || println("Simulation completed for $(day_ctr[:day]) days.")
    end

    #######################

    # simulatio history series for plotting: arrays NOT dataframes
    series = Dict(loc=>Dict(:cum=>cumhistmx[loc], :new=>newhistmx[loc]) for loc in locales)

    # sum agegrps to total for all conditions
    hist_total_agegrps!(series, locales)

    for loc in locales
        add_totinfected_series!(series, loc)
    end

    @show idxtime, sprtime, trtime, histtime

    return alldict, series
end



################################################################################
#  Build and update daily history series
################################################################################

@views function do_history!(locales, popdat, cumhist, newhist, agegrp_idx)
    thisday = day_ctr[:day]
    for loc in locales
        dat = popdat[loc]  # source
        cumdat = cumhist[loc]   # sink
        newdat = newhist[loc]   # sink

        #
        # cumulative data
        #
        @inbounds for age in agegrps
            # get the source data: status
            dat_age = dat[agegrp_idx[loc][age]]
            status_today = countsarr(dat_age.status, unexposed:dead)    # outcomes for thisday

            # get the source data: conditions in (nil, mild, sick, severe)
            filt_infectious = findall(dat_age.status .== infectious)
            if size(filt_infectious, 1) > 0
                sick_today = countsarr(dat_age.cond[filt_infectious], infectious_cases)  # ditto
            else   # there can be days when no one is infected
                sick_today = Dict()
            end

            # insert into sink: cum
            for i in statuses  # 1:4
                cumdat[thisday, map2series[i][age]] = get(status_today, i, 0)
            end

            for i in infectious_cases # 5:8
                cumdat[thisday, map2series[i][age]] = get(sick_today, i, 0)
            end

            for i in all_conds
                if thisday == 1
                    newdat[thisday, map2series[i][age]] = get(status_today, i, 0)
                else  # on all other days
                    newdat[thisday, map2series[i][age]] = (
                        cumdat[thisday, map2series[i][age]] 
                        - cumdat[thisday - 1, map2series[i][age]]
                        )         
                end
            end
        end # for age in agegrps

    end # for loc in locales

end # function


"""
    Count how many times each value of input array is found in an array
    of comparison values.
    Slightly faster than StatsBase: counts (2x) or countmap (4x).
    Requires integer values in a continuous range.
"""
function countsarr(arr, compare_vals)
    vals_range = minimum(compare_vals):maximum(compare_vals)
    ret = zeros(Int, length(vals_range))
    ret = OffsetVector(ret, vals_range)  # enables indexing 5:8, etc.
    @inbounds for i in arr
        ret[i] += 1
    end
    return ret
end


function hist_total_agegrps!(series, locales)
    for loc in locales
        for kind in [:cum, :new]
            for cond in all_conds
                series[loc][kind][:,map2series[cond][totalcol]] = sum(series[loc][kind][:,map2series[cond][agegrps]],dims=2)
            end
        end
    end
end


function review_history(histmx)
    for i in 1:size(histmx, 3)
        println("   *** Day $i ***")
        display(hcat(histmx[:,:,i], [:Unexposed, :Infectious, :Recovered, :Dead, :Nil, :Mild, :Sick, :Severe]))
        print("Press enter or type q and enter to quit> "); resp = chomp(readline())
        if resp == "q"; break; end
    end
end


# a single locale that already has both new and cum series
function add_totinfected_series!(series, locale)
    if !(haskey(series[locale], :cum) && haskey(series[locale], :new))
        error("locale series must contain both :cum and :new series")
        return
    end
    # for new
    @views begin
        n = size(series[locale][:new],1)
        series[locale][:new] = hcat(series[locale][:new], zeros(T_int[], n, 6))
        series[locale][:new][:,map2series.totinfected] = ( (series[locale][:new][:,map2series.unexposed] .< 0 ) .*
                                                          abs.(series[locale][:new][:,map2series.unexposed]) ) 
        # for cum
        series[locale][:cum] = hcat(series[locale][:cum], zeros(T_int[], n, 6))
        cumsum!(series[locale][:cum][:,map2series.totinfected], series[locale][:new][:,map2series.totinfected], dims=1)  
    end
    return
end



#####################################################################################
#  other functions used in simulation
#####################################################################################


function empty_all_caches!()
    # empty tracking queues
    !isempty(spreadq) && (deleteat!(spreadq, 1:length(spreadq)))   
    !isempty(transq) && (deleteat!(transq, 1:length(transq)))   
    !isempty(tntq) && (deleteat!(tntq, 1:length(tntq)))   
    !isempty(r0q) && (deleteat!(r0q, 1:length(r0q)))  
    cleanup_stash(sim_stash) 
end


"""
    optfindall(p, X, maxlen=0)

Returns indices to X where p, a filter, is true.
Filters should be anonymous functions.
For maxlen=0, the length of the temporary vector is length(x).
For maxlen=n, the length of the temporary vector is n.
For maxlen=0.x, the length of temporary vector is 0.x * length(x) and
x should be in (0.0, 1.0).
"""
function optfindall(p, X, maxlen=1)
    if maxlen==1
        out = Vector{Int}(undef, length(X))
    elseif isa(maxlen, Int)
        out = Vector{Int}(undef, maxlen)
    else
        out = Vector{Int}(undef, floor(Int, maxlen * length(X)))
    end
    ind = 0
    @inbounds for (i, x) in pairs(X)
        if p(x)
            out[ind+=1] = i
        end
    end
    resize!(out, ind)
    return out
end


#######################################################################################
#  probability
#######################################################################################


# discrete integer histogram
function bucket(x; vals)
    if isempty(vals)
        vals = range(minimum(x), stop = maximum(x))
    end
    [count(x .== i) for i in vals]
end


# range counts to discretize PDF of continuous outcomes
function histo(x)
    big = ceil(maximum(x))
    bins = Int(big)
    sm = floor(minimum(x))
    ret = zeros(T_int[], bins)
    binbounds = collect(1:bins)
    @inbounds for i = 1:bins
        n = count(x -> i-1 < x <= i,x)
        ret[i] = T_int[](n)
    end
    return ret, binbounds
end


"""
Returns continuous value that represents gamma outcome for a given
approximate center point (scale value of gamma).  We can interpret this
as a funny sort of probability or as a number outcome from a gamma
distributed sample.
1.2 provides a good shape with long tail right and big clump left
"""
function gamma_prob(target; shape=1.0)
    @assert 0.0 <= target <= 99.0 "target must be between 0.0 and 99.0"
    dgamma = Gamma(shape,target)
    pr = rand(dgamma, 1)[1] / 100.0
end


"""
Returns a single number of successes for a
sampled outcome of cnt tries with the input pr of success.
"""
function binomial_one_sample(cnt, pr)::T_int[]
    return rand.(Binomial.(cnt, pr))
end



"""
    categorical_sim(prs::Vector{Float64}, do_assert=true)
    categorical_sim(prs::Vector{Float64}, n::Int, do_assert=true)

Approximates sampling from a categorical distribution.
prs is an array of floats that must sum to 1.0.
do_assert determines if an assert tests this sum. For a single trial, this runs
in 10% of the time of rand(Categorical(prs)). For multiple trials, the 
second method runs in less than 50% of the time.

The second method generates results for n trials. 
The assert test is done only once if do_assert is true.
"""
function categorical_sim(prs, do_assert=true)
    do_assert && @assert isapprox(sum(prs), 1.0)
    x = rand()
    cumpr = 0.0
    i = 0
    for pr in prs
        cumpr += pr
        i += 1
        if x <= cumpr 
            break
        end
    end
    i
end

function categorical_sim(prs, n::Int, do_assert=true)
    do_assert && @assert isapprox(sum(prs), 1.0)
    ret = Vector{Int}(undef, n)
    
    @inbounds for i in 1:n
        ret[i] = categorical_sim(prs, false)
    end
    ret
end



#############################################################
#  other convenience functions
#############################################################


# to access a column of a TypedTable using a variable that holds the symbol
getcol = TypedTables.getproperty


function printsp(xs...)
    for x in xs
       print(x," ")
    end
   println()
end

sparsify!(x, eps=1e-8) = x[abs.(x) .< eps] .= 0.0;


#############################################################
#  experiments
#############################################################


####################################################################################
#   convenience functions for reading and inputting population statistics
#                in the population data matrices
####################################################################################

mutable struct actions
   tests::Array{Array{Int,1},1}  # [[column index, value]]
   cmps::Array{Function, 1}     # must have same number of elements as tests
   todo::Array{Array{Int,1},1}
   setters::Array{Function, 1}  # must have same number of elements as todo
end


mutable struct filts 
   tests::Array{Array{Int,1},1}
   cmps::Array{Function, 1}     # must have same number of elements as tests
end


function make_sick!(dat; cnt, fromage, tocond, tosickday=1)

    @assert size(cnt, 1) == size(fromage, 1)

    filt_unexp = optfindall(==(unexposed), dat.status, 1) # must be unexposed

    for i in 1:size(fromage, 1)  # by target age groups

        filt_age = dat.agegrp[filt_unexp] .== fromage[i] # age of the unexposed
        rowrange = 1:cnt[i]
        filt_all = filt_unexp[filt_age][rowrange]
        cols = [col_status, col_cond, col_sickday]

        dat.status[filt_all] .= infectious
        dat.cond[filt_all] .= tocond
        dat.sickday[filt_all] .= tosickday
    end
end


function incr(a,b)
    a .+= b
end

function setval(a,b)
    a .= b
end

