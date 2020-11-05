
####################################################################################
#   simulation runner
####################################################################################


function run_a_sim(n_days, locales; runcases=[], spreadcases=[], showr0 = true, silent=true, set_int_type=Int64,
            geofilename="../data/geo2data.csv", 
            dtfilename="../parameters/dec_tree_all_25.yml",
            spfilename="../parameters/spread_params_ilm.yml")

    empty_all_caches!() # from previous runs

    T_int[] = set_int_type # update the global type of ints with the input value

    # access input data and pre-allocate storage
    alldict = setup(n_days, locales; geofilename=geofilename, 
                    dectreefilename=dtfilename, spfilename=spfilename)

        dt_dict = alldict["dt_dict"]  # decision trees for transition
        popdat = alldict["dat"]["popdat"]
        agegrp_idx = alldict["dat"]["agegrp_idx"]
        cumhistmx = alldict["dat"]["cumhistmx"]
        newhistmx = alldict["dat"]["newhistmx"]
        geodf = alldict["geo"]
        spread_params = alldict["sp"]

        env = initialize_sim_env(geodf; spread_params...)

    # start the day counter at zero
    reset!(ctr, :day)  # return and reset key to 0 :day leftover from prior runs

    locales = locales   # force local scope to be visible in the loop


    # check age distribution
    # dat = popdat[first(locales)]
    # agedist = countmap(dat[:,cpop_agegrp])
    # println(agedist)
    # println(sum(values(agedist)))

    ######################
    # simulation loop
    ######################
    sprtime = 0
    trtime = 0
    simtime = 0
    histtime = 0
    for i = 1:n_days
        inc!(ctr, :day)  # increment the simulation day counter
        silent || println("simulation day: ", ctr[:day])
        for loc in locales     # @inbounds

            density_factor = geodf[geodf[:fips] .== loc, :density_factor][]
            for case in runcases
                # case(loc, popdat, isolatedmx, testmx, env)   
                case(loc, popdat, [], [], env)   
            end
            sprtime += @elapsed  spread!(popdat, loc, spreadcases, env, density_factor)  # sptime += @elapsed 
            trtime += @elapsed  transition!(popdat, loc, dt_dict)                       # trtime += @elapsed 

            # r0 displayed every 10 days
            if showr0 && (mod(ctr[:day],10) == 0)   # do we ever want to do this by locale -- maybe
                current_r0 = r0_sim(age_dist, popdat, loc, dt_dict, env, density_factor)
                println("day $(ctr[:day]), locale $loc: rt = $current_r0")
            end

        end

        histtime += @elapsed do_history!(locales, popdat, cumhistmx, newhistmx, agegrp_idx)

    end
    silent || println("Simulation completed for $(ctr[:day]) days.")
    #######################

    # simulatio history series for plotting: arrays NOT dataframes
    series = Dict(loc=>Dict(:cum=>cumhistmx[loc], :new=>newhistmx[loc]) for loc in locales)

    # sum agegrps to total for all conditions
    hist_total_agegrps!(series, locales)

    for loc in locales
        add_totinfected_series!(series, loc)
    end

    @show sprtime, trtime, histtime

    return alldict, env, series
end



################################################################################
#  Build and update daily history series
################################################################################


# const map2series = (unexposed=1:6, infectious=7:12, recovered=13:18, dead=19:24, 
#                     nil=25:30, mild=31:36, sick=37:42, severe=43:48, totinfected=49:54)

function do_history!(locales, popdat, cumhist, newhist, agegrp_idx)
    thisday = ctr[:day]
    for loc in locales
        dat = popdat[loc]  # source
        cumdat = cumhist[loc]   # sink
        newdat = newhist[loc]   # sink

        #
        # cumulative data
        #
        @views for age in agegrps
            # get the source data: status
            dat_age = dat[agegrp_idx[loc][age]]
            status_today = countsarr(dat_age.status, unexposed:dead)    # outcomes for thisday

            # get the source data: conditions in (nil, mild, sick, severe)
            filt_infectious = findall(dat_age.status .== infectious)
            if size(filt_infectious, 1) > 0
                sick_today = countsarr(dat_age.cond[filt_infectious], nil:severe)  # ditto
            else   # there can be days when no one is infected
                sick_today = Dict()
            end

            # insert into sink: cum
            for i in unexposed:dead  # 1:4
                cumdat[thisday, map2series[i][age]] = get(status_today, i, 0)
            end

            for i in nil:severe # 5:8
                cumdat[thisday, map2series[i][age]] = get(sick_today, i, 0)
            end

            for i in [unexposed, infectious, recovered, dead, nil, mild, sick, severe]
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
    Allows counting appearances of all values appearing in an array.
    Slightly faster than StatsBase: counts (2x) or countmap (4x).
    Requires integer values in a continuous range.
"""
function countsarr(input, vals)
    ret = zeros(Int, size(vals,1))
    ret = OffsetVector(ret, vals)  # enables indexing 5:8, etc.
    @inbounds for i in input
        ret[i] += 1
    end
    return ret
end


function hist_total_agegrps!(series, locales)
    for loc in locales
        for kind in [:cum, :new]
            for cond in conditions
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
Returns an array of row indices that satisfy the filters a named tuple of related columns.
The related columns can be a TypedTable or a columntable, as created by
Tables.columntable().

Filters must be an array of (symbol, comparison function, value) where symbol is a column reference 
in the named tuple.
   
example: 

```julia
(:status, ==, 5)
```      

"""
function sq_query(dat, filters::Array{Popquery, 1})
    filts = copy(filters)
    # ret = Array{Int,1}
    # q = popfirst!(filts)
    # idx1 = @elapsed ret = map(x -> q.op(x, q.val), getproperty(dat, q.col))   # findall()
    # @show size(ret)

    # idxloop = @elapsed for q in filts
    #     ret = map(x -> q.op(x, q.val), getproperty(dat, q.col)[ret])   # findall()
    #     @show size(ret)
    # end

    ret = trues(size(dat,1))

    # ret .& X for X in [map(x -> q.op(x, q.val), getproperty(dat, q.col)) for q in filts]    # reduce(X -> X .& ret,    ) 

    # reduce(&, [BitArray(map(x -> q.op(x, q.val),getproperty(locdat, q.col))) for q in filts]; init = trues(size(locdat,1)))
    # there is a place for a foreach() somewhere in that mess

    # this works but is only a bit faster than multiple findall's
    for q in filters
        ret .&= map(x -> q.op(x, q.val), getproperty(dat, q.col))
    end

    return findall(ret)
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


function categorical_sample(probvec, trials)::Array{T_int[],1}
    x = rand(Categorical(probvec), trials)
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


function update!(dat, cnt, actions::actions)  

    filt = falses(size(dat,1))   # put this in the env so only do once

    filt[:] = actions.cmps[1].(dat[:, actions.tests[1][1]], actions.tests[1][2])
    for i in 2:length(actions.tests)
        filt[:] .&= actions.cmps[i].(dat[:, actions.tests[i][1]], actions.tests[i][2])
    end

    rowsel = cnt == 0 ? (:) : 1:cnt  # (:) selects all matches

    for i = 1:length(actions.todo)
        sel = view(dat, filt, actions.todo[i][1])

        try
            sel[rowsel] =  actions.setters[i](sel[rowsel], actions.todo[i][2])
        catch
            @warn("no match or too many updates")
        end
    end

end


function make_sick!(dat; cnt, fromage, tocond, tolag=1)

    @assert size(cnt, 1) == size(fromage, 1)

    filt_unexp = findall(dat.status .== unexposed) # must be unexposed

    for i in 1:size(fromage, 1)  # by target age groups

        filt_age = dat.agegrp[filt_unexp] .== fromage[i] # age of the unexposed
        rowrange = 1:cnt[i]
        filt_all = filt_unexp[filt_age][rowrange]
        cols = [cpop_status, cpop_cond, cpop_lag]

        dat.status[filt_all] .= infectious
        dat.cond[filt_all] .= tocond
        dat.lag[filt_all] .= tolag
    end
end


# function make_sick!(dat; cnt, fromage, tocond, tolag=1)

#     @assert size(cnt, 1) == size(fromage, 1)

#     filt_unexp = findall((dat[:,cpop_status] .== unexposed)) # must be unexposed

#     for i in 1:size(fromage, 1)  # by target age groups

#         filt_age = dat[filt_unexp, cpop_agegrp] .== fromage[i] # age of the unexposed
#         rowrange = 1:cnt[i]
#         filt_all = filt_unexp[filt_age][rowrange]
#         cols = [cpop_status, cpop_cond, cpop_lag]

#         dat[filt_all, cols] .= [infectious tocond tolag]
#     end
# end



function change_sick!(dat; cnt, fromcond, fromage, fromlag, tests=[], tocond)
    cs_actions = actions([[cpop_cond, fromcond], [cpop_agegrp, fromage],
                              [cpop_lag, fromlag], [cpop_status, infectious]], # tests
                      [==, ==, ==, ==],  # cmps
                      [[cpop_cond, tocond]], # todo
                      [setval])  # setters

    update!(dat, cnt, cs_actions)
end


# this isn't going to work
function bump_sick!(dat; cnt, fromcond, fromage, fromlag, tests=[])


    bs_actions = actions([[cpop_status, infectious]], # tests
                      [==],  # cmps
                      [[cpop_lag, 1]], # todo
                      [incr])  # setters

    if fromcond != 0
        push!(bs_actions.tests, [cpop_cond, fromcond])
        push!(bs_actions.cmps, ==)
    end
    if fromage != 0
        push!(bs_actions.tests, [cpop_agegrp, fromage])
        push!(bs_actions.cmps, ==)
    end
    if fromlag != 0
        push!(bs_actions.tests, [cpop_lag, fromlag])
        push!(bs_actions.cmps, ==)
    end


    update!(dat, cnt, bs_actions)
end

# this isn't going to work
function make_dead!(dat; cnt, fromage, fromlag, fromcond, tests=[])
    md_actions = actions([[cpop_cond, fromcond], [cpop_agegrp, fromage],
                                [cpop_lag, fromlag], [cpop_status, infectious]], # tests
                          [==, ==, ==, ==],  # cmps
                          [[cpop_status, dead]], # todo
                          [setval])  # setters

    update!(dat, cnt, md_actions)
end

# this isn't going to work
function make_recovered!(dat; cnt, fromage, fromlag, fromcond, tests=[])
    mr_actions = actions([[cpop_cond, fromcond], [cpop_agegrp, fromage],
                                [cpop_lag, fromlag], [cpop_status, infectious]], # tests
                          [==, ==, ==, ==],  # cmps
                          [[cpop_status, recovered]], # todo
                          [setval])  # setters

    update!(dat, cnt, mr_actions)
end


function incr(a,b)
    a .+= b
end

function setval(a,b)
    a .= b
end

