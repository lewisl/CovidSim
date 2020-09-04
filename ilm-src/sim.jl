
####################################################################################
#   simulation runner
####################################################################################


function run_a_sim(n_days, locales; runcases=[], spreadcases=[], showr0 = true, silent=true, set_int_type=Int64,
            geofilename="../data/geo2data.csv", 
            dtfilename="../parameters/dec_tree_all_25.yml",
            spfilename="../parameters/spread_params.yml")

    empty_all_caches!() # from previous runs

    T_int[] = set_int_type # update the global type of ints with the input value

    # access input data and pre-allocate storage
    alldict = setup(n_days, locales; geofilename=geofilename, 
                    dectreefilename=dtfilename, spfilename=spfilename)

        dt = alldict["dt"]  # decision trees for transition
        all_decpoints = alldict["decpoints"]
        openmx = alldict["dat"]["openmx"]
        agegrp_idx = alldict["dat"]["agegrp_idx"]
        cumhistmx = alldict["dat"]["cumhistmx"]
        newhistmx = alldict["dat"]["newhistmx"]
        # isolatedmx = alldict["dat"]["isolatedmx"]
        # testmx = alldict["dat"]["testmx"]
        geodf = alldict["geo"]
        spread_params = alldict["sp"]

        env = initialize_sim_env(geodf; spread_params...)

        # initial data for building data series of simulation outcomes
        starting_unexposed = [sum(openmx[loc][:,cpop_status]) for loc in locales]
        # starting_unexposed = reduce(hcat, [grab(unexposed, agegrps, 1, loc, openmx) for loc in locales])
        starting_unexposed = (size(locales,1) == 1 ? Dict(locales[1]=>starting_unexposed[1]) : 
            Dict(locales[i]=>starting_unexposed[i] for i in 1:size(locales,1)))

    # start the day counter at zero
    reset!(ctr, :day)  # return and reset key to 0 :day leftover from prior runs

    locales = locales   # force local scope to be visible in the loop

    ######################
    # simulation loop
    ######################
    sptime = 0
    trtime = 0
    for i = 1:n_days
        inc!(ctr, :day)  # increment the simulation day counter
        silent || println("simulation day: ", ctr[:day])
        @inbounds for loc in locales
            density_factor = geodf[geodf[:fips] .== loc, :density_factor][]
            for case in runcases
                # case(loc, openmx, isolatedmx, testmx, env)   
                case(loc, openmx, [], [], env)   
            end
            if isempty(spreadcases)
                sptime += @elapsed spread!(openmx, loc, env,  density_factor)
            else
                sptime += @elapsed spread!(openmx, loc, spreadcases, env,  density_factor)
            end
            trtime += @elapsed transition!(dt, all_decpoints, loc, openmx, agegrp_idx)   # transition infectious cases "in the open"
        end
        # transition!(dt, all_decpoints, isolatedmx)  # transition infectious cases isolation
        # transition!(dt, all_decpoints, testmx) # transition infectious cases in test and trace

        # r0 displayed every 10 days
        if showr0 && (mod(ctr[:day],10) == 0)   # do we ever want to do this by locale -- maybe
            current_r0 = sim_r0(env, dt, all_decpoints)
            println("at day $(ctr[:day]) r0 = $current_r0")
        end

        do_history!(locales, opendat=openmx, cumhist=cumhistmx, newhist=newhistmx, 
            starting_unexposed=starting_unexposed)

    end
    silent || println("Simulation completed for $(ctr[:day]) days.")
    #######################

    # "history" series for plotting: NOT dataframes, but arrays
    series = Dict(loc=>Dict(:cum=>cumhistmx[loc], :new=>newhistmx[loc]) for loc in locales)

    # for loc in locales
    #     add_totinfected_series!(series, loc)
    # end

    @show sptime, trtime

    return alldict, env, series
end



################################################################################
#  Build and update daily history series
################################################################################


# const map2series = (unexposed=1:6, infectious=7:12, recovered=13:18, dead=19:24, 
#                     nil=25:30, mild=31:36, sick=37:42, severe=43:48, totinfected=49:54)

function do_history!(locales; opendat, cumhist, newhist, starting_unexposed)
    thisday = ctr[:day]
    for loc in locales
        dat = opendat[loc]  # source
        cumdat = cumhist[loc]   # sink
        newdat = newhist[loc]   # sink

        statusday = countmap(dat[:, cpop_status])    # outcomes for thisday
        sickday = countmap(dat[dat[:,cpop_status] .== infectious, cpop_cond])  # ditto

        for i in unexposed:dead  # 1:4
            cumdat[thisday, map2series[i][totalcol]] = get(statusday, i, 0)
        end

        for i in nil:severe # 5:8
            cumdat[thisday, map2series[i][totalcol]] = get(sickday, i, 0)
        end

        # new today
            # dead = today:dead - yesterday:dead
            # recovered = today:recovered - yesterday:recovered
            # infectious = today:infectious - yesterday:infectious + newdead + newrecovered
            # little bit harder for nil, mild, sick, severe
            # TODO getting it working for each agegrp
        for i in [dead, recovered]
            if thisday == 1
                newdat[thisday, map2series[i][totalcol]] = get(statusday, i, 0)
            else  # on all other days
                newdat[thisday, map2series[i][totalcol]] = (
                    cumdat[thisday, map2series[i][totalcol]] 
                    - cumdat[thisday - 1, map2series[i][totalcol]]
                    )         
            end
        end
        if thisday == 1
            newdat[thisday, map2series[infectious][totalcol]] = get(statusday, infectious, 0)
        else
            newdat[thisday, map2series[infectious][totalcol]] = (
                cumdat[thisday, map2series[infectious][totalcol]]
                - cumdat[thisday - 1, map2series[infectious][totalcol]]
                + newdat[thisday, map2series[recovered][totalcol]]
                + newdat[thisday, map2series[dead][totalcol]]
                )
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


# a single locale, either cumulative or new
# function make_series(histmx)
#     s = zeros(T_int[], size(histmx,3), prod(size(histmx)[1:2]))
#     for i in 1:size(histmx, 3)
#         @views s[i, :] = reduce(vcat,[histmx[j, :, i] for j in 1:size(histmx,1)])'
#     end
#     return s
# end

# a single locale that already has both new and cum series
function add_totinfected_series!(series, locale)
    if !(haskey(series[locale], :cum) && haskey(series[locale], :new))
        error("locale series must contain both :cum and :new series")
        return
    end
    # for new
    n = size(series[locale][:new],1)
    series[locale][:new] = hcat(series[locale][:new], zeros(T_int[], n, 6))
    series[locale][:new][:,map2series.totinfected] = ( (series[locale][:new][:,map2series.unexposed] .< 0 ) .*
                                                      abs.(series[locale][:new][:,map2series.unexposed]) ) 
    # for cum
    series[locale][:cum] = hcat(series[locale][:cum], zeros(T_int[], n, 6))
    @views cumsum!(series[locale][:cum][:,map2series.totinfected], series[locale][:new][:,map2series.totinfected], dims=1)  
    return
end



#####################################################################################
#  other functions used in simulation
#####################################################################################

# returns a single r0 value
function sim_r0(env, dt, all_decpoints)  # named args must be provided by caller
    # captures current population condition 
    pct_unexposed = sum(env.simple_accessible[1,:]) / sum(env.simple_accessible)
    sa_pct = [pct_unexposed,(1-pct_unexposed)/2.0,(1-pct_unexposed)/2.0]   

    # if social_distancing case with split population
    if haskey(spread_stash, :case_cf) || haskey(spread_stash, :case_tf)
        compliance = env.sd_compliance
        cf = spread_stash[:case_cf]; tf = spread_stash[:case_tf]
        r0_comply = r0_sim(compliance = compliance, cf=cf, tf=tf, dt=dt, decpoints=all_decpoints, sa_pct=sa_pct, env=env).r0

        cf = spread_stash[:default_cf]; tf = spread_stash[:default_tf]
        r0_nocomply = r0_sim(compliance=(1.0 .- compliance), cf=cf, tf=tf, dt=dt, decpoints=all_decpoints,
                             sa_pct=sa_pct, env=env).r0

        # this works if all compliance values are the same; approximate otherwise
        current_r0 = round(mean(compliance) * r0_comply + (1.0-mean(compliance)) * r0_nocomply, digits=2)
    else
        cf =  env.contact_factors
        tf = env.touch_factors     
        current_r0 = round(r0_sim(cf=cf, tf=tf, dt=dt, decpoints=all_decpoints, sa_pct=sa_pct, env=env).r0, digits=2)   
    end
    return current_r0
end


function empty_all_caches!()
    # empty tracking queues
    !isempty(spreadq) && (deleteat!(spreadq, 1:length(spreadq)))   
    !isempty(transq) && (deleteat!(transq, 1:length(transq)))   
    !isempty(tntq) && (deleteat!(tntq, 1:length(tntq)))   
    !isempty(r0q) && (deleteat!(r0q, 1:length(r0q)))  
    cleanup_stash(spread_stash) 
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

    filt_unexp = findall((dat[:,cpop_status] .== unexposed)) # must be unexposed

    for i in 1:size(fromage, 1)  # by target age groups

        filt_age = dat[filt_unexp, cpop_agegrp] .== fromage[i] # age of the unexposed
        rowrange = 1:cnt[i]
        filt_all = filt_unexp[filt_age][rowrange]
        cols = [cpop_status, cpop_cond, cpop_lag]

        dat[filt_all, cols] .= [infectious tocond tolag]
    end


end


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

