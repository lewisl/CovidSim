



####################################################################################
#   simulation runner
####################################################################################


function run_a_sim(n_days, locales; runcases=[], spreadcases=[], showr0 = true, silent=true, set_int_type=Int64,
            geofilename="../data/geo2data.csv", 
            dtfilename="../parameters/dec_tree_all_25.yml",
            spfilename="../parameters/spread_params_grp.yml")
#=
    see cases.jl for runcases and spreadcases
=#

    empty_all_caches!() # from previous runs

    T_int[] = set_int_type # update the global type of ints with the input value

    # access input data and pre-allocate storage
    alldict = setup(n_days; geofilename=geofilename, 
                    dectreefilename=dtfilename, spfilename=spfilename)

        dt_dict = alldict["dt"]  # decision trees for transition
        openmx = alldict["dat"]["openmx"]
        cumhistmx = alldict["dat"]["cumhistmx"]
        newhistmx = alldict["dat"]["newhistmx"]
        isolatedmx = alldict["dat"]["isolatedmx"]
        testmx = alldict["dat"]["testmx"]
        geodata = alldict["geo"]
        spread_params = alldict["sp"]
        fips_locs = alldict["fips_locs"]

        env = initialize_sim_env(geodata; spread_params...)
        density_factors = Dict(loc => 
            geodata[geodata[:, fips] .== loc, density_fac][1] for loc in locales)

        # initial data for building data series of simulation outcomes
        starting_unexposed = reduce(hcat, [grab(unexposed, agegrps, 1, loc, openmx) for loc in locales])
        starting_unexposed = (size(locales,1) == 1 ? Dict(locales[1]=>starting_unexposed) : 
            Dict(locales[i]=>starting_unexposed[i,:] for i in 1:size(locales,1)))

    # start the day counter at zero
    reset!(ctr, :day)  # return and reset key :day leftover from prior runs

    locales = locales   # force local scope to be visible in the loop


    # check age distribution
    # dat = openmx[first(locales)]
    # println(dat[1, unexposed, :])
    # println(sum(dat[1, unexposed, :]))


    ######################
    # simulation loop
    ######################
    for i = 1:n_days
        inc!(ctr, :day)  # increment the simulation day counter
        silent || println("simulation day: ", ctr[:day])
        @inbounds for loc in locales
            density_factor = density_factors[loc]
            for case in runcases
                case(loc, openmx, isolatedmx, testmx, env)   
            end
            spread!(loc, spreadcases, openmx, env,  density_factor)
            transition!(openmx, dt_dict, loc)   # transition infectious cases "in the open"
        end
        transition!(isolatedmx, dt_dict)  # transition infectious cases isolation
        transition!(testmx, dt_dict) # transition infectious cases in test and trace

        # r0 displayed every 10 days
        if showr0 && (mod(ctr[:day],10) == 0)   # do we ever want to do this by locale -- maybe
            current_r0 = sim_r0(env, dt_dict)
            println("at day $(ctr[:day]) r0 = $current_r0")
        end

        # println("day $(ctr[:day]) all locales ", keys(isolatedmx))
        do_history!(locales, opendat=openmx, cumhist=cumhistmx, newhist=newhistmx, 
            starting_unexposed=starting_unexposed)
    end
    silent || println("Simulation completed for $(ctr[:day]) days.")
    #######################

    # "history" series for plotting: NOT dataframes, but arrays
    series = Dict(loc=>Dict(:cum=>make_series(cumhistmx[loc]), :new=>make_series(newhistmx[loc])) for loc in locales)
    for loc in locales
        add_totinfected_series!(series, loc)
    end

    return alldict, env, series
end



################################################################################
#  Build and update daily history series
################################################################################

function do_history!(locales; opendat, cumhist, newhist, starting_unexposed)
    # capture a snapshot of the end-of-day population matrix
    thisday = ctr[:day]
    if thisday == 1
        @inbounds for locale in locales
            zerobase = zeros(T_int[], size(newhist[locale])[1:2])
            zerobase[1,agegrps] .+= starting_unexposed[locale]
            zerobase[1,totalcol] = sum(starting_unexposed[locale])

            @views cumhist[locale][:, agegrps, thisday] = reshape(sum(opendat[locale],dims=1), n_conditions, n_agegrps)
            @views cumhist[locale][:, totalcol, thisday] = sum(cumhist[locale][:, agegrps, thisday], dims=2)
            newhist[locale][:,:,thisday] = cumhist[locale][:,:, thisday] .- zerobase
        end
    else  # on all other days...
        @inbounds for locale in locales
            cumhist[locale][:,agegrps, thisday] = reshape(sum(opendat[locale],dims=1), n_conditions, n_agegrps)
            @views cumhist[locale][:,totalcol, thisday] = sum(cumhist[locale][:, agegrps, thisday], dims=2) # @views 
            @views newhist[locale][:,:,thisday] = cumhist[locale][:,:,thisday] .- cumhist[locale][:,:,thisday-1] # @views 
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
function make_series(histmx)  # this function just changes the shape and indexing of the data
    s = zeros(T_int[], size(histmx,3), prod(size(histmx)[1:2])) # days, conds * (agegrps + 1)
    for day in 1:size(histmx, 3)
        @views s[day, :] = reduce(vcat,[histmx[cond, :, day] for cond in 1:size(histmx,1)])'
    end
    return s
end

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


####################################################################################
#   convenience functions for reading and inputting population statistics
#                in the population data matrices
####################################################################################

"""
    function grab(condition, agegrp, lag, locale, dat)

Inputs condition, agegrp and lag can be single or multiple (array or range).
Only one locale can be accessed. Caller should loop over locales to retrieve
data from multiple locales.

The returned array has dimensions (lag, condition, agegrp).
"""
function grab(condition, agegrp, lag, locale, dat)
    # @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single Int or NamedTuple"

    return dat[locale][lag, condition, agegrp]
end


"""
    function input!(val, condition, agegrp, lag, locale, dat)

Input val can be an array. Its dimensions must match the inputs for lag, condition, agegrp.
Only one locale can be provided.

ex: if size(val) is (25, 4, 5) then length(lag) must = 25, length(condition) must = 4, and length(agegrp) must = 5.

Inputs overwrite existing data at the referenced location of the target population matrix dat.
"""
function input!(val, condition, agegrp, lag, locale, dat)
    # @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single Int or NamedTuple"
    current = grab(condition, agegrp, lag, locale, dat)
    dat[locale][lag, condition, agegrp] = val
end


"""
    function plus!(val, condition, agegrp, lag, locale, dat)

Input val can be an array. Its dimensions must match the inputs for lag, condition, agegrp.
Only one locale can be provided.

ex: if size(val) is (25, 4, 5) then length(lag) must = 25, length(condition) must = 4, and length(agegrp) must = 5.

Inputs are added to the existing data at the referenced location of the target population matrix dat.
"""
function plus!(val, condition, agegrp, lag, locale, dat) 
    # @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single Int or NamedTuple"

    dat[locale][lag, condition, agegrp] += val    # T(val)
end


"""
    function minus!(val, condition, agegrp, lag, locale, dat)

Input val can be an array. Its dimensions must match the inputs for lag, condition, agegrp.
Only one locale can be provided.

ex: if size(val) is (25, 4, 5) then length(lag) must = 25, length(condition) must = 4, and length(agegrp) must = 5.

If subtraction from the existing data would result in negative values at the referenced locations of the target population matrix dat,
an error will be raised. The population matrix must contain positive integer values.
"""
function minus!(val, condition, agegrp, lag, locale, dat)
    # @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single Int or NamedTuple"
    current = grab(condition, agegrp, lag, locale, dat)
    @assert sum(val) <= sum(current) "subtracting > than existing: day $(ctr[:day]) loc $locale lag $lag cond $condition agegrp $agegrp"
    dat[locale][lag, condition, agegrp] -= val
end


#####################################################################################
#  other functions used in simulation
#####################################################################################

# returns a single r0 value
function sim_r0(env, dt_dict)  # named args must be provided by caller
    # captures current population condition 
    pct_unexposed = sum(env.simple_accessible[1,:]) / sum(env.simple_accessible)
    sa_pct = [pct_unexposed,(1-pct_unexposed)/2.0,(1-pct_unexposed)/2.0]   

    # if social_distancing case with split population
    if haskey(spread_stash, :case_cf) || haskey(spread_stash, :case_tf)
        compliance = env.sd_compliance
        cf = spread_stash[:case_cf]; tf = spread_stash[:case_tf]
        r0_comply = r0_sim(compliance = compliance, cf=cf, tf=tf, dt_dict=dt_dict, sa_pct=sa_pct, env=env).r0

        cf = spread_stash[:default_cf]; tf = spread_stash[:default_tf]
        r0_nocomply = r0_sim(compliance=(1.0 .- compliance), cf=cf, tf=tf, dt_dict=dt_dict, sa_pct=sa_pct, env=env).r0

        # this works if all compliance values are the same; approximate otherwise
        current_r0 = round(mean(compliance) * r0_comply + (1.0-mean(compliance)) * r0_nocomply, digits=2)
    else
        cf =  env.contact_factors
        tf = env.touch_factors     
        current_r0 = round(r0_sim(cf=cf, tf=tf, dt_dict=dt_dict, sa_pct=sa_pct, env=env).r0, digits=2)   
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

"""
    categorical_sim(prs::Vector{Float64}, do_assert=true)
    categorical_sim(prs::Vector{Float64}, n::Int, do_assert=true)

Approximates sampling from a categorical distribution.
prs is an array of floats that must sum to 1.0.
do_assert determines if an assert tests this sum. For a single trial, this runs
in 10% of the time of rand(Categorical(prs)). For multiple trial, the 
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


function printsp(xs...)
    for x in xs
       print(x," ")
    end
   println()
end

sparsify!(x, eps=1e-8) = x[abs.(x) .< eps] .= 0.0;
