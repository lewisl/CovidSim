
# pre-allocate large arrays, accessed and modified frequently
# hold complex parameter sets
mutable struct Env
    geodata::Array{Any, 2}
    spreaders::Array{Int64,3} # laglim,4,5
    all_accessible::Array{Int64,3} # laglim,6,5
    numcontacts::Array{Int64,3} # laglim,4,5
    simple_accessible::Array{Int64,2} # 6,5
    numtouched::Array{Int64,2} # laglim,5
    lag_contacts::Array{Int,1} # laglim,
    riskmx::Array{Float64,2} # laglim,5
    contact_factors::Array{Float64,2}  # 4,5 parameters for spread!
    touch_factors::Array{Float64,2}  #  6,5  parameters for spread!
    send_risk_by_lag::Array{Float64,1}  # laglim,  parameters for spread!
    recv_risk_by_age::Array{Float64,1}  # 5,  parameters for spread!
    sd_compliance::Array{Float64,2} # (6,5) social_distancing compliance unexp,recov,nil:severe by age

    # constructor with keyword arguments and type compatible fillins--not suitable as defaults, see initialize_sim_env
    function Env(;          geodata=[0 "" ], # geodata
                            spreaders=zeros(Int,0,0,0),   # semicolon for all keyword (named) arguments)
                            all_accessible=zeros(Int,0,0,0),
                            numcontacts=zeros(Int,0,0,0),
                            simple_accessible=zeros(Int,0,0),
                            numtouched=zeros(Int,0,0),
                            lag_contacts=zeros(Int,laglim),
                            riskmx=zeros(Float64,0,0),
                            contact_factors=zeros(Float64,0,0),
                            touch_factors=zeros(Float64,0,0),
                            send_risk_by_lag=zeros(Float64,laglim),
                            recv_risk_by_age=zeros(Float64,5),
                            sd_compliance=ones(Float64,6,5)        )
        return new(geodata, spreaders, all_accessible, numcontacts, simple_accessible,
                   numtouched, lag_contacts, riskmx, contact_factors,
                   touch_factors, send_risk_by_lag, recv_risk_by_age, sd_compliance)
    end
end

# switches and values used to control cases that run during simulation
const RunControl = Dict{Symbol,Any}()

# control constants
const age_dist = [0.251, 0.271,   0.255,   0.184,   0.039]
const laglim = 25
const lags = 1:laglim   # rows

# geo data: fips,county,city,state,sizecat,pop,density
const fips = 1
const county = 2
const city = 3
const state = 4
const sizecat = 5
const popsize = 6
const density = 7
const anchor = 8
const restrict = 9
const density_fac = 10

# population centers sizecats
const major = 1  # 20
const large = 2  # 50
const medium = 3
const small = 4
const smaller = 5
const rural = 6

# condition_outcome columns and stat series columns
const unexposed = 1
const infectious = 2
const recovered = 3
const dead = 4
const nil = 5
const mild = 6
const sick = 7
const severe = 8
const totinfected = 9
const travelers = 10
const isolated = 11

# columns of history series: first 5 cols are agegrps, 6th is total
const map2series = (unexposed=1:6, infectious=7:12, recovered=13:18, dead=19:24, 
                    nil=25:30, mild=31:36, sick=37:42, severe=43:48, totinfected=49:54)
const total = 6

const conditions = [unexposed, infectious, recovered, dead, nil, mild, sick, severe]
const condnames = Dict(1=>"unexposed", 2=>"infectious", 3=>"recovered", 4=>"dead",
                       5=>"nil", 6=>"mild", 7=>"sick", 8=>"severe", 9=>"totinfected")
const infectious_cases = [nil, mild, sick, severe]
const transition_cases = [recovered, nil, mild, sick, severe, dead]

# transition_prob_rows
const to_recovered = 1
const to_nil = 2
const to_mild = 3
const to_sick = 4
const to_severe = 5
const to_dead = 6

# agegrp channels at dimension 3
const a1 = 1 # 0-19
const a2 = 2 # 20-39
const a3 = 3 # 40-59
const a4 = 4 # 60-79
const a5 = 5 # 80+
const agegrps = 1:5
const ages = length(agegrps)


# traveling constants
const travprobs = [1.0, 2.0, 3.0, 3.0, 0.4] # by age group


####################################################################################
#   simulation runner
####################################################################################


function run_a_sim(n_days, locales; runcases=[], spreadcases=[],showr0 = true, silent=true,
            geofilename="../data/geo2data.csv", 
            dtfilename="../parameters/dec_tree_all_25.csv",
            spfilename="../parameters/spread_params.toml")
#=
    see cases.jl for runcases and spreadcases
=#
    !isempty(spreadq) && (deleteat!(spreadq, 1:length(spreadq)))   # empty it
    !isempty(transq) && (deleteat!(transq, 1:length(transq)))   # empty it
    !isempty(tntq) && (deleteat!(tntq, 1:length(tntq)))   # empty it

    # access input data and pre-allocate storage
    alldict = setup(n_days; geofilename=geofilename, 
                    dectreefilename=dtfilename, spfilename=spfilename)

        dt = alldict["dt"]  # decision trees for transition
        openmx = alldict["dat"]["openmx"]
        cumhistmx = alldict["dat"]["cumhistmx"]
        newhistmx = alldict["dat"]["newhistmx"]
        isolatedmx = alldict["dat"]["isolatedmx"]
        testmx = alldict["dat"]["testmx"]
        geodata = alldict["geo"]
        spread_params = alldict["sp"]
        fips_locs = alldict["fips_locs"]

    env = initialize_sim_env(geodata; spread_params...)

    # start the day counter at zero
    reset!(ctr, :day)  # return and reset key :day leftover from prior runs

    # initial data for building data series of simulation outcomes
    starting_unexposed = reduce(hcat, [grab(unexposed, agegrps, 1, loc, dat=openmx) for loc in locales])
    starting_unexposed = (size(locales,1) == 1 ? Dict(locales[1]=>starting_unexposed) : 
        Dict(locales[i]=>starting_unexposed[i,:] for i in 1:size(locales,1)))

    locales = locales   # force local scope to be visible in the loop
    for i = 1:n_days
        inc!(ctr, :day)  # increment the simulation day counter
        silent || println("simulation day: ", ctr[:day])
        for loc in locales
            density_factor = geodata[geodata[:, fips] .== loc, density_fac][1]
            for case in runcases
                case(loc; opendat=openmx, isodat=isolatedmx, testdat=testmx, env=env)
            end
            spread!(loc, density_factor, dat=openmx, env=env, spreadcases=spreadcases)
            transition!(dt, loc, dat=openmx)   # transition all infectious cases "in the open"
        end
        transition!(dt, dat=isolatedmx)  # transition all infectious cases / locales in isolation
        transition!(dt, dat=testmx) # transition tracked test to stay in sync with openmx
        if showr0 && (mod(ctr[:day],10) == 0)   # do we ever want to do this by locale -- maybe
            current_r0 = sim_r0(env=env, dt=dt)
            println("at day $(ctr[:day]) r0 = $current_r0")
        end

        # println("day $(ctr[:day]) all locales ", keys(isolatedmx))
        do_history!(locales, opendat=openmx, cumhist=cumhistmx, newhist=newhistmx, starting_unexposed=starting_unexposed)
    end
    println("Simulation completed for $(ctr[:day]) days.")

    # "history" series for plotting: NOT dataframes, but arrays
    series = Dict(loc=>Dict(:cum=>make_series(cumhistmx[loc]), :new=>make_series(newhistmx[loc])) for loc in locales)
    for loc in locales
        add_totinfected_series!(series, loc)
    end

    return alldict, env, series
end


function do_history!(locales; opendat, cumhist, newhist, starting_unexposed)
    # capture a snapshot of the end-of-day population matrix
    thisday = ctr[:day]
    if thisday == 1
        for locale in locales
            zerobase = zeros(Int, size(newhist[locale])[1:2])
            zerobase[1,1:5] .+= starting_unexposed[locale]
            zerobase[1,6] = sum(starting_unexposed[locale])

            cumhist[locale][:, 1:5, thisday] = reshape(sum(opendat[locale],dims=1), 8,5)
            cumhist[locale][:, 6, thisday] = sum(cumhist[locale][:, 1:5, thisday], dims=2)
            newhist[locale][:,:,thisday] = cumhist[locale][:,:, thisday] .- zerobase
        end
    else  # on all other days...
        for locale in locales
            cumhist[locale][:,1:5, thisday] = reshape(sum(opendat[locale],dims=1), 8,5)
            @views cumhist[locale][:,6, thisday] = sum(cumhist[locale][:, 1:5, thisday], dims=2) # @views 
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
function make_series(histmx)
    s = zeros(Int, size(histmx,3), prod(size(histmx)[1:2]))
    for i in 1:size(histmx, 3)
        s[i, :] = reduce(vcat,[histmx[j, :, i] for j in 1:size(histmx,1)])'
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
    series[locale][:new] = hcat(series[locale][:new], zeros(Int,n,6))
    series[locale][:new][:,map2series.totinfected] = ( (series[locale][:new][:,map2series.unexposed] .< 0 ) .*
                                                      abs.(series[locale][:new][:,map2series.unexposed]) ) 
    # for cum
    series[locale][:cum] = hcat(series[locale][:cum], zeros(Int,n,6))
    @views cumsum!(series[locale][:cum][:,map2series.totinfected], series[locale][:new][:,map2series.totinfected], dims=1)  
    return
end


function sim_r0(;env=env, dt=dt)
    # captures current population condition 
    pct_unexposed = sum(env.simple_accessible[1,:]) / sum(env.simple_accessible)
    sa_pct = [pct_unexposed,(1-pct_unexposed)/2.0,(1-pct_unexposed)/2.0]   

    # if social_distancing case with split population
    if haskey(spread_stash, :case_cf) || haskey(spread_stash, :case_tf)
        compliance = env.sd_compliance
        cf = spread_stash[:case_cf]; tf = spread_stash[:case_tf]
        r0_comply = r0_sim(compliance = compliance, cf=cf, tf=tf, dt=dt, sa_pct=sa_pct, env=env).r0

        cf = spread_stash[:default_cf]; tf = spread_stash[:default_tf]
        r0_nocomply = r0_sim(compliance=(1.0 .- compliance), cf=cf, tf=tf, dt=dt, sa_pct=sa_pct, env=env).r0

        # this works if all compliance values are the same; approximate otherwise
        current_r0 = round(mean(compliance) * r0_comply + (1.0-mean(compliance)) * r0_nocomply, digits=2)
    else
        cf =  env.contact_factors
        tf = env.touch_factors     
        current_r0 = round(r0_sim(cf=cf, tf=tf, dt=dt, sa_pct=sa_pct, env=env).r0, digits=2)   
    end
    return current_r0
end


function initialize_sim_env(geodata; contact_factors, touch_factors, send_risk, recv_risk)
 
    # initialize the simulation Env
    ret =   Env(geodata=geodata,
                spreaders=zeros(Int64,laglim,4,5),
                all_accessible=zeros(Int64,laglim,6,5),
                numcontacts=zeros(Int64,laglim,4,5),
                simple_accessible=zeros(Int64,6,5),
                numtouched=zeros(Int64,laglim,5),
                lag_contacts=zeros(Int64,laglim),
                riskmx = zeros(Float64,laglim,5),
                contact_factors = contact_factors,
                touch_factors = touch_factors,
                send_risk_by_lag = send_risk,
                recv_risk_by_age = recv_risk,
                sd_compliance = zeros(6,5))
    ret.riskmx = send_risk_by_recv_risk(ret.send_risk_by_lag, ret.recv_risk_by_age)

    # contact_factors and touch_factors look like:
    #=
        contact_factors = 
                [ 1    1.8    1.8     1.5     1.0;    # nil
                  1    1.7    1.7     1.4     0.9;    # mild
                0.7    1.0    1.0     0.7     0.5;   # sick
                0.5    0.8    0.8     0.5     0.3]  # severe

      # agegrp    1     2      3       4       5

        touch_factors = 
                [.55    .62     .58     .4    .35;    # unexposed
                 .55    .62     .58     .4    .35;    # recovered
                 .55    .62     .58     .4    .35;    # nil
                 .55    .6      .5      .35   .28;    # mild
                 .28   .35      .28     .18   .18;    # sick
                 .18   .18      .18     .18   .18]   # severe
    =#

    return ret
end


#######################################################################################
#  probability
#######################################################################################


# discrete integer histogram
function bucket(x; vals=[])
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
    ret = zeros(Int, bins)
    binbounds = collect(1:bins)
    for i = 1:bins
        n = count(x -> i-1 < x <= i,x)
        ret[i] = n
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
function binomial_one_sample(cnt, pr)::Int
    return rand.(Binomial.(cnt, pr))
end


function categorical_sample(probvec, trials)
    x = rand(Categorical(probvec), trials)
end


####################################################################################
#   convenience functions for reading and inputting population statistics
#                in the simulation data matrices
####################################################################################

"""
    function grab(condition, agegrp, lag, locale; dat=openmx)

Inputs condition, agegrp and lag can be single or multiple (array or range).
Only one locale can be accessed. Caller should loop over locales to retrieve
data from multiple locales.

The returned array has dimensions (lag, condition, agegrp).
"""
function grab(condition, agegrp, lag, locale; dat=openmx)
    @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single Int or NamedTuple"
    return dat[locale][lag, condition, agegrp]
end


"""
    function input!(val, condition, agegrp, lag, locale; dat=openmx)

Input val can be an array. Its dimensions must match the inputs for lag, condition, agegrp.
Only one locale can be provided.

ex: if size(val) is (25, 4, 5) then length(lag) must = 25, length(condition) must = 4, and length(agegrp) must = 5.

Inputs overwrite existing data at the referenced location of the target population matrix dat.
"""
function input!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single Int or NamedTuple"
    current = grab(condition, agegrp, lag, locale; dat=dat)
    dat[locale][lag, condition, agegrp] = val
end


"""
    function plus!(val, condition, agegrp, lag, locale; dat=openmx)

Input val can be an array. Its dimensions must match the inputs for lag, condition, agegrp.
Only one locale can be provided.

ex: if size(val) is (25, 4, 5) then length(lag) must = 25, length(condition) must = 4, and length(agegrp) must = 5.

Inputs are added to the existing data at the referenced location of the target population matrix dat.
"""
function plus!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single Int or NamedTuple"
    dat[locale][lag, condition, agegrp] += val
end


"""
    function minus!(val, condition, agegrp, lag, locale; dat=openmx)

Input val can be an array. Its dimensions must match the inputs for lag, condition, agegrp.
Only one locale can be provided.

ex: if size(val) is (25, 4, 5) then length(lag) must = 25, length(condition) must = 4, and length(agegrp) must = 5.

If subtraction from the existing data would result in negative values at the referenced locations of the target population matrix dat,
an error will be raised. The population matrix must contain positive integer values.
"""
function minus!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single Int or NamedTuple"
    current = grab(condition, agegrp, lag, locale; dat=dat)
    @assert sum(val) <= sum(current) "subtracting > than existing: day $(ctr[:day]) loc $locale lag $lag cond $condition agegrp $agegrp"
    dat[locale][lag, condition, agegrp] -= val
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
