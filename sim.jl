
# pre-allocate large arrays, accessed and modified frequently
# hold complex parameter sets
mutable struct Env
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
    function Env(;          spreaders=zeros(Int,0,0,0),   # semicolon for all keyword (named) arguments)
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
        return new(spreaders, all_accessible, numcontacts, simple_accessible,
                   numtouched, lag_contacts, riskmx, contact_factors,
                   touch_factors, send_risk_by_lag, recv_risk_by_age, sd_compliance)
    end

end

#  stash for temporary values changed during simulation cases
#      to change just once and then get the originals back
#      it is the users responsibility to get rid of stuff
const stash = Dict{Symbol, Array}()

# control constants
const age_dist = [0.251, 0.271,   0.255,   0.184,   0.039]
const laglim = 25
const lags = 1:laglim   # rows

# geo data: id,fips,county,city,state,sizecat,pop,density
const id = 1
const fips = 2
const county = 3
const city = 4
const state = 5
const sizecat = 6
const popsize = 7
const density = 8
const anchor = 9
const restrict = 10
const density_fac = 11

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
const travelers = 9
const isolated = 10

const conditions = [unexposed, infectious, recovered, dead, nil, mild, sick, severe]
const condnames = Dict(1=>"unexposed",2=>"infectious",3=>"recovered", 4=>"dead",
                       5=>"nil",6=>"mild",7=>"sick",8=>"severe")
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


"""
Map a condition index from rows in the data matrix to
indices for the transition probabilities:

```
               unexposed  infectious  recovered  dead   nil  mild  sick  severe
data rows         1          2            3        4     5     6     7     8
transition pr    -1         -1            1        6     2     3     4     5
```

Tranition pr indices that return -1 are not used and will raise an error.

- Use with text literal in code as mapcond2tran.nil => 2
- Use with variables that stand for the data rows as mapcond2tran[nil] => 2
"""
const mapcond2tran = (unexposed=-1, infectious=-1, recovered=1, dead=6, nil=2, mild=3, sick=4, severe=5)


####################################################################################
#   simulation runner
####################################################################################


function run_a_sim(geofilename, n_days, locales; runcases=[], spreadcases=[],
                   dtfilename = "dec_tree_all.csv", showr0 = true, silent=true)
#=
    see cases.jl for runcases and spreadcases
=#
    !isempty(spreadq) && (deleteat!(spreadq, 1:length(spreadq)))   # empty it

    # access input data and pre-allocate storage
    alldict = setup(geofilename, n_days; dectreefilename=dtfilename, geolim=15)
        dt = alldict["dt"]  # decision trees for transition
        openmx = alldict["dat"]["openmx"]
        cumhistmx = alldict["dat"]["cumhistmx"]
        newhistmx = alldict["dat"]["newhistmx"]
        isolatedmx = alldict["dat"]["isolatedmx"]
        geodata = alldict["geo"]
    dseries = build_series(locales)   # should use numgeo here?

    # simulation environment: pre-allocate arrays and initialize parameters for spread!
    env = initialize_sim_env()

    # start the counter at zero
    reset!(ctr, :day)  # remove key :day leftover from prior runs

    # initial data for building data series of simulation outcomes
    starting_unexposed = reduce(hcat, [grab(unexposed, agegrps, 1, i, dat=openmx) for i in locales])
    starting_unexposed = (size(locales,1) == 1 ? Dict(locales[1]=>starting_unexposed) : 
        Dict(locales[i]=>starting_unexposed[i,:] for i in size(locales,1)))

    locales = locales   # force local scope to be visible in the loop
    for i = 1:n_days
        inc!(ctr, :day)  # update the simulation day counter
        @debug "\n\n Start day $(ctr[:day])"
        silent || println("simulation day: ", ctr[:day])
        for locale in locales
            for case in runcases
                case(locale, opendat=openmx, isodat=isolatedmx, env=env)
            end
            density_factor = geodata[locale,density_fac]
            spread!(locale, density_factor, dat=openmx, env=env, spreadcases=spreadcases)
            transition!(dt, locale, dat=openmx)   # transition all infectious cases "in the open"
            transition!(dt, locale, dat=isolatedmx)  # transition all infectious cases in isolation
            if showr0 && (mod(ctr[:day],10) == 0)
                current_r0 = sim_r0(env=env)
                println("at day $(ctr[:day]) r0 = $current_r0")
            end
        end
        do_history!(locales, opendat=openmx, cumhist=cumhistmx, newhist=newhistmx, starting_unexposed=starting_unexposed)
    end
    println("Simulation completed for $(ctr[:day]) days.")

    series = Dict(loc=>Dict(:cum=>make_series(cumhistmx[loc]), :new=>make_series(newhistmx[loc])) for loc in locales)

    return alldict, env, series
end


function do_history!(locales; opendat, cumhist, newhist, starting_unexposed)
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
            cumhist[locale][:,1:5, thisday] = sum(opendat[locale],dims=1)[1,:,:]
            cumhist[locale][:,6, thisday] = sum(cumhist[locale][:, 1:5, thisday], dims=2)
            newhist[locale][:,:,thisday] = cumhist[locale][:,:,thisday] .- cumhist[locale][:,:,thisday-1]
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


# a single locale, either cumulative or new matrix
function make_series(histmx)
    s = zeros(Int, size(histmx,3), prod(size(histmx)[1:2]))
    for i in 1:size(histmx, 3)
        s[i, :] = reduce(vcat,[histmx[j, :, i] for j in 1:size(histmx,1)])'
    end
    return s
end



function sim_r0(;env=env)
    # captures current population condition 
    pct_unexposed = sum(env.simple_accessible[1,:]) / sum(env.simple_accessible)
    sa_pct = [pct_unexposed,(1-pct_unexposed)/2.0,(1-pct_unexposed)/2.0]   

    # if social_distancing case with split population
    if haskey(spread_stash, :case_cf) || haskey(spread_stash, :case_tf)
        compliance = env.sd_compliance
        cf = spread_stash[:case_cf]; tf = spread_stash[:case_tf]
        r0_comply = r0_sim(compliance = compliance, cf=cf, tf=tf, sa_pct=sa_pct).r0

        cf = spread_stash[:default_cf]; tf = spread_stash[:default_tf]
        r0_nocomply = r0_sim(compliance=(1.0 .- compliance), cf=cf, tf=tf, sa_pct=sa_pct).r0

        # this works if all compliance values are the same; approximate otherwise
        current_r0 = round(mean(compliance) * r0_comply + (1.0-mean(compliance)) * r0_nocomply, digits=2)
    else
        cf =  env.contact_factors
        tf = env.touch_factors     
        current_r0 = round(r0_sim(cf=cf, tf=tf, sa_pct=sa_pct).r0, digits=2)   
    end
    return current_r0
end


function initialize_sim_env()
    @assert laglim >= 19 "laglim must be >= 19--got $laglim"
    # expand send_risk_by_lag length to match laglim
    filln = laglim - 18
    mid = fill(0.3, filln)
    last3 = [0.3,0.3,0.1]
    first15 = [0.0, .3, .7, .8, .9, .9, .8, .7, .6, .5, .4, .3, .3, 0.3, 0.3]
    send_risk = [first15..., mid..., last3...]

    # initialize the simulation Env
    ret =   Env(spreaders=zeros(Int64,laglim,4,5),
                all_accessible=zeros(Int64,laglim,6,5),
                numcontacts=zeros(Int64,laglim,4,5),
                simple_accessible=zeros(Int64,6,5),
                numtouched=zeros(Int64,laglim,5),
                lag_contacts=zeros(Int64,laglim),
                riskmx = zeros(Float64,laglim,5),
                contact_factors =       [ 1    1.8    1.8     1.5    1;    # nil
                                          1    1.7    1.7     1.4   0.9;    # mild
                                        0.7    1.0    1.0     0.7   0.5;   # sick
                                        0.5    0.8    0.8     0.5   0.2],  # severe
                              # agegrp    1     2      3       4     5
                touch_factors =         [.55    .62     .58     .4   .35;    # unexposed
                                         .55    .62     .58     .4   .35;    # recovered
                                         .55    .62     .58     .4   .35;    # nil
                                         .55    .6      .5     .35   .28;    # mild
                                         .28   .35     .28    .18  .18;    # sick
                                         .18   .18     .18    .18  .18],   # severe
                send_risk_by_lag = send_risk,
                recv_risk_by_age = [.1, .4, .4, .50, .55],
                sd_compliance = zeros(6,5))
    ret.riskmx = send_risk_by_recv_risk(ret.send_risk_by_lag, ret.recv_risk_by_age)
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


function grab(condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    return dat[locale][lag, condition, agegrp]
end


function input!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    dat[locale][lag, condition, agegrp] = val
end


function plus!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    dat[locale][lag, condition, agegrp] += val
end


function minus!(val, condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    dat[locale][lag, condition, agegrp] -= val
end


function total!(condition, agegrp, lag, locale; dat=openmx)
    @assert length(locale) == 1 "locale must be a scalar"
    sum(dat[locale][lag, condition, agegrp])
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
