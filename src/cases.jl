####################################################################
# cases.jl
#       pass these cases to run_a_sim in kwarg runcases as a list
#       cases = [CovidSim.<funcname>, CovidSim.<funcname>, ...]  
#       then run_a_sim(geofilename, n_days, locales; runcases=cases)
####################################################################


"""
Generate seeding cases.
inputs: day, cnt, lag, cond, agegrp
Two of the inputs may refer to multiple items and must match in number of items.

Returns a function that can be used in runcases input to run_a_sim.
"""
function seed_case_gen(day, cnt, lag, cond, agegrp) # these args go into the returned seed! case
    function scase(locale; opendat=openmx, isodat=isolatedmx, env=env)  # args must match runcases loop in run_a_sim
        seed!(day, cnt, lag, cond, agegrp, locale, dat=opendat)
    end
end

# some generated seed! cases-->these are global (in code)
# seed_6_12 = seed_case_gen(8, [0,6,6,0,0], 5, nil, agegrps)
# seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 5, nil, agegrps)


# some isolation cases
function isolate_case_1(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day] == 15
        isolate!(.25,[unexposed, nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.70,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    elseif ctr[:day] == 23
        isolate!(.50,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.70,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_1(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 120
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(1.0,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function isolate_case_2(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day] == 15
        isolate!(.40,[unexposed, nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    elseif ctr[:day] == 23
        isolate!(.60,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_2(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 69
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(1.0,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_2b(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 84
        unisolate!(.6,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(.6,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end


function isolate_case_3(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day] == 40
        isolate!(.40,[unexposed, nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    elseif ctr[:day] == 50
        isolate!(.60,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_3(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 80
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(1.0,[mild,sick,severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

"""
- define a cases as mycase=Spreadcase(15,cf_array,tf_array,compliance_array_or_float)
- pass these cases to run_a_sim in kwarg spreadcases as a list--they'll be run in function spread!
scases = [mycase, case_2, ...]  then run_a_sim(geofilename, n_days, locales; spreadcases=scases)
- cases above can be combined with these passing in both runcases and spreadcases
"""
struct Spreadcase
    day::Int
    cf::Array{Float64,2}  # (4,5)
    tf::Array{Float64,2}  # (6,5)
    compliance::Union{Float64, Array{Float64,2}}
end


function spread_case_setter(cases=[]; env=env)
    for case in cases
        c = case(env=env)
        if c.day == ctr[:day]
            # before the case starts--ignore it
            # after the case--it's already in effect--nothing to change
            if iszero(c.compliance)  # signal to shutdown cases and restore defaults
                # restore defaults for spread!  
                default_env = initialize_sim_env()
                env.sd_compliance = default_env.sd_compliance
                env.contact_factors = default_env.contact_factors
                env.touch_factors = default_env.touch_factors
                delete!(spread_stash, :case_cf)
                delete!(spread_stash, :case_tf)
            else
            # set contact_factors, touch_factors to the case values
                # copy c.cf to spread_stash
                # copy c.cf to env.contact_factors->active contact_factors  => when running spreadsteps
                if !haskey(spread_stash, :default_cf)  # should only ever happen once for an entire simulation
                    spread_stash[:default_cf] = copy(env.contact_factors)
                end
                if !haskey(spread_stash, :default_tf)
                    spread_stash[:default_tf] = copy(env.touch_factors)
                end

                spread_stash[:case_cf] = copy(c.cf)  # shouldn't need copy, but it's safer
                spread_stash[:case_tf] = copy(c.tf)  #             "

            # set the compliance note: compliance is the same for spreaders and accessible
                    # it varies by agegrp and condition if desired
                # check all compliance values in [0.0, 1.0]
                @assert c.compliance .>= 0.0 "compliance values must be positive"
                @assert c.compliance .<= 1.0 "compliance values must be in [0.0,1.0]"
                env.sd_compliance .= copy(c.compliance)  # TODO do we need to copy? takes 2x time
            end # if for current day case
        end  # if test for today
    end # case for loop
end  # function case_setter


function spread_case_runner(density_factor, all_unexposed; env=env)
    spread_stash[:spreaders] = copy(env.spreaders)  # stash today's spreaders--isolated from env
    spread_stash[:simple_accessible] = copy(env.simple_accessible) # stash today's accessible--isolated from env
    newinfected = []  # capture infected for comply and nocomply groups
    for i in [:comply,:nocomply]
        if i == :comply  # split the spreaders and accessible, set the case factors
            env.spreaders[:]= round.(Int,permutedims(permutedims(copy(spread_stash[:spreaders]),[2,3,1]) .*
                                       env.sd_compliance[3:6,:], [3,1,2]))
            env.simple_accessible[:]= round.(Int,copy(spread_stash[:simple_accessible]) .*
                                             env.sd_compliance)
            env.contact_factors = copy(spread_stash[:case_cf])
            env.touch_factors = copy(spread_stash[:case_tf])
        else  # i == :nocomply other split of spreaders and accessible, restore default factors
            env.spreaders[:]= round.(Int, permutedims(permutedims(copy(spread_stash[:spreaders]),[2,3,1]) .*
                                        (1.0 .- env.sd_compliance[3:6,:]), [3,1,2]))
            env.simple_accessible[:]= round.(Int, copy(spread_stash[:simple_accessible]) .*
                                             (1.0 .- env.sd_compliance))
            # set the default contact_factors and touch_factors
            env.contact_factors = copy(spread_stash[:default_cf])
            env.touch_factors = copy(spread_stash[:default_tf])
        end  # if
        push!(newinfected, spreadsteps(density_factor, all_unexposed, env=env))
        if i == :comply
            spread_stash[:comply_contacts] = copy(env.numcontacts)
            spread_stash[:comply_touched] = copy(env.numtouched)
            spread_stash[:comply_spreaders] = copy(env.spreaders)
        end
    end  # for loop
    # total values for comply + nocomply
    newinfected = newinfected[1] .+ newinfected[2]
end


function sd_gen(;start=45, comply=.7, cf=(.2, 1.6), tf=(.18,.7))
    function sd_mod(;env=env)
        sd_mod = Spreadcase(start,
                        round.(shifter(env.contact_factors, cf...),digits=2),
                        round.(shifter(env.touch_factors, tf...),digits=2),
                        comply)
    end
end

# copy beyond the comment and run in the REPL, use as input
#
# mod_45 = sd_gen()  # with defaults
# mod_90 = sd_gen(start=90,cf=(.2,1.5), tf=(.18,.6),comply=.85)
# str_45 = sd_gen(start=45, comply=.90, cf=(.2,1.0), tf=(.18,.3))
# str_55 = sd_gen(start=55, comply=.95, cf=(.2,1.0), tf=(.18,.3))
# zer = sd_gen(start=90, comply=0.0)


function test_and_trace(;tc_day=400, sensitivity=.98, specificity=0.95, infect_prior=0.5,
    test_pct=.70, test_delay=1, env=env, dat=openmx, generations = 3,
    locale)   # must supply

    map2access = (unexposed= 1, infectious=-1, recovered= 2, dead=-1, nil= 3, mild=  4, sick= -1, severe= -1)
    thisday = ctr[:day]
    to_test = Dict(i=>zeros(Int, 4, 5) for i in 1:generations)
    postests = Dict(i=>zeros(Int, 4, 5) for i in 1:generations)
    postouched = Dict(i=>zeros(Int, 4, 5) for i in 1:generations)

    # who gets tested?  #  before sum (24,4,5); after (4,5)
    avail_to_test = sum(grab([unexposed, recovered, nil, mild], agegrps, lags, locale, dat=dat), dims=1)[1,:,:] 
    # get density_factor for locale

    for gen in 1:3
        if gen == 1
            to_test[gen] = avail_to_test
        else
            to_test[gen] = postouched[gen-1]
        end

        postests[gen] = simtests(to_test[gen]; tc_day=tc_day, sensitivity=sensitivity, 
                             specificity=specificity, infect_prior=0.5)

        # trace 
        poscontacts[gen] = pos_contacts(postests[gen], env=env)  # TODO need to pass in locale's density_factor
        postouched[gen] = pos_touched(poscontacts[gen], env=env)

        # isolate the initial positives and their positive contacts
            # repeat for 2 more generations or stop if touches exhausted
        lag = 5
        for a in agegrps
            for cond in [unexposed, recovered, nil, mild]
                isolate_by!(postest[gen][map2access[cond],a], cond, a, lag, 
                    (locale=locale, startday=thisday), opendat=dat, isodat=isolatedmx )
            end
        end
        # what to do about lag?  we assume that the person and the tester don't know how long
            # the person has had the disease. Simulation knows, but we dropped the information because 
            # the testers can't know it and the required quarantine will be the same.  
            # When we isolate, we still need to transition people. TODO--go back and keep the lag information.  
            # For now, assume lag = 5 because 14 days gets to a decision point.
    end
    
end


function simtests(to_test; tc_day, sensitivity, specificity, infect_prior, env=env)

    # randomly distribute the tests across disease conditions by age group and condition
    avail_by_age = sum(to_test,dims=1)    #(1,5)
    avail_by_age_pct = vec(avail_by_age ./ sum(avail_by_age))
    avail_by_cond_pct = to_test ./ avail_by_age
    x = categorical_sample(avail_by_age_pct, tc_day)
    dist_tests = [count(x .== i) for i in agegrps]
    apply_tests = round.(Int, avail_by_cond_pct .* reshape(dist_tests, 1, agegrps))

    # test results  # Bayesian probs. better but public tests don't use Bayes interpretation
        # specificity tests apply to actual true not-infected: unexposed recovered
        false_pos = rand.(Binomial.(apply_tests[1:2,:], 1.0 - specificity))
        # sensitivity tests apply to actual true infected: nil mild
        true_pos = rand.(Binomial.(apply_tests[3:4,:], sensitivity))

    test_results = [false_pos; true_pos]  # (4,5)
end


function bayes(sensitivity, specificity, pr_pop)
    pr_pos_given_pos_test = (
                       sensitivity * pr_pop /
       (sensitivity * pr_pop + (1.0 - specificity) * (1.0 - pr_pop))
    )
end


function pos_contacts(pos_tests, density_factor=1.0; env=env)   # should we use the env matrices?   maybe not--fix it later

    # variables from env
    # spreaders = env.spreaders  # laglim,4,5
    # numcontacts = env.numcontacts  # the result

    all_accessible = env.all_accessible
    contactcnt = copy(env.contact_factors[1,:])  # we want only the row for the nil condition
    contact_factors = repeat(contactcnt',4,1)
    numcontacts = zeros(Int, 4,5)

    if sum(env.simple_accessible) == 0   # this can happen with a social distancing case with 100% compliance
        numcontacts[:] .= 0
        return
    end

    pos_conds, pos_ages = size(pos_tests)


    for agegrp in 1:pos_ages
        for cond in 1:pos_conds
            scale = contact_factors[cond, agegrp]
            poscount = pos_tests[cond, agegrp]
            if poscount == 0
                numcontacts[cond, agegrp] = 0
                continue
            end
            dgamma = Gamma(1.0, density_factor * scale)  #shape, scale
            x = rand(dgamma,poscount)
            numcontacts[cond, agegrp] = round(Int,sum(x))
        end
    end

    # correct over contacting: this can happen with high scale towards the high point of infection
    oc_ratio = sum(numcontacts) / sum(all_accessible)
    if oc_ratio > 1.0
        println(ctr[:day],"warning: overcontact ratio ", oc_ratio)
        numcontacts[:] = round.(1.0/oc_ratio .* numcontacts)
    end

    return numcontacts
end


function pos_touched(poscontacts; env=env)
#=
    - who has been touched by those who tested positive? 
    - poscontacts (4, 5)  cond x agegrp: unexposed, recovered, nil, mild
=#
    # map to access maps conditions to the rows of simple_accessible and touch_factors
    map2access = (unexposed= 1, infectious=-1, recovered= 2, dead=-1, nil= 3, mild=  4, sick= 5, severe= 6)

    simple_accessible = env.simple_accessible
    apply_tf =     view(env.touch_factors,map2access[unexposed]:map2access[mild], agegrps)
    numtouched =        zeros(Int, 4,5)         # result 

    # println("apply_tf = ")
    # display(apply_tf)

    if sum(simple_accessible) == 0   
        numtouched[:] .= 0
        return
    end

    track_accessible = view(simple_accessible, map2access.unexposed:map2access.mild, :)  # (4, 5)

    println("track_accessible = ")
    display(track_accessible)
    # s_a_pct is dist. of accessible by agegrp and uexposed, recovered, nil, mild  (4,5)
    t_a_pct = round.(reshape(track_accessible ./ sum(track_accessible), 20), digits=5) # % for each cell
    if !isapprox(sum(t_a_pct), 1.0)
        t_a_pct = t_a_pct ./ sum(t_a_pct) # normalize so sums to 1.0
    end

    # println("size(t_a_pct) = ", size(t_a_pct))
    # println("t_a_pct = ", t_a_pct)


    dcat = Categorical(t_a_pct) # categorical distribution by agegrp and unexposed, recovered, nil, mild

    for cond in [unexposed, recovered, nil, mild]
        for a in agegrps 
            cond_age = poscontacts[map2access[cond],a]
            if cond_age == 0
                numtouched[cond, a] = 0
            else  # probabilistically see who of the accessible is significantly touched
                x = rand(dcat, cond_age) # probabistically distribute contacts for an age across accessible by unexposed|recovered|nil|mild, agegrp
                peeps = reshape([count(x .== i) for i in 1:20], 4,5)# (5,) distribute across all 3 groups, but only take unexposed
                # println("peeps = ")
                # display(peeps)
                cnt = binomial_one_sample.(peeps, apply_tf) 
                # println("cond $cond, age $a cnt = ")
                # display(cnt)
                # numtouched[map2access[cond], a] += sum(clamp.(cnt, 0, track_accessible))
                numtouched .+= clamp.(cnt, 0, track_accessible)
                # println("numtouched = ")
                # display(numtouched)
            end
        end
    end
    return numtouched
end