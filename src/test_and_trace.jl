"""
Quar_Loc is a type alias for 

```NamedTuple{(:locale, :start_date),Tuple{Int64,Int64}}```

example:  (locale=53038, start_date=50)

Used as a locale in an isolatedmx 3-d array. Quar_Loc locale provides
the geographic locale and the date that a cohort of infected folks
entered quarantine. Locale values are US Census FIPS codes at the county level.
"""
const Quar_Loc = NamedTuple{(:locale, :start_date),Tuple{Int64,Int64}}


const map2access = (unexposed= 1, infectious=-1, recovered= 2, dead=-1, 
                  nil= 3, mild=  4, sick= -1, severe= -1)


"""
Generate test_and_trace case.
required inputs: policyday
keyword arguments: tc_perday, sensitivity, specificity, infect_prior, test_pct, test_delay
required keywork arguments: env, opendat, isodat, locale

Returns a function to use in runcases input to run_a_sim.
"""
function t_n_t_case_gen(start_day, end_day;         # these args go into the returned t_n_t_case
    tc_perday=400, sensitivity=.90, specificity=0.90, infect_prior=0.5, test_pct=.95,  # optional
    test_delay=1, generations=3, qdays=15) 
    # args match runcases loop in run_a_sim
    function scase(locale; opendat, isodat, env)  # case loop in run_a_sim provides required args
        t_n_t_case(start_day, end_day; 
                   env=env, opendat=opendat, isodat=isodat, locale=locale,  # from case loop
                   tc_perday=tc_perday, sensitivity=sensitivity, specificity=specificity, 
                   infect_prior=infect_prior, test_pct=test_pct,  
                   test_delay=test_delay, generations=generations, qdays=qdays)
    end
end


function t_n_t_case(start_date, end_date; 
                env, opendat, isodat, locale,     # from case loop
                tc_perday=1000, sensitivity=.95, specificity=0.90, infect_prior=0.5, test_pct=.95,  
                test_delay=1, generations=3, qdays=15)

    # do we need to UNquarantine anyone today?  this finds all dated quarantines until they are gone
    for k in keys(isodat)
        if typeof(k) <: Quar_Loc # then k is a Quar_Loc = (locale=locale, start_date=thisday)
            if k.locale == locale
                if k.start_date + qdays == ctr[:day] # end of quarantine is today
                   allout = t_n_t_unquarantine(k, opendat=opendat, isodat=isodat, env=env)
                   delete!(isodat, k)  
                   push!(tntq, (day=ctr[:day], unq=sum(allout)))
                end
            end
        end
    end

    test_and_trace(start_date, end_date; 
        env=env, opendat=opendat, isodat=isodat, locale=locale,   # from case loop
        tc_perday=tc_perday, sensitivity=sensitivity, specificity=specificity, # optional
        infect_prior=infect_prior, test_pct=test_pct, test_delay=test_delay, 
        generations=generations, qdays=qdays)
end


function test_and_trace(start_date, end_date; 
    env, opendat, isodat, locale,  # from case loop
    tc_perday=1000, sensitivity=.95, specificity=0.95, infect_prior=0.5, test_pct=.95, # optional
    test_delay=1, generations=3, qdays=15)
    
    thisday = ctr[:day]
    test_conds = [unexposed, recovered, nil, mild]

    if start_date <= thisday < end_date

        # who gets tested?  # (25,4,5)
        avail_to_test = grab(test_conds, agegrps, lags, locale, dat=opendat)

        if sum(avail_to_test) == 0
            println("on $thisday no one to test")
            return
        end

        to_test = Dict(i=>zeros(Int, laglim, 4, 5) for i in 1:generations) # TODO put these in a struct
        postests = Dict(i=>zeros(Int, laglim, 4, 5) for i in 1:generations)
        poscontacts = Dict(i=>zeros(Int, laglim, 4, 5) for i in 1:generations)
        postouched = Dict(i=>zeros(Int, laglim, 4, 5) for i in 1:generations)

        # create new tracking locale if doesn't exist
        qloc = (locale=locale, start_date=thisday)
        if !haskey(isodat, qloc)
            isodat[qloc] = zeros(Int, laglim, length(conditions), length(agegrps))
        end
        
        density_factor = env.geodata[env.geodata[:, fips] .== locale, density_fac][1]
        conducted = 0 # per gen
        perday_conducted = 0 # per day

        for gen in 1:generations
            if gen == 1
                to_test[gen] = avail_to_test
            else
                to_test[gen] = postouched[gen-1]
                tc_perday -= conducted
            end

            # test
            postests[gen], conducted = simtests(to_test[gen]; tc_perday=tc_perday, sensitivity=sensitivity, 
                            specificity=specificity, infect_prior=0.05, test_pct=test_pct, env=env)
            perday_conducted += conducted

            # trace  
            # contacts
                target_cf = repeat(view(env.contact_factors,1,:)',4,1) # use unexposed for all rows
            poscontacts[gen] = how_many_contacts!(poscontacts[gen], 
                                    postests[gen],  # equivalent to spreaders in spread
                                    avail_to_test,
                                    target_cf, 
                                    density_factor, env=env)  
            # touched--consquential contact that we count
                target_tf = view(env.touch_factors,map2access[unexposed]:map2access[mild], agegrps)
            postouched[gen] = how_many_touched!(postouched[gen], poscontacts[gen], 
                                        avail_to_test, test_conds, 
                                        target_tf, env=env, kind=:trace)
            # isolate
            t_n_t_quarantine(postests[gen], qloc::Quar_Loc; opendat=opendat, isodat=isodat, env=env)

        end  # for gen 

        push!(tntq,(day=thisday, avail=reduce(+, map(sum,values(to_test))), 
                    conducted=perday_conducted, postests=reduce(+, map(sum,values(postests))), 
                    poscontacts=reduce(+, map(sum,values(poscontacts))), 
                    postouched=reduce(+, map(sum,values(postouched)))))

    end  # if start_day
end


function simtests(to_test; tc_perday=1000, sensitivity=.9, specificity=.9, infect_prior=.05, 
    test_pct=.95, env=env)

    # distribute the tests across disease conditions by age group and condition
           # we could do probabilistically but the probs are very small and the whole thing
           # is somewhat artificial: we can capture "randomness" by randomly ignoring x% of the results

    if tc_perday <= 0
        println("got here")
        return zeros(Int,laglim, 4, length(agegrps)), 0
    end

    today_tests = rand(Binomial(tc_perday, test_pct), 1)[1]

    # println("today $(ctr[:day]) tc_perday $tc_perday  today tests $today_tests")

    tst_pct = to_test ./ sum(to_test)
    tst_pct[isnan.(tst_pct)] .= 0.0  # eliminate NaNs

    dist_tests = round.(Int, (tst_pct .- 1e-5) .* today_tests) # tst_pct .- 1e-5

    # test results  

        # specificity tests apply to actual true not-infected: unexposed recovered
        @views false_pos = rand.(Binomial.(dist_tests[:, 1:2, :], 1.0 - specificity))

        # sensitivity tests apply to actual true infected: nil, mild
        @views true_pos = rand.(Binomial.(dist_tests[:, 3:4, :], sensitivity))
        false_neg = dist_tests[:, 3:4,:] - true_pos   # TODO should report this for curiosity

        pos_results = cat(false_pos,true_pos,dims=2)  # (25,4,5)

        # println(" day $(ctr[:day])  ")
        # println(" False Pos  quarantined even though not sick: ", sum(false_pos) ,", ", 
        #         round(sum(false_pos) / sum(dist_tests[:, 1:2, :]), digits=4))
        # println(" True Pos  quarantined because actually sick: ", sum(true_pos) ,", ", 
        #         round(sum(true_pos) / sum(dist_tests[:, 3:4, :]), digits=4))
        # println(" False Neg  not quarantined even though sick: ", sum(false_neg) ,", ",
        #         round(sum(false_neg) / sum(dist_tests[:, 3:4,:]), digits=4))
        # println(" Total positive results, total tests conducted ", sum(pos_results), ", ", 
            # sum(dist_tests))

    return pos_results, sum(dist_tests)
end


function bayes(sensitivity, specificity, pr_pop) 
    pr_pos_given_pos_test = (
                       sensitivity * pr_pop /
       (sensitivity * pr_pop + (1.0 - specificity) * (1.0 - pr_pop))
    ) # Bayesian probs. better but public tests don't use Bayes interpretation!
end


function t_n_t_quarantine(postests, qloc::Quar_Loc; opendat, isodat, env)
    test_conds = [unexposed, recovered, nil, mild] 
    isolate_by!(postests, test_conds, agegrps, lags, qloc, opendat=opendat, isodat=isodat)
end


function t_n_t_unquarantine(qloc::Quar_Loc; opendat, isodat, env)
    ret_conds = [unexposed, recovered, dead, nil, mild, sick, severe] 
    allout = grab(ret_conds, agegrps, lags, qloc, dat=isodat)

    # println("day $(ctr[:day]) unquaranting this many ", sum(allout))

    unisolate_by!(allout, ret_conds, agegrps, lags, qloc; 
                  opendat=opendat, isodat=isodat, mode=:plus) # do plus! becasuse we delete the qloc
    return allout
end