"""
Generate test_and_trace case.
required inputs: policyday
keyword arguments: tc_perday, sensitivity, specificity, infect_prior, test_pct, test_delay
required keywork arguments: env, opendat, isodat, locale

Returns a function that can be used in runcases input to run_a_sim.
"""
function t_n_t_case_gen(start_day, end_day;         # these args go into the returned t_n_t_case
    tc_perday=400, sensitivity=.90, specificity=0.90, infect_prior=0.5, test_pct=.70,  # optional
    test_delay=1, generations=3, qdays=15) 
    # args match runcases loop in run_a_sim
    function scase(locale; opendat, isodat, env)  # case loop in run_a_sim provides required args
        t_n_t_case(start_day, end_day; env=env, opendat=opendat, isodat=isodat, locale=locale,
                   tc_perday=tc_perday, sensitivity=sensitivity, specificity=specificity, 
                   infect_prior=infect_prior, test_pct=test_pct,  
                   test_delay=test_delay, generations=generations, qdays=qdays)
    end
end


function t_n_t_case(start_date, end_date; env, opendat, isodat, locale,     # required keyword args
    tc_perday=400, sensitivity=.95, specificity=0.90, infect_prior=0.5, test_pct=.70,  # optional
    test_delay=1, generations=3, qdays=15)

    # do we need to unquarantine anyone today?  this finds all dated isodats until they are gone
    for k in keys(isodat)
        if typeof(k) <: NamedTumple
            if k.locale == locale
                if k.start_date + qdays == ctr[:day] # end of quarantine is today
                   t_n_t_unquarantine(locale, k.start_date, opendat=opendat, isodat=isodat, env=env)
                   delete!(isodat, (locale=locale, start_date=k.start_date))
                end
            end
        end
    end

    # time to test again   # first line of args are all required
    test_and_trace(start_date, end_date; env=env, opendat=opendat, isodat=isodat, locale=locale,   
        tc_perday=400, sensitivity=.98, specificity=0.95, infect_prior=0.5, test_pct=.70, # optional
        test_delay=1, generations=3, qdays=15)

end


function test_and_trace(start_date, end_date; env, opendat, isodat, locale,   # required keyword args
    tc_perday=400, sensitivity=.98, specificity=0.95, infect_prior=0.5, test_pct=.70, # optional
    test_delay=1, generations = 3, qdays = 15)
    
    thisday = ctr[:day]
    map2access = (unexposed= 1, infectious=-1, recovered= 2, dead=-1, 
                  nil= 3, mild=  4, sick= -1, severe= -1)

    if start_date <= thisday < end_date

        to_test = Dict(i=>zeros(Int, 4, 5) for i in 1:generations)
        postests = Dict(i=>zeros(Int, 4, 5) for i in 1:generations)
        poscontacts = Dict(i=>zeros(Int, 4, 5) for i in 1:generations)
        postouched = Dict(i=>zeros(Int, 4, 5) for i in 1:generations)

        # who gets tested?  #  before sum (24,4,5); after (4,5)
        avail_to_test = sum(grab([unexposed, recovered, nil, mild], agegrps, 
                                 lags, locale, dat=dat), dims=1)[1,:,:] 
        
        density_factor = geodata[geodata[fips, :] .== locale, density_fac]

        for gen in 1:3
            if gen == 1
                to_test[gen] = avail_to_test
            else
                to_test[gen] = postouched[gen-1]
            end

            postests[gen] = simtests(to_test[gen]; tc_perday=tc_perday, sensitivity=sensitivity, 
                                 specificity=specificity, infect_prior=0.5)

            # trace   # note: reshape is a view so no need to reshape back to original shape
            poscontacts[gen] = how_many_contacts!(reshape(poscontacts[gen],1,4,5), reshape(postests[gen],1,4,5), 
                               repeat(view(env.contact_factors,1,:)',4,1), density_factor, env=env)  

            # target_accessible  (4, 5)
                target_tf = view(env.touch_factors,map2access[unexposed]:map2access[mild], agegrps)
                target_accessible = view(env.simple_accessible, map2access.unexposed:map2access.mild, :)  # (4, 5)
            postouched[gen] = how_many_touched!(postouched[gen], reshape(poscontacts[gen],1,4,5), 
                                        target_accessible, [unexposed, recovered, nil, mild], 
                                        target_tf, env=env, kind=:trace)

            # isolate the initial positives and their positive contacts
                # repeat for 2 more generations or stop if touches exhausted
            lag = 5
            for a in agegrps
                for cond in [unexposed, recovered, nil, mild]
                    isolate_by!(postest[gen][map2access[cond],a], cond, a, lag, 
                        (locale=locale, start_date=thisday), opendat=dat, isodat=isolatedmx )
                end
            end
            # what to do about lag?  we assume that the person and the tester don't know how long
                # the person has had the disease. Simulation knows, but we dropped the 
                # information because the testers can't know it and the required quarantine 
                # will be the same.  TODO--go back and keep the lag information. 
                # When we isolate, we still need to transition people.  
                # For now, assume lag = 5 because 14 days gets to a decision point.
        end  # for gen 
    end  # if start_day
end


function simtests(to_test; tc_perday, sensitivity, specificity, infect_prior, env=env)

    # randomly distribute the tests across disease conditions by age group and condition
    avail_by_age = sum(to_test,dims=1)    #(1,5)
    avail_by_age_pct = vec(avail_by_age ./ sum(avail_by_age))
    avail_by_cond_pct = to_test ./ avail_by_age
    x = categorical_sample(avail_by_age_pct, tc_perday)
    dist_tests = [count(x .== i) for i in agegrps]
    apply_tests = round.(Int, avail_by_cond_pct .* reshape(dist_tests, 1, agegrps))

    # test results  # Bayesian probs. better but public tests don't use Bayes interpretation
        # specificity tests apply to actual true not-infected: unexposed recovered
        false_pos = rand.(Binomial.(apply_tests[1:2,:], 1.0 - specificity))
        # sensitivity tests apply to actual true infected: nil mild
        true_pos = rand.(Binomial.(apply_tests[3:4,:], sensitivity))
        false_neg = apply_tests[3:4,:] - true_pos   # TODO should report this for curiosity

    test_results = [false_pos; true_pos]  # (4,5)
end


function bayes(sensitivity, specificity, pr_pop)
    pr_pos_given_pos_test = (
                       sensitivity * pr_pop /
       (sensitivity * pr_pop + (1.0 - specificity) * (1.0 - pr_pop))
    )
end



function t_n_t_unquarantine(loc, start_date, opendat, isodat, env)
    getconds = [unexposed, recovered, nil, mild, sick, severe] # last 2 will turn up zero
    getfrom = (locale=loc, start_date=start_date)
    allin = grab(getconds, agegrps, lags, getfrom, dat=isodat)

    # only need to do plus! becasuse we delete this isolatedmx value
    for a in agegrps
        for cond in [unexposed, recovered, nil, mild]
            unisolate_by!(num, cond, a, lag, (locale=loc, start_date=start_date);
                         opendat=opendat, isodat=isodat)
        end
    end
end