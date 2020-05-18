"""
Quar_Loc is a type alias for 

```NamedTuple{(:locale, :start_date),Tuple{Int64,Int64}}```

example:  (locale=53038, start_date=50)

Used as a locale in an isolatedmx 3-d array. Quar_Loc locale provides
the geographic locale and the simulation ordinal date that a cohort of infected folks
entered quarantine. Locale values are US Census FIPS codes at the county level.
"""
const Quar_Loc = NamedTuple{(:locale, :start_date),Tuple{Int64,Int64}}


const map2access = (unexposed= 1, infectious=-1, recovered=2, dead=-1, 
                  nil= 3, mild=  4, sick= -1, severe= -1)


# cache values needed during test and track
const tnt_stash = Dict{Quar_Loc, Array}()


"""
    function t_n_t_case_gen(start_day, end_day;         
        tc_perday=400, sensitivity=.90, specificity=0.90, infect_prior=0.5, 
        test_pct=.95,  test_delay=1, generations=3, qdays=15) 

Returns a function to use in runcases input to function ```run_a_sim```.
"""
function t_n_t_case_gen(start_day, end_day;         # these args go into the returned t_n_t_case
    tc_perday=400, sensitivity=.90, specificity=0.90, infect_prior=0.5, test_pct=.95,  # optional
    q_comply=0.8, c_comply=0.9, breakout_pct=.3, test_delay=1, generations=3, qdays=15) 
    # args match runcases loop in run_a_sim
    function scase(locale; opendat, isodat, env)  # case loop in run_a_sim provides required args
        t_n_t_case(start_day, end_day; 
                   env=env, opendat=opendat, isodat=isodat, locale=locale,  # from case loop
                   tc_perday=tc_perday, sensitivity=sensitivity, specificity=specificity, 
                   infect_prior=infect_prior, test_pct=test_pct, q_comply=q_comply, c_comply=c_comply,
				   breakout_pct=breakout_pct, test_delay=test_delay, generations=generations, qdays=qdays)
    end
end


function t_n_t_case(start_date, end_date; 
                env, opendat, isodat, locale,     # from case loop
                tc_perday=1000, sensitivity=.95, specificity=0.90, infect_prior=0.5, test_pct=.95,  
                q_comply=0.8, c_comply=0.9, breakout_pct=.3, test_delay=1, generations=3, qdays=15)

    thisday = ctr[:day]
    ret_conds = [unexposed, recovered, nil, mild, sick, severe] 

    # do we need to UNquarantine anyone today?  
        # finds each dated quarantines until they are gone
        # check for breakouts from quarantine
    for k in keys(isodat)
        if typeof(k) <: Quar_Loc # then k is a Quar_Loc = (locale=locale, start_date=thisday)
            if k.locale == locale
                if k.start_date < thisday < k.start_date + qdays
                    day_of_q = thisday - k.start_date
                    cnt = get(tnt_stash[k], day_of_q, 0) # are there breakouts?
                    if cnt > 0
                        # println("  GOT HERE:  unquarantine breakouts  $cnt ")
                        t_n_t_unquarantine(cnt_2_array(cnt, isodat[k]), k, opendat=opendat, 
                                           isodat=isodat, env=env)
                        push!(tntq, (day=thisday, breakout=cnt))
                    end
                    # end
                elseif k.start_date + qdays == thisday # end of quarantine is today
                    cnt = grab(ret_conds, agegrps, lags, k, dat=isodat)
                    t_n_t_unquarantine(cnt, k, opendat=opendat, isodat=isodat, env=env)
                    push!(tntq, (day=ctr[:day], unquarantine=sum(cnt)))
                    delete!(isodat, k)  # remove dated locale
                    delete!(tnt_stash, k)  # remove stash for dated locale
                end
            end
        end
    end

    test_and_trace(start_date, end_date; 
        env=env, opendat=opendat, isodat=isodat, locale=locale,   # from case loop
        tc_perday=tc_perday, sensitivity=sensitivity, specificity=specificity, # optional
        infect_prior=infect_prior, test_pct=test_pct, q_comply=q_comply, c_comply=c_comply, 
		breakout_pct=breakout_pct, test_delay=test_delay, generations=generations, qdays=qdays)
end


function test_and_trace(start_date, end_date; 
    env, opendat, isodat, locale,  # from case loop
    tc_perday=1000, sensitivity=.95, specificity=0.95, infect_prior=0.5, test_pct=.95, # optional
    q_comply=0.8, c_comply=0.9, breakout_pct=.3, test_delay=1, generations=3, qdays=15)
    
    thisday = ctr[:day]
    test_conds = [unexposed, recovered, nil, mild]

    if thisday  == 1    # TODO when is the right time?  what is the right cleanup?
        cleanup_stash(tnt_stash)
    end


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

        # create new tracking locale and stash if doesn't exist
        qloc = (locale=locale, start_date=thisday)
        if !haskey(isodat, qloc)
            isodat[qloc] = zeros(Int, laglim, length(conditions), length(agegrps))
        end
        if !haskey(tnt_stash, qloc)
            tnt_stash[qloc] = zeros(Int, qdays-1)  # there is no breakout on the last day--everyone's out
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

            # trace contacts
                target_cf = repeat(view(env.contact_factors,1,:)',4,1) # use unexposed for all rows
            poscontacts[gen] = how_many_contacts!(poscontacts[gen], 
                                    postests[gen],  # equivalent to spreaders in spread
                                    avail_to_test,
                                    target_cf, 
                                    density_factor, env=env)  								
			poscontacts[gen] = round.(Int, c_comply .* poscontacts[gen])

            # contacts lead to consquential touches that we count
                target_tf = view(env.touch_factors,map2access[unexposed]:map2access[mild], agegrps)
            postouched[gen] = how_many_touched!(postouched[gen], poscontacts[gen], 
                                        avail_to_test, test_conds, 
                                        target_tf, env=env, kind=:trace)

            # isolate positives and cache breakout
            put_in = round.(Int, q_comply .* postests[gen])
            t_n_t_quarantine(put_in, qloc::Quar_Loc; opendat=opendat, isodat=isodat, env=env)
            if breakout_pct != 0.0  # future breakouts from this quarantine cohort
                tnt_stash[qloc] .+= Int.(breakout(breakout_pct, postests[gen]))

                # println("  tnt_stash $qloc ", tnt_stash[qloc])


            end
        end  # for gen 

        # statistics
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

    if tc_perday <= 0  # earlier generations of test and trace used up the available tests today
        return zeros(Int,laglim, 4, length(agegrps)), 0
    end

    # today_tests = rand(Binomial(tc_perday, test_pct), 1)[1]

    # println("today $(ctr[:day]) tc_perday $tc_perday  today tests $today_tests")

    alloc_pct = reshape(to_test ./ sum(to_test), length(to_test))
    alloc_pct[isnan.(alloc_pct)] .= 0.0  # eliminate NaNs (underflow)
    @assert isapprox(sum(alloc_pct), 1.0) "probabilities must sum to 1.0"

    # dist_tests = round.(Int, (alloc_pct .- 1e-5) .* today_tests) 
    dist_tests = zeros(Int, size(to_test))
    dcat = Categorical(alloc_pct)
    x = rand(dcat, round(Int, test_pct * tc_perday))
    dist_tests[:] = reshape([count(x .== i) for i in 1:length(dcat.p)], size(dist_tests))
    dist_tests[:] = clamp.(dist_tests, 0, to_test)

    sum(dist_tests) > sum(to_test) && (@warn "Happy Day: more tests than people to test")

    # test results  

        # specificity tests apply to actual true not-infected: unexposed recovered
        false_pos = rand.(Binomial.(dist_tests[:, 1:2, :], 1.0 - specificity)) # @views 

        # sensitivity tests apply to actual true infected: nil, mild
        true_pos = rand.(Binomial.(dist_tests[:, 3:4, :], sensitivity))  # @views 
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


function t_n_t_unquarantine(cnt, qloc::Quar_Loc; opendat, isodat, env)
    ret_conds = [unexposed, recovered, nil, mild, sick, severe] 
    unisolate_by!(cnt, ret_conds, agegrps, lags, qloc; 
                  opendat=opendat, isodat=isodat, mode=:plus) # delete the qloc when unq all
end


function breakout(breakout_pct, inq; 
    breakout_dist = [0.0, 0.0, 0.0, 0.0, 0.0, 0.03, 0.03, 
                    0.03, 0.03, 0.03, 0.03, 0.05, 0.05, 0.05, 0.67])

    m = breakout_pct / sum(breakout_dist[1:14])
    breakout_dist[1:14] .*= m
    breakout_dist[15] = 1.0 - sum(breakout_dist[1:14])

    dcat = Categorical(breakout_dist)
    outs = rand(dcat, sum(inq))
    outs = [count(outs .== i) for i in 1:length(breakout_dist)][1:length(breakout_dist)-1]
end



"""
    function cnt_2_array(cnt, pop_mat)

Distribute a count of people to be unquarantined across
the cells of a population matrix, limited to the number of
quarantined people in each cell of the quarantine population
matrix.
"""
function cnt_2_array(cnt, pop_mat; ret_conds=[nil, mild, sick, severe, unexposed, recovered])
    (ls, _, as) = axes(pop_mat)
    map2unq = (unexposed=1, infectious=-1, recovered=2,  dead=-1, nil=3, mild=4, sick=5, severe=6)
    new_mat = zeros(Int, length(ls), length(ret_conds), length(as))
    for l in ls, c in ret_conds, a in as
        if cnt <= 0; break; end
        put = clamp(cnt,0,pop_mat[l,c,a])
        new_mat[l,map2unq[c],a] = put
        cnt -= put
    end
    return new_mat
end