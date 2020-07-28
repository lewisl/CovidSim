#############
# spread.jl
#############

"""
Stash for temporary values changed during simulation cases
- to change just once and then get the originals back
- it is the users responsibility to empty the stash
- there may be (will be) side effects if you don't empty the stash between simulations
"""
const spread_stash = Dict{Symbol, Array}()

"""
How far do the infectious people spread the virus to
previously unexposed people, by agegrp?  For a single locale...
"""
function spread!(locale, density_factor = [1.0]; spreadcases=[], dat=openmx, env=env)

    if ctr[:day]  == 1    # TODO when is the right time?  what is the right cleanup?
        cleanup_stash(spread_stash)
    end

    # variables from env
    spreaders =         env.spreaders
    all_accessible =    env.all_accessible
    simple_accessible = env.simple_accessible
    contacts =          env.contacts
    touched =           env.touched
    lag1 =              1

    # set function scope for variables modified in loop--> this is the result
    newinfected = zeros(T_int[], 5) # by agegrp

    # how many spreaders  TODO grab their condition.  Separate probs by condition
    # spreaders[:] = grab(infectious_cases, agegrps, lags, locale, dat=dat) # laglim x 4 x 5 lag x cond x agegrp
    spreaders[:] = grab(infectious_cases, agegrps, lags, locale, dat=dat) # laglim x 4 x 5 lag x cond x agegrp

    if sum(spreaders) == T_int[](0)
        return
    end

    all_accessible[:] = grab([unexposed,recovered, nil, mild, sick, severe],agegrps,lags, locale, dat=dat)
    simple_accessible[:] = sum(all_accessible, dims=1)[1,:,:] # sum all the lags result (6,5)  @views 

    all_unexposed = grab(unexposed, agegrps, lag1, locale, dat=dat)  # (5, ) agegrp for lag 1

    # set and run spreadcases or just run spreadsteps
    spread_case_setter(spreadcases, env=env)  # bounces out right away if empty
    if iszero(env.sd_compliance) || isone(env.sd_compliance)  # no compliance or full compliance--no split, just run once
        newinfected = spreadsteps(density_factor, all_unexposed, env=env)
    else # we need to split people into complying and noncomplying
        newinfected = spread_case_runner(density_factor, all_unexposed, env=env)
    end  # no active case or active case

    # move the people from unexposed:agegrp to infectious:agegrp and nil
    plus!.(newinfected, infectious, agegrps, lag1, locale, dat=dat)
    plus!.(newinfected, nil, agegrps, lag1, locale, dat=dat)
    minus!.(newinfected, unexposed, agegrps, lag1, locale, dat=dat)

    push!(spreadq, (day=ctr[:day], locale=locale,
                spreaders = sum(spreaders) + sum(get(spread_stash, :comply_spreaders, 0)),
                contacts = sum(contacts) + sum(get(spread_stash, :comply_contacts, 0)),
                touched = sum(touched) + sum(get(spread_stash, :comply_touched, 0)),
                accessible = sum(all_accessible),
                unexposed=sum(grab(unexposed, agegrps, lag1, locale, dat=dat)),
                infected=sum(newinfected)))

    return
end


"""
function spreadsteps(density_factor, all_unexposed; env=env)

Perform the 3 fundamental steps that spread the virus to unexposed people.
Returns newinfected.
"""
function spreadsteps(density_factor, all_unexposed; env=env)

    how_many_contacts!(density_factor, env=env)

    how_many_touched!(;env=env)

    newinfected = how_many_infected(all_unexposed, env=env)    # (5,)  # how many people become infected?
end


# this method used spread; full method used by test_and_trace
function how_many_contacts!(density_factor=1.0; env=env)
    how_many_contacts!(env.contacts, env.spreaders, env.simple_accessible, env.contact_factors, density_factor; env=env)
end


"""
How many contacts do spreaders attempt to reach?  This is based on the characteristics of the
spreaders.
"""
function how_many_contacts!(contacts, spreaders, target_accessible, contact_factors, 
    density_factor; env=env)
    #=  how_many_contacts--assumes anyone not isolated is equally likely to be touched.
        how_many_touched corrects this by allocating contacts in each cell by population proportion.
        Spreaders is usually small compared to all_accessible but there is a check.
    =#

    if sum(env.simple_accessible) == T_int[](0)   # this can happen with a social distancing case with 100% compliance
        contacts[:] .= T_int[](0)
        return
    end

    # how many people are contacted by each spreader?  Think of this as reaching out...
        # contacts is the potential number of people contacted by a spreader in each
        # cell by lag (laglim), infectious cond (4), and agegrp(5)
    sp_lags, sp_conds, sp_ages = size(spreaders)
    for agegrp in 1:sp_ages
        for cond in 1:sp_conds
            @inbounds for lag in 1:sp_lags
                scale = contact_factors[cond, agegrp]
                spcount = spreaders[lag, cond, agegrp]
                if spcount == T_int[](0)
                    contacts[lag, cond, agegrp] = T_int[](0)
                    continue
                end
                dgamma = Gamma(1.0, density_factor * scale)  #shape, scale
                x = rand(dgamma,spcount)
                contacts[lag, cond, agegrp] = round(T_int[], sum(x))
            end
        end
    end

    # correct over contacting: may rarely happen with high scale around the high point of infection
    oc_ratio = sum(contacts) / sum(env.all_accessible)
    if oc_ratio > 1.0
        @warn "day: $(ctr[:day]) warning: overcontact ratio > 1.0: $oc_ratio"
        contacts[:] = round.(T_int[], 1.0/oc_ratio .* contacts)
    end

    return contacts
end

 
"""
How many contacts result in consequential touches? Each cell in contacts is the number of
contacts indexed by the characteristics of the *spreader*. Map each cell of contacts to
the accessible population (susceptible and NOT) proportionally to the cells of the accessible. 
Then, filter through binomial distribution for a "successful" touch based on characteristics 
of the recipient.

An alternative method is used for test_and_trace, which prepares its inputs in the caller.
"""
function how_many_touched!(;env=env)

    map2touch = (unexposed= 1, infectious=3, recovered=2, dead=-1, 
                 nil= -1, mild= -1, sick= -1, severe= -1)
    # consolidate the cells of the contacts made by spreaders by summing across agegrp and condition
    # reduce the number of cells 500 / (4 * 5) → 25
    env.lag_contacts[:] = sum(env.contacts,dims=(2,3))[:,:,1] #  @views (laglim, ) contacts by lag after sum by cond, agegrp
    contacts = reshape(env.lag_contacts, laglim, 1, 1)  # laglim x 4 x 5: lag x cond x agegrp; cond in {nil, mild, sick, severe}
    touched = env.touched
    touched .= T_int[](0)
    target_accessible = env.simple_accessible  # (6,5) 6 conds: unexp, recovered, nil, mild, sick, severe by agegrps
    target_tf = env.touch_factors  # view(env.touch_factors, map2touch.unexposed, :)  w/o view (6,5)

    totaccessible = sum(target_accessible)
    peeps = zeros(T_int[], size(target_accessible))

    # this can happen with a social distancing or quarantine case with 100% compliance
        # or early/late in epidemic when there are no spreaders to make contacts
    if totaccessible == T_int[](0)  || sum(contacts) == T_int[](0)
        return touched  # already initialized to zeros by caller
    end

    # map to access maps conditions to the rows of simple_accessible and touch_factors
    map2access = (unexposed= 1, infectious=-1, recovered= 2, dead=-1,
                  nil= 3, mild=  4, sick= 5, severe= 6)
    # map to top 3 rows of simple_accessible
    map2touch = (unexposed= 1, infectious=3, recovered=2, dead=-1,
                 nil= -1, mild= -1, sick= -1, severe= -1)

    # t_a_pct is dist. of accessible by agegrp and target conds (15,)
    t_a_pct = round.(reshape(target_accessible ./ totaccessible, length(target_accessible)), digits=10) # % for each cell

    if !isapprox(sum(t_a_pct), 1.0)
        @debug @warn "how_many_touched: sum(t_a_pct) $(sum(t_a_pct)) not equal 1.0: normalized to sum to 1.0"
        t_a_pct[:] = t_a_pct ./ sum(t_a_pct) # normalize to sum to 1.0
    end

    dcat = Categorical(t_a_pct) # categorical distribution by cells of target_accessible

    for l in eachindex(contacts)
        subgroup = contacts[l] # source of contacts, broken down by shape of spreaders
        if subgroup == 0
            continue
        else
            x = rand(dcat, subgroup)  # distribute cell across all accessible cells
            peeps[:] = reshape([count(x .== i) for i in 1:length(dcat.p)], size(target_accessible))  
            cnt = binomial_one_sample.(peeps, target_tf)  # in each cell of accessible: successful touch
            touched[l, :, :] .= cnt
        end
    end

    return touched
end


# method used by test_and_trace
function how_many_touched!(touched, contacts, target_accessible, target_conds,
                           target_tf; env=env)

    totaccessible = sum(target_accessible)
    peeps = zeros(T_int[], size(target_accessible))

    # this can happen with a social distancing or quarantine case with 100% compliance
        # or early/late in epidemic when there are no spreaders to make contacts
    if totaccessible == 0  || sum(contacts) == 0
        return touched  # already initialized to zeros by caller
    end

    # map to access maps conditions to the rows of simple_accessible and touch_factors
    map2access = (unexposed= 1, infectious=-1, recovered= 2, dead=-1,
                  nil= 3, mild=  4, sick= 5, severe= 6)
    # map to top 3 rows of simple_accessible
    map2touch = (unexposed= 1, infectious=3, recovered=2, dead=-1,
                 nil= -1, mild= -1, sick= -1, severe= -1)

    # t_a_pct is dist. of accessible by agegrp and target conds (15,)
    t_a_pct = round.(reshape(target_accessible ./ totaccessible, length(target_accessible)), digits=10) # % for each cell

    if !isapprox(sum(t_a_pct), 1.0)
        @debug @warn "how_many_touched: sum(t_a_pct) $(sum(t_a_pct)) not equal 1.0: normalized to sum to 1.0"
        t_a_pct[:] = t_a_pct ./ sum(t_a_pct) # normalize to sum to 1.0
    end

    dcat = Categorical(t_a_pct) # categorical distribution by dims of target_accessible

    for l in lags  
        for cond in target_conds
            @inbounds for a in agegrps     # agegrps 
                subgroup = contacts[l, map2access[cond], a]   # [l, map2access[cond], a]
                if subgroup == T_int[](0)
                    continue
                else   # distribute 1 contact subgroup to all cells of touched
                    x = rand(dcat, subgroup) # across all the touch target groups
                    peeps[:] = reshape([count(x .== i) for i in 1:length(dcat.p)], size(target_accessible))
                    @inbounds for l in lags  # consolidate the touched
                        @views cnt = binomial_one_sample.(peeps[l,:,:], target_tf)  
                        touched[l,:,:] .+= cnt  # total the results from each contact cell
                    end
                end
            end
        end 
    end

    return touched
end


"""
function how_many_infected(all_unexposed; env=env)

This function is the last step of spreadsteps. Based on previous calculations of contacts
and touched:
- multiply the transmissibility of the spreader times the transmissibility of the touched
    by lag for spreaders and by agegrp for the touched
- use the infection factor in a binomial sample:  was the contact "successful" in causing infection?
- test to be sure we don't exceed the unexposed and reduce touches to 95% of unexposed by agegrp

returns newinfected
"""
function how_many_infected(all_unexposed; env=env)

    # only unexposed (= susceptible) can become infected
    @views touched_by_lag_age = env.touched[:, unexposed, :]  # (laglim,5)

    newinfected = zeros(T_int[], length(agegrps))  # (5,)

    if sum(env.touched) == T_int[](0)   # might happen with a social distancing or quarantine case with 100% compliance
        return newinfected
    end

    for age in agegrps
        for lag in lags
            newsick = binomial_one_sample(touched_by_lag_age[lag, age], env.riskmx[lag, age])  # draws, probability
            newinfected[age] += newsick
        end
        # newsick = clamp(newinfected[age], T_int[](0), floor(T_int[],.95 * all_unexposed[age]))
    end

    return newinfected  # (length of agegrps, ) (only condition is nil, assumed lag = 1 => first day infected)
end


function send_risk_by_recv_risk(send_risk, recv_risk)
    recv_risk' .* send_risk  # (laglim, agegrps)
end


function cleanup_stash(stash)
    for k in keys(stash)
        delete!(stash, k)
    end
end


function r0_sim(;env=env, sa_pct=[1.0,0.0,0.0], density_factor=1.0, dt=[], decpoints=[], cf=[], tf=[],
                compliance=[], shift_contact=(), shift_touch=(), disp=false)
    # factor_source must be one of: r0env, or env of current simulation
    # setup separate environment
    r0env = initialize_sim_env(env.geodata; contact_factors=env.contact_factors, touch_factors=env.touch_factors,
                               send_risk=env.send_risk_by_lag, recv_risk=env.recv_risk_by_age);
    r0mx = data_dict(1; lags=laglim, conds=length(conditions), agegrps=ages)  # single locale
    locale = 1
    population = convert(T_int[], 2_000_000)
    setup_unexposed!(r0mx, population, locale)

    # setup data
    all_unexposed = grab(unexposed, agegrps, 1, locale, dat=r0mx)  # (5, ) agegrp for lag 1
    track_infected = zeros(T_int[], 5)
    track_contacts = zeros(T_int[], laglim, 4, 5)
    track_touched = zeros(T_int[], laglim, 6, 5)

    r0env.all_accessible[:] = grab([unexposed,recovered, nil, mild, sick, severe], agegrps, lags, locale, dat=r0mx)  #   laglim x 6 x 5  lag x cond by agegrp
    r0env.simple_accessible[:] = sum(r0env.all_accessible, dims=1)[1,:,:] # sum all the lags result (6,5)
    if !isempty(compliance)
        r0env.simple_accessible[:] = round.(T_int[], compliance .* r0env.simple_accessible)
    end

    if sa_pct[1] != 1.0
        sa_pct = [sa_pct[1],sa_pct[2],sa_pct[3], fill(sa_pct[3]./4.0, 3)...]
        res = [r0env.simple_accessible[1,:] .* i for i in sa_pct]
        sanew = zeros(T_int[], 6, 5)
        for i in 1:6
           sanew[i,:] .= res[i]
        end
        r0env.simple_accessible[:] = round.(T_int[], sanew)
    end

    age_relative = round.(T_int[], age_dist ./ minimum(age_dist))
    starting_spreaders = ones(T_int[], laglim,4,5)
    for i in 1:5
        starting_spreaders[:,:,i] .= age_relative[i]
    end
    if !isempty(dt)
        starting_spreaders[2:laglim, :, :] .= T_int[](0)
        starting_spreaders .*= T_int[](20)
    else
        starting_spreaders[1,:,:] .= T_int[](0);
    end
    r0env.spreaders = copy(starting_spreaders)
    input!(starting_spreaders,infectious_cases,agegrps,lags,locale,dat=r0mx)

    # parameters that drive r0
    !isempty(cf) && (r0env.contact_factors = deepcopy(cf))
    !isempty(tf) && (r0env.touch_factors = deepcopy(tf))
    isempty(shift_contact)  || (r0env.contact_factors =shifter(r0env.contact_factors, shift_contact...))
    isempty(shift_touch) || (r0env.touch_factors = shifter(r0env.touch_factors, shift_touch...))

    stopat = !isempty(dt) ? laglim : 1

    for i = 1:stopat
        disp && println("test day = $i, spreaders = $(sum(r0env.spreaders))")

        track_contacts .+= how_many_contacts!(density_factor, env=r0env)
        track_touched .+= how_many_touched!(env=r0env)
        track_infected .+= how_many_infected(all_unexposed, env=r0env)

        if !isempty(dt)  # optionally transition
            transition!(dt, all_decpoints, locale, dat=r0mx)
            r0env.spreaders = copy(grab(infectious_cases,agegrps,lags,locale,dat=r0mx))
        end
    end

    if !isempty(dt)
        tot_spreaders = sum(starting_spreaders)
    else
        tot_spreaders = round.(T_int[], sum(starting_spreaders) / 24)
    end
    tot_contacts = sum(track_contacts)
    tot_touched = sum(track_touched)
    tot_infected = sum(track_infected)
    r0 = tot_infected / tot_spreaders
    contact_ratio = tot_contacts / tot_spreaders  
    touch_ratio = tot_touched / tot_spreaders

    if disp
        contact_factors = round.(r0env.contact_factors, digits=3)
        touch_factors = round.(r0env.touch_factors, digits=3)
            println("r0 = $r0  contact_ratio=$contact_ratio  touch_ratio=$touch_ratio")
            println("spreaders = $tot_spreaders, contacts = $tot_contacts, touched = $tot_touched, infected = $tot_infected")
            print("contact_factors ")
            display(contact_factors)
            print("touch_factors ")
            display(touch_factors)
    end
    
    ret = (day=ctr[:day], r0=r0, spreaders=tot_spreaders, contacts=tot_contacts, touched=tot_touched, infected=tot_infected,
           contact_ratio=contact_ratio,  touch_ratio=touch_ratio)
    push!(r0q, ret)
    return ret
end


function r0_table(n=6, cfstart = 0.9, tfstart = 0.3; env=env, dt=dt)
    tbl = zeros(n+1,n+1)
    cfiter = [cfstart + (i-1) * .1 for i=1:n]
    tfiter = [tfstart + (i-1) * 0.05 for i=1:n]
    for (j,cf) in enumerate(cfiter)
        for (i,tf) = enumerate(tfiter)
            tbl[i+1,j+1] = r0_sim(env=env, dt=dt, decpoints=decpoints, shift_contact=(0.2,cf), shift_touch=(.18,tf)).r0
        end
    end
    tbl[1, 2:n+1] .= cfiter
    tbl[2:n+1, 1] .= tfiter
    tbl[:] = round.(tbl, digits=2)
    display(tbl)
    return tbl
end

#=
approximate r0 values from model
using default age distribution
model selects a c_f based on age and infectious case
model selects a t_f based on age and condition (includes unexposed and recovered)
r0 depends on the selection of both c_f and t_f
Note: simulation uses samples so generated values will vary

           c_f
  tf       1.1   1.2   1.3   1.4   1.5   1.6   1.7   1.8   1.9   2.0
           ----------------------------------------------------------
     0.18 | 0.38| 0.38 | 0.42 | 0.46 | 0.49 | 0.51 | 0.55 | 0.57 | 0.59 | 0.64
     0.23 | 0.47| 0.47 | 0.49 | 0.55 | 0.64 | 0.65 | 0.68 | 0.68 | 0.73 | 0.77
     0.28 | 0.53| 0.61 | 0.62 | 0.65 | 0.69 | 0.73 | 0.79 | 0.82 | 0.83 | 0.88
     0.33 | 0.61| 0.66 | 0.7  | 0.79 | 0.8  | 0.83 | 0.9  | 0.95 | 0.99 | 1.04
     0.38 | 0.7 | 0.74 | 0.85 | 0.84 | 0.94 | 0.98 | 1.04 | 1.08 | 1.11 | 1.17
     0.43 | 0.8 | 0.85 | 0.89 | 0.93 | 1.03 | 1.11 | 1.16 | 1.2  | 1.27 | 1.34
     0.48 | 0.88| 0.91 | 0.99 | 1.03 | 1.16 | 1.23 | 1.26 | 1.32 | 1.42 | 1.47
     0.53 | 0.97| 1.06 | 1.08 | 1.18 | 1.26 | 1.27 | 1.42 | 1.47 | 1.52 | 1.61
     0.58 | 1.01| 1.09 | 1.17 | 1.25 | 1.33 | 1.43 | 1.52 | 1.52 | 1.68 | 1.76
     0.63 | 1.11| 1.2  | 1.25 | 1.38 | 1.42 | 1.5  | 1.65 | 1.75 | 1.78 | 1.95


=#|
