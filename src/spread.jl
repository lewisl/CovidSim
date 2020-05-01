#############
# spread.jl
#############


#  stash for temporary values changed during simulation cases
#      to change just once and then get the originals back
#      it is the users responsibility to empty the stash
const spread_stash = Dict{Symbol, Array}()

"""
How far do the infectious people spread the virus to
previously unexposed people, by agegrp?  For a single locale...
"""
function spread!(locale, density_factor = [1.0]; spreadcases=[], dat=openmx, env=env)

    if ctr[:day]  == 1    # TODO when is the right time?  what is the right cleanup?
        cleanup_spread_stash()
    end

    # variables from env
    spreaders =         env.spreaders
    all_accessible =    env.all_accessible
    simple_accessible = env.simple_accessible
    numcontacts =       env.numcontacts
    numtouched =        env.numtouched

    # set function scope for variables modified in loop
    newinfected = zeros(Int, 5) # by agegrp

    # how many spreaders  TODO grab their condition.  Separate probs by condition
    spreaders[:] = grab(infectious_cases, agegrps, lags, locale, dat=dat) # laglim x 4 x 5 lag x cond x agegrp

    if sum(spreaders) == 0
        return
    end

    all_accessible[:] = grab([unexposed,recovered, nil, mild, sick, severe],agegrps,lags, locale, dat=dat)  #   laglim x 6 x 5  lag x cond by agegrp
    simple_accessible[:] = sum(all_accessible, dims=1)[1,:,:] # sum all the lags result (6,5)

    all_unexposed = grab(unexposed, agegrps, 1, locale, dat=dat)  # (5, ) agegrp for lag 1

    # set and run spreadcases
    spread_case_setter(spreadcases, env=env)  # bounces out right away if empty
    if iszero(env.sd_compliance) || isone(env.sd_compliance)  # no compliance or full compliance--no split, just run once
        newinfected = spreadsteps(density_factor, all_unexposed, env=env)
    else # we need to split people into complying and noncomplying
        newinfected = spread_case_runner(density_factor, all_unexposed, env=env)
    end  # no active case or active case

    lag1 = 1

    # test if newinfected > unexposed
    # this hasn't triggered for a long time
    # for agegrp in agegrps
    #     if newinfected[agegrp] > grab(unexposed, agegrp, lag1, locale, dat=dat)
    #         @warn "big problem: infected exceeds unexposed in spread!"
    #     end
    # end

    # move the people from unexposed:agegrp to infectious:agegrp and nil
    plus!.(newinfected, infectious, agegrps, lag1, locale, dat=dat)
    plus!.(newinfected, nil, agegrps, lag1, locale, dat=dat)
    minus!.(newinfected, unexposed, agegrps, lag1, locale, dat=dat)

    push!(spreadq, (day=ctr[:day], locale=locale,
                spreaders = sum(spreaders) + sum(get(spread_stash, :comply_spreaders, 0)),
                contacts = sum(numcontacts) + sum(get(spread_stash, :comply_contacts, 0)),
                touched = sum(numtouched) + sum(get(spread_stash, :comply_touched, 0)),
                accessible = sum(all_accessible),
                unexposed=sum(grab(unexposed, agegrps, lag1, locale, dat=dat)),
                infected=sum(newinfected)))

    return
end


function spreadsteps(density_factor, all_unexposed; env=env)
    how_many_contacts!(density_factor, env=env) # how many people are contacted based on spreader? updates env.numcontacts
    how_many_touched!(env=env)   # lag x agegrp # how many are touched based on recipient and contacts?
    newinfected = how_many_infected(all_unexposed, env=env)    # (5,)  # how many people become infected?
end


"""
How many contacts do spreaders attempt to reach?  This is based on the characteristics of the
spreaders.
"""
function how_many_contacts!(density_factor=1.0; env=env)
    #=  This originally ignores the conditions of the touched--assumes they are all equally likely to be touched
        how_many_touched corrects this.
        We assume spreaders is small compared to all_accessible. At some point this might not be true:
        how_many_touched also handles this.
    =#

    # variables from env
    spreaders = env.spreaders  # laglim,4,5
    numcontacts = env.numcontacts  # the result
    contact_factors = env.contact_factors
    all_accessible = env.all_accessible

    if sum(env.simple_accessible) == 0   # this can happen with a social distancing case with 100% compliance
        numcontacts[:] .= 0
        return
    end

    sp_lags, sp_conds, sp_ages = size(spreaders)

    # how many people are contacted by each spreader?  Think of this as reaching out...
        # numcontacts is the potential number of people contacted by a spreader in each
        # cell by lag (laglim), infectious cond (4), and agegrp(5)
    for agegrp in 1:sp_ages
        for cond in 1:sp_conds
            for lag in 1:sp_lags
                scale = contact_factors[cond, agegrp]
                spcount = spreaders[lag, cond, agegrp]
                if spcount == 0
                    numcontacts[lag, cond, agegrp] = 0
                    continue
                end
                dgamma = Gamma(1.0, density_factor * scale)  #shape, scale
                x = rand(dgamma,spcount)
                numcontacts[lag, cond, agegrp] = round(Int,sum(x))
            end
        end
    end

    # correct over contacting: this can happen with high scale towards the high point of infection
    oc_ratio = sum(numcontacts) / sum(all_accessible)
    if oc_ratio > 1.0
        println(ctr[:day],"warning: overcontact ratio ", oc_ratio)
        numcontacts[:] = round.(1.0/oc_ratio .* numcontacts)
    end
end


"""
For potential contacts by spreaders reaching out, how many of the accessible (susceptible and NOT)
are consequentially "touched" by a spreader? This is based on the characteristics of the
recipient. This also splits the recipients by the characteristics of the spreaders and the recipients.
"""
function how_many_touched!(; env=env)
#=
    - who is accessible: start with all to capture effect of "herd immunity", when it arises
    - all_accessible laglim x 6 x 5  lag x cond x agegrp, includes unexposed, recovered, nil, mild, sick, severe
    - numcontacts laglim x 4 x 5  lag x cond x agegrp, includes nil, mild, sick, severe = the infectious

    There is a lot of setup before we get to business here.
=#
    # variables from simulation environment
    numcontacts =       env.numcontacts         # outcome from how_many_contacts!
    simple_accessible = env.simple_accessible
    lag_contacts =      env.lag_contacts
    touch_factors =     env.touch_factors
    numtouched =        env.numtouched          # result by lag of spreader and agegrp or recipient

    if sum(simple_accessible) == 0   # this can happen with a social distancing case with 100% compliance
        numtouched[:] .= 0
        return
    end

    # map to access maps conditions to the rows of simple_accessible and touch_factors
    map2access = (unexposed= 1, infectious=-1, recovered= 2, dead=-1, nil= 3, mild=  4, sick= 5, severe= 6)

    lag_contacts[:] = sum(numcontacts,dims=(2,3))[:,:,1] # (laglim, ) contacts by lag after sum by cond, agegrp
    totaccessible = sum(simple_accessible)

    # sum accessible to unexposed, recovered, infectious by agegrps (don't use nil, mild, sick, severe conditions)
            # unexposed, recovered, sum of(nil, mild, sick, severe) = infectious  6,5 with last 3 rows all zeros
    simple_accessible[:] = [simple_accessible[1:2,:]; sum(simple_accessible[3:6,:],dims=1); zeros(Int,3,5)];  # (6, 5)
    # s_a_pct is dist. of accessible by agegrp and uexposed, recovered, infectious (15,)
    s_a_pct = round.(reshape(simple_accessible[1:3,:] ./ totaccessible, 15), digits=5) # % for each cell
    if !isapprox(sum(s_a_pct), 1.0)
        s_a_pct = s_a_pct ./ sum(s_a_pct) # normalize so sums to 1.0
    end

    # we only use this to makes sure we don't touch more people than there are
    unexp_by_contact_lag = zeros(Int, laglim,5) # lags, agegrps
    lagpct = zeros(laglim)
    lagpct[:] = lag_contacts ./ (sum(lag_contacts) + 1e-8)
    if sum(lagpct) == 0.0
        lagpct[:] = fill(1.0/laglim, laglim)
    end
    @assert isapprox(sum(lagpct), 1.0) "pct must sum to 1.0; got $(sum(lagpct))"
    for i in agegrps, j in lags
          unexp_by_contact_lag[j,i] = round(Int,simple_accessible[map2access.unexposed, i] * lagpct[j])
    end

    # now to business: who gets touched in unexposed by agegrp?
    #=
        - folks in each cell of numcontacts touch a sample of the accessible reduced by the touch factor
        - draw a categorical sample for each cell to distribute them across the contact categories
        - we consider the folks who get contacts in infectious and recovered, because that happens--reduces
           the number of unexposed who are touched
        - we throw out the folks in infectious and recovered and keep those in unexposed
        - we should care about not touching more than the number of unexposed?   (or accessible?)
        - env.touch_factors is (5,6): agegrps by unexposed, recovered, nil, sick, mild, severe
    =#

    # map to top 3 rows of simple_accessible
    mapi = (unexposed= 1, infectious=3, recovered=2, dead=-1, nil= -1, mild= -1, sick= -1, severe= -1)

    dcat = Categorical(s_a_pct) # categorical distribution by agegrp and unexposed, recovered, infectious
    numtouched[:] = zeros(Int, laglim, length(agegrps)) # (laglim,5)

    # loop over numcontacts lag vector
    for lag in lags
        lc = lag_contacts[lag]
        if lc == 0
            numtouched[lag,:] .= 0
            continue
        end
        x = rand(dcat, lc) # (length(lc),) probabistically distribute contacts for a lag across accessible by unexposed|recovered|infectious, agegrp
        peeps = reshape([count(x .== i) for i in 1:15], 3,5)[1,:]  # (5,) distribute across all 3 groups, but only take unexposed
        for a in agegrps # probabilistically see who of the accessible is significantly touched
            cnt = binomial_one_sample(peeps[a], touch_factors[mapi.unexposed, a])
            numtouched[lag, a] = clamp(cnt, 0, unexp_by_contact_lag[lag, a])
        end
    end
end


function how_many_infected(all_unexposed; env=env)
    #=
        - multiply the transmissibility of the spreader times the transmissibility of the contacts
                by lag for spreaders and by agegrp for the contacts
        - use the infection factor in a binomial sample:  was the contact "successful" in causing infection?
        - we'll test to be sure we don't exceed the unexposed and reduce touches to 80% of unexposed by agegrp
    =#

    # variables from env
    touched_by_lag_age = env.numtouched  # (laglim,5)  all_unexposed (5,)

    newinfected = zeros(Int, length(agegrps))  # (5,)

    if sum(env.numtouched) == 0   # this can happen with a social distancing case with 100% compliance
        return newinfected
    end

    # setup risk table   now part of initialization: caller must update if send or recv risk has changed
    # env.riskmx[:] = send_risk_by_recv_risk(env.send_risk_by_lag, env.recv_risk_by_age)  # (laglim,5)

    for age in agegrps
        for lag in lags
            newsick = binomial_one_sample(touched_by_lag_age[lag, age], env.riskmx[lag, age])  # draws, probability
            newsick = clamp(newsick, 0, floor(Int,.9 * all_unexposed[age]))
            newinfected[age] += newsick
        end
    end

    @debug "\n newly infected: $newinfected  \n"

    return newinfected  # (length of agegrps, ) (only condition is nil, assumed lag = 1 => first day infected)
end


function send_risk_by_recv_risk(send_risk, recv_risk)
    recv_risk' .* send_risk  # (laglim, agegrps)
end


function cleanup_spread_stash()
    for k in keys(spread_stash)
        delete!(spread_stash, k)
    end
end


function r0_sim(;sa_pct=[1.0,0.0,0.0], density_factor=1.0, dt=[], cf=[], tf=[],
                compliance=[], shift_contact=(), shift_touch=(), pri=false)
    # factor_source must be one of: r0env, or env of current simulation
    # setup separate environment
    r0env = initialize_sim_env();
    r0mx = data_dict(1)  # single locale
    locale = 1
    popsize = 1_000_000
    for agegrp in agegrps
        r0mx[locale][1, unexposed, agegrp] = floor(Int,age_dist[agegrp] * popsize)
    end
    track_infected = zeros(Int, 5)

    # setup data
    all_unexposed = grab(unexposed, agegrps, 1, locale, dat=r0mx)  # (5, ) agegrp for lag 1

    r0env.all_accessible[:] = grab([unexposed,recovered, nil, mild, sick, severe], agegrps, lags, locale, dat=r0mx)  #   laglim x 6 x 5  lag x cond by agegrp
    r0env.simple_accessible[:] = sum(r0env.all_accessible, dims=1)[1,:,:] # sum all the lags result (6,5)
    if !isempty(compliance)
        r0env.simple_accessible[:] = round.(compliance .* r0env.simple_accessible)
    end

    if sa_pct[1] != 1.0
        sa_pct = [sa_pct[1],sa_pct[2],sa_pct[3], fill(sa_pct[3]./4.0, 3)...]
        res = [r0env.simple_accessible[1,:] .* i for i in sa_pct]
        sanew = zeros(6,5)
        for i in 1:6
           sanew[i,:] .= res[i]
        end
        r0env.simple_accessible[:] = round.(Int, sanew)
    end

    spreaders = ones(Int,25,4,5)
    spreaders[1,:,:] .= 0;
    r0env.spreaders = spreaders

    # parameters that drive r0
    !isempty(cf) && (r0env.contact_factors = deepcopy(cf))
    !isempty(tf) && (r0env.touch_factors = deepcopy(tf))
    isempty(shift_contact)  || (r0env.contact_factors =shifter(r0env.contact_factors, shift_contact...))
    isempty(shift_touch) || (r0env.touch_factors = shifter(r0env.touch_factors, shift_touch...))

    for i = 1:laglim
        isempty(dt) || (transition!(dt, locale,dat=r0mx)) # optionally start with a transition

        how_many_contacts!(density_factor, env=r0env)

        how_many_touched!(env=r0env)

        track_infected .+= how_many_infected(all_unexposed, env=r0env)
    end

    # for last iteration, but close across most of them
    tot_spreaders = sum(r0env.spreaders)
    tot_contacts = sum(r0env.numcontacts)
    tot_touched = sum(r0env.numtouched)
    tot_infected = sum(track_infected)
    r0 = tot_infected / tot_spreaders

    if pri
        contact_factors = round.(r0env.contact_factors, digits=3)
        touch_factors = round.(r0env.touch_factors, digits=3)

            println("r0 = $r0")
            println("spreaders = $tot_spreaders, contacts = $tot_contacts, touched = $tot_touched, infected = $tot_infected")
            print("contact_factors ")
            display(contact_factors)
            print("touch_factors ")
            display(touch_factors) 
    end
    return(r0=r0,spreaders=tot_spreaders,contacts=tot_contacts,touched=tot_touched,infected=tot_infected)
end


function r0_table(n=6, cfstart = 0.9, tfstart = 0.3)
    tbl = zeros(n+1,n+1)
    cfiter = [cfstart + (i-1) * .1 for i=1:n]
    tfiter = [tfstart + (i-1) * 0.05 for i=1:n]
    for (j,cf) in enumerate(cfiter)
        for (i,tf) = enumerate(tfiter)
            tbl[i+1,j+1] = r0_sim(shift_contact=(0.2,cf), shift_touch=(.18,tf)).r0
        end
    end
    tbl[1, 2:n+1] .= cfiter
    tbl[2:n+1, 1] .= tfiter
    tbl[:] = round.(tbl, digits=2)
    display(tbl)
    return tbl
end

#=
approximate r0 values from model (simulation uses samples)
           max c_f
            1.1   1.2   1.3   1.4   1.5   1.6   1.7   1.8
max t_f    -----------------------------------------------
     0.25 | 0.57  0.58  0.63  0.73  0.72  0.79  0.74  0.87
     0.3  | 0.59  0.7   0.76  0.74  0.91  1.0   0.95  1.01
     0.35 | 0.62  0.82  0.8   0.99  0.93  1.11  1.11  1.12
     0.4  | 0.82  0.88  0.98  1.08  0.99  1.31  1.38  1.33
     0.45 | 0.94  1.02  1.1   1.12  1.23  1.35  1.39  1.52
     0.5  | 1.03  1.14  1.12  1.24  1.44  1.45  1.5   1.58
     0.55 | 1.07  1.19  1.24  1.52  1.55  1.46  1.67  1.7
     0.6  | 1.18  1.33  1.45  1.44  1.53  1.82  1.78  2.07
=#