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
    newinfected = zeros(Int, 5)

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

    # test if newinfected > unexposed
    lag = 1
    for agegrp in agegrps
        if newinfected[agegrp] > grab(unexposed, agegrp, lag, locale, dat=dat)
            @warn "big problem: infected exceeds unexposed in spread!"
        end
    end

    # move the people from unexposed:agegrp to infectious:agegrp and nil
    plus!.(newinfected, infectious, agegrps, lag, locale, dat=dat)
    plus!.(newinfected, nil, agegrps, lag, locale, dat=dat)
    minus!.(newinfected, unexposed, agegrps, lag, locale, dat=dat)

    push!(spreadq, (day=ctr[:day], locale=locale,
                spreaders = sum(spreaders) + sum(get(spread_stash, :comply_spreaders, 0)),
                contacts = sum(numcontacts) + sum(get(spread_stash, :comply_contacts, 0)),
                touched = sum(numtouched) + sum(get(spread_stash, :comply_touched, 0)),
                accessible = sum(all_accessible),
                unexposed=sum(grab(unexposed, agegrps, lag, locale, dat=dat)),
                infected=sum(newinfected)))
    # add to stats queue for today
    queuestats(sum(newinfected), locale, spreadstat) # sum(5 agegroups), nil is the default, single locale

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


function spread_sanity(x; shape= 1.0,scale=1.6,pr=0.6) # defaults are before social distancing or isolation
    send_risk_by_lag = [.1, .3, .7, .8, .9, .9, .8, .7, .6, .5, .3, .1, .1, 0.05, 0.05, 0.05, 0, 0, 0]
    recv_risk_by_age = [.1, .4, .4, .50, .55]    
    risk = mean(recv_risk_by_age' .* send_risk_by_lag) 
    dgamma = Gamma(shape, scale)
    contact(x) = round(Int,sum(rand(dgamma,x)))
    touch(x) = rand(Binomial(x, pr))
    infect(x) = round(Int, risk .* x)
    samp = [reduce(+,[infect(touch(contact(x))) for i in 1:14]) for i in 1:50] # sample of newinfected
    mean_contacts = round(Int, mean([contact(x) for i in 1:50]))
    mean_touches = round(Int, mean([touch(mean_contacts) for i in 1:50]))
    mean_infected = round(Int, mean(samp))
    return (r0=mean(samp) / x, mean_contacts=mean_contacts, mean_touches=mean_touches, mean_infected=mean_infected)
end

#=
    R0 tables
    arithmetic mean of risk in the model is .118

tables by pr
for pr = 0.25    risk by scale
    risk   scale->
     0.0   1.0   1.2   1.4   1.6   1.8   2.0
     0.15  0.45  0.54  0.63  0.72  0.81  0.9
     0.25  0.75  0.89  1.05  1.2   1.35  1.5
     0.35  1.06  1.26  1.47  1.68  1.9   2.09
     0.45  1.35  1.63  1.89  2.16  2.45  2.69
     0.55  1.66  1.98  2.32  2.63  2.98  3.29
     0.65  1.96  2.33  2.73  3.13  3.53  3.9

for pr = 0.35    risk by scale
    risk   scale->
     0.0   1.0   1.2   1.4   1.6   1.8   2.0
     0.15  0.63  0.75  0.89  1.01  1.14  1.27
     0.25  1.05  1.26  1.47  1.68  1.89  2.1
     0.35  1.47  1.76  2.06  2.35  2.66  2.93
     0.45  1.9   2.26  2.66  3.04  3.42  3.79
     0.55  2.31  2.78  3.25  3.67  4.16  4.64
     0.65  2.73  3.29  3.81  4.34  4.96  5.46

for pr = 0.45    risk by scale
    risk   scale->
     0.0   1.0   1.2   1.4   1.6   1.8   2.0
     0.15  0.81  0.97  1.13  1.29  1.46  1.61
     0.25  1.35  1.62  1.89  2.16  2.42  2.69
     0.35  1.89  2.26  2.64  3.02  3.38  3.77
     0.45  2.43  2.92  3.4   3.9   4.37  4.85
     0.55  2.95  3.55  4.15  4.72  5.35  5.94
     0.65  3.51  4.22  4.89  5.65  6.34  7.0


for pr = 0.45    risk by scale
    risk   scale->
     0.0   1.0   1.2   1.4   1.6   1.8   2.0
     0.15  0.99  1.19  1.38  1.58  1.78  1.98
     0.25  1.65  1.98  2.32  2.64  2.98  3.29
     0.35  2.32  2.77  3.23  3.7   4.16  4.63
     0.45  2.97  3.54  4.15  4.75  5.35  5.95
     0.55  3.63  4.34  5.09  5.82  6.56  7.27
     0.65  4.3   5.15  6.0   6.82  7.73  8.61

for pr = 0.65    risk by scale
    risk   scale->
     0.0   1.0   1.2   1.4   1.6   1.8    2.0
     0.15  1.17  1.4   1.64  1.87  2.1    2.35
     0.25  1.95  2.34  2.72  3.12  3.5    3.91
     0.35  2.72  3.29  3.82  4.38  4.93   5.47
     0.45  3.51  4.22  4.93  5.59  6.33   7.04
     0.55  4.29  5.17  6.02  6.88  7.72   8.58
     0.65  5.05  6.09  7.08  8.13  9.18  10.18

for pr = 0.75    risk by scale
    risk   scale->
     0.0   1.0   1.2   1.4   1.6    1.8    2.0
     0.15  1.35  1.62  1.88  2.16   2.43   2.7
     0.25  2.25  2.72  3.15  3.62   4.04   4.5
     0.35  3.15  3.77  4.41  5.02   5.66   6.31
     0.45  4.06  4.85  5.65  6.5    7.29   8.08
     0.55  4.97  5.94  6.95  7.92   8.95   9.91
     0.65  5.86  7.04  8.2   9.36  10.54  11.67


for pr = 0.85    risk by scale
    risk   scale->
     0.0   1.0   1.2   1.4    1.6    1.8    2.0
     0.15  1.52  1.84  2.14   2.46   2.75   3.06
     0.25  2.54  3.05  3.58   4.08   4.6    5.11
     0.35  3.57  4.3   4.98   5.71   6.43   7.15
     0.45  4.58  5.5   6.43   7.35   8.26   9.16
     0.55  5.6   6.74  7.83   8.97  10.11  11.22
     0.65  6.62  7.95  9.27  10.57  11.92  13.26


tables by risk
for risk = .15
    pr     scale->
     0.0   1.0   1.2   1.4   1.6   1.8   2.0   2.2   2.5
     0.15  0.27  0.33  0.38  0.43  0.48  0.54  0.59  0.67
     0.25  0.45  0.54  0.63  0.72  0.81  0.9   1.0   1.12
     0.35  0.63  0.75  0.89  1.01  1.13  1.26  1.38  1.58
     0.45  0.81  0.97  1.13  1.3   1.46  1.62  1.78  2.03
     0.55  0.99  1.18  1.39  1.58  1.78  1.98  2.17  2.47
     0.65  1.17  1.4   1.65  1.87  2.11  2.34  2.58  2.92
     0.75  1.35  1.62  1.89  2.16  2.44  2.69  2.97  3.37
     0.85  1.53  1.84  2.14  2.44  2.75  3.06  3.37  3.82

for risk = .15
    pr     scale->
     0.0   1.0   1.2   1.4   1.6   1.8   2.0   2.2   2.5
     0.15  0.45  0.54  0.63  0.72  0.81  0.9   0.99  1.12
     0.25  0.75  0.9   1.04  1.2   1.35  1.49  1.66  1.88
     0.35  1.05  1.27  1.47  1.67  1.88  2.1   2.31  2.63
     0.45  1.35  1.61  1.89  2.16  2.44  2.69  2.96  3.39
     0.55  1.65  1.98  2.31  2.64  2.98  3.3   3.62  4.13
     0.65  1.95  2.33  2.72  3.13  3.51  3.89  4.29  4.88
     0.75  2.25  2.7   3.14  3.61  4.06  4.5   4.94  5.63
     0.85  2.55  3.06  3.58  4.08  4.58  5.1   5.6   6.37

=#
