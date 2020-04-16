#############
# spread.jl
#############


#  stash for temporary values changed during simulation cases
#      to change just once and then get the originals back
#      it is the users responsibility to get rid of stuff
const spread_stash = Dict{Symbol, Array}()

"""
How far do the infectious people spread the virus to
previously unexposed people, by agegrp?  For a single locale...
"""
function spread!(locale, density_factor = [1.0]; spreadcases=[], dat=openmx, env=env)

    if ctr[:day]  == 1    # TODO when is the right time?  what is the right cleanup?
        cleanup_spread_cases()
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
    spreaders[:] = grab(infectious_cases, agegrps, lags, locale, dat=dat) # 19 x 4 x 5 lag x cond x agegrp

    if sum(spreaders) == 0
        return
    end

    all_accessible[:] = grab([unexposed,recovered, nil, mild, sick, severe],agegrps,lags, locale, dat=dat)  #   19 x 6 x 5  lag x cond by agegrp
    simple_accessible[:] = sum(all_accessible, dims=1)[1,:,:] # sum all the lags result (6,5)
    all_unexposed = grab(unexposed, agegrps, 1, locale, dat=dat)  # (5, ) agegrp for lag 1

    # set and run spreadcases
    case_setter(spreadcases, env=env)  # bounces out right away if empty
    if iszero(env.sd_compliance) || isone(env.sd_compliance)  # no compliance or full compliance--no split, just run once
        newinfected = spreadsteps(density_factor, all_unexposed, env=env)
    else # TODO this wants to be a case_runner function
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
            end
        end  # for loop
        # total values for comply + nocomply
        newinfected = newinfected[1] .+ newinfected[2]
        spreaders = spread_stash[:spreaders]
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

    push!(dayq, (day=ctr[:day], locale=locale, spreaders = sum(spreaders), 
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
How many contacts do spreaders attempt to make?  This is based on the characteristics of the
spreaders.
"""
function how_many_contacts!(density_factor=1.0; env=env)
    #=  This originally ignores the conditions of the touched--assumes they are all equally likely to be touched
        how_many_touched corrects this.
        We assume spreaders is small compared to all_accessible. At some point this might not be true:
        how_many_touched also handles this.
    =#

    # variables from env
    spreaders = env.spreaders  # 19,4,5
    numcontacts = env.numcontacts
    contact_factors = env.contact_factors
    all_accessible = env.all_accessible

    sp_lags, sp_conds, sp_ages = size(spreaders)
    # numcontacts = zeros(Int, sp_lags, sp_conds, sp_ages)  # 19 x 4 x 5 lag x cond x agegrp

    # how many people are contacted by each spreader?  Think of this as reaching out...
        # numcontacts is the potential number of people contacted by a spreader in each
        # cell by lag (19), infectious cond (4), and agegrp(5)
    for agegrp in 1:sp_ages
        for cond in 1:sp_conds
            for lag in 1:sp_lags
                scale = contact_factors[cond, agegrp]
                spcount = spreaders[lag, cond, agegrp]
                lag_reduce = .038 * 19  # reduce contact scale by lag between 1.0 and 0.3 -- WHY
                dgamma = Gamma(1.0, lag_reduce * density_factor * scale)  #shape, scale
                x = rand(dgamma,spcount)

                if isempty(x)
                    x=0
                end
                numcontacts[lag, cond, agegrp] = round(sum(x)) 
            end
        end
    end

    # correct over contacting
    oc_ratio = sum(numcontacts) / sum(all_accessible)
    if oc_ratio > 1.0
        println(ctr[:day]," overcontact ratio ", oc_ratio)
        numcontacts[:] = round.(1.0/oc_ratio .* numcontacts)
    end
    
    # return numcontacts
end


"""
For potential contacts by spreaders reaching out, how many of the accessible (susceptible and NOT)
are actuallly "touched" by a spreader? This is loosely based on the characteristics of the
receiver.
"""
function how_many_touched!(; env=env)
#=
    - who is accessible: start with all to capture effect of "herd immunity", when it arises
    - all_accessible 19 x 6 x 5  lag x cond x agegrp, includes unexposed, recovered, nil, mild, sick, severe
    - numcontacts 19 x 4 x 5  lag x cond x agegrp, includes nil, mild, sick, severe = the infectious

    There is a lot of setup before we get to business here.
=#
    # variables from simulation environment
    numcontacts =       env.numcontacts
    simple_accessible = env.simple_accessible
    lag_contacts =      env.lag_contacts
    touch_factors =     env.touch_factors
    numtouched =        env.numtouched

    # map to access maps conditions to the rows of simple_accessible and touch_factors
    map2access = (unexposed= 1, infectious=-1, recovered= 2, dead=-1, nil= 3, mild=  4, sick= 5, severe= 6)

    lag_contacts[:] = sum(numcontacts,dims=(2,3))[:,:,1] # (19, ) contacts by lag after sum by cond, agegrp
    totaccessible = sum(simple_accessible)

    # sum accessible to unexposed, recovered, infectious by agegrps (don't use nil, mild, sick, severe conditions)
            # unexposed, recovered, sum of(nil, mild, sick, severe) = infectious  6,5 with last 3 rows all zeros
    simple_accessible[:] = [simple_accessible[1:2,:]; sum(simple_accessible[3:6,:],dims=1); zeros(Int,3,5)];  # (6, 5)
    # s_a_pct is dist. of accessible by agegrp and uexposed, recovered, infectious (15,)
    s_a_pct = round.(reshape(simple_accessible[1:3,:] ./ totaccessible, 15), digits=5) # % for each cell
    if !isapprox(sum(s_a_pct), 1.0, atol=1e-8) 
        s_a_pct = s_a_pct ./ sum(s_a_pct) # normalize so sums to 1.0
    end

    # we only use this to makes sure we don't touch more people than there are
    unexp_by_contact_lag = zeros(Int, 19,5) # lags, agegrps
    lagpct = zeros(19)
    lagpct[:] = lag_contacts ./ (sum(lag_contacts) + 1e-8)
    if sum(lagpct) == 0.0
        lagpct[:] = fill(1.0/19.0, 19)
    end
    @assert isapprox(sum(lagpct), 1.0, atol=1e-4) "pct must sum to 1.0; got $(sum(lagpct))"
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
    numtouched[:] = zeros(Int, 19, length(agegrps)) # (19,5)
    # loop over numcontacts lag vector
    for lag in lags
        lc = lag_contacts[lag]
        x = rand(dcat, lc) # (15,) probabistically distribute contacts for a lag across accessible by unexposed|recovered|infectious, agegrp
        peeps = reshape([count(x .== i) for i in 1:15], 3,5)[1,:]  # (5,) distribute across all 3 groups, but only take unexposed
        for a in agegrps # probabilistically see who of the accessible is significantly touched
            cnt = binomial_one_sample(peeps[a], touch_factors[map2access.unexposed, a])
            numtouched[lag, a] = clamp(cnt, 0, unexp_by_contact_lag[lag, a])
        end
    end

    # return numtouched  # (19,5)
end


function how_many_infected(all_unexposed; env=env)
    # transmissibility by agegrp of recipient
    #=
        - multiply the transmissibility of the spreader times the transmissibility of the contacts
                by lag for spreaders and by agegrp for the contacts
        - use the infection factor in a binomial sample:  was the contact "successful" in causing infection?
        - we'll test to be sure we don't exceed the unexposed and reduce touches to 80% of unexposed by agegrp
    =#

    # touched_by_lag_age (19,5)     all_unexposed (5,)

    # variables from env
    touched_by_lag_age = env.numtouched

    # setup risk table
    env.riskmx[:] = send_risk_by_recv_risk(env.send_risk_by_lag, env.recv_risk_by_age)  # (19,5)

    newinfected = zeros(Int, length(agegrps))  # (5,)
    for age in agegrps
        for lag in lags
            newsick = binomial_one_sample(touched_by_lag_age[lag, age], env.riskmx[lag, age])
            newsick = clamp(newsick, 0, floor(Int,.8 * all_unexposed[age]))
            newinfected[age] += newsick
        end
    end

    @debug "\n newly infected: $newinfected  \n"

    return newinfected  # (num of agegrps, ) (only condition is nil, assumed lag = 1 => first day infected)
end


function send_risk_by_recv_risk(send_risk, recv_risk)
    repeat(recv_risk',outer=19) .* send_risk
end


function cleanup_spread_cases()
    for k in [:oldcf, :casecf, :old_comp, :new_comp, :simple_accessible,
        :spreaders, :comply_infected, :noncomply_infected]
        delete!(spread_stash, k)
    end
end