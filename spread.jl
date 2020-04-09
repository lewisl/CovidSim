#############
# spread.jl
#############


"""
How far do the infectious people spread the virus to
previously unexposed people, by agegrp?  For a single locale...
"""
function spread!(locale; dat=openmx, env=env)


    # how many spreaders  TODO grab their condition.  Separate probs by condition
    env.spreaders[:] = grab(infectious_cases, agegrps, lags, locale, dat=dat) # 19 x 4 x 5 lag x cond x agegrp
    @debug "  all the spreaders $(sum(env.spreaders))"

    if sum(env.spreaders) == 0
        return
    end

    env.all_accessible[:] = grab([unexposed,recovered, nil, mild, sick, severe],agegrps,lags, locale, dat=dat)  #   19 x 6 x 5  lag x cond by agegrp
    all_unexposed = grab(unexposed, agegrps, 1, locale, dat=dat)  # (5, ) agegrp for lag 1

    # how many people are contacted based on characteristics of spreader
    env.numcontacts[:] = how_many_contacts(env.spreaders, contact_factors, env=env)
    @debug "  all the contacts $(sum(env.numcontacts))"

    # how many people are touched based on characteristics of recipient and potential contacts?
    numtouched = how_many_touched(env.numcontacts, env.all_accessible, env=env)   # cond x agegrp
    @debug "  all the touched $(sum(numtouched))"

    newinfected = how_many_infected(numtouched, all_unexposed, env=env)    # x agegrp (only condition is nil, assumed lag = 1)
    @debug "  all the newly infected $(sum(newinfected))"

    lag = 1
    # test if newinfected > unexposed
    for agegrp in agegrps
        if newinfected[agegrp] > grab(unexposed, agegrp, lag, locale, dat=dat)
            println("big problem: infected exceeds unexposed")
        end
    end

    # move the people from unexposed:agegrp to infectious:agegrp and nil
    plus!.(newinfected, infectious, agegrps, lag, locale, dat=dat)
    plus!.(newinfected, nil, agegrps, lag, locale, dat=dat)
    minus!.(newinfected, unexposed, agegrps, lag, locale, dat=dat)

    push!(bugq, (day=ctr[:day], locale=locale, spreaders = sum(env.spreaders), contacts = sum(env.numcontacts),
                    touched = sum(numtouched),
                    unexposed=sum(grab(unexposed, agegrps, lag, locale, dat=dat)),
                    infected=sum(newinfected)))
    # add to stats queue for today
    queuestats(sum(newinfected), locale, spreadstat) # sum(5 agegroups), nil is the default, single locale

    return
end


# TODO calculate density_factor in setup, per locale
# TODO fix logistic shift and scale
test_density = rand((5000:3_000_000),20)  # use US Census data
function minmax(x)
    x_max = maximum(x, dims=1)
    x_min = minimum(x, dims=1)
    minmax_density = (x .- x_min) ./ (x_max .- x_min .+ 1e-08)
end
scale_minmax(x, newmin, newmax) = x .* (newmax - newmin) .+ newmin


# contact factors for the spreaders
            # agegrp     1     2      3       4     5
const contact_factors = [1     2.5    2.5     1.5   1;  # nil
                         1     2.5    2.5     1.5   1;  # mild
                         0.7   1.0    1.0     0.7   0.5;  # sick
                         0.5   0.8    0.8     0.5  0.2  # severe
                        ]


"""
How many contacts do spreaders attempt to make?  This is based on the characteristics of the
spreaders.
"""
function how_many_contacts(spreaders, contact_factors, density_factor=1.3; scale=6, env=env)
    #=  This originally ignores the conditions of the touched--assumes they are all equally likely to be touched
        how_many_touched corrects this.
        We assume spreaders is small compared to all_accessible. At some point this might not be true:
        how_many_touched also handles this.
    =#
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
                dgamma = Gamma(1.2, density_factor * scale)  #shape, scale
                x = round.(Int,rand(dgamma,spcount))

                if isempty(x)
                else
                    env.numcontacts[lag, cond, agegrp] = sum(x)
                end
            end
        end
    end

    return env.numcontacts
end


"""
For potential contacts by spreaders reaching out, how many of the accessible (susceptible and NOT)
are actuallly "touched" by a spreader? This is loosely based on the characteristics of the
receiver.
"""
function how_many_touched(numcontacts, all_accessible; env=env)
#=
    - who is accessible: start with all to capture effect of "herd immunity", when it arises
    - all_accessible 19 x 6 x 5  lag x cond x agegrp, includes unexposed, recovered, nil, mild, sick, severe
    - numcontacts 19 x 4 x 5  lag x cond x agegrp, includes nil, mild, sick, severe = the infectious

    There is a lot of setup before we get to business here.
=#

    # parameters for accessibility of the accessible--not who gets sick, just who is touched!
    map2access = (unexposed= 1, infectious=-1, recovered= 2, dead=-1, nil= 3, mild=  4, sick= 5, severe= 6)
    low, low2, low3 = .6, .5, .3
    lowest, lowest2, lowest3, vlowest = .35, .28, .18, .15
    highest, high, high2 = .9, .8, .7

    access_table = zeros(count(x->x>0, map2access),length(agegrps))  #     6 x 5    cond x agegrp
    access_table[map2access.unexposed, :] .= [low, highest, high, low, lowest]  # only using this one for now
    access_table[map2access.recovered, :] .= [low, highest, high, low, lowest]  # row goes across agegrps
    access_table[map2access.nil, :] .= [low, highest, high, low, lowest]
    access_table[map2access.mild, :] .= [low, highest, high2, low2, lowest2]
    access_table[map2access.sick, :] .= [lowest2, lowest, lowest2, lowest3, lowest3]
    access_table[map2access.severe,:] .= vlowest
    # not varying access pct for lag of infectious states because we ignore increased viral load from repeat exposures

    env.lag_contacts[:] = sum(numcontacts,dims=(2,3))[:,:,1] # (19, ) contacts by lag after sum by cond, agegrp
    # totcontacts = sum(env.lag_contacts)
    totaccessible = sum(all_accessible)


    # simplify accessible to unexposed, recovered, infectious by agegrps
    env.simple_accessible[:] = sum(all_accessible, dims=1)[1,:,:] # sum all the lags result (6,5)
    # next: unexposed, recovered, sum of(nil, mild, sick, severe) = infectious
    env.simple_accessible[:] = [env.simple_accessible[1:2,:]; sum(env.simple_accessible[3:6,:],dims=1); zeros(Int,3,5)];  # (6, 5)

    sptime = @elapsed begin
        s_a_split_by_lag = zeros(Int, 19,5)

        pct = zeros(19)
        pct[:] = env.lag_contacts ./ (sum(env.lag_contacts) + 1e-8)

        if sum(pct) == 0.0
            pct[:] = fill(1.0/19.0, 19)
        end


        @assert isapprox(sum(pct), 1.0, atol=1e-4) "pct must sum to 1.0 $(sum(pct))"
        for i in agegrps, j in lags
              s_a_split_by_lag[j,i] = round(Int,env.simple_accessible[map2access.unexposed, i] * pct[j])
        end


        s_a_pct = round.(reshape(env.simple_accessible[1:3,:] ./ totaccessible, 15), digits=3) # % for each cell
        if !isapprox(sum(s_a_pct), 1.0, atol=1e-8)
            s_a_pct = s_a_pct ./ sum(s_a_pct) # normalize so sums to 1.0
        end
    end

    # now to business: who gets touched in unexposed by agegrp?

    #=
        - folks in each cell of numcontacts touch a sample of the accessible reduced by the accessibility factor
        - draw a categorical sample for each cell to distribute them across the contact categories
        - we consider the folks who get contacts in infectious and recovered, because that happens--reduces
           the number of unexposed who are touched
        - we throw out the folks in infectious and recovered and keep those in unexposed
        - we should care about not touching more than the number of accessible
    =#

    mapi = (unexposed= 1, infectious=3, recovered=2, dead=-1, nil= -1, mild= -1, sick= -1, severe= -1)

    dcat = Categorical(s_a_pct) # categorical distribution by accessible pct
    touched_by_age_cond = zeros(Int, 19, length(agegrps)) # (19,5)
    # loop over numcontacts lag vector
    for lag in lags
        lc = env.lag_contacts[lag]
        x = rand(dcat, lc) # probabistically distribute contacts for 1 lag across accessible by cond, agegrp

        peeps = reshape([count(x .== i) for i in 1:15], 3,5)[1,:]  # (5,) after distributing across accessible,
                # use only the first row for unexposed by agegrp

        for a in agegrps # probabilistically see who of the accessible is "willing" to be touched
            cnt = binomial_one_sample(peeps[a], access_table[map2access.unexposed, a])
            touched_by_age_cond[lag,a] = clamp(cnt, 0, s_a_split_by_lag[lag,a])
        end
    end

    return touched_by_age_cond  # (19,5)
end


function how_many_infected(touched_by_age_cond, all_unexposed; env=env)
    # transmissibility by agegrp of recipient
    #=
        - multiply the transmissibility of the spreader times the transmissibility of the contacts
                by lag for spreaders and by agegrp for the contacts
        - use the infection factor in a binomial sample:  was the contact "successful" in causing infection?
        - we'll test to be sure we don't exceed the unexposed and reduce touches to 80% of unexposed by agegrp
    =#

    # touched_by_age_cond (19,5)     all_unexposed (5,)

    # setup risk table
    env.riskmx[:] = send_risk_by_recv_risk(send_risk_by_lag, recv_risk_by_age)  # (19,5)

    newinfected = zeros(Int, length(agegrps))  # (5,)
    for age in agegrps
        for lag in lags
            newsick = binomial_one_sample(touched_by_age_cond[lag, age], env.riskmx[lag, age])
            newsick = clamp(newsick, 0, floor(Int,.8 * all_unexposed[age]))
            newinfected[age] += newsick
        end
    end

    @debug "\n newly infected: $newinfected  \n"

    return newinfected
end
