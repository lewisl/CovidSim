#############
# spread.jl
#############

"""
Stash for temporary values changed during simulation cases
- to change just once and then get the originals back
- it is the users responsibility to empty the stash
- there may be (will be) side effects if you don't empty the stash between simulations
"""
const spread_stash = Dict{Symbol, Any}()


struct Spreadcase
    day::Int
    cf::Tuple{Float64,Float64}  # (4,5)
    tf::Tuple{Float64,Float64}  # (6,5)
    comply::Float64
end


function sd_gen(;start=45, comply=.7, cf=(.2, 1.6), tf=(.18,.7))
    Spreadcase(start, cf, tf, comply)
end



"""
How far do the infectious people spread the virus to
previously unexposed people, by agegrp?  For a single locale...

The alternative method processes spreadcases for social distancing. If comply percentage
is not 1.0 or 0.0, the population is split into complying and non-complying.
"""
function spread!(dat, locale::Int, spreadcases, env, density_factor::Float64 = 1.0)

    locdat = locale == 0 ? dat : dat[locale]

    # spread_idx = findall(locdat[:, cpop_status] .== infectious) # index of all spreaders
    spread_idx = findall(locdat.status .== infectious) # index of all spreaders
    n_spreaders = size(spread_idx, 1);
    # contactable_idx = findall(locdat[:, cpop_status] .!= dead) # index of all potential contacts
    contactable_idx = findall(locdat.status .!= dead) # index of all potential contacts

    # println("finding the contactable took ", @elapsed findall(locdat.status .!= dead))

    n_contactable = size(contactable_idx, 1)

    do_case = get(spread_stash, :do_case, false) # we've never had a case or we shut down the previous case

    if  isempty(spreadcases) && !do_case  # no cases left and do_case is false 

        n_contacts, n_touched, n_newly_infected = ( 
            _spread!(locdat, spread_idx, 
                contactable_idx, env.contact_factors, env.touch_factors, env.riskmx, env.shape,
                density_factor)
            )

    else  # a case may start today OR we have an active case

        n_contacts, n_touched, n_newly_infected = (
            spread_cases(locdat, spreadcases, spread_idx, contactable_idx, env, density_factor)
            )
    end

    push!(spreadq,
        (day=ctr[:day], locale=locale, spreaders=n_spreaders, contacts=n_contacts,
            touched=n_touched, infected=n_newly_infected)
        )

    return n_spreaders, n_contacts, n_touched, n_newly_infected
end


function spread_cases(locdat, spreadcases, spread_idx, contactable_idx, env, density_factor)

    for (i,case) in enumerate(spreadcases)
        if case.day == ctr[:day] # there is a case that starts today!
            # case = i
            if case.comply == 0.0  # cancel spreadcase and restore defaults
                spread_stash[:do_case] = false  # run spread! without cases
            elseif 0.0 < case.comply <= 1.0    # use case inputs 
                spread_stash[:do_case] = true
                spread_stash[:comply] = case.comply
                spread_stash[:cf] = shifter(env.contact_factors, case.cf...)
                spread_stash[:tf] = shifter(env.touch_factors, case.tf...)
            else
                @error "spreadcase comply value must be >= 0.0 and <= 1.0: got $(case.comply)"
            end
            # pop this case:  we don't need to look at it again
            deleteat!(spreadcases, i)
            break   # we can only have one active case at a time; new cases replaces prior case; do one per day
        end
    end # if we go through loop w/o finding a case today, then nothing changes in spread_stash

    do_case = get(spread_stash, :do_case, false)

    if do_case
        n_spreaders = size(spread_idx, 1);
        n_contactable = size(contactable_idx, 1)            

        if spread_stash[:comply] == 1.0 # we have a case that applies to everyone using the case parameters

            n_contacts, n_touched, n_newly_infected = _spread!(locdat, spread_idx, 
                            contactable_idx, spread_stash[:cf], spread_stash[:tf], 
                            env.riskmx, env.shape, density_factor)
                

        elseif 0.0 < spread_stash[:comply] < 1.0  # split the population into comply and nocomply for 0.0 < comply < 1.0: 

            n_spreaders_comply = round(Int, spread_stash[:comply] * n_spreaders)
            n_contactable_comply = round(Int, spread_stash[:comply] * n_contactable)

            n_contacts = n_touched = n_newly_infected = 0

            for pass in [:comply, :nocomply]
                if pass == :comply  # use case input factors
                    pass_cf = spread_stash[:cf]
                    pass_tf = spread_stash[:tf]

                    nsp = n_spreaders_comply  
                    ncon = n_contactable_comply            
                elseif pass == :nocomply # use default factors
                    pass_cf = env.contact_factors
                    pass_tf = env.touch_factors

                    nsp = n_spreaders - n_spreaders_comply  
                    ncon = n_contactable - n_contactable_comply            
                end

                pass_spread_idx = sample(spread_idx, nsp, replace=false)
                pass_contactable_idx = sample(contactable_idx, ncon, replace=false)

                n_contacts, n_touched, n_newly_infected = .+((n_contacts, n_touched, n_newly_infected),
                            _spread!(locdat, pass_spread_idx, pass_contactable_idx, pass_cf, pass_tf, env.riskmx, env.shape, density_factor))
            end
        end
    else  # a case was zero'ed out today--this is same as running spread! -> that might happen next day if no more cnt_accessible
        n_contacts, n_touched, n_newly_infected = _spread!(locdat, spread_idx, contactable_idx, 
                        env.contact_factors, env.touch_factors, env.riskmx, env.shape,
                        density_factor)
    end

    return n_contacts, n_touched, n_newly_infected
    
end


function _spread!(locdat, spread_idx, contactable_idx, contact_factors, touch_factors, riskmx, shape, density_factor)

    # how many contacts?

    spreaders_to_contacts = zeros(Int, size(spread_idx,1), 2) # second column for lag of the spreader

    @inbounds for i in 1:size(spread_idx, 1)  # for each spreader  # size(spreaders_to_contacts, 1)
        p = spread_idx[i]

        thiscond = locdat.cond[p] - 4  # map 5-8 to 1-4
        thisagegrp = locdat.agegrp[p]
        thislag = locdat.lag[p]

        scale = density_factor * contact_factors[thiscond, thisagegrp]

        spreaders_to_contacts[i, 1] = round(Int,rand(Gamma(shape, scale))) # cnt of contacts for 1 spreader
        spreaders_to_contacts[i, 2] =  thislag
    end

    # n_contacts = sum(spreaders_to_contacts.nc)     # sum(spreaders_to_contacts[:,1])   sum(spreaders_to_contacts.nc) 
    n_contacts = sum(spreaders_to_contacts[:,1])  

    n_contactable = size(contactable_idx, 1)


    # assign the contacts 
    n_target_contacts = min(n_contacts, n_contactable)

    contact_people = sample(contactable_idx, n_target_contacts, replace=false)


    # which contacts are consequential touches? which touched get infected?
    n_touched = 0
    n_newly_infected = 0

    stop = 0
    @inbounds for i in 1:size(spread_idx,1)  # nc=numContacts, lag=lag of spreader

        nc = spreaders_to_contacts[i,1]
        lag = spreaders_to_contacts[i,2]

        start = stop + 1; stop = stop + nc
        # stop = stop > n_target_contacts ? n_target_contacts : stop  
        # @assert stop <= n_target_contacts

        @inbounds @views for person in contact_people[start:stop]
            # person's characteristics
            status = locdat.status[person]  # TODO below crap needs to be fixed
            agegrp = locdat.agegrp[person]
            cond = locdat.cond[person]
            lookup = if status == unexposed
                        1  # row 1
                     elseif status == recovered
                        2
                     else
                        cond - 2  # 5:8 - 2 -> 3:6
                     end

            # touch outcome
            touched = rand(Binomial(1, touch_factors[lookup, agegrp]))
            n_touched += touched

            # infection outcome
            if (touched == 1) && (status == unexposed)    # (characteristic == unexposed)
                prob = riskmx[lag, agegrp]
                newly_infected = rand(Binomial(1, prob))
                if newly_infected == 1
                    locdat.cond[person] = nil # nil === asymptomatic or pre-symptomatic
                    locdat.status[person] = infectious
                    # lag remains zero because person was unexposed; transition! function updates lag
                end
                n_newly_infected += newly_infected
            end
        end
    end

    return n_contacts, n_touched, n_newly_infected
end


function send_risk_by_recv_risk(send_risk, recv_risk)
    recv_risk' .* send_risk  # (laglim, agegrps)
end


function cleanup_stash(stash)
    for k in keys(stash)
        delete!(stash, k)
    end
end


function r0_sim(age_dist, dat, locale::Int, dt_dict, env, density_factor, pop=1_000_000; scale=10)

    # setup a fake locale or use the current locale in the simulation
    @views if locale == 0 # simulate with fake locale with pop people
        # create population
        r0pop = pop_data(pop, age_dist=age_dist,intype=Int16,cols="all")

        # create spreaders: ages by age_dist, cond is nil
        age_relative = round.(T_int[], age_dist ./ minimum(age_dist))
        age_relative .*= scale # update with scale
        cnt_spreaders = sum(age_relative)
        for i in agegrps
            idx = findfirst(x->x==i, r0pop.agegrp)
            for j = 1:age_relative[i]
                r0pop.status[idx] = infectious
                r0pop.cond[idx] = nil
                r0pop.lag[idx] = 1
                idx += 1
            end
        end
    else  # simulate r0 at the current state of locale being simulated
        r0pop = deepcopy(dat[locale])
        ignore_idx = findall(r0pop.status .== infectious)
        r0pop.status[ignore_idx] .= recovered # can't catch what they already have; won't spread for calc of r0

        cnt_accessible = count(r0pop.status .== unexposed)
        age_relative = round.(T_int[], age_dist ./ minimum(age_dist)) # counts by agegrp
        scale = set_by_level(cnt_accessible)
        age_relative .*= scale # update with scale
        cnt_spreaders = sum(age_relative)

        for i in agegrps
            idx = findall((r0pop.agegrp .== i) .& (r0pop.status .== unexposed))
            for j = 1:age_relative[i]
                p = idx[j]
                r0pop.status[p] = infectious
                r0pop.cond[p] = nil
                r0pop.lag[p] = 1
            end
        end        
    end

    ret = [0,0,0,0] # n_spreaders, n_contacts, n_touched, n_newly_infected 
    for i = 1:laglim                                          
        ret[:] .+= spread!(r0pop, 0, [], env, density_factor)  # dat, locale, spreadcases, env, density_factor::Float64 = 1.0)

        transition!(r0pop, 0, dt_dict) 

        # eliminate the new spreaders so we only track the original spreaders
        newsick_idx = findall(r0pop.lag .== 1)

        r0pop.status[newsick_idx] .= unexposed
        r0pop.lag[newsick_idx] .= 0
    end

    r0 =  ret[4] / cnt_spreaders   # n_newly_infected / cnt_spreaders
    return r0
end


function set_by_level(x, levels=[[1, 300_000], [5, 500_000], [10, 10_000_000_000]])
    ret = 0
    for lvl in levels
        if x <= lvl[2]
            ret = lvl[1]
            break
        end
    end
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
  tf       1.1   1.2      1.3   1.4     1.5   1.6    1.7    1.8    1.9    2.0
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
