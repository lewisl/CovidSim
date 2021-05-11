################################
# spread.jl for ilm model
################################

"""
Stash for temporary values changed during simulation cases
- to change just once and then get the originals back
- it is the users responsibility to empty the stash
- there may be (will be) side effects if you don't empty the stash between simulations
"""
const sim_stash = Dict{Symbol, Any}()


struct Spreadcase
    day::Int
    cf::Tuple{Float64,Float64}  # (4,5)  # contact factors
    tf::Tuple{Float64,Float64}  # (6,5)  # touch factors
    comply::Float64             # compliance percentage
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
@inline function spread!(locdat, infect_idx, contactable_idx, spreadcases, 
    spreadparams, density_factor)

    contact_factors = spreadparams.contact_factors
    touch_factors = spreadparams.touch_factors
    riskmx = spreadparams.riskmx
    shape = spreadparams.shape

    # assign contacts, do touches, do new infections
    @inbounds @fastmath for p in infect_idx      # p is the person who is the spreader

        # spreader's characteristics

        # if in_quarantine(locdat, p, 0.0)  # person is in quarantine--can't spread
        #     return  # or continue
        # end
        spreadercond = locdat.cond[p]  # map 5-8 to 1-4
        spreaderagegrp = locdat.agegrp[p]
        spreadersickday = locdat.sickday[p]

        scale = density_factor * contact_factors[spreaderagegrp][condnames[spreadercond]]

        nc = round(Int,rand(Gamma(shape, scale))) # number of contacts for 1 spreader
        # n_contacts += nc
                                # TODO we could keep track of contacts for contact tracing
        @inbounds @fastmath for contact in sample(contactable_idx, nc, replace=false) # people can get contacted more than once
            # contacts's characteristics
            # if in_quarantine(locdat, contact, 0.0)  # person is in quarantine--can't spread
            #     return
            # end
            contactstatus = locdat.status[contact]  
            contactagegrp = locdat.agegrp[contact]
            contactcond = locdat.cond[contact]
            contactlookup = if contactstatus == unexposed   # set lookup row in touch_factors
                        "unexposed"  # row 1
                     elseif contactstatus == recovered
                        "recovered"  # row 2
                     else
                        condnames[contactcond]  # text names of conds 5:8 - 2 -> rows 3:6
                     end

            if contactstatus == unexposed  # only condition that can get infected   TODO: handle reinfection of recovered
                # touch outcome
                touched = rand(Binomial(1, touch_factors[contactagegrp][contactlookup]))
                # n_touched += touched

                # infection outcome
                if (touched == 1) && (contactstatus == unexposed)    # TODO some recovered people will become susceptible again
                    prob = riskmx[spreadersickday, contactagegrp]            # TODO also vaccinated people will have partially unsusceptible
                    newly_infected = rand(Binomial(1, prob))
                    if newly_infected == 1
                        locdat.cond[contact] = nil # nil === asymptomatic or pre-symptomatic
                        locdat.status[contact] = infectious
                        locdat.sickday[contact] = 1
                        # NOT ANY MORE sickday remains zero because person was unexposed; transition! function updates sickday
                    end
                    # n_newly_infected += newly_infected
                end
            end
        end
    end

    return # n_contacts, n_touched, n_newly_infected
end


@inline function old_spread!(locdat, infect_idx, contactable_idx, spreadcases, spreadparams, 
                         density_factor::Float64 = 1.0)

    # locdat = locale == 0 ? dat : dat[locale]

    # filttime = @elapsed
    # begin # indices of all spreaders; indices of all who could be contacted
        # infect_idx = findall(locdat.status .== infectious)    # optfindall(==(infectious), locdat.status, 0.5)   # must use parens around 1st comparison for operator precedence
        # contactable_idx = findall(locdat.status .!= dead)     # optfindall(!=(dead), locdat.status, 1)
        # shuffle!(infect_idx); shuffle!(contactable_idx)

        # n_spreaders = size(infect_idx, 1)
        # n_contactable = size(contactable_idx, 1)
    # end
    # @show filttime

    do_case = get(sim_stash, :do_case, false) # we've never had a case or we shut down the previous case

    if isempty(spreadcases) && !do_case  # no cases left and do_case is false 

        _spread!(locdat, infect_idx, contactable_idx,
                spreadparams.contact_factors, spreadparams.touch_factors, 
                spreadparams.riskmx, spreadparams.shape, density_factor)  # n_contacts, n_touched, n_newly_infected = 
            
    else  # a case may start today OR we have an active case

        spread_cases(locdat, spreadcases,           # n_contacts, n_touched, n_newly_infected = 
                infect_idx, contactable_idx, spreadparams, density_factor)
            
    end

    # push!(spreadq,   (day=day_ctr[:day], locale=locale, spreaders=n_spreaders, contacts=n_contacts,
                    #   touched=n_touched, infected=n_newly_infected))

    return #n_spreaders, n_contacts, n_touched, n_newly_infected
end


function spread_cases(locdat, spreadcases, infect_idx, contactable_idx, spreadparams, density_factor)

    for (i,case) in enumerate(spreadcases)
        if case.day == day_ctr[:day] # there is a case that starts today!
            if case.comply == 0.0  # cancel spreadcase and restore defaults
                sim_stash[:do_case] = false  # run spread! without cases
            elseif 0.0 < case.comply <= 1.0    # use case inputs 
                sim_stash[:do_case] = true
                sim_stash[:comply] = case.comply
                sim_stash[:cf] = shifter(spreadparams.contact_factors, case.cf...)
                sim_stash[:tf] = shifter(spreadparams.touch_factors, case.tf...)
            else
                @error "spreadcase comply value must be >= 0.0 and <= 1.0: got $(case.comply)"
            end
            
            deleteat!(spreadcases, i) # pop this case:  we don't need to look at it again
            break   # we can only have one active case at a time; new cases replaces prior case; do one per day
        end
    end # if we go through loop w/o finding a case today, then nothing changes in sim_stash

    do_case = get(sim_stash, :do_case, false)

    if do_case
        n_spreaders = size(infect_idx, 1);
        n_contactable = size(contactable_idx, 1)            

        if sim_stash[:comply] == 1.0 # we have a case that applies to everyone using the case parameters

             _spread!(locdat, infect_idx,                                      # n_contacts, n_touched, n_newly_infected =
                            contactable_idx, sim_stash[:cf], sim_stash[:tf], 
                            spreadparams.riskmx, spreadparams.shape, density_factor)
                

        elseif 0.0 < sim_stash[:comply] < 1.0  # split the population into comply and nocomply for 0.0 < comply < 1.0: 
            @views begin
                n_spreaders_comply = round(Int, sim_stash[:comply] * n_spreaders)
                n_contactable_comply = round(Int, sim_stash[:comply] * n_contactable)
                infect_idx_split = Dict()
                contactable_idx_split = Dict()

                shuffle!(infect_idx); shuffle!(contactable_idx)  # enable random split of each
                infect_idx_split[:comply] = infect_idx[1:n_spreaders_comply]
                infect_idx_split[:nocomply] = infect_idx[(n_spreaders_comply + 1):n_spreaders]

                contactable_idx_split[:comply] = contactable_idx[1:n_contactable_comply]
                contactable_idx_split[:nocomply] = contactable_idx[(n_contactable_comply + 1):n_contactable]


                # n_contacts = n_touched = n_newly_infected = 0

                for pass in [:comply, :nocomply]
                    if pass == :comply  # use case input factors
                        pass_cf = sim_stash[:cf]
                        pass_tf = sim_stash[:tf]

                        nsp = n_spreaders_comply  
                        ncon = n_contactable_comply            
                    elseif pass == :nocomply # use default factors
                        pass_cf = spreadparams.contact_factors
                        pass_tf = spreadparams.touch_factors

                        nsp = n_spreaders - n_spreaders_comply  
                        ncon = n_contactable - n_contactable_comply            
                    end

                    pass_infect_idx = infect_idx_split[pass]
                    pass_contactable_idx = contactable_idx_split[pass] 

                    # n_contacts, n_touched, n_newly_infected = .+((n_contacts, n_touched, n_newly_infected),
                    _spread!(locdat, pass_infect_idx, pass_contactable_idx, pass_cf, pass_tf, spreadparams.riskmx, spreadparams.shape, density_factor)   #)
                end
            end # begin block
        end
    else  # a case was zero'ed out today--this is same as running spread! -> that might happen next day if no more cnt_accessible
        _spread!(locdat, infect_idx, contactable_idx,                                    # n_contacts, n_touched, n_newly_infected =
                spreadparams.contact_factors, spreadparams.touch_factors, spreadparams.riskmx, 
                spreadparams.shape, density_factor)
    end

    return #n_contacts, n_touched, n_newly_infected
end





function send_risk_by_recv_risk(send_risk, recv_risk)
    recv_risk' .* send_risk  # (sickdaylim, agegrps)
end


function cleanup_stash(stash)
    for k in keys(stash)
        delete!(stash, k)
    end
end


function r0_sim(age_dist, dat, locale::Int, dt_dict, spreadparams, density_factor, pop=1_000_000; scale=10)

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
                r0pop.sickday[idx] = 1
                idx += 1
            end
        end
    else  # simulate r0 at the current state of locale being simulated
        r0pop = deepcopy(dat[locale])
        ignore_idx = optfindall(==(infectious), r0pop.status, 0.5)
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
                r0pop.sickday[p] = 1
            end
        end        
    end

    ret = [0,0,0,0] # n_spreaders, n_contacts, n_touched, n_newly_infected 
    for i = 1:sickdaylim                                          
        ret[:] .+= spread!(r0pop, 0, [], spreadparams, density_factor)  

        transition!(r0pop, 0, dt_dict) 

        # eliminate the new spreaders so we only track the original spreaders
        newsick_idx = findall(r0pop.sickday .== 1)

        r0pop.status[newsick_idx] .= unexposed
        r0pop.sickday[newsick_idx] .= 0
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


function r0_table(n=6, cfstart = 0.9, tfstart = 0.3; spreadparams=spreadparams, dt=dt)
    tbl = zeros(n+1,n+1)
    cfiter = [cfstart + (i-1) * .1 for i=1:n]
    tfiter = [tfstart + (i-1) * 0.05 for i=1:n]
    for (j,cf) in enumerate(cfiter)
        for (i,tf) = enumerate(tfiter)
            tbl[i+1,j+1] = r0_sim(spreadparams=spreadparams, dt=dt, decpoints=decpoints, shift_contact=(0.2,cf), shift_touch=(.18,tf)).r0
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
