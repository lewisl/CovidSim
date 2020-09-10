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
function spread!(dat, locale::Int, env, density_factor::Float64 = 1.0)  # no spreadcases
    riskmx = env.riskmx
    locdat = dat[locale]
        
    spread_idx = findall(locdat[:, cpop_status] .== infectious) # index of all spreaders
    n_spreaders = size(spread_idx, 1);
    contactable_idx = findall(locdat[:, cpop_status] .!= dead) # index of all potential contacts
    n_contactable = size(contactable_idx, 1)

    n_spreaders, n_contacts, n_touched, n_newly_infected =  _spread!(locdat, spread_idx, 
                contactable_idx, env.contact_factors, env.touch_factors, env.riskmx, 
                density_factor, env.shape)

    push!(spreadq,
            (day=ctr[:day], locale=locale, spreaders=n_spreaders, contacts=n_contacts,
                touched=n_touched, infected=n_newly_infected)
            )
    
    return n_spreaders, n_contacts, n_touched, n_newly_infected
end


function spread!(dat, locale::Int, spreadcases::Array{Spreadcase, 1}, env, density_factor::Float64 = 1.0)
    riskmx = env.riskmx
    locdat = dat[locale]

    case = Spreadcase(1, (0,0), (0,0), 0.0)
    for i in spreadcases
        if i.day == ctr[:day] # there is a case that starts today!
            case = i
            if case.comply == 0.0  # use defaults for everyone
                spread_stash[:do_case] = false  # run spread! without cases
            else    # use case inputs 
                spread_stash[:do_case] = true
                spread_stash[:comply] = case.comply
                spread_stash[:cf] = shifter(env.contact_factors, case.cf...)
                spread_stash[:tf] = shifter(env.touch_factors, case.tf...)
            end
            break  # found a case that applies today--skip any later cases
        end
    end # if we go through loop w/o finding a case today, then either:
        # nothing is assigned in spread_stash OR nothing changes in spread_stash


    if !(get(spread_stash, :do_case, false)) # we've never had a case or we shut down the previous case
        # use defaults with no cases
        n_spreaders, n_contacts, n_touched, n_newly_infected = spread!(dat, locale, env, density_factor) 

    elseif spread_stash[:comply] == 1.0 # we have a case that applies to all with case parameters
        spread_idx = findall(locdat[:, cpop_status] .== infectious)
        n_spreaders = size(spread_idx, 1);
        contactable_idx = findall(locdat[:, cpop_status] .!= dead)
        n_contactable = size(contactable_idx, 1)

        n_spreaders, n_contacts, n_touched, n_newly_infected = _spread!(locdat, spread_idx, 
                        contactable_idx, spread_stash[:cf], spread_stash[:tf], 
                        env.riskmx, density_factor, env.shape)

    else # split the population into comply and nocomply for 0.0 < comply < 1.0: 
        spread_idx = findall(locdat[:, cpop_status] .== infectious)
        n_spreaders = size(spread_idx, 1)
        n_spreaders_comply = round(Int, spread_stash[:comply] * n_spreaders)

        contactable_idx = findall(locdat[:, cpop_status] .!= dead)
        n_contactable = size(spread_idx, 1)
        n_contactable_comply = round(Int, spread_stash[:comply] * n_contactable)

        n_newly_infected = 0

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

            n_spreaders, n_contacts, n_touched, n_newly_infected = .+((n_spreaders, n_contacts, n_touched, n_newly_infected),
                        _spread!(locdat, pass_spread_idx, pass_contactable_idx, pass_cf, pass_tf, env.riskmx, density_factor, env.shape))
        end
        
    end

    push!(spreadq,
            (day=ctr[:day], locale=locale, spreaders=n_spreaders, contacts=n_contacts,
                touched=n_touched, infected=n_newly_infected)
            )

    return n_spreaders, n_contacts, n_touched, n_newly_infected

end


function _spread!(locdat, spread_idx, contactable_idx, contact_factors, touch_factors, riskmx, density_factor, shape)
    # how many spreaders
    n_spreaders = size(spread_idx,1)

    # how many contacts?
    spreaders_to_contacts = zeros(Int, size(spread_idx,1),2) # second column for lag of the spreader
    for i in 1:size(spreaders_to_contacts, 1)  # for each spreader
        scale = density_factor * contact_factors[locdat[spread_idx[i], cpop_cond]-4, locdat[spread_idx[i], cpop_agegrp]]
        spreaders_to_contacts[i, 1] = round(Int,rand(Gamma(shape, scale))) # cnt of contacts for 1 spreader
        spreaders_to_contacts[i, 2] = locdat[spread_idx[i], cpop_lag]   # lag of this spreader 
    end
    n_contacts = sum(spreaders_to_contacts[:,1])
    n_contactable = size(contactable_idx, 1)

    # assign the contacts 
    n_target_contacts = min(n_contacts, n_contactable)
    contact_people = sample(contactable_idx, n_target_contacts, replace=false) # specific people contacted

    # which contacts are consequential touches? which touched get infected?
    n_touched = 0
    n_newly_infected = 0

    stop = 0
    for (nc, lag) in eachrow(spreaders_to_contacts)  # nc=numContacts, lag=lag of spreader
        start = stop + 1; stop = stop + nc

        stop = stop > n_target_contacts ? n_target_contacts : stop  

        for person in contact_people[start:stop]
            # person's characteristics
            status = locdat[person, cpop_status]  # TODO below crap needs to be fixed
            agegrp = locdat[person, cpop_agegrp]
            characteristic =  status in [1,3] ? [1,0,2][status] : locdat[person, cpop_cond]-2 # max(0,ilmat[person, cpop_cond]-2
            @debug characteristic < 1 && error("bad characteristic value")

            # touch outcome
            touched = rand(Binomial(1, touch_factors[characteristic, agegrp]))
            n_touched += touched

            # infection outcome
            if touched == 1 && characteristic == unexposed
                prob = riskmx[lag, agegrp]
                newly_infected = rand(Binomial(1, prob))
                if newly_infected == 1
                    locdat[person, cpop_cond] = nil # nil === asymptomatic or pre-symptomatic
                    locdat[person, cpop_status] = infectious
                end
                n_newly_infected += newly_infected
            end
        end
    end

    # println("ratio of infected to spreaders: $(n_newly_infected/n_spreaders)")

    return n_spreaders, n_contacts, n_touched, n_newly_infected
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
    r0mx = data_dict(1; lags=laglim, conds=length(conditions), agegrps=n_agegrps)  # single locale
    locale = 1
    population = convert(T_int[], 2_000_000)
    setup_unexposed!(r0mx, population, locale)

    # setup data
    all_unexposed = grab(unexposed, agegrps, 1, locale, r0mx)  # (5, ) agegrp for lag 1
    track_infected = zeros(T_int[], 5)
    track_contacts = zeros(T_int[], laglim, 4, 5)
    track_touched = zeros(T_int[], laglim, 6, 5)

    r0env.all_accessible[:] = grab([unexposed,recovered, nil, mild, sick, severe], agegrps, lags, locale, r0mx)  #   laglim x 6 x 5  lag x cond by agegrp
    r0env.simple_accessible[:] = sum(r0env.all_accessible, dims=1)[1,:,:] # sum all the lags result (6,5)
    if !isempty(compliance)
        r0env.simple_accessible[:] = round.(T_int[], compliance .* r0env.simple_accessible)
    end

    if sa_pct[1] != 1.0
        sa_pct = [sa_pct[1],sa_pct[2],sa_pct[3], fill(sa_pct[3]./4.0, 3)...]
        res = [r0env.simple_accessible[1,:] .* i for i in sa_pct]
        sanew = zeros(T_int[], 6, 5)
        @inbounds for i in 1:6
           sanew[i,:] .= round.(Int,res[i])
        end
        r0env.simple_accessible[:] = round.(T_int[], sanew)
    end

    age_relative = round.(T_int[], age_dist ./ minimum(age_dist))
    r0env.spreaders[:] = ones(T_int[], laglim, 4, agegrps)
    @inbounds for i in 1:5
        r0env.spreaders[:,:,i] .= age_relative[i]
    end
    if !isempty(dt)
        r0env.spreaders[2:laglim, :, :] .= T_int[](0)
        r0env.spreaders .*= T_int[](20)
        tot_spreaders = sum(r0env.spreaders)
    else
        r0env.spreaders[1,:,:] .= T_int[](0);
        tot_spreaders = round.(T_int[], sum(r0env.spreaders) / (laglim - 1))
    end

    input!(r0env.spreaders,infectious_cases,agegrps,lags,locale,r0mx)

    # parameters that drive r0
    !isempty(cf) && (r0env.contact_factors[:] = deepcopy(cf))
    !isempty(tf) && (r0env.touch_factors[:] = deepcopy(tf))
    isempty(shift_contact)  || (r0env.contact_factors[:] =shifter(r0env.contact_factors, shift_contact...))
    isempty(shift_touch) || (r0env.touch_factors[:] = shifter(r0env.touch_factors, shift_touch...))

    stopat = !isempty(dt) ? laglim : 1

    for i = 1:stopat
        disp && println("test day = $i, spreaders = $(sum(r0env.spreaders))")

        track_contacts .+= how_many_contacts!(r0env, density_factor)
        track_touched .+= how_many_touched!(r0env)
        track_infected .+= how_many_infected(all_unexposed, r0env)

        if !isempty(dt)  # optionally transition
            transition!(dt, decpoints, locale, r0mx)
            r0env.spreaders[:] = grab(infectious_cases,agegrps,lags,locale, r0mx)
        end
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
