################################
# spread.jl for ilm model
#    social distancing cases
#    spreading the infection
################################


#######################################################################
# social distancing cases
#      struct to hold parameters for defining the case
#      implement the case: 
#           - set social distance compliance for each person
#           - define the contact_factors and touch_factors for the case
#######################################################################          
# mod_90 = sd_gen(start=90,cf=(.2,1.5), tf=(.18,.6),comply=.85)
# str_45 = sd_gen(start=45, comply=.90, cf=(.2,1.0), tf=(.18,.3))
# str_55 = sd_gen(start=55, comply=.95, cf=(.2,1.0), tf=(.18,.3))

Base.@kwdef struct Spreadcase                 # Base.@kwdef -> use keyword arguments in constructor
    day::Int
    cfdelta::Tuple{Float64,Float64}  
    tfdelta::Tuple{Float64,Float64}  
    comply::Float64             # compliance fraction
    cfcase::Dict{Int64, Dict{String, Float64}}
    tfcase::Dict{Int64, Dict{String, Float64}}
end

function sd_gen(;startday::Int, comply::Float64, cf::Tuple{Float64, Float64},
    tf::Tuple{Float64, Float64}, name::Union{String, Symbol}, agegrps=[])
    function scase(locale, dat, spreadparams, sdcases, ages)   # scase(locale, dat, spreadparams, sdcases)
        s_d_seed!(dat, sdcases, startday, comply, cf, tf, name, agegrps, locale, spreadparams, ages)
    end
end


@inline function s_d_seed!(dat, sdcases, startday, comply, cf, tf, name, agegrps, locale, spreadparams, ages)
    @assert 0.0 <= comply <= 1.0  "comply must be floating point in 0.0 to 1.0 inclusive"
    
    if startday == day_ctr[:day]
        name = Symbol(name)
        locdat = dat[locale]

        if comply == 0.0  # magic signal: if comply is zero turn off this case for the selected agegrps
            cancel_sd_case!(locdat, sdcases, name, agegrps, ages)
            return
        end

        # create the Spreadcase in sdcases
        sdcases[name] = Spreadcase(
                        day     = startday,   
                        cfdelta = cf,         
                        tfdelta = tf,         
                        comply  = comply,     
                        cfcase  = shifter(spreadparams.contact_factors, cf...),  
                        tfcase  = shifter(spreadparams.touch_factors, tf...)     
                        )

        # load the s_d_comply column of the population table
        # filter1 is everyone who is unexposed, recovered or sick: nil or mild
        filter1 = findall(((locdat.status .== 1) .| (locdat.status .== 3)) .| ((locdat.cond .== 5) .| (locdat.cond .== 6)))
        if (comply == 1.0)   # include everyone in filter1
            complyfilter = filter1
        else
            complyfilter = sample(filter1, round(Int, comply*length(filter1)), replace=false)
        end
        
        if isempty(agegrps)   # include all agegrps
            locdat.s_d_comply[complyfilter] .= name
        else
            byage_idx = intersect(complyfilter, union((ages[i] for i in agegrps)...))
            locdat.s_d_comply[byage_idx] .= name
        end
    end
end


function cancel_sd_case!(locdat, sdcases, name, agegrps, ages)
    # filter on who is in this case now
    incase_idx = findall(locdat.s_d_comply .== name)

    if isempty(agegrps)   # include all agegrps
        locdat.s_d_comply[incase_idx] .= :none
        delete!(sdcases, name)  # there is no one left...
    else  # only turn it off for some agegrps
        byage_idx = intersect(incase_idx, union((ages[i] for i in agegrps)...))
        locdat.s_d_comply[byage_idx] .= :none
    end    

end


@inline @inbounds @fastmath function numcontacts(density_factor, shape, agegrp, cond, contact_factors)::Int 
    scale = density_factor * contact_factors[agegrp][condnames[cond]]
    round(Int,rand(Gamma(shape, scale)))
end

@inline @inbounds @fastmath function numcontacts(density_factor, shape, agegrp, cond, acase::Spreadcase)::Int
    scale = density_factor * acase.cfcase[agegrp][condnames[cond]]  
    round(Int,rand(Gamma(shape, scale)))
end

@inline @inbounds @fastmath function istouched(agegrp, lookup, touch_factors)::Int
    rand(Binomial(1, touch_factors[agegrp][lookup]))    
end

@inline @inbounds @fastmath function istouched(agegrp, lookup, acase::Spreadcase)::Int
    rand(Binomial(1, acase.tfcase[agegrp][lookup]))  
end



"""
How far do the infectious people spread the virus to
previously unexposed people, by agegrp?  For a single locale...
"""
@inline function spread!(locdat, infect_idx, contactable_idx, sdcases, spreadparams, density_factor)

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
        spreadercond = locdat.cond[p]  
        spreaderagegrp = locdat.agegrp[p]
        spreadersickday = locdat.sickday[p]
        spreadersdcomply = locdat.s_d_comply[p]

        nc =    if spreadersdcomply == :none
                    numcontacts(density_factor, shape, spreaderagegrp, spreadercond, contact_factors) 
                else
                    numcontacts(density_factor, shape, spreaderagegrp, spreadercond, sdcases[spreadersdcomply])
                end
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
            contactcomply = locdat.s_d_comply[contact]
            contactlookup = if contactstatus == unexposed   # set lookup row in touch_factors
                        "unexposed"  # row 1
                     elseif contactstatus == recovered
                        "recovered"  # row 2
                     else
                        condnames[contactcond]  # text names of conds 5:8 - 2 -> rows 3:6
                     end

            if contactstatus == unexposed  # only condition that can get infected   TODO: handle reinfection of recovered
                touched =   if contactcomply == :none
                                istouched(contactagegrp, contactlookup, touch_factors)
                            else
                                istouched(contactagegrp, contactlookup, sdcases[contactcomply])
                            end

                # infection outcome
                if (touched == 1) && (contactstatus == unexposed)    # TODO some recovered people will become susceptible again
                    prob = riskmx[spreadersickday, contactagegrp]            # TODO also vaccinated people will have partially unsusceptible
                    newly_infected = rand(Binomial(1, prob))
                    if newly_infected == 1
                        locdat.cond[contact] = nil # nil === asymptomatic or pre-symptomatic
                        locdat.status[contact] = infectious
                        locdat.sickday[contact] = 1
                    end
                end  # if (touched ...)
            end  # if contactstatus
        end  # for contact in sample(...)
    end  # for p in infect_idx

    return # n_contacts, n_touched, n_newly_infected
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
