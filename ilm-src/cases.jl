####################################################################
# cases.jl
#       pass these cases to run_a_sim in kwarg runcases as a list
#       cases = [CovidSim.<funcname>, CovidSim.<funcname>, ...]  
#       then run_a_sim(geofilename, n_days, locales; runcases=cases)
####################################################################

####################################################################
# 
# - define a cases as mycase=Spreadcase(15,cf_array,tf_array,compliance_array_or_float)
# - pass these cases to run_a_sim in kwarg spreadcases as a list--they'll be run in function spread!
#        scases = [mycase, case_2, ...]  then run_a_sim(geofilename, n_days, locales; spreadcases=scases)
# - cases above can be combined with these passing in both runcases and spreadcases
# 
####################################################################

####################################################################
# seeding cases
####################################################################

"""
Generate seeding cases.
inputs: day, cnt, sickday, cond, agegrp
Two of the inputs may refer to multiple items and must match in number of items.

Returns a function that can be used in runcases input to run_a_sim.
"""
function seed_case_gen(day, cnt, sickday, cond, agegrp) # these args go into the returned seed! case
    function scase(locale, opendat, spreaddict)  # args must match runcases loop in run_a_sim
        seed!(day, cnt, sickday, cond, agegrp, locale, opendat)
    end
end


"""
    seed!(day, cnt, sickday, conds, agegrps, locale, dat)

This is the action function that implements a seeding case.
"""
function seed!(day, cnt, sickday, conds, agegrps, locale, dat)
    @assert length(sickday) == 1 "input only one sickday value"
    # @warn "Seeding is for testing and may result in case counts out of balance"
    if day == day_ctr[:day]
        println("*** seed day $(day_ctr[:day]) locale $locale....")
        for loc in locale
            for cond in conds
                @assert (cond in [nil, mild, sick, severe]) "Seed cases must have conditions of nil, mild, sick, or severe" 
                make_sick!(dat[loc]; cnt=cnt, fromage=agegrps, tocond=nil, tosickday=sickday)
            end
        end
    end
end


# some generated seed! cases-->these are global (in code)
# seed_6_12 = seed_case_gen(8, [0,6,6,0,0], 5, nil, agegrps)
# seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 5, nil, agegrps)


####################################################################
# isolation cases
####################################################################




####################################################################
# spread cases
####################################################################

        # all in spread.jl




# copy beyond the comment and run in the REPL, use as input
#
# mod_45 = sd_gen()  # with defaults
# mod_90 = sd_gen(start=90,cf=(.2,1.5), tf=(.18,.6),comply=.85)
# str_45 = sd_gen(start=45, comply=.90, cf=(.2,1.0), tf=(.18,.3))
# str_55 = sd_gen(start=55, comply=.95, cf=(.2,1.0), tf=(.18,.3))
# zer = sd_gen(start=90, comply=0.0)


