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
    function scase(locale, opendat, isodat, testdat, env)  # args must match runcases loop in run_a_sim
        seed!(day, cnt, sickday, cond, agegrp, locale, opendat)
    end
end


# some generated seed! cases-->these are global (in code)
# seed_6_12 = seed_case_gen(8, [0,6,6,0,0], 5, nil, agegrps)
# seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 5, nil, agegrps)


####################################################################
# isolation cases
####################################################################

function isolate_case_1(locale; opendat, isodat, testdat, env)
    if day_ctr[:day] == 15
        isolate!(.25,[unexposed, nil],agegrps,1,locale, opendat, isodat)
        isolate!(.70,[mild,sick, severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    elseif day_ctr[:day] == 23
        isolate!(.50,[unexposed,nil],agegrps,1,locale, opendat, isodat)
        isolate!(.70,[mild,sick, severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    end
end

function unisolate_case_1(locale; opendat, isodat, testdat, env)
    if day_ctr[:day]  == 120
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale, opendat, isodat)
        unisolate!(1.0,[mild,sick, severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    end
end

function isolate_case_2(locale; opendat, isodat, testdat, env)
    if day_ctr[:day] == 15
        isolate!(.40,[unexposed, nil],agegrps,1,locale, opendat, isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    elseif day_ctr[:day] == 23
        isolate!(.60,[unexposed,nil],agegrps,1,locale, opendat, isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    end
end

function unisolate_case_2(locale; opendat, isodat, testdat, env)
    if day_ctr[:day]  == 69
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale, opendat, isodat)
        unisolate!(1.0,[mild,sick, severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    end
end

function unisolate_case_2b(locale; opendat, isodat, testdat, env)
    if day_ctr[:day]  == 84
        unisolate!(.6,[unexposed,nil],agegrps,1,locale, opendat, isodat)
        unisolate!(.6,[mild,sick, severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    end
end


function isolate_case_3(locale; opendat, isodat, testdat, env)
    if day_ctr[:day] == 40
        isolate!(.40,[unexposed, nil],agegrps,1,locale, opendat, isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    elseif day_ctr[:day] == 50
        isolate!(.60,[unexposed,nil],agegrps,1,locale, opendat, isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    end
end

function unisolate_case_3(locale; opendat, isodat, testdat, env)
    if day_ctr[:day]  == 80
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale, opendat, isodat)
        unisolate!(1.0,[mild,sick,severe],agegrps,1:sickdaylim,locale, opendat, isodat)
    end
end


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


