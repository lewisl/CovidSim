####################################################################
# cases.jl
#       pass these cases to run_a_sim in kwarg runcases as a list
#       cases = [CovidSim.<funcname>, CovidSim.<funcname>, ...]  
#       then run_a_sim(geofilename, n_days, locales; runcases=cases)
####################################################################


"""
Generate seeding cases.
inputs: day, cnt, lag, cond, agegrp
Two of the inputs may refer to multiple items and must match in number of items.

Returns a function that can be used in runcases input to run_a_sim.
"""
function seed_case_gen(day, cnt, lag, cond, agegrp) # these args go into the returned seed! case
    function scase(locale; opendat=openmx, isodat=isolatedmx, env=env)  # args must match runcases loop in run_a_sim
        seed!(day, cnt, lag, cond, agegrp, locale, dat=opendat)
    end
end

# some generated seed! cases-->these are global (in code)
# seed_6_12 = seed_case_gen(8, [0,6,6,0,0], 5, nil, agegrps)
# seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 5, nil, agegrps)


# some isolation cases
function isolate_case_1(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day] == 15
        isolate!(.25,[unexposed, nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.70,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    elseif ctr[:day] == 23
        isolate!(.50,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.70,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_1(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 120
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(1.0,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function isolate_case_2(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day] == 15
        isolate!(.40,[unexposed, nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    elseif ctr[:day] == 23
        isolate!(.60,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_2(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 69
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(1.0,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_2b(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 84
        unisolate!(.6,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(.6,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end


function isolate_case_3(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day] == 40
        isolate!(.40,[unexposed, nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    elseif ctr[:day] == 50
        isolate!(.60,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_3(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 80
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(1.0,[mild,sick,severe],agegrps,1:laglim,locale; opendat=opendat, isodat=isodat)
    end
end

# define a cases as mycase=Spreadcase(15,cf_array,tf_array,compliance_array_or_float)
# pass these cases to run_a_sim in kwarg spreadcases as a list--they'll be run in function spread!
# scases = [mycase, case_2, ...]  then run_a_sim(geofilename, n_days, locales; spreadcases=scases)
# cases above can be combined with these passing in both runcases and spreadcases
struct Spreadcase
    day::Int
    cf::Array{Float64,2}  # (4,5)
    tf::Array{Float64,2}  # (6,5)
    compliance::Union{Float64, Array{Float64,2}}
end

function spread_case_setter(cases=[]; env=env)
    for case in cases
        c = case(env=env)
        if c.day == ctr[:day]
            # before the case starts--ignore it
            # after the case--it's already in effect--nothing to change
            if iszero(c.compliance)  # signal to shutdown cases and restore defaults
                # restore defaults for spread!  
                default_env = initialize_sim_env()
                env.sd_compliance = default_env.sd_compliance
                env.contact_factors = default_env.contact_factors
                env.touch_factors = default_env.touch_factors
                delete!(spread_stash, :case_cf)
                delete!(spread_stash, :case_tf)
            else
            # set contact_factors, touch_factors to the case values
                # copy c.cf to spread_stash
                # copy c.cf to env.contact_factors->active contact_factors  => when running spreadsteps
                if !haskey(spread_stash, :default_cf)  # should only ever happen once for an entire simulation
                    spread_stash[:default_cf] = copy(env.contact_factors)
                end
                if !haskey(spread_stash, :default_tf)
                    spread_stash[:default_tf] = copy(env.touch_factors)
                end

                spread_stash[:case_cf] = copy(c.cf)  # shouldn't need copy, but it's safer
                spread_stash[:case_tf] = copy(c.tf)  #             "

            # set the compliance note: compliance is the same for spreaders and accessible
                    # it varies by agegrp and condition if desired
                # check all compliance values in [0.0, 1.0]
                @assert c.compliance .>= 0.0 "compliance values must be positive"
                @assert c.compliance .<= 1.0 "compliance values must be in [0.0,1.0]"
                env.sd_compliance .= copy(c.compliance)  # TODO do we need to copy? takes 2x time
            end # if for current day case
        end  # if test for today
    end # case for loop
end  # function case_setter


function spread_case_runner(density_factor, all_unexposed; env=env)
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
            spread_stash[:comply_spreaders] = copy(env.spreaders)
        end
    end  # for loop
    # total values for comply + nocomply
    newinfected = newinfected[1] .+ newinfected[2]
end


function sd_gen(;start=45, comply=.7, cf=(.2, 1.6), tf=(.18,.7))
    function sd_mod(;env=env)
        sd_mod = Spreadcase(start,
                        round.(shifter(env.contact_factors, cf...),digits=2),
                        round.(shifter(env.touch_factors, tf...),digits=2),
                        comply)
    end
end

# copy beyond the comment and run in the REPL, use as input
#
# mod_45 = sd_gen()  # with defaults
# mod_90 = sd_gen(start=90,cf=(.2,1.5), tf=(.18,.6),comply=.85)
# str_45 = sd_gen(start=45, comply=.90, cf=(.2,1.0), tf=(.18,.3))
# str_55 = sd_gen(start=55, comply=.95, cf=(.2,1.0), tf=(.18,.3))
# zer = sd_gen(start=90, comply=0.0)
