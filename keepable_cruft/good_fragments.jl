using DelimitedFiles


mutable struct Row_indices
    start::Int64
    age_brackets::Int64
    sickday_brackets::Int64
    idx::Union{Array{UnitRange{Int64}, 1}, UnitRange{Int64}}
end
    
    # alternate creator method to build array of indices
    Row_indices(start, ages, sickdays) = Row_indices(
            start,
            ages,
            sickdays,
            create_row_indices(start, ages, sickdays)
        )


function build_indices(ages, sickdays)

    # access as scalars: row indices of population stats for each locale (column)
    id = 1
    size_cat = 2
    popsize =  3

    # population conditions
    unexposed = Row_indices(4,1,1,4:8) # active or cumulative
    exposed = Row_indices(9, ages,sickdays)  # aging over 18 days; means not_infectious
    infectious = Row_indices(104, ages, sickdays) # aging over 18 days
    isolated = Row_indices(sickdaylim9, ages, sickdays) # aging over 18 days
    recovered = Row_indices(674, 1, 1, 674:678) # active or cumulative

    # population outcomes
    nil = Row_indices(294, ages, sickdays) # isn't this the same as exposed and not infectious?
    mild = Row_indices(389, ages, sickdays)
    sick = Row_indices(484, ages, sickdays)
    severe = Row_indices(579, ages, sickdays)
    dead = Row_indices(679,1,1,679:683) # active or cumulative

    travel = Row_indices(684,1,1,684:688)

    return Dict("id" => id, "size_cat" => size_cat, "popsize" => popsize, "unexposed" => unexposed, 
            "exposed" => exposed, "infectious" => infectious, "isolated" => isolated, "recovered" => recovered, 
            "nil" => nil, "mild" => mild, "sick" => sick, "severe" => severe, "dead" => dead, "travel" => travel)
end

function create_row_indices(start, age_brackets, sickday_brackets)
    # simulation data indices: rows are features,  columns are locales
    ret = UnitRange{Int64}[]
    f = start
    for i = 1:age_brackets
        stop = f + sickday_brackets-1
        push!(ret, f:stop)
        f = stop + 1
    end
    return ret
end

function create_row_indices(obj::Row_indices)  
    ret = UnitRange{Int64}[]
    f = obj.start
    for i = 1:obj.age_brackets
        stop = f + obj.sickday_brackets - 1
        push!(ret, f:stop)
        f = stop + 1
    end
    return ret
end

function setup(symdata_filename, age_brackets, sickday_brackets)
    symdata = readdlm(symdata_filename, ',')
    pop =  build_indices(5, sickdaylim)  
    println("sick indices: ", pop["sick"].idx)
    println("locale 1 popsize: ", symdata[pop["popsize"],1])
    println("locale 2 unexposed age bracket 1: ", symdata[pop["unexposed"].idx[1][1],2])
    return (symdata, pop)

end


function daystep()  # one day of simulation

    # some people travel in; distribution of them are infectious


    # evolve residents across all sickdays
            # start from last sickday, by age group:
                    # severe distribute to sick, severe or die remain at sickday 18
                    # sick distribute to severe, mild remain at sickday 18
                    # mild distribute to sick, recovered (no nil at this point) at sickday 18
                    # nil distribute to recovered 
            # for earlier sickdays:
                    # severe distribute to sick, severe or die => advance one sickday
                    # sick distribute to sick, severe, mild => advance one sickday
                    # mild distribute to sick, recovered => advance one sickday
                    # nil distribute to nil, mild, recovered => advance one sickday

            # go dead => remove from previous condition; remove from resident, remove from totalpop:  make sure we can count # dead today (or more than yesterday)
            # go recovered => remove from previous condition; remove from exposed? remove from infectious        


    # infectious ones contact distribution of residents across age&condition&sickdays


    # of those contacted, a distribution become exposed & nil at sickday 0


    # some residents travel out
    

    # gather stats

end

####################################################################################
#   convenience functions for reading and inputting population statistics
####################################################################################

# single age, single sickday, one or more locales
# example: grab("exposed", 1, 1, 1:3)
function grab(item::String, age::Int, sickday::Int, locale; dat=symdata, popidx=pop)
    return dat[popidx[item].idx[age][sickday], locale]
end

# single age, multiple sickdays, one or more locales
# example: grab("exposed", 1, 1:10, 1:3)  for sickday 0 to 9 days
function grab(item::String, age::Int, sickday::UnitRange{Int64}, locale; dat=symdata, popidx=pop)
    return dat[popidx[item].idx[age][sickday[1]]:popidx[item].idx[age][sickday[end]], locale]
end 

# 1 or more values to single age, single sickday, one more more locales
function input!(val, item::String, age::Int, sickday::Int, locale; dat=symdata, popidx=pop)
    dat[popidx[item].idx[age][sickday], locale] .= val
end

# 1 or more values to single age, multiple sickdays, one more more locales
function input!(val, item::String, age::Int, sickday::UnitRange{Int64}, locale; dat=symdata, popidx=pop)
    dat[popidx[item].idx[age][sickday[1]]:popidx[item].idx[age][sickday[end]], locale] .= val
end


###############################################################################
function update_infectious!(locale, dat=openmx) # by single locale
    for agegrp in agegrps
        tot = total!([nil, mild, sick, severe],agegrp,:,locale) # sum across cases and sickdays per locale and agegroup
        input!(tot, infectious, agegrp, 1, locale) # update the infectious total for the locale and agegroup
    end
end


    # a not ok way to estimate R0:  need each subsequent day of infected for the originating group; this is only 1 day
    # starters = sum(spreaders)
    # newlyinfected = sum(byage)
    # onestep_r0 = newlyinfected / starters
    # @debug "Spreaders $starters to newly infected: $newlyinfected for r0: $onestep_r0"



#################################################################################
# not such a good way to update an ILM pop matrix
#################################################################################


function update!(dat; cnt=cnt, tests=[[cpop_status, ==, 1], [cpop_agegrp, ==, 3]], 
                    todo=[[cpop_cond, 5], [cpop_sickday, 1]])

    tcnt = length(tests)

    tests = [[tst[1], tst[2], tst[3]] for tst in tests]

    @show typeof(tests)
    @show typeof(todo)

    truthtests = falses(tcnt)  # allocate once
    @show typeof(truthtests)

    did = 0
    n_rows = size(dat, 1)
    if cnt == 0
        cnt = n_rows
    end

    @inbounds for i = 1:n_rows  
        if did < cnt
            @simd for j in 1:tcnt
                truthtests[j] = tests[j][2](dat[i, tests[j][1]], tests[j][3])
            end

            if all(truthtests) 

                @inbounds @simd for act in todo
                    dat[i, act[1]] = act[2]
                end

                did += 1
            end
        else
            break
        end
    end

end

# this is an alternative for filtering the ilm pop matrix
    prefilt = [actions.cmps[i].(dat[:, actions.tests[i][1]], actions.tests[i][2]) for i in 1:length(actions.tests)]
    filt[:] = .&(
                 prefilt...
                )
