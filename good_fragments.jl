using DelimitedFiles


mutable struct Row_indices
    start::Int64
    age_brackets::Int64
    lag_brackets::Int64
    idx::Union{Array{UnitRange{Int64}, 1}, UnitRange{Int64}}
end
    
    # alternate creator method to build array of indices
    Row_indices(start, ages, lags) = Row_indices(
            start,
            ages,
            lags,
            create_row_indices(start, ages, lags)
        )


function build_indices(ages, lags)

    # access as scalars: row indices of population stats for each locale (column)
    id = 1
    size_cat = 2
    popsize =  3

    # population conditions
    unexposed = Row_indices(4,1,1,4:8) # active or cumulative
    exposed = Row_indices(9, ages,lags)  # aging over 18 days; means not_infectious
    infectious = Row_indices(104, ages, lags) # aging over 18 days
    isolated = Row_indices(199, ages, lags) # aging over 18 days
    recovered = Row_indices(674, 1, 1, 674:678) # active or cumulative

    # population outcomes
    nil = Row_indices(294, ages, lags) # isn't this the same as exposed and not infectious?
    mild = Row_indices(389, ages, lags)
    sick = Row_indices(484, ages, lags)
    severe = Row_indices(579, ages, lags)
    dead = Row_indices(679,1,1,679:683) # active or cumulative

    travel = Row_indices(684,1,1,684:688)

    return Dict("id" => id, "size_cat" => size_cat, "popsize" => popsize, "unexposed" => unexposed, 
            "exposed" => exposed, "infectious" => infectious, "isolated" => isolated, "recovered" => recovered, 
            "nil" => nil, "mild" => mild, "sick" => sick, "severe" => severe, "dead" => dead, "travel" => travel)
end

function create_row_indices(start, age_brackets, lag_brackets)
    # simulation data indices: rows are features,  columns are locales
    ret = UnitRange{Int64}[]
    f = start
    for i = 1:age_brackets
        stop = f + lag_brackets-1
        push!(ret, f:stop)
        f = stop + 1
    end
    return ret
end

function create_row_indices(obj::Row_indices)  
    ret = UnitRange{Int64}[]
    f = obj.start
    for i = 1:obj.age_brackets
        stop = f + obj.lag_brackets - 1
        push!(ret, f:stop)
        f = stop + 1
    end
    return ret
end

function setup(symdata_filename, age_brackets, lag_brackets)
    symdata = readdlm(symdata_filename, ',')
    pop =  build_indices(5, 19)  
    println("sick indices: ", pop["sick"].idx)
    println("locale 1 popsize: ", symdata[pop["popsize"],1])
    println("locale 2 unexposed age bracket 1: ", symdata[pop["unexposed"].idx[1][1],2])
    return (symdata, pop)

end


function daystep()  # one day of simulation

    # some people travel in; distribution of them are infectious


    # evolve residents across all lags
            # start from last lag, by age group:
                    # severe distribute to sick, severe or die remain at lag 18
                    # sick distribute to severe, mild remain at lag 18
                    # mild distribute to sick, recovered (no nil at this point) at lag 18
                    # nil distribute to recovered 
            # for earlier lags:
                    # severe distribute to sick, severe or die => advance one lag
                    # sick distribute to sick, severe, mild => advance one lag
                    # mild distribute to sick, recovered => advance one lag
                    # nil distribute to nil, mild, recovered => advance one lag

            # go dead => remove from previous condition; remove from resident, remove from totalpop:  make sure we can count # dead today (or more than yesterday)
            # go recovered => remove from previous condition; remove from exposed? remove from infectious        


    # infectious ones contact distribution of residents across age&condition&lags


    # of those contacted, a distribution become exposed & nil at lag 0


    # some residents travel out
    

    # gather stats

end

####################################################################################
#   convenience functions for reading and inputting population statistics
####################################################################################

# single age, single lag, one or more locales
# example: grab("exposed", 1, 1, 1:3)
function grab(item::String, age::Int, lag::Int, locale; dat=symdata, popidx=pop)
    return dat[popidx[item].idx[age][lag], locale]
end

# single age, multiple lags, one or more locales
# example: grab("exposed", 1, 1:10, 1:3)  for lag 0 to 9 days
function grab(item::String, age::Int, lag::UnitRange{Int64}, locale; dat=symdata, popidx=pop)
    return dat[popidx[item].idx[age][lag[1]]:popidx[item].idx[age][lag[end]], locale]
end 

# 1 or more values to single age, single lag, one more more locales
function input!(val, item::String, age::Int, lag::Int, locale; dat=symdata, popidx=pop)
    dat[popidx[item].idx[age][lag], locale] .= val
end

# 1 or more values to single age, multiple lags, one more more locales
function input!(val, item::String, age::Int, lag::UnitRange{Int64}, locale; dat=symdata, popidx=pop)
    dat[popidx[item].idx[age][lag[1]]:popidx[item].idx[age][lag[end]], locale] .= val
end


###############################################################################
function update_infectious!(locale, dat=openmx) # by single locale
    for agegrp in agegrps
        tot = total!([nil, mild, sick, severe],agegrp,:,locale) # sum across cases and lags per locale and agegroup
        input!(tot, infectious, agegrp, 1, locale) # update the infectious total for the locale and agegroup
    end
end


    # a not ok way to estimate R0:  need each subsequent day of infected for the originating group; this is only 1 day
    # starters = sum(spreaders)
    # newlyinfected = sum(byage)
    # onestep_r0 = newlyinfected / starters
    # @debug "Spreaders $starters to newly infected: $newlyinfected for r0: $onestep_r0"