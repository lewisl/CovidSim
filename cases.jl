
# pass these cases to run_a_sim in kwarg runcases as a list
# cases = [CovidSim.<funcname>, CovidSim.<funcname>, ...]  then run_a_sim(geofilename, n_days, locales; runcases=cases)
function isolate_case_1(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day] == 15
        isolate!(.25,[unexposed, nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.70,[mild,sick, severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat)
    elseif ctr[:day] == 23
        isolate!(.50,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.70,[mild,sick, severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_1(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 69
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(1.0,[mild,sick, severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat) 
    end
end

function isolate_case_2(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day] == 15
        isolate!(.40,[unexposed, nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat)
    elseif ctr[:day] == 23
        isolate!(.60,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_2(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 69
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(1.0,[mild,sick, severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat) 
    end
end

function unisolate_case_2b(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 84
        unisolate!(.6,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(.6,[mild,sick, severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat) 
    end
end


function isolate_case_3(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day] == 40
        isolate!(.40,[unexposed, nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat)
    elseif ctr[:day] == 50
        isolate!(.60,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        isolate!(.75,[mild,sick, severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat)
    end
end

function unisolate_case_3(locale; opendat=openmx, isodat=isolatedmx, env=env)
    if ctr[:day]  == 80
        unisolate!(1.0,[unexposed,nil],agegrps,1,locale; opendat=opendat, isodat=isodat)
        unisolate!(1.0,[mild,sick,severe],agegrps,1:19,locale; opendat=opendat, isodat=isodat) 
    end
end


# pass these cases to run_a_sim in kwarg spreadcases as a list--they'll be run in function spread!
# scases = [CovidSim.<funcname>, CovidSim.<funcname>, ...]  then run_a_sim(geofilename, n_days, locales; spreadcases=scases)
# cases above can be combined with these passing in both runcases and spreadcases

function social_distance_case_1(density_factor, all_unexposed; env=env)
    startday = 40; endday = 150

    if ctr[:day] >= startday  && ctr[:day] <= endday

        # copy the old values
        if !haskey(spread_stash, :old_comp)
            spread_stash[:old_comp] = copy(env.sd_compliance) # copy isolates stashed array from changes
        end
        if !haskey(spread_stash, :oldcf)
            spread_stash[:oldcf] = copy(env.contact_factors) # copy isolates stashed array from changes
        end

        # do the math once; retrieve thereafter
        env.sd_compliance[:] = if !haskey(spread_stash, :new_comp)
                                    spread_stash[:new_comp] = copy(spread_stash[:old_comp]) .* 0.95
                                else
                                    copy(spread_stash[:new_comp])  # env.sd_compliance is not changed during the case
                                end  # 75 to 85 ns

        # change the contact factors, just once: 
        env.contact_factors[:] = if !haskey(spread_stash, :casecf)  # squash the contact_factors, stash, return
                                    oldmin = minimum(env.contact_factors)
                                    oldmax = maximum(env.contact_factors)
                                    newmin = oldmin * .9
                                    newmax = oldmax * .6 > newmin ? oldmax * .6 : newmin * 1.1
                                    spread_stash[:casecf] = round.(shifter(env.contact_factors, newmin, newmax),digits=2)
                                    copy(spread_stash[:casecf])  # copy is needed isolate the stash
                                    # print("new cf  ");println(spread_stash[:casecf]) 
                                    # print("  old cf "); println(spread_stash[:oldcf])
                                else  
                                    copy(spread_stash[:casecf]) # get the case contact_factors Do we need to copy each day?
                                end

        # println(" compliant sd_compliance ",env.sd_compliance[3:6,:])

        # get spreaders for the beginning of a new spread! day--before splitting into compliant and noncompliant
        spread_stash[:all_spr] = copy(env.spreaders)  # stash today's spreaders--isolated from env

        # run spread! for the complying
        env.spreaders[:]= round.(Int,permutedims(permutedims(copy(spread_stash[:all_spr]),[2,3,1]) .* env.sd_compliance[3:6,:], [3,1,2]))
        how_many_contacts!(density_factor, env=env)
        how_many_touched!(env=env)
        spread_stash[:comply_infected] = how_many_infected(all_unexposed, env=env)

        # println(" compliant spreaders  $(sum(env.spreaders))   ")
        # println(" compliant contacts  $(sum(env.numcontacts))   ")
        # println(" compliant touched   $(sum(env.numtouched))  ")
        # println(" compliant infected   $(sum(spread_stash[:comply_infected]))       ")

        # run spread! for the non-complying
        # set spreaders to non-complying
        # println(" noncompliant sd_compliance ",env.sd_compliance[3:6,:])  -- need the copy?
        env.spreaders[:] = round.(Int, permutedims(permutedims(copy(spread_stash[:all_spr]),[2,3,1]) .* (1.0 .- env.sd_compliance[3:6,:]), [3,1,2]))

        # println(" split of spreaders: total: ", sum(spread_stash[:all_spr]), 
        #         " comply ", sum(round.(Int,permutedims(permutedims(spread_stash[:all_spr],[2,3,1]) .* env.sd_compliance[3:6,:], [3,1,2]))),
        #         " noncomply ", sum(round.(Int, permutedims(permutedims(copy(spread_stash[:all_spr]),[2,3,1]) .* (1.0 .- env.sd_compliance[3:6,:]), [3,1,2])))         )

        # retrieve old contact_factors for non-complying
        env.contact_factors[:] = copy(spread_stash[:oldcf]) # copy  isolates the stashed copy
        # run the spread
        how_many_contacts!(density_factor, env=env)
        how_many_touched!(env=env)
        spread_stash[:noncomply_infected] = how_many_infected(all_unexposed, env=env)

        # println(" noncompliant spreaders  $(sum(env.spreaders))   ")
        # println(" noncompliant contacts  $(sum(env.numcontacts))   ")
        # println(" noncompliant touched   $(sum(env.numtouched))  ")
        # println(" noncompliant infected   $(sum(spread_stash[:noncomply_infected]))       ")

        newinfected = spread_stash[:noncomply_infected] .+ spread_stash[:comply_infected]
    else # before and after the case days
        if ctr[:day] == endday + 1
            # reset contact_factors and sd_compliance back to "normal"  this will stick
            env.contact_factors[:] = copy(spread_stash[:oldcf])
            env.sd_compliance[:] = copy(spread_stash[:old_comp])
        end
        how_many_contacts!(density_factor, env=env)
        how_many_touched!(env=env)
        newinfected = how_many_infected(all_unexposed, env=env)
    end

    return newinfected
end



# function unsocial_distance_case_1(density_factor; env=env)
#     if ctr[:day] == 75
#         oldmin = minimum(env.contact_factors)
#         oldmax = maximum(env.contact_factors)
#         newmin = oldmin / .9
#         newmax = oldmax / .85 > newmin ? oldmax / .85 : newmin / 1.1
#         env.contact_factors[:] = shifter(env.contact_factors, newmin, newmax)
#     end
# end

function social_distance_case_2(density_factor; env=env)
    if ctr[:day] == 40
        oldmin = minimum(env.touch_factors)
        oldmax = maximum(env.touch_factors)
        newmin = oldmin * .9
        newmax = oldmax * .8 > newmin ? oldmax * .8 : newmin * 1.4
        env.touch_factors[:] = shifter(env.touch_factors, newmin, newmax)
    end
end

function unsocial_distance_case_2(density_factor; env=env)
    if ctr[:day] == 75
        oldmin = minimum(env.touch_factors)
        oldmax = maximum(env.touch_factors)
        newmin = oldmin / .9
        newmax = oldmax / .85 > newmin ? oldmax / .85 : newmin / 1.4
        env.touch_factors[:] = shifter(env.touch_factors, newmin, newmax)
    end
end