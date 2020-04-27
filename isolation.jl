##########################
# isolation.jl
##########################


"""
Place people who are isolating into the isolated queue: isolatedq

You can enter a percentage (as a fraction in [0.0, 1.0]), an array
of percentages, a number, or an array of numbers.

Use a dot after the function name to apply the same pct or number
input to several conditions, agegrps, lags, or locales.

Use a dot after the function name to apply an array: one or more
of agegrp, cond, or locale must have the same number of elements as the input.
"""
function isolate!(pct::Float64,cond,agegrp,lag,locale; opendat=openmx, isodat=isolatedmx)
    for c in cond
        for age in agegrp
            for l in lag
                isolate_by!(pct::Float64,c,age,l,locale; opendat=opendat, isodat=isodat)
            end
        end
    end
end

function isolate_by!(pct::Float64,cond,agegrp,lag,locale; opendat=openmx, isodat=isolatedmx)
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, lag, locale, dat=opendat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, 0, available)  # limit to max
    cnt < scnt && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
end


function isolate_by!(num::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    @assert num > 0 "num must be greater than zero"
    available = grab(cond, agegrp, lag, locale, dat=opendat)  # max
    cnt = clamp(num, 0, available)  # limit to max
    cnt < num && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing
end  # this one works


function _isolate!(cnt::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    minus!(cnt, cond, agegrp, lag, locale, dat=opendat)  # move out 
    update_infectious!(locale, dat=opendat)
    plus!(cnt, cond, agegrp, lag, locale, dat=isodat)  # move in
    update_infectious!(locale, dat=isodat)
    cnt !== 0 && queuestats(cnt=cnt, locale=locale, conds=cond, agegrp=agegrp, event=:isolate)
    return nothing  # this one works!
end


function unisolate!(pct::Float64,cond,agegrp,lag,locale; opendat=openmx, isodat=isolatedmx)
    for c in cond
        for age in agegrp
            for l in lag
                unisolate_by!(pct::Float64,c,age,l,locale; opendat=opendat, isodat=isodat)
            end
        end
    end
end


function unisolate_by!(pct::Float64,cond,agegrp,lag,locale; opendat = openmx, isodat=isolatedmx)
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, lag, locale, dat=isodat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, 0, available)  # limit to max
    cnt < scnt && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing  # this one works!
end


function unisolate_by!(num::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    @assert num > 0 "num must be greater than zero"
    available = grab(cond, agegrp, lag, locale, dat=isodat)  # max
    cnt = clamp(num, 0, available)  # limit to max
    cnt < num && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing
end  # this one works


function _unisolate!(cnt::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    minus!(cnt, cond, agegrp, lag, locale, dat=isodat)
    update_infectious!(locale, dat=isodat)
    plus!(cnt, cond, agegrp, lag, locale, dat=opendat)
    update_infectious!(locale, dat=opendat)
    cnt !== 0 && queuestats(cnt=-cnt, locale=locale, conds=cond, agegrp=agegrp, event=:isolate)
    return nothing
end