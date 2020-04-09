##########################
# isolation.jl
##########################



"""
For each simulation day entry in queue isolatedq:
- Remove people from open environment in array openmx.
- Add people to the isolated environment in array isolatedmx.
"""
function isolate_move!(;opendat=openmx, isodat=isolatedmx) # use only keyword args
    while !isempty(isolatedq)
        g = dequeue!(isolatedq)
        cnt = clamp(g.cnt, 0, grab(g.cond, g.agegrp, g.lag, g.locale, dat=opendat))
        cnt < g.cnt && (@warn "You tried to isolate more people than were in the category: proceeding with existing people.")
        minus!(cnt,g.cond, g.agegrp, g.lag, g.locale, dat=opendat)
        plus!(cnt, g.cond, g.agegrp, g.lag, g.locale, dat=isodat)
    end
end  # this works


function unisolate_move!(;opendat=openmx, isodat=isolatedmx)
    while !isempty(isolatedq)
        g = dequeue!(isolatedq)
        @assert g.cnt <= 0 "Unisolate assumes that cnt is negative: found positive"
        cnt = -g.cnt  # switch it to positve so the operations make sense
        cnt = clamp(cnt, 0, grab(g.cond, g.agegrp, g.lag, g.locale, dat=isodat))
        cnt < -g.cnt && (@warn "You tried to unisolate more people than were in the category: proceeding with existing people.")
        minus!(cnt,g.cond, g.agegrp, g.lag, g.locale, dat=isodat)
        plus!(cnt, g.cond, g.agegrp, g.lag, g.locale, dat=opendat)
    end
end  # this works

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
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, lag, locale, dat=opendat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, 0, available)  # limit to max
    cnt < scnt && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
end


function isolate!(num::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    @assert num > 0 "num must be greater than zero"
    available = grab(cond, agegrp, lag, locale, dat=opendat)  # max
    cnt = clamp(num, 0, available)  # limit to max
    cnt < num && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _isolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing
end  # this one works


function _isolate!(cnt::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    minus!(cnt, cond, agegrp, lag, locale, dat=opendat)  # move out 
    plus!(cnt, cond, agegrp, lag, locale, dat=isodat)  # move in
    cnt !== 0 && queuestats(cnt, locale, isolatestat)
    return nothing  # this one works!
end


function unisolate!(pct::Float64,cond,agegrp,lag,locale; opendat = openmx, isodat=isolatedmx)
    @assert 0.0 <= pct <= 1.0 "pct must be between 0.0 and 1.0"
    available = grab(cond, agegrp, lag, locale, dat=isodat)  # max
    scnt = binomial_one_sample(available, pct)  # sample
    cnt = clamp(scnt, 0, available)  # limit to max
    cnt < scnt && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing  # this one works!
end


function unisolate!(num::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    @assert num > 0 "num must be greater than zero"
    available = grab(cond, agegrp, lag, locale, dat=isodat)  # max
    cnt = clamp(num, 0, available)  # limit to max
    cnt < num && (@warn "You tried to isolate more people than were in the category: proceeding with available.")
    _unisolate!(cnt, cond, agegrp, lag, locale; opendat=opendat, isodat=isodat)
    return nothing
end  # this one works


function _unisolate!(cnt::Int64, cond, agegrp, lag, locale; opendat=openmx, isodat=isolatedmx)
    minus!(cnt, cond, agegrp, lag, locale, dat=isodat)
    plus!(cnt, cond, agegrp, lag, locale, dat=opendat)
    cnt !== 0 && queuestats(-cnt, locale, isolatestat)  # this doesn't work for unisolate
    return nothing
end