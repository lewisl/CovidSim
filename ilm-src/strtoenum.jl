function strtoenum(enumgrp::Type{<:Enum{T}}, str::String) where {T<:Integer}
    dict = Base.Enums.namemap(enumgrp)
    found = -1
    srch = Symbol(str)
    for (k,v) in dict
        if v == srch
            found = k
            break
        end
    end
    if found != -1
        return enumgrp(found)
    end
end


function strtoenum2(enumgrp::Type{<:Enum{T}}, str::String) where {T<:Integer}
    srch = Symbol(str)
    found = -1
    for val in instances(enumgrp)
        if srch == Symbol(val)
            found = val
            break
        end
    end
    if found != -1
        return found
    end
end


function strtoenum2(enumgrp::Type{<:Enum{T}}, sym::Symbol) where {T<:Integer}
    srch = sym
    found = -1
    for val in instances(enumgrp)
        if srch == Symbol(val)
            found = val
            break
        end
    end
    if found != -1
        return found
    end
end

# here's how to do it!
#=
julia> inst = instances(condition)
(nil, mild, sick, severe)

julia> syms = Symbol.(inst)
(:nil, :mild, :sick, :severe)

julia> condsym = Dict(zip(inst, syms))
Dict{condition, Symbol} with 4 entries:
  sick   => :sick
  mild   => :mild
  severe => :severe
  nil    => :nil

julia> symcond = Dict(zip(syms, inst))
Dict{Symbol, condition} with 4 entries:
  :sick   => sick
  :mild   => mild
  :nil    => nil
  :severe => severe

julia> @btime $symcond[:severe]
  5.831 ns (0 allocations: 0 bytes)
severe::condition = 8  
=#