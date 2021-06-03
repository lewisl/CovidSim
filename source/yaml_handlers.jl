# functions to handle YAML file custom tags

function array_enum(obj)
    # for tag !array_enum
    status_cond_map.(obj)
end
    
function status_cond_map(x)
    # for tags !status  !cond
    val = if x == 0
                condition(0)
          elseif x <= Int(typemax(status))
                status(x)
          elseif x <= Int(typemax(condition))
                condition(x)
          else
                @assert false "Invalid status or condition $x"
          end
end

function agegroup(x)
    # for tag !agegrp
    agegrp(x)
end