
#############################################################
# dec_tree.jl
# decision tree for transition
#############################################################


#=
{                                           # the enclosing dict
    { 1:                                    # each agegrp             
        { 5:                                # transition day for an agegrp                      
            { 5:                            # current condition on this transition day          
                { probs: [0.4, 0.5, 0.1]    # probs and outcomes for the next transition    
                    outcomes: [5, 6, 7]  
                }         
            }
        }    
    } 
}    
=#

function setup_dt(dtfilename)
    trees = YAML.load_file(dtfilename)

    newdict = (
    Dict(agegrp(k1) =>         
        Dict(k2 =>             
            Dict(condition(k3) => 
                Dict(Symbol(k4) => v4 for (k4, v4) in v3) 
                                            for (k3, v3) in v2)
                                                for (k2, v2) in v1)
                                                    for (k1,v1) in trees)
          )



    # Convert values in :outcomes to Enum condition or status
    for (k1,v1) in newdict
        for (k2,v2) in v1
            for (k3,v3) in v2
                for (k4, v4) in v3
                    if k4 == :outcomes
                        outs = [i > Int(dead) ? condition(i) : status(i)  for i in v4]
                        newdict[k1][k2][k3][k4] = outs
                    end
                end
            end
        end
    end

    return newdict
end


function display_tree(tree)
    for agegrp in keys(tree)
        agetree = tree[agegrp]
        println("agegrp: ", agegrp, " =>")
        for sickday in keys(agetree)
            sickdaytree = agetree[sickday]
            println("    sickday: ", sickday, " =>")
            for fromcond in keys(sickdaytree)
                condtree = sickdaytree[fromcond]
                println("        fromcond: ", fromcond, " =>")
                print("            probs: => ")
                println(condtree[:probs])
                #
                print("            outcomes: => ")
                println(condtree[:outcomes])
                #
                # println("            branches: =>")
                # for branch in keys(condtree["branches"])
                #     println("                ", condtree["branches"][branch])   
                # end
            end  # for fromcond
        end  # for sickday
    end   # for agegrp     
end


function sanitycheck(dectree)
    for i in agegrps
        seqs = CovidSim_ilm.walksequence(dectree[i])
        probs, allpr = CovidSim_ilm.verifyprobs(seqs)
        println("for agegroup ", i)
        for p in pairs(probs)
            println("    ",p)
        end
        println("    Prob total: ",allpr)
    end
end


"""
Find all sequences of conditions by transition date and current condition through to new conditions
for a single agegrp.
"""
function walksequence(dt_by_age)
    # find the top nodes
    dt_by_age = sort(dt_by_age)
    breakdays = collect(keys(dt_by_age))
    k1 = first(breakdays)
    todo = [] # array of node sequences 
    done = [] # ditto

    # gather the outcomes at the first breakday for the starting conditions
    # no transition has happened yet: these are initial conditions
    for fromcond in keys(dt_by_age[k1])
        for i in 1:length(dt_by_age[k1][fromcond][:outcomes])
            outcome = dt_by_age[k1][fromcond][:outcomes][i]
            prob = dt_by_age[k1][fromcond][:probs][i]
            push!(todo, [(sickday=k1, fromcond=fromcond, tocond=outcome, prob=prob)])
        end
    end


    # build sequences from top to terminal states: recovered or dead
    while !isempty(todo)
        seq = popfirst!(todo)  
        lastnode = seq[end]
        # @show lastnode
        breakday, fromcond, tocond = lastnode
        # tocond = condition(tocond)
        nxtidx = findfirst(isequal(breakday), breakdays) + 1
        for brk in breakdays[nxtidx:end]
            if tocond in keys(dt_by_age[brk])   # keys are the fromcond at the next break day so previous tocond == current fromcond
                for i in 1:length(dt_by_age[brk][tocond][:outcomes])
                    outcome = dt_by_age[brk][tocond][:outcomes][i]
                    prob = dt_by_age[brk][tocond][:probs][i]
                    if (outcome == dead) | (outcome == recovered)  # terminal node reached--no more nodes to add
                        push!(done, vcat(seq, (sickday=brk, fromcond=tocond, tocond=outcome, prob=prob)))
                    else  # not at a terminal outcome: still more nodes to add
                        push!(todo, vcat(seq, (sickday=brk, fromcond=tocond, tocond=outcome, prob=prob)))
                    end
                end
                break # we found the tocond as a matching fromcond
            end
        end
    end

    return done
end


function verifyprobs(seqs)
    ret = Dict(dead=>0.0, recovered=>0.0)
    allpr = 0.0

    for seq in seqs
        pr = mapreduce(x->getindex(x,:prob), *, seq)
        outcome = last(seq).tocond
        ret[outcome] += pr
        allpr += pr
    end
    return ret, allpr
end




# new trees look like this to replace the tuple as a node identifier with nested dictionaries
# "branches" dict is no longer included

 #=
instead of nodes looking like this:
    (9, 7) =>
        "probs" =>
            [0.85, 0.15]
        "outcomes" =>
            [7, 8]
        "branches" =>
            Dict("tocond" => 7, "pr" = > 0.85, "next" => (3, 3))
            Dict("tocond" => 8, "pr" = > 0.15, "next" => (3, 4))
    (9, 5) =>
        "probs" =>
            [0.8, 0.2]
        "outcomes" =>
            [3, 7]
        "branches" =>
            Dict("tocond" => 3, "pr" = > 0.8, "next" => (0, 0))
            Dict("tocond" => 7, "pr" = > 0.2, "next" => (3, 3))

a node looks like this:
    9 =>
        7 =>
            "probs" =>
                [0.85, 0.15]
            "outcomes" =>
                [7, 8]
            "branches" =>
                Dict("tocond" => 7, "pr" = > 0.85, "next" => (3, 3))
                Dict("tocond" => 8, "pr" = > 0.15, "next" => (3, 4))
        5 =>
            "probs" =>
                [0.8, 0.2]
            "outcomes" =>
                [3, 7]
            "branches" =>
                Dict("tocond" => 3, "pr" = > 0.8, "next" => (0, 0))
                Dict("tocond" => 7, "pr" = > 0.2, "next" => (3, 3))
=#

#  what a tree looks like for 5 agegrps
#= 
agegrp: 5 =>
    sickday: 25 =>
        fromcond: 7 =>
            probs: => [0.682, 0.318]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.682)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.318)
        fromcond: 8 =>
            probs: => [0.676, 0.324]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.676)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.324)
    sickday: 19 =>
        fromcond: 8 =>
            probs: => [0.49, 0.24, 0.27]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.49)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.24)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.27)
    sickday: 14 =>
        fromcond: 6 =>
            probs: => [0.7, 0.3]
            outcomes: => [3, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.7)
                Dict{Any, Any}("tocond" => 7, "next" => (25, 7), "pr" => 0.3)
        fromcond: 7 =>
            probs: => [0.7, 0.1, 0.2]
            outcomes: => [3, 7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.7)
                Dict{Any, Any}("tocond" => 7, "next" => (25, 7), "pr" => 0.1)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.2)
        fromcond: 8 =>
            probs: => [0.12, 0.67, 0.21]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.12)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.67)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.21)
    sickday: 9 =>
        fromcond: 5 =>
            probs: => [0.5, 0.5]
            outcomes: => [3, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.5)
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.5)
        fromcond: 6 =>
            probs: => [0.4, 0.6]
            outcomes: => [6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 6, "next" => (14, 6), "pr" => 0.4)
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.6)
        fromcond: 7 =>
            probs: => [0.6, 0.4]
            outcomes: => [7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.6)
                Dict{Any, Any}("tocond" => 8, "next" => (14, 8), "pr" => 0.4)
    sickday: 5 =>
        fromcond: 5 =>
            probs: => [0.1, 0.5, 0.4]
            outcomes: => [5, 6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 5, "next" => (9, 5), "pr" => 0.1)
                Dict{Any, Any}("tocond" => 6, "next" => (9, 6), "pr" => 0.5)
                Dict{Any, Any}("tocond" => 7, "next" => (9, 7), "pr" => 0.4)
agegrp: 4 =>
    sickday: 25 =>
        fromcond: 7 =>
            probs: => [0.76, 0.24]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.76)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.24)
        fromcond: 8 =>
            probs: => [0.688, 0.312]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.688)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.312)
    sickday: 19 =>
        fromcond: 8 =>
            probs: => [0.81, 0.13, 0.06]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.81)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.13)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.06)
    sickday: 14 =>
        fromcond: 6 =>
            probs: => [1.0]
            outcomes: => [3]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 1.0)
        fromcond: 7 =>
            probs: => [0.8, 0.1, 0.1]
            outcomes: => [3, 7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.8)
                Dict{Any, Any}("tocond" => 7, "next" => (25, 7), "pr" => 0.1)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.1)
        fromcond: 8 =>
            probs: => [0.165, 0.715, 0.12]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.165)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.715)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.12)
    sickday: 9 =>
        fromcond: 5 =>
            probs: => [0.62, 0.38]
            outcomes: => [3, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.62)
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.38)
        fromcond: 6 =>
            probs: => [1.0]
            outcomes: => [6]
            branches: =>
                Dict{Any, Any}("tocond" => 6, "next" => (14, 6), "pr" => 1.0)
        fromcond: 7 =>
            probs: => [0.78, 0.22]
            outcomes: => [7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.78)
                Dict{Any, Any}("tocond" => 8, "next" => (14, 8), "pr" => 0.22)
    sickday: 5 =>
        fromcond: 5 =>
            probs: => [0.15, 0.6, 0.25]
            outcomes: => [5, 6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 5, "next" => (9, 5), "pr" => 0.15)
                Dict{Any, Any}("tocond" => 6, "next" => (9, 6), "pr" => 0.6)
                Dict{Any, Any}("tocond" => 7, "next" => (9, 7), "pr" => 0.25)
agegrp: 2 =>
    sickday: 25 =>
        fromcond: 7 =>
            probs: => [0.964, 0.036]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.964)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.036)
        fromcond: 8 =>
            probs: => [0.964, 0.036]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.964)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.036)
    sickday: 19 =>
        fromcond: 8 =>
            probs: => [0.922, 0.072, 0.006]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.922)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.072)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.006)
    sickday: 14 =>
        fromcond: 6 =>
            probs: => [1.0]
            outcomes: => [3]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 1.0)
        fromcond: 7 =>
            probs: => [0.83, 0.1, 0.07]
            outcomes: => [3, 7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.83)
                Dict{Any, Any}("tocond" => 7, "next" => (25, 7), "pr" => 0.1)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.07)
        fromcond: 8 =>
            probs: => [0.474, 0.514, 0.012]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.474)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.514)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.012)
    sickday: 9 =>
        fromcond: 5 =>
            probs: => [0.85, 0.15]
            outcomes: => [3, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.85)
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.15)
        fromcond: 6 =>
            probs: => [1.0]
            outcomes: => [6]
            branches: =>
                Dict{Any, Any}("tocond" => 6, "next" => (14, 6), "pr" => 1.0)
        fromcond: 7 =>
            probs: => [0.9, 0.1]
            outcomes: => [7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.9)
                Dict{Any, Any}("tocond" => 8, "next" => (14, 8), "pr" => 0.1)
    sickday: 5 =>
        fromcond: 5 =>
            probs: => [0.2, 0.7, 0.1]
            outcomes: => [5, 6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 5, "next" => (9, 5), "pr" => 0.2)
                Dict{Any, Any}("tocond" => 6, "next" => (9, 6), "pr" => 0.7)
                Dict{Any, Any}("tocond" => 7, "next" => (9, 7), "pr" => 0.1)
agegrp: 3 =>
    sickday: 25 =>
        fromcond: 7 =>
            probs: => [0.958, 0.042]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.958)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.042)
        fromcond: 8 =>
            probs: => [0.958, 0.042]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.958)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.042)
    sickday: 19 =>
        fromcond: 8 =>
            probs: => [0.856, 0.126, 0.018]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.856)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.126)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.018)
    sickday: 14 =>
        fromcond: 6 =>
            probs: => [0.9, 0.1]
            outcomes: => [3, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.9)
                Dict{Any, Any}("tocond" => 7, "next" => (25, 7), "pr" => 0.1)
        fromcond: 7 =>
            probs: => [0.85, 0.14, 0.01]
            outcomes: => [3, 7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.85)
                Dict{Any, Any}("tocond" => 7, "next" => (25, 7), "pr" => 0.14)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.01)
        fromcond: 8 =>
            probs: => [0.776, 0.206, 0.018]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.776)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.206)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.018)
    sickday: 9 =>
        fromcond: 5 =>
            probs: => [0.9, 0.1]
            outcomes: => [3, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.9)
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.1)
        fromcond: 6 =>
            probs: => [1.0]
            outcomes: => [6]
            branches: =>
                Dict{Any, Any}("tocond" => 6, "next" => (14, 6), "pr" => 1.0)
        fromcond: 7 =>
            probs: => [0.9, 0.1]
            outcomes: => [7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.9)
                Dict{Any, Any}("tocond" => 8, "next" => (14, 8), "pr" => 0.1)
    sickday: 5 =>
        fromcond: 5 =>
            probs: => [0.2, 0.7, 0.1]
            outcomes: => [5, 6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 5, "next" => (9, 5), "pr" => 0.2)
                Dict{Any, Any}("tocond" => 6, "next" => (9, 6), "pr" => 0.7)
                Dict{Any, Any}("tocond" => 7, "next" => (9, 7), "pr" => 0.1)
agegrp: 1 =>
    sickday: 25 =>
        fromcond: 7 =>
            probs: => [0.976, 0.024]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.976)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.024)
        fromcond: 8 =>
            probs: => [0.91, 0.09]
            outcomes: => [3, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.91)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.09)
    sickday: 19 =>
        fromcond: 8 =>
            probs: => [0.891, 0.106, 0.003]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.891)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.106)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.003)
    sickday: 14 =>
        fromcond: 6 =>
            probs: => [1.0]
            outcomes: => [3]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 1.0)
        fromcond: 7 =>
            probs: => [0.85, 0.12, 0.03]
            outcomes: => [3, 7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.85)
                Dict{Any, Any}("tocond" => 7, "next" => (25, 7), "pr" => 0.12)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.03)
        fromcond: 8 =>
            probs: => [0.692, 0.302, 0.006]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.692)
                Dict{Any, Any}("tocond" => 8, "next" => (19, 8), "pr" => 0.302)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.006)
    sickday: 9 =>
        fromcond: 5 =>
            probs: => [0.9, 0.1]
            outcomes: => [3, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.9)
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.1)
        fromcond: 6 =>
            probs: => [1.0]
            outcomes: => [6]
            branches: =>
                Dict{Any, Any}("tocond" => 6, "next" => (14, 6), "pr" => 1.0)
        fromcond: 7 =>
            probs: => [0.95, 0.05]
            outcomes: => [7, 8]
            branches: =>
                Dict{Any, Any}("tocond" => 7, "next" => (14, 7), "pr" => 0.95)
                Dict{Any, Any}("tocond" => 8, "next" => (14, 8), "pr" => 0.05)
    sickday: 5 =>
        fromcond: 5 =>
            probs: => [0.4, 0.5, 0.1]
            outcomes: => [5, 6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 5, "next" => (9, 5), "pr" => 0.4)
                Dict{Any, Any}("tocond" => 6, "next" => (9, 6), "pr" => 0.5)
                Dict{Any, Any}("tocond" => 7, "next" => (9, 7), "pr" => 0.1)
=#


# what a dec_points dict looks like:

#= type Dict{Int64, Array{Tuple{Int64, Int64}, 1}}   
    5 => 
        [1, 1]
    9 => 
        [2, 1]
        [2, 2]
        [2, 3]
    14 => 
        [3, 2]
        [3, 3]
        [3, 4]
    19 => 
        [4, 4]
    25 => 
        [5, 3]
        [5, 4)
=#

# what the yaml parameter file looks like for a single agregroup
#    when it is loaded and processed, arrays will be changed to tuples

#=
1:                                          # agegrp
  [5,5]:                                      # node is [sickdayday effective, from condition]
    - {tocond: 5, next: [9, 5], pr: 0.4}
    - {tocond: 6, next: [9, 6], pr: 0.5}
    - {tocond: 7, next: [9, 7], pr: 0.1}
  [9,5]:                                      # node
    - {tocond: 3, next: [0,0], pr: 0.9}         # branch. node [0,0] denotes recovered
    - {tocond: 7, next: [14,7], pr: 0.1}        # branch
  [9,6]:                                      # node
    - {tocond: 6, next: [14, 6], pr: 1.0}       # branch
  [9,7]:
    - {tocond: 7, next: [14, 7], pr: 0.95}
    - {tocond: 8, next: [14, 8], pr: 0.05}
  [14,6]:
    - {tocond: 3, next: [0,0], pr: 1.0}
  [14,7]:
    - {tocond: 3, next: [0,0], pr: 0.85}
    - {tocond: 7, next: [25, 7], pr: 0.12}
    - {tocond: 8, next: [19, 8], pr: 0.03}
  [14,8]:
    - {tocond: 3, next: [0,0], pr: 0.692}
    - {tocond: 8, next: [19,8], pr: 0.302}
    - {tocond: 4, next: [0,5], pr: 0.006}     # branch. node [0,5] denotes dead
  [19,8]:
    - {tocond: 3, next: [0,0], pr: 0.891}
    - {tocond: 8, next: [25,8], pr: 0.106}
    - {tocond: 4, next: [0,5], pr: 0.003}
  [25,7]:
    - {tocond: 3, next: [0,0], pr: 0.976}
    - {tocond: 4, next: [0,5], pr: 0.024}
  [25,8]:
    - {tocond: 3, next: [0,0], pr: 0.91}
    - {tocond: 4, next: [0,5], pr: 0.09}
2:                                           # start of next agegrp
=#