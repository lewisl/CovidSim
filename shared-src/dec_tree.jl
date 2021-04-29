
#############################################################
# dec_tree.jl
# decision tree for transition
#############################################################


function setup_dt(dtfilename)
    trees = YAML.load_file(dtfilename)
    # next: change 2nd level keys from 2 item array{Int} [9, 5] to Tuple{Int, Int} (9,5)
    trees = Dict(i => Dict(Tuple(k)=>trees[i][k] for k in keys(trees[i])) for i in keys(trees))

    # next: change the type of next node item from array{Int} [25, 8] to Tuple{Int, Int} (25, 8)
    for agegrp in agegrps
        for (k,v) in trees[agegrp]
           for item in v
                item["next"] = Tuple(item["next"])
            end
        end
    end


    # pre-calculate the array of probabilities for all branches at a node
    # pre-calculate the array of outcome conditions ("tocond") for all branches at a node

    newdict = Dict()
    for agegrp in agegrps
        newdict[agegrp] = Dict()
        for node in keys(trees[agegrp])  # node is (lag, fromcond)
            lag       = node[1]
            fromcond  = node[2]
                probs = [branch["pr"] for branch in trees[agegrp][node]]
                outcomes = [branch["tocond"] for branch in trees[agegrp][node]]
                branches = [branch for branch in trees[agegrp][node]]
            if haskey(newdict[agegrp], lag)
                newdict[agegrp][lag][fromcond] = Dict("probs"=>probs, "outcomes"=>outcomes, "branches"=>branches)
            else
                newdict[agegrp][lag]=Dict()
                newdict[agegrp][lag][fromcond] = Dict("probs"=>probs, "outcomes"=>outcomes, "branches"=>branches)
            end
        end
    end
    newdict = Dict(i=>sort(newdict[i], rev=true) for i in agegrps) 

    lags_by_age = Dict{Int,Array{Int,1}}()  # empty
    fromconds_by_age = Dict{Int,Array{Int,1}}()  # empty
    for agegrp in agegrps
        lags_by_age[agegrp] = [k[1] for k in collect(keys(trees[agegrp]))]
        fromconds_by_age[agegrp] = [k[2] for k in collect(keys(trees[agegrp]))]
    end

    decpoints = Dict{Int,Array{Int, 1}}()
    for i in agegrps
        decpoints[i] = unique([k[1] for k in keys(trees[i])])
    end

    return Dict("dt"=>newdict, "decpoints"=>decpoints, "lags"=>lags_by_age, "fromconds"=>fromconds_by_age)
end


function display_tree(tree)
    for agegrp in keys(tree)
        agetree = tree[agegrp]
        println("agegrp: ", agegrp, " =>")
        for lag in keys(agetree)
            lagtree = agetree[lag]
            println("    lag: ", lag, " =>")
            for fromcond in keys(lagtree)
                condtree = lagtree[fromcond]
                println("        fromcond: ", fromcond, " =>")
                print("            probs: => ")
                println(condtree["probs"])
                #
                print("            outcomes: => ")
                println(condtree["outcomes"])
                #
                println("            branches: =>")
                for branch in keys(condtree["branches"])
                    println("                ", condtree["branches"][branch])   
                end
            end  # for fromcond
        end  # for lag
    end   # for agegrp     
end


function walktree(dt, top)
    done = []
    todo = [[top]]
    while !isempty(todo)
        currentpath = popfirst!(todo)
        endnode = currentpath[end]
        for br in dt[endnode]["branches"]
            # if br.next[1] == 0
            if br["next"][1] == 0
                push!(done, vcat(currentpath, [br["next"]]))  # append without modifying currentpath
            else
                push!(todo, vcat(currentpath, [br["next"]]))   
            end
        end
    end
    return done
end


function sanity_test_all(trees)
    tbl = zeros(length(trees),4)
    for (i, tree) in trees
        paths = walktree(tree, (5,5))
        res = sanity_test(paths, tree)
        tbl[i, :] .= [i, res.total, res.recovered, res.dead]
    end
    return tbl
end

function sanity_test_all(dtfname::String)
    trees, decpoints = setup_dt(dtfname)
    sanity_test_all(trees)
end

# TODO: check that probs of all branches at a node add to one

function sanity_test(paths, tree)
    probs = []
    outcomes = []
    deadpr = 0.0
    recoveredpr = 0.0
    for path in paths
        prs = get_the_probs(path, tree)
        res = prs[1]
        prs = prs[2]
        push!(probs,(res, prod(prs)))
    end
    for item in probs
        if item[1] == "recovered"
            recoveredpr += item[2]
        else
            deadpr += item[2]
        end
    end
    return (recovered=recoveredpr, dead=deadpr, total=recoveredpr+deadpr, probs=probs)
end

function get_the_probs(path, tree)
    probs = []
    for cnt in 1:length(path)-1
        it1, it2 = path[cnt], path[cnt+1]
        node = tree[it1]
        for br in node["branches"]
            if br["next"] == it2
                push!(probs, br["pr"])
            end
        end
    end
    if path[end] == (0,0)
        probs = ["recovered", probs]
    elseif path[end] == (0,5)
        probs = ["dead", probs]
    else
        error("didn't work")
    end
    return probs
end


# new trees look like this to replace the tuple as a node identifier with nested dictionaries

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
    lag: 25 =>
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
    lag: 19 =>
        fromcond: 8 =>
            probs: => [0.49, 0.24, 0.27]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.49)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.24)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.27)
    lag: 14 =>
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
    lag: 9 =>
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
    lag: 5 =>
        fromcond: 5 =>
            probs: => [0.1, 0.5, 0.4]
            outcomes: => [5, 6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 5, "next" => (9, 5), "pr" => 0.1)
                Dict{Any, Any}("tocond" => 6, "next" => (9, 6), "pr" => 0.5)
                Dict{Any, Any}("tocond" => 7, "next" => (9, 7), "pr" => 0.4)
agegrp: 4 =>
    lag: 25 =>
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
    lag: 19 =>
        fromcond: 8 =>
            probs: => [0.81, 0.13, 0.06]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.81)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.13)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.06)
    lag: 14 =>
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
    lag: 9 =>
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
    lag: 5 =>
        fromcond: 5 =>
            probs: => [0.15, 0.6, 0.25]
            outcomes: => [5, 6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 5, "next" => (9, 5), "pr" => 0.15)
                Dict{Any, Any}("tocond" => 6, "next" => (9, 6), "pr" => 0.6)
                Dict{Any, Any}("tocond" => 7, "next" => (9, 7), "pr" => 0.25)
agegrp: 2 =>
    lag: 25 =>
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
    lag: 19 =>
        fromcond: 8 =>
            probs: => [0.922, 0.072, 0.006]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.922)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.072)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.006)
    lag: 14 =>
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
    lag: 9 =>
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
    lag: 5 =>
        fromcond: 5 =>
            probs: => [0.2, 0.7, 0.1]
            outcomes: => [5, 6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 5, "next" => (9, 5), "pr" => 0.2)
                Dict{Any, Any}("tocond" => 6, "next" => (9, 6), "pr" => 0.7)
                Dict{Any, Any}("tocond" => 7, "next" => (9, 7), "pr" => 0.1)
agegrp: 3 =>
    lag: 25 =>
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
    lag: 19 =>
        fromcond: 8 =>
            probs: => [0.856, 0.126, 0.018]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.856)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.126)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.018)
    lag: 14 =>
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
    lag: 9 =>
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
    lag: 5 =>
        fromcond: 5 =>
            probs: => [0.2, 0.7, 0.1]
            outcomes: => [5, 6, 7]
            branches: =>
                Dict{Any, Any}("tocond" => 5, "next" => (9, 5), "pr" => 0.2)
                Dict{Any, Any}("tocond" => 6, "next" => (9, 6), "pr" => 0.7)
                Dict{Any, Any}("tocond" => 7, "next" => (9, 7), "pr" => 0.1)
agegrp: 1 =>
    lag: 25 =>
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
    lag: 19 =>
        fromcond: 8 =>
            probs: => [0.891, 0.106, 0.003]
            outcomes: => [3, 8, 4]
            branches: =>
                Dict{Any, Any}("tocond" => 3, "next" => (0, 0), "pr" => 0.891)
                Dict{Any, Any}("tocond" => 8, "next" => (25, 8), "pr" => 0.106)
                Dict{Any, Any}("tocond" => 4, "next" => (0, 5), "pr" => 0.003)
    lag: 14 =>
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
    lag: 9 =>
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
    lag: 5 =>
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
  [5,5]:                                      # node is [lagday effective, from condition]
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