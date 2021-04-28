
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
        for node in keys(trees[agegrp])
            a = node[1]
            b = node[2]
                probs = [branch["pr"] for branch in trees[agegrp][node]]
                outcomes = [branch["tocond"] for branch in trees[agegrp][node]]
                branches = [branch for branch in trees[agegrp][node]]
            if haskey(newdict[agegrp], a)
                newdict[agegrp][a][b] = Dict("probs"=>probs, "outcomes"=>outcomes, "branches"=>branches)
            else
                newdict[agegrp][a]=Dict()
                newdict[agegrp][a][b] = Dict("probs"=>probs, "outcomes"=>outcomes, "branches"=>branches)
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
    for k in keys(sort(tree))
        println(k)
        for item in tree[k]
            println("   ", item)
        end
    end
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


#  what a tree looks like for a single agegrp:

#= type Dict{Array{Int64, 1}, Array{Dict, 1}}  
(5, 5) =>
    "probs"
        [0.2, 0.65, 0.15]
    "outcomes"
        [5, 6, 7]
    "branches"
        Dict("tocond" => 5, "pr" = > 0.2, "next" => (2, 1))
        Dict("tocond" => 6, "pr" = > 0.65, "next" => (2, 2))
        Dict("tocond" => 7, "pr" = > 0.15, "next" => (2, 3))
(9, 5) =>
    "probs"
        [0.8, 0.2]
    "outcomes"
        [3, 7]
    "branches"
        Dict("tocond" => 3, "pr" = > 0.8, "next" => (0, 0))
        Dict("tocond" => 7, "pr" = > 0.2, "next" => (3, 3))
(9, 6) =>
    "probs"
        [1.0]
    "outcomes"
        [6]
    "branches"
        Dict("tocond" => 6, "pr" = > 1.0, "next" => (3, 2))
(9, 7) =>
    "probs"
        [0.85, 0.15]
    "outcomes"
        [7, 8]
    "branches"
        Dict("tocond" => 7, "pr" = > 0.85, "next" => (3, 3))
        Dict("tocond" => 8, "pr" = > 0.15, "next" => (3, 4))
(14, 6) =>
    "probs"
        [1.0]
    "outcomes"
        [3]
    "branches"
        Dict("tocond" => 3, "pr" = > 1.0, "next" => (0, 0))
(14, 7) =>
    "probs"
        [0.8, 0.1, 0.1]
    "outcomes"
        [3, 7, 8]
    "branches"
        Dict("tocond" => 3, "pr" = > 0.8, "next" => (0, 0))
        Dict("tocond" => 7, "pr" = > 0.1, "next" => (5, 3))
        Dict("tocond" => 8, "pr" = > 0.1, "next" => (4, 4))
(14, 8) =>
    "probs"
        [0.45, 0.5, 0.05]
    "outcomes"
        [3, 8, 4]
    "branches"
        Dict("tocond" => 3, "pr" = > 0.45, "next" => (0, 0))
        Dict("tocond" => 8, "pr" = > 0.5, "next" => (4, 4))
        Dict("tocond" => 4, "pr" = > 0.05, "next" => (0, 5))
(19, 8) =>
    "probs"
        [0.85, 0.1, 0.05]
    "outcomes"
        [3, 8, 4]
    "branches"
        Dict("tocond" => 3, "pr" = > 0.85, "next" => (0, 0))
        Dict("tocond" => 8, "pr" = > 0.1, "next" => (5, 4))
        Dict("tocond" => 4, "pr" = > 0.05, "next" => (0, 5))
(25, 7) =>
    "probs"
        [0.9, 0.1]
    "outcomes"
        [3, 4]
    "branches"
        Dict("tocond" => 3, "pr" = > 0.9, "next" => (0, 0))
        Dict("tocond" => 4, "pr" = > 0.1, "next" => (0, 5))
(25, 8) =>
    "probs"
        [0.6, 0.4]
    "outcomes"
        [3, 4]
    "branches"
        Dict("tocond" => 3, "pr" = > 0.6, "next" => (0, 0))
        Dict("tocond" => 4, "pr" = > 0.4, "next" => (0, 5))
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