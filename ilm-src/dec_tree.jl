
#############################################################
# dec_tree.jl
# decision tree for transition
#############################################################

#= TODO
- add sanity check that all branches at a node have same from condition

=#



using DelimitedFiles

# columns of dec_tree_csv files 
const agegrp_col = 1
const node_col = 2
const lag_col = 3
const from_col = 4
const to_col = 5
const prob_col = 6
const next_col = 7


function setup_dt_old(dtfname)
    arr = read_dectree_file(dtfname)
    dectrees = create_dec_trees(arr)
    all_decpoints = reduce(merge, dectrees[agegrp].dec_points for agegrp in agegrps)
    return dectrees, all_decpoints
end


function setup_dt(dtfilename)
    treedict = YAML.load_file(dtfilename)
    dt = Dict(i=>sort(treedict[i], rev=true) for i in agegrps) 

    # pre-calculate the array of probabilities for all branches at a node
    # pre-calculate the array of outcome conditions ("tocond") for all branches at a node

    for agegrp in agegrps
        tree = dt[agegrp]
        for node in keys(tree)
            probs = [branch["pr"] for branch in tree[node]]
            outcomes = [branch["tocond"] for branch in tree[node]]
            branches = [branch for branch in tree[node]]
            tree[node] = Dict("probs"=>probs, "outcomes"=>outcomes, "branches"=>branches)
        end
    end

    lags_by_age = Dict{Int,Array{Int,1}}()  # empty
    fromconds_by_age = Dict{Int,Array{Int,1}}()  # empty
    for agegrp in agegrps
        lags_by_age[agegrp] = [k[1] for k in collect(keys(dt[agegrp]))]
        fromconds_by_age[agegrp] = [k[2] for k in collect(keys(dt[agegrp]))]
    end

    decpoints = Dict{Int,Array{Int, 1}}()
    for i in agegrps
        decpoints[i] = unique([k[1] for k in keys(dt[i])])
    end

    return Dict("dt"=>dt, "decpoints"=>decpoints, "lags"=>lags_by_age, "fromconds"=>fromconds_by_age)
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
        paths = walktree(tree,[5,5])
        res = sanity_test(paths, tree)
        tbl[i, :] .= [i, res.total, res.recovered, res.dead]
    end
    return tbl
end

function sanity_test_all(dtfname::String)
    trees, decpoints = setup_dt(dtfname)
    sanity_test_all(trees)
end


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
    if path[end] == [0,0]
        probs = ("recovered", probs)
    elseif path[end] == [0,5]
        probs = ("dead", probs)
    end
    return probs
end


#  what a tree looks like:

#= type Dict{Array{Int64, 1}, Array{Dict, 1}}  
[5, 5] =>
   Dict("tocond" => 5, "pr" = > 0.2, "next" => [2, 1])
   Dict("tocond" => 6, "pr" = > 0.65, "next" => [2, 2])
   Dict("tocond" => 7, "pr" = > 0.15, "next" => [2, 3])
[9, 5] =>
   Dict("tocond" => 3, "pr" = > 0.8, "next" => [0, 0])
   Dict("tocond" => 7, "pr" = > 0.2, "next" => [3, 3])
[9, 6] =>
   Dict("tocond" => 6, "pr" = > 1.0, "next" => [3, 2])
[9, 7] =>
   Dict("tocond" => 7, "pr" = > 0.85, "next" => [3, 3])
   Dict("tocond" => 8, "pr" = > 0.15, "next" => [3, 4])
[14, 6] =>
   Dict("tocond" => 3, "pr" = > 1.0, "next" => [0, 0])
[14, 7] =>
   Dict("tocond" => 3, "pr" = > 0.8, "next" => [0, 0])
   Dict("tocond" => 7, "pr" = > 0.1, "next" => [5, 3])
   Dict("tocond" => 8, "pr" = > 0.1, "next" => [4, 4])
[14, 8] =>
   Dict("tocond" => 3, "pr" = > 0.45, "next" => [0, 0])
   Dict("tocond" => 8, "pr" = > 0.5, "next" => [4, 4])
   Dict("tocond" => 4, "pr" = > 0.05, "next" => [0, 5])
[19, 8] =>
   Dict("tocond" => 3, "pr" = > 0.85, "next" => [0, 0])
   Dict("tocond" => 8, "pr" = > 0.1, "next" => [5, 4])
   Dict("tocond" => 4, "pr" = > 0.05, "next" => [0, 5])
[25, 7] =>
   Dict("tocond" => 3, "pr" = > 0.9, "next" => [0, 0])
   Dict("tocond" => 4, "pr" = > 0.1, "next" => [0, 5])
[25, 8] =>
   Dict("tocond" => 3, "pr" = > 0.6, "next" => [0, 0])
   Dict("tocond" => 4, "pr" = > 0.4, "next" => [0, 5])
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
