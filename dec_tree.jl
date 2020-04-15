
#=
dec_tree.jl
decision tree for transition
=#

using DelimitedFiles


struct Branch 
    fromcond::Int 
    tocond::Int
    pr::Float64
    next::Tuple{Int,Int}
    fromcondname::String
    tocondname::String
end


function setup_dt(dtfname; node_starts_filename = "dec_tree_starts.csv")
    arr = read_dectree_file(dtfname)
    dectrees = create_node_dict(arr)
    nodestarts = build_nodestarts(node_starts_filename)
    return (dt=dectrees, starts=nodestarts)
end

function read_dectree_file(fname)
    # return the data at index 1, not the header at index 2
    arr=readdlm(fname, ',', header=true, comments=true, comment_char='#')[1]
end


function create_node_dict(arr::Array) # this wants to be recursive--another time...
    arrdectrees = []  # array of decision trees (dicts)
    ages = unique!(arr[:,1])
    for agegrp in ages
        xx = view(arr, arr[:,1] .== agegrp, :) # view on one agegrp
        dectree = Dict{Tuple, Array}()  # Dict of nodes.  each node is an array of branches 
        nodelist = unique!(xx[:,2])
        for strnode in nodelist
            node = eval(Meta.parse(strnode))  # convert the string to a tuple
            rr = xx[xx[:,2].==strnode,:]  # rows for the branches of this node
            arrbranches = []   # array of branches
            for br in 1:size(rr,1)  # make a branch
                # field values
                fromcond = eval(Symbol(rstrip(rr[br,3])))  # convert a string to a variable name, which holds an Int
                tocond = eval(Symbol(rstrip(rr[br,4])))
                pr = rr[br,5]    # float
                next = eval(Meta.parse(rr[br,6]))  # convert a string to a tuple
                fromcondname = condnames[fromcond]  # lookup text name for the numeric index
                tocondname = condnames[tocond]    

                newbr = Branch(fromcond, tocond, pr, next, fromcondname, tocondname)
                push!(arrbranches, newbr)        
            end
            dectree[node] = arrbranches  # put all branches for this node in the dict for the agegrp
        end
        push!(arrdectrees, dectree) # done with all of the nodes, put the tree in the array
    end
    return arrdectrees
end


function build_nodestarts(fname)
    arr = readdlm(fname, ',', header=true, comments=true, comment_char='#')[1]
    ns = Dict{Int, Array{Tuple{Int,Int},1}}()
    for r in eachrow(arr)
        day = r[1]
        tpl = eval(Meta.parse(r[2]))
        @assert day in (1:19) "day must be in 1:19"
        @assert typeof(tpl) <: Tuple{Int, Int} "Tuple input as text must convert to a tuple of 2 ints"
        if haskey(ns, day)
            push!(ns[day], tpl)    
        else
            ns[day] = [tpl]
        end
    end
    return ns
end


function display_tree(tree)
    for k in keys(sort(tree))
        println(k)
        for item in tree[k]
            println("   ", item)
        end
    end
end

# built this by eyeballing the correct definition of the tree structure;
#   a really smart person would write code to follow the branches in a tree
paths2 =    [
               [(1,1),(2,1),(0,0)],
               [(1,1),(2,1),(3,3),(0,0)],
               [(1,1),(2,1),(3,3),(4,4),(0,0)],
               [(1,1),(2,1),(3,3),(4,4),(0,5)],
               [(1,1),(2,2),(3,2),(0,0)],
               [(1,1),(2,3),(3,3),(0,0)],
               [(1,1),(2,3),(3,3),(4,4),(0,0)],
               [(1,1),(2,3),(3,3),(4,4),(0,5)],
               [(1,1),(2,3),(3,4),(0,0)],
               [(1,1),(2,3),(3,4),(0,5)],
               [(1,1),(2,3),(3,4),(4,4),(0,0)],
               [(1,1),(2,3),(3,4),(4,4),(0,5)],
            ]


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

function get_the_probs(seq, tree)
    probs = []
    for cnt in 1:length(seq)-1
        it1, it2 = seq[cnt], seq[cnt+1]
        node = tree[it1]
        for br in node
            if br.next == it2
                push!(probs, br.pr)
            end
        end
    end
    if seq[end] == (0,0)
        probs = ("recovered", probs)
    elseif seq[end] == (0,5)
        probs = ("dead", probs)
    end
    return probs
end

# hard coded  TODO:  create a version that follows the branches
function sanity_test_1(tree; atol=1e-3)
    r = Float64[]  # recovered
    d = Float64[]  # dead

    # (1,1).1 * (2,2) * (3,2)                   recovered
    rec = tree[(1,1)][1].pr * tree[(2,2)][1].pr * tree[(3,2)][1].pr
    push!(r,rec)

    # (1,1).2 * (2,3).1 * (3,3).1               recovered
    rec = tree[(1,1)][2].pr * tree[(2,3)][1].pr * tree[(3,3)][1].pr
    push!(r,rec)

    # (1,1).2 * (2,3).1 * (3,3).2 * (4,4).1     recovered, total with same endpoint below
    rec1 = tree[(1,1)][2].pr * tree[(2,3)][1].pr * tree[(3,3)][2].pr * tree[(4,4)][1].pr

    # (1,1).2 * (2,3).2  * (3,4).2 * (4,4).1    recovered
    rec2 = tree[(1,1)][2].pr * tree[(2,3)][2].pr * tree[(3,4)][2].pr * tree[(4,4)][1].pr
    push!(r,rec1 + rec2)

    # (1,1).2 * (2,3).2  * (3,4).1              recovered
    rec = tree[(1,1)][2].pr * tree[(2,3)][2].pr * tree[(3,4)][1].pr
    push!(r,rec)

    # (1,1).2 * (2,3).2  * (3,4).2 * (4,4).2    dead, total with same endpoint belowq
    dead1 = tree[(1,1)][2].pr * tree[(2,3)][1].pr * tree[(3,3)][2].pr * tree[(4,4)][2].pr

    # (1,1).2 * (2,3).1 * (3,3).2 * (4,4).2     dead
    dead2 = tree[(1,1)][2].pr * tree[(2,3)][1].pr * tree[(3,3)][2].pr * tree[(4,4)][2].pr
    push!(d,dead1 + dead2)

    # (1,1).2 * (2,3).2  * (3,4).3              dead
    dead = tree[(1,1)][2].pr * tree[(2,3)][2].pr * tree[(3,4)][3].pr
    push!(d,dead)

    totprob = sum(r) + sum(d)
    tstresult = isapprox(totprob, 1.0, atol=atol)
    println("Total probability: ", totprob)
    println("Dead probability:  ", sum(d))
    println("Rec probability:   ", sum(r))
    if tstresult
        println("Passed: probabilities of all outcomes sum to 1.0")
    else
        println("Failed: probabilities of all outcomes do not sum to 1.0")
    end

    return r, d
end

#=
    mutable struct step
        node::Tuple
        pr::Float64
        tocond::Int
        path::Dict
    end


    function find_all_leaves(tree)
        leaves = []
        for (node,br_arr) in tree
            for br in br_arr
                if br.next[1] == 0
                    push!(leaves, (node, br.tocondname, br.pr))
                end
            end
        end
        return leaves
    end

    function find_parents_of_leaves(leaves, tree)
        parents = Dict()
        for leaf in leaves
            lfnode = leaf[1]
            for (node,br_arr) in tree
                for br in br_arr
                    if br.next == lfnode
                        haskey(parents, (lfnode, node)) || (parents[(lfnode, node)] = br.pr)
                    end
                end
            end
        end
        return parents
    end

    function to_root(p_of_l, tree)
        paths = []
        for k in keys(p_of_l)
            kl = k[2]
            for (node,br_arr) in tree
                for br in br_arr
                    if br == (1,1)  # already at root
                    end
                    if br.next == kl
                        push!(paths, ((kl, br.pr), k))
                    end
                end
            end
        end
        return paths
    end
=#


#=
    function across(node, tree, outcomes = Dict())
        println("enter across ", node)

        across = [tree[node][i].next for i in 1:length(tree[node])]
        outcomes[node] = [tree[node][i].pr for i in 1:length(tree[node])]
        println("across ", across)
        for nextnode in across
            outcomes = Dict
            br = tree[node]
            println("across loop ", node, " ", nextnode)
            res = br.pr * back[1]
            push!(outcome, (back[2], res))
            println(" before end of across loop ", outcome)
        end

    end

    function down(node, tree)
        println("enter down ", node)

        if node == (0,0)
            return (1.0, "recovered")
        elseif node == (0,5)
            return (1.0, "dead")
        elseif length(tree[node]) > 1
             across(node, tree)
        elseif tree[node][1].next == (0,0)
            br = tree[node][1]
            return (br.pr, br.tocondname)
        else
            println("down else ", tree[node][1])
            br = tree[node][1]
            return br.pr * down(br.next, tree)
        end

    end
=#




# function dectree_sanity(tree, top=(1,1); atol=1e-3)
#     numbranches = Dict()
#     for k in keys(tree)
#         numbranches[k] = length(tree[k])
#     end
#     results = []
#     node = (1,1)
#     prevnode = (1,1)
#     prod = 1.0
#     while prevnode != (0,0)
#         numbranches = length(tree[node])
#         i = 1
#         while i <= numbranches
#             prod *= prod * tree[node][i].pr 
#             nxt = tree[node][i].next
#             outcome = tree[node][i].tocondname
#             if nxt == (0,0)
#                 push!(results, (outcome, prod))
#                 prod= 1
#             end
#             i += 1
#         end

#             prevnode = node
#             node = nxt



#     end
#     return numbranches
# end


#=
    (1, 1)
       CovidSim.Branch(5, 5, 0.2, (2, 1), "nil", "nil")
       CovidSim.Branch(5, 6, 0.65, (2, 2), "nil", "mild")
       CovidSim.Branch(5, 7, 0.15, (2, 3), "nil", "sick")
    (2, 1)
       CovidSim.Branch(5, 3, 0.8, (0, 0), "nil", "recovered")
       CovidSim.Branch(5, 7, 0.2, (3, 3), "nil", "sick")
    (2, 2)
       CovidSim.Branch(6, 6, 1.0, (3, 2), "mild", "mild")
    (2, 3)
       CovidSim.Branch(7, 7, 0.85, (3, 3), "sick", "sick")
       CovidSim.Branch(7, 8, 0.15, (3, 4), "sick", "severe")
    (3, 2)
       CovidSim.Branch(6, 3, 1.0, (0, 0), "mild", "recovered")
    (3, 3)
       CovidSim.Branch(7, 3, 0.85, (0, 0), "sick", "recovered")
       CovidSim.Branch(7, 8, 0.15, (4, 4), "sick", "severe")
    (3, 4)
       CovidSim.Branch(8, 3, 0.45, (0, 0), "severe", "recovered")
       CovidSim.Branch(8, 8, 0.5, (4, 4), "severe", "severe")
       CovidSim.Branch(8, 4, 0.05, (0, 0), "severe", "dead")
    (4, 4)
       CovidSim.Branch(8, 3, 0.9, (0, 0), "severe", "recovered")
       CovidSim.Branch(8, 4, 0.1, (0, 0), "severe", "dead")
=#


#=
    function dectree_sanity(tree; atol=1e-3)
        allnodes = keys(tree)
        # find every node that has at least one leaf
        hasleaves = [has_leaf(n,tree) for n in allnodes if !(isnothing(has_leaf(n,tree)))]

        # find every parent of each node that has leaves
        leafparents = []
        for node in hasleaves
            println(node[1], " ", parentsof(node[1], tree))
        end
    end

    function has_leaf(node, tree)  # a node is an array of branches
        leaves = []
        for br in tree[node]
            if br.next == (0,0)
                push!(leaves, br.tocondname)
            end
        end
        return isempty(leaves) ? nothing  : (node, leaves)
    end    

    function parentsof(node, tree)
        ret = []
        nlvl = node[1]
        uppers = [nn for nn in keys(tree) if nn[1] == nlvl -1]
        for nn in uppers
            for br in tree[nn] # an array of branches
                if br.next == node
                    push!(ret, nn)
                end
            end
        end
        return ret
    end
=#


#=
    function dectree_sanity(tree, node=(1,1), paths=Dict(); atol=1e-3)
        for br in tree[node]
            if br.next == (0,0)
                return step((0,0), br.pr, br.tocond)
            else
                paths[node] = step(br.next, br.pr, 0, dectree_sanity(tree, node))
                return dectree_sanity(tree, node=br.next, )

    end
=#


#= 
    function dectree_sanity(tree; atol=1e-3)
        remaining = Dict()
        for ntup in keys(tree)
            for br in tree[ntup]
                lvl = ntup[1]
                if  haskey(remaining, lvl)
                    push!(remaining[lvl], (ntup, br.next, br.pr))
                else
                    remaining[lvl] = [(ntup, br.next, br.pr)]
                end
            end
        end
        outs = 0
        prods = Dict()
        for lvl = 1:2    # fix constant later
            println(remaining[lvl])
            for fac in remaining[lvl]
                println(fac)
                if fac[2] != (0,0)
                    if  haskey(prods, outs)
                        push!(prods[outs], fac[3])
                    else
                        outs += 1
                        prods[outs] = [fac[3]]
                    end
                else
                    if  haskey(prods, outs)
                        push!(prods[outs], fac[3])
                    else

                        prods[outs] = [fac[3]]
                    end      
                    
                end
            end
            println(prods)
        end
        return remaining,prods
    end
=#

# function path_outcomes(paths, tree)
#     outcomes = []
#     println("no. of paths ", length(paths))
#     for path in paths
#         display(path)
#         pr = 1.0
#         shortstop = length(path) -1
#         for i in 1:shortstop
#             # [(1,1),(2,1),(0,0)],
#             if i == shortstop
#                 node = tree[path[i]]
#                 for br in node
#                     if br.next[1] == 0
#                         pr *= br.pr 
#                     end
#                 end
#                 endnode = path[i+1]
#                 if endnode == (0,0)
#                     # pr *= br.pr
#                     push!(outcomes,("recovered", round(pr, digits=4)))
#                     # println(outcomes)
#                     break
#                 elseif endnode == (0,5)
#                     # pr *= br.pr
#                     push!(outcomes, ("dead", round(pr,digits=4)))
#                     # println(outcomes)
#                     break
#                 end
#             end
#             for br in tree[path[i]]   # DON'T WALK THROUGH; JUST GET THE ONE WE NEED
#                 # println(br)
#                 # @assert false
#                 if br.next == path[i+1]
#                     pr *= br.pr 
#                 end # if
#             end # br loop
#         end # path node loop
#     end # path loop
#     return outcomes
# end # function


# function sanity_test_paths(paths, tree)
#     outcomes = path_outcomes(paths, tree)
#     deadtotalpr = 0.0
#     recoveredtotalpr =  0.0
#     for item in outcomes
#         if item[1] == "recovered"
#             recoveredtotalpr += item[2]
#         else
#             deadtotalpr += item[2]
#         end
#     end
#     return (deadpr=deadtotalpr, recoveredpr = recoveredtotalpr, outcomes=outcomes)
# end