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


function setup_dt(dtfname; node_starts_filename = "nodestarts.csv")
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


# hard coded  TODO:  create a version that follows the branches
function sanity_test(tree; atol=1e-3)
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
what the result looks like:
Dict key 
    value::Array of Branch

(2, 3)
    Any[Branch(7, 7, 0.7, (3, 3), "sick"), Branch(7, 8, 0.7, (3, 4), "sick", severe")]
(3, 2)
    Any[Branch(6, 3, 1.0, (0, 0), "mild", "recovered")]
(3, 3)
    Any[Branch(7, 3, 0.7, (0, 0), "sick", recovered"), Branch(7, 8, 0.3, (4, 4), "sick", "severe")]
(2, 2)
    Any[Branch(6, 6, 1.0, (3, 2), "mild", "mild")]
(1, 1)
    Any[Branch(5, 6, 0.7, (2, 2), "nil", "mild"), Branch(5, 7, 0.3, (2, 3), "nil", "sick")]
(4, 4)
    Any[Branch(8, 3, 0.8, (0, 0), "severe", "recovered"), Branch(8, 4, 0.2, (0, 0), "severe", "dead")]
(3, 4)
    Any[Branch(8, 3, 0.2, (0, 0), "severe", "recovered"), Branch(8, 8, 0.7, (4, 4), "severe", "severe"), 
        Branch(8, 4, 0.1, (0, 0), "severe", "dead")]

=#