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


function read_dectree_file(fname)
    # return the data at index 1, not the header at index 2
    arr=readdlm(fname, ',', header=true, comments=true, comment_char='#')[1]
end


function create_node_dict(arr::Array)
    arrdectrees = []  # array of decision trees (dicts)
    ages = unique!(arr[:,1])
    for agegrp in ages
        xx = view(arr, arr[:,1].== agegrp, :) # view on one agegrp
        dectree = Dict{Tuple, Array}()  # Dict of nodes.  each node is an array of branches 
        nodelist = unique!(xx[:,2])
        for strnode in nodelist
            node = eval(Meta.parse(strnode)) 
            rr = xx[xx[:,2].==strnode,:]  # rows for the branches of this node
            arrbranches = []   # array of branches
            for br in 1:size(rr,1)  # make a branch
                fromcond = eval(Symbol(rstrip(rr[br,3])))
                tocond = eval(Symbol(rstrip(rr[br,4])))
                pr = rr[br,5]
                next = eval(Meta.parse(rr[br,6]))
                fromcondname = condnames[fromcond]
                tocondname = condnames[tocond]    
                newbr = Branch(fromcond, tocond, pr, next, fromcondname, tocondname)
                push!(arrbranches, newbr)        
            end
            dectree[node] = arrbranches
        end
        push!(arrdectrees, dectree)
    end
    return arrdectrees
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
    println("rec:   ", rec)
    push!(r,rec)
    # (1,1).2 * (2,3).1 * (3,3).1               recovered
    rec = tree[(1,1)][2].pr * tree[(2,3)][1].pr * tree[(3,3)][1].pr
    println("rec:   ", rec)
    push!(r,rec)
    # (1,1).2 * (2,3).1 * (3,3).2 * (4,4).1     
    rec1 = tree[(1,1)][2].pr * tree[(2,3)][1].pr * tree[(3,3)][2].pr * tree[(4,4)][1].pr


    # (1,1).2 * (2,3).2  * (3,4).2 * (4,4).1    recovered
    rec2 = tree[(1,1)][2].pr * tree[(2,3)][2].pr * tree[(3,4)][2].pr * tree[(4,4)][1].pr
    println("rec:   ", rec1 + rec2)
    push!(r,rec1 + rec2)

    # (1,1).2 * (2,3).2  * (3,4).1              recovered
    rec = tree[(1,1)][2].pr * tree[(2,3)][2].pr * tree[(3,4)][1].pr
    println("rec:   ", rec)
    push!(r,rec)


    # (1,1).2 * (2,3).2  * (3,4).2 * (4,4).2    
    dead1 = tree[(1,1)][2].pr * tree[(2,3)][1].pr * tree[(3,3)][2].pr * tree[(4,4)][2].pr

    # (1,1).2 * (2,3).1 * (3,3).2 * (4,4).2     dead
    dead2 = tree[(1,1)][2].pr * tree[(2,3)][1].pr * tree[(3,3)][2].pr * tree[(4,4)][2].pr
    println("dead:  ", dead2+dead2)
    push!(d,dead1 + dead2)

    # (1,1).2 * (2,3).2  * (3,4).3              dead
    dead = tree[(1,1)][2].pr * tree[(2,3)][2].pr * tree[(3,4)][3].pr
    println("dead:  ", dead)
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
    and value::Array of Branch

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