#=
tests.jl
decision tree for transition

=#

using DelimitedFiles

# condition_outcome columns
const unexposed = 1
const infectious = 2
const recovered = 3
const dead = 4
const nil = 5
const mild = 6
const sick = 7
const severe = 8

const condnames = Dict(1=>"unexposed",2=>"infectious",3=>"recovered", 4=>"dead",
                       5=>"nil",6=>"mild",7=>"sick",8=>"severe")


struct Branch 
    fromcond::Int 
    tocond::Int
    pr::Float64
    next::Tuple{Int,Int}
    fromcondname::String
    tocondname::String
end

function create_node_dict(dat::Array)
    nrows = size(dat,1)
    i = 0    
    arrdectrees = []  # array of decision trees
    arrnode = []   # reused array of nodes to force function scope
    dectree = Dict{Tuple, Array}()  # reused Dict arrays of nodes to force function scope
    prev_agegrp = 0
    prev_node = (0,0)
    while i < nrows
        i += 1
                @debug println("  i: ", i)
        agegrp = dat[i,1]
        node = eval(Meta.parse(dat[i,2])) # node test
                @debug println("   agegrp: ", agegrp, " prev_agegrp ", prev_agegrp)
        if agegrp != prev_agegrp
                    @debug "  Finishing dectree, starting new"
            if i != 1 
                push!(arrdectrees, dectree) # done with old
            end
            prev_agegrp = agegrp
            dectree = Dict{Tuple, Array}() # start new
        end
                @debug println("      node: ", node, "   prev_node: ", prev_node, " ", prev_node == node)
        if node != prev_node
                    @debug "  Finishing node, starting new"
            if i != 1
                        @debug println("      ", prev_node)
                        @debug println("      ", arrnode)
                dectree[prev_node] = arrnode  # done with the old
                        @debug println("      ", keys(dectree))
            end
            prev_node = node
            arrnode = []  # new
        end

        # assign struct field values
            fromcond = eval(Symbol(rstrip(dat[i,3])))
            tocond = eval(Symbol(rstrip(dat[i,4])))
            pr = dat[i,5]
            next = eval(Meta.parse(dat[i,6]))
            fromcondname = condnames[fromcond]
            tocondname = condnames[tocond]
                    @debug println("         branch: ", cond, " ", pr, " ", next, " ", condname)
        br = Branch(fromcond, tocond, pr, next, fromcondname, tocondname)
                    @debug println("         branch obj: ", br)
        push!(arrnode, br)  
    end

    # for end of file
        @debug println("    ", prev_node)
        dectree[prev_node] = arrnode
        push!(arrdectrees, dectree) 

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

function read_dectree_file(fname)
    arr=readdlm(fname, ',', header=true, comments=true, comment_char='#')
end


#=
what the result looks like:
Dict key and value

(2, 3)
Any[Branch(7, 7, 0.7, (3, 3), "sick"), Branch(7, 8, 0.7, (3, 4), "severe")]
(3, 2)
Any[Branch(6, 3, 1.0, (0, 0), "recovered")]
(3, 3)
Any[Branch(7, 3, 0.7, (0, 0), "recovered"), Branch(7, 8, 0.3, (4, 4), "severe")]
(2, 2)
Any[Branch(6, 6, 1.0, (3, 2), "mild")]
(1, 1)
Any[Branch(5, 6, 0.7, (2, 2), "mild"), Branch(5, 7, 0.3, (2, 3), "sick")]
(4, 4)
Any[Branch(8, 3, 0.8, (0, 0), "recovered"), Branch(8, 4, 0.2, (0, 0), "dead")]
(3, 4)
Any[Branch(8, 3, 0.2, (0, 0), "recovered"), Branch(8, 8, 0.7, (4, 4), "severe"), Branch(8, 4, 0.1, (0, 0), "dead")]

=#