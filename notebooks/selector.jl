# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# %%
using TypedTables
using StatsBase
using BenchmarkTools
using Random
using PrettyPrint

# %%
function create_peeps(n)
   ret = Table(
        pp = collect(1:n),
        touches = zeros(Int,n),
        contacts = zeros(Int, n),
        status = zeros(Int, n)
        )  
end

# %%
function static_setup(pop, n_touches, n_unexp)
    for i in eachindex(pop)
        pop.touches[i] = n_touches
    end
    pop.status[1:n_unexp] .= 1
    start_at = n_unexp + 1
    percat = round(Int, (length(pop)-start_at)/3)
    for j in 2:3
        n = start_at + percat - 1
        pop.status[start_at:n] .= j
        start_at += percat
    end
    pop.status[start_at:end] .= 4
    shuffle!(pop.status)
    return Vector{Int}(undef, n_touches)
end


# %%
function run_touches!(pop, touchtrack, quiet=true)
    quiet || println("starting...")
    touchcount = 0
    alive = findall(pop.status .!= 4) # the living can't contact the dead and visa versa
    sick = findall(pop.status == 2)
    for pass in [1,2]
        # pass 1: do all alive with only 2 tries
        # pass 2: do sick people for all remaining contacts
        if pass == 1
            group = alive
        elseif pass == 2
            group = sick
        end
        for p in group
            touchtrack[:] .= 0
            n = 0
            tries = 0
            maxtries = pass == 1 ? 2 : pop.touches[p]
            while pop.contacts[p] < pop.touches[p]
                
                for mycontact in sample(alive, maxtries)
                    tries += 1

                    # rejected contacts
                    if p == mycontact # can't contact yourself
                        continue
                    end
                    if pop.contacts[mycontact] >= pop.touches[mycontact] # can't contact a person over his/her limit
                        continue
                    end
                    if mycontact in touchtrack # can't contact the same person twice
                        continue
                    end

                    quiet || println(p, " => ", mycontact)
                    touchtrack[n+=1] = mycontact
                    pop.contacts[p] += 1
                    pop.contacts[mycontact] += 1
                    touchcount += 1

                    if tries > maxtries
                        break
                    end
                    
                end
                
            end
        end
    end
    return touchcount
end    

# %% [markdown]
# ## Start at next cell

# %%
limit = 5
num = 100_000
num_unexp = 25000
tab = create_peeps(num)
mytouches= static_setup(tab, limit, num_unexp)
tab

# %%
countmap(tab.status)

# %%
@time cnt = run_touches!(tab, mytouches)
tab

# %%
cnt

# %%
countmap(tab.contacts)

# %%
countmap(tab.contacts[tab.status .== 2])

# %%
countmap(tab.contacts[tab.status .== 1])

# %%
countmap(tab.contacts[tab.status .== 3])

# %%
countmap(tab.contacts[tab.status .== 4])

# %% [markdown]
# ## Different ways to randomize the contacts
# - shuffle! the vector of indices
# - randperm: since the indices are basically 1:n: the challenge is that we have holes for excluded rows
# - rand: select from the index vector
# - sample: select n from the index vector

# %%
@btime randperm(500);

# %%
v = collect(1:500)
@btime shuffle!(v);

# %%
@btime for i in 1:100
    sample($v, 5);
end

# %%
@btime rand($v, 500);  # this is the same as shuffle--not in place so the result has to be allocated

# %%
@btime for i in 1:100
    rand($v,5)
end

# %%
collect(1:5:500) # starting index for 5 element ranges

# %%
# this is the winner, starting with shuffle!(vec)


@btime for i in 1:5:500
    (@view $v[i:i+4])[1]
end  

# %%
vec()
