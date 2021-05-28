# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
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
        contact_limit = zeros(Int,n),
        contacts = zeros(Int, n),
        status = zeros(Int, n)
        )  
end

# %%
function static_setup(pop, limit, n_unexp)
    for i in eachindex(pop)
        pop.contact_limit[i] = rand(1:limit)
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
end


# %%

function run_touches!(pop, alive_idx, limit)

    overall_count = 0

    # alive = findall(pop.status .!= 4) # the living can't contact the dead and visa versa
    overall_limit = sum(pop.contact_limit[alive_idx])
    shuffle!(alive_idx) # go through the indices in random order
    next_alive = 0
    n_alive = length(alive_idx)

    @inbounds for i in 1:n_alive   # p is the person reaching out to make contacts
        p = alive_idx[i]
        if pop.contacts[p] >= pop.contact_limit[p]
            continue                # this person is done making contacts--received their limit
        end
            
        while abs(overall_limit - overall_count) > 2    # made all the contacts we could
            next_alive = next_alive < n_alive ? next_alive + 1 : i + 1 # contacts below "i" are exhausted
            mycontact = alive_idx[next_alive]      # mycontact is the index of the person on the receiving end

            # rejected contacts
            if p == mycontact # can't contact yourself
                continue
            elseif pop.contacts[mycontact] >= pop.contact_limit[mycontact] # can't exceed contact's limit
                continue
            end  # no else needed-> just drop through

            pop.contacts[p] += 1
            pop.contacts[mycontact] += 1
            overall_count += 2  # each contact increases the total contacts by 2, 1 for each person

            if pop.contacts[p] >= pop.contact_limit[p]  # the contact initiator can't exceed his/her limit
                break                # why break instead of continue: we need a new 'p'
            end

        end   # while abs(...)
    end   # for i in 1:n_alive
    return overall_count
end    


# %% [markdown]
# ## Start at next cell

# %%
limit = 5
num = 100_000
num_unexp = 45_000
tab = create_peeps(num)
static_setup(tab, limit, num_unexp)
tab

# %%
alive = findall(tab.status .!= 4);

# %%
countmap(tab.status)

# %%
@time counts = run_touches!(tab, alive, limit)
tab

# %%
180 * .0034

# %%
(counts, sum(tab[alive].contacts))

# %%
countmap(tab.contacts[tab.status .== 2])

# %%
sum(values(ans))

# %%
countmap(tab.contacts[tab.status .== 1])

# %%
sum(values(ans))

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
#
# # %%
# @btime randperm(500);
#
# # %%
# v = collect(1:500)
# @btime shuffle!(v);
#
# # %%
# @btime for i in 1:100
#     sample($v, 5);
# end
#
# # %%
# @btime rand($v, 500);  # this is the same as shuffle--not in place so the result has to be allocated
#
# # %%
# @btime for i in 1:100
#     rand($v,5)
# end
#
# # %%
# collect(1:5:500) # starting index for 5 element ranges
#
# # %%
# # this is the winner, starting with shuffle!(vec)
# # well, not so much. when you run a lot more and save the view to get all of its contacts,
# #    there are lots'o' allocations
#
# @btime for i in 1:5:500
#     (@view $v[i:i+4])[1]
# end  
#
# # %%
# # bit of a failed experiment...  ...probably could be improved by taking only one sample from the view--
# # but then why do it?
#
# # %%
#
# function alt_touches!(pop, limit, quiet=true)
#     quiet || println("starting...")
#     overallcount = 0
#     contact_tracker = Vector{Int}(undef, limit)
#
#     alive = findall(pop.status .!= 4) # the living can't contact the dead and visa versa
#     shuffle!(alive) # go through the indices in random order
#     next_alive = 0
#     n_alive = length(alive)
#
#     sick = findall(pop.status .== 2)
#
#     for p in sick
#         contact_tracker[:] .= 0
#         cnt_contacts = 0  # index for contact_track and count of valid contacts for each p
#         lim_contacts = pop.contact_limit[p]
#             
#         for i in 1:limit:n_alive-limit 
#
#             guys = @view alive[i:i+limit-1]
#
#             for j in 1:limit
#                 mycontact = guys[j]
#
#                 # rejected contacts
#                 if p == mycontact # can't contact yourself
#                     continue
#                 elseif pop.contacts[mycontact] >= pop.contact_limit[mycontact] # can't exceed contact's limit
#                     continue
#                 # elseif pop.contacts[p] >= pop.contact_limit[p]
#                 #     continue
#                 elseif mycontact in contact_tracker # can't contact the same person twice
#                     continue
#                 end  # no else needed-> just drop through
#
#                 quiet || println(p, " => ", mycontact)
#                 cnt_contacts += 1  # increment number of valid contacts
#                 if cnt_contacts >= lim_contacts
#                     break  # end the while loop; get the next p in sick
#                 else
#                     contact_tracker[cnt_contacts] = mycontact
#                     pop.contacts[p] += 1
#                     pop.contacts[mycontact] += 1
#                     overallcount += 1
#                     # @assert pop.contacts[p] <= pop.contact_limit[p] "spreader $p exceeded contact limit at $next_alive"
#                     # @assert pop.contacts[mycontact] <= pop.contact_limit[mycontact] "contact $mycontact exceeded contact limit at $next_alive"
#                 end
#             end # for j
#         end # for i
#     end # for p in sick
# #     end
#     return overallcount
# end    
#
# # %%
