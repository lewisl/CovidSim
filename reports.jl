# <md>
# Let's load the simulation code and do a first run to get things started
#
# <codecell>
using CovidSim
seed_1_24 = seed_case_gen(1, [0,12,12,0,0], 3, nil, agegrps)
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 3, nil, agegrps)
alldict, dseries, starting_unexposed = run_a_sim("geo2data.csv",180,11, dtfilename="dec_tree_all_25.csv",
       silent=true,spreadcases=[],
       runcases=[seed_1_6]);
geo = alldict["geo"];
cumplot(dseries,11, geo=geo)
# <md>
# The plot we see is for Phoenix (actually, all of Maricopa county). With no social distancing
# or isolation in effect, 58% of all of the people have contracted Coronavirus and, of these,
# about 1% have died. We started this with 6 people who arrived in Phoenix who were asymptomatic, but
# on the 4th day of their COVID-19 infection.  These 6 people spread the disease to over half of
# Phoenix. The graph shows 180 days which we can assume starts around February 1, 2020--so this graph
# ends on just about August 1.
#
# - The blue line plots the daily number of people who have not been exposed.
# - The orange line shows the *active*  cases for those who are infected. This
# is not what you see in most plots.  This plot shows those who are currently infected,
# which means the number of people who have died or recovered at a given point are not shown.
# What you are used to seeing in most reports is *total infected*, which is the total number
# of people who have ever been infected.
# - The green line shows the cumulative total number of people who have recovered as of each day.
# As the number of people who have recovered increases, there are fewer people who can become infected.
# This is why the orange line goes down and the infection "burns itself out."
# - This comes at the cost of the 27,159 people who have died.
#
# This depiction is dramatically more severe than anything we have observed anywhere in the world because...
#   ...everywhere people have isolated themselves or practiced social distancing, which have reduced the number
# of people who have gotten the disease. Without isolation, social distancing, or other methods
# (such as wearing *really* effective masks and replacing or cleaning them frequently) of preventing
# transmission of the virus, something like this is what would happen.
