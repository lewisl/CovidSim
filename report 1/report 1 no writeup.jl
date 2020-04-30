# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.4.0
#     language: julia
#     name: julia-1.4
# ---

# %% hide_input=true
using CovidSim

# %% hide_input=true
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 5, nil, agegrps)
seed_6_12 = seed_case_gen(8, [0,6,6,0,0], 5, nil, agegrps)

# %% hide_input=true
alldict, evn, series = run_a_sim(180,11, silent=true,
        spreadcases=[], 
        runcases=[seed_1_6, seed_6_12]);
geo = alldict["geo"];

# %% hide_input=true slideshow={"slide_type": "slide"}
cumplot(series,11,geo=geo)

# %% hide_input=true slideshow={"slide_type": "slide"}
cumplot(series, 11, [infectious, dead], geo=geo)

# %% hide_input=true slideshow={"slide_type": "slide"}
dayplot(spreadq,[:spreaders, :contacts, :touched])

# %% slideshow={"slide_type": "slide"}
dayplot(spreadq)

# %% [markdown] slideshow={"slide_type": "slide"}
# #### Social Distancing

# %% [markdown]
# Here is a case when social distancing starts on day 60 with 75% compliance and is rather strong in terms of contacts and touches.

# %%
str_60 = sd_gen(start=60, comply=.75, cf=(.3,1.2), tf=(.18,.4));
alldict, env, series = run_a_sim(180,11, silent=true,
    spreadcases=[str_60], 
    runcases=[seed_1_6, seed_6_12]);

# %% slideshow={"slide_type": "slide"}
cumplot(series, 11, geo=geo)

# %% slideshow={"slide_type": "slide"}
cumplot(series, 11, [infectious, dead],geo=geo)

# %%
76670/2605750

# %%
str_50 = sd_gen(start=50, comply=.75, cf=(.3,1.2), tf=(.18,.4));
alldict, env, series = run_a_sim(180, 11, silent=true,
    spreadcases=[str_50],
    runcases=[seed_1_6, seed_6_12]);

# %% slideshow={"slide_type": "slide"}
cumplot(series, 11, geo=geo)

# %% [markdown]
# We cut not quite another 30% of infections compared to starting isolation 10 days later.  With exponential growth, 
# when you start mitigation practices makes a big difference. 
#
# Looking at this the authorities could say, "We are in great shape; we have this virus licked." 
#
# Let's see what happens if we open up to where we were before after 30 days of social distancing.

# %%
open_all = sd_gen(start=80, comply=0.0, cf=(.2,1.8), tf=(.18,.62)); 
alldict, env, series = run_a_sim(180,10,silent=true,
    spreadcases=[str_50,open_all],
    runcases=[seed_1_6, seed_6_12]);  

# %% slideshow={"slide_type": "slide"}
cumplot(series, 10, geo=geo)

# %% [markdown] slideshow={"slide_type": "slide"}
# It's kind of like nothing happened because of social distancing. We stopped the rise early; then we started the rise.
# We ended up in nearly the same place, except everything started later and the simulation ends before
# the infection curve burns itself out.
#
# |              | Infections | Deaths|
# |--------------|---------|----------|
# |No restrictions | 2,605,750 | 27,199 |
# | Restrict. on/off | 2,522,693 | 23,244 |
#
# Because Maricopa County has a large population of over 4 million, 18 people coming in with infections takes a while to get the exponential growth curve "rolling".  What would be more typical is a double bump.  Social distancing starts after significant growth is observed; then it brings the curve down; the people cheer; social distancing ends. 
#
# Let's see what happens in a smaller place like Omaha, Nebraska and we won't open up totally, which does get us back 
# to what looks like no social distancing at all.

# %%
close = sd_gen(start=50, comply=.75, cf=(.3,1.2), tf=(.18,.4)); open = sd_gen(start=90, comply=.75, cf=(.2,1.5), tf=(.18,.5));
alldict, env, series = run_a_sim(180,10, silent=true, 
    spreadcases=[close, open], 
    runcases=[seed_1_6, seed_6_12]);

# %% slideshow={"slide_type": "slide"}
cumplot(series, 10,[infectious, dead], geo=geo)

# %%
alldict, env, series = run_a_sim(180, 10, silent=true,
    spreadcases=[],
    runcases=[seed_1_6, seed_6_12]);

# %% slideshow={"slide_type": "slide"}
cumplot(series, 10,[infectious, dead, recovered], geo=geo)

# %%
