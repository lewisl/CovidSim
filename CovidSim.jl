# TODO
    # make sure contacted is less than accessible
    # more info
        # get fatality rate by age and co-morbidity CDC, Italian NIH
        # by agegroup, hospitalization %, ICU admission %, fatality %
        # UW virology, expansion of deaths by state on log chart
    # probably need to raise effective death rate for people in 60-80 and 80+  per Vox article
    # need to reduce effective infection rates
    # look for speed improvement: preallocation


module CovidSim

export                  # functions
    run_a_sim,
    seed!,
    transition!,
    spread!,
    how_many_contacted,
    how_many_touched,
    how_many_infected,
    cumplot,
    newplot,
    grab,
    input!,
    plus!,
    minus!,
    total!,
    reviewbugs,
    showq

export                  # control constants
    bugq,
    age_dist,
    lags


export      # constants for geo data
    id,
    state,
    size_cat,
    popsize,
    major,
    large,
    medium,
    small,
    smaller,
    rural

export              # constants for indices to data tables
    unexposed,
    infectious,
    recovered,
    dead,
    nil,
    mild,
    sick,
    severe,
    conditions,
    condnames,
    infectious_cases,
    transition_cases,
    mapcond2tran,
    series_colnames,
    to_recovered,
    to_nil,
    to_mild,
    to_severe,
    to_sick,
    to_dead,
    a1,
    a2,
    a3,
    a4,
    a5,
    agegrps,
    ages,
    recv_risk_by_age

export            # queues for tracking
    travelq,   
    isolatedq,
    newstatq
    


using DelimitedFiles
using DataStructures
using DataFrames
using Random
using Distributions
using StatsBase
using Printf
using PyPlot
import Seaborn

include("dec_tree.jl")
include("setup.jl")
include("sim.jl")



end # module CovidSim