# TODO
    # sanity check for spread mean outcomes
    # transition loop should auto-recognize where the branches should run
    # more info
        # get fatality rate by age and co-morbidity CDC, Italian NIH
        # by agegroup, hospitalization %, ICU admission %, fatality %
        # UW virology, expansion of deaths by state on log chart
    # probably need to raise effective death rate for people in 60-80 and 80+  per Vox article

# Done


module CovidSim

using DelimitedFiles
using DataStructures
using DataFrames
using Random
using Distributions
using StatsBase
using Printf
using Plots

include("dec_tree.jl")
include("setup.jl")
include("sim.jl")
include("tracking.jl")
include("isolation.jl")
include("spread.jl")


export                  # functions for simulation
    run_a_sim,
    seed!,
    transition!,
    spread!,
    how_many_contacts,
    how_many_touched,
    how_many_infected,
    isolate!,
    _isolate!,
    unisolate!,
    _unisolate!,
    cumplot,
    newplot,
    grab,
    input!,
    plus!,
    minus!,
    total!,
    Env


export                  # functions for setup
    build_data,
    setup

export                  # functions for tracking
    reviewdays,
    showq

export            # queues for tracking
    travelq,
    isolatedq,
    newstatq,
    dayq

export                  # functions for decision trees
    setup_dt,
    read_dectree_file,
    create_node_dict,
    display_tree,
    sanity_test,
    get_the_probs

export                  # control constants
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



end # module CovidSim
