# TODO
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
using StatsPlots
using Debugger
using Dates

include("dec_tree.jl")
include("setup.jl")
include("sim.jl")
include("tracking.jl")
include("transition.jl")
include("spread.jl")
include("cases.jl")
include("johns_hopkins_data.jl")


export                  # functions for simulation
    run_a_sim,
    seed!,
    transition!,
    spread!,
    how_many_contacts!,
    how_many_touched,
    how_many_infected,
    isolate!,
    unisolate!,
    grab,
    input!,
    plus!,
    minus!,
    total!,
    Env,
    Spreadcase,
    sd_gen,
    seed_case_gen,
    case_setter,
    initialize_sim_env,
    r0_sim

export                  # functions for setup
    build_data,
    setup

export                  # functions for tracking
    reviewdays,
    showq,
    cumplot,
    newplot,
    dayplot,
    dayanimate2,
    review_history,
    make_series

export            # queues for tracking
    travelq,
    isolatedq,
    newstatq,
    spreadq,
    day2df,
    map2series

export                  # functions for decision trees
    setup_dt,
    read_dectree_file,
    create_node_dict,
    display_tree,
    sanity_test,
    get_the_probs

export                 # functions for accessing data from Johns Hopkins
    get_real_data,
    loc2df,
    read_actual

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
