# TODO
    # get rid of all uses of integer literals for lag, agegrps, conditions
        # in dectrees, for instance
    # more info
        # get fatality rate by age and co-morbidity CDC, Italian NIH
        # by agegroup, hospitalization %, ICU admission %, fatality %
        # UW virology, expansion of deaths by state on log chart
    # fix all the travel functions to latest APIs
    # put in an inflection measure

# TODO for individual level model
    # get rid of indirection for lowlevel population updates
    # add transq to ilm and test

# Done
    # do unquarantine for ilm



__precompile__(true)

module CovidSim_ilm

# required
using DelimitedFiles
using DataStructures
using DataFrames
using CSV
using Random
using Distributions
using StatsBase
using Printf
using Plots
using PlotThemes
using StatsPlots
using Debugger
using Dates
using YAML
using TypedTables
using OffsetArrays
using SplitApplyCombine


# order matters for these includes!
include("dec_tree.jl")
include("setup.jl")
include("tracking.jl")
include("cases.jl")
include("test_and_trace.jl")
include("transition.jl")
include("spread.jl")
include("sim.jl")
include("../shared/johns_hopkins_data.jl")

# functions for simulation
export    
    setup,              
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
    sumit,
    SimEnv,
    initialize_sim_env,
    r0_sim,
    set_by_level

# functions for cases
export
    test_and_trace,     
    Spreadcase,
    sd_gen,
    seed_case_gen,
    t_n_t_case_gen,
    case_setter,
    bayes,
    tntq,
    shifter

# functions for setup
export                  
    build_data,
    setup

# functions for tracking
export                  
    reviewdays,
    showq,
    cumplot,
    newplot,
    dayplot,
    dayanimate2,
    review_history,
    make_series,
    virus_outcome

# queues and caches (variables) for tracking
export            
    travelq,
    spreadq,
    transq,
    day2df,
    map2series,
    ctr,
    sim_stash

# functions for decision trees
export                  
    setup_dt,
    read_dectree_file,
    create_node_dict,
    display_tree,
    sanity_test,
    get_the_probs

# functions for accessing data from Johns Hopkins
export                 
    get_real_data,
    loc2df,
    read_actual

# control constants
export                  
    age_dist,
    lags,
    laglim

# constants for geo data
export      
    fips,
    state,
    size_cat,
    popsize,
    major,
    large,
    medium,
    small,
    smaller,
    rural

# constants for indices to population matrix
export              
    unexposed,
    infectious,
    recovered,
    dead,
    nil,
    mild,
    sick,
    severe,
    totinfected,
    conditions,
    condnames,
    infectious_cases,
    transition_cases,
    map2series,
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
    n_agegrps,
    recv_risk_by_age,
    col_status,
    col_agegrp,
    col_cond,
    col_lag,
    col_cluster,
    col_recov_day,
    col_dead_day,
    col_vax,
    col_vax_day,
    col_test,
    col_test_day,
    col_quar,
    col_quar_day,
    totalcol


###########################################################################
# module global constants (except in Julia things aren't really constant!)
###########################################################################

# datatype constants
const T_int = Ref(Int64)  # this creates a reference type accessed or modified in functions as T_int[]


################################################################
# constants for data structure indices
################################################################

# control constants
const age_dist = [0.251, 0.271,   0.255,   0.184,   0.039]
const laglim = 25
const lags = 1:laglim   # rows

# geo data: fips,county,city,state,sizecat,pop,density
const fips = 1
const county = 2
const city = 3
const state = 4
const sizecat = 5
const popsize = 6
const density = 7
const anchor = 8
const restrict = 9
const density_fac = 10

# population centers sizecats
const major = 1  # 20
const large = 2  # 50
const medium = 3
const small = 4
const smaller = 5
const rural = 6

# stats series/dataframe columns
const unexposed         = 1  # note that 1:8 are also values for status and cond
const infectious        = 2
const recovered         = 3
const dead              = 4
const nil               = 5
const mild              = 6
const sick              = 7
const severe            = 8
const totinfected       = 9
const travelers         = 10
const isolated          = 11

# columns of population matrix not used with TypedTables, but these are still correct
const col_status       = 1
const col_agegrp       = 2
const col_cond         = 3
const col_lag          = 4
const col_cluster      = 5
const col_recov_day    = 6
const col_dead_day     = 7
const col_susceptible  = 8
const col_vax          = 9
const col_vax_day      = 10
const col_test         = 11
const col_test_day     = 12
const col_quar         = 13
const col_quar_day     = 14

# columns of history series: first 5 cols are agegrps, 6th is total
const map2series = (unexposed=1:6, infectious=7:12, recovered=13:18, dead=19:24, 
                    nil=25:30, mild=31:36, sick=37:42, severe=43:48, totinfected=49:54)
const totalcol = 6

const conditions = [unexposed, infectious, recovered, dead, nil, mild, sick, severe]
const n_conditions = length(conditions)
const condnames = Dict(1=>"unexposed", 2=>"infectious", 3=>"recovered", 4=>"dead",
                       5=>"nil", 6=>"mild", 7=>"sick", 8=>"severe", 9=>"totinfected")
const infectious_cases = [nil, mild, sick, severe]
const transition_cases = [recovered, nil, mild, sick, severe, dead]

# transition_prob_rows
const to_recovered = 1
const to_nil = 2
const to_mild = 3
const to_sick = 4
const to_severe = 5
const to_dead = 6

# agegrp channels at dimension 3
const a1 = 1 # 0-19
const a2 = 2 # 20-39
const a3 = 3 # 40-59
const a4 = 4 # 60-79
const a5 = 5 # 80+
const agegrps = 1:5
const n_agegrps = length(agegrps)


# struct to hold current status indices
mutable struct Current_idx
   data::Array{Int,1}
   last::Int
   Current_idx(n) = new(zeros(Int,n),0)
end


# traveling constants
const travprobs = [1.0, 2.0, 3.0, 3.0, 0.4] # by age group


end # module CovidSim
