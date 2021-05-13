# TODO
    # get rid of all uses of integer literals for sickday, agegrps, conditions
        # in dectrees, for instance
    # more info
        # get fatality rate by age and co-morbidity CDC, Italian NIH
        # by agegroup, hospitalization %, ICU admission %, fatality %
        # UW virology, expansion of deaths by state on log chart
    # fix all the travel functions to latest APIs
    # put in an inflection measure

# TODO for individual level model
    # rewrite r0 simulator
    # add transq to ilm and test
    # vaccination with 3 vaccines / 1 or 2 shots
    # extend to one year
    # change meaningful constants to _nil_, _severe_, _dead_ etc.???

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
using Dates
using YAML
using TypedTables
using OffsetArrays
using SplitApplyCombine


# order matters for these includes!
include("../shared-src/dec_tree.jl")
include("setup.jl")
include("../shared-src/tracking.jl")
include("cases.jl")
include("test_and_trace.jl")
include("transition.jl")
include("spread.jl")
include("sim.jl")
include("../shared-src/johns_hopkins_data.jl")

# functions for simulation
export    
    setup,              
    run_a_sim,
    day_ctr,
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
    map2series

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
    sickdays,
    sickdaylim

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
    notsick,
    nil,
    mild,
    sick,
    severe,
    totinfected,
    statuses,
    conditions,
    allconds,
    condnames,
    infectious_cases,
    transition_cases,
    map2series,
    series_colnames,
    age0_19,
    age20_39,
    age40_59, 
    age60_79, 
    age80_up, 
    agegrps,
    n_agegrps,
    recv_risk,
    totalcol


###########################################################################
# module global constants (except in Julia things aren't really constant!)
###########################################################################

# datatype constants
const T_int = Ref(Int64)  # this creates a reference type accessed or modified in functions as T_int[]


"""
- use incr!(day_ctr, :day) for day of the simulation:  creates and adds 1
- use reset!(day_ctr, :day) to remove :day and return its current value, set it to 0
- use day_ctr[:day] to return current value of day
"""
const day_ctr = counter(Symbol) # from package DataStructures


################################################################
# constants for data structure indices
################################################################

# control constants
const age_dist = [0.251, 0.271, 0.255, 0.184, 0.039]
const sickdaylim = 25
const sickdays = 1:sickdaylim   # rows

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

# status
const unexposed         = 1  
const infectious        = 2
const recovered         = 3
const dead              = 4

# agegrp 
const age0_19           = 1 
const age20_39          = 2 
const age40_59          = 3 
const age60_79          = 4 
const age80_up          = 5 

# condition
const notsick           = 0
const nil               = 5
const mild              = 6
const sick              = 7
const severe            = 8

const statuses          = [unexposed,infectious,recovered,dead]
const conditions        = [nil,mild,sick,severe]
const all_conds         = [unexposed,infectious,recovered,dead,nil,mild,sick,severe]
const infectious_cases  = [nil, mild, sick, severe]
const transition_cases  = [recovered, nil, mild, sick, severe, dead]
const agegrps           = [age0_19,age20_39,age40_59,age60_79,age80_up]
const n_agegrps         = length(agegrps)
const condnames         = Dict(0=>"notsick", 1=>"unexposed", 2=>"infectious", 3=>"recovered", 4=>"dead",
                                5=>"nil", 6=>"mild", 7=>"sick", 8=>"severe", 9=>"totinfected")

const totinfected       = 9
const travelers         = 10
const isolated          = 11

# columns of history series: first 5 cols are agegrps, 6th is total
const map2series = (unexposed=1:6, infectious=7:12, recovered=13:18, dead=19:24, 
                    nil=25:30, mild=31:36, sick=37:42, severe=43:48, totinfected=49:54)
const totalcol = 6


# traveling constants
const travprobs = [1.0, 2.0, 3.0, 3.0, 0.4] # by age group

end # module CovidSim
