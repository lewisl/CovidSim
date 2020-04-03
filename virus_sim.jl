#=
Virus simulation converted from Python to Julia
=#

using Printf
using Random

abstract type Virus end

    mutable struct SimpleVirus <: Virus
        max_birth_prob::Float64
        clear_prob::Float64
    end

    mutable struct ResistantVirus <: Virus
        max_birth_prob::Float64
        clear_prob::Float64
        resistances::Dict{String,Bool}
        mut_prob::Float64
    end


# methods for virus types

function does_clear(v::Virus)
      return (rand() < v.clear_prob ? true : false)
end    


function reproduce(v::SimpleVirus, pop_density)
    if rand() < v.max_birth_prob * (1 - pop_density)
        vnew = SimpleVirus(v.max_birth_prob, v.clear_prob)
        return vnew
    else
        return 0  # caller has to test for this!
    end
end


function is_resistant_to(v::ResistantVirus, drug)
    return get(v.resistances, drug, false)  # actual value or false as default value
end


function reproduce(v::ResistantVirus, pop_density, active_drugs)
    drnum = 0  # set before the loop because we need it after the loop ends
    for (drnum, dr) in enumerate(active_drugs)
        if !is_resistant_to(v, dr)
            drnum = 0
            break
        end
    end
    isresistant = drnum == length(active_drugs) ? true : false
    if isresistant && rand() < v.max_birth_prob * (1 - pop_density)
        restraits = copy(v.resistances) # don't change the source list!
        for dr in keys(restraits)
            if rand() < v.mut_prob
                restraits[dr] = !restraits[dr]
            end
        end
        vnew = ResistantVirus(v.max_birth_prob, v.clear_prob, restraits, v.mut_prob)
        return vnew
    else
        return 0  # caller has to test for this!
    end
end

# types for patient

abstract type Patient end

    mutable struct SimplePatient <: Patient
        viruses::Vector{SimpleVirus}
        max_pop::Float64
    end

    mutable struct TreatedPatient <: Patient
        viruses::Vector{ResistantVirus}
        max_pop::Float64
        drugs::Vector{String}
    end
    # constructor for assuming no drugs have been prescribed
    TreatedPatient(viruses, max_pop) = TreatedPatient(viruses,max_pop,[])

# methods for patient types
function add_prescription!(tp::TreatedPatient, newdrug)
    if !in(newdrug, tp.drugs)
        push!(tp.drugs, newdrug)
    end
end

function get_total_pop(p::Patient)
    return length(p.viruses)
end


function get_resistant_pop(tp::TreatedPatient, druglist) # viruses in patient resistant to list
    respop = 0
    for v in tp.viruses
        v_drugs = v.resistances
        fully_resistant = true
        for test_drug in druglist
            if !in(test_drug,keys(v_drugs)) # if test drug not in virus's list, then virus can't resist it
                fully_resistant = false
                break
            end
            if !v_drugs[test_drug]    # if resistance value is false
                fully_resistant = false
                break
            end 
        end
        if fully_resistant
            respop += 1
        end
    end
    return respop
end


function update!(sp::SimplePatient)  
    killviruses = SimpleVirus[]
    addviruses = SimpleVirus[]
    for v in sp.viruses
        if does_clear(v)
            push!(killviruses, v) # can't delete things in the collection we loop over
        end
    end
    for v in killviruses
        deleteat!(sp.viruses,findfirst(sp.viruses,v))
    end
    popdensity = length(sp.viruses) / sp.max_pop
    for v in sp.viruses
        childv = reproduce(v, popdensity)
        if isa(childv, SimpleVirus)  # make sure we got a child virus not 0
            push!(addviruses, childv)
        end
    end
    for a in addviruses
        push!(sp.viruses, a)
    end
    return length(sp.viruses)
end

function update!(tp::TreatedPatient)  
    killviruses = ResistantVirus[]
    addviruses = ResistantVirus[]
    for v in tp.viruses
        if does_clear(v)
            push!(killviruses, v) # can't delete things in the collection we loop over
        end
    end
    for v in killviruses
        deleteat!(tp.viruses,findfirst(isequal(v),tp.viruses))
    end
    popdensity = length(tp.viruses) / tp.max_pop
    for v in tp.viruses
        childv = reproduce(v, popdensity, tp.drugs)
        if isa(childv, ResistantVirus) # make sure we got a child virus not 0
            push!(addviruses, childv)
        end
    end
    for a in addviruses
        push!(tp.viruses, a)
    end
    return length(tp.viruses)
end


# # ======================== simulation runs ========================================

using PyPlot
Random.seed!(1234)

"""
    For each of numTrials trials, instantiates a patient, runs a simulation for
    150 timesteps, adds guttagonol, and runs the simulation for an additional
    150 timesteps.  At the end plots the average virus population size
    (for both the total virus population and the guttagonol-resistant virus
    population) as a function of time.

    numViruses: number of ResistantVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: maximum clearance probability (a float between 0-1)
    resistances: a dictionary of drugs that each ResistantVirus is resistant to
                 (e.g., {'guttagonol': False})
    mutProb: mutation probability for each ResistantVirus particle
             (a float between 0-1). 
    numTrials: number of simulation runs to execute (an integer)
    
"""
function simulationWithDrug(numViruses::Int64, maxPop::Int64, maxBirthProb::Float64, 
                            clearProb::Float64, resistances::Dict{String,Bool}, 
                            mutProb::Float64, numTrials::Int64)
    # example: simulationWithDrug(1, 20, 1.0, 0.0, Dict("guttagonol"=> true), 1.0, 5)

    timesteps = 150
    cycles = 2
    viruses = [ResistantVirus(maxBirthProb,clearProb, resistances, mutProb) for i in 1:numViruses]
    vpops = [0 for i in 1:(cycles * timesteps)]
    resistpops = [0 for i in 1:(cycles * timesteps)]

    for i in 1:numTrials
        pat = TreatedPatient(viruses, maxPop)
        for t in 1:timesteps
            update!(pat)
            vpops[t] += get_total_pop(pat)
            resistpops[t] += get_resistant_pop(pat, ["guttagonol"])
        end
        add_prescription!(pat, "guttagonol")
        for t in timesteps+1:cycles*timesteps
            update!(pat)
            vpops[t] += get_total_pop(pat)
            resistpops[t] += get_resistant_pop(pat, ["guttagonol"])
        end
    end

    numTrials = float(numTrials)
    vpops = [i / numTrials for i in vpops]
    resistpops = [i / numTrials for i in resistpops]
    figure()
    plot(vpops)
    plot(resistpops)
    xlabel("Simulation Time Steps")
    ylabel("Avg. Virus Population")
    title("Simulated Avg. Virus Population")
    legend(("Virus Pop","Resistant Pop"), loc=4)
    # show()
end

# # ================= simulationDelayedTreatment ========================================

"""
    Runs numTrials simulations to show the relationship between delayed
    treatment and patient outcome using a histogram.

    Histograms of final total virus populations are displayed for delays of 300,
    150, 75, 0 timesteps (followed by an additional 150 timesteps of
    simulation).

    numTrials: number of simulation runs to execute (an integer)
"""
function simulationDelayedTreatment(numTrials::Int64)

    numViruses = 100  # 100
    maxPop = 1000  # 1000
    maxBirthProb = 0.1  # 0.1
    clearProb = 0.05  # 0.05
    input_resistances = Dict("guttagonol"=> false)   # false
    mutProb = 0.005
    timesteps_after = 150
    divisor = float(numTrials) # to compute averages across numtrials

    # set different experimental conditions before beginning trials
    for (cond_cnt,timesteps_before) in enumerate([300, 150, 75, 0])  # [300, 150, 75, 0]
        # static across all trials
        totalpops = [0 for i in 1:(timesteps_before + timesteps_after)]
        resistpops = [0 for i in 1:(timesteps_before + timesteps_after)]
        totalpop_pertrial = []
        for i in 1:numTrials
            input_viruses = [ResistantVirus(maxBirthProb,clearProb, input_resistances, mutProb) 
                            for i in 1:numViruses]
            pat = TreatedPatient(input_viruses, maxPop)
            for t in 1:timesteps_before
                update!(pat)
                totalpops[t] += get_total_pop(pat)
                resistpops[t] += get_resistant_pop(pat, ["guttagonol"])
            end
            # if rand() < .2 # do nothing for 20% uncooperative patients (trials)
            # else
            #     add_prescription!(pat, "guttagonol")
            # end
            add_prescription!(pat, "guttagonol")
            for t in (timesteps_before + 1):(timesteps_before + timesteps_after)
                update!(pat)
                totalpops[t] += get_total_pop(pat)
                resistpops[t] += get_resistant_pop(pat, ["guttagonol"])
            end
            push!(totalpop_pertrial, get_total_pop(pat))
        end
        totalpops = [i / divisor for i in totalpops]
        resistpops = [i / divisor for i in resistpops]

        figure(1)
        subplot(2,2,cond_cnt)
        title(@sprintf("Delay %0i ", timesteps_before))
        plot(totalpops)
        plot(resistpops)
        # xlabel("Time Steps")
        ylabel("Avg. Virus Pop.")
        # legend(("Virus Pop","Resistant Pop"), loc=4)

        figure(2)
        subplot(2,2,cond_cnt)
        title(@sprintf("Delay %0i ", timesteps_before))
        plt[:hist](totalpop_pertrial,bins=50)
    end
    figure(1)
    suptitle("Simulated Virus Pop. Over Time")

    figure(2)
    suptitle("Total Virus Pop. Dist.")

    print("\nEnd of program.  Press enter to close charts.> "); chomp(readline())
    close("all")
end

# # ================= simulationTwoDrugDelayedTreatment ====================================

"""
    Runs numTrials simulations to show the relationship between delayed
    treatment with a drug cocktail and patient outcome using a histogram.

    Histograms of final total virus populations are displayed for delays of 300,
    150, 75, 0 timesteps (preceded by and followed by an additional 150 timesteps of
    simulation).

    numTrials: number of simulation runs to execute (an integer)
"""
function simulationTwoDrugDelayedTreatment(numTrials::Int64)
    # simulationWithDrug(1, 20, 1.0, 0.0, {"guttagonol": True}, 1.0, 5)
    # def simulationWithDrug(numViruses, maxPop, maxBirthProb, clearProb, input_resistances,
    #                    mutProb, numTrials):

    numViruses = 100  # 100
    maxPop = 1000  # 1000
    maxBirthProb = 0.1  # 0.1
    clearProb = 0.05  # 0.05
    input_resistances = Dict("guttagonol"=> false, "grimpex"=> false)   # false
    mutProb = 0.005  # 0.005
    timesteps_one = 150
    # timesteps_two set below to vary
    timesteps_three = 150
    divisor = float(numTrials) # to compute averages across numtrials

    # set different experimental conditions before beginning trials
    for (cond_cnt,timesteps_two) in enumerate([300, 150, 75, 0])  # [300, 150, 75, 0]
        # static across all trials
        totalpops = [0 for i in 1:(timesteps_one + timesteps_two + timesteps_three)]
        resistpops = [0 for i in 1:(timesteps_one + timesteps_two + timesteps_three)]
        totalpop_pertrial = []
        for i in 1:numTrials
            input_viruses = [ResistantVirus(maxBirthProb,clearProb, input_resistances, mutProb) 
                            for i in 1:numViruses]
            pat = TreatedPatient(input_viruses, maxPop)
            # first phase
            for t in 1:timesteps_one
                update!(pat)
                totalpops[t] += get_total_pop(pat)
                resistpops[t] += get_resistant_pop(pat, ["guttagonol", "grimpex"])
            end
            # second phase
            add_prescription!(pat, "guttagonol")
            for t in (timesteps_one + 1):(timesteps_one + timesteps_two)
                update!(pat)
                totalpops[t] += get_total_pop(pat)
                resistpops[t] += get_resistant_pop(pat, ["guttagonol", "grimpex"])
            end
            # third phase
            add_prescription!(pat, "grimpex")
            for t in (timesteps_one + timesteps_two + 1):(timesteps_one + timesteps_two + timesteps_three)
                update!(pat)
                totalpops[t] += get_total_pop(pat)
                resistpops[t] += get_resistant_pop(pat, ["guttagonol", "grimpex"])
            end
            push!(totalpop_pertrial, get_total_pop(pat))
        end
        totalpops = [i / divisor for i in totalpops]
        resistpops = [i / divisor for i in resistpops]

        figure(1)
        subplot(2,2,cond_cnt)
        title(@sprintf("Delay %0i ", timesteps_two))
        plot(totalpops)
        plot(resistpops)
        # xlabel("Time Steps")
        ylabel("Avg. Virus Pop.")
        # legend(("Virus Pop","Resistant Pop"), loc=4)

        figure(2)
        subplot(2,2,cond_cnt)
        title(@sprintf("Delay %0i ", timesteps_two))
        plt[:hist](totalpop_pertrial,bins=50)
    end
    figure(1)
    suptitle("Simulated Virus Pop. Over Time")

    figure(2)
    suptitle("Total Virus Pop. Dist.")

    print("\nEnd of program.  Press enter to close charts.> "); chomp(readline())
    close("all")
end


# # ============================= some tests ===========================

function test_resist()
    virus1 = ResistantVirus(1.0, 0.0, Dict("drug1"=> true), 0.0)
    virus2 = ResistantVirus(1.0, 0.0, Dict("drug1"=> false, "drug2"=> true), 0.0)
    virus3 = ResistantVirus(1.0, 0.0, Dict("drug1"=> true, "drug2"=> true), 0.0)
    patient = TreatedPatient([virus1, virus2, virus3], 100)
    @assert get_total_pop(patient) == 3
    @assert get_resistant_pop(patient,["drug1"]) == 2
    @assert get_resistant_pop(patient,["drug2"]) == 2
    @assert get_resistant_pop(patient,["drug1","drug2"]) == 1
    @assert get_resistant_pop(patient,["drug3"]) == 0
    @assert get_resistant_pop(patient,["drug1", "drug3"]) == 0
    @assert get_resistant_pop(patient,["drug1","drug2", "drug3"]) == 0
end

function test_siminputs()
    virus1 = ResistantVirus(1.0, 0.0, Dict("guttagonol"=> true), 1.0)
    patient = TreatedPatient([virus1], 20)
    @assert get_resistant_pop(patient, ["guttagonol"]) == 1
end