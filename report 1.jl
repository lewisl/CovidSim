using Plots



#=      pauser
        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end
=#

        






function run_report_1()

    seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 5, nil, agegrps)
    seed_6_12 = seed_case_gen(8, [0,6,6,0,0], 5, nil, agegrps)

    alldict, dseries, starting_unexposed = run_a_sim("geo2data.csv",180,11,
           dtfilename="dec_tree_all_25.csv", 
           silent=true,spreadcases=[],
           runcases=[seed_1_6, seed_6_12]);
    geo = alldict["geo"];

    print("\n\n**********\n")
    println("Phoenix with no social distancing or isolation...")
    print("**********\n\n")

    cumplot(dseries,11,geo=geo)

        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end

    print("\n\n**********\n")
    println("Let's look more closely at the infectious and the dead...")
    print("**********\n\n")

    cumplot(dseries,11,geo=geo, [:Infectious, :Dead])

        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end

    print("\n\n**********\n")
    println("How does the disease actually spread?")
    print("**********\n\n")

    dayplot(spreadq,[:spreaders, :contacts, :touched])

        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end

    dayplot(spreadq)

        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end


    print("\n\n**********\n")
    println("Let's do some strong social distancing with 75% compliance, starting on day 60...")
    print("**********\n\n")

    str_60 = sd_gen(start=60, comply=.75, cf=(.3,1.2), tf=(.18,.4));
    alldict, dseries, starting_unexposed = run_a_sim("geo2data.csv",180,11,
           dtfilename="dec_tree_all_25.csv", 
           silent=true,spreadcases=[str_60],
           runcases=[seed_1_6, seed_6_12]);
    cumplot(dseries, 11, geo=geo)

        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end

    print("\n\n**********\n")
    println("Focus on infectious and dead...")
    print("**********\n\n")

    cumplot(dseries, 11, [:Infectious, :Dead],geo=geo)

        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end

    76670/2605750


    print("\n\n**********\n")
    println("Let's start social distancing 10 days sooner...")
    print("**********\n\n")

    str_50 = sd_gen(start=50, comply=.75, cf=(.3,1.2), tf=(.18,.4));
    alldict, dseries, starting_unexposed = run_a_sim("geo2data.csv",180,11,
           dtfilename="dec_tree_all_25.csv", 
           silent=true,spreadcases=[str_50],
           runcases=[seed_1_6, seed_6_12]);
    cumplot(dseries, 11, geo=geo)

        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end

    print("\n\n**********\n")
    println("We're on top of it; Let's open up and party...")
    print("**********\n\n")

    open_80 = sd_gen(start=80, comply=1.0, cf=(.2,1.8), tf=(.18,.62));
    open_all = sd_gen(start=80, comply=0.0, cf=(.2,1.8), tf=(.18,.62)); # 0% compliance is a signal to end social distancing
    alldict, dseries, starting_unexposed = run_a_sim("geo2data.csv",180,11,
           dtfilename="dec_tree_all_25.csv", 
           silent=true,spreadcases=[str_50,open_all],  # strong social distancing, then open
           runcases=[seed_1_6, seed_6_12]);
    cumplot(dseries, 11, geo=geo)

        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end

    print("\n\n**********\n")
    println("Phoenix is too big.  Let's go to Omaha...")
    print("**********\n\n")

    close = sd_gen(start=50, comply=.75, cf=(.3,1.2), tf=(.18,.4));
    open = sd_gen(start=90, comply=.75, cf=(.2,1.5), tf=(.18,.5
            ));
    alldict, dseries, starting_unexposed = run_a_sim("geo2data.csv",180,14,
           dtfilename="dec_tree_all_25.csv", 
           silent=true,spreadcases=[close, open],
           runcases=[seed_1_6, seed_6_12]);
    cumplot(dseries, 14,[:Infectious, :Dead, :Recovered], geo=geo)

        print("press enter to continue or enter q to quit> "); resp = chomp(readline())
            if resp == "q"
                return
            end

    print("\n\n**********\n")
    println("Here's Omaha with no social distancing or isolation...")
    print("**********\n\n")

    alldict, dseries, starting_unexposed = run_a_sim("geo2data.csv",180,14,
           dtfilename="dec_tree_all_25.csv", 
           silent=true,spreadcases=[],
           runcases=[seed_1_6, seed_6_12]);
    cumplot(dseries, 14,[:Infectious, :Dead, :Recovered], geo=geo)
end
