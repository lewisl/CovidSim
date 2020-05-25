

#function daystep()  # one day of simulation

    # DONE loop across days of the simulation
        # function hooks for:...
        # DONE seeding events--especially inbound international travel
        # outbreaks--when we don't have pockets implemented
        # rules--that restrict movement and affect touches and infection contacts  DEFER


    # need rules for how many people travel across locales; then rules for travel restrictions
         # DONE = travelin!: some people travel in; distribution of them are infectious
        # pull item from queue:  add to destination: agegrp/condition, agegrp/outcome
        #                        remove from source: agegrp/condition, agegrp/outcome
        # TODO add to log

    # some people will isolate: need to come up with rules for isolation, leaks from isolation
        #   DONE = isolate_queue! put people into isolation by agegroup, locale, condition, lag
                    # and isolate_move!=>remove from open, unisolate_queue!, unisolate_move!
                    # add isolation people to queue
                    # unisolate:  remove people from isolatedmx and put them back into openmx
        # DONE transition people in isolation across conditions
        

    # DONE = transition! distribute residents across all lags and conditions
        # openmx
        # isolatedmx
        # loop across all locales to transition all
            # start from last lag, by age group:
                    # severe distribute to sick, severe or die remain at lag 18
                    # sick distribute to severe, mild remain at lag 18
                    # mild distribute to sick, recovered (no nil at this point) at lag 18
                    # nil distribute to nil, mild, sick, severe at lags 2-6
            # DONE for earlier lags:  - part of decision trees for transition
                    # severe distribute to sick, severe or die => advance one lag
                    # sick distribute to sick, severe, mild => advance one lag
                    # mild distribute to sick, recovered => advance one lag
                    # nil distribute to nil, mild, recovered => advance one lag

            # go dead => remove from previous condition; remove from resident, remove from totalpop:  make sure we can count # dead today (or more than yesterday)
            # DONE go recovered => remove from previous condition; remove from exposed? remove from infectious

            # TODO log new dead, new recovered  to change state



    # DONE:  spread from infectious people to unexposed|recovered people(not isolated)
            #DONE contact distribution of residents across age&condition&lags
            # DONE: handle partial (treating as full) immunity of recovered
            # DONE = spread!  touchers -> num_touched -> split by age/condition
            # DONE in spread!  of those contacted, a distribution become exposed & nil at lag 0

            # DONE log new infectious


    # DONE = travelout! some residents travel out
        # choose distribution of people traveling by age and condition
        # add to queue

    # TODO summarize current state values for unexposed, infectious, dead, recovered
        # add to series using change_state, summarize function

#end

# TODO Build locale data in a dataframe
#   cols: id, state, locale, pop, density, travelprobs, gender split, age dist (maybe)
#
#
#
#
#
