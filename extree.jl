#=
this is the sample:
(1, 1)
   CovidSim.Branch(5, 5, 0.2, (2, 1), "nil", "nil")
   CovidSim.Branch(5, 6, 0.65, (2, 2), "nil", "mild")
   CovidSim.Branch(5, 7, 0.15, (2, 3), "nil", "sick")
(2, 1)
   CovidSim.Branch(5, 3, 0.8, (0, 0), "nil", "recovered")
   CovidSim.Branch(5, 7, 0.2, (3, 3), "nil", "sick")
(2, 2)
   CovidSim.Branch(6, 6, 1.0, (3, 2), "mild", "mild")
(2, 3)
   CovidSim.Branch(7, 7, 0.85, (3, 3), "sick", "sick")
   CovidSim.Branch(7, 8, 0.15, (3, 4), "sick", "severe")
(3, 2)
   CovidSim.Branch(6, 3, 1.0, (0, 0), "mild", "recovered")
(3, 3)
   CovidSim.Branch(7, 3, 0.8, (0, 0), "sick", "recovered")
   CovidSim.Branch(7, 7, 0.1, (5, 3), "sick", "sick")
   CovidSim.Branch(7, 8, 0.1, (4, 4), "sick", "severe")
(3, 4)
   CovidSim.Branch(8, 3, 0.45, (0, 0), "severe", "recovered")
   CovidSim.Branch(8, 8, 0.5, (4, 4), "severe", "severe")
   CovidSim.Branch(8, 4, 0.05, (0, 5), "severe", "dead")
(4, 4)
   CovidSim.Branch(8, 3, 0.85, (0, 0), "severe", "recovered")
   CovidSim.Branch(8, 8, 0.1, (5, 4), "severe", "severe")
   CovidSim.Branch(8, 4, 0.05, (0, 5), "severe", "dead")
(5, 3)
   CovidSim.Branch(7, 3, 0.9, (0, 0), "sick", "recovered")
   CovidSim.Branch(7, 4, 0.1, (0, 5), "sick", "dead")
(5, 4)
   CovidSim.Branch(8, 3, 0.6, (0, 0), "severe", "recovered")
   CovidSim.Branch(8, 4, 0.4, (0, 5), "severe", "dead")
=#

function walktree(dt, top)
    done = []
    not_done = [[top]]
    cnt = 1

    while !isempty(not_done)
        cnt += 1
        currentpath = popfirst!(not_done)
        endnode = currentpath[end]
        for br in dt[endnode]
            if br.next[1] == 0
                push!(done, [currentpath..., br.next])
            else
                push!(not_done, [currentpath..., br.next])   # append without modifying currentpath
            end
        end
    end
    return done
end