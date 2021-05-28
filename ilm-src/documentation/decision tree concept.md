#### Decision Tree Concept

In an SEIR simulation of the COVID-19 outbreak (Susceptible, Exposed, Infected, Removed (Recovery or Death), I use a decision tree to transition people who have gotten infected from nil (asymptomatic), to mild, to sick, to severe, to recovered, (sadly) to dead. The simulated folks don't go through all conditions.  After a number of days of being infectious they reach decision nodes every few days that probabilistically send them to a different condition over the course of 25 days. The leaf nodes of the tree are outcomes of either dead or recovered.  It's essential that all the probabilities of all leaves total to 1.0.  Each age group has its own tree so I sanity check all of them.

I thought I *should* be able to do this recursively, but my brain is more loopy than recursive.  I had done the text book recursion cases of factorial and following the branches of a binary tree, with the classic paths down left and right branches. But, the illness transition decision trees are more ragged.  They can have any number of branches (though always 1 to 3) and each branch can have different depths (longest is 6).  For 2 days I beat myself up, got impatient, and gave up.  I built the paths by hand (the human brain can see the paths almost immediately) and then fed the paths as input to the sanity check. I put the problem down to get more important stuff done.

Tonight, I took it on again.  In half an hour in < 20 lines of code I solved it. It's not recursive, but it is conceptually recursive (I did learn something...).  The trickiest thing wasn't the recursive-lite approach. I realized I need to append to the current path without modifying the current path. In ML-like languages we'd do this with "cons", which creates a pair.  In ML, one abhors mutating a variable so you can only call with a new value (the argument) or return a new value. And recursion performs the loop. In Julia, we love to mutate data structures in place for efficiency.  There is append!, but not append. But, it's easy to use vcat create a new array that is an append--or "cons" of the old array (the head) and its new tail. A while loop is a bit like recursion:  we stop when the conditions have been met that we are done, not by a fixed number of iterations. The while loop runs 17 times for a tree with 22 leaves--and thus 22 branches. Each time through we finish at least one branch of a given depth.

```julia
function walktree(dt, top)
    done = []
    todo = [[top]]

    while !isempty(todo)
        currentpath = popfirst!(todo)
        endnode = currentpath[end]
        for br in dt[endnode]
            if br.next[1] == 0
                push!(done, vcat(currentpath, br.next))  # append without modifying currentpath
            else
                push!(todo, vcat(currentpath, br.next))   
            end
        end
    end
    return done
end
```

Instead of pushing things onto the stack, I push them onto the todo list. I pop an item from the todo list and there are 2 outcomes: either the next node down a branch is a leaf and that path is done; or the next node will lead to another and the path is one more node longer and gets pushed onto the todo list.

This is conceptually like recursion. Each addition to a path is tested:  In one case, the end condition is met and we push it onto done (in recursion we'd return a static result up the call chain and NOT make another recursive call). In the other case, we have not reached a leaf and we push the path extended by one node onto todo (in recursion, we'd call ourself again with the extended, incomplete branch).

It's sort of inverted recursion because I extend the path by a node before deciding what to do with it. In recursion, the completion of the paths comes up as the recursive function returns static results instead of making another call. But, it's short, simple and relatively obvious because at each step there are only 2 choices.  I am sure it's not very efficient and there are lots of allocations, but recursion wouldn't be that efficient and there would be lots of allocations pushed onto the stack. But, the trees are small--see the one below.

It's a bit tricky because a tree is a bit of a gnarly nested data structure: 

- a tree is a dict of nodes;
- each node is an array of branches (largest number of branches is 3);
- a branch is a struct, which includes the next node or a sentinel node that says "leaf"
- leaf nodes are encoded as either (0,0) or (0,5): in textbook examples there is a boolean "isleaf"--either way a simple test tells us whether we are at the end of a path or not.

Now, I can quickly tell for all the trees if the branch endpoints all add up to 1 and get the expected value proportions for recovering or dying.  This was really satisfying, if slightly morbid.

Here is what a tree looks like:

```julia
#=
(1, 1)  # node and indented array of branches
   CovidSim.Branch(5, 5, 0.2, (2, 1), "nil", "nil")  # from, to, probability, next node, from name, to name
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

```


Here is the output for a valid collection of 5 trees, one for each age group. The leaf node probabilities all sum to one (3rd column is total probs for "recover" and 4th column is total probs for "dead").

```
5×4 Array{Float64,2}:
 1.0  1.0  0.999511  0.000488522
 2.0  1.0  0.999332  0.000668336
 3.0  1.0  0.99616   0.0038404
 4.0  1.0  0.980863  0.0191366
 5.0  1.0  0.850242  0.149758
```