### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ ef025bd0-787f-496d-9284-282daeb2c190
begin
	# using CovidSim_group
	push!(LOAD_PATH, joinpath(homedir(), "Dropbox/Online Coursework/Covid/ilm-src")) # using CovidSim_ilm

	try
		using Revise
	catch e 
		@warn(e.msg)
	end

end

# ╔═╡ a8021925-f200-42b5-b65a-b1ef9c0cc435
using TypedTables

# ╔═╡ d407275f-1e48-4e0a-a803-aa7152e533cb
using CovidSim_ilm

# ╔═╡ 5a3536dc-b402-4ab9-b741-d64723b9f766
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# ╔═╡ bac12bbe-2957-42c9-99a7-61f0eeae0d53
result_dict, env, series = run_a_sim(180, 38015, showr0=false, silent=true, spreadcases=[], runcases=[seed_1_6]);

# ╔═╡ df7ef52f-893c-4219-9fc8-eb8c096abff8
cumplot(series, 38015)

# ╔═╡ ee7b2803-27c1-4464-83ae-efbad2aaa8b6
65012/95626

# ╔═╡ 742d5523-218d-4571-aada-6470fd805488
keys(result_dict)


# ╔═╡ Cell order:
# ╠═ef025bd0-787f-496d-9284-282daeb2c190
# ╠═a8021925-f200-42b5-b65a-b1ef9c0cc435
# ╠═d407275f-1e48-4e0a-a803-aa7152e533cb
# ╠═5a3536dc-b402-4ab9-b741-d64723b9f766
# ╠═bac12bbe-2957-42c9-99a7-61f0eeae0d53
# ╠═df7ef52f-893c-4219-9fc8-eb8c096abff8
# ╠═ee7b2803-27c1-4464-83ae-efbad2aaa8b6
# ╠═742d5523-218d-4571-aada-6470fd805488
