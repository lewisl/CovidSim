# Seattle

plot(1:n, series[seattle.fips][:cum][1+adjdays:n+adjdays,map2series.totinfected[6]], label="Simulation",dpi=150, size=(400,260), tickfontsize=6)
plot!(1:n, tricities[:,:sea_infected], label="Reported Cases",dpi=150, size=(400,260), tickfontsize=6)
title!("King County Simulation vs. Reported Cases", titlefontsize=8)
xlabel!("Days: Feb. 1 to April 30", guidefontsize=8)



# New York

plot(1:122, vcat(zeros(abs(adjdays)),series[newyork.fips][:cum][1:n+adjdays,map2series.totinfected[6]]), label="Simulation",dpi=150, size=(400,260), tickfontsize=6)
# plot(1:122, series[newyork.fips][:cum][1+adjdays:n+adjdays,map2series.totinfected[6]], label="Simulation",dpi=150, size=(400,260), tickfontsize=6)
plot!(1:122, tricities[:,:nyc_infected], label="Reported Cases",dpi=150, size=(400,260), tickfontsize=7)
title!("New York City Simulation vs. Reported Cases", titlefontsize=7)
xlabel!("Days: Feb. 1 to April 30", guidefontsize=8)
