### R0 Interpretation

R0 in the model is *not* a direct input as it is in some differential equation models. Instead, it is the result of the spread! function, which models interaction of people by group (agegroup, disease condition) and the probability of virus transmission.  So, it falls out from the simulated interaction, much as R0 is in the real world can't be observed directly, but is the result of the transmissibility of the virus and the social patterns of people's interactions.

In the function ```how_many_contacts!```, the contact_factors determine how many people are contacted by infected people we refer to as "spreaders." The touch_factors determine whether a contact is "consequential" and distributes the contacts across the composition of the population.  This includes receiving people who are not susceptible, either because they have become (partially) immune through recovery or because they are already infected. This reduces the transmission of the virus because the spreader has a lower probability of contacting another person who is susceptible. This is the social dynamic of R0: as the body of people who acquire immunity by recovery or by vaccination (currently not modeled) grows larger, the effective R0 goes down because fewer people get the disease.

The combination of spreaders making contacts and recipients receptive to a touch is based on pseudo-random sampling from 2 distributions:
- The contact_factor is the scale input to a gamma distribution sample which produces a distribution of number of contacts. A larger contact factor results in wider dispersion of the sample output--e.g., a longer tail to the right while the mode is close to the left--typically between 0 and 2.
- The touch_factor is the "probability of success" input to a pseudo-random sample from the binomial distribution that determines the success of the touch.

The biological transmissibility of the disease, whether by respiratory droplets or aerosols or acquisition from physical surfaces, is still not fully understood. The model simplistically uses a multiplicative probability based on the number of days the spreader has had the disease and the age group of the recipient ("touched" person).

The inputs below are judgment inputs that are "sanity checked" to produce R0 values in the ranges that epidemiologists have calcuated (or assumed?) based on observations of the early stages of the spread of the Coronavirus in several countries. The tendency is that younger people, who believe that they are not infected, even though in the simulation--and in reality--they may be infectious asymptomatic or with very mild symptoms, will make more contacts. Younger recipients without observable symptoms will, likewise, be more accessible.  Older people and people in the disease categories of sick and severe will make and accept a lower rate of contact.

While R0 is very hard to "observe" in the wild, it is a useful diagnostic to interpret how aggressive a model is in spreading the disease. The model provides an R0 simulator, which calculates the R0 from a stage 1 cohort of spreaders; traces the spreaders through all 25 days of their disease; and simply adds up how many people get the virus from this stage 1 cohort.  Those infected "disappear" out of the simulations so that their spreading is not added to the spreading of the initial cohort. Why a cohort?  Why not just 1 person? We want to have a dispersion of people across age groups and disease conditions because we don't know what the characteristics of a single "patient 0" would be. Thus, we calculate an average across the entire range of characteristics in roughly the distribution that they occur.




Here are the probabilities input to a binomial sample
for the defaults in the model.

rows: age group of the recipients
columns: lag of the spreaders

 0.0   0.0   0.0   0.0   0.0
 0.03  0.12  0.12  0.15  0.165
 0.07  0.28  0.28  0.35  0.385
 0.08  0.32  0.32  0.4   0.44
 0.09  0.36  0.36  0.45  0.495
 0.09  0.36  0.36  0.45  0.495
 0.08  0.32  0.32  0.4   0.44
 0.07  0.28  0.28  0.35  0.385
 0.06  0.24  0.24  0.3   0.33
 0.05  0.2   0.2   0.25  0.275
 0.04  0.16  0.16  0.2   0.22
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.03  0.12  0.12  0.15  0.165
 0.01  0.04  0.04  0.05  0.055
