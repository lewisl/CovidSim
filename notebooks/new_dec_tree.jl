# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# # New approach to decision trees for transition

# %%
using CovidSim_ilm

# %%
using StatsBase
using DelimitedFiles
using Distributions
using PrettyPrint
using JSON2

# %%
cd(joinpath(homedir(),"Dropbox/Online Coursework/Covid/ilm-src"))

# %%
ilmat = readdlm("../data/ilmtestdata.csv", ',',Int, header=true)[1]
ilmat = repeat(ilmat, 10_000)
refresh = copy(ilmat)

# %% [markdown]
# ### Current YAML Approach as of 4/28/2021

# %%
dectreefilename="../parameters/dec_tree_all_25.yml"

# %%
dectree_dict = setup_dt(dectreefilename)

# %%
dectree = dectree_dict["dt"]

# %%
display_tree(dectree)

# %% [markdown]
# ## Experiments with YAML 

# %%
dectreefilename="../parameters/dec_tree_all_25.csv"

# %%
using JSON
using Pkg.TOML
using YAML

# %%
dt, all_decpoints = CovidSim.setup_dt_old(dectreefilename)

# %%
newdict1 = Dict(string(k) => dt[1].tree[k] for k in keys(dt[1].tree))

# %%
keys(dt[1].tree)

# %%
JSON.json(dt[1].tree)

# %%
JSON.json(dt)

# %%
TOML.print(newdict1)

# %% jupyter={"source_hidden": true} tags=[]
# this is json from the YAML
jsontxt = """
[
  {
    "nodes": [
      {
        "branches": [
          {
            "pr": 1.0, 
            "tocond": 7, 
            "id": "to 7", 
            "next": "(0, 0)"
          }, 
          {
            "pr": 0.85, 
            "tocond": 7, 
            "id": "to 7", 
            "next": "(0, 0)"
          }, 
          {
            "pr": 0.15, 
            "tocond": 8, 
            "id": "to 8", 
            "next": "(4, 4)"
          }
        ], 
        "id": "(14, 6)"
      }
    ], 
    "agegrp": 1
  }, 
  {
    "nodes": [
      {
        "branches": [
          {
            "pr": 1.0, 
            "tocond": 7, 
            "id": "to 7", 
            "next": "(0, 0)"
          }
        ], 
        "id": "(14, 6)"
      }, 
      {
        "branches": [
          {
            "pr": 0.85, 
            "tocond": 7, 
            "id": "to 7", 
            "next": "(0, 0)"
          }, 
          {
            "pr": 0.15, 
            "tocond": 8, 
            "id": "to 8", 
            "next": "(4, 4)"
          }
        ], 
        "id": "(14, 7)"
      }
    ], 
    "agegrp": 2
  }
]
"""

# %%
tst = JSON.parse(jsontxt)

# %%
keys(tst[2])

# %%
tst[2]["nodes"]

# %%
tst[2]["agegrp"]

# %%
# this is hand authored json
jsonhand = """
{
"1":{
      "(5,5)":[{"tocond":5,"pr":0.2,"next":[9,5]},
                {"tocond":6,"pr":0.65,"next":[9,6]},
                {"tocond":7,"pr":0.15,"next":[9,7]}],
      "(9,5)":[{"tocond":3,"pr":0.8,"next":[0,0]},
                {"tocond":7,"pr":0.2,"next":[14,7]}],
      "(9,6)":[{"tocond":6,"pr":1.0,"next":[14,6]}],
      "(9,7)":[{"tocond":7,"pr":0.85,"next":[14,7]},
                {"tocond":8,"pr":0.15,"next":[14,8]}],
      "(14,6)":[{"tocond":3,"pr":1.0,"next":[0,0]}],
      "(14,7)":[{"tocond":3,"pr":0.85,"next":[0,0]},
                {"tocond":8,"pr":0.15,"next":[19,8]}],
      "(14,8)":[{"tocond":3,"pr":0.45,"next":[0,0]},
                {"tocond":8,"pr":0.5,"next":[19,8]},
                {"tocond":4,"pr":0.05,"next":[0,5]}],
      "(19,8)":[{"tocond":3,"pr":0.9,"next":[0,0]},
                {"tocond":4,"pr":0.1,"next":[0,5]}]
    },
 "2":{
      "(5,5)":[{"tocond":5,"pr":0.2,"next":[9,5]},
                {"tocond":6,"pr":0.65,"next":[9,6]},
                {"tocond":7,"pr":0.15,"next":[9,7]}],
      "(9,5)":[{"tocond":3,"pr":0.8,"next":[0,0]},
                {"tocond":7,"pr":0.2,"next":[14,7]}],
      "(9,6)":[{"tocond":6,"pr":1.0,"next":[14,6]}],
      "(9,7)":[{"tocond":7,"pr":0.85,"next":[14,7]},
                {"tocond":8,"pr":0.15,"next":[14,8]}],
      "(14,6)":[{"tocond":3,"pr":1.0,"next":[0,0]}],
      "(14,7)":[{"tocond":3,"pr":0.8,"next":[10,0]},
                {"tocond":8,"pr":0.2,"next":[19,8]}],
      "(14,8)":[{"tocond":3,"pr":0.35,"next":[0,0]},
                {"tocond":8,"pr":0.55,"next":[19,8]},
                {"tocond":4,"pr":0.1,"next":[0,5]}],
      "(19,8)":[{"tocond":3,"pr":0.85,"next":[0,0]},
                {"tocond":4,"pr":0.15,"next":[0,5]}]
      },
 "3":{
      "(5,5)":[{"tocond":5,"pr":0.2,"next":[9,5]},
                {"tocond":6,"pr":0.6,"next":[9,6]},
                {"tocond":7,"pr":0.2,"next":[9,7]}],
      "(9,5)":[{"tocond":3,"pr":0.7,"next":[0,0]},
                {"tocond":7,"pr":0.3,"next":[14,7]}],
      "(9,6)":[{"tocond":6,"pr":1.0,"next":[14,6]}],
      "(9,7)":[{"tocond":7,"pr":0.85,"next":[14,7]},
                {"tocond":8,"pr":0.15,"next":[14,8]}],
      "(14,6)":[{"tocond":3,"pr":1.0,"next":[0,0]}],
      "(14,7)":[{"tocond":3,"pr":0.8,"next":[10,0]},
                {"tocond":8,"pr":0.2,"next":[19,8]}],
      "(14,8)":[{"tocond":3,"pr":0.35,"next":[0,0]},
                {"tocond":8,"pr":0.55,"next":[19,8]},
                {"tocond":4,"pr":0.1,"next":[0,5]}],
      "(19,8)":[{"tocond":3,"pr":0.85,"next":[0,0]},
                {"tocond":4,"pr":0.15,"next":[0,5]}]
      },
 "4":{
      "(5,5)":[{"tocond":5,"pr":0.2,"next":[9,5]},
                {"tocond":6,"pr":0.65,"next":[9,6]},
                {"tocond":7,"pr":0.15,"next":[9,7]}],
      "(9,5)":[{"tocond":3,"pr":0.8,"next":[0,0]},
                {"tocond":7,"pr":0.2,"next":[14,7]}],
      "(9,6)":[{"tocond":6,"pr":1.0,"next":[14,6]}],
      "(9,7)":[{"tocond":7,"pr":0.85,"next":[14,7]},
                {"tocond":8,"pr":0.15,"next":[14,8]}],
      "(14,6)":[{"tocond":3,"pr":1.0,"next":[0,0]}],
      "(14,7)":[{"tocond":3,"pr":0.8,"next":[10,0]},
                {"tocond":8,"pr":0.2,"next":[19,8]}],
      "(14,8)":[{"tocond":3,"pr":0.35,"next":[0,0]},
                {"tocond":8,"pr":0.55,"next":[19,8]},
                {"tocond":4,"pr":0.1,"next":[0,5]}],
      "(19,8)":[{"tocond":3,"pr":0.85,"next":[0,0]},
                {"tocond":4,"pr":0.15,"next":[0,5]}]
      },
 "5":{
      "(5,5)":[{"tocond":5,"pr":0.2,"next":[9,5]},
                {"tocond":6,"pr":0.65,"next":[9,6]},
                {"tocond":7,"pr":0.15,"next":[9,7]}],
      "(9,5)":[{"tocond":3,"pr":0.8,"next":[0,0]},
                {"tocond":7,"pr":0.2,"next":[14,7]}],
      "(9,6)":[{"tocond":6,"pr":1.0,"next":[14,6]}],
      "(9,7)":[{"tocond":7,"pr":0.85,"next":[14,7]},
                {"tocond":8,"pr":0.15,"next":[14,8]}],
      "(14,6)":[{"tocond":3,"pr":1.0,"next":[0,0]}],
      "(14,7)":[{"tocond":3,"pr":0.8,"next":[10,0]},
                {"tocond":8,"pr":0.2,"next":[19,8]}],
      "(14,8)":[{"tocond":3,"pr":0.35,"next":[0,0]},
                {"tocond":8,"pr":0.55,"next":[19,8]},
                {"tocond":4,"pr":0.1,"next":[0,5]}],
      "(19,8)":[{"tocond":3,"pr":0.85,"next":[0,0]},
                {"tocond":4,"pr":0.15,"next":[0,5]}]
      }

}
"""

# %% jupyter={"outputs_hidden": true} tags=[]
hand = JSON.parse(jsonhand)

# %%
hand["1"]

# %%
hand["1"]["(5,5)"]

# %%
using YAML

# %%
YAML.write_file("../parameters/test-output.yml", hand)

# %% jupyter={"source_hidden": true} tags=[]
newyaml = """
4:
  (9,6):
    - tocond: 6
      next:
        - 14
        - 6
      pr: 1.0
  (9,5):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.8
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.2
  (14,6):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 1.0
  (9,7):
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.85
    - tocond: 8
      next:
        - 14
        - 8
      pr: 0.15
  (5,5):
    - tocond: 5
      next:
        - 9
        - 5
      pr: 0.2
    - tocond: 6
      next:
        - 9
        - 6
      pr: 0.65
    - tocond: 7
      next:
        - 9
        - 7
      pr: 0.15
  (14,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.35
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.55
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.1
  (14,7):
    - tocond: 3
      next:
        - 10
        - 0
      pr: 0.8
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.2
  (19,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.85
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.15
1:
  (9,6):
    - tocond: 6
      next:
        - 14
        - 6
      pr: 1.0
  (9,5):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.8
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.2
  (14,6):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 1.0
  (9,7):
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.85
    - tocond: 8
      next:
        - 14
        - 8
      pr: 0.15
  (5,5):
    - tocond: 5
      next:
        - 9
        - 5
      pr: 0.2
    - tocond: 6
      next:
        - 9
        - 6
      pr: 0.65
    - tocond: 7
      next:
        - 9
        - 7
      pr: 0.15
  (14,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.45
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.5
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.05
  (14,7):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.85
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.15
  (19,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.9
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.1
5:
  (9,6):
    - tocond: 6
      next:
        - 14
        - 6
      pr: 1.0
  (9,5):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.8
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.2
  (14,6):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 1.0
  (9,7):
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.85
    - tocond: 8
      next:
        - 14
        - 8
      pr: 0.15
  (5,5):
    - tocond: 5
      next:
        - 9
        - 5
      pr: 0.2
    - tocond: 6
      next:
        - 9
        - 6
      pr: 0.65
    - tocond: 7
      next:
        - 9
        - 7
      pr: 0.15
  (14,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.35
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.55
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.1
  (14,7):
    - tocond: 3
      next:
        - 10
        - 0
      pr: 0.8
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.2
  (19,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.85
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.15
2:
  (9,6):
    - tocond: 6
      next:
        - 14
        - 6
      pr: 1.0
  (9,5):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.8
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.2
  (14,6):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 1.0
  (9,7):
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.85
    - tocond: 8
      next:
        - 14
        - 8
      pr: 0.15
  (5,5):
    - tocond: 5
      next:
        - 9
        - 5
      pr: 0.2
    - tocond: 6
      next:
        - 9
        - 6
      pr: 0.65
    - tocond: 7
      next:
        - 9
        - 7
      pr: 0.15
  (14,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.35
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.55
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.1
  (14,7):
    - tocond: 3
      next:
        - 10
        - 0
      pr: 0.8
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.2
  (19,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.85
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.15
3:
  (9,6):
    - tocond: 6
      next:
        - 14
        - 6
      pr: 1.0
  (9,5):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.7
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.3
  (14,6):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 1.0
  (9,7):
    - tocond: 7
      next:
        - 14
        - 7
      pr: 0.85
    - tocond: 8
      next:
        - 14
        - 8
      pr: 0.15
  (5,5):
    - tocond: 5
      next:
        - 9
        - 5
      pr: 0.2
    - tocond: 6
      next:
        - 9
        - 6
      pr: 0.6
    - tocond: 7
      next:
        - 9
        - 7
      pr: 0.2
  (14,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.35
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.55
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.1
  (14,7):
    - tocond: 3
      next:
        - 10
        - 0
      pr: 0.8
    - tocond: 8
      next:
        - 19
        - 8
      pr: 0.2
  (19,8):
    - tocond: 3
      next:
        - 0
        - 0
      pr: 0.85
    - tocond: 4
      next:
        - 0
        - 5
      pr: 0.15
"""

# %%
roundtrip_yaml = YAML.load(newyaml)

# %%
roundtrip_yaml[1]

# %%
roundtrip_yaml[1]["(9,7)"]

# %% jupyter={"source_hidden": true} tags=[]
dense_yaml = """
4:
  (9,6):
    - {tocond: 6, next: [14,6], pr: 1.0}
  (9,5):
    - {tocond: 3, next: [0,0], pr: 0.8}
    - {tocond: 7, next: [14,7], pr: 0.2}
  (14,6):
    - {tocond: 3, next: [0,0], pr: 1.0}
  (9,7):
    - {tocond: 7, next: [14, 7], pr: 0.85}
    - {tocond: 8, next: [14, 8], pr: 0.15}
  (5,5):
    - {tocond: 5, next: [9, 5], pr: 0.2}
    - {tocond: 6, next: [9, 6], pr: 0.65}
    - {tocond: 7, next: [9, 7], pr: 0.15}
  (14,8):
    - {tocond: 3, next: [0,0], pr: 0.35}
    - {tocond: 8, next: [19, 8], pr: 0.55}
    - {tocond: 4, next: [0,5], pr: 0.1}
  (14,7):
    - {tocond: 3, next: [0, 0], pr: 0.8}
    - {tocond: 8, next: [19, 8], pr: 0.2}
  (19,8):
    - {tocond: 3, next: [0,0], pr: 0.85}
    - {tocond: 4, next: [0,5], pr: 0.15}
1:
  (9,6):
    - {tocond: 6, next: [14, 6], pr: 1.0}
  (9,5):
    - {tocond: 3, next: [0,0], pr: 0.8}
    - {tocond: 7, next: [14, 7], pr: 0.2}
  (14,6):
    - {tocond: 3, next: [0,0], pr: 1.0}
  (9,7):
    - {tocond: 7, next: [14, 7], pr: 0.85}
    - {tocond: 8, next: [14, 8], pr: 0.15}
  (5,5):
    - {tocond: 5, next: [9, 5], pr: 0.2}
    - {tocond: 6, next: [9, 6], pr: 0.65}
    - {tocond: 7, next: [9, 7], pr: 0.15}
  (14,8):
    - {tocond: 3, next: [0,0], pr: 0.45}
    - {tocond: 8, next: [19, 8], pr: 0.5}
    - {tocond: 4, next: [0,5], pr: 0.05}
  (14,7):
    - {tocond: 3, next: [0,0], pr: 0.85}
    - {tocond: 8, next: [19, 8], pr: 0.15}
  (19,8):
    - {tocond: 3, next: [0,0], pr: 0.9}
    - {tocond: 4, next: [0,5], pr: 0.1}
5:
  (9,6):
    - {tocond: 6, next: [14,6], pr: 1.0}
  (9,5):
    - {tocond: 3, next: [0,0], pr: 0.8}
    - {tocond: 7, next: [14,7], pr: 0.2}
  (14,6):
    - {tocond: 3, next: [0,0], pr: 1.0}
  (9,7):
    - {tocond: 7, next: [14, 7], pr: 0.85}
    - {tocond: 8, next: [14, 8], pr: 0.15}
  (5,5):
    - {tocond: 5, next: [9, 5], pr: 0.2}
    - {tocond: 6, next: [9, 6], pr: 0.65}
    - {tocond: 7, next: [9, 7], pr: 0.15}
  (14,8):
    - {tocond: 3, next: [0,0], pr: 0.35}
    - {tocond: 8, next: [19, 8], pr: 0.55}
    - {tocond: 4, next: [0,5], pr: 0.1}
  (14,7):
    - {tocond: 3, next: [0, 0], pr: 0.8}
    - {tocond: 8, next: [19, 8], pr: 0.2}
  (19,8):
    - {tocond: 3, next: [0,0], pr: 0.85}
    - {tocond: 4, next: [0,5], pr: 0.15}
2:
  (9,6):
    - {tocond: 6, next: [14,6], pr: 1.0}
  (9,5):
    - {tocond: 3, next: [0,0], pr: 0.8}
    - {tocond: 7, next: [14, 7], pr: 0.2}
  (14,6):
    - {tocond: 3, next: [0,0], pr: 1.0}
  (9,7):
    - {tocond: 7, next: [14, 7], pr: 0.85}
    - {tocond: 8, next: [14, 7], pr: 0.15}
  (5,5):
    - {tocond: 5, next: [9, 5], pr: 0.2}
    - {tocond: 6, next: [9, 6], pr: 0.65}
    - {tocond: 7, next: [9, 7], pr: 0.15}
  (14,8):
    - {tocond: 3, next: [0,0], pr: 0.35}
    - {tocond: 8, next: [19, 8], pr: 0.55}
    - {tocond: 4, next: [0,5], pr: 0.1}
  (14,7):
    - {tocond: 3, next: [0, 0], pr: 0.8}
    - {tocond: 8, next: [19, 8], pr: 0.2}
  (19,8):
    - {tocond: 3, next: [0,0], pr: 0.85}
    - {tocond: 4, next: [0,5], pr: 0.15}
3:
  (9,6):
    - {tocond: 6, next: [14,6], pr: 1.0}
  (9,5):
    - {tocond: 3, next: [0,0], pr: 0.8}
    - {tocond: 7, next: [14,7], pr: 0.2}
  (14,6):
    - {tocond: 3, next: [0,0], pr: 1.0}
  (9,7):
    - {tocond: 7, next: [14, 7], pr: 0.85}
    - {tocond: 8, next: [14, 8], pr: 0.15}
  (5,5):
    - {tocond: 5, next: [9, 5], pr: 0.2}
    - {tocond: 6, next: [9, 6], pr: 0.65}
    - {tocond: 7, next: [9, 7], pr: 0.15}
  (14,8):
    - {tocond: 3, next: [0,0], pr: 0.35}
    - {tocond: 8, next: [19, 8], pr: 0.55}
    - {tocond: 4, next: [0,5], pr: 0.1}
  (14,7):
    - {tocond: 3, next: [0, 0], pr: 0.8}
    - {tocond: 8, next: [19, 8], pr: 0.2}
  (19,8):
    - {tocond: 3, next: [0,0], pr: 0.85}
    - {tocond: 4, next: [0,5], pr: 0.15}
"""

# %%
dense_literal_yaml = YAML.load(dense_yaml)

# %% jupyter={"source_hidden": true} tags=[]
hand_yaml = """
3:
  (5,5):
    - {tocond: 5, next: [9, 5], pr: 0.2}
    - {tocond: 6, next: [9, 6], pr: 0.7}
    - {tocond: 7, next: [9, 7], pr: 0.1}
  (9,5):
    - {tocond: 3, next: [0,0], pr: 0.85}
    - {tocond: 7, next: [14, 7], pr: 0.15}
  (9,6):
    - {tocond: 6, next: [14,6], pr: 1.0}
  (9,7):
    - {tocond: 7, next: [14, 7], pr: 0.9}
    - {tocond: 8, next: [14, 8], pr: 0.1}
  (14,6):
    - {tocond: 3, next: [0,0], pr: 1.0}
  (14,7):
    - {tocond: 3, next: [0, 0], pr: 0.83}
    - {tocond: 8, next: [19, 7], pr: 0.1}
    - {tocond: 8, next: [19, 8], pr: 0.07}
  (14,8):
    - {tocond: 3, next: [0,0], pr: 0.474}
    - {tocond: 8, next: [19, 8], pr: 0.514}
    - {tocond: 4, next: [0,5], pr: 0.012}
  (19,8):
    - {tocond: 3, next: [0,0], pr: 0.922}
    - {tocond: 8, next: [0,5], pr: 0.072}
    - {tocond: 4, next: [0,5], pr: 0.006}
  (25,7):
    - {tocond: 3, next: [0,0], pr: 0.964}
    - {tocond: 4, next: [0,5], pr: 0.036}
  (25,8):
    - {tocond: 3, next: [0,0], pr: 0.964}
    - {tocond: 4, next: [0,5], pr: 0.036}
"""

# %%
hand_test = YAML.load(hand_yaml)

# %%
dense_test = YAML.load_file("../parameters/dec_tree_all_25.yml")

# %%
decpoints = Dict{Int,Array{Int, 1}}()
for i in 1:5
    decpoints[i] = unique([k[1] for k in keys(dense_test[i])])
end
decpoints

# %%
display_tree(dense_test[1])

# %%
dense_test[4]

# %%
length(dense_test)

# %%
CovidSim.sanity_test_all(dense_test)

# %%
std_file = "../parameters/dec_tree_all_25.yml"

# %% jupyter={"outputs_hidden": true}
treedict, decpoints = setup_dt(std_file)

# %%
treedict

# %%
treedict[1][[9,5]]

# %%
decpoints

# %%
