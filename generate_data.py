import json
import simplesimdb as simplesim
import magneticfielddb as mag
import numpy as np


# select magnetic field and R_0
params = mag.select( "COMPASS/compass_1X.json")
params["R_0"]=545.0
tau = 0.0 # 1.0

inputfile = {
    "magnetic_field" :
    {
        "curvmode" : "toroidal",
        "input" : "params",
        "params" : params
    },
    "grid" :
    {
        "n" : 4,
        "Nx" : 192, # 12*16 or 18*16 (max is about 3*312 on 1 GPU)
        "Ny" : 336, # 21*16 or 31*16 (max is about 3*576 on 1 GPU)
        "Nz" : 32,
        "scaleR" : [1.45,1.25], # 2.7
        "scaleZ" : [2.475,2.25]   # 4.725
    },
    "advection":
    {
        "slope-limiter" : "none"
    },
    "boundary" :
    {
        "wall" :
        {
            "type" : "sol_pfr",
            "boundary" : [1.15, 0.97],
            "alpha" : [0.10,0.10],
            "penalization" : 1e-2, # 1 seems to prevent good timestep in adaptive
            "modify-B" : True,
            "penalize-rhs" : True,
            "nwall" : 0.2,
            "uwall" : 0.0
        },
        "sheath" :
        {
            #"type" : "wall",
            "type" : "insulating",
            #"type" : "bohm",
            "boundary" : 7/32, # large boundary seems to stabilize
            "alpha" : 7/32-2/32, # should also be large
            "penalization" : 5.0, # larger runs better
            "penalize-rhs" : True,
            "coordinate" : "s",
            "max_angle" : 4
        },
        "bc" :
        {
            "density" : ["NEU", "NEU"],
            "velocity": ["NEU", "NEU"],
            "potential":["DIR", "DIR"],
            "aparallel":["NEU", "NEU"]
        }
    },
    "flags" : [],
    "init" :
    {
        "type" : "fields",
        "density" :
        {
            "type" : "ne",
            "ntilde":
            {
                "type" : "turbulence",
                "amplitude" : 1e-4,
                "revolutions" : 1,
                "parallel" : "gaussian",
                "sigma_z" : 0.5
            },
            "profile":
            {
                "type" : "aligned",
                "npeak" : 8.5,
                "nsep" : 1.0,
                "background" : 0.2
            },
            "damping":
            {
                "type" : "alignedPFR",
                "alpha" : [0.1,0.03],
                "boundary" : [1.15, 0.97]
            }
        },
        "potential": { "type" : "zero_pol"},
        "velocity":
        {
            #"type" : "zero",
            "type" : "ui",
            "profile" : "linear_cs"
        },
        "aparallel": { "type" : "zero"}
    },
    "source" :
    {
        "minne" : 0.2,  # minne seems to rescue the sol ue instability
        "minrate" : 1.0, #
        "minalpha" : 0.05, # smaller seems better
        #"minne" : 0.,  # without sheath no need for density forcing
        "type" : "influx",
        "rate" : 1e-4,
        "profile":
        {
            "type" : "aligned",
            "npeak" : 1.0,
            "nsep" : 1.0/8.5,
            "background" : 0.0,
        },
        "ntilde" : {"type" : "zero"},
        "damping":
        {
            "type" : "alignedX",
            "alpha" : 0.2,
            "boundary" : 0.55,
        }
    },
    "regularization" :
    {
        "order" : 2,
        "direction" : "forward",
        "nu_perp_n" : 1e-4,
        "nu_perp_u" : 1e-4,
        "nu_parallel_n" : 1e3
    },
    "elliptic":
    {
        "stages": 5,
        "eps_pol" : [1e-6, 0.5, 0.5, 0.5, 0.5], # 0.5 prevents outliers on stage 0
        "eps_gamma" : 1e-8, # gamma inverts much faster than pol
        "eps_ampere" : 1e-8, # ampere inverts much faster than pol
        "direction" : "forward", # centered can make oscillations
        "jumpfactor" : 1.
    },
    "FCI":
    {
        "refine" : [6,6], # 12 may be better for conservation than 6
        "rk4eps" : 1e-6,
        "periodify" : False,
        "bc" : "along_field",
        "interpolation-method" : "linear-nearest"
    },
    "physical":
    {
        "mu" : -0.00027244371074816386,
        "tau" : tau,
        "beta" : 1e-4,
        "resistivity" : 1e-4,
        "epsilon_D" : 4.1458919332419e-05,
        "viscosity" : "value",
        "nu_parallel" : [ 3700, 114]

    },
    "timestepper":
    {
        #"type" : "multistep-imex",
        #"tableau" : "ImEx-TVB-3-3",
        ##"dt" : 0.1,
        #"dt" : 0.02,
        #"eps_time" : 1e-8
        #"type" : "adaptive",
        "tableau" : "Bogacki-Shampine-4-2-3",
        "type" : "adaptive",
        "rtol": 1e-5,
        "atol" : 1e-6,
        "output-mode" : "deltaT",
        "deltaT" : 100,
        "reject-limit" : 10
    },
    "output" :
    {
        "type" : "netcdf",
        "itstp" : 200,
        "maxout": 1000,
        "compression": [1,1],
        "equations":[
            "Basic",
            "Mass-conserv",
            "Energy-theorem",
            "Toroidal-momentum",
            "Parallel-momentum",
            "Zonal-Flow-Energy",
            "COCE"
        ]
    }
}

m = simplesim.Manager( directory="data", executable="./submit_job.sh", filetype="nc")

eps_map = { # resistivity -> [ epsilon_D , source_rate, deltaT, previous name]# [1ms]
    3e-4: [2.3936318213431424e-05, 2e-4,    100], # 14667
    1e-4: [4.1458919332419e-05,    1e-4,    100], # 19303
    3e-5: [7.569328442972544e-05,  0.5e-4,  100], # 26082
    1e-5: [0.00013110461385575586, 0.35e-4, 150], # 34326
    3e-6: [0.00023936318237861086, 0.30e-4, 200], # 46382
    1e-6: [0.0004145891932469323,  0.25e-4, 200]  # 61042
}
for key in  eps_map:
    inputfile["physical"]["resistivity"] = key
    inputfile["physical"]["epsilon_D"] =  eps_map[key][0]
    inputfile["source"]["rate"] =  eps_map[key][1]
    inputfile["timestepper"]["deltaT"] =  eps_map[key][2]
    print( inputfile["physical"])
    for i in range(0,6) : # set number of sims here
        if m.exists( inputfile,i) :
            print( "Simulation already run ", m.outfile( inputfile, i))
        else:
            print( "Run Simulation ", m.outfile( inputfile, i))
            m.create( inputfile, i, error="display")

testfile = inputfile
testfile["grid"]["Nx"] = 48
testfile["grid"]["Ny"] = 80
testfile["flags"] = ["symmetric"]

testfile["output"]["type"] = "glfw"
testfile["output"]["itstp"] = 1
testfile["timestepper"]["deltaT"] = 1
testfile["output"]["window"] = {"rows":6, "reduction" : 4, "width" :200, "height" : 400}
with open( "test.json", 'w') as f:
    inputstring = json.dump( inputfile, f,
        sort_keys=True, ensure_ascii=True, indent=4)


