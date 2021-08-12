import json
import simplesimdb as simplesim
import magneticfielddb as mag
import numpy as np


# select magnetic field and R_0
params = mag.select( "COMPASS/compass_1X.json")
params["R_0"]=545.0
# select input
inputfile = {
    "magnetic_field" :
    {
        "curvmode" : "toroidal",
        "input" : "params",
        "params" : params
    },
    "grid" :
    {
        "n" : 3,
        "Nx" : 192, # 12*16 or 18*16 (max is about 312 on 1 GPU)
        "Ny" : 336, # 21*16 or 31*16 (max is about 576 on 1 GPU)
        "Nz" : 32,
        "scaleR" : [1.45,1.25], # 2.7
        "scaleZ" : [2.4,2.25]   # 4.65
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
            "boundary" : 3/16, # large boundary seems to stabilize
            "alpha" : 3/16-1/32, # should also be large
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
    "timestepper":
    {
        "type" : "adaptive",
        "tableau" : "Bogacki-Shampine-4-2-3",
        "rtol": 1e-5,
        "atol" : 1e-9,
        "output-mode" : "equidistant",
        "Tend" : 1e5,
        "reject-limit" : 10
    },
    "regularization" :
    {
        "order" : 2,
        "direction" : "forward",
        "nu_perp_n" : 5e-3,
        "nu_perp_u" : 5e-3,
    },
    "advection" :
    {
        "slope-limiter": "none"
    },
    "elliptic":
    {
        "stages": 3,
        "eps_pol" : [1e-6, 100, 100],
        "eps_gamma" : 1e-7,
        "eps_ampere" : 1e-7,
        "direction" : "forward",
        "jumpfactor" : 1.0
    },
    "FCI":
    {
        "refine" : [5,5],
        "rk4eps" : 1e-6,
        "periodify" : False,
        "bc" : "along_field",
        "interpolation-method" : "linear"
    },
    "physical":
    {
        "mu" : -0.00027244371074816386,
        "tau" : 1.0,
        "beta" : 1e-4,
        "resistivity" : 1e-4,
        "epsilon_D" : 4.1458919332419e-05,
        "viscosity" : "braginskii"
    },
    "output" :
    {
        "type" : "netcdf",
        "inner_loop" : 50,
        "itstp" : 50,
        "maxout": 2000,
        "compression": [2,2]
    }
}

#m = simplesim.Manager( directory="data", executable="./submit_job.sh", filetype="nc")
#
#m.create( inputfile, 0, error="display")

testfile = inputfile
testfile["grid"]["Nx"] = 48
testfile["grid"]["Ny"] = 80
testfile["flags"] = ["symmetric"]

testfile["output"]["type"] = "glfw"
testfile["output"]["inner_loop"] = 2
testfile["output"]["itstp"] = 1
testfile["output"]["window"] = {"rows":6, "reduction" : 4, "width" :200, "height" : 400}
with open( "test.json", 'w') as f:
    inputstring = json.dump( inputfile, f,
        sort_keys=True, ensure_ascii=True, indent=4)


