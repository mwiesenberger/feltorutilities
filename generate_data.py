import simplesimdb as simplesim
import feltorutilities as fp
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
        "Nx" : 192, # 12*16 (max is about 312 on 1 GPU)
        "Ny" : 352, # 22*16 (max is about 576 on 1 GPU)
        "Nz" : 32,
        "scaleR" : [1.45,1.25], # 2.7
        "scaleZ" : [2.6,2.25]   # 4.85
    },
    "boundary" :
    {
        "wall" :
        {
            "type" : "sol_pfr",
            "boundary" : [1.15, 0.95],
            "alpha" : [0.1,0.05],
            "penalization" : 1.0,
            "modify-B" : True,
            "penalize-rhs" : False
        },
        "sheath" :
        {
            "type" : "insulating",
            "boundary" : 0.25,
            "alpha" : 0.15,
            "penalization" : 1.0,
            "penalize-rhs" : True
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
    "init" : { "type" : "zero" },
    "source" :
    {
        "type" : "influx",
        "rate" : 2e-3,
        "ntilde":
        {
            "type" : "turbulence",
            "amplitude" : 0.01,
            "rk4eps" : 1e-6,
            "refine" : [5,5],
            "revolutions" : 1,
            "sigma_z" : 0.5
        },
        "profile":
        {
            "type" : "aligned",
            "amplitude" : 1.0
        },
        "damping":
        {
            "type" : "aligned",
            "alpha" : 0.2,
            "boundary" : 1.0
        }
    },
    "timestepper":
    {
        "tableau" : "TVB-3-3",
        "dt" : 1e-2
    },
    "regularization" :
    {
        "order" : 2,
        "direction" : "centered",
        "nu_perp_n" : 2e-5,
        "nu_perp_u" : 2e-5,
    },
    "elliptic":
    {
        "stages": 3,
        "eps_pol" : [1e-6, 100, 100],
        "eps_gamma" : 1e-7,
        "eps_ampere" : 1e-7,
        "direction" : "centered",
        "jumpfactor" : 1.0
    },
    "FCI":
    {
        "refine" : [2,2],
        "rk4eps" : 1e-6,
        "periodify" : True,
        "bc" : "along_field"
    },
    "physical":
    {
        "mu" : -0.00027244371074816386,
        "tau" : 1.0,
        "beta" : 5e-5,
        "resistivity" : 5e-5,
        "epsilon_D" : 8.291783866749523e-05,
        "viscosity" : "braginskii"
    },
    "output" :
    {
        "type" : "netcdf",
        "inner_loop" : 5,
        "itstp" : 500,
        "maxout": 50,
        "compression": [2,2]
    }
}

m = simplesim.Manager( directory="data", executable="submit_job.sh", filetype="nc")
#
m.create( inputfile, 0)

#test = simplesim.Repeater( "./feltor.sh", "test.json", "test.nc")
#testfile = inputfile
#testfile["grid"]["Nx"] = 48
#testfile["grid"]["Ny"] = 88
#testfile["flags"] = ["symmetric"]
#
#testfile["output"]["type"] = "glfw"
#testfile["output"]["inner_loop"] = 2
#testfile["output"]["itstp"] = 1
#testfile["output"]["window"] = {"rows":6, "reduction" : 4, "width" :200, "height" : 400}
#test.run( testfile)


