import json

import magneticfielddb as mag
import simplesimdb as simplesim

# select magnetic field and R_0
params = mag.select( "COMPASS/compass_1X.json")
params["R_0"]=545.0

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
        "direction" : "centered", # forward makes a weird instability
        "nu_perp_n" : 1e-3,
        "nu_perp_u" : 1e-3,
        "nu_parallel_n" : 5e2
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
        "refine" : [12,12], # 12 may be better for conservation than 6
        "rk4eps" : 1e-6,
        "periodify" : False,
        "bc" : "along_field",
        "interpolation-method" : "linear-nearest"
    },
    "physical":
    {
        "mu" : -0.00027244371074816386,
        "tau" : 1.0,
        "beta" : 1e-4,
        "resistivity" : 1e-4,
        "epsilon_D" : 4.1458919332419e-05,
        "viscosity" : "value",
        "nu_parallel" : [ 3700, 114]

    },
    "timestepper":
    {
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
        "itstp" : 100,
        "maxout": 100,
        "compression": [2,2],
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

# nu_perp = 1e-3 is at the edge of the CFL condition for the velocity upwind correction
#Tabular format for simulations:
simulations = [
    # tau(0), resistivity(1), epsilon_D(2) , source_rate(3), deltaT(4), nu_e(5), nu_i(6), nu_p(7), nu_n(8), direction]# [1ms]
    #[0, 3e-4, 2.3936318213431424e-05, 1.6e-4,  100, 1233,  114, 1e-3,  5e2, "forward"], # 14700
    #[1, 3e-4, 2.3936318213431424e-05, 2e-4,    100, 1233,  114, 2e-3, 10e2, "forward"], # 14700
    #[0, 1e-4, 4.1458919332419e-05,    0.7e-4,  100, 3700,  114, 1e-3,  5e2, "forward"], # 19300
    [1, 1e-4, 4.1458919332419e-05,    1e-4,    100, 3700,  114, 1e-3,  5e2, "forward"], # 19300
    #[0, 3e-5, 7.569328442972544e-05,  0.5e-4,  100, 3700,  114, 2e-3, 10e2, "forward"], # 26100
    #[1, 3e-5, 7.569328442972544e-05,  0.5e-4,  100, 3700,  381, 1e-3,  5e2, "forward"], # 26100
    #[0, 1e-5, 0.00013110461385575586, 0.35e-4, 150, 3700,  114, 1e-3,  5e2, "forward"], # 34300
    #[1, 1e-5, 0.00013110461385575586, 0.35e-4, 150, 3700, 1143, 1.5e-3, 7.5e2, "forward"], # 34300
    #[0, 3e-6, 0.00023936318237861086, 0.30e-4, 200, 3700,  114, 1e-3,  5e2, "forward"], # 46400
    #[1, 3e-6, 0.00023936318237861086, 0.30e-4, 200, 3700, 3700, 1e-3, 10e2,"centered"], # 46400
    #[0, 1e-6, 0.0004145891932469323,  0.25e-4, 200, 3700,  114, 1e-3,  5e2], # 61000
    #[1, 1e-6, 0.0004145891932469323,  0.25e-4, 200, 3700, 3700, 1e-3, 10e2,"centered"]  # 61000
]
for values in  simulations:
    inputfile["physical"]["tau"] = values[0]
    inputfile["physical"]["resistivity"] = values[1]
    inputfile["physical"]["epsilon_D"] =  values[2]
    inputfile["source"]["rate"] =  values[3]
    inputfile["timestepper"]["deltaT"] =  values[4]
    inputfile["physical"]["nu_parallel"] =  [values[5], values[6]]
    inputfile["regularization"]["nu_perp_n"] = values[7]
    inputfile["regularization"]["nu_perp_u"] = values[7]
    inputfile["regularization"]["nu_parallel_n"] = values[8]
    inputfile["regularization"]["direction"] = values[9]
    print( inputfile["physical"])
    for i in range(0,1) : # set number of sims here
        if m.exists( inputfile,i) :
            print( "Simulation already run ", m.outfile( inputfile, i))
        else:
            print( "Run Simulation ", m.outfile( inputfile, i))
            m.create( inputfile, i, error="display")

#testfile = inputfile
#testfile["grid"]["Nx"] = 48
#testfile["grid"]["Ny"] = 80
#testfile["flags"] = ["symmetric"]
#
#testfile["output"]["type"] = "glfw"
#testfile["output"]["itstp"] = 1
#testfile["timestepper"]["deltaT"] = 1
#testfile["output"]["window"] = {"rows":6, "reduction" : 4, "width" :200, "height" : 400}
#with open( "test.json", 'w') as f:
#    inputstring = json.dump( inputfile, f,
#        sort_keys=True, ensure_ascii=True, indent=4)
