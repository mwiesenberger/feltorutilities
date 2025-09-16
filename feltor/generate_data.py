import json

import simplesimdb as simplesim

##########################################
# These simulations reproduce the data in
# "Wiesenberger and Held, Plasma Phys. Control. Fusion 66 065003 (2024)"
##########################################

with open( "compass_1X-input.json", "r" ) as f:
    inputfile = json.load(f)

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
