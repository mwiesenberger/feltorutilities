import json

import simplesimdb as simplesim

with open( "compass_1X-input.json", "r" ) as f:
    inputfile = json.load(f)

m = simplesim.Manager( directory="data", executable="./submit_job.sh", filetype="nc")

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


