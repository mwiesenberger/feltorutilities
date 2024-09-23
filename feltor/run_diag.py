# It is safe to run several instances of this script as simplesimdb recognises which netcdf files exist
# and feltordiag creates the netcdf file at program start almost instantly

import json

import simplesimdb as simplesim

data = simplesim.Manager( directory="/m100_scratch/userexternal/mwiesenb/feltor3d", filetype = "nc")
diag = simplesim.Manager(
    executable="./feltordiag.sh",
    directory="/m100_scratch/userexternal/mwiesenb/diag",
    filetype = "nc")

config = {
        # for feltordiag
        "n": 3,
        "Npsi": 64,
        "Neta": 640,
        "Kphi": 10,
        "fx_0" : 0.125,
        "fsa" : "toroidal-average", # or "convoluted-toroidal-average"
        "diagnostics":
        [
            "fsa",
            # "fsa2d",
            # "cta2d",
            # "cta2dX",
            # "fluc2d", # these can also be computed from cta2d, fsa2d
            "ifs",
            "std_fsa",
            "ifs_lcfs",
            "ifs_norm"
        ],
        # for interpolate in 3d
        "fine-grid-factor" : 2,
        "cta-interpolation" : "dg",
        "x-grid-interpolation" : "dg",
        "time-reduction-factor" : 10
    }

with open( "config.json", 'w') as f:
    json.dump( config, f,
        sort_keys=True, ensure_ascii=True, indent=4)

content = data.table()
for pp in content :
    for i in range( data.count(pp)):
        print( "Data ", data.outfile(pp), i)
        diag.create( pp,i,error="display",stdout="display")
