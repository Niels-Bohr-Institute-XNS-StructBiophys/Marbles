import numpy as np
import sys
import os

col = int(sys.argv[1])

mine = np.loadtxt("out_check")
ref  = np.loadtxt("../bead_modeling_flat/out_check")
thresh = 1e-12

if( col == -1 ):
    mine_f = mine#[:,0]
    ref_f = ref#[:,0]
else:
    mine_f = mine[:,col]
    ref_f = ref[:,col]

err = []

check = ( len(mine_f) == len(ref_f) )
if( not check ):
    print("## LEN PROBLEM!! ##")
    print( f"REF: {len(ref_f)}, MINE: {len(mine_f)}" )
    exit()
else:
    print("Same len confirmed! ", len(ref_f), len(mine_f))

dim = mine_f.shape[0]
count = 0
count2 = 0
tr = thresh * ref_f.max()
tm = thresh * mine_f.max()

for i in range( dim ):

    if( ref_f[i] == 0 and mine_f[i] == 0 ):
        rel_err = 0
    elif( ref_f[i] < tr and mine_f[i] < tm ):
        rel_err = 0
        count2 += 1
    else:
        rel_err = np.fabs( (mine_f[i] - ref_f[i])/ref_f[i] )

    err.append( rel_err )
    if( rel_err > 5e-4 ):
        print( f"CHECK! Rel err: {rel_err},  {mine_f[i]} != {ref_f[i]} at {i}" )
        count += 1

print( "Average rel err: ", np.mean(err) )
print( "Percentage of errors: ",  count / dim )
print( "Negligible numbers: ", count2 / dim )
